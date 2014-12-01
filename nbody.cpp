#include <cstdlib>
#include <ctime>
#include <cmath>
#include <csignal>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <chrono>

#include <xmmintrin.h>

#include "bitmap_image.hpp"

typedef float v4sf __attribute__((vector_size (16)));

struct NBodySettings
{
    unsigned int n_particles;
    unsigned int n_steps;
    unsigned int img_width;
    unsigned int img_height;
    unsigned int plot_every;
    float        time_step;
    float        max_initspeed;
    float        max_initmass;
    float        max_initcharge;
    float        min_initspeed;
    float        min_initmass;
    float        min_initcharge;

    NBodySettings()
    : n_particles(1000),
      n_steps(1000),
      img_width(500),
      img_height(500),
      plot_every(10),
      time_step(1.0f),
      max_initspeed(30.0f),
      max_initmass(1000.0f),
      max_initcharge(1000.0f),
      min_initspeed(1e-6f), //cannot be zero
      min_initmass(0.0f),
      min_initcharge(-1000.0f)
    {
    }

    friend std::ostream & operator<<(std::ostream & os, const NBodySettings & s);
};

void
run_simulation(const NBodySettings & s);

void
save_image(float * x,
           float * y,
           float * m,
           float * q,
           const NBodySettings & s,
           unsigned int seq);

void
signal_handler(int signum);

void
init_array(float * arr,
           float min,
           float max,
           unsigned int n);

std::ostream &
operator<<(std::ostream & os, const NBodySettings & s);

void
print_status(unsigned int step, const NBodySettings & s);

void
print_help(const char * runcmnd);

volatile bool g_interrupted = false;

int
main(int argc, char *argv[])
{
    NBodySettings s;

    srand(time(NULL));
    signal(SIGINT, signal_handler);

    if(argc < 3) {
        print_help(*argv);
        return EXIT_FAILURE;
    }

    // Skip execution command
    --argc, ++argv;

    // Required arguments
    s.n_particles = std::atoi(*argv++), --argc;
    s.n_steps     = std::atoi(*argv++), --argc;

    // Optional arguments
    if(0 < argc)
        s.img_width      = std::atoi(*argv++), --argc;
    if(0 < argc)
        s.img_height     = std::atoi(*argv++), --argc;
    if(0 < argc)
        s.plot_every     = std::atoi(*argv++), --argc;
    if(0 < argc)
        s.time_step      = std::atof(*argv++), --argc;
    if(0 < argc)
        s.max_initspeed  = std::atof(*argv++), --argc;
    if(0 < argc)
        s.max_initmass   = std::atof(*argv++), --argc;
    if(0 < argc)
        s.max_initcharge = std::atof(*argv++), --argc;
    if(0 < argc)
        s.min_initspeed  = std::atof(*argv++), --argc;
    if(0 < argc)
        s.min_initmass   = std::atof(*argv++), --argc;
    if(0 < argc)
        s.min_initcharge = std::atof(*argv++), --argc;

    run_simulation(s);

    return EXIT_SUCCESS;
}

void
run_simulation(const NBodySettings & s)
{
    std::cout << "Running simulation with:" << std::endl;
    std::cout << s;

    float * x  = new float[s.n_particles];
    float * y  = new float[s.n_particles];
    float * xn = new float[s.n_particles];
    float * yn = new float[s.n_particles];
    float * vx = new float[s.n_particles];
    float * vy = new float[s.n_particles];
    float * m  = new float[s.n_particles];
    float * q  = new float[s.n_particles];

    init_array(x,                0,      s.img_width, s.n_particles);
    init_array(y,                0,     s.img_height, s.n_particles);
    init_array(vx, s.min_initspeed,  s.max_initspeed, s.n_particles);
    init_array(vy, s.min_initspeed,  s.max_initspeed, s.n_particles);
    init_array(m,   s.min_initmass,   s.max_initmass, s.n_particles);
    init_array(q, s.min_initcharge, s.max_initcharge, s.n_particles);

    auto begin = std::chrono::steady_clock::now();

    const float half = 0.5f;
    const float dt = s.time_step;
    const float dt2 = dt * dt;

    const v4sf vhalf = { half, half, half, half };
    const v4sf vdt = { dt, dt, dt, dt };
    const v4sf vdt2 = vdt * vdt;

    unsigned int step;
    for (step = 0; step < s.n_steps && !g_interrupted; ++step) {
#       pragma omp parallel for schedule(static) default(none) \
            firstprivate(x, y, xn, yn, vx, vy, m, q) shared(s)
        for (unsigned int i = 0; i < s.n_particles - 3; i += 4) {
            v4sf ax = { 0.0f, 0.0f, 0.0f, 0.0f };
            v4sf ay = { 0.0f, 0.0f, 0.0f, 0.0f };

            const v4sf xi = *(v4sf *)(x + i);
            const v4sf yi = *(v4sf *)(y + i);
            const v4sf mi = *(v4sf *)(m + i);
            const v4sf qi = *(v4sf *)(q + i);

            for (unsigned int j = 0; j < s.n_particles; ++j) {
                const v4sf xj = { x[j], x[j], x[j], x[j] };
                const v4sf yj = { y[j], y[j], y[j], y[j] };
                const v4sf mj = { m[j], m[j], m[j], m[j] };
                const v4sf qj = { q[j], q[j], q[j], q[j] };

                const v4sf dx = xj - xi;
                const v4sf dy = yj - yi;

                const v4sf invr = __builtin_ia32_rsqrtps( dx * dx + dy * dy + vhalf );
                const v4sf coef = (mj - qi * qj / mi) * invr * invr * invr;

                /* accumulate the acceleration from gravitational attraction */
                ax += coef * dx;
                ay += coef * dy;
            }

            const v4sf vxi = *(v4sf *)(vx + i);
            const v4sf vyi = *(v4sf *)(vy + i);

            /* update position of particle "i" */
            *(v4sf *)(xn + i) = xi + vxi * vdt + vhalf * ax * vdt2;
            *(v4sf *)(yn + i) = yi + vyi * vdt + vhalf * ay * vdt2;

            /* update velocity of particle "i" */
            *(v4sf *)(vx + i) = vxi + ax * vdt;
            *(v4sf *)(vy + i) = vyi + ay * vdt;
        }

        for (unsigned int i = s.n_particles & ~0x3; i < s.n_particles; ++i) {
            float ax = 0.0f;
            float ay = 0.0f;

            for (unsigned int j = 0; j < s.n_particles; ++j) {
                const float dx = x[j] - x[i];
                const float dy = y[j] - y[i];
                const float invr = 1.0f / std::sqrt(dx * dx + dy * dy + half);
                const float coef = (m[j] - q[i] * q[j] / m[i]) * invr * invr * invr;

                /* accumulate the acceleration from gravitational attraction */
                ax += coef * dx;
                ay += coef * dy;
            }

            /* update position of particle "i" */
            xn[i] = x[i] + vx[i] * dt + half * ax * dt2;
            yn[i] = y[i] + vy[i] * dt + half * ay * dt2;

            /* update velocity of particle "i" */
            vx[i] += ax * dt;
            vy[i] += ay * dt;
        }

#ifdef VISUAL
#       pragma omp parallel for schedule(static) default(none) \
            firstprivate(xn, yn, vx, vy) shared(s)
        for (unsigned int i = 0; i < s.n_particles; ++i) {
            if (xn[i] < 0) {
                xn[i] = -xn[i];
                vx[i] = 0.5f * std::fabs(vx[i]);
            } else if (xn[i] > s.img_width - 1) {
                xn[i] = 2 * (s.img_width - 1) - xn[i];
                vx[i] = -0.5f * std::fabs(vx[i]);
            }

            if (yn[i] < 0) {
                yn[i] = -yn[i];
                vy[i] = 0.5f * std::fabs(vy[i]);
            } else if (yn[i] > s.img_height - 1) {
                yn[i] = 2 * (s.img_height - 1) - yn[i];
                vy[i] = -0.5f * std::fabs(vy[i]);
            }
        }

        if (step % s.plot_every == 0) {
            save_image(x, y, m, q, s, step / s.plot_every);
        }

        print_status(step, s);
#endif

        std::swap(x, xn);
        std::swap(y, yn);
    }

#ifdef VISUAL
    if (step % s.plot_every == 0) {
        save_image(x, y, m, q, s, step / s.plot_every);
    }
    print_status(step, s);
    std::cout << std::endl;
#else
    std::cout << "Processed: " << step << " steps" << std::endl;
#endif

    auto end = std::chrono::steady_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
    std::cout << "Elapsed time: "
              << std::fixed << std::setprecision(3) << elapsed / 1000.0
              << " s" << std::endl;

    delete [] x;
    delete [] y;
    delete [] xn;
    delete [] yn;
    delete [] vx;
    delete [] vy;
    delete [] m;
    delete [] q;
}

void
save_image(float * x,
           float * y,
           float * m,
           float * q,
           const NBodySettings & s,
           unsigned int seq)
{
    bitmap_image image(s.img_width, s.img_height);
    image_drawer drawer(image);

    float max_val = std::max(std::max(std::fabs(s.max_initcharge),
                                      std::fabs(s.min_initcharge)),
                             s.max_initmass);
    float norm = 255.0f / max_val;

    image.set_all_channels(0);

    for (unsigned int p = 0; p < s.n_particles; ++p) {
        float xr = std::round(x[p]);
        float yr = std::round(y[p]);

        if (0 <= xr && xr <= s.img_width  - 1
         && 0 <= yr && yr <= s.img_height - 1) {
            if (q[p] >= 0.0f) {
                drawer.pen_color(std::round(q[p] * norm),
                                 std::round(m[p] * norm),
                                 0);
            } else if (q[p] < 0.0f) {
                drawer.pen_color(0,
                                 std::round(m[p] * norm),
                                 std::round(-q[p] * norm));
            }

            drawer.plot_pixel(xr, yr);
        }
    }

    std::ostringstream oss;
    oss << "output_" << std::setfill('0') << std::setw(5) << seq << ".bmp";
    image.save_image(oss.str());
}

void
signal_handler(int signum)
{
    g_interrupted = signum; //to silence (set but) not used warning
}

void
init_array(float * arr,
           float min,
           float max,
           unsigned int n)
{
    if (min == max) {
        std::fill(arr, arr + n, min);
    } else {
        if(max < min)
            std::swap(min, max);

        const float norm = (max - min) / RAND_MAX;

        for (unsigned int i = 0; i < n; ++i) {
            arr[i] = min + norm * rand();
        }
    }
}

std::ostream &
operator<<(std::ostream & os, const NBodySettings & s)
{
    os << " n_particles=    " << s.n_particles    << std::endl
       << " n_steps=        " << s.n_steps        << std::endl;
#ifdef VISUAL
    os << " img_width=      " << s.img_width      << std::endl
       << " img_height=     " << s.img_height     << std::endl
       << " plot_every=     " << s.plot_every     << std::endl
       << " time_step=      " << s.time_step      << std::endl
       << " max_initspeed=  " << s.max_initspeed  << std::endl
       << " max_initmass=   " << s.max_initmass   << std::endl
       << " max_initcharge= " << s.max_initcharge << std::endl
       << " min_initspeed=  " << s.min_initspeed  << std::endl
       << " min_initmass=   " << s.min_initmass   << std::endl
       << " min_initcharge= " << s.min_initcharge << std::endl;
#endif
    return os;
}

void
print_status(unsigned int step, const NBodySettings & s)
{

    std::cout << '\r'
              << "step: " << step << '/' << s.n_steps << ' '
              << "plotted: " << 1 + step / s.plot_every << '/' << 1 + s.n_steps / s.plot_every
              << std::flush;
}


void
print_help(const char * runcmnd)
{
    std::cout << "Usage: " << runcmnd
              << " #particles #steps [width height plot_every time_step max_initspeed max_initmass max_initcharge min_initspeed min_initmass min_initcharge]" << std::endl
              << " defaults:" << std::endl
              << NBodySettings(); //print default values
}
