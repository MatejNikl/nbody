#include <cstdlib>
#include <ctime>
#include <cmath>
#include <csignal>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <chrono>

#include "bitmap_image.hpp"

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

    /*
    vx[0] = vy[0] = 0;
    x[0] =  s.img_width / 2;
    y[0] = s.img_height / 2;
    m[0] = -1e6;
    */

    auto begin = std::chrono::steady_clock::now();


    float dt = s.time_step;
    unsigned int step;

    for (step = 0; step < s.n_steps && !g_interrupted; ++step) {
        for (unsigned int i = 0; i + 3 < s.n_particles; ++i) {
            float ax_0 = 0.0f;
            float ax_1 = 0.0f;
            float ax_2 = 0.0f;
            float ax_3 = 0.0f;

            float ay_0 = 0.0f;
            float ay_1 = 0.0f;
            float ay_2 = 0.0f;
            float ay_3 = 0.0f;

            for (unsigned int j = 0; j < s.n_particles; ++j) {
                float dx_0 = x[j] - x[i + 0];
                float dx_1 = x[j] - x[i + 1];
                float dx_2 = x[j] - x[i + 2];
                float dx_3 = x[j] - x[i + 3];

                float dy_0 = y[j] - y[i + 0];
                float dy_1 = y[j] - y[i + 1];
                float dy_2 = y[j] - y[i + 2];
                float dy_3 = y[j] - y[i + 3];

                float invr_0 = 1.0f / std::sqrt(dx_0 * dx_0 + dy_0 * dy_0 + 0.5f);
                float invr_1 = 1.0f / std::sqrt(dx_1 * dx_1 + dy_1 * dy_1 + 0.5f);
                float invr_2 = 1.0f / std::sqrt(dx_2 * dx_2 + dy_2 * dy_2 + 0.5f);
                float invr_3 = 1.0f / std::sqrt(dx_3 * dx_3 + dy_3 * dy_3 + 0.5f);

                float coef_0 = (m[j] - q[i + 0] * q[j] / m[i + 0]) * invr_0 * invr_0 * invr_0;
                float coef_1 = (m[j] - q[i + 1] * q[j] / m[i + 1]) * invr_1 * invr_1 * invr_1;
                float coef_2 = (m[j] - q[i + 2] * q[j] / m[i + 2]) * invr_2 * invr_2 * invr_2;
                float coef_3 = (m[j] - q[i + 3] * q[j] / m[i + 3]) * invr_3 * invr_3 * invr_3;

                ax_0 += coef_0 * dx_0; /* accumulate the acceleration from gravitational attraction */
                ax_1 += coef_1 * dx_1;
                ax_2 += coef_2 * dx_2;
                ax_3 += coef_3 * dx_3;

                ay_0 += coef_0 * dy_0;
                ay_1 += coef_1 * dy_1;
                ay_2 += coef_2 * dy_2;
                ay_3 += coef_3 * dy_3;
            }

            xn[i + 0] = x[i + 0] + vx[i + 0] * dt + 0.5f * ax_0 * dt * dt; /* update position of particle "i" */
            xn[i + 1] = x[i + 1] + vx[i + 1] * dt + 0.5f * ax_1 * dt * dt;
            xn[i + 2] = x[i + 2] + vx[i + 2] * dt + 0.5f * ax_2 * dt * dt;
            xn[i + 3] = x[i + 3] + vx[i + 3] * dt + 0.5f * ax_3 * dt * dt;

            yn[i + 0] = y[i + 0] + vy[i + 0] * dt + 0.5f * ay_0 * dt * dt;
            yn[i + 1] = y[i + 1] + vy[i + 1] * dt + 0.5f * ay_1 * dt * dt;
            yn[i + 2] = y[i + 2] + vy[i + 2] * dt + 0.5f * ay_2 * dt * dt;
            yn[i + 3] = y[i + 3] + vy[i + 3] * dt + 0.5f * ay_3 * dt * dt;

            vx[i + 0] += ax_0 * dt; /* update velocity of particle "i" */
            vx[i + 1] += ax_1 * dt;
            vx[i + 2] += ax_2 * dt;
            vx[i + 3] += ax_3 * dt;

            vy[i + 0] += ay_0 * dt;
            vy[i + 1] += ay_1 * dt;
            vy[i + 2] += ay_2 * dt;
            vy[i + 3] += ay_3 * dt;

#ifdef VISUAL
            if (xn[i] < 0) {
                xn[i] = -xn[i];
                vx[i] = 0.5f * std::fabs(vx[i]);
            } else if (xn[i] > s.img_width) {
                xn[i] = 2 * s.img_width - xn[i];
                vx[i] = -0.5f * std::fabs(vx[i]);
            }

            if (yn[i] < 0) {
                yn[i] = -yn[i];
                vy[i] = 0.5f * std::fabs(vy[i]);
            } else if (yn[i] > s.img_height) {
                yn[i] = 2 * s.img_height - yn[i];
                vy[i] = -0.5f * std::fabs(vy[i]);
            }
#endif
        }

        std::swap(x, xn);
        std::swap(y, yn);

#ifdef VISUAL
        std::cout << '\r'
                  << "step: " << 1 + step << '/' << s.n_steps << ' '
                  << "plotted: " << std::ceil(((float) step) / s.plot_every) << '/' << std::ceil(((float) s.n_steps) / s.plot_every)
                  << std::flush;

        if (step % s.plot_every == 0) {
            save_image(x, y, m, q, s, step / s.plot_every);
        }
#endif
    }

#ifdef VISUAL
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
}

void
save_image(float * x,
           float * y,
           float * m,
           float * q,
           const NBodySettings & s,
           unsigned int seq)
{
    static const int pen_width = 2;
    bitmap_image image(s.img_width + pen_width * 2, s.img_height + pen_width * 2);
    image_drawer drawer(image);

    float max_val = std::max(std::max(std::fabs(s.max_initcharge),
                                      std::fabs(s.min_initcharge)),
                             s.max_initmass);
    float norm = 255.0f / max_val;

    image.set_all_channels(0);
    drawer.pen_width(pen_width);

    for (unsigned int p = 0; p < s.n_particles; ++p) {

        if (0 <= x[p] && x[p] <= s.img_width
         && 0 <= y[p] && y[p] <= s.img_height) {
            if (q[p] >= 0.0f) {
                drawer.pen_color(std::round(q[p] * norm),
                                 std::round(m[p] * norm),
                                 0);
            } else if (q[p] < 0.0f) {
                drawer.pen_color(0,
                                 std::round(m[p] * norm),
                                 std::round(q[p] * norm));
            }

            drawer.plot_pen_pixel(std::round(x[p] + pen_width),
                                  std::round(y[p] + pen_width));
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
print_help(const char * runcmnd)
{
    std::cout << "Usage: " << runcmnd
              << " #particles #steps [width height plot_every time_step max_initspeed max_initmass max_initcharge min_initspeed min_initmass min_initcharge]" << std::endl
              << " defaults:" << std::endl
              << NBodySettings(); //print default values
}
