#include <cstdlib>
#include <ctime>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <limits>

#include "bitmap_image.hpp"

struct NBodySettings
{
    unsigned int n_particles;
    unsigned int n_steps;
    unsigned int img_width;
    unsigned int img_height;
    float        time_step;
    float        max_initspeed;
    float        max_initmass;
    float        max_initcharge;
    float        min_initspeed;
    float        min_initmass;
    float        min_initcharge;

    NBodySettings()
    : n_particles(1000u),
      n_steps(1000),
      img_width(500),
      img_height(500),
      time_step(1.0f),
      max_initspeed(30.0f),
      max_initmass(1000.0f),
      max_initcharge(1000.0f),
      min_initspeed(0.0f),
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
           const NBodySettings & s,
           unsigned int seq);

void
init_array(float * arr,
           float min,
           float max,
           unsigned int n);

void
print_help(const char * runcmnd);



int main(int argc, char *argv[])
{
    NBodySettings s;

    srand(time(NULL));

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

    vx[0] = vy[0] = 0;
    x[0] =  s.img_width / 2;
    y[0] = s.img_height / 2;
    m[0] = -1e6;

    float dt = s.time_step;
    for (unsigned int step = 0; step < s.n_steps; ++step) {
        for (unsigned int i = 0; i < s.n_particles; ++i) {
            float ax = 0.0f;
            float ay = 0.0f;

            for (unsigned int j = 0; j < s.n_particles; ++j) {
                float dx = x[j] - x[i];
                float dy = y[j] - y[i];
                float invr = 1.0f / std::sqrt(dx * dx + dy * dy + 1.0f);
                float coef = (m[j] - q[i] * q[j] / m[i]) * invr * invr * invr;

                ax += coef * dx; /* accumulate the acceleration from gravitational attraction */
                ay += coef * dy;
            }

            xn[i] = x[i] + vx[i] * dt + 0.5f * ax * dt * dt; /* update position of particle "i" */
            yn[i] = y[i] + vy[i] * dt + 0.5f * ay * dt * dt;
            vx[i] += ax * dt; /* update velocity of particle "i" */
            vy[i] += ay * dt;

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
        }

        std::swap(x, xn);
        std::swap(y, yn);

#ifdef VISUAL
        std::cout << '\r' << "step: " << step + 1 << '/' << s.n_steps << std::flush;
        save_image(x, y, s, step);
#endif
    }

#ifdef VISUAL
    std::cout << std::endl;
#endif

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
           const NBodySettings & s,
           unsigned int seq)
{
    static bitmap_image image(s.img_width, s.img_height);
    static image_drawer drawer(image);

    image.set_all_channels(0);
    drawer.pen_color(255, 255, 255);

    for (unsigned int p = 0; p < s.n_particles; ++p) {
        if (0 <= x[p] && x[p] <= s.img_width
         && 0 <= y[p] && y[p] <= s.img_height) {
            drawer.plot_pixel(x[p], y[p]);
        }
    }

    std::ostringstream oss;
    oss << "output_" << std::setfill('0') << std::setw(5) << seq << ".bmp";
    image.save_image(oss.str());
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
       << " n_steps=        " << s.n_steps        << std::endl
       << " img_width=      " << s.img_width      << std::endl
       << " img_height=     " << s.img_height     << std::endl
       << " time_step=      " << s.time_step      << std::endl
       << " max_initspeed=  " << s.max_initspeed  << std::endl
       << " max_initmass=   " << s.max_initmass   << std::endl
       << " max_initcharge= " << s.max_initcharge << std::endl
       << " min_initspeed=  " << s.min_initspeed  << std::endl
       << " min_initmass=   " << s.min_initmass   << std::endl
       << " min_initcharge= " << s.min_initcharge << std::endl;
#ifndef VISUAL
    os << "(non-visual mode)" << std::endl;
#endif
    return os;
}

void
print_help(const char * runcmnd)
{
    std::cout << "Usage: " << runcmnd
              << " #particles #steps [width height time_step min_initspeed min_initmass min_initcharge min_initspeed min_initmass min_initcharge]" << std::endl
              << " defaults:" << std::endl
              << NBodySettings(); //print default values
}
