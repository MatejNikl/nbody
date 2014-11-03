#include <cstdlib>
#include <ctime>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <limits>

#include "bitmap_image.hpp"

void
run_simulation(unsigned int n_particles,
               unsigned int n_steps,
               unsigned int img_width,
               unsigned int img_height,
               float time_step,
               float max_initspeed,
               float max_initmass);

void
save_image(float * x,
           float * y,
           unsigned int n_particles,
           unsigned int width,
           unsigned int height,
           unsigned int seq);

void
init_array(float * arr,
           float min,
           float max,
           unsigned int n);

void
print_help(const char * runcmd);

int main(int argc, char *argv[])
{
    unsigned int n_particles;
    unsigned int n_steps;
    unsigned int img_width     = 500;
    unsigned int img_height    = 500;
    float        time_step     = 1.0f;
    float        max_initspeed = 30.0f;
    float        max_initmass  = 1000.0f;

    srand(time(NULL));

    if(argc < 3) {
        print_help(*argv);
        return EXIT_FAILURE;
    }

    // Skip execution command
    --argc, ++argv;

    // Required arguments
    n_particles = std::atoi(*argv++), --argc;
    n_steps     = std::atoi(*argv++), --argc;

    // Optional arguments
    if(0 < argc)
        img_width     = std::atoi(*argv++), --argc;
    if(0 < argc)
        img_height    = std::atoi(*argv++), --argc;
    if(0 < argc)
        time_step     = std::atof(*argv++), --argc;
    if(0 < argc)
        max_initspeed = std::atof(*argv++), --argc;
    if(0 < argc)
        max_initmass  = std::atof(*argv++), --argc;

    run_simulation(n_particles,
                   n_steps,
                   img_width,
                   img_height,
                   time_step,
                   max_initspeed,
                   max_initmass);

    return EXIT_SUCCESS;
}

void
run_simulation(unsigned int n_particles,
               unsigned int n_steps,
               unsigned int img_width,
               unsigned int img_height,
               float time_step,
               float max_initspeed,
               float max_initmass)
{
    std::cout << "Running simulation with:"          << std::endl
              << " n_particles=   " << n_particles   << std::endl
              << " n_steps=       " << n_steps       << std::endl
              << " img_width=     " << img_width     << std::endl
              << " img_height=    " << img_height    << std::endl
              << " time_step=     " << time_step     << std::endl
              << " max_initspeed= " << max_initspeed << std::endl
              << " max_initmass=  " << max_initmass  << std::endl;

    float * x  = new float[n_particles];
    float * y  = new float[n_particles];
    float * xn = new float[n_particles];
    float * yn = new float[n_particles];
    float * vx = new float[n_particles];
    float * vy = new float[n_particles];
    float * m  = new float[n_particles];

    init_array(x,  0,              img_width,     n_particles);
    init_array(y,  0,              img_height,    n_particles);
    init_array(vx, -max_initspeed, max_initspeed, n_particles);
    init_array(vy, -max_initspeed, max_initspeed, n_particles);
    init_array(m,  0,              max_initmass,  n_particles);

    vx[0] = vy[0] = 0;
    x[0] =  img_width / 2;
    y[0] = img_height / 2;
    m[0] = -1e6;

    for (unsigned int s = 0; s < n_steps; ++s) {
        for (unsigned int i = 0; i < n_particles; ++i) {
            float ax = 0.0f;
            float ay = 0.0f;

            for (unsigned int j = 0; j < n_particles; ++j) {
                float dx = x[j] - x[i];
                float dy = y[j] - y[i];
                float invr = 1.0f / sqrt(dx * dx + dy * dy + 1e-03f);
                float coef = m[j] * invr * invr * invr;

                ax += coef * dx; /* accumulate the acceleration from gravitational attraction */
                ay += coef * dy;
            }

            xn[i] = x[i] + vx[i] * time_step + 0.5f * ax * time_step * time_step; /* update position of particle "i" */
            yn[i] = y[i] + vy[i] * time_step + 0.5f * ay * time_step * time_step;
            vx[i] += ax * time_step; /* update velocity of particle "i" */
            vy[i] += ay * time_step;

            if (xn[i] < 0) {
                xn[i] = -xn[i];
                vx[i] = 0.5f * std::fabs(vx[i]);
            } else if (xn[i] > img_width) {
                xn[i] = 2 * img_width - xn[i];
                vx[i] = -0.5f * std::fabs(vx[i]);
            }

            if (yn[i] < 0) {
                yn[i] = -yn[i];
                vy[i] = 0.5f * std::fabs(vy[i]);
            } else if (yn[i] > img_height) {
                yn[i] = 2 * img_height - yn[i];
                vy[i] = -0.5f * std::fabs(vy[i]);
            }
        }

        std::swap(x, xn);
        std::swap(y, yn);

#ifdef VISUAL
        std::cout << '\r' << "step: " << s + 1 << '/' << n_steps << std::flush;
        save_image(x, y, n_particles, img_width, img_height, s);
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
           unsigned int n_particles,
           unsigned int width,
           unsigned int height,
           unsigned int seq)
{
    static const unsigned int pen_width = 1;
    static bitmap_image image(width + 2 * pen_width, height + 2 * pen_width);
    static image_drawer drawer(image);

    image.set_all_channels(0);
    drawer.pen_width(pen_width);
    drawer.pen_color(255, 255, 255);

    for (unsigned int p = 0; p < n_particles; ++p) {
        if (0 <= x[p] && x[p] <= width
         && 0 <= y[p] && y[p] <= height) {
            drawer.plot_pen_pixel(x[p] + pen_width, y[p] + pen_width);
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

void
print_help(const char * runcmd)
{
    std::cout << "Usage: " << runcmd
              << " #particles #steps [width height time_step max_initspeed max_initmass]" << std::endl
              << " (default width=500 height=500 time_step=1 max_initspeed=30 max_initmass=1000)" << std::endl;
}
