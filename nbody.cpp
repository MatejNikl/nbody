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

    float * x    = new float[n_particles];
    float * y    = new float[n_particles];
    float * xnew = new float[n_particles];
    float * ynew = new float[n_particles];
    float * xvel = new float[n_particles];
    float * yvel = new float[n_particles];
    float * mass = new float[n_particles];

    init_array(x,    0,              img_width,     n_particles);
    init_array(y,    0,              img_height,    n_particles);
    init_array(xvel, -max_initspeed, max_initspeed, n_particles);
    init_array(yvel, -max_initspeed, max_initspeed, n_particles);
    init_array(mass, 0,              max_initmass,  n_particles);

    xvel[0] = yvel[0] = 0;
    x[0] =  img_width / 2;
    y[0] = img_height / 2;
    mass[0] = -1e6;

    for (unsigned int s = 0; s < n_steps; ++s) {
        for (unsigned int i = 0; i < n_particles; ++i) {
            float ax = 0.0f;
            float ay = 0.0f;

            for (unsigned int j = 0; j < n_particles; ++j) {
                float dx = x[j] - x[i];
                float dy = y[j] - y[i];
                float invr = 1.0/sqrt(dx*dx + dy*dy + 10);
                float f = mass[j]*invr*invr*invr;
                ax += f*dx; /* accumulate the acceleration from gravitational attraction */
                ay += f*dy;
            }

            xnew[i] = x[i] + time_step*xvel[i] + 0.5f*time_step*time_step*ax; /* update position of particle "i" */
            ynew[i] = y[i] + time_step*yvel[i] + 0.5f*time_step*time_step*ay;
            xvel[i] += time_step*ax; /* update velocity of particle "i" */
            yvel[i] += time_step*ay;

            if (xnew[i] < 0.0f) {
                xnew[i] = 0.0f;
                xvel[i] = std::fabs(xvel[i]) / 2;
            } else if (xnew[i] > img_width) {
                xnew[i] = img_width;
                xvel[i] = -std::fabs(xvel[i]) / 2;
            }

            if (ynew[i] < 0.0f) {
                ynew[i] = 0.0f;
                yvel[i] = std::fabs(yvel[i]) / 2;
            } else if (ynew[i] > img_height) {
                ynew[i] = img_height;
                yvel[i] = -std::fabs(yvel[i]) / 2;
            }
        }

        std::swap(x, xnew);
        std::swap(y, ynew);

#ifdef VISUAL
        std::cout << '\r' << "step: " << s + 1 << '/' << n_steps << std::flush;
        save_image(x, y, n_particles, img_width, img_height, s);
#endif
    }

#ifdef VISUAL
    std::cout << std::endl;
#endif

    delete [] x;
    delete [] xnew;
    delete [] y;
    delete [] ynew;
    delete [] xvel;
    delete [] yvel;
    delete [] mass;
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
