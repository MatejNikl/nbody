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
               float        time_step,
               unsigned int max_initspeed,
               unsigned int img_width,
               unsigned int img_height);
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
    float        time_step     = 1;
    unsigned int max_initspeed = 30;
    unsigned int img_width     = 500;
    unsigned int img_height    = 500;


#ifndef VISUAL
    if (argc < 3) {
        print_help(argv[0]);
        return 1;
    }

    n_particles = atoi(argv[1]);
    n_steps     = atoi(argv[2]);
#else
    if (argc < 7) {
        print_help(argv[0]);
        return 1;
    }

    n_particles   = std::atoi(argv[1]);
    n_steps       = std::atoi(argv[2]);
    time_step     = std::atof(argv[3]);
    max_initspeed = std::atoi(argv[4]);
    img_width     = std::atoi(argv[5]);
    img_height    = std::atoi(argv[6]);
#endif

    run_simulation(n_particles,
                   n_steps,
                   time_step,
                   max_initspeed,
                   img_width,
                   img_height);

    return 0;
}

void
run_simulation(unsigned int n_particles,
               unsigned int n_steps,
               float        time_step,
               unsigned int max_initspeed,
               unsigned int img_width,
               unsigned int img_height)
{
    float * x    = new float[n_particles];
    float * xnew = new float[n_particles];
    float * y    = new float[n_particles];
    float * ynew = new float[n_particles];
    float * xvel = new float[n_particles];
    float * yvel = new float[n_particles];
    float * mass = new float[n_particles];

    srand(time(NULL));

    init_array(x,    0,                        img_width,     n_particles);
    init_array(y,    0,                        img_height,    n_particles);
    init_array(xvel, -((float) max_initspeed), max_initspeed, n_particles);
    init_array(yvel, -((float) max_initspeed), max_initspeed, n_particles);
    init_array(mass, 0,                        1000,          n_particles);

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
    if (min < max) {
        float norm = RAND_MAX / (max - min);

        for (unsigned int i = 0; i < n; ++i) {
            arr[i] = (rand() / norm) + min;
        }
    } else {
        std::fill(arr, arr + n, min);
    }
}

void
print_help(const char * runcmd)
{
    std::cout << "Usage: " << runcmd << " #particles #steps";
#ifdef VISUAL
    std::cout << " time_step max_initspeed width height";
#endif
    std::cout << std::endl;
}

