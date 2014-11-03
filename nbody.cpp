#include <cstdlib>
#include <ctime>
#include <cmath>
#include <iostream>

#include "bitmap_image.hpp"

void
run_simulation(unsigned int n_particles,
                unsigned int n_steps,
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
    unsigned int max_initspeed = 30;
    unsigned int img_width     = 500;
    unsigned int img_height    = 500;


#ifndef VISUAL
    if (3 > argc) {
        print_help(argv[0]);
        return 1;
    }

    n_particles = atoi(argv[1]);
    n_steps     = atoi(argv[2]);
#else
    if (6 > argc) {
        print_help(argv[0]);
        return 1;
    }

    n_particles   = atoi(argv[1]);
    n_steps       = atoi(argv[2]);
    max_initspeed = atoi(argv[3]);
    img_width     = atoi(argv[4]);
    img_height    = atoi(argv[5]);
#endif

    run_simulation(n_particles,
                    n_steps,
                    max_initspeed,
                    img_width,
                    img_height);

    return 0;
}

void
run_simulation(unsigned int n_particles,
                unsigned int n_steps,
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

    init_array(x,    0,                        img_width,     n_particles);
    init_array(y,    0,                        img_height,    n_particles);
    init_array(xvel, -((float) max_initspeed), max_initspeed, n_particles);
    init_array(yvel, -((float) max_initspeed), max_initspeed, n_particles);
    init_array(mass, 0,                        100,           n_particles);

    for (unsigned int s = 0; s < n_steps; ++s) {
        for (unsigned int p = 0; p < n_particles; ++p) {
            xnew[p] = x[p] + xvel[p];
            ynew[p] = y[p] + yvel[p];
        }

        std::swap(x, xnew);
        std::swap(y, ynew);

#ifdef VISUAL
        std::cout << '\r' << "step: " << s + 1 << '/' << n_steps;
        save_image(x, y, n_particles, img_width, img_height, s);
        char a;
        std::cin >> a;
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
    static bitmap_image image(width, height);
    static image_drawer drawer(image);

    // set background to black
    image.set_all_channels(0);
    drawer.pen_width(3);
    drawer.pen_color(255, 255, 255);

    for (unsigned int p = 0; p < n_particles; ++p) {
        if (0 <= x[p] && x[p] <= width
         && 0 <= y[p] && y[p] <= height) {
            drawer.plot_pixel(x[p], y[p]);
        }
    }

    image.save_image("output.bmp");
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
    std::cout << " max_initspeed width height";
#endif
    std::cout << std::endl;
}

