#include <cmath>

#include <string>
#include <iostream>

#include <dlfcn.h>
#include <xmmintrin.h>

#include "simulators.h"

typedef float v4sf __attribute__((vector_size (16)));

simulator_t
get_simulator(
    const char* name
    )
{
    void* h;
    simulator_t sim;

    std::string fname( "simulator_" );
    fname += name;

    /* GCC is whining about data pointer being converted
       to function pointer ... */
#   pragma GCC diagnostic push
#   pragma GCC diagnostic ignored "-Wpedantic"
    if( !(h = dlopen( NULL, RTLD_LAZY | RTLD_GLOBAL )) ||
        !(sim = (simulator_t)dlsym( h, fname.c_str() )) )
#   pragma GCC diagnostic pop
    {
        std::cerr << dlerror() << std::endl;
        return NULL;
    }

    return sim;
}

extern "C" unsigned int
simulator_naive(
    const simulator_conf_t* conf,
    const simulator_data_t* data
    )
{
    LOAD_CONF( conf );
    LOAD_DATA( data );

    const float half = 0.5f;
    const float dt2 = dt * dt;

    unsigned int step;
    for (step = 0; step < n_steps; ++step) {
        for (unsigned int i = 0; i < n_particles; ++i) {
            float ax = 0.0f;
            float ay = 0.0f;

            for (unsigned int j = 0; j < n_particles; ++j) {
                float dx = x[j] - x[i];
                float dy = y[j] - y[i];
                float invr = 1.0f / std::sqrt(dx * dx + dy * dy + half);
                float coef = (m[j] - q[i] * q[j] / m[i]) * invr * invr * invr;

                ax += coef * dx; /* accumulate the acceleration from gravitational attraction */
                ay += coef * dy;
            }

            xn[i] = x[i] + vx[i] * dt + half * ax * dt2; /* update position of particle "i" */
            yn[i] = y[i] + vy[i] * dt + half * ay * dt2;
            vx[i] += ax * dt; /* update velocity of particle "i" */
            vy[i] += ay * dt;
        }

#ifdef VISUAL
        for (unsigned int i = 0; i < n_particles; ++i) {
            if (xn[i] < 0) {
                xn[i] = -xn[i];
                vx[i] = 0.5f * std::fabs(vx[i]);
            } else if (xn[i] > img_width - 1) {
                xn[i] = 2 * (img_width - 1) - xn[i];
                vx[i] = -0.5f * std::fabs(vx[i]);
            }

            if (yn[i] < 0) {
                yn[i] = -yn[i];
                vy[i] = 0.5f * std::fabs(vy[i]);
            } else if (yn[i] > img_height - 1) {
                yn[i] = 2 * (img_height - 1) - yn[i];
                vy[i] = -0.5f * std::fabs(vy[i]);
            }
        }
#endif

        if( !(*cb)(cb_arg, step, xn, yn, vx, vy) )
            break;

        std::swap(x, xn);
        std::swap(y, yn);
    }

    return step;
}
