#include <cmath>
#include <cstring>

#include <string>
#include <iostream>

#include <dlfcn.h>
#include <immintrin.h>

#include "simulators.h"

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
simulator_unroll_2_0(
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
        for (unsigned int i = 0; i < n_particles - 1; i += 2) {
            float ax_0 = 0.0f;
            float ax_1 = 0.0f;

            float ay_0 = 0.0f;
            float ay_1 = 0.0f;

            for (unsigned int j = 0; j < n_particles; ++j) {
                float dx_0 = x[j] - x[i + 0];
                float dx_1 = x[j] - x[i + 1];

                float dy_0 = y[j] - y[i + 0];
                float dy_1 = y[j] - y[i + 1];

                float invr_0 = 1.0f / std::sqrt(dx_0 * dx_0 + dy_0 * dy_0 + half);
                float invr_1 = 1.0f / std::sqrt(dx_1 * dx_1 + dy_1 * dy_1 + half);

                float coef_0 = (m[j] - q[i + 0] * q[j] / m[i + 0]) * invr_0 * invr_0 * invr_0;
                float coef_1 = (m[j] - q[i + 1] * q[j] / m[i + 1]) * invr_1 * invr_1 * invr_1;

                ax_0 += coef_0 * dx_0; /* accumulate the acceleration from gravitational attraction */
                ax_1 += coef_1 * dx_1;

                ay_0 += coef_0 * dy_0;
                ay_1 += coef_1 * dy_1;
            }

            xn[i + 0] = x[i + 0] + vx[i + 0] * dt + half * ax_0 * dt2; /* update position of particle "i" */
            xn[i + 1] = x[i + 1] + vx[i + 1] * dt + half * ax_1 * dt2;

            yn[i + 0] = y[i + 0] + vy[i + 0] * dt + half * ay_0 * dt2;
            yn[i + 1] = y[i + 1] + vy[i + 1] * dt + half * ay_1 * dt2;

            vx[i + 0] += ax_0 * dt; /* update velocity of particle "i" */
            vx[i + 1] += ax_1 * dt;

            vy[i + 0] += ay_0 * dt;
            vy[i + 1] += ay_1 * dt;
        }

        if (n_particles & 0x1) {
            const unsigned int end = n_particles - 1;

            float ax = 0.0f;
            float ay = 0.0f;

            for (unsigned int j = 0; j < n_particles; ++j) {
                const float dx = x[j] - x[end];
                const float dy = y[j] - y[end];
                const float invr = 1.0f / std::sqrt(dx * dx + dy * dy + half);
                const float coef = (m[j] - q[end] * q[j] / m[end]) * invr * invr * invr;

                /* accumulate the acceleration from gravitational attraction */
                ax += coef * dx;
                ay += coef * dy;
            }

            /* update position of particle "end" */
            xn[end] = x[end] + vx[end] * dt + half * ax * dt2;
            yn[end] = y[end] + vy[end] * dt + half * ay * dt2;

            /* update velocity of particle "end" */
            vx[end] += ax * dt;
            vy[end] += ay * dt;
        }

        if( !(*cb)(cb_arg, step) )
            break;

        std::swap(x, xn);
        std::swap(y, yn);
    }

    return step;
}

extern "C" unsigned int
simulator_unroll_4_0(
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
        for (unsigned int i = 0; i < n_particles - 3; i += 4) {
            float ax_0 = 0.0f;
            float ax_1 = 0.0f;
            float ax_2 = 0.0f;
            float ax_3 = 0.0f;

            float ay_0 = 0.0f;
            float ay_1 = 0.0f;
            float ay_2 = 0.0f;
            float ay_3 = 0.0f;

            for (unsigned int j = 0; j < n_particles; ++j) {
                float dx_0 = x[j] - x[i + 0];
                float dx_1 = x[j] - x[i + 1];
                float dx_2 = x[j] - x[i + 2];
                float dx_3 = x[j] - x[i + 3];

                float dy_0 = y[j] - y[i + 0];
                float dy_1 = y[j] - y[i + 1];
                float dy_2 = y[j] - y[i + 2];
                float dy_3 = y[j] - y[i + 3];

                float invr_0 = 1.0f / std::sqrt(dx_0 * dx_0 + dy_0 * dy_0 + half);
                float invr_1 = 1.0f / std::sqrt(dx_1 * dx_1 + dy_1 * dy_1 + half);
                float invr_2 = 1.0f / std::sqrt(dx_2 * dx_2 + dy_2 * dy_2 + half);
                float invr_3 = 1.0f / std::sqrt(dx_3 * dx_3 + dy_3 * dy_3 + half);

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

            xn[i + 0] = x[i + 0] + vx[i + 0] * dt + half * ax_0 * dt2; /* update position of particle "i" */
            xn[i + 1] = x[i + 1] + vx[i + 1] * dt + half * ax_1 * dt2;
            xn[i + 2] = x[i + 2] + vx[i + 2] * dt + half * ax_2 * dt2;
            xn[i + 3] = x[i + 3] + vx[i + 3] * dt + half * ax_3 * dt2;

            yn[i + 0] = y[i + 0] + vy[i + 0] * dt + half * ay_0 * dt2;
            yn[i + 1] = y[i + 1] + vy[i + 1] * dt + half * ay_1 * dt2;
            yn[i + 2] = y[i + 2] + vy[i + 2] * dt + half * ay_2 * dt2;
            yn[i + 3] = y[i + 3] + vy[i + 3] * dt + half * ay_3 * dt2;

            vx[i + 0] += ax_0 * dt; /* update velocity of particle "i" */
            vx[i + 1] += ax_1 * dt;
            vx[i + 2] += ax_2 * dt;
            vx[i + 3] += ax_3 * dt;

            vy[i + 0] += ay_0 * dt;
            vy[i + 1] += ay_1 * dt;
            vy[i + 2] += ay_2 * dt;
            vy[i + 3] += ay_3 * dt;
        }

        for (unsigned int i = n_particles & ~0x3u; i < n_particles; ++i) {
            float ax = 0.0f;
            float ay = 0.0f;

            for (unsigned int j = 0; j < n_particles; ++j) {
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

        if( !(*cb)(cb_arg, step) )
            break;

        std::swap(x, xn);
        std::swap(y, yn);
    }

    return step;
}


extern "C" unsigned int
simulator_unroll_8_0(
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
        for (unsigned int i = 0; i < n_particles - 7; i += 8) {
            float ax_0 = 0.0f;
            float ax_1 = 0.0f;
            float ax_2 = 0.0f;
            float ax_3 = 0.0f;
            float ax_4 = 0.0f;
            float ax_5 = 0.0f;
            float ax_6 = 0.0f;
            float ax_7 = 0.0f;

            float ay_0 = 0.0f;
            float ay_1 = 0.0f;
            float ay_2 = 0.0f;
            float ay_3 = 0.0f;
            float ay_4 = 0.0f;
            float ay_5 = 0.0f;
            float ay_6 = 0.0f;
            float ay_7 = 0.0f;

            for (unsigned int j = 0; j < n_particles; ++j) {
                float dx_0 = x[j] - x[i + 0];
                float dx_1 = x[j] - x[i + 1];
                float dx_2 = x[j] - x[i + 2];
                float dx_3 = x[j] - x[i + 3];
                float dx_4 = x[j] - x[i + 4];
                float dx_5 = x[j] - x[i + 5];
                float dx_6 = x[j] - x[i + 6];
                float dx_7 = x[j] - x[i + 7];

                float dy_0 = y[j] - y[i + 0];
                float dy_1 = y[j] - y[i + 1];
                float dy_2 = y[j] - y[i + 2];
                float dy_3 = y[j] - y[i + 3];
                float dy_4 = y[j] - y[i + 4];
                float dy_5 = y[j] - y[i + 5];
                float dy_6 = y[j] - y[i + 6];
                float dy_7 = y[j] - y[i + 7];

                float invr_0 = 1.0f / std::sqrt(dx_0 * dx_0 + dy_0 * dy_0 + half);
                float invr_1 = 1.0f / std::sqrt(dx_1 * dx_1 + dy_1 * dy_1 + half);
                float invr_2 = 1.0f / std::sqrt(dx_2 * dx_2 + dy_2 * dy_2 + half);
                float invr_3 = 1.0f / std::sqrt(dx_3 * dx_3 + dy_3 * dy_3 + half);
                float invr_4 = 1.0f / std::sqrt(dx_4 * dx_4 + dy_4 * dy_4 + half);
                float invr_5 = 1.0f / std::sqrt(dx_5 * dx_5 + dy_5 * dy_5 + half);
                float invr_6 = 1.0f / std::sqrt(dx_6 * dx_6 + dy_6 * dy_6 + half);
                float invr_7 = 1.0f / std::sqrt(dx_7 * dx_7 + dy_7 * dy_7 + half);

                float coef_0 = (m[j] - q[i + 0] * q[j] / m[i + 0]) * invr_0 * invr_0 * invr_0;
                float coef_1 = (m[j] - q[i + 1] * q[j] / m[i + 1]) * invr_1 * invr_1 * invr_1;
                float coef_2 = (m[j] - q[i + 2] * q[j] / m[i + 2]) * invr_2 * invr_2 * invr_2;
                float coef_3 = (m[j] - q[i + 3] * q[j] / m[i + 3]) * invr_3 * invr_3 * invr_3;
                float coef_4 = (m[j] - q[i + 4] * q[j] / m[i + 4]) * invr_4 * invr_4 * invr_4;
                float coef_5 = (m[j] - q[i + 5] * q[j] / m[i + 5]) * invr_5 * invr_5 * invr_5;
                float coef_6 = (m[j] - q[i + 6] * q[j] / m[i + 6]) * invr_6 * invr_6 * invr_6;
                float coef_7 = (m[j] - q[i + 7] * q[j] / m[i + 7]) * invr_7 * invr_7 * invr_7;

                ax_0 += coef_0 * dx_0; /* accumulate the acceleration from gravitational attraction */
                ax_1 += coef_1 * dx_1;
                ax_2 += coef_2 * dx_2;
                ax_3 += coef_3 * dx_3;
                ax_4 += coef_4 * dx_4;
                ax_5 += coef_5 * dx_5;
                ax_6 += coef_6 * dx_6;
                ax_7 += coef_7 * dx_7;

                ay_0 += coef_0 * dy_0;
                ay_1 += coef_1 * dy_1;
                ay_2 += coef_2 * dy_2;
                ay_3 += coef_3 * dy_3;
                ay_4 += coef_4 * dy_4;
                ay_5 += coef_5 * dy_5;
                ay_6 += coef_6 * dy_6;
                ay_7 += coef_7 * dy_7;
            }

            xn[i + 0] = x[i + 0] + vx[i + 0] * dt + half * ax_0 * dt2; /* update position of particle "i" */
            xn[i + 1] = x[i + 1] + vx[i + 1] * dt + half * ax_1 * dt2;
            xn[i + 2] = x[i + 2] + vx[i + 2] * dt + half * ax_2 * dt2;
            xn[i + 3] = x[i + 3] + vx[i + 3] * dt + half * ax_3 * dt2;
            xn[i + 4] = x[i + 4] + vx[i + 4] * dt + half * ax_4 * dt2;
            xn[i + 5] = x[i + 5] + vx[i + 5] * dt + half * ax_5 * dt2;
            xn[i + 6] = x[i + 6] + vx[i + 6] * dt + half * ax_6 * dt2;
            xn[i + 7] = x[i + 7] + vx[i + 7] * dt + half * ax_7 * dt2;

            yn[i + 0] = y[i + 0] + vy[i + 0] * dt + half * ay_0 * dt2;
            yn[i + 1] = y[i + 1] + vy[i + 1] * dt + half * ay_1 * dt2;
            yn[i + 2] = y[i + 2] + vy[i + 2] * dt + half * ay_2 * dt2;
            yn[i + 3] = y[i + 3] + vy[i + 3] * dt + half * ay_3 * dt2;
            yn[i + 4] = y[i + 4] + vy[i + 4] * dt + half * ay_4 * dt2;
            yn[i + 5] = y[i + 5] + vy[i + 5] * dt + half * ay_5 * dt2;
            yn[i + 6] = y[i + 6] + vy[i + 6] * dt + half * ay_6 * dt2;
            yn[i + 7] = y[i + 7] + vy[i + 7] * dt + half * ay_7 * dt2;

            vx[i + 0] += ax_0 * dt; /* update velocity of particle "i" */
            vx[i + 1] += ax_1 * dt;
            vx[i + 2] += ax_2 * dt;
            vx[i + 3] += ax_3 * dt;
            vx[i + 4] += ax_4 * dt;
            vx[i + 5] += ax_5 * dt;
            vx[i + 6] += ax_6 * dt;
            vx[i + 7] += ax_7 * dt;

            vy[i + 0] += ay_0 * dt;
            vy[i + 1] += ay_1 * dt;
            vy[i + 2] += ay_2 * dt;
            vy[i + 3] += ay_3 * dt;
            vy[i + 4] += ay_4 * dt;
            vy[i + 5] += ay_5 * dt;
            vy[i + 6] += ay_6 * dt;
            vy[i + 7] += ay_7 * dt;
        }

        for (unsigned int i = n_particles & ~0x7u; i < n_particles; ++i) {
            float ax = 0.0f;
            float ay = 0.0f;

            for (unsigned int j = 0; j < n_particles; ++j) {
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

        if( !(*cb)(cb_arg, step) )
            break;

        std::swap(x, xn);
        std::swap(y, yn);
    }

    return step;
}

extern "C" unsigned int
simulator_unroll_0_2(
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

            for (unsigned int j = 0; j < n_particles - 1; j += 2) {
                float dx_0 = x[j + 0] - x[i];
                float dx_1 = x[j + 1] - x[i];

                float dy_0 = y[j + 0] - y[i];
                float dy_1 = y[j + 1] - y[i];

                float invr_0 = 1.0f / std::sqrt(dx_0 * dx_0 + dy_0 * dy_0 + half);
                float invr_1 = 1.0f / std::sqrt(dx_1 * dx_1 + dy_1 * dy_1 + half);

                float coef_0 = (m[j + 0] - q[i] * q[j + 0] / m[i]) * invr_0 * invr_0 * invr_0;
                float coef_1 = (m[j + 1] - q[i] * q[j + 1] / m[i]) * invr_1 * invr_1 * invr_1;

                ax += coef_0 * dx_0; /* accumulate the acceleration from gravitational attraction */
                ax += coef_1 * dx_1;

                ay += coef_0 * dy_0;
                ay += coef_1 * dy_1;
            }

            if (n_particles & 0x1) {
                const unsigned int end = n_particles - 1;

                float dx = x[end] - x[i];
                float dy = y[end] - y[i];

                float invr = 1.0f / std::sqrt(dx * dx + dy * dy + half);
                float coef = (m[end] - q[i] * q[end] / m[i]) * invr * invr * invr;

                ax += coef * dx; /* accumulate the acceleration from gravitational attraction */
                ay += coef * dy;
            }

            xn[i] = x[i] + vx[i] * dt + half * ax * dt2; /* update position of particle "i" */
            yn[i] = y[i] + vy[i] * dt + half * ay * dt2;

            vx[i] += ax * dt; /* update velocity of particle "i" */
            vy[i] += ay * dt;
        }

        if( !(*cb)(cb_arg, step) )
            break;

        std::swap(x, xn);
        std::swap(y, yn);
    }

    return step;
}

extern "C" unsigned int
simulator_unroll_0_4(
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

            for (unsigned int j = 0; j < n_particles - 3; j += 4) {
                float dx_0 = x[j + 0] - x[i];
                float dx_1 = x[j + 1] - x[i];
                float dx_2 = x[j + 2] - x[i];
                float dx_3 = x[j + 3] - x[i];

                float dy_0 = y[j + 0] - y[i];
                float dy_1 = y[j + 1] - y[i];
                float dy_2 = y[j + 2] - y[i];
                float dy_3 = y[j + 3] - y[i];

                float invr_0 = 1.0f / std::sqrt(dx_0 * dx_0 + dy_0 * dy_0 + half);
                float invr_1 = 1.0f / std::sqrt(dx_1 * dx_1 + dy_1 * dy_1 + half);
                float invr_2 = 1.0f / std::sqrt(dx_2 * dx_2 + dy_2 * dy_2 + half);
                float invr_3 = 1.0f / std::sqrt(dx_3 * dx_3 + dy_3 * dy_3 + half);

                float coef_0 = (m[j + 0] - q[i] * q[j + 0] / m[i]) * invr_0 * invr_0 * invr_0;
                float coef_1 = (m[j + 1] - q[i] * q[j + 1] / m[i]) * invr_1 * invr_1 * invr_1;
                float coef_2 = (m[j + 2] - q[i] * q[j + 2] / m[i]) * invr_2 * invr_2 * invr_2;
                float coef_3 = (m[j + 3] - q[i] * q[j + 3] / m[i]) * invr_3 * invr_3 * invr_3;

                ax += coef_0 * dx_0; /* accumulate the acceleration from gravitational attraction */
                ax += coef_1 * dx_1;
                ax += coef_2 * dx_2;
                ax += coef_3 * dx_3;

                ay += coef_0 * dy_0;
                ay += coef_1 * dy_1;
                ay += coef_2 * dy_2;
                ay += coef_3 * dy_3;
            }

            for (unsigned int j = n_particles & ~0x3u; j < n_particles; ++j) {
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

        if( !(*cb)(cb_arg, step) )
            break;

        std::swap(x, xn);
        std::swap(y, yn);
    }

    return step;
}


extern "C" unsigned int
simulator_unroll_0_8(
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

            for (unsigned int j = 0; j < n_particles - 7; j += 8) {
                float dx_0 = x[j + 0] - x[i];
                float dx_1 = x[j + 1] - x[i];
                float dx_2 = x[j + 2] - x[i];
                float dx_3 = x[j + 3] - x[i];
                float dx_4 = x[j + 4] - x[i];
                float dx_5 = x[j + 5] - x[i];
                float dx_6 = x[j + 6] - x[i];
                float dx_7 = x[j + 7] - x[i];

                float dy_0 = y[j + 0] - y[i];
                float dy_1 = y[j + 1] - y[i];
                float dy_2 = y[j + 2] - y[i];
                float dy_3 = y[j + 3] - y[i];
                float dy_4 = y[j + 4] - y[i];
                float dy_5 = y[j + 5] - y[i];
                float dy_6 = y[j + 6] - y[i];
                float dy_7 = y[j + 7] - y[i];

                float invr_0 = 1.0f / std::sqrt(dx_0 * dx_0 + dy_0 * dy_0 + half);
                float invr_1 = 1.0f / std::sqrt(dx_1 * dx_1 + dy_1 * dy_1 + half);
                float invr_2 = 1.0f / std::sqrt(dx_2 * dx_2 + dy_2 * dy_2 + half);
                float invr_3 = 1.0f / std::sqrt(dx_3 * dx_3 + dy_3 * dy_3 + half);
                float invr_4 = 1.0f / std::sqrt(dx_4 * dx_4 + dy_4 * dy_4 + half);
                float invr_5 = 1.0f / std::sqrt(dx_5 * dx_5 + dy_5 * dy_5 + half);
                float invr_6 = 1.0f / std::sqrt(dx_6 * dx_6 + dy_6 * dy_6 + half);
                float invr_7 = 1.0f / std::sqrt(dx_7 * dx_7 + dy_7 * dy_7 + half);

                float coef_0 = (m[j + 0] - q[i] * q[j + 0] / m[i]) * invr_0 * invr_0 * invr_0;
                float coef_1 = (m[j + 1] - q[i] * q[j + 1] / m[i]) * invr_1 * invr_1 * invr_1;
                float coef_2 = (m[j + 2] - q[i] * q[j + 2] / m[i]) * invr_2 * invr_2 * invr_2;
                float coef_3 = (m[j + 3] - q[i] * q[j + 3] / m[i]) * invr_3 * invr_3 * invr_3;
                float coef_4 = (m[j + 4] - q[i] * q[j + 4] / m[i]) * invr_4 * invr_4 * invr_4;
                float coef_5 = (m[j + 5] - q[i] * q[j + 5] / m[i]) * invr_5 * invr_5 * invr_5;
                float coef_6 = (m[j + 6] - q[i] * q[j + 6] / m[i]) * invr_6 * invr_6 * invr_6;
                float coef_7 = (m[j + 7] - q[i] * q[j + 7] / m[i]) * invr_7 * invr_7 * invr_7;

                ax += coef_0 * dx_0; /* accumulate the acceleration from gravitational attraction */
                ax += coef_1 * dx_1;
                ax += coef_2 * dx_2;
                ax += coef_3 * dx_3;
                ax += coef_4 * dx_4;
                ax += coef_5 * dx_5;
                ax += coef_6 * dx_6;
                ax += coef_7 * dx_7;

                ay += coef_0 * dy_0;
                ay += coef_1 * dy_1;
                ay += coef_2 * dy_2;
                ay += coef_3 * dy_3;
                ay += coef_4 * dy_4;
                ay += coef_5 * dy_5;
                ay += coef_6 * dy_6;
                ay += coef_7 * dy_7;
            }

            for (unsigned int j = n_particles & ~0x7u; j < n_particles; ++j) {
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

        if( !(*cb)(cb_arg, step) )
            break;

        std::swap(x, xn);
        std::swap(y, yn);
    }

    return step;
}

extern "C" unsigned int
simulator_unroll_2_2(
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
        for (unsigned int i = 0; i < n_particles - 1; i += 2) {
            float ax_0 = 0.0f;
            float ax_1 = 0.0f;

            float ay_0 = 0.0f;
            float ay_1 = 0.0f;

            for (unsigned int j = 0; j < n_particles - 1; j += 2) {
                float dx_0_0 = x[j + 0] - x[i + 0];
                float dx_0_1 = x[j + 1] - x[i + 0];
                float dy_0_0 = y[j + 0] - y[i + 0];
                float dy_0_1 = y[j + 1] - y[i + 0];

                float dx_1_0 = x[j + 0] - x[i + 1];
                float dx_1_1 = x[j + 1] - x[i + 1];
                float dy_1_0 = y[j + 0] - y[i + 1];
                float dy_1_1 = y[j + 1] - y[i + 1];

                float invr_0_0 = 1.0f / std::sqrt(dx_0_0 * dx_0_0 + dy_0_0 * dy_0_0 + half);
                float invr_0_1 = 1.0f / std::sqrt(dx_0_1 * dx_0_1 + dy_0_1 * dy_0_1 + half);
                float invr_1_0 = 1.0f / std::sqrt(dx_1_0 * dx_1_0 + dy_1_0 * dy_1_0 + half);
                float invr_1_1 = 1.0f / std::sqrt(dx_1_1 * dx_1_1 + dy_1_1 * dy_1_1 + half);

                float coef_0_0 = (m[j + 0] - q[i + 0] * q[j + 0] / m[i + 0]) * invr_0_0 * invr_0_0 * invr_0_0;
                float coef_0_1 = (m[j + 1] - q[i + 0] * q[j + 1] / m[i + 0]) * invr_0_1 * invr_0_1 * invr_0_1;
                float coef_1_0 = (m[j + 0] - q[i + 1] * q[j + 0] / m[i + 1]) * invr_1_0 * invr_1_0 * invr_1_0;
                float coef_1_1 = (m[j + 1] - q[i + 1] * q[j + 1] / m[i + 1]) * invr_1_1 * invr_1_1 * invr_1_1;

                ax_0 += coef_0_0 * dx_0_0; /* accumulate the acceleration from gravitational attraction */
                ax_0 += coef_0_1 * dx_0_1;
                ax_1 += coef_1_0 * dx_1_0;
                ax_1 += coef_1_1 * dx_1_1;

                ay_0 += coef_0_0 * dy_0_0; /* accumulate the acceleration from gravitational attraction */
                ay_0 += coef_0_1 * dy_0_1;
                ay_1 += coef_1_0 * dy_1_0;
                ay_1 += coef_1_1 * dy_1_1;
            }

            if (n_particles & 0x1) {
                const unsigned int end = n_particles - 1;

                float dx_0 = x[end] - x[i + 0];
                float dx_1 = x[end] - x[i + 1];

                float dy_0 = y[end] - y[i + 0];
                float dy_1 = y[end] - y[i + 1];

                float invr_0 = 1.0f / std::sqrt(dx_0 * dx_0 + dy_0 * dy_0 + half);
                float invr_1 = 1.0f / std::sqrt(dx_1 * dx_1 + dy_1 * dy_1 + half);

                float coef_0 = (m[end] - q[i + 0] * q[end] / m[i + 0]) * invr_0 * invr_0 * invr_0;
                float coef_1 = (m[end] - q[i + 1] * q[end] / m[i + 1]) * invr_1 * invr_1 * invr_1;

                ax_0 += coef_0 * dx_0; /* accumulate the acceleration from gravitational attraction */
                ax_1 += coef_1 * dx_1;

                ay_0 += coef_0 * dy_0;
                ay_1 += coef_1 * dy_1;
            }

            xn[i + 0] = x[i + 0] + vx[i + 0] * dt + half * ax_0 * dt2; /* update position of particle "i" */
            xn[i + 1] = x[i + 1] + vx[i + 1] * dt + half * ax_1 * dt2;

            yn[i + 0] = y[i + 0] + vy[i + 0] * dt + half * ay_0 * dt2;
            yn[i + 1] = y[i + 1] + vy[i + 1] * dt + half * ay_1 * dt2;

            vx[i + 0] += ax_0 * dt; /* update velocity of particle "i" */
            vx[i + 1] += ax_1 * dt;

            vy[i + 0] += ay_0 * dt;
            vy[i + 1] += ay_1 * dt;
        }

        if (n_particles & 0x1) {
            float ax = 0.0f;
            float ay = 0.0f;
            const unsigned int end = n_particles - 1;

            for (unsigned int j = 0; j < n_particles; ++j) {
                const float dx = x[j] - x[end];
                const float dy = y[j] - y[end];
                const float invr = 1.0f / std::sqrt(dx * dx + dy * dy + half);
                const float coef = (m[j] - q[end] * q[j] / m[end]) * invr * invr * invr;

                /* accumulate the acceleration from gravitational attraction */
                ax += coef * dx;
                ay += coef * dy;
            }

            /* update position of particle "end" */
            xn[end] = x[end] + vx[end] * dt + half * ax * dt2;
            yn[end] = y[end] + vy[end] * dt + half * ay * dt2;

            /* update velocity of particle "end" */
            vx[end] += ax * dt;
            vy[end] += ay * dt;
        }

        if( !(*cb)(cb_arg, step) )
            break;

        std::swap(x, xn);
        std::swap(y, yn);
    }

    return step;
}

extern "C" unsigned int
simulator_vec_4_0_4(
    const simulator_conf_t* conf,
    const simulator_data_t* data
    )
{
    LOAD_CONF( conf );
    LOAD_DATA( data );

    const float half = 0.5f;
    const float dt2 = dt * dt;

#   ifdef ACCURATE_VEC
    const __v4sf vone  = { 1.0f, 1.0f, 1.0f, 1.0f };
#   endif

    const __v4sf vhalf = { half, half, half, half };
    const __v4sf vdt = { dt, dt, dt, dt };
    const __v4sf vdt2 = vdt * vdt;

    unsigned int step;
    for (step = 0; step < n_steps; ++step) {
        for (unsigned int i = 0; i < n_particles - 3; i += 4) {
            __v4sf ax = { 0.0f, 0.0f, 0.0f, 0.0f };
            __v4sf ay = ax;

            const __v4sf xi = *(__v4sf *)(x + i);
            const __v4sf yi = *(__v4sf *)(y + i);
            const __v4sf mi = *(__v4sf *)(m + i);
            const __v4sf qi = *(__v4sf *)(q + i);

            for (unsigned int j = 0; j < n_particles; ++j) {
                const __v4sf xj = { x[j], x[j], x[j], x[j] };
                const __v4sf yj = { y[j], y[j], y[j], y[j] };
                const __v4sf mj = { m[j], m[j], m[j], m[j] };
                const __v4sf qj = { q[j], q[j], q[j], q[j] };

                const __v4sf dx = xj - xi;
                const __v4sf dy = yj - yi;

#               ifdef ACCURATE_VEC
                const __v4sf invr = vone / _mm_sqrt_ps( dx * dx + dy * dy + vhalf );
#               else
                const __v4sf invr = _mm_rsqrt_ps( dx * dx + dy * dy + vhalf );
#               endif
                const __v4sf coef = (mj - qi * qj / mi) * invr * invr * invr;

                /* accumulate the acceleration from gravitational attraction */
                ax += coef * dx;
                ay += coef * dy;
            }

            const __v4sf vxi = *(__v4sf *)(vx + i);
            const __v4sf vyi = *(__v4sf *)(vy + i);

            /* update position of particle "i" */
            *(__v4sf *)(xn + i) = xi + vxi * vdt + vhalf * ax * vdt2;
            *(__v4sf *)(yn + i) = yi + vyi * vdt + vhalf * ay * vdt2;

            /* update velocity of particle "i" */
            *(__v4sf *)(vx + i) = vxi + ax * vdt;
            *(__v4sf *)(vy + i) = vyi + ay * vdt;
        }

        for (unsigned int i = n_particles & ~0x3u; i < n_particles; ++i) {
            float ax = 0.0f;
            float ay = 0.0f;

            for (unsigned int j = 0; j < n_particles; ++j) {
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

        if( !(*cb)(cb_arg, step) )
            break;

        std::swap(x, xn);
        std::swap(y, yn);
    }

    return step;
}

extern "C" unsigned int
simulator_vec_8_0_4(
    const simulator_conf_t* conf,
    const simulator_data_t* data
    )
{
    LOAD_CONF( conf );
    LOAD_DATA( data );

    const float half = 0.5f;
    const float dt2 = dt * dt;

#   ifdef ACCURATE_VEC
    const __v4sf vone  = { 1.0f, 1.0f, 1.0f, 1.0f };
#   endif

    const __v4sf vhalf = { half, half, half, half };
    const __v4sf vdt = { dt, dt, dt, dt };
    const __v4sf vdt2 = vdt * vdt;

    unsigned int step;
    for (step = 0; step < n_steps; ++step) {
        for (unsigned int i = 0; i < n_particles - 7; i += 8) {
            __v4sf ax_0 = { 0.0f, 0.0f, 0.0f, 0.0f };
            __v4sf ay_0 = ax_0;
            __v4sf ax_4 = ax_0;
            __v4sf ay_4 = ax_0;

            const __v4sf xi_0 = *(__v4sf *)(x + i + 0);
            const __v4sf yi_0 = *(__v4sf *)(y + i + 0);
            const __v4sf mi_0 = *(__v4sf *)(m + i + 0);
            const __v4sf qi_0 = *(__v4sf *)(q + i + 0);

            const __v4sf xi_4 = *(__v4sf *)(x + i + 4);
            const __v4sf yi_4 = *(__v4sf *)(y + i + 4);
            const __v4sf mi_4 = *(__v4sf *)(m + i + 4);
            const __v4sf qi_4 = *(__v4sf *)(q + i + 4);

            for (unsigned int j = 0; j < n_particles; ++j) {
                const __v4sf xj = { x[j], x[j], x[j], x[j] };
                const __v4sf yj = { y[j], y[j], y[j], y[j] };
                const __v4sf mj = { m[j], m[j], m[j], m[j] };
                const __v4sf qj = { q[j], q[j], q[j], q[j] };

                const __v4sf dx_0 = xj - xi_0;
                const __v4sf dy_0 = yj - yi_0;
                const __v4sf dx_4 = xj - xi_4;
                const __v4sf dy_4 = yj - yi_4;

#               ifdef ACCURATE_VEC
                const __v4sf invr_0 = vone / _mm_sqrt_ps( dx_0 * dx_0 + dy_0 * dy_0 + vhalf );
                const __v4sf invr_4 = vone / _mm_sqrt_ps( dx_4 * dx_4 + dy_4 * dy_4 + vhalf );
#               else
                const __v4sf invr_0 = _mm_rsqrt_ps( dx_0 * dx_0 + dy_0 * dy_0 + vhalf );
                const __v4sf invr_4 = _mm_rsqrt_ps( dx_4 * dx_4 + dy_4 * dy_4 + vhalf );
#               endif
                const __v4sf coef_0 = (mj - qi_0 * qj / mi_0) * invr_0 * invr_0 * invr_0;
                const __v4sf coef_4 = (mj - qi_4 * qj / mi_4) * invr_4 * invr_4 * invr_4;

                /* accumulate the acceleration from gravitational attraction */
                ax_0 += coef_0 * dx_0;
                ay_0 += coef_0 * dy_0;
                ax_4 += coef_4 * dx_4;
                ay_4 += coef_4 * dy_4;
            }

            const __v4sf vxi_0 = *(__v4sf *)(vx + i + 0);
            const __v4sf vyi_0 = *(__v4sf *)(vy + i + 0);
            const __v4sf vxi_4 = *(__v4sf *)(vx + i + 4);
            const __v4sf vyi_4 = *(__v4sf *)(vy + i + 4);

            /* update position of particle "i" */
            *(__v4sf *)(xn + i + 0) = xi_0 + vxi_0 * vdt + vhalf * ax_0 * vdt2;
            *(__v4sf *)(yn + i + 0) = yi_0 + vyi_0 * vdt + vhalf * ay_0 * vdt2;
            *(__v4sf *)(xn + i + 4) = xi_4 + vxi_4 * vdt + vhalf * ax_4 * vdt2;
            *(__v4sf *)(yn + i + 4) = yi_4 + vyi_4 * vdt + vhalf * ay_4 * vdt2;

            /* update velocity of particle "i" */
            *(__v4sf *)(vx + i + 0) = vxi_0 + ax_0 * vdt;
            *(__v4sf *)(vy + i + 0) = vyi_0 + ay_0 * vdt;
            *(__v4sf *)(vx + i + 4) = vxi_4 + ax_4 * vdt;
            *(__v4sf *)(vy + i + 4) = vyi_4 + ay_4 * vdt;
        }

        for (unsigned int i = n_particles & ~0x7u; i < n_particles; ++i) {
            float ax = 0.0f;
            float ay = 0.0f;

            for (unsigned int j = 0; j < n_particles; ++j) {
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

        if( !(*cb)(cb_arg, step) )
            break;

        std::swap(x, xn);
        std::swap(y, yn);
    }

    return step;
}

#ifdef __AVX__
extern "C" unsigned int
simulator_vec_8_0_8(
    const simulator_conf_t* conf,
    const simulator_data_t* data
    )
{
    LOAD_CONF( conf );
    LOAD_DATA( data );

    const float half = 0.5f;
    const float dt2 = dt * dt;

#   ifdef ACCURATE_VEC
    const __v8sf vone  = { 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f };
#   endif
    const __v8sf vhalf = { half, half, half, half, half, half, half, half };
    const __v8sf vdt = { dt, dt, dt, dt, dt, dt, dt, dt };
    const __v8sf vdt2 = vdt * vdt;

    unsigned int step;
    for (step = 0; step < n_steps; ++step) {
        for (unsigned int i = 0; i < n_particles - 7; i += 8) {
            __v8sf ax = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
            __v8sf ay = ax;

            const __v8sf xi = *(__v8sf *)(x + i);
            const __v8sf yi = *(__v8sf *)(y + i);
            const __v8sf mi = *(__v8sf *)(m + i);
            const __v8sf qi = *(__v8sf *)(q + i);

            for (unsigned int j = 0; j < n_particles; ++j) {
                const __v8sf xj = { x[j], x[j], x[j], x[j], x[j], x[j], x[j], x[j] };
                const __v8sf yj = { y[j], y[j], y[j], y[j], y[j], y[j], y[j], y[j] };
                const __v8sf mj = { m[j], m[j], m[j], m[j], m[j], m[j], m[j], m[j] };
                const __v8sf qj = { q[j], q[j], q[j], q[j], q[j], q[j], q[j], q[j] };

                const __v8sf dx = xj - xi;
                const __v8sf dy = yj - yi;

#               ifdef ACCURATE_VEC
                const __v8sf invr = vone / _mm256_sqrt_ps( dx * dx + dy * dy + vhalf );
#               else
                const __v8sf invr = _mm256_rsqrt_ps( dx * dx + dy * dy + vhalf );
#               endif
                const __v8sf coef = (mj - qi * qj / mi) * invr * invr * invr;

                /* accumulate the acceleration from gravitational attraction */
                ax += coef * dx;
                ay += coef * dy;
            }

            const __v8sf vxi = *(__v8sf *)(vx + i);
            const __v8sf vyi = *(__v8sf *)(vy + i);

            /* update position of particle "i" */
            *(__v8sf *)(xn + i) = xi + vxi * vdt + vhalf * ax * vdt2;
            *(__v8sf *)(yn + i) = yi + vyi * vdt + vhalf * ay * vdt2;

            /* update velocity of particle "i" */
            *(__v8sf *)(vx + i) = vxi + ax * vdt;
            *(__v8sf *)(vy + i) = vyi + ay * vdt;
        }

        for (unsigned int i = n_particles & ~0x7u; i < n_particles; ++i) {
            float ax = 0.0f;
            float ay = 0.0f;

            for (unsigned int j = 0; j < n_particles; ++j) {
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

        if( !(*cb)(cb_arg, step) )
            break;

        std::swap(x, xn);
        std::swap(y, yn);
    }

    return step;
}
#endif

extern "C" unsigned int
simulator_vec_12_0_4(
    const simulator_conf_t* conf,
    const simulator_data_t* data
    )
{
    LOAD_CONF( conf );
    LOAD_DATA( data );

    const float half = 0.5f;
    const float dt2 = dt * dt;

#   ifdef ACCURATE_VEC
    const __v4sf vone  = { 1.0f, 1.0f, 1.0f, 1.0f };
#   endif

    const __v4sf vhalf = { half, half, half, half };
    const __v4sf vdt = { dt, dt, dt, dt };
    const __v4sf vdt2 = vdt * vdt;

    unsigned int step;
    for (step = 0; step < n_steps; ++step) {
        unsigned int i;
        for (i = 0; i < n_particles - 11; i += 12) {
            __v4sf ax_0 = { 0.0f, 0.0f, 0.0f, 0.0f };
            __v4sf ay_0 = ax_0;
            __v4sf ax_4 = ax_0;
            __v4sf ay_4 = ax_0;
            __v4sf ax_8 = ax_0;
            __v4sf ay_8 = ax_0;

            const __v4sf xi_0 = *(__v4sf *)(x + i + 0);
            const __v4sf yi_0 = *(__v4sf *)(y + i + 0);
            const __v4sf mi_0 = *(__v4sf *)(m + i + 0);
            const __v4sf qi_0 = *(__v4sf *)(q + i + 0);

            const __v4sf xi_4 = *(__v4sf *)(x + i + 4);
            const __v4sf yi_4 = *(__v4sf *)(y + i + 4);
            const __v4sf mi_4 = *(__v4sf *)(m + i + 4);
            const __v4sf qi_4 = *(__v4sf *)(q + i + 4);

            const __v4sf xi_8 = *(__v4sf *)(x + i + 8);
            const __v4sf yi_8 = *(__v4sf *)(y + i + 8);
            const __v4sf mi_8 = *(__v4sf *)(m + i + 8);
            const __v4sf qi_8 = *(__v4sf *)(q + i + 8);

            for (unsigned int j = 0; j < n_particles; ++j) {
                const __v4sf xj = { x[j], x[j], x[j], x[j] };
                const __v4sf yj = { y[j], y[j], y[j], y[j] };
                const __v4sf mj = { m[j], m[j], m[j], m[j] };
                const __v4sf qj = { q[j], q[j], q[j], q[j] };

                const __v4sf dx_0 = xj - xi_0;
                const __v4sf dy_0 = yj - yi_0;
                const __v4sf dx_4 = xj - xi_4;
                const __v4sf dy_4 = yj - yi_4;
                const __v4sf dx_8 = xj - xi_8;
                const __v4sf dy_8 = yj - yi_8;

#               ifdef ACCURATE_VEC
                const __v4sf invr_0 = vone / _mm_sqrt_ps( dx_0 * dx_0 + dy_0 * dy_0 + vhalf );
                const __v4sf invr_4 = vone / _mm_sqrt_ps( dx_4 * dx_4 + dy_4 * dy_4 + vhalf );
                const __v4sf invr_8 = vone / _mm_sqrt_ps( dx_8 * dx_8 + dy_8 * dy_8 + vhalf );
#               else
                const __v4sf invr_0 = _mm_rsqrt_ps( dx_0 * dx_0 + dy_0 * dy_0 + vhalf );
                const __v4sf invr_4 = _mm_rsqrt_ps( dx_4 * dx_4 + dy_4 * dy_4 + vhalf );
                const __v4sf invr_8 = _mm_rsqrt_ps( dx_8 * dx_8 + dy_8 * dy_8 + vhalf );
#               endif
                const __v4sf coef_0 = (mj - qi_0 * qj / mi_0) * invr_0 * invr_0 * invr_0;
                const __v4sf coef_4 = (mj - qi_4 * qj / mi_4) * invr_4 * invr_4 * invr_4;
                const __v4sf coef_8 = (mj - qi_8 * qj / mi_8) * invr_8 * invr_8 * invr_8;

                /* accumulate the acceleration from gravitational attraction */
                ax_0 += coef_0 * dx_0;
                ay_0 += coef_0 * dy_0;
                ax_4 += coef_4 * dx_4;
                ay_4 += coef_4 * dy_4;
                ax_8 += coef_8 * dx_8;
                ay_8 += coef_8 * dy_8;
            }

            const __v4sf vxi_0 = *(__v4sf *)(vx + i + 0);
            const __v4sf vyi_0 = *(__v4sf *)(vy + i + 0);
            const __v4sf vxi_4 = *(__v4sf *)(vx + i + 4);
            const __v4sf vyi_4 = *(__v4sf *)(vy + i + 4);
            const __v4sf vxi_8 = *(__v4sf *)(vx + i + 8);
            const __v4sf vyi_8 = *(__v4sf *)(vy + i + 8);

            /* update position of particle "i" */
            *(__v4sf *)(xn + i + 0) = xi_0 + vxi_0 * vdt + vhalf * ax_0 * vdt2;
            *(__v4sf *)(yn + i + 0) = yi_0 + vyi_0 * vdt + vhalf * ay_0 * vdt2;
            *(__v4sf *)(xn + i + 4) = xi_4 + vxi_4 * vdt + vhalf * ax_4 * vdt2;
            *(__v4sf *)(yn + i + 4) = yi_4 + vyi_4 * vdt + vhalf * ay_4 * vdt2;
            *(__v4sf *)(xn + i + 8) = xi_8 + vxi_8 * vdt + vhalf * ax_8 * vdt2;
            *(__v4sf *)(yn + i + 8) = yi_8 + vyi_8 * vdt + vhalf * ay_8 * vdt2;

            /* update velocity of particle "i" */
            *(__v4sf *)(vx + i + 0) = vxi_0 + ax_0 * vdt;
            *(__v4sf *)(vy + i + 0) = vyi_0 + ay_0 * vdt;
            *(__v4sf *)(vx + i + 4) = vxi_4 + ax_4 * vdt;
            *(__v4sf *)(vy + i + 4) = vyi_4 + ay_4 * vdt;
            *(__v4sf *)(vx + i + 8) = vxi_8 + ax_8 * vdt;
            *(__v4sf *)(vy + i + 8) = vyi_8 + ay_8 * vdt;
        }

        for ( ; i < n_particles; ++i) {
            float ax = 0.0f;
            float ay = 0.0f;

            for (unsigned int j = 0; j < n_particles; ++j) {
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

        if( !(*cb)(cb_arg, step) )
            break;

        std::swap(x, xn);
        std::swap(y, yn);
    }

    return step;
}

extern "C" unsigned int
simulator_vec_16_0_4(
    const simulator_conf_t* conf,
    const simulator_data_t* data
    )
{
    LOAD_CONF( conf );
    LOAD_DATA( data );

    const float half = 0.5f;
    const float dt2 = dt * dt;

#   ifdef ACCURATE_VEC
    const __v4sf vone  = { 1.0f, 1.0f, 1.0f, 1.0f };
#   endif

    const __v4sf vhalf = { half, half, half, half };
    const __v4sf vdt = { dt, dt, dt, dt };
    const __v4sf vdt2 = vdt * vdt;

    unsigned int step;
    for (step = 0; step < n_steps; ++step) {
        for (unsigned int i = 0; i < n_particles - 15; i += 16) {
            __v4sf ax_0 = { 0.0f, 0.0f, 0.0f, 0.0f };
            __v4sf ay_0 = ax_0;
            __v4sf ax_4 = ax_0;
            __v4sf ay_4 = ax_0;
            __v4sf ax_8 = ax_0;
            __v4sf ay_8 = ax_0;
            __v4sf ax_12 = ax_0;
            __v4sf ay_12 = ax_0;

            const __v4sf xi_0 = *(__v4sf *)(x + i + 0);
            const __v4sf yi_0 = *(__v4sf *)(y + i + 0);
            const __v4sf mi_0 = *(__v4sf *)(m + i + 0);
            const __v4sf qi_0 = *(__v4sf *)(q + i + 0);

            const __v4sf xi_4 = *(__v4sf *)(x + i + 4);
            const __v4sf yi_4 = *(__v4sf *)(y + i + 4);
            const __v4sf mi_4 = *(__v4sf *)(m + i + 4);
            const __v4sf qi_4 = *(__v4sf *)(q + i + 4);

            const __v4sf xi_8 = *(__v4sf *)(x + i + 8);
            const __v4sf yi_8 = *(__v4sf *)(y + i + 8);
            const __v4sf mi_8 = *(__v4sf *)(m + i + 8);
            const __v4sf qi_8 = *(__v4sf *)(q + i + 8);

            const __v4sf xi_12 = *(__v4sf *)(x + i + 12);
            const __v4sf yi_12 = *(__v4sf *)(y + i + 12);
            const __v4sf mi_12 = *(__v4sf *)(m + i + 12);
            const __v4sf qi_12 = *(__v4sf *)(q + i + 12);

            for (unsigned int j = 0; j < n_particles; ++j) {
                const __v4sf xj = { x[j], x[j], x[j], x[j] };
                const __v4sf yj = { y[j], y[j], y[j], y[j] };
                const __v4sf mj = { m[j], m[j], m[j], m[j] };
                const __v4sf qj = { q[j], q[j], q[j], q[j] };

                const __v4sf dx_0 = xj - xi_0;
                const __v4sf dy_0 = yj - yi_0;
                const __v4sf dx_4 = xj - xi_4;
                const __v4sf dy_4 = yj - yi_4;
                const __v4sf dx_8 = xj - xi_8;
                const __v4sf dy_8 = yj - yi_8;
                const __v4sf dx_12 = xj - xi_12;
                const __v4sf dy_12 = yj - yi_12;

#               ifdef ACCURATE_VEC
                const __v4sf invr_0 = vone / _mm_sqrt_ps( dx_0 * dx_0 + dy_0 * dy_0 + vhalf );
                const __v4sf invr_4 = vone / _mm_sqrt_ps( dx_4 * dx_4 + dy_4 * dy_4 + vhalf );
                const __v4sf invr_8 = vone / _mm_sqrt_ps( dx_8 * dx_8 + dy_8 * dy_8 + vhalf );
                const __v4sf invr_12 = vone / _mm_sqrt_ps( dx_12 * dx_12 + dy_12 * dy_12 + vhalf );
#               else
                const __v4sf invr_0 = _mm_rsqrt_ps( dx_0 * dx_0 + dy_0 * dy_0 + vhalf );
                const __v4sf invr_4 = _mm_rsqrt_ps( dx_4 * dx_4 + dy_4 * dy_4 + vhalf );
                const __v4sf invr_8 = _mm_rsqrt_ps( dx_8 * dx_8 + dy_8 * dy_8 + vhalf );
                const __v4sf invr_12 = _mm_rsqrt_ps( dx_12 * dx_12 + dy_12 * dy_12 + vhalf );
#               endif
                const __v4sf coef_0 = (mj - qi_0 * qj / mi_0) * invr_0 * invr_0 * invr_0;
                const __v4sf coef_4 = (mj - qi_4 * qj / mi_4) * invr_4 * invr_4 * invr_4;
                const __v4sf coef_8 = (mj - qi_8 * qj / mi_8) * invr_8 * invr_8 * invr_8;
                const __v4sf coef_12 = (mj - qi_12 * qj / mi_12) * invr_12 * invr_12 * invr_12;

                /* accumulate the acceleration from gravitational attraction */
                ax_0 += coef_0 * dx_0;
                ay_0 += coef_0 * dy_0;
                ax_4 += coef_4 * dx_4;
                ay_4 += coef_4 * dy_4;
                ax_8 += coef_8 * dx_8;
                ay_8 += coef_8 * dy_8;
                ax_12 += coef_12 * dx_12;
                ay_12 += coef_12 * dy_12;
            }

            const __v4sf vxi_0 = *(__v4sf *)(vx + i + 0);
            const __v4sf vyi_0 = *(__v4sf *)(vy + i + 0);
            const __v4sf vxi_4 = *(__v4sf *)(vx + i + 4);
            const __v4sf vyi_4 = *(__v4sf *)(vy + i + 4);
            const __v4sf vxi_8 = *(__v4sf *)(vx + i + 8);
            const __v4sf vyi_8 = *(__v4sf *)(vy + i + 8);
            const __v4sf vxi_12 = *(__v4sf *)(vx + i + 12);
            const __v4sf vyi_12 = *(__v4sf *)(vy + i + 12);

            /* update position of particle "i" */
            *(__v4sf *)(xn + i + 0) = xi_0 + vxi_0 * vdt + vhalf * ax_0 * vdt2;
            *(__v4sf *)(yn + i + 0) = yi_0 + vyi_0 * vdt + vhalf * ay_0 * vdt2;
            *(__v4sf *)(xn + i + 4) = xi_4 + vxi_4 * vdt + vhalf * ax_4 * vdt2;
            *(__v4sf *)(yn + i + 4) = yi_4 + vyi_4 * vdt + vhalf * ay_4 * vdt2;
            *(__v4sf *)(xn + i + 8) = xi_8 + vxi_8 * vdt + vhalf * ax_8 * vdt2;
            *(__v4sf *)(yn + i + 8) = yi_8 + vyi_8 * vdt + vhalf * ay_8 * vdt2;
            *(__v4sf *)(xn + i + 12) = xi_12 + vxi_12 * vdt + vhalf * ax_12 * vdt2;
            *(__v4sf *)(yn + i + 12) = yi_12 + vyi_12 * vdt + vhalf * ay_12 * vdt2;

            /* update velocity of particle "i" */
            *(__v4sf *)(vx + i + 0) = vxi_0 + ax_0 * vdt;
            *(__v4sf *)(vy + i + 0) = vyi_0 + ay_0 * vdt;
            *(__v4sf *)(vx + i + 4) = vxi_4 + ax_4 * vdt;
            *(__v4sf *)(vy + i + 4) = vyi_4 + ay_4 * vdt;
            *(__v4sf *)(vx + i + 8) = vxi_8 + ax_8 * vdt;
            *(__v4sf *)(vy + i + 8) = vyi_8 + ay_8 * vdt;
            *(__v4sf *)(vx + i + 12) = vxi_12 + ax_12 * vdt;
            *(__v4sf *)(vy + i + 12) = vyi_12 + ay_12 * vdt;
        }

        for (unsigned int i = n_particles & ~0xfu; i < n_particles; ++i) {
            float ax = 0.0f;
            float ay = 0.0f;

            for (unsigned int j = 0; j < n_particles; ++j) {
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

        if( !(*cb)(cb_arg, step) )
            break;

        std::swap(x, xn);
        std::swap(y, yn);
    }

    return step;
}

#ifdef __AVX__
extern "C" unsigned int
simulator_vec_16_0_8(
    const simulator_conf_t* conf,
    const simulator_data_t* data
    )
{
    LOAD_CONF( conf );
    LOAD_DATA( data );

    const float half = 0.5f;
    const float dt2 = dt * dt;

#   ifdef ACCURATE_VEC
    const __v8sf vone  = { 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f };
#   endif
    const __v8sf vhalf = { half, half, half, half, half, half, half, half };
    const __v8sf vdt = { dt, dt, dt, dt, dt, dt, dt, dt };
    const __v8sf vdt2 = vdt * vdt;

    unsigned int step;
    for (step = 0; step < n_steps; ++step) {
        for (unsigned int i = 0; i < n_particles - 15; i += 16) {
            __v8sf ax_0 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
            __v8sf ay_0 = ax_0;
            __v8sf ax_8 = ax_0;
            __v8sf ay_8 = ax_0;

            const __v8sf xi_0 = *(__v8sf *)(x + i + 0);
            const __v8sf yi_0 = *(__v8sf *)(y + i + 0);
            const __v8sf mi_0 = *(__v8sf *)(m + i + 0);
            const __v8sf qi_0 = *(__v8sf *)(q + i + 0);

            const __v8sf xi_8 = *(__v8sf *)(x + i + 8);
            const __v8sf yi_8 = *(__v8sf *)(y + i + 8);
            const __v8sf mi_8 = *(__v8sf *)(m + i + 8);
            const __v8sf qi_8 = *(__v8sf *)(q + i + 8);

            for (unsigned int j = 0; j < n_particles; ++j) {
                const __v8sf xj = { x[j], x[j], x[j], x[j], x[j], x[j], x[j], x[j] };
                const __v8sf yj = { y[j], y[j], y[j], y[j], y[j], y[j], y[j], y[j] };
                const __v8sf mj = { m[j], m[j], m[j], m[j], m[j], m[j], m[j], m[j] };
                const __v8sf qj = { q[j], q[j], q[j], q[j], q[j], q[j], q[j], q[j] };

                const __v8sf dx_0 = xj - xi_0;
                const __v8sf dy_0 = yj - yi_0;
                const __v8sf dx_8 = xj - xi_8;
                const __v8sf dy_8 = yj - yi_8;

#               ifdef ACCURATE_VEC
                const __v8sf invr_0 = vone / _mm256_sqrt_ps( dx_0 * dx_0 + dy_0 * dy_0 + vhalf );
                const __v8sf invr_8 = vone / _mm256_sqrt_ps( dx_8 * dx_8 + dy_8 * dy_8 + vhalf );
#               else
                const __v8sf invr_0 = _mm256_rsqrt_ps( dx_0 * dx_0 + dy_0 * dy_0 + vhalf );
                const __v8sf invr_8 = _mm256_rsqrt_ps( dx_8 * dx_8 + dy_8 * dy_8 + vhalf );
#               endif
                const __v8sf coef_0 = (mj - qi_0 * qj / mi_0) * invr_0 * invr_0 * invr_0;
                const __v8sf coef_8 = (mj - qi_8 * qj / mi_8) * invr_8 * invr_8 * invr_8;

                /* accumulate the acceleration from gravitational attraction */
                ax_0 += coef_0 * dx_0;
                ay_0 += coef_0 * dy_0;
                ax_8 += coef_8 * dx_8;
                ay_8 += coef_8 * dy_8;
            }

            const __v8sf vxi_0 = *(__v8sf *)(vx + i + 0);
            const __v8sf vyi_0 = *(__v8sf *)(vy + i + 0);
            const __v8sf vxi_8 = *(__v8sf *)(vx + i + 8);
            const __v8sf vyi_8 = *(__v8sf *)(vy + i + 8);

            /* update position of particle "i" */
            *(__v8sf *)(xn + i + 0) = xi_0 + vxi_0 * vdt + vhalf * ax_0 * vdt2;
            *(__v8sf *)(yn + i + 0) = yi_0 + vyi_0 * vdt + vhalf * ay_0 * vdt2;
            *(__v8sf *)(xn + i + 8) = xi_8 + vxi_8 * vdt + vhalf * ax_8 * vdt2;
            *(__v8sf *)(yn + i + 8) = yi_8 + vyi_8 * vdt + vhalf * ay_8 * vdt2;

            /* update velocity of particle "i" */
            *(__v8sf *)(vx + i + 0) = vxi_0 + ax_0 * vdt;
            *(__v8sf *)(vy + i + 0) = vyi_0 + ay_0 * vdt;
            *(__v8sf *)(vx + i + 8) = vxi_8 + ax_8 * vdt;
            *(__v8sf *)(vy + i + 8) = vyi_8 + ay_8 * vdt;
        }

        for (unsigned int i = n_particles & ~0xfu; i < n_particles; ++i) {
            float ax = 0.0f;
            float ay = 0.0f;

            for (unsigned int j = 0; j < n_particles; ++j) {
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

        if( !(*cb)(cb_arg, step) )
            break;

        std::swap(x, xn);
        std::swap(y, yn);
    }

    return step;
}
#endif


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

        if( !(*cb)(cb_arg, step) )
            break;

        std::swap(x, xn);
        std::swap(y, yn);
    }

    return step;
}

