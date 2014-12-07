#ifndef NBODY_SIMULATORS_H
#define NBODY_SIMULATORS_H

typedef struct
{
    float* x;
    float* y;
    float* xn;
    float* yn;

    float* vx;
    float* vy;

    const float* m;
    const float* q;
} simulator_data_t;

#define LOAD_DATA( data )                       \
    float* x = (data)->x;                       \
    float* y = (data)->y;                       \
    float* xn = (data)->xn;                     \
    float* yn = (data)->yn;                     \
    float* const vx = (data)->vx;               \
    float* const vy = (data)->vy;               \
    const float* const m = (data)->m;           \
    const float* const q = (data)->q

typedef bool (*simulator_callback_t)(
    void* arg,
    unsigned int step,
    float* x,
    float* y,
    float* vx,
    float* vy
    );

typedef struct
{
    unsigned int n_particles;
    unsigned int n_steps;
    float dt;
    simulator_callback_t cb;
    void* cb_arg;
} simulator_conf_t;

#define LOAD_CONF( conf )                                       \
    const unsigned int n_particles = conf->n_particles;         \
    const unsigned int n_steps = conf->n_steps;                 \
    const float dt = conf->dt;                                  \
    const simulator_callback_t cb = conf->cb;                   \
    void* const cb_arg = conf->cb_arg;

typedef unsigned int (*simulator_t)(
    const simulator_conf_t* conf,
    const simulator_data_t* data
    );

simulator_t get_simulator(const char* name);

#endif /* !NBODY_SIMULATORS_H */
