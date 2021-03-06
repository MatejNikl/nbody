#ifndef NBODY_NBODYSIM_H
#define NBODY_NBODYSIM_H

#include <iostream>
#include <string>
#include <vector>

#include <immintrin.h>

#include "aligned_allocator.hpp"
#include "simulators.h"

class NBodySim {
    friend std::ostream& operator<<(std::ostream& os, const NBodySim& s);

public:
    static void register_signal(int signum);

    NBodySim(unsigned int n_particles,
             unsigned int n_steps,
             const std::string& simulator);
    NBodySim(std::istream & s);

    void load_default_settings();
    bool load_settings(std::istream & s);

    bool load_particles(std::istream & s);
    bool dump_particles(std::ostream & s) const;

    void initialize();
    void run_simulation();
    void save_image(unsigned int seq) const;
    void print_status(unsigned int step) const;

private:
    struct CONF_KEYS {
        static const std::string N_PARTICLES;
        static const std::string N_STEPS;
        static const std::string SIMULATOR;
        static const std::string VISUAL;
        static const std::string IMG_WIDTH;
        static const std::string IMG_HEIGHT;
        static const std::string PLOT_EVERY;
        static const std::string TIME_STEP;
        static const std::string MAX_INITSPEED;
        static const std::string MAX_INITMASS;
        static const std::string MAX_INITCHARGE;
        static const std::string MIN_INITSPEED;
        static const std::string MIN_INITMASS;
        static const std::string MIN_INITCHARGE;
        static const std::string SEED;
        static const std::string IMG_PREFIX;
        static const std::string DUMP_FILE;
    };

    template <class ForwardIterator>
    static void init_array(ForwardIterator first,
                           ForwardIterator last,
                           float min,
                           float max);
    static void signal_handler(int signum);

    static bool simulator_cb(void* arg,
                             unsigned int step);
    static bool simulator_cb_visual(void* arg,
                                    unsigned int step);

    unsigned int m_n_particles;
    unsigned int m_n_steps;
    unsigned int m_img_width;
    unsigned int m_img_height;
    unsigned int m_plot_every;
    unsigned int m_seed;

    bool m_visual;

    float m_time_step;
    float m_max_initspeed;
    float m_max_initmass;
    float m_max_initcharge;
    float m_min_initspeed;
    float m_min_initmass;
    float m_min_initcharge;

#ifdef __AVX__
    typedef __v8sf align_type;
#else
    typedef __v4sf align_type;
#endif
    std::vector<float, aligned_allocator<float, sizeof(align_type)> > m_x;
    std::vector<float, aligned_allocator<float, sizeof(align_type)> > m_y;
    std::vector<float, aligned_allocator<float, sizeof(align_type)> > m_xn;
    std::vector<float, aligned_allocator<float, sizeof(align_type)> > m_yn;
    std::vector<float, aligned_allocator<float, sizeof(align_type)> > m_vx;
    std::vector<float, aligned_allocator<float, sizeof(align_type)> > m_vy;
    std::vector<float, aligned_allocator<float, sizeof(align_type)> > m_m;
    std::vector<float, aligned_allocator<float, sizeof(align_type)> > m_q;

    std::string m_img_prefix;
    std::string m_dumpfile;

    std::string m_simulator;
    simulator_t m_simfun;

    static volatile bool interrupted;
};

#endif //NBODY_NBODYSIM_H
