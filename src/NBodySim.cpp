#include <chrono>
#include <cmath>
#include <csignal>
#include <cstdlib>
#include <cstring> // strerror()
#include <ctime>
#include <iomanip>
#include <limits>
#include <sstream>
#include <xmmintrin.h>

#include "NBodySim.h"
#include "bitmap_image.hpp"

typedef float v4sf __attribute__((vector_size (16)));

volatile bool NBodySim::interrupted = false;

const std::string NBodySim::CONF_KEYS::N_PARTICLES    = "n_particles";
const std::string NBodySim::CONF_KEYS::N_STEPS        = "n_steps";
const std::string NBodySim::CONF_KEYS::IMG_WIDTH      = "img_width";
const std::string NBodySim::CONF_KEYS::IMG_HEIGHT     = "img_height";
const std::string NBodySim::CONF_KEYS::PLOT_EVERY     = "plot_every";
const std::string NBodySim::CONF_KEYS::TIME_STEP      = "time_step";
const std::string NBodySim::CONF_KEYS::MAX_INITSPEED  = "max_initspeed";
const std::string NBodySim::CONF_KEYS::MAX_INITMASS   = "max_initmass";
const std::string NBodySim::CONF_KEYS::MAX_INITCHARGE = "max_initcharge";
const std::string NBodySim::CONF_KEYS::MIN_INITSPEED  = "min_initspeed";
const std::string NBodySim::CONF_KEYS::MIN_INITMASS   = "min_initmass";
const std::string NBodySim::CONF_KEYS::MIN_INITCHARGE = "min_initcharge";
const std::string NBodySim::CONF_KEYS::SEED           = "rand_seed";
const std::string NBodySim::CONF_KEYS::IMG_PREFIX     = "img_prefix";
const std::string NBodySim::CONF_KEYS::DUMP_FILE      = "dump_file";


template <class ForwardIterator>
void
NBodySim::init_array(ForwardIterator first,
                     ForwardIterator last,
                     float min,
                     float max)
{
    if (min == max) {
        std::fill(first, last, min);
    } else {
        if(max < min)
            std::swap(min, max);

        const float norm = (max - min) / RAND_MAX;

        while (first != last) {
            *first = min + norm * rand();
            ++first;
        }
    }
}

void
NBodySim::signal_handler(int signum)
{
    NBodySim::interrupted = signum;
}


NBodySim::NBodySim(unsigned int n_particles,
                   unsigned int n_steps,
                   unsigned int seed)
: m_n_particles(n_particles),
  m_n_steps(n_steps),
  m_img_width(800),
  m_img_height(600),
  m_plot_every(10),
  m_seed(seed),
  m_time_step(0.001f),
  m_max_initspeed(0.0f),
  m_max_initmass(1000.0f),
  m_max_initcharge(1000.0f),
  m_min_initspeed(0.0f),
  m_min_initmass(1e-6f), //cannot be zero
  m_min_initcharge(-1000.0f)
{
    init_arrays();
}

const char * trimchars = " \t\n\r\f\v\"";

static inline std::string &
rtrim(std::string& s, const char * t = trimchars)
{
    s.erase(s.find_last_not_of(t) + 1);
    return s;
}

static inline std::string &
ltrim(std::string& s, const char * t = trimchars)
{
    s.erase(0, s.find_first_not_of(t));
    return s;
}

static inline std::string &
trim(std::string & s, const char * t = trimchars)
{
    return ltrim(rtrim(s, t), t);
}

NBodySim::NBodySim(std::istream & s)
: NBodySim()
{
    load_settings(s);
    init_arrays();
}

std::istream &
NBodySim::load_settings(std::istream & s)
{
    std::string line;
    try {
        while (std::getline(s, line)) {
            if (!line.empty() && line[0] == '#') continue; //skip comments

            std::istringstream iss(line);
            std::string key;

            if (std::getline(iss, key, '=')) {
                trim(key);
                std::string value;
                if (std::getline(iss, value)) {
                    trim(value);

                    try {
                        if (key == CONF_KEYS::N_PARTICLES) {
                            m_n_particles = std::stoul(value);
                        } else if (key == CONF_KEYS::N_STEPS) {
                            m_n_steps = std::stoul(value);
                        } else if (key == CONF_KEYS::IMG_WIDTH) {
                            m_img_width = std::stoul(value);
                        } else if (key == CONF_KEYS::IMG_HEIGHT) {
                            m_img_height = std::stoul(value);
                        } else if (key == CONF_KEYS::PLOT_EVERY) {
                            m_plot_every = std::stoul(value);
                        } else if (key == CONF_KEYS::TIME_STEP) {
                            m_time_step = std::stof(value);
                        } else if (key == CONF_KEYS::MAX_INITSPEED) {
                            m_max_initspeed = std::stof(value);
                        } else if (key == CONF_KEYS::MAX_INITMASS) {
                            m_max_initmass = std::stof(value);
                        } else if (key == CONF_KEYS::MAX_INITCHARGE) {
                            m_max_initcharge = std::stof(value);
                        } else if (key == CONF_KEYS::MIN_INITSPEED) {
                            m_min_initspeed = std::stof(value);
                        } else if (key == CONF_KEYS::MIN_INITMASS) {
                            m_min_initmass = std::stof(value);
                        } else if (key == CONF_KEYS::MIN_INITCHARGE) {
                            m_min_initcharge = std::stof(value);
                        } else if (key == CONF_KEYS::SEED) {
                            m_seed = std::stoul(value);
                        } else if (key == CONF_KEYS::IMG_PREFIX) {
                            m_img_prefix = value;
                        } else if (key == CONF_KEYS::DUMP_FILE) {
                            m_dumpfile = value;
                        } else {
                            std::cerr << "Ignoring unknown key: " << key << " = " << value << std::endl;
                        }
                    } catch (std::exception & ex) {
                        std::cout << "error" << std::endl;
                        std::cerr << "Could not use " << key << " = " << value
                                  << " pair, because " << value << " is not a valid value." << std::endl;
                    }
                }
            }
        }
    } catch (std::exception & ex) {
        if (errno != 0) throw ex; // do not throw on EOF
    }

    return s;
}

std::istream &
NBodySim::load_particles(std::istream & s)
{
    std::string line;
    try {
        for (unsigned int i = 0; i < m_n_particles && std::getline(s, line); ++i) {
            if (!line.empty() && line[0] == '#') continue; //skip comments

            std::istringstream iss(line);
            iss >> m_x[i]
                >> m_y[i]
                >> m_xn[i]
                >> m_yn[i]
                >> m_vx[i]
                >> m_vy[i]
                >> m_m[i]
                >> m_q[i];
        }
    } catch (std::exception & ex) {
        if (errno != 0) throw ex; // do not throw on EOF
    }

    return s;
}

std::ostream &
NBodySim::dump_particles(std::ostream & s) const
{
    s << "#x\ty\txnext\tynext\tvelx\tvely\tmass\tcharge" << std::endl;

    for (unsigned int i = 0; i < m_n_particles; ++i) {
        s << m_x[i] << '\t'
          << m_y[i] << '\t'
          << m_xn[i] << '\t'
          << m_yn[i] << '\t'
          << m_vx[i] << '\t'
          << m_vy[i] << '\t'
          << m_m[i] << '\t'
          << m_q[i] << '\t'
          << std::endl;
    }

    return s;
}

void
NBodySim::init_arrays()
{
    m_x.resize(m_n_particles);
    m_y.resize(m_n_particles);
    m_xn.resize(m_n_particles);
    m_yn.resize(m_n_particles);
    m_vx.resize(m_n_particles);
    m_vy.resize(m_n_particles);
    m_m.resize(m_n_particles);
    m_q.resize(m_n_particles);

    srand(m_seed);

    if (m_n_particles != 0) {
        NBodySim::init_array(m_x.begin(),  m_x.end(),  0,                m_img_width - 1);
        NBodySim::init_array(m_y.begin(),  m_y.end(),  0,                m_img_height - 1);
        NBodySim::init_array(m_vx.begin(), m_vx.end(), m_min_initspeed,  m_max_initspeed);
        NBodySim::init_array(m_vy.begin(), m_vy.end(), m_min_initspeed,  m_max_initspeed);
        NBodySim::init_array(m_m.begin(),  m_m.end(),  m_min_initmass,   m_max_initmass);
        NBodySim::init_array(m_q.begin(),  m_q.end(),  m_min_initcharge, m_max_initcharge);
    }
}

void
NBodySim::run_simulation()
{
    std::cout << "Running simulation with:" << std::endl;
    std::cout << *this;

#ifdef VISUAL
    const unsigned int img_height = m_img_height;
    const unsigned int img_width  = m_img_width;
    const unsigned int plot_every = m_plot_every;
#endif
    const unsigned int n_particles = m_n_particles;

    float * x  = m_x.data();
    float * y  = m_y.data();
    float * xn = m_xn.data();
    float * yn = m_yn.data();
    float * const vx = m_vx.data();
    float * const vy = m_vy.data();
    const float * const m = m_m.data();
    const float * const q = m_q.data();

    const float half = 0.5f;
    const float dt = m_time_step;
    const float dt2 = dt * dt;

    const v4sf vhalf = { half, half, half, half };
    const v4sf vdt = { dt, dt, dt, dt };
    const v4sf vdt2 = vdt * vdt;

    auto begin = std::chrono::steady_clock::now();

    unsigned int step;
    for (step = 0; step < m_n_steps && !NBodySim::interrupted; ++step) {
#       pragma omp parallel for schedule(static) default(none) \
            firstprivate(x, y, xn, yn, vx, vy, m, q, n_particles)
        for (unsigned int i = 0; i < n_particles - 3; i += 4) {
            v4sf ax = { 0.0f, 0.0f, 0.0f, 0.0f };
            v4sf ay = { 0.0f, 0.0f, 0.0f, 0.0f };

            const v4sf xi = *(v4sf *)(x + i);
            const v4sf yi = *(v4sf *)(y + i);
            const v4sf mi = *(v4sf *)(m + i);
            const v4sf qi = *(v4sf *)(q + i);

            for (unsigned int j = 0; j < n_particles; ++j) {
                const v4sf xj = { x[j], x[j], x[j], x[j] };
                const v4sf yj = { y[j], y[j], y[j], y[j] };
                const v4sf mj = { m[j], m[j], m[j], m[j] };
                const v4sf qj = { q[j], q[j], q[j], q[j] };

                const v4sf dx = xj - xi;
                const v4sf dy = yj - yi;

                const v4sf invr = __builtin_ia32_rsqrtps( dx * dx + dy * dy + vhalf );
                const v4sf coef = (mj - qi * qj / mi) * invr * invr * invr;

                /* accumulate the acceleration from gravitational attraction */
                ax += coef * dx;
                ay += coef * dy;
            }

            const v4sf vxi = *(v4sf *)(vx + i);
            const v4sf vyi = *(v4sf *)(vy + i);

            /* update position of particle "i" */
            *(v4sf *)(xn + i) = xi + vxi * vdt + vhalf * ax * vdt2;
            *(v4sf *)(yn + i) = yi + vyi * vdt + vhalf * ay * vdt2;

            /* update velocity of particle "i" */
            *(v4sf *)(vx + i) = vxi + ax * vdt;
            *(v4sf *)(vy + i) = vyi + ay * vdt;
        }

        for (unsigned int i = n_particles & ~0x3; i < n_particles; ++i) {
            float ax = 0.0f;
            float ay = 0.0f;

            for (unsigned int j = 0; j < n_particles; ++j) {
                const float dx = x[j] - x[i];
                const float dy = y[j] - y[i];
                const float invr = 1.0f / std::sqrt(dx * dx + dy * dy + half);
                const float coef = (m[j] - q[i] * q[j] / m[i])
                                                        * invr * invr * invr;

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

#ifdef VISUAL
#       pragma omp parallel for schedule(static) default(none) \
            firstprivate(xn, yn, vx, vy, n_particles, img_height, img_width, plot_every)
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

        if (step % plot_every == 0) {
            save_image(step / plot_every);
        }

        print_status(step);
#endif

        std::swap(x, xn);
        std::swap(y, yn);
        std::swap(m_x, m_xn);
        std::swap(m_y, m_yn);
    }

#ifdef VISUAL
    if (step % m_plot_every == 0) {
        save_image(step / m_plot_every);
    }
    print_status(step);
    std::cout << std::endl;
#else
    std::cout << "Processed: " << step << " steps" << std::endl;
#endif

    auto end = std::chrono::steady_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
    std::cout << "Elapsed time: "
              << std::fixed << std::setprecision(3) << elapsed / 1000.0
              << " s" << std::endl;


    if (!m_dumpfile.empty()) {
        if (m_dumpfile == "-") {
            std::cout << "Dumping particles to stdout:" << std::endl;
            dump_particles(std::cout);
        } else {
            std::ofstream f;
            f.exceptions(std::ios::failbit | std::ios::badbit);

            try {
                std::cout << "Dumping particles to file '" << m_dumpfile << "'..." << std::flush;
                f.open(m_dumpfile);
                dump_particles(f);
                f.close();
                std::cout << "success" << std::endl;
            } catch (std::exception & ex) {
                std::cout << "error" << std::endl;
                std::cerr << "Reason: " << strerror(errno) << std::endl;
            }
        }
    }
}

void
NBodySim::save_image(unsigned int seq) const
{
    bitmap_image image(m_img_width, m_img_height);
    image_drawer drawer(image);

    float max_val = std::max(std::max(std::fabs(m_max_initcharge),
                                      std::fabs(m_min_initcharge)),
                             m_max_initmass);
    float norm = 255.0f / max_val;

    image.set_all_channels(0);

    for (unsigned int p = 0; p < m_n_particles; ++p) {
        float xr = std::round(m_x[p]);
        float yr = std::round(m_y[p]);

        if (0 <= xr && xr <= m_img_width  - 1
         && 0 <= yr && yr <= m_img_height - 1) {
            if (m_q[p] >= 0.0f) {
                drawer.pen_color(std::round(m_q[p] * norm),
                                 std::round(m_m[p] * norm),
                                 0);
            } else if (m_q[p] < 0.0f) {
                drawer.pen_color(0,
                                 std::round(m_m[p] * norm),
                                 std::round(-m_q[p] * norm));
            }

            drawer.plot_pixel(xr, yr);
        }
    }

    std::ostringstream oss;
    oss << m_img_prefix << std::setfill('0') << std::setw(5) << seq << ".bmp";
    image.save_image(oss.str());
}


void
NBodySim::print_status(unsigned int step) const
{
    std::cout << '\r'
              << "step: " << step << '/' << m_n_steps << ' '
              << "plotted: " << 1 + step / m_plot_every << '/' << 1 + m_n_steps / m_plot_every
              << std::flush;
}

void
NBodySim::register_signal(int signum)
{
    signal(signum, NBodySim::signal_handler);
}

template <class T>
static inline void
print_aligned(std::ostream & os,
                   const std::string & s,
                   const T & val)
{
#ifdef VISUAL
    static const int w = 25;
#else
    static const int w = 17;
#endif
    static const char f = ' ';
    os << s << '=' << std::setfill(f) << std::setw(w - s.length()) << val << std::endl;
}

std::ostream &
operator<<(std::ostream & os, const NBodySim & s)
{
       print_aligned(os, NBodySim::CONF_KEYS::N_PARTICLES,    s.m_n_particles);
       print_aligned(os, NBodySim::CONF_KEYS::N_STEPS,        s.m_n_steps);
#ifdef VISUAL
       print_aligned(os, NBodySim::CONF_KEYS::IMG_WIDTH,      s.m_img_width);
       print_aligned(os, NBodySim::CONF_KEYS::IMG_HEIGHT,     s.m_img_height);
       print_aligned(os, NBodySim::CONF_KEYS::PLOT_EVERY,     s.m_plot_every);
       print_aligned(os, NBodySim::CONF_KEYS::TIME_STEP,      s.m_time_step);
       print_aligned(os, NBodySim::CONF_KEYS::MAX_INITSPEED,  s.m_max_initspeed);
       print_aligned(os, NBodySim::CONF_KEYS::MAX_INITMASS,   s.m_max_initmass);
       print_aligned(os, NBodySim::CONF_KEYS::MAX_INITCHARGE, s.m_max_initcharge);
       print_aligned(os, NBodySim::CONF_KEYS::MIN_INITSPEED,  s.m_min_initspeed);
       print_aligned(os, NBodySim::CONF_KEYS::MIN_INITMASS,   s.m_min_initmass);
       print_aligned(os, NBodySim::CONF_KEYS::MIN_INITCHARGE, s.m_min_initcharge);
       print_aligned(os, NBodySim::CONF_KEYS::IMG_PREFIX,     s.m_img_prefix);
#endif
       print_aligned(os, NBodySim::CONF_KEYS::DUMP_FILE,      s.m_dumpfile);
       print_aligned(os, NBodySim::CONF_KEYS::SEED,           s.m_seed);
    return os;
}
