#include <chrono>
#include <csignal>
#include <cstdlib>
#include <cstring> // strerror()
#include <ctime>

#include <stdexcept>
#include <iomanip>
#include <limits>
#include <sstream>

#include "NBodySim.h"
#include "bitmap_image.hpp"

/*************************************************************************/
/* Trimming functions                                                    */
/*************************************************************************/
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

/*************************************************************************/
/* Printing functions                                                    */
/*************************************************************************/
template <class T>
static inline void
print_aligned(std::ostream & os,
              const std::string & s,
              const T & val)
{
    static const int w = 25;
    static const char f = ' ';
    os << s << '=' << std::setfill(f) << std::setw(w - s.length()) << val << std::endl;
}

std::ostream &
operator<<(std::ostream & os, const NBodySim & s)
{
    print_aligned(os, NBodySim::CONF_KEYS::N_PARTICLES,    s.m_n_particles);
    print_aligned(os, NBodySim::CONF_KEYS::N_STEPS,        s.m_n_steps);
    print_aligned(os, NBodySim::CONF_KEYS::SIMULATOR,      s.m_simulator);
    print_aligned(os, NBodySim::CONF_KEYS::VISUAL,         s.m_visual);

    if (s.m_visual) {
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
    }

    print_aligned(os, NBodySim::CONF_KEYS::DUMP_FILE,      s.m_dumpfile);
    print_aligned(os, NBodySim::CONF_KEYS::SEED,           s.m_seed);
    return os;
}

/*************************************************************************/
/* NBodySim                                                              */
/*************************************************************************/
volatile bool NBodySim::interrupted = false;

void
NBodySim::register_signal(int signum)
{
    signal(signum, NBodySim::signal_handler);
}

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

bool
NBodySim::simulator_cb(void* arg,
                       unsigned int)
{
    NBodySim* self = (NBodySim*)arg;

    std::swap(self->m_x, self->m_xn);
    std::swap(self->m_y, self->m_yn);
    return !NBodySim::interrupted;
}

bool
NBodySim::simulator_cb_visual(void* arg,
                              unsigned int step)
{
    NBodySim* self = (NBodySim*)arg;

    const unsigned int n_particles = self->m_n_particles;
    const unsigned int img_width = self->m_img_width;
    const unsigned int img_height = self->m_img_height;

    float * x = self->m_x.data();
    float * y = self->m_y.data();
    float * vx = self->m_vx.data();
    float * vy = self->m_vy.data();

    for (unsigned int i = 0; i < n_particles; ++i) {
        if (x[i] < 0) {
            x[i] = -x[i];
            vx[i] = 0.5f * std::fabs(vx[i]);
        } else if (x[i] > img_width - 1) {
            x[i] = 2 * (img_width - 1) - x[i];
            vx[i] = -0.5f * std::fabs(vx[i]);
        }

        if (y[i] < 0) {
            y[i] = -y[i];
            vy[i] = 0.5f * std::fabs(vy[i]);
        } else if (y[i] > img_height - 1) {
            y[i] = 2 * (img_height - 1) - y[i];
            vy[i] = -0.5f * std::fabs(vy[i]);
        }
    }

    if (step % self->m_plot_every == 0)
        self->save_image(step / self->m_plot_every);

    self->print_status(step);
    return simulator_cb(self, step);
}

NBodySim::NBodySim(unsigned int n_particles,
                   unsigned int n_steps,
                   const std::string& simulator)
{
    load_default_settings();
    m_n_particles = n_particles;
    m_n_steps = n_steps;
    m_simulator = simulator;

    initialize();
}

NBodySim::NBodySim(std::istream & s)
{
    load_default_settings();
    load_settings(s);

    initialize();
}

void
NBodySim::load_default_settings()
{
    m_n_particles = 0;
    m_n_steps = 0;
    m_img_width = 800;
    m_img_height = 600;
    m_plot_every = 10;
    m_seed = time(nullptr);

    m_visual = false;

    m_time_step = 0.001f;
    m_max_initspeed = 0.0f;
    m_max_initmass = 1000.0f;
    m_max_initcharge = 1000.0f;
    m_min_initspeed = 0.0f;
    m_min_initmass = 1e-6f; //cannot be zero
    m_min_initcharge = -1000.0f;

    m_simulator = "naive";
}

bool
NBodySim::load_settings(std::istream & s)
{
    std::string line;

    try {
        while (std::getline(s, line)) {
            if(line.empty() || line[0] == '#')
                // skip empty lines and comments
                continue;

            std::istringstream iss(line);
            std::string key, value;

            if(!std::getline(iss, key, '=') ||
               !std::getline(iss, value))
                continue;

            trim(key);
            trim(value);

            try {
                if (key == CONF_KEYS::N_PARTICLES) {
                    m_n_particles = std::stoul(value);
                } else if (key == CONF_KEYS::N_STEPS) {
                    m_n_steps = std::stoul(value);
                } else if (key == CONF_KEYS::SIMULATOR) {
                    m_simulator = value;
                } else if (key == CONF_KEYS::VISUAL) {
                    m_visual = std::stoul(value);
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
    } catch (std::exception & ex) {
        if (!s.eof())
            throw;
    }

    return s || s.eof();
}

bool
NBodySim::load_particles(std::istream & s)
{
    unsigned int i = 0;
    std::string line;

    try {
        while (i < m_n_particles && std::getline(s, line)) {
            if(line.empty() || line[0] == '#')
                // skip empty lines and comments
                continue;

            std::istringstream iss(line);
            iss >> m_x[i]
                >> m_y[i]
                >> m_xn[i]
                >> m_yn[i]
                >> m_vx[i]
                >> m_vy[i]
                >> m_m[i]
                >> m_q[i];
            ++i;
        }
    } catch (std::exception & ex) {
        if (!s.eof())
            throw;
    }

    return s || s.eof();
}

bool
NBodySim::dump_particles(std::ostream & s) const
{
    s << "#x\ty\txnext\tynext\tvelx\tvely\tmass\tcharge" << std::endl;

    s << std::setprecision(10) << std::fixed;

    for (unsigned int i = 0; i < m_n_particles; ++i) {
        s << m_x[i]  << '\t'
          << m_y[i]  << '\t'
          << m_xn[i] << '\t'
          << m_yn[i] << '\t'
          << m_vx[i] << '\t'
          << m_vy[i] << '\t'
          << m_m[i]  << '\t'
          << m_q[i]  << '\t'
          << std::endl;
    }

    return s;
}

void
NBodySim::initialize()
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

    if( !(m_simfun = get_simulator( m_simulator.c_str() )) )
        throw std::runtime_error( "Simulator not found" );

    if (0 < m_n_particles) {
        init_array(m_x.begin(),  m_x.end(),  0,                m_img_width - 1);
        init_array(m_y.begin(),  m_y.end(),  0,                m_img_height - 1);
        init_array(m_vx.begin(), m_vx.end(), m_min_initspeed,  m_max_initspeed);
        init_array(m_vy.begin(), m_vy.end(), m_min_initspeed,  m_max_initspeed);
        init_array(m_m.begin(),  m_m.end(),  m_min_initmass,   m_max_initmass);
        init_array(m_q.begin(),  m_q.end(),  m_min_initcharge, m_max_initcharge);
    }
}

void
NBodySim::run_simulation()
{
    std::cout << "Running simulation with:" << std::endl;
    std::cout << *this;

    simulator_conf_t conf;
    memset( &conf, 0, sizeof(conf) );
    conf.n_particles = m_n_particles;
    conf.n_steps = m_n_steps;
    conf.img_width = m_img_width;
    conf.img_height = m_img_height;
    conf.dt = m_time_step;
    conf.cb = (m_visual ? simulator_cb_visual : simulator_cb);
    conf.cb_arg = this;

    simulator_data_t data;
    memset( &data, 0, sizeof(data) );
    data.x = m_x.data();
    data.y = m_y.data();
    data.xn = m_xn.data();
    data.yn = m_yn.data();
    data.vx = m_vx.data();
    data.vy = m_vy.data();
    data.m = m_m.data();
    data.q = m_q.data();

    auto begin = std::chrono::steady_clock::now();
    unsigned int step = (*m_simfun)(&conf, &data);

    if( m_visual ) {
        if (step % m_plot_every == 0) {
            save_image(step / m_plot_every);
        }
        print_status(step);
        std::cout << std::endl;
    } else {
        std::cout << "Processed: " << step << " steps" << std::endl;
    }

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

/*************************************************************************/
/* NBodySim::CONF_KEYS                                                   */
/*************************************************************************/
const std::string NBodySim::CONF_KEYS::N_PARTICLES    = "n_particles";
const std::string NBodySim::CONF_KEYS::N_STEPS        = "n_steps";
const std::string NBodySim::CONF_KEYS::SIMULATOR      = "simulator";
const std::string NBodySim::CONF_KEYS::VISUAL         = "visual";
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
