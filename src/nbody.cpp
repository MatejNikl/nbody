#include <chrono>
#include <cmath>
#include <csignal>
#include <cstring> // strerror()
#include <fstream>
#include <iostream>
#include <sstream>
#include <unistd.h>

#include "NBodySim.h"


int
write_default_config(const std::string & fn);

void
print_help(const char * runcmnd);

bool
file_exists(const std::string & file);

bool
is_integer_only(const char * a,
                unsigned int & res);

int
main(int argc, char *argv[])
{
    NBodySim * s = nullptr;


    if (argc == 2 || argc == 3) {
        unsigned int n_particles;
        unsigned int n_steps;

        if (argc == 3 && is_integer_only(argv[1], n_particles)
                      && is_integer_only(argv[2], n_steps)) {
            s = new NBodySim(n_particles, n_steps);
        } else {
            if (!file_exists(argv[1])) return write_default_config(argv[1]);

            std::ifstream f;
            f.exceptions(std::ios::failbit | std::ios::badbit);

            try {
                std::cout << "Loading settings from config file '" << argv[1] << "'..." << std::flush;
                f.open(argv[1]);
                s = new NBodySim(f);
                f.close();
                std::cout << "success" << std::endl;

                if (argc == 3) {
                    std::cout << "Loading particles from file '" << argv[2] << "'..." << std::flush;
                    f.open(argv[2]);
                    s->load_particles(f);
                    f.close();
                    std::cout << "success" << std::endl;
                }
            } catch (const std::ifstream::failure & e) {
                std::cout << "error" << std::endl;
                std::cerr << "Reason: " << strerror(errno) << std::endl;
                return EXIT_FAILURE;
            }
        }
    } else {
        print_help(*argv);
        return EXIT_FAILURE;
    }

    NBodySim::register_signal(SIGINT);

    s->run_simulation();
    delete s;

    return EXIT_SUCCESS;
}

int
write_default_config(const std::string & fn)
{
    if (file_exists(fn)) {
        std::cout << "File '" << fn << "' already exists." << std::endl
                  << "Press Enter to continue, Ctrl+C to cancel..." << std::flush;
        if (std::cin.get() != '\n') return EXIT_FAILURE;
    }

    std::ofstream f;
    f.exceptions(std::ios::failbit | std::ios::badbit);

    try {
        std::cout << "Writing default config to '" << fn << "'..." << std::flush;
        f.open(fn);
        f << NBodySim(0, 0);
        f.close();
        std::cout << "success" << std::endl;
    } catch (const std::ifstream::failure & e) {
        std::cout << "error" << std::endl;
        std::cerr << "Reason: " << strerror(errno) << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

void
print_help(const char * runcmnd)
{
    std::cout << "Usage:" << std::endl
              << runcmnd << " #particles #steps" << std::endl
              << runcmnd << " non-existant file (writes default one)" << std::endl
              << runcmnd << " settings_file [particles_file]" << std::endl;
}

bool file_exists(const std::string & file)
{
    return access(file.c_str(), F_OK) != -1;
}

bool
is_integer_only(const char * a,
           unsigned int & res)
{
    char * end;
    res = std::strtoul(a, &end, 10);
    return *end == '\0';
}
