#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <vector>
#include "bsp_util.h"


int main(int argc, char* argv[])
{

    std::vector< std::string > args;
    std::string inputfn, outputfn;
    bool help = false;
    bool israndomsphere = true;
    bool ismesh2bsp = false;
    bool ismill = false;
    int resolution = 30;

    try {
        for (int i = 1; i < argc; i++) {
            if (strcmp("--randomsphere", argv[i]) == 0) {
                israndomsphere = true;
            }
            else if (strcmp("--mesh2bsp", argv[i]) == 0) {
                israndomsphere = false;
                ismesh2bsp = true;
            }
            else if (strcmp("--mill", argv[i]) == 0) {
                israndomsphere = false;
                ismill = true;
            }
            else if (strcmp("--millresolution", argv[i]) == 0) {
                if (++i >= argc) {
                    std::cerr << "Missing mill resolution argument!" << std::endl;
                    return -1;
                }
                resolution = atoi(argv[i]);
            }
            else if (strcmp("--output", argv[i]) == 0) {
                if (++i >= argc) {
                    std::cerr << "missing output argument" << std::endl;
                    return -1;
                }

                outputfn = argv[i];
            }
            else if (strcmp("--input", argv[i]) == 0) {
                if (++i >= argc) {
                    std::cerr << "missing input argument" << std::endl;
                    return -1;
                }

                inputfn = argv[i];
            }
            else
                args.push_back(argv[i]);
        }

    }
    catch (const std::exception& e) {
        std::cout << "Error: " << e.what() << std::endl;
        help = true;
    }

    if (!israndomsphere) {
        if (inputfn.empty() || outputfn.empty()) {
            help = true;
        }
    }

    if (!args.empty() || help) {
        std::cout << "--randomsphere" << std::endl;
        std::cout << "--mesh2bsp --input inputfilepath.obj --output resultfilepath.obj" << std::endl;
        std::cout << "--mill --millresolution 30 --input inputfilepath.obj --output resultfilepath.obj" << std::endl;


        std::cout << "current argc == " << argc << '\n';

        for (int ndx{}; ndx != argc; ++ndx)
            std::cout << "argv[" << ndx << "] == " << std::quoted(argv[ndx]) << '\n';
        std::cout << "argv[" << argc << "] == "
            << static_cast<void*>(argv[argc]) << '\n';

        return -1;
    }

    if (israndomsphere) {
        random_sphere(outputfn, 1000);
    }
    else {
        Mesh sh;
        sh.load(inputfn);
        mesh2bsp(sh, ismill, true, outputfn.c_str());
    }

    return 0;
}