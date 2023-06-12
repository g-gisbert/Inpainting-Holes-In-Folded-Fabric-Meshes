#include "application.hpp"


int main(int argc, char** argv) {

    if (argc != 2) {
        std::cerr << "Wrong number of arguments. The only argument is the data path." << std::endl;
        return EXIT_FAILURE;
    }
    std::string path(argv[1]);

    Application::init(path);

    return EXIT_SUCCESS;
}