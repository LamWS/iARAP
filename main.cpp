#include <iostream>
#include "mesh.h"
#include "Optimizer.h"
#include "path_config.h"

int main(int argc, char *argv[]) {
    std::cout << argv[0] << std::endl;
    if (argc < 2) {
        std::cout << "Usage: " << argv[0] << " <input mesh>" << std::endl;
        return 1;
    }
    Mesh mesh;
    mesh.read_tetramesh(std::string(argv[1]));
    Optimizer optimizer(mesh, output_path);
    optimizer.solve(mesh, 0.01, 1000);
    return 0;
}
