#include <iostream>
#include "mesh.h"
#include "Optimizer.h"
int main() {
    Mesh mesh;
    mesh.read_tetramesh("../data/cube5.msh");
    Optimizer optimizer(mesh);
    optimizer.solve(mesh, 0.01, 1000);
    return 0;
}
