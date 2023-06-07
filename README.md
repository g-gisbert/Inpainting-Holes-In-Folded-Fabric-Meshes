![Teaser](teaser.png)
# Inpainting Holes in Folded Fabric Meshes

This repository contains a research prototype implementation of the paper Inpainting Holes in Folded
Fabric Meshes, SMI 2023. 

The method proposes to fill holes in triangle mesh surfaces representing fabrics. The software allows
the user to provide a mesh with a hole and to fill it.

## Build Instructions

The project requires CGAL to be installed on your computer. It also relies on geometry-central, 
polyscope and libIGL as well as a C++17 compiler.
You can build the project as follows:
```
git clone https://github.com/g-gisbert/Inpainting-Holes-In-Folded-Fabric-Meshes
cd Inpainting-Holes-In-Folded-Fabric-Meshes
mkdir build
cd build
cmake ..
make -j8
```

## Run Instructions


