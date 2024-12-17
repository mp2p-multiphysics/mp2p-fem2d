#ifndef CONTAINER_TYPEDEF
#define CONTAINER_TYPEDEF
#include <unordered_map>
#include <vector>

namespace FEM2D
{

// nested int vectors
typedef std::vector<int> VectorInt;
typedef std::vector<VectorInt> VectorInt2D;

// nested double vectors
typedef std::vector<double> VectorDouble;
typedef std::vector<VectorDouble> VectorDouble2D;
typedef std::vector<VectorDouble2D> VectorDouble3D;
typedef std::vector<VectorDouble3D> VectorDouble4D;

// typedef of unordered maps
typedef std::unordered_map<int, int> MapIntInt;

}

#endif
