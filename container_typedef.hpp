#ifndef CONTAINER_TYPEDEF
#define CONTAINER_TYPEDEF
#include <unordered_map>
#include <vector>

// typedef of vectors of numbers
typedef std::vector<double> VectorDouble;
typedef std::vector<int> VectorInt;

// typedef of unordered maps
typedef std::unordered_map<int, int> MapIntInt;

// nested vectors for integration
typedef std::vector<double> Vector1D;
typedef std::vector<Vector1D> Vector2D;
typedef std::vector<Vector2D> Vector3D;
typedef std::vector<Vector3D> Vector4D;

// nested map and vectors for integration
typedef std::unordered_map<int, std::unordered_map<int, Vector1D>> MapIntIntVector1D;
typedef std::unordered_map<int, std::unordered_map<int, Vector2D>> MapIntIntVector2D;
typedef std::unordered_map<int, std::unordered_map<int, Vector3D>> MapIntIntVector3D;

#endif
