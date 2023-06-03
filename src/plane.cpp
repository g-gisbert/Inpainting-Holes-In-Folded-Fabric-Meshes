#include "plane.hpp"


Plane::Plane(Vector3 _n) : n(_n), shift(0.) {}
Plane::Plane(Vector3 _n, double _shift) : n(_n), shift(_shift) {}

Vector3 Plane::projection(Vector3 p) const {
    return p - dot(n, p) / norm(n) * n + n*shift;
}
