#include "sphere.hpp"


bool Sphere::isInside(const Vector3& p) noexcept {
    return norm(p - center) < r;
}


Vector3 Sphere::closestSurfacePoint(const Vector3& p) noexcept {
    return (p - center).normalize() * r + center;
}