#include "cylinder.hpp"

bool Cylinder::isInside(const Vector3& p) noexcept {
    const Vector3 pos = p - center;
    Vector3 projection = dot(pos, up) * up + center;
    return norm(p-projection) < r && norm(center-projection) < h_2;
}

Vector3 Cylinder::closestSurfacePoint(const Vector3 &p) noexcept {
    const Vector3 pos = p - center;
    Vector3 projection = dot(pos, up) * up + center;
    return (p-projection).normalize() * r + projection;
}
