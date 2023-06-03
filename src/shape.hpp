#ifndef INC_2DCONTOUR_SHAPE_HPP
#define INC_2DCONTOUR_SHAPE_HPP


class Shape {
public:
    virtual ~Shape() = default;
    virtual bool isInside(const Vector3& p) noexcept = 0;
    virtual Vector3 closestSurfacePoint(const Vector3& p) noexcept = 0;

protected:
    static int count;
};


#endif //INC_2DCONTOUR_SHAPE_HPP
