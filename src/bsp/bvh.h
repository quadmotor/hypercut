#pragma once


#include <Eigen/Core>
#include <Eigen/Geometry>
#include "pmp/surface_mesh.h"


// adapt from famous
// https://github.com/wjakob/instant-meshes

typedef double Float;

typedef Eigen::Matrix<Float, 2, 1>                              Vector2f;
typedef Eigen::Matrix<Float, 3, 1>                              Vector3f;

struct Ray {
    Vector3f o, d;
    Float mint, maxt;

    Ray(const Vector3f& o, const Vector3f& d) :
        o(o), d(d), mint(0), maxt(std::numeric_limits<Float>::infinity()) { }

    Ray(const Vector3f& o, const Vector3f& d, Float mint, Float maxt) :
        o(o), d(d), mint(mint), maxt(maxt) { }

    Vector3f operator()(Float t) const { return o + t * d; }
};

struct AABB {
    Vector3f min, max;

    AABB() { clear(); }

    AABB(const Vector3f& min, const Vector3f& max) : min(min), max(max) {}

    void clear() {
        const Float inf = std::numeric_limits<Float>::infinity();
        min.setConstant(inf);
        max.setConstant(-inf);
    }

    void expandBy(const Vector3f& p) {
        min = min.cwiseMin(p);
        max = max.cwiseMax(p);
    }

    void expandBy(const AABB& aabb) {
        min = min.cwiseMin(aabb.min);
        max = max.cwiseMax(aabb.max);
    }

    bool contains(const Vector3f& p) {
        return (p.array() >= min.array()).all() &&
            (p.array() <= max.array()).all();
    }

    bool rayIntersect(const Ray& ray) const {
        Float nearT = -std::numeric_limits<Float>::infinity();
        Float farT = std::numeric_limits<Float>::infinity();

        for (int i = 0; i < 3; i++) {
            Float origin = ray.o[i];
            Float minVal = min[i], maxVal = max[i];

            if (ray.d[i] == 0) {
                if (origin < minVal || origin > maxVal)
                    return false;
            }
            else {
                Float t1 = (minVal - origin) / ray.d[i];
                Float t2 = (maxVal - origin) / ray.d[i];

                if (t1 > t2)
                    std::swap(t1, t2);

                nearT = std::max(t1, nearT);
                farT = std::min(t2, farT);

                if (!(nearT <= farT))
                    return false;
            }
        }

        return ray.mint <= farT && nearT <= ray.maxt;
    }

    Float squaredDistanceTo(const Vector3f& p) const {
        Float result = 0;
        for (int i = 0; i < 3; ++i) {
            Float value = 0;
            if (p[i] < min[i])
                value = min[i] - p[i];
            else if (p[i] > max[i])
                value = p[i] - max[i];
            result += value * value;
        }
        return result;
    }

    int largestAxis() const {
        Vector3f extents = max - min;

        if (extents[0] >= extents[1] && extents[0] >= extents[2])
            return 0;
        else if (extents[1] >= extents[0] && extents[1] >= extents[2])
            return 1;
        else
            return 2;
    }

    Float surfaceArea() const {
        Vector3f d = max - min;
        return 2.0f * (d[0] * d[1] + d[0] * d[2] + d[1] * d[2]);
    }

    Vector3f center() const {
        return 0.5f * (min + max);
    }

    static AABB merge(const AABB& aabb1, const AABB& aabb2) {
        return AABB(aabb1.min.cwiseMin(aabb2.min), aabb1.max.cwiseMax(aabb2.max));
    }
};




/* BVH node in 32 bytes */
struct BVHNode {
    union {
        struct {
            unsigned flag : 1;
            uint32_t size : 31;
            uint32_t start;
        } leaf;

        struct {
            uint32_t unused;
            uint32_t rightChild;
        } inner;
    };
    AABB aabb;

    inline bool isLeaf() const {
        return leaf.flag == 1;
    }

    inline bool isInner() const {
        return leaf.flag == 0;
    }

    inline bool isUnused() const {
        return inner.unused == 0 && inner.rightChild == 0;
    }

    inline uint32_t start() const {
        return leaf.start;
    }

    inline uint32_t end() const {
        return leaf.start + leaf.size;
    }
};

class BVH {
    friend struct BVHBuildTask;
    /* Cost values for BVH surface area heuristic */
    enum { T_aabb = 1, T_tri = 1 };
public:
    BVH(const pmp::SurfaceMesh &msh, const AABB& aabb);

    ~BVH();

    Float diskRadius() const { return mDiskRadius; }

    void build();

    bool rayIntersect(Ray ray) const;

    bool rayIntersect(Ray ray, uint32_t& idx, Float& t, Vector2f* uv = nullptr) const;


protected:
    bool rayIntersectTri(const Ray& ray, uint32_t i, Float& t, Vector2f& uv) const;
    bool rayIntersectDisk(const Ray& ray, uint32_t i, Float& t) const;
    void refitBoundingBoxes(uint32_t node_idx = 0);
    std::pair<Float, uint32_t> statistics(uint32_t node_idx = 0) const;

protected:
    std::vector<BVHNode> mNodes;
    uint32_t* mIndices;
    Float mDiskRadius;
    const pmp::SurfaceMesh& msh;
};


