#include "bvh.h"


struct Bins {
    static const int BIN_COUNT = 8;
    Bins() { memset(counts, 0, sizeof(uint32_t) * BIN_COUNT); }
    uint32_t counts[BIN_COUNT];
    AABB bounds[BIN_COUNT];
};


struct BVHBuildTask
{

static void execute_serially(BVH& bvh, uint32_t node_idx, uint32_t* start, uint32_t* end, uint32_t* temp) {
    uint32_t size = end - start;
    BVHNode& node = bvh.mNodes[node_idx];

    auto F = [&](const int lid, const int fid) {
        auto fh = bvh.msh.halfedge(pmp::Face(fid));
        if (lid == 0) {
            return bvh.msh.to_vertex(fh).idx();
        }
        else if (1 == lid) {
            return bvh.msh.to_vertex(bvh.msh.next_halfedge(fh)).idx();
        }
        else if (2 == lid) {
            return bvh.msh.to_vertex(bvh.msh.next_halfedge(bvh.msh.next_halfedge(fh))).idx();
        }
        else {
            assert(false);
            return pmp::Vertex().idx();
        }
    };

    auto V = [&](int axis, int nid) {
        auto p = bvh.msh.position(pmp::Vertex(nid));
        return p[axis];
    };

    auto Vcol = [&](const int nid) {
        auto p = bvh.msh.position(pmp::Vertex(nid));
        return Vector3f(p[0], p[1], p[2]);
    };


    Float best_cost = BVH::T_tri * size;
    int64_t best_index = -1, best_axis = -1;
    float* left_areas = (float*)temp;

    for (int axis = 0; axis < 3; ++axis) {
        {
            std::sort(start, end, [&](uint32_t f1, uint32_t f2) {
                return
                    (V(axis, F(0, f1)) + V(axis, F(1, f1)) + V(axis, F(2, f1))) <
                    (V(axis, F(0, f2)) + V(axis, F(1, f2)) + V(axis, F(2, f2)));
                });
        }

        AABB aabb;
        for (uint32_t i = 0; i < size; ++i) {
            uint32_t f = *(start + i);
            {
                aabb.expandBy(Vcol(F(0, f)));
                aabb.expandBy(Vcol(F(1, f)));
                aabb.expandBy(Vcol(F(2, f)));
            }
            left_areas[i] = (float)aabb.surfaceArea();
        }
        if (axis == 0)
            node.aabb = aabb;

        aabb.clear();

        Float tri_factor = BVH::T_tri / node.aabb.surfaceArea();
        for (uint32_t i = size - 1; i >= 1; --i) {
            uint32_t f = *(start + i);
            {
                aabb.expandBy(Vcol(F(0, f)));
                aabb.expandBy(Vcol(F(1, f)));
                aabb.expandBy(Vcol(F(2, f)));
            }

            float left_area = left_areas[i - 1];
            float right_area = aabb.surfaceArea();
            uint32_t prims_left = i;
            uint32_t prims_right = size - i;

            Float sah_cost = 2.0f * BVH::T_aabb +
                tri_factor * (prims_left * left_area +
                    prims_right * right_area);
            if (sah_cost < best_cost) {
                best_cost = sah_cost;
                best_index = i;
                best_axis = axis;
            }
        }
    }

    if (best_index == -1) {
        /* Splitting does not reduce the cost, make a leaf */
        node.leaf.flag = 1;
        node.leaf.start = start - bvh.mIndices;
        node.leaf.size = size;
        return;
    }


    {
        std::sort(start, end, [&](uint32_t f1, uint32_t f2) {
            return
                (V(best_axis, F(0, f1)) + V(best_axis, F(1, f1)) + V(best_axis, F(2, f1))) <
                (V(best_axis, F(0, f2)) + V(best_axis, F(1, f2)) + V(best_axis, F(2, f2)));
            });
    }

    uint32_t left_count = best_index;
    int node_idx_left = node_idx + 1;
    int node_idx_right = node_idx + 2 * left_count;
    node.inner.rightChild = node_idx_right;
    node.inner.unused = 0;

    execute_serially(bvh, node_idx_left, start, start + left_count, temp);
    execute_serially(bvh, node_idx_right, start + left_count, end, temp + left_count);
}

};


BVH::BVH(const pmp::SurfaceMesh& msh, const AABB& aabb)
    : mIndices(nullptr), mDiskRadius(0.f), msh(msh)
{
    {
        const auto nf = msh.n_faces();
        mNodes.resize(2 * nf);
        memset(mNodes.data(), 0, sizeof(BVHNode) * mNodes.size());
        mNodes[0].aabb = aabb;
        mIndices = new uint32_t[nf];
    }
}


bool BVH::rayIntersect(Ray ray, uint32_t& idx, Float& t, Vector2f* uv) const {
    if (mNodes.empty())
        return false;

    uint32_t node_idx = 0, stack[64];
    uint32_t stack_idx = 0;
    bool hit = false;
    t = std::numeric_limits<Float>::infinity();

    {
        while (true) {
            const BVHNode& node = mNodes[node_idx];

            if (!node.aabb.rayIntersect(ray)) {
                if (stack_idx == 0)
                    break;
                node_idx = stack[--stack_idx];
                continue;
            }

            if (node.isInner()) {
                stack[stack_idx++] = node.inner.rightChild;
                node_idx++;
                assert(stack_idx < 64);
            }
            else {
                Float _t;
                Vector2f _uv;
                for (uint32_t i = node.start(), end = node.end(); i < end; ++i) {
                    if (rayIntersectTri(ray, mIndices[i], _t, _uv)) {
                        idx = mIndices[i];
                        t = ray.maxt = _t;
                        hit = true;
                        if (uv)
                            *uv = _uv;
                    }
                }
                if (stack_idx == 0)
                    break;
                node_idx = stack[--stack_idx];
                continue;
            }
        }
    }


    return hit;
}

auto top(pmp::Point& p)
{
    return Vector3f(p[0], p[1], p[2]);
}

bool BVH::rayIntersectTri(const Ray& ray, uint32_t i, Float& t, Vector2f& uv) const {
    auto fh = msh.halfedge(pmp::Face(i));
    auto v0 = msh.to_vertex(fh);
    auto v1 = msh.to_vertex(msh.next_halfedge(fh));
    auto v2 = msh.to_vertex(msh.next_halfedge(msh.next_halfedge(fh)));
    
    auto pa = msh.position(v0);
    auto pb = msh.position(v1);
    auto pc = msh.position(v2);

    Vector3f p0 = top(pa);
    Vector3f p1 = top(pb);
    Vector3f p2 = top(pc);


    Vector3f edge1 = p1 - p0, edge2 = p2 - p0;
    Vector3f pvec = ray.d.cross(edge2);

    Float det = edge1.dot(pvec);
    if (det == 0.0f)
        return false;
    Float inv_det = 1.0f / det;

    Vector3f tvec = ray.o - p0;
    Float u = tvec.dot(pvec) * inv_det;
    if (u < 0.0f || u > 1.0f)
        return false;

    Vector3f qvec = tvec.cross(edge1);
    Float v = ray.d.dot(qvec) * inv_det;

    if (v < 0.0f || u + v > 1.0f)
        return false;

    Float tempT = edge2.dot(qvec) * inv_det;
    if (tempT < ray.mint || tempT > ray.maxt)
        return false;

    t = tempT;
    uv << u, v;
    return true;
}



BVH::~BVH() {
    delete[] mIndices;
}


void BVH::build()
{
    const auto total_size = msh.n_faces();

    uint32_t* temp = new uint32_t[total_size];

    BVHBuildTask::execute_serially(*this, 0u, mIndices, mIndices + total_size, temp);
    delete[] temp;
}


void buildscanline(pmp::SurfaceMesh& surface) {};
