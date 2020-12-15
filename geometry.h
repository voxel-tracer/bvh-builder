#pragma once

// following are needed to properly compile vec3.h
//#define __host__
//#define __device__
#include "vec3.h"

struct Bounds3 {
    vec3 pMin;
    vec3 pMax;

    Bounds3() {
        float minNum = std::numeric_limits<float>::lowest();
        float maxNum = std::numeric_limits<float>::max();
        pMin = vec3(maxNum, maxNum, maxNum);
        pMax = vec3(minNum, minNum, minNum);
    }

    vec3 Diagonal() const { return pMax - pMin; }

    vec3 Center() const { return (pMax + pMin) / 2; }

    int MaximumExtent() const {
        vec3 d = Diagonal();
        if (d.x() > d.y() && d.x() > d.z())
            return 0;
        else if (d.y() > d.z())
            return 1;
        else
            return 2;
    }

    float SurfaceArea() const {
        vec3 d = Diagonal();
        return 2 * (d.x() * d.y() + d.x() * d.z() + d.y() * d.z());
    }

    bool MedianSplits(Bounds3 bounds) const {
        int dim = MaximumExtent();
        float pos = Center()[dim];
        return bounds.pMin[dim] < pos && bounds.pMax[dim] > pos;
    }

    Bounds3 leftSplit() const {
        int dim = MaximumExtent();
        float pos = Center()[dim];

        Bounds3 left = *this;
        left.pMax[dim] = pos;

        return left;
    }

    Bounds3 rightSplit() const {
        int dim = MaximumExtent();
        float pos = Center()[dim];

        Bounds3 right = *this;
        right.pMin[dim] = pos;

        return right;
    }

    bool Contains(Bounds3 node) const {
        //for (int a = 0; a < 3; a++) {
        //    if (node.pMax[a] < pMin[a] || node.pMin[a] > pMax[a])
        //        return false;
        //}
        //return true;
        return node.pMin.x() >= pMin.x() && node.pMin.y() >= pMin.y() && node.pMin.z() >= pMin.z() &&
            node.pMax.x() <= pMax.x() && node.pMax.y() <= pMax.y() && node.pMax.z() <= pMax.z();
    }

    bool IsValid() const {
        for (int a = 0; a < 3; a++)
            if (pMin[a] > pMax[a]) return false;
        return true;
    }
};

#ifdef GEOMETRY_IMPLEMENTATION

Bounds3 Union(const Bounds3& b1, const Bounds3& b2) {
    Bounds3 ret;
    ret.pMin = min(b1.pMin, b2.pMin);
    ret.pMax = max(b1.pMax, b2.pMax);
    return ret;
}

Bounds3 Union(const Bounds3& b, const vec3& p) {
    Bounds3 ret;
    ret.pMin = min(b.pMin, p);
    ret.pMax = max(b.pMax, p);
    return ret;
}

Bounds3 Intersect(const Bounds3& b1, const Bounds3& b2) {
    Bounds3 ret = b1;
    // _min must be in [node._min, node._max]
    ret.pMin = max(ret.pMin, b2.pMin);
    ret.pMin = min(ret.pMin, b2.pMax);
    // _max must be in [node._min, node._max]
    ret.pMax = max(ret.pMax, b2.pMin);
    ret.pMax = min(ret.pMax, b2.pMax);
    return ret;
}

bool operator!=(const Bounds3& b1, const Bounds3& b2) {
    return b1.pMin != b2.pMin || b1.pMax != b2.pMax;
}
#else

Bounds3 Union(const Bounds3& b1, const Bounds3& b2);
Bounds3 Union(const Bounds3& b, const vec3& p);
Bounds3 Intersect(const Bounds3& b1, const Bounds3& b2);
bool operator!=(const Bounds3& b1, const Bounds3& b2);

#endif // GEOMETRY_IMPLEMENTATION

struct Triangle {
    int id;
    vec3 v[3];
    float texCoords[6];
    unsigned char meshID;
    Bounds3 bounds;

    int s; // number of splits
    int bestIdx; // most important node that this triangle intersects

    Triangle() {}
    Triangle(const Triangle& tri, const Bounds3 bounds) : bounds(bounds), id(tri.id) {
        for (int i = 0; i < 3; i++)
            v[i] = tri.v[i];
        meshID = tri.meshID;
        bestIdx = tri.bestIdx;
        for (auto i = 0; i < 6; i++)
            texCoords[i] = tri.texCoords[i];
    }
    Triangle(int id, vec3 v0, vec3 v1, vec3 v2, float tc[6], unsigned char mID): id(id) {
        v[0] = v0;
        v[1] = v1;
        v[2] = v2;
        meshID = mID;

        for (auto i = 0; i < 6; i++)
            texCoords[i] = tc[i];

        bounds = Union(bounds, v0);
        bounds = Union(bounds, v1);
        bounds = Union(bounds, v2);
    }
};

// file structs compatible with my current renderer, may remove this in the future
struct LinearTriangle {
    LinearTriangle() {}
    LinearTriangle(const Triangle& t) {
        v[0] = t.v[0];
        v[1] = t.v[1];
        v[2] = t.v[2];
        meshID = t.meshID;
        for (auto i = 0; i < 6; i++)
            texCoords[i] = t.texCoords[i];
    }

    vec3 v[3];
    float texCoords[6];
    unsigned char meshID;
};

// this is duplicated with bvh.cpp, but it's fine for now
struct LinearBVHNode {
    __host__ __device__ LinearBVHNode() {}

    Bounds3 bounds;
    union {
        int primitivesOffset;   // leaf
        int firstChildOffset;  // interior
    };
    uint16_t nPrimitives;   // 0 -> interior node
    uint8_t axis;           // interior node: xyz
    uint8_t pad[1];         // ensure 32 bytes total size
};
