#pragma once

// following are needed to properly compile vec3.h
#define __host__
#define __device__
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

    int MaximumExtent() const {
        vec3 d = Diagonal();
        if (d.x() > d.y() && d.x() > d.z())
            return 0;
        else if (d.y() > d.z())
            return 1;
        else
            return 2;
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
#else

Bounds3 Union(const Bounds3& b1, const Bounds3& b2);
Bounds3 Union(const Bounds3& b, const vec3& p);

#endif // GEOMETRY_IMPLEMENTATION

struct Triangle {
    vec3 v[3];
    float texCoords[6];
    unsigned char meshID;
    Bounds3 bounds;

    int s; // number of splits
    int bestIdx; // most important node that this triangle intersects

    Triangle() {}
    Triangle(const Triangle& tri, const Bounds3 bounds) : bounds(bounds) {
        for (int i = 0; i < 3; i++)
            v[i] = tri.v[i];
        meshID = tri.meshID;
        bestIdx = tri.bestIdx;
        for (auto i = 0; i < 6; i++)
            texCoords[i] = tri.texCoords[i];
    }
    Triangle(vec3 v0, vec3 v1, vec3 v2, float tc[6], unsigned char mID) {
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