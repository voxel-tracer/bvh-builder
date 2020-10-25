
#include <iostream>
#include <time.h>
#include <string>

#define TINYOBJLOADER_IMPLEMENTATION 
#include "tiny_obj_loader.h"

// following are needed to properly compile vec3.h
#define __host__
#define __device__

//#define SAH_BVH

#include "vec3.h"

struct mat3x3 {
    vec3 rows[3];

    vec3 mul(const vec3& v) const {
        return vec3(
            dot(rows[0], v),
            dot(rows[1], v),
            dot(rows[2], v)
        );
    }
};

const mat3x3 xUp = {
    vec3(0,-1,0),
    vec3(1,0,0),
    vec3(0,0,1)
};

const mat3x3 yUp = {
    vec3(1,0,0),
    vec3(0,1,0),
    vec3(0,0,1)
};

const mat3x3 zUp = {
    vec3(1,0,0),
    vec3(0,0,1),
    vec3(0,-1,0)
};

struct aabb {
    vec3 _min;
    vec3 _max;

    aabb() : _min(vec3(INFINITY, INFINITY, INFINITY)), _max(vec3(-INFINITY, -INFINITY, -INFINITY)) {}
    aabb(const aabb& node): _min(node._min), _max(node._max) {}

    void grow(vec3 v) {
        _min = min(_min, v);
        _max = max(_max, v);
    }

    void grow(const aabb &b) {
        _min = min(_min, b._min);
        _max = max(_max, b._max);
    }

    void intersect(const aabb& node) {
        _min = max(_min, node._min);
        _max = min(_max, node._max);
    }

    vec3 centroid() const {
        return (_max + _min) * 0.5;
    }

    unsigned int split_axis() const { return max_component(_max - _min); }

    // return true if this aabb's split plane splits the passed node
    bool doesItSplit(const aabb& node) const {
        int axis = split_axis();
        float middle = centroid()[axis];
        return middle > node._min[axis] && middle < node._max[axis];
    }

    vec3 size() const { return _max - _min; }

    aabb leftSplit() const {
        int axis = split_axis();
        float splitPos = centroid()[axis];

        aabb left(*this);
        left._max[axis] = splitPos;

        return left;
    }

    aabb rightSplit() const {
        int axis = split_axis();
        float splitPos = centroid()[axis];

        aabb right(*this);
        right._min[axis] = splitPos;

        return right;
    }

    bool contains(const aabb& node) const {
        return node._min.x() >= _min.x() && node._min.y() >= _min.y() && node._min.z() >= _min.z() &&
            node._max.x() <= _max.x() && node._max.y() <= _max.y() && node._max.z() <= _max.z();
    }
};

struct btriangle {
    vec3 v[3];
    float texCoords[6];
    unsigned char meshID;
    aabb bounds;

    int s; // number of splits
    int bestIdx; // most important node that this triangle intersects

    btriangle() {}
    btriangle(const btriangle& tri, const aabb& node) {
        v[0] = tri.v[0];
        v[1] = tri.v[1];
        v[2] = tri.v[2];
        meshID = tri.meshID;
        bestIdx = tri.bestIdx;

        for (auto i = 0; i < 6; i++)
            texCoords[i] = tri.texCoords[i];

        bounds = node;
    }

    btriangle(vec3 v0, vec3 v1, vec3 v2, float tc[6], unsigned char mID) {
        v[0] = v0;
        v[1] = v1;
        v[2] = v2;
        meshID = mID;

        for (auto i = 0; i < 6; i++)
            texCoords[i] = tc[i];

        bounds.grow(v0);
        bounds.grow(v1);
        bounds.grow(v2);
    }
};

struct triangle {
    triangle() {}
    triangle(const btriangle &t) {
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

#ifdef SAH_BVH
struct sah_aabb {
    vec3 bmin;
    vec3 bmax;

    sah_aabb(): bmin(vec3(FLT_MAX, FLT_MAX, FLT_MAX)), bmax(vec3(-FLT_MAX, -FLT_MAX, -FLT_MAX)) {}
    sah_aabb(const triangle& t) {
        bmin = min(t.v[0], min(t.v[1], t.v[2]));
        bmax = max(t.v[0], max(t.v[1], t.v[2]));
    }
    sah_aabb(vec3 min, vec3 max) :bmin(min), bmax(max) {}

    unsigned int split_axis() const { return max_component(bmax - bmin); }

    float area() const {
        vec3 size = bmax - bmin;
        return size.x() * size.y() * size.z();
    }
};

sah_aabb merge(const sah_aabb& a, const sah_aabb& b) {
    return sah_aabb(min(a.bmin, b.bmin), max(a.bmax, b.bmax));
}

struct sah_bvh_node {
    vec3 min;
    vec3 max;

    uint32_t start;
    uint32_t length;

    sah_bvh_node* left;
    sah_bvh_node* right;

    sah_bvh_node(const sah_aabb& aabb, sah_bvh_node* left, sah_bvh_node* right) : min(aabb.bmin), max(aabb.bmax), start(0), length(0), left(left), right(right) {}
    sah_bvh_node(const sah_aabb& aabb, uint32_t start, uint32_t length) : min(aabb.bmin), max(aabb.bmax), start(start), length(length), left(NULL), right(NULL) {}
};
#endif

struct scene {
    btriangle* tris;
    int numTris;

    aabb* bvh;
    int bvh_size;

    vec3 bMin;
    vec3 bMax;

    void release() {
        delete[] tris;
        delete[] bvh;
    }
};

int center_x_compare(const void* a, const void* b) {
    float xa = ((btriangle*)a)->bounds.centroid().x();
    float xb = ((btriangle*)b)->bounds.centroid().x();

    if (xa < xb) return -1;
    else if (xb < xa) return 1;
    return 0;
}

int center_y_compare(const void* a, const void* b) {
    float ya = ((btriangle*)a)->bounds.centroid().y();
    float yb = ((btriangle*)b)->bounds.centroid().y();

    if (ya < yb) return -1;
    else if (yb < ya) return 1;
    return 0;
}

int center_z_compare(const void* a, const void* b) {
    float za = ((btriangle*)a)->bounds.centroid().z();
    float zb = ((btriangle*)b)->bounds.centroid().z();

    if (za < zb) return -1;
    else if (zb < za) return 1;
    return 0;
}

int min(int a, int b) {
    return a < b ? a : b;
}

int max(int a, int b) {
    return a > b ? a : b;
}

uint64_t build_bvh(aabb* nodes, int idx, btriangle* l, int n, int m, int numPrimitivesPerLeaf) {
    aabb node;
    aabb centroid_bounds;
    for (int i = 0; i < m; i++) {
        node.grow(l[i].bounds);
        centroid_bounds.grow(l[i].bounds.centroid());
    }
    nodes[idx] = node;

    if (m > numPrimitivesPerLeaf) {
        const unsigned int axis = centroid_bounds.split_axis();
        if (axis == 0)
            qsort(l, m, sizeof(btriangle), center_x_compare);
        else if (axis == 1)
            qsort(l, m, sizeof(btriangle), center_y_compare);
        else
            qsort(l, m, sizeof(btriangle), center_z_compare);

        // split the primitives such that at most n/2 are on the left of the split and the rest are on the right
        // given we have m primitives, left will get min(n/2, m) and right gets max(0, m - n/2)
        return 1 + 
            build_bvh(nodes, idx * 2, l, n / 2, min(n / 2, m), numPrimitivesPerLeaf) +
            build_bvh(nodes, idx * 2 + 1, l + n / 2, n / 2, max(0, m - (n / 2)), numPrimitivesPerLeaf);
    } else {
        return 1;
    }
}

float computeEmptySurfaceArea(const btriangle &tri) {
    // compute aabb surface area
    vec3 boundsExtent = tri.bounds.size();
    float Aaabb = 2 * (boundsExtent.x() * boundsExtent.y() + boundsExtent.x() * boundsExtent.z() + boundsExtent.y() * boundsExtent.z());

    // compute ideal surface area
    vec3 d1 = tri.v[1] - tri.v[0];
    vec3 d2 = tri.v[2] - tri.v[0];
    vec3 cr = cross(d1, d2);
    float Aideal = abs(cr.x()) + abs(cr.y()) + abs(cr.z());

    float value = (Aaabb - Aideal);
    return value;
}

float computeNodeImportance(const aabb &triBounds, int triIdx, const aabb* bvh, int bvhSize, int &bestIdx) {
    int level = log2f(bvhSize);
    int best = level;
    // we know that a bvh of size bvhSize has half its nodes just for the last level (one node per triangle bounds)
    // => index of first triangle's bvh node = bvhSize / 2
    int idx = bvhSize/2 + triIdx;

    bestIdx = -1;
    while (idx > 0) {
        if (bvh[idx].doesItSplit(triBounds)) {
            bestIdx = idx;
            best = level;
        }
        level--;
        idx >>= 1;
    }

    if (bestIdx == -1)
        std::cerr << " triangle " << triIdx << " doesn't have a best node" << std::endl;

    return powf(2, -best);
}

float computePriority(const btriangle &tri, int triIdx, const aabb* bvh, int bvhSize, int &bestIdx) {
    float importance = computeNodeImportance(tri.bounds, triIdx, bvh, bvhSize, bestIdx);
    float emptyArea = computeEmptySurfaceArea(tri);

    float X = 2.0f;
    float Y = 1.0f / 3.0f;
    float priority = powf(powf(X, importance) * emptyArea, Y);
    return priority;
}

float sumf(const float* priority, int size) {
    float s = 0;
    for (int t = 0; t < size; t++)
        s += priority[t];
    return s;
}

int sumi(const float* priority, float D, int size) {
    int s = 0;
    for (int t = 0; t < size; t++)
        s += floor(D * priority[t]);
    return s;
}

float computeSplitD(const float* priority, int size, int Sbudget) {
    // start by estimating Dmin, Dmax
    float Dmin = Sbudget / sumf(priority, size);
    int Smin = sumi(priority, Dmin, size);
    float Dmax = Dmin * Sbudget / Smin;
    int Smax = sumi(priority, Dmax, size);
    std::cerr << "Dmin = " << Dmin << " (" << Smin << "). Dmax = " << Dmax << " (" << Smax << ")" << std::endl;
    // use bissection and iterate 6 times
    while (true) {
        float Davg = (Dmin + Dmax) / 2;
        int Savg = sumi(priority, Davg, size);
        if (Savg == Smin || Savg == Smax) break;

        if (Savg > Sbudget) {
            Smax = Savg;
            Dmax = Davg;
        }
        else {
            Smin = Savg;
            Dmin = Davg;
        }
        std::cerr << "Dmin = " << Dmin << " (" << Smin << "). Dmax = " << Dmax << " (" << Smax << ")" << std::endl;
    }

    return Dmin;
}

int computeNumSplits(btriangle* tris, int size, const aabb* bvh, int bvhSize, int Sbudget) {
    // compute and store split priority for each triangle
    float* priority = new float[size];
    for (int t = 0; t < size; t++)
        priority[t] = computePriority(tris[t], t, bvh, bvhSize, tris[t].bestIdx);
    // compute split D
    float D = computeSplitD(priority, size, Sbudget);
    // compute num splits for each triangle
    int numSplits = 0;
    for (int t = 0; t < size; t++) {
        tris[t].s = floor(D * priority[t]);
        numSplits += tris[t].s;
    }

    delete[] priority;

    return numSplits;
}

bool nodeSplitIntersects(const aabb& node, const aabb& bounds) {
    int dim = node.split_axis();
    float pos = node.centroid()[dim];
    return bounds._min[dim] < pos && bounds._max[dim] > pos;
}

aabb findBestNode(const aabb& node, const aabb& bounds) {
    aabb cur = node;
    while (!nodeSplitIntersects(cur, bounds)) {
        aabb left = cur.leftSplit();
        if (left.contains(bounds))
            cur = left;
        else
            cur = cur.rightSplit();
    }
    return cur;
}

float clamp(float v, float vmin, float vmax) {
    return fminf(vmax, fmaxf(vmin, v));
    //return min(vmax, max(vmin, v));
}

void split(const btriangle& tri, const aabb& bestNode, aabb& left, aabb& right) {
    int dim = bestNode.split_axis();
    float pos = bestNode.centroid()[dim];

    vec3 v1 = tri.v[2];
    for (int i = 0; i < 3; i++) {
        vec3 v0 = v1;
        v1 = tri.v[i]; 
        float v1p = v1[dim];
        float v0p = v0[dim];

        if (v0p <= pos) left.grow(v0);
        if (v0p >= pos) right.grow(v0);
        // check if edge intersects the plane and store it in both sides
        if ((v0p < pos && v1p > pos) || (v0p > pos && v1p < pos)) {
            vec3 t = lerp(v0, v1, clamp((pos - v0p) / (v1p - v0p), 0.0f, 1.0f));
            left.grow(t);
            right.grow(t);
        }
    }
    // intersect original bounds in case tri is a split already
    left._max[dim] = pos;
    right._min[dim] = pos;
    left.intersect(tri.bounds);
    right.intersect(tri.bounds);

    // validate that left and right bounds are inside triangle bounds and their union produces triangle bounds
    if (!tri.bounds.contains(left))
        std::cerr << "tri.bounds !contains left" << std::endl;
    if (!tri.bounds.contains(right))
        std::cerr << "tri.bounds !contains right" << std::endl;
    aabb bounds(left);
    bounds.grow(right);
    if (tri.bounds._min != bounds._min || tri.bounds._max != bounds._max)
        std::cerr << "tri.bounds != union(left, right)" << std::endl;
}

void updateSplitCount(int s, const aabb& left, const aabb& right, int &sleft, int &sright) {
    float wa = left.size()[left.split_axis()];
    float wb = right.size()[right.split_axis()];

    sleft = floor((s - 1) * wa / (wa + wb) + 0.5f);
    sright = s - 1 - sleft;
}

aabb* build_bvh(btriangle* l, unsigned int numLeaves, int& bvhSize, bool splitTriangles, btriangle** tris, int &numTrisWithSplits) {
    std::cout << "numLeaves: " << numLeaves << std::endl;
    // number of leaves that is a power of 2, this is the max width of a complete binary tree
    const int pow2NumLeaves = (int)powf(2.0f, ceilf(log2f(numLeaves)));
    std::cout << "pow2NumLeaves: " << pow2NumLeaves << std::endl;
    // total number of nodes in the tree
    bvhSize = pow2NumLeaves * 2;
    std::cout << "bvh_size: " << bvhSize << std::endl;
    // allocate enough nodes to hold the whole tree, even if some of the nodes will remain unused
    aabb* bvh = new aabb[bvhSize];
    uint64_t numNodes = build_bvh(bvh, 1, l, pow2NumLeaves, numLeaves, 1);
    std::cerr << "num internal nodes = " << numNodes << std::endl;

    // should we split ? and can we actually split ?
    if (!splitTriangles || pow2NumLeaves == numLeaves) {
        return bvh; // do not split
    }
    // start splitting the triangles
    std::cerr << "start splitting triangles. We will add " << (pow2NumLeaves - numLeaves) << " extra splits" << std::endl;
    // compute for each triangle its num splits and most important node idx
    int numSplits = computeNumSplits(l, numLeaves, bvh, bvhSize, pow2NumLeaves - numLeaves);
    std::cerr << "adding an additional " << numSplits << " splits" << std::endl;
    // allocate enough triangles to store all triangles including the splits
    int total = numLeaves + numSplits;
    btriangle* l2 = new btriangle[total];
    for (int i = 0; i < numLeaves; i++) {
        l2[i] = l[i];
    }

    // actually split all triangles according to their numSplit (.s)
    int nextSplit = numLeaves;
    for (int i = 0; i < total; i++) {
        btriangle& tri = l2[i];
        // keep splitting this triangle until we run out of splits
        while (tri.s > 0) {
            //std::cerr << "splitting triangle " << i << " with " << tri.s << " splits" << std::endl;
            // find most important node that splits this triangle, starting from its bestIdx
            aabb node = findBestNode(bvh[tri.bestIdx], tri.bounds);
            aabb left, right;
            split(tri, node, left, right);
            // replace current triangle with left split
            tri.bounds = left;
            // add new split to the end of the array
            btriangle split(tri, right);
            // update split count for left and right splits
            int sleft, sright;
            updateSplitCount(tri.s, left, right, sleft, sright);
            tri.s = sleft;
            split.s = sright;
            // store new split at end of array
            l2[nextSplit++] = split;

            if (left.size()[left.split_axis()] < 0.00001f)
                std::cerr << " left bound empty" << std::endl;
            if (right.size()[right.split_axis()] < 0.00001f)
                std::cerr << " right bound empty" << std::endl;
        }
    }

    if (nextSplit != total) std::cerr << "nextSplit < total : " << nextSplit << " < " << total << std::endl;

    // now rebuild the bvh using all triangles including the splits
    std::cerr << "rebuilding bvh with splits" << std::endl;
    numNodes = build_bvh(bvh, 1, l2, pow2NumLeaves, total, 1);
    std::cerr << "num internal nodes = " << numNodes << std::endl;

    *tris = l2;
    numTrisWithSplits = total;

    //* tris = l2;
    //numTrisWithSplits = numLeaves;
    return bvh;
}

#ifdef SAH_BVH
sah_bvh_node* build_sah_bvh(triangle* tris, int n, uint32_t index, int triStart, sah_aabb *boxes, float* left_area, float* right_area, int &numInternals, int &numLeaves, uint32_t &maxIndex) {
    maxIndex = max(maxIndex, index);

    sah_aabb main_box;
    for (auto i = 0; i < n; i++) {
        sah_aabb box(tris[i]);
        main_box = merge(main_box, box);
    }

    if (n <= 5) {
        return new sah_bvh_node(main_box, triStart, n);
    }

    // find splitting axis
    int axis = main_box.split_axis();
    if (axis == 0)
        qsort(tris, n, sizeof(triangle), bmin_x_compare);
    else if (axis == 1)
        qsort(tris, n, sizeof(triangle), bmin_y_compare);
    else
        qsort(tris, n, sizeof(triangle), bmin_z_compare);

    // collect bounds of all triangles
    for (auto i = 0; i < n; i++) {
        boxes[i] = sah_aabb(tris[i]);
    }

    // compute left_area[i] = area of bounding box of all triangles [0, i]
    sah_aabb left_box;
    for (auto i = 0; i < n; i++) {
        left_box = merge(left_box, sah_aabb(tris[i]));
        left_area[i] = left_box.area();
    }

    // compute right_area[i] = area of bounding box of all triangles [i, n[
    sah_aabb right_box;
    for (auto i = n - 1; i > 0; i--) {
        right_box = merge(right_box, sah_aabb(tris[i]));
        right_area[i] = right_box.area();
    }

    // find partition with minimal SAH
    float minSAH = FLT_MAX;
    int minSAHidx = 0;
    for (auto i = 0; i < n - 1; i++) {
        float SAH = i * left_area[i] + (n - i - 1) * right_area[i + 1];
        if (SAH < minSAH) {
            minSAHidx = i;
            minSAH = SAH;
        }
    }

    float mainSAH = n * main_box.area();
    if (mainSAH <= minSAH) {
        return new sah_bvh_node(main_box, triStart, n);
    }
        
    // best split is [0, minSAHidx], ]minSAHidx, n[
    sah_bvh_node* left;
    if (minSAHidx == 0) {// leaf node
        left = new sah_bvh_node(boxes[0], triStart, 1);
        numLeaves++;
    } else {
        left = build_sah_bvh(tris, minSAHidx + 1, index * 2, triStart, boxes, left_area, right_area, numInternals, numLeaves, maxIndex);
        numInternals++;
    }

    sah_bvh_node* right;
    if (minSAHidx == n - 2) {// leaf right node
        right = new sah_bvh_node(boxes[n - 1], triStart + n - 1, 1);
        numLeaves++;
    } else {
        right = build_sah_bvh(tris + minSAHidx + 1, n - minSAHidx - 1, index * 2 + 1, triStart + minSAHidx + 1, boxes, left_area, right_area, numInternals, numLeaves, maxIndex);
        numInternals++;
    }
    return new sah_bvh_node(main_box, left, right);
}

void clean(sah_bvh_node* node) {
    if (node->left) clean(node->left);
    if (node->right) clean(node->right);
    delete node;
}

void build_sah_bvh(triangle *tris, int n) {
    sah_aabb* boxes = new sah_aabb[n];
    float* left_area = new float[n];
    float* right_area = new float[n];

    int numInternals = 0;
    int numLeaves = 0;
    uint32_t maxIndex = 0;

    std::cerr << "Building BVH using SAH algorithm..." << std::endl;
    sah_bvh_node* root = build_sah_bvh(tris, n, 1, 0, boxes, left_area, right_area, numInternals, numLeaves, maxIndex);
    std::cout << " max index " << maxIndex << std::endl;
    std::cout << " numInternals " << numInternals << std::endl;
    std::cout << " numLeaves " << numLeaves << std::endl;

    clean(root);
}
#endif

// center scene around origin
scene initScene(const std::vector<btriangle> &tris, float scale, bool centerAndScale, bool splitTriangles) {
    scene sc;

    int size = tris.size();
    sc.numTris = size;
    sc.tris = new btriangle[sc.numTris];

    for (int i = 0; i < size; i++) {
        sc.tris[i] = tris[i];
    }

    // center scene around origin
    aabb bounds;
    for (int i = 0; i < size; i++) {
        bounds.grow(sc.tris[i].bounds);
    }
    vec3 ctr = bounds.centroid(); // this is the model center
    ctr[1] = bounds._min[1]; // make sure we can put the floor at y = 0

    std::cerr << " original model bounds:" << std::endl;
    std::cerr << "  min: " << bounds._min << std::endl;
    std::cerr << "  max: " << bounds._max << std::endl;

    if (centerAndScale) {
        // find max size across all axes
        const float maxSize = max(bounds._max - bounds._min);
        if (scale == 0) scale = maxSize;
        // we want to normalize the model so that its maxSize is scale and centered around ctr
        // for each vertex v:
        //  v = v- ctr // ctr is new origin
        //  v = v / maxSize // scale model to fit in a bbox with maxSize 1
        //  v = v * scale // scale model so that maxSize = scale
        // => v = (v - ctr) * scale/maxSize
        aabb new_bounds;
        vec3 vs[3];
        for (int t = 0; t < size; t++) {
            for (int i = 0; i < 3; i++) {
                vs[i] = (sc.tris[t].v[i] - ctr) * scale / maxSize;
            }
            sc.tris[t] = btriangle(vs[0], vs[1], vs[2], sc.tris[t].texCoords, sc.tris[t].meshID);
            new_bounds.grow(sc.tris[t].bounds);
        }

        bounds = new_bounds;
        std::cerr << " updated model bounds:" << std::endl;
        std::cerr << "  min: " << bounds._min << std::endl;
        std::cerr << "  max: " << bounds._max << std::endl;
    }

    sc.bMin = bounds._min;
    sc.bMax = bounds._max;

    // build bvh
#ifdef SAH_BVH
    build_sah_bvh(sc.tris, size);
#else
    btriangle* trisWithSplits;
    int numTrisWithSplits;
    sc.bvh = build_bvh(sc.tris, size, sc.bvh_size, splitTriangles, &trisWithSplits, numTrisWithSplits);
    // replace tris with trisWithSplits
    delete[] sc.tris;
    sc.tris = trisWithSplits;
    sc.numTris = numTrisWithSplits;
#endif // SAH_BVH

    return sc;
}

bool loadFromObj(const std::string& filepath, const mat3x3& mat, std::vector<btriangle> &tris, unsigned char meshID) {
    tinyobj::attrib_t attrib;
    std::vector<tinyobj::shape_t> shapes;
    std::vector<tinyobj::material_t> materials;

    std::string warn;
    std::string err;

    bool ret = tinyobj::LoadObj(&attrib, &shapes, &materials, &warn, &err, filepath.c_str());

    if (!warn.empty())
        std::cerr << warn << std::endl;
    if (!err.empty())
        std::cerr << err << std::endl;

    if (!ret)
        return false;

    std::cerr << " num vertices " << attrib.vertices.size() << std::endl;
    if (!materials.empty())
        std::cerr << " materials size " << materials.size() << std::endl;
    if (!attrib.texcoords.empty())
        std::cerr << " texcoord size " << attrib.texcoords.size() << std::endl;
    if (!attrib.colors.empty())
        std::cerr << " colors size " << attrib.colors.size() << std::endl;
    if (!attrib.normals.empty())
        std::cerr << " normals size " << attrib.normals.size() << std::endl;

    // loop over shapes and copy all triangles to vertices vector
    uint32_t triIdx = 0;
    vec3 verts[3];
    float tc[6];
    for (auto s = 0; s < shapes.size(); s++) {
        // loop over faces
        size_t index_offset = 0;
        for (auto f = 0; f < shapes[s].mesh.num_face_vertices.size(); f++, triIdx++) {
            int fv = shapes[s].mesh.num_face_vertices[f];
            if (fv != 3)
                std::cerr << "face " << f << " of shape " << s << " has " << fv << " vertices" << std::endl;

            // loop over vertices in the face
            for (auto v = 0; v < 3; v++) {
                // access to vertex
                tinyobj::index_t idx = shapes[s].mesh.indices[index_offset + v];
                tinyobj::real_t vx = attrib.vertices[3 * idx.vertex_index + 0];
                tinyobj::real_t vy = attrib.vertices[3 * idx.vertex_index + 1];
                tinyobj::real_t vz = attrib.vertices[3 * idx.vertex_index + 2];
                tinyobj::real_t tx = attrib.texcoords[2 * idx.texcoord_index + 0];
                tinyobj::real_t ty = attrib.texcoords[2 * idx.texcoord_index + 1];

                vec3 tri(vx, vy, vz);
                verts[v] = mat.mul(tri);
                tc[v * 2 + 0] = tx;
                tc[v * 2 + 1] = ty;
            }
            tris.push_back(btriangle(verts[0], verts[1], verts[2], tc, meshID));
            index_offset += 3;
        }
    }

    return true;
}

// 00.04: no more triangle.center
void save(const std::string output, const scene& sc, int numPrimitivesPerLeaf) {
    // convert btriangle to triangle
    triangle* tris = new triangle[sc.numTris];
    for (int i = 0; i < sc.numTris; i++) {
        tris[i] = triangle(sc.tris[i]);
    }

    std::fstream out(output, std::ios::out | std::ios::binary);
    const char* HEADER = "BVH_00.04";
    out.write(HEADER, strlen(HEADER) + 1);
    out.write((char*)&sc.numTris, sizeof(int));
    out.write((char*)tris, sizeof(triangle) * sc.numTris);
    out.write((char*)&sc.bvh_size, sizeof(int));
    out.write((char*)sc.bvh, sizeof(aabb) * sc.bvh_size);
    out.write((char*)&sc.bMin, sizeof(vec3));
    out.write((char*)&sc.bMax, sizeof(vec3));
    out.write((char*)&numPrimitivesPerLeaf, sizeof(int));
    out.close();

    delete[] tris;
}

int main() {
    float scale = 100.0f;
    mat3x3 mat = yUp;
    bool splitTriangles = true;

    //TODO include scale in the transformation mat, so we can scale models separately
    std::string basePath = "C:\\Users\\adene\\models\\glsl-assets\\staircase\\";
    std::vector<btriangle> tris;
    bool success = true;

    success = success && loadFromObj(basePath + "Black.obj", yUp, tris, 0);
    success = success && loadFromObj(basePath + "Brass.obj", yUp, tris, 1);
    success = success && loadFromObj(basePath + "BrushedAluminium.obj", yUp, tris, 2);
    success = success && loadFromObj(basePath + "Candles.obj", yUp, tris, 3);
    success = success && loadFromObj(basePath + "ChairSeat.obj", yUp, tris, 4);
    success = success && loadFromObj(basePath + "Glass.obj", yUp, tris, 5);
    success = success && loadFromObj(basePath + "Gold.obj", yUp, tris, 6);
    success = success && loadFromObj(basePath + "Lampshade.obj", yUp, tris, 7);
    success = success && loadFromObj(basePath + "MagnoliaPaint.obj", yUp, tris, 8);
    success = success && loadFromObj(basePath + "Painting1.obj", yUp, tris, 9);
    success = success && loadFromObj(basePath + "Painting2.obj", yUp, tris, 10);
    success = success && loadFromObj(basePath + "Painting3.obj", yUp, tris, 11);
    success = success && loadFromObj(basePath + "StainlessSteel.obj", yUp, tris, 12);
    success = success && loadFromObj(basePath + "Wallpaper.obj", yUp, tris, 13);
    success = success && loadFromObj(basePath + "WhitePaint.obj", yUp, tris, 14);
    success = success && loadFromObj(basePath + "WhitePlastic.obj", yUp, tris, 15);
    success = success && loadFromObj(basePath + "WoodChair.obj", yUp, tris, 16);
    success = success && loadFromObj(basePath + "WoodFloor.obj", yUp, tris, 17);
    success = success && loadFromObj(basePath + "WoodLamp.obj", yUp, tris, 18);
    success = success && loadFromObj(basePath + "WoodStairs.obj", yUp, tris, 19);
    
    //success = success && loadFromObj("D:\\models\\obj\\cube\\cube.obj", yUp, tris, 0);
    if (!success) {
        std::cerr << "Failed to load obj files" << std::endl;
        return -1;
    }
    std::cerr << "read " << tris.size() << " triangles" << std::endl;
    scene s = initScene(tris, scale, false, splitTriangles); // passing scale=0 disables scaling the model

#ifndef SAH_BVH
    //save("D:\\models\\obj\\cube.bvh", s, numPrimitivesPerLeaf);
    save("C:\\Users\\adene\\models\\BVH\\staircase.bvh", s, 1);
#endif // !SAH_BVH

    delete[] s.tris;
}
