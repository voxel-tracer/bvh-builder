
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

struct triangle {
    triangle() {}
    triangle(vec3 v0, vec3 v1, vec3 v2, float tc[6], unsigned char mID) {
        v[0] = v0;
        v[1] = v1;
        v[2] = v2;
        meshID = mID;

        for (auto i = 0; i < 6; i++)
            texCoords[i] = tc[i];
    }

    vec3 center() const {
        //return (v[0] + v[1] + v[2]) / 3;
        return (bounds_min() + bounds_max()) / 2;
    }

    vec3 bounds_min() const {
        return min(v[0], min(v[1], v[2]));
    }

    vec3 bounds_max() const {
        return max(v[0], max(v[1], v[2]));
    }

    vec3 v[3];
    float texCoords[6];
    unsigned char meshID;
};

/**
non-leaf nodes represent bounding box of all spheres that are inside it
leaf nodes represent 2 spheres
*/
struct bvh_node {
    bvh_node() {}
    bvh_node(const vec3& min, const vec3& max) :min(min), max(max) {}

    unsigned int split_axis() const { return max_component(max - min); }

    void markLeaf() {
        min[0] += 1000;
    }

    float volume() const {
        float volume = 1.0f;
        vec3 extent = max - min;
        for (size_t i = 0; i < 3; i++)
            if (extent[i] > 0.0f) volume *= extent[i];
        return volume;
    }

    vec3 min;
    vec3 max;

private:
    void swap(int idx) {
        float tmp = min[idx];
        min[idx] = max[idx];
        max[idx] = tmp;
    }
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
    triangle* tris;
    int numTris;

    bvh_node* bvh;
    int bvh_size;
    int numLevels;

    vec3 bMin;
    vec3 bMax;

    void release() {
        delete[] tris;
        delete[] bvh;
    }
};

vec3 minof(const triangle* l, int n) {
    vec3 m(INFINITY, INFINITY, INFINITY);
    for (int t = 0; t < n; t++) // for each triangle
        m = min(m, l[t].bounds_min());
    return m;
}

vec3 center_minof(const triangle* l, int n) {
    vec3 m(INFINITY, INFINITY, INFINITY);
    for (int t = 0; t < n; t++) // for each triangle
        m = min(m, l[t].center());
    return m;
}

vec3 maxof(const triangle* l, int n) {
    vec3 m(-INFINITY, -INFINITY, -INFINITY);
    for (int t = 0; t < n; t++) // for each triangle
        m = max(m, l[t].bounds_max());
    return m;
}

vec3 center_maxof(const triangle* l, int n) {
    vec3 m(-INFINITY, -INFINITY, -INFINITY);
    for (int t = 0; t < n; t++) // for each triangle
        m = max(m, l[t].center());
    return m;
}

int bmin_x_compare(const void* a, const void* b) {
    float xa = ((triangle*)a)->bounds_min().x();
    float xb = ((triangle*)b)->bounds_min().x();

    if (xa < xb) return -1;
    else if (xb < xa) return 1;
    return 0;
}

int center_x_compare(const void* a, const void* b) {
    float xa = ((triangle*)a)->center().x();
    float xb = ((triangle*)b)->center().x();

    if (xa < xb) return -1;
    else if (xb < xa) return 1;
    return 0;
}

int bmin_y_compare(const void* a, const void* b) {
    float ya = ((triangle*)a)->bounds_min().y();
    float yb = ((triangle*)b)->bounds_min().y();

    if (ya < yb) return -1;
    else if (yb < ya) return 1;
    return 0;
}

int center_y_compare(const void* a, const void* b) {
    float ya = ((triangle*)a)->center().y();
    float yb = ((triangle*)b)->center().y();

    if (ya < yb) return -1;
    else if (yb < ya) return 1;
    return 0;
}

int bmin_z_compare(const void* a, const void* b) {
    float za = ((triangle*)a)->bounds_min().z();
    float zb = ((triangle*)b)->bounds_min().z();

    if (za < zb) return -1;
    else if (zb < za) return 1;
    return 0;
}

int center_z_compare(const void* a, const void* b) {
    float za = ((triangle*)a)->center().z();
    float zb = ((triangle*)b)->center().z();

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

struct bvh_result {
    int numNodes = 0;
    float cost = 0.0f;
    bvh_node* node;
};

bvh_result build_bvh(bvh_node* nodes, int idx, triangle* l, int n, int m, float It) {
    bvh_node node = bvh_node(minof(l, m), maxof(l, m));
    if ((node.min[0] >= 500) || ((node.min[0] + 1000) < 500)) std::cerr << "node.min.x " << node.min[0] << std::endl;

    nodes[idx] = node;
    bvh_result result = {
        1,
        m * It, // no split cost
        nodes + idx
    };

    if (m > 1) {
        bvh_node centroid_bounds(center_minof(l, m), center_maxof(l, m));
        const unsigned int axis = centroid_bounds.split_axis();
        if (axis == 0)
            qsort(l, m, sizeof(triangle), center_x_compare);
        else if (axis == 1)
            qsort(l, m, sizeof(triangle), center_y_compare);
        else
            qsort(l, m, sizeof(triangle), center_z_compare);

        // split the primitives such that at most n/2 are on the left of the split and the rest are on the right
        // given we have m primitives, left will get min(n/2, m) and right gets max(0, m - n/2)
        bvh_result left = build_bvh(nodes, idx * 2, l, n / 2, min(n / 2, m), It);
        bvh_result right = build_bvh(nodes, idx * 2 + 1, l + n / 2, n / 2, max(0, m - (n / 2)), It);

        // compute SAH split cost = 2 * cost of traversing each child node + cost of each node scaled by the probability to hit the node
        float splitCost = 2 + left.cost * left.node->volume() / result.node->volume() + right.cost * right.node->volume() / result.node->volume();
        if (result.cost <= splitCost) { // do not split
            result.node->markLeaf();
        } else { // split
            result.numNodes += left.numNodes;
            result.numNodes += right.numNodes;
            result.cost = splitCost;
        }
    } else if (m > 0) {
        result.node->markLeaf();
    }

    return result;
}

bvh_node* build_bvh(triangle* l, unsigned int numPrimitives, int& bvh_size, int &numLevels, float It) {
    const int numLeaves = numPrimitives;
    std::cout << "numLeaves: " << numLeaves << std::endl;
    // number of leaves that is a power of 2, this is the max width of a complete binary tree
    numLevels = ceilf(log2f(numLeaves));
    const int pow2NumLeaves = (int)powf(2.0f, numLevels);
    std::cout << "pow2NumLeaves: " << pow2NumLeaves << std::endl;
    // total number of nodes in the tree
    bvh_size = pow2NumLeaves * 2;
    std::cout << "bvh_size: " << bvh_size << std::endl;
    // allocate enough nodes to hold the whole tree, even if some of the nodes will remain unused
    bvh_node* nodes = new bvh_node[bvh_size];
    bvh_result result = build_bvh(nodes, 1, l, pow2NumLeaves, numPrimitives, It);
    std::cout << "num internal nodes created = " << result.numNodes << std::endl;

    return nodes;
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
scene initScene(const std::vector<triangle> &tris, float scale, bool centerAndScale, float It) {
    scene sc;

    // copy triangles to sc
    sc.numTris = tris.size();
    sc.tris = new triangle[sc.numTris];
    for (int i = 0; i < sc.numTris; i++) {
        sc.tris[i] = tris[i];
    }

    // center scene around origin
    vec3 mn = minof(sc.tris, sc.numTris);
    vec3 mx = maxof(sc.tris, sc.numTris);
    vec3 ctr = (mx + mn) / 2; // this is the model center
    ctr[1] = mn[1]; // make sure we can put the floor at y = 0

    std::cerr << " original model bounds:" << std::endl;
    std::cerr << "  min: " << mn << std::endl;
    std::cerr << "  max: " << mx << std::endl;

    if (centerAndScale) {
        // find max size across all axes
        const float maxSize = max(mx - mn);
        if (scale == 0) scale = maxSize;
        // we want to normalize the model so that its maxSize is scale and centered around ctr
        // for each vertex v:
        //  v = v- ctr // ctr is new origin
        //  v = v / maxSize // scale model to fit in a bbox with maxSize 1
        //  v = v * scale // scale model so that maxSize = scale
        // => v = (v - ctr) * scale/maxSize
        for (int t = 0; t < sc.numTris; t++) {
            for (int i = 0; i < 3; i++) {
                vec3 v = sc.tris[t].v[i];
                sc.tris[t].v[i] = (v - ctr) * scale / maxSize;
            }
        }

        mn = minof(sc.tris, sc.numTris);
        mx = maxof(sc.tris, sc.numTris);
        std::cerr << " updated model bounds:" << std::endl;
        std::cerr << "  min: " << mn << std::endl;
        std::cerr << "  max: " << mx << std::endl;
    }

    sc.bMin = mn;
    sc.bMax = mx;

    // build bvh
#ifdef SAH_BVH
    build_sah_bvh(sc.tris, size);
#else
    sc.bvh = build_bvh(sc.tris, sc.numTris, sc.bvh_size, sc.numLevels, It);
#endif // SAH_BVH

    return sc;
}

bool loadFromObj(const std::string& filepath, const mat3x3& mat, std::vector<triangle> &tris, unsigned char meshID) {
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

    //std::cerr << " num vertices " << attrib.vertices.size() << std::endl;
    //if (!materials.empty())
    //    std::cerr << " materials size " << materials.size() << std::endl;
    //if (!attrib.texcoords.empty())
    //    std::cerr << " texcoord size " << attrib.texcoords.size() << std::endl;
    //if (!attrib.colors.empty())
    //    std::cerr << " colors size " << attrib.colors.size() << std::endl;
    //if (!attrib.normals.empty())
    //    std::cerr << " normals size " << attrib.normals.size() << std::endl;

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
            tris.push_back(triangle(verts[0], verts[1], verts[2], tc, meshID));
            index_offset += 3;
        }
    }

    return true;
}

// 00.04: no more triangle.center
// 00.05: SAH based node merging, - numPrimitivesPerLeaf, + numLevels
// 00.06: add 500 to bvh_node.min[0] if its a leaf
void save(const std::string output, const scene& sc) {
    std::fstream out(output, std::ios::out | std::ios::binary);
    const char* HEADER = "BVH_00.06";
    out.write(HEADER, strlen(HEADER) + 1);
    out.write((char*)&sc.numTris, sizeof(int));
    out.write((char*)sc.tris, sizeof(triangle) * sc.numTris);
    out.write((char*)&sc.bvh_size, sizeof(int));
    out.write((char*)sc.bvh, sizeof(bvh_node) * sc.bvh_size);
    out.write((char*)&sc.bMin, sizeof(vec3));
    out.write((char*)&sc.bMax, sizeof(vec3));
    out.write((char*)&sc.numLevels, sizeof(int));
    out.close();
}

int main(int argc, char** argv) {
    float scale = 100.0f;
    mat3x3 mat = yUp;
    int numPrimitivesPerLeaf = 1;
    float It = 2.0f;

    if (argc > 1) It = strtof(argv[1], NULL);
    std::cerr << "building bvh with It = " << It << std::endl;

    //TODO include scale in the transformation mat, so we can scale models separately
    std::string basePath = "C:\\Users\\adene\\models\\glsl-assets\\staircase\\";
    std::vector<triangle> tris;
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
    scene s = initScene(tris, scale, false, It);

#ifndef SAH_BVH
    //save("D:\\models\\obj\\cube.bvh", s);
    save("C:\\Users\\adene\\models\\BVH\\staircase.bvh", s);
#endif // !SAH_BVH

    delete[] s.tris;
}
