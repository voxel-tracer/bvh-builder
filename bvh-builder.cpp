#include <iostream>
#include <time.h>
#include <string>
#include <algorithm>
#include <vector>

#define TINYOBJLOADER_IMPLEMENTATION 
#include "tiny_obj_loader.h"

#define GEOMETRY_IMPLEMENTATION
#define __host__
#define __device__
#include "geometry.h"

#include "bvh.h"
#include "bvh_io.h"


bool loadFromObj(const std::string& filepath, std::vector<std::shared_ptr<Triangle>> &tris, unsigned char meshID) {
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
    vec3 verts[3];
    float tc[6];
    for (auto s = 0; s < shapes.size(); s++) {
        // loop over faces
        size_t index_offset = 0;
        for (auto f = 0; f < shapes[s].mesh.num_face_vertices.size(); f++) {
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

                verts[v] = vec3(vx, vy, vz);
                tc[v * 2 + 0] = tx;
                tc[v * 2 + 1] = ty;
            }
            tris.push_back(std::make_shared<Triangle>(tris.size(), verts[0], verts[1], verts[2], tc, meshID));
            index_offset += 3;
        }
    }

    return true;
}

struct TriSplitInfo {
    float priority;
    Bounds3 bestNode;
    int s;
};

float computeNodeImportance(std::shared_ptr<Triangle> tri, const std::vector<LinearBVHNode>& nodes, TriSplitInfo *info, int triIdx) {
    // go down the BVH tree and stop at first node that the triangle intersects its median
    int idx = 0;
    int level = 1;
    while (true) {
        LinearBVHNode node = nodes[idx];
        if (node.nPrimitives > 0)
            break; // we reached the triangle's leaf node, we can't go further
        if (node.bounds.MedianSplits(tri->bounds))
            break; // we found a node that splits the triangle

        // go to left/right node depending on which side the triangle falls in
        if (nodes[node.firstChildOffset].bounds.Contains(tri->bounds))
            idx = node.firstChildOffset; // go left
        else if (nodes[node.firstChildOffset + 1].bounds.Contains(tri->bounds))
            idx = node.firstChildOffset +1; // go right
        else {
            //LinearBVHNode left = nodes[node.firstChildOffset];
            //LinearBVHNode right = nodes[node.firstChildOffset + 1];
            //std::cerr << "ERROR bestNode doesn't contain triangle " << triIdx << "(" << tri->id << ")" << std::endl;
            break;
        }
        level++; // go down one level
    }

    info->bestNode = nodes[idx].bounds;
    //if (!info->bestNode.Contains(tri->bounds))
    //    std::cerr << "ERROR bestNode doesn't contain its triangle" << std::endl;

    return powf(2, -level);
}

float computeEmptySurfaceArea(std::shared_ptr<Triangle> tri) {
    // compute aabb surface area
    float Aaabb = tri->bounds.SurfaceArea();

    // compute ideal surface area
    vec3 d1 = tri->v[1] - tri->v[0];
    vec3 d2 = tri->v[2] - tri->v[0];
    vec3 cr = cross(d1, d2);
    float Aideal = abs(cr.x()) + abs(cr.y()) + abs(cr.z());

    return (Aaabb - Aideal);
}

void computePriorityAndBestNode(std::shared_ptr<Triangle> tri, const std::vector<LinearBVHNode>& nodes, TriSplitInfo* splitInfo, int triIdx) {
    float importance = computeNodeImportance(tri, nodes, splitInfo, triIdx);
    float emptyArea = computeEmptySurfaceArea(tri);

    float X = 2.0f;
    float Y = 1.0f / 3.0f;
    splitInfo->priority = powf(powf(X, importance) * emptyArea, Y);
}

float sumf(const std::vector<TriSplitInfo> &splitInfos) {
    float s = 0;
    for (int t = 0; t < splitInfos.size(); t++)
        s += splitInfos[t].priority;
    return s;
}

int sumi(const std::vector<TriSplitInfo>& splitInfos, float D) {
    int s = 0;
    for (int t = 0; t < splitInfos.size(); t++)
        s += floor(D * splitInfos[t].priority);
    return s;
}

float computeSplitD(std::vector<TriSplitInfo>& splitInfos, int Sbudget) {
    // start by estimating Dmin, Dmax
    float Dmin = Sbudget / sumf(splitInfos);
    int Smin = sumi(splitInfos, Dmin);
    float Dmax = Dmin * Sbudget / Smin;
    int Smax = sumi(splitInfos, Dmax);
    //std::cerr << "Dmin = " << Dmin << " (" << Smin << "). Dmax = " << Dmax << " (" << Smax << ")" << std::endl;
    // use bissection and iterate 6 times
    while (true) {
        float Davg = (Dmin + Dmax) / 2;
        int Savg = sumi(splitInfos, Davg);
        if (Savg == Smin || Savg == Smax) break;

        if (Savg > Sbudget) {
            Smax = Savg;
            Dmax = Davg;
        }
        else {
            Smin = Savg;
            Dmin = Davg;
        }
        //std::cerr << "Dmin = " << Dmin << " (" << Smin << "). Dmax = " << Dmax << " (" << Smax << ")" << std::endl;
    }

    return Dmin;
}

int computeNumSplits(
        const std::vector<std::shared_ptr<Triangle>> &primitives,
        const std::vector<LinearBVHNode> &nodes, 
        std::vector<TriSplitInfo> &splitInfos, 
        int Sbudget) {
    // compute and store split priority for each triangle
    for (int t = 0; t < primitives.size(); t++)
        computePriorityAndBestNode(primitives[t], nodes, &splitInfos[t], t);
    // compute split D
    float D = computeSplitD(splitInfos, Sbudget);
    // compute num splits for each triangle
    int numSplits = 0;
    for (int t = 0; t < primitives.size(); t++) {
        splitInfos[t].s = floor(D * splitInfos[t].priority);
        numSplits += splitInfos[t].s;
    }

    return numSplits;
}

Bounds3 findBestNode(Bounds3 node, std::shared_ptr<Triangle> tri) {
    Bounds3 cur = node;
    while (!cur.MedianSplits(tri->bounds)) {
        Bounds3 left = cur.leftSplit();
        if (left.Contains(tri->bounds))
            cur = left;
        else
            cur = cur.rightSplit();
    }
    return cur;
}

float clamp(float v, float vmin, float vmax) {
    return fminf(vmax, fmaxf(vmin, v));
}

void split(std::shared_ptr<Triangle> tri, Bounds3 bestNode, Bounds3 &left, Bounds3 &right) {
    int dim = bestNode.MaximumExtent();
    float pos = bestNode.Center()[dim];

    vec3 v1 = tri->v[2];
    for (int i = 0; i < 3; i++) {
        vec3 v0 = v1;
        v1 = tri->v[i];
        float v1p = v1[dim];
        float v0p = v0[dim];

        if (v0p <= pos) left = Union(left, v0);
        if (v0p >= pos) right = Union(right, v0);
        // check if edge intersects the plane and store it in both sides
        if ((v0p < pos && v1p > pos) || (v0p > pos && v1p < pos)) {
            vec3 t = lerp(v0, v1, clamp((pos - v0p) / (v1p - v0p), 0.0f, 1.0f));
            left = Union(left, t);
            right = Union(right, t);
        }
    }
    // intersect original bounds in case tri is a split already
    left.pMax[dim] = pos;
    right.pMin[dim] = pos;
    left = Intersect(left, tri->bounds);
    right = Intersect(right, tri->bounds);

    // validate that left and right bounds are inside triangle bounds and their union produces triangle bounds
    if (!tri->bounds.Contains(left))
        std::cerr << "tri.bounds !contains left" << std::endl;
    if (!tri->bounds.Contains(right))
        std::cerr << "tri.bounds !contains right" << std::endl;
    Bounds3 bounds = Union(left, right);
    if (tri->bounds.pMin != bounds.pMin || tri->bounds.pMax != bounds.pMax)
        std::cerr << "tri.bounds != union(left, right)" << std::endl;
}

void updateSplitCount(int s, Bounds3 left, Bounds3 right, int& sleft, int& sright) {
    float wa = left.Diagonal()[left.MaximumExtent()];
    float wb = right.Diagonal()[right.MaximumExtent()];

    sleft = floor((s - 1) * wa / (wa + wb) + 0.5f);
    sright = s - 1 - sleft;
}

BVHAccel* split(BVHAccel* accel, float Sexcess) {
    // start splitting the triangles
    int Sbudget = floor(accel->primitives.size() * Sexcess);
    std::cerr << "start splitting triangles. We will add at most " << Sbudget << " extra splits" << std::endl;

    // compute for each triangle its num splits and most important node idx
    std::vector<TriSplitInfo> splitInfos(accel->primitives.size());
    int numSplits = computeNumSplits(accel->primitives, accel->nodes, splitInfos, Sbudget);
    std::cerr << "adding an additional " << numSplits << " splits" << std::endl;

    // allocate enough triangles to store all triangles including the splits
    int total = accel->primitives.size() + numSplits;
    std::vector<std::shared_ptr<Triangle>> primitives(accel->primitives);

    // actually split all triangles according to their numSplit (.s)
    int nextSplit = accel->primitives.size();

    primitives.resize(total);
    splitInfos.resize(total);
    //int skipSplits = 0;
    //int leftToSplit = numSplits;
    //bool doneSplitting = false;
    for (int i = 0; /*!doneSplitting &&*/ i < total; i++) {
        std::shared_ptr<Triangle> tri = primitives[i];
        TriSplitInfo &info = splitInfos[i];
        // keep splitting this triangle until we run out of splits
        while (/*!doneSplitting &&*/ info.s > 0) {
            //if (--skipSplits > 0) {
            //    tri.s = 0;
            //    continue;
            //}
            //std::cerr << "splitting triangle " << i << " with " << info.s << " splits" << std::endl;
            // find most important node that splits this triangle, starting from its bestIdx
            Bounds3 node = findBestNode(info.bestNode, tri);
            Bounds3 left, right;
            split(tri, node, left, right);
            // validate splits
            if (!left.IsValid())
                std::cerr << "left split invalid" << std::endl;
            if (!right.IsValid())
                std::cerr << "right split invalid" << std::endl;
            if (left.Diagonal()[left.MaximumExtent()] < 0.00001f)
                std::cerr << " left bound empty" << std::endl;
            if (right.Diagonal()[right.MaximumExtent()] < 0.00001f)
                std::cerr << " right bound empty" << std::endl;

            // compute split count for left and right splits
            int sleft, sright;
            updateSplitCount(info.s, left, right, sleft, sright);

            // replace current triangle with left split
            tri->bounds = left;
            info.s = sleft;

            // add new split to the end of the array
            primitives[nextSplit] = std::make_shared<Triangle>(*tri, right);
            splitInfos[nextSplit] = { info.priority, info.bestNode, sright };
            nextSplit++;

            //if (--leftToSplit == 0)
            //    doneSplitting = true;
        }
    }

    delete accel;
    return new BVHAccel(primitives);
}

int main(int argc, char** argv) {
    std::string basePath = "C:\\Users\\adene\\models\\glsl-assets\\staircase\\";
    std::vector<std::shared_ptr<Triangle>> tris;
    bool success = true;
    bool splitTriangles = true;

    int Sexcess = 0;
    if (argc > 1) {
        Sexcess = strtol(argv[1], NULL, 10);
    }

    success = success && loadFromObj(basePath + "Black.obj", tris, 0);
    success = success && loadFromObj(basePath + "Brass.obj", tris, 1);
    success = success && loadFromObj(basePath + "BrushedAluminium.obj", tris, 2);
    success = success && loadFromObj(basePath + "Candles.obj", tris, 3);
    success = success && loadFromObj(basePath + "ChairSeat.obj", tris, 4);
    success = success && loadFromObj(basePath + "Glass.obj", tris, 5);
    success = success && loadFromObj(basePath + "Gold.obj", tris, 6);
    success = success && loadFromObj(basePath + "Lampshade.obj", tris, 7);
    success = success && loadFromObj(basePath + "MagnoliaPaint.obj", tris, 8);
    success = success && loadFromObj(basePath + "Painting1.obj", tris, 9);
    success = success && loadFromObj(basePath + "Painting2.obj", tris, 10);
    success = success && loadFromObj(basePath + "Painting3.obj", tris, 11);
    success = success && loadFromObj(basePath + "StainlessSteel.obj", tris, 12);
    success = success && loadFromObj(basePath + "Wallpaper.obj", tris, 13);
    success = success && loadFromObj(basePath + "WhitePaint.obj", tris, 14);
    success = success && loadFromObj(basePath + "WhitePlastic.obj", tris, 15);
    success = success && loadFromObj(basePath + "WoodChair.obj", tris, 16);
    success = success && loadFromObj(basePath + "WoodFloor.obj", tris, 17);
    success = success && loadFromObj(basePath + "WoodLamp.obj", tris, 18);
    success = success && loadFromObj(basePath + "WoodStairs.obj", tris, 19);
    
    if (!success) {
        std::cerr << "Failed to load obj files" << std::endl;
        return -1;
    }
    std::cerr << "read " << tris.size() << " triangles" << std::endl;
    time_t start = clock();
    BVHAccel *accel = new BVHAccel(tris);
    std::cerr << "BVH build took " << ((float)(clock() - start)) / CLOCKS_PER_SEC << " seconds" << std::endl;

    if (splitTriangles && Sexcess > 0) {
        start = clock();
        accel = split(accel, Sexcess / 100.0f);
        std::cerr << "Splitting and rebuilding the BVH took " << ((float)(clock() - start)) / CLOCKS_PER_SEC << " seconds" << std::endl;
    }

    saveBVH("C:\\Users\\adene\\models\\BVH\\staircase.bvh", accel->primitives, accel->nodes, 1);
    delete accel;
}
