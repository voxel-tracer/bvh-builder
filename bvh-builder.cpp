#include <iostream>
#include <time.h>
#include <string>
#include <algorithm>
#include <vector>

#define __host__
#define __device__
#define GEOMETRY_IMPLEMENTATION
#include "geometry.h"

#include "trisplit.h"

#define TINYOBJLOADER_IMPLEMENTATION 
#include "tiny_obj_loader.h"


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

BVHAccel* split(BVHAccel* accel, BVHAccel::SplitMethod splitMethod, float internalCost, float Sexcess) {
    std::vector<std::shared_ptr<Triangle>> splitPrimitives;
    split(accel->primitives, accel->nodes, Sexcess, splitPrimitives);

    delete accel;
    return new BVHAccel(splitPrimitives, 1, splitMethod, internalCost);
}

int main(int argc, char** argv) {
    std::string basePath = "C:\\Users\\adene\\models\\glsl-assets\\staircase\\";
    std::vector<std::shared_ptr<Triangle>> tris;
    bool success = true;
    bool splitTriangles = true;
    BVHAccel::SplitMethod splitMethod = BVHAccel::SplitMethod::SAH;
    float internalCost = 1.0f;

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
    BVHAccel *accel = new BVHAccel(tris, 1, splitMethod, internalCost);
    std::cerr << "BVH build took " << ((float)(clock() - start)) / CLOCKS_PER_SEC << " seconds" << std::endl;

    if (splitTriangles && Sexcess > 0) {
        start = clock();
        accel = split(accel, splitMethod, internalCost, Sexcess / 100.0f);
        std::cerr << "Splitting and rebuilding the BVH took " << ((float)(clock() - start)) / CLOCKS_PER_SEC << " seconds" << std::endl;
    }

    saveBVH("C:\\Users\\adene\\models\\BVH\\staircase.bvh", accel->primitives, accel->nodes, 1);
    delete accel;
}
