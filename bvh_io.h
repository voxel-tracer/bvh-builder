#pragma once

#include <vector>
#include <memory>
#include <fstream>

#include "geometry.h"

// 0.05 LinearBVHNode instead of Bounds3, no scene Bounds3
// 0.06 children nodes are stored in successive offsets

void saveBVH(std::string output, 
             std::vector<std::shared_ptr<Triangle>> &primitives,
             std::vector<LinearBVHNode> &nodes,
             int maxPrimsInNode) {
    std::fstream out(output, std::ios::out | std::ios::binary);
    // start with header
    const char* HEADER = "BVH_00.06";
    out.write(HEADER, strlen(HEADER) + 1);

    // convert Triangle to triangle and write them to disk
    int numPrimitives = primitives.size();
    LinearTriangle* tris = new LinearTriangle[numPrimitives];
    for (int i = 0; i < numPrimitives; i++)
        tris[i] = LinearTriangle(*(primitives[i]));
    out.write((char*)&numPrimitives, sizeof(int));
    out.write((char*)tris, sizeof(LinearTriangle) * numPrimitives);
    delete[] tris;

    // save linear nodes
    int totalNodes = nodes.size();
    LinearBVHNode* nodesArray = new LinearBVHNode[totalNodes];
    for (int i = 0; i < totalNodes; i++)
        nodesArray[i] = nodes[i];
    out.write((char*)&totalNodes, sizeof(int));
    out.write((char*)nodesArray, sizeof(LinearBVHNode) * totalNodes);
    delete[] nodesArray;

    // write scene bounds
    out.write((char*)&maxPrimsInNode, sizeof(int));
    out.close();
}

bool loadBVH(std::string input, 
             std::vector<LinearTriangle> &primitives,
             std::vector<LinearBVHNode> &nodes,
             int *maxPrimsInNode) {
    std::fstream in(input, std::ios::in | std::ios::binary);

    const char* HEADER = "BVH_00.06";
    int headerLen = strlen(HEADER) + 1;
    char* header = new char[headerLen];
    in.read(header, headerLen);
    if (strcmp(HEADER, header) != 0) {
        std::cerr << "invalid header " << header << std::endl;
        return false;
    }

    // load triangles
    int numTriangles;
    in.read((char*)&numTriangles, sizeof(int));
    LinearTriangle* tris = new LinearTriangle[numTriangles];
    in.read((char*)tris, sizeof(LinearTriangle) * numTriangles);
    primitives.reserve(numTriangles);
    for (int i = 0; i < numTriangles; i++)
        primitives.push_back(tris[i]);
    delete[] tris;

    // load BVH nodes
    int numNodes;
    in.read((char*)&numNodes, sizeof(int));
    LinearBVHNode* rNodes = new LinearBVHNode[numNodes];
    in.read((char*)rNodes, sizeof(LinearBVHNode) * numNodes);
    nodes.reserve(numNodes);
    for (int i = 0; i < numNodes; i++)
        nodes.push_back(rNodes[i]);
    delete[] rNodes;

    in.read((char*)maxPrimsInNode, sizeof(int));

    return true;
}
