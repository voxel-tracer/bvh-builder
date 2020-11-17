#include "bvh.h"
#include "geometry.h"

#include <algorithm>
#include <fstream>

// BVHAccel Local Declarations
struct BVHPrimitiveInfo {
    BVHPrimitiveInfo() : primitiveNumber(0) {}
    BVHPrimitiveInfo(size_t primitiveNumber, const Bounds3& bounds)
        : primitiveNumber(primitiveNumber), bounds(bounds), centroid(.5f * bounds.pMin + .5f * bounds.pMax) { }
    size_t primitiveNumber;
    Bounds3 bounds;
    vec3 centroid;
};

struct BVHBuildNode {
    void InitLeaf(int first, int n, const Bounds3& b) {
        firstPrimOffset = first;
        nPrimitives = n;
        bounds = b;
        children[0] = children[1] = nullptr;
    }

    void InitInterior(int axis, BVHBuildNode* c0, BVHBuildNode* c1) {
        children[0] = c0;
        children[1] = c1;
        bounds = Union(c0->bounds, c1->bounds);
        splitAxis = axis;
        nPrimitives = 0;
    }

    Bounds3 bounds;
    BVHBuildNode* children[2];
    uint32_t splitAxis, firstPrimOffset, nPrimitives;
};

struct LinearBVHNode {
    Bounds3 bounds;
    union {
        int primitivesOffset;   // leaf
        int secondChildOffset;  // interior
    };
    uint16_t nPrimitives;   // 0 -> interior node
    uint8_t axis;           // interior node: xyz
    uint8_t pad[1];         // ensure 32 bytes total size
};

// file structs compatible with my current renderer, may remove this in the future
struct FileTriangle {
    FileTriangle() {}
    FileTriangle(const Triangle& t) {
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

BVHAccel::BVHAccel(const std::string out, const std::vector<std::shared_ptr<Triangle>>& p, int maxPrimsInNode, SplitMethod splitMethod)
    : maxPrimsInNode(maxPrimsInNode), splitMethod(splitMethod), primitives(p) {
    if (primitives.size() == 0)
        return;
    // Build BVH from primitives
    // Initialize primitiveInfo array from primitives
    std::vector<BVHPrimitiveInfo> primitiveInfo(primitives.size());
    for (size_t i = 0; i < primitives.size(); ++i)
        primitiveInfo[i] = { i, primitives[i]->bounds };
    // Build BVH tree for primitives using PrimitiveInfo
    int totalNodes = 0;
    std::vector<std::shared_ptr<Triangle>> orderedPrims;
    BVHBuildNode* root;
    // we only support SplitMethod::EqualCounts
    root = recursiveBuild(primitiveInfo, 0, primitives.size(), &totalNodes, orderedPrims);
    primitives.swap(orderedPrims);
    // Compute representation of depth-first traversal of BVH tree
    // node at index idx has its children stored in idx*2, idx*2+1. This allows the traversal to deduce indices without explicitely storing them
    totalNodes = maxOffset(root, 1);
    nodes = new LinearBVHNode[totalNodes];
    flattenBVHTree(root, 1);
    save(out, totalNodes);
}

BVHAccel::~BVHAccel() { delete[] nodes; }

BVHBuildNode* BVHAccel::recursiveBuild(
    std::vector<BVHPrimitiveInfo>& primitiveInfo,
    int start, int end, int* totalNodes,
    std::vector<std::shared_ptr<Triangle>>& orderedPrims) {
    BVHBuildNode* node = new BVHBuildNode();
    (*totalNodes)++;
    // compute bounds of all primitives in BVH node
    Bounds3 bounds;
    for (int i = start; i < end; i++)
        bounds = Union(bounds, primitiveInfo[i].bounds);
    int nPrimitives = end - start;
    if (nPrimitives == 1) {
        // create leaf BVHBuildNode
        int firstPrimOffset = orderedPrims.size();
        for (int i = start; i < end; i++) {
            int primNum = primitiveInfo[i].primitiveNumber;
            orderedPrims.push_back(primitives[primNum]);
        }
        node->InitLeaf(firstPrimOffset, nPrimitives, bounds);
        return node;
    } else {
        // Compute bound of primitive centroids, chose split dimension dim
        Bounds3 centroidBounds;
        for (int i = start; i < end; i++)
            centroidBounds = Union(centroidBounds, primitiveInfo[i].centroid);
        int dim = centroidBounds.MaximumExtent();
        // Partition primitives into two sets and build children
        int mid = (start + end) / 2;
        if (centroidBounds.pMax[dim] == centroidBounds.pMin[dim]) {
            // create leaf BVHBuildNode
            int firstPrimOffset = orderedPrims.size();
            for (int i = start; i < end; i++) {
                int primNum = primitiveInfo[i].primitiveNumber;
                orderedPrims.push_back(primitives[primNum]);
            }
            node->InitLeaf(firstPrimOffset, nPrimitives, bounds);
            return node;
        } else {
            // Partition primitives based on SplitMethod, we only support EqualCounts for now
            // Partition primitives into equal sized subsets
            std::nth_element(&primitiveInfo[start], &primitiveInfo[mid], &primitiveInfo[end - 1] + 1,
                [dim](const BVHPrimitiveInfo& a, const BVHPrimitiveInfo& b) {
                    return a.centroid[dim] < b.centroid[dim];
                });
            node->InitInterior(dim,
                recursiveBuild(primitiveInfo, start, mid, totalNodes, orderedPrims),
                recursiveBuild(primitiveInfo, mid, end, totalNodes, orderedPrims));
        }
    }
    return node;
}

int BVHAccel::maxOffset(BVHBuildNode* node, int offset) {
    if (node->nPrimitives > 0) {
        return offset; // Leaf node
    } else {
        return std::max(
            maxOffset(node->children[0], offset * 2),
            maxOffset(node->children[1], offset * 2 + 1));
    }
}

void BVHAccel::flattenBVHTree(BVHBuildNode* node, int offset) {
    LinearBVHNode* linearNode = &nodes[offset];
    linearNode->bounds = node->bounds;
    if (node->nPrimitives > 0) {
        linearNode->primitivesOffset = node->firstPrimOffset;
        linearNode->nPrimitives = node->nPrimitives;
    } else {
        // Create interior flattened BVH node
        linearNode->axis = node->splitAxis;
        linearNode->nPrimitives = 0;
        flattenBVHTree(node->children[0], offset*2);
        flattenBVHTree(node->children[1], offset * 2 + 1);
    }
}

// 0.05 LinearBVHNode instead of Bounds3
void BVHAccel::save(std::string output, int totalNodes) {
    std::fstream out(output, std::ios::out | std::ios::binary);
    // start with header
    const char* HEADER = "BVH_00.05";
    out.write(HEADER, strlen(HEADER) + 1);

    // convert Triangle to triangle and write them to disk
    int numPrimitives = primitives.size();
    FileTriangle* tris = new FileTriangle[numPrimitives];
    for (int i = 0; i < numPrimitives; i++)
        tris[i] = FileTriangle(*(primitives[i]));
    out.write((char*)&numPrimitives, sizeof(int));
    out.write((char*)tris, sizeof(FileTriangle) * numPrimitives);
    delete[] tris;

    // save linear nodes
    out.write((char*)&totalNodes, sizeof(int));
    out.write((char*)nodes, sizeof(LinearBVHNode) * totalNodes);
    // write scene bounds
    out.write((char*)&nodes[1].bounds.pMin, sizeof(vec3));
    out.write((char*)&nodes[1].bounds.pMax, sizeof(vec3));
    out.write((char*)&maxPrimsInNode, sizeof(int));
    out.close();

}