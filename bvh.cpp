#include "bvh.h"

#define __host__
#define __device__
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

void BVHAccel::computeQuality(const BVHBuildNode* node, std::vector<int> &bins, int *total, float binWidth) {
    if (node->nPrimitives == 0) {
        BVHBuildNode* left = node->children[0];
        BVHBuildNode* right = node->children[1];
        // compute overlap of children nodes and the ratio of its surface area vs this node's surface area
        float qa = Intersect(left->bounds, right->bounds).SurfaceArea() / node->bounds.SurfaceArea();
        int bin = floor(qa / binWidth);
        bins[bin]++;
        (*total)++;
        computeQuality(node->children[0], bins, total, binWidth);
        computeQuality(node->children[1], bins, total, binWidth);
    }
}

BVHAccel::BVHAccel(const std::vector<std::shared_ptr<Triangle>>& p, int maxPrimsInNode, SplitMethod splitMethod)
    : maxPrimsInNode(maxPrimsInNode), splitMethod(splitMethod), primitives(p) {
    if (primitives.empty()) return;
    // Build BVH from primitives

    // Initialize primitiveInfo array from primitives
    std::vector<BVHPrimitiveInfo> primitiveInfo(primitives.size());
    for (size_t i = 0; i < primitives.size(); ++i)
        primitiveInfo[i] = { i, primitives[i]->bounds };

    // Build BVH tree for primitives using PrimitiveInfo
    int totalNodes = 0;
    std::vector<std::shared_ptr<Triangle>> orderedPrims;
    orderedPrims.reserve(primitives.size());
    BVHBuildNode* root;
    // we only support SplitMethod::EqualCounts
    root = recursiveBuild(primitiveInfo, 0, primitives.size(), &totalNodes, orderedPrims);
    primitives.swap(orderedPrims);
    primitiveInfo.resize(0);
    // compute estimated node quality
    int numBins = 10;
    int total = 0;
    float binWidth = 1.0f / numBins;
    std::vector<int> bins(numBins);
    computeQuality(root, bins, &total, binWidth);
    std::cerr << "BVH created with " << totalNodes << " nodes for " << (int)primitives.size() << std::endl;
    // compute and display histogram values
    std::cerr << "BVH quality histogram:" << std::endl;
    for (int b = 0; b < numBins; b++) {
        float percent = bins[b] * 100.0f/ total;
        std::cerr << " [" << (b * binWidth) << ", " << ((b + 1) * binWidth) << "]\t = " << percent << std::endl;
    }

    // Compute representation of depth-first traversal of BVH tree
    nodes.resize(totalNodes);
    int offset = 1;
    flattenBVHTree(root, 0, &offset);
}

BVHAccel::~BVHAccel() { }

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
            // Partition primitives based on SplitMethod
            switch (splitMethod) {
            case SplitMethod::EqualCounts:
                // Partition primitives into equal sized subsets
                std::nth_element(&primitiveInfo[start], &primitiveInfo[mid], &primitiveInfo[end - 1] + 1,
                    [dim](const BVHPrimitiveInfo& a, const BVHPrimitiveInfo& b) {
                        return a.centroid[dim] < b.centroid[dim];
                    });
                break;
            case SplitMethod::SAH:
                // Partition primitives using approximate SAH
                if (nPrimitives < 4) {
                    // Partition primitives into equal sized subsets
                    mid = (start + end) / 2;
                    std::nth_element(&primitiveInfo[start], &primitiveInfo[mid], &primitiveInfo[end - 1] + 1,
                        [dim](const BVHPrimitiveInfo& a, const BVHPrimitiveInfo& b) {
                            return a.centroid[dim] < b.centroid[dim];
                        });
                } else {
                    // Allocate BucketInfo for SAH partition buckets
                    constexpr int nBuckets = 12;
                    struct BucketInfo {
                        int count = 0;
                        Bounds3 bounds;
                    };
                    BucketInfo buckets[nBuckets];
                    // Initialize BucketInfo for SAH partition buckets
                    for (int i = start; i < end; ++i) {
                        int b = nBuckets * centroidBounds.Offset(primitiveInfo[i].centroid)[dim];
                        if (b == nBuckets) b = nBuckets - 1;
                        buckets[b].count++;
                        buckets[b].bounds = Union(buckets[b].bounds, primitiveInfo[i].bounds);
                    }
                    // Compute costs for splitting after each bucket
                    float cost[nBuckets - 1];
                    for (int i = 0; i < nBuckets-1; ++i) {
                        Bounds3 b0, b1;
                        int count0 = 0, count1 = 0;
                        for (int j = 0; j <= i; ++j) {
                            b0 = Union(b0, buckets[j].bounds);
                            count0 += buckets[j].count;
                        }
                        for (int j = i+1; j <= nBuckets; ++j) {
                            b1 = Union(b0, buckets[j].bounds);
                            count1 += buckets[j].count;
                        }
                        cost[i] = .125f * (count0 * b0.SurfaceArea() + count1 * b1.SurfaceArea()) / bounds.SurfaceArea();
                    }
                    // Find bucket to split that minimizes SAH metric
                    float minCost = cost[0];
                    int minCostSplitBucket = 0;
                    for (int i = 1; i < nBuckets-1; ++i) {
                        if (cost[i] < minCost) {
                            minCost = cost[i];
                            minCostSplitBucket = i;
                        }
                    }
                    // Either create leaf or split primitives at selected SAH bucket
                    float leafCost = nPrimitives;
                    if (nPrimitives > maxPrimsInNode || minCost < leafCost) {
                        BVHPrimitiveInfo* pmid = std::partition(&primitiveInfo[start], &primitiveInfo[end - 1] + 1,
                            [=](const BVHPrimitiveInfo& pi) {
                                int b = nBuckets * centroidBounds.Offset(pi.centroid)[dim];
                                if (b == nBuckets) b = nBuckets - 1;
                                return b <= minCostSplitBucket;
                            });
                        mid = pmid - &primitiveInfo[0];
                    } else {
                        // Create leaf BVHBuildNode
                        int firstPrimOffset = orderedPrims.size();
                        for (int i = start; i < end; i++) {
                            int primNum = primitiveInfo[i].primitiveNumber;
                            orderedPrims.push_back(primitives[primNum]);
                        }
                        node->InitLeaf(firstPrimOffset, nPrimitives, bounds);
                        return node;
                    }
                }
                break;
            }

            node->InitInterior(dim,
                recursiveBuild(primitiveInfo, start, mid, totalNodes, orderedPrims),
                recursiveBuild(primitiveInfo, mid, end, totalNodes, orderedPrims));
        }
    }
    return node;
}

void BVHAccel::flattenBVHTree(BVHBuildNode* node, int offset, int* firstChildOffset) {
    // store node at offset
    // and its children, if any, at firstChildOffset and firstChildOffset+1
    LinearBVHNode* linearNode = &nodes[offset];
    linearNode->bounds = node->bounds;
    if (node->nPrimitives > 0) {
        linearNode->primitivesOffset = node->firstPrimOffset;
        linearNode->nPrimitives = node->nPrimitives;
    }
    else {
        // Create interior flattened BVH node
        linearNode->axis = node->splitAxis;
        linearNode->nPrimitives = 0;
        linearNode->firstChildOffset = (*firstChildOffset);
        (*firstChildOffset) += 2; // reserve space for both children
        flattenBVHTree(node->children[0], linearNode->firstChildOffset, firstChildOffset);
        flattenBVHTree(node->children[1], linearNode->firstChildOffset+1, firstChildOffset);
    }
}
