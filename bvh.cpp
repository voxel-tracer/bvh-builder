#include "bvh.h"

#define __host__
#define __device__
#include "geometry.h"

#include <algorithm>
#include <fstream>

void split(const vec3* triVerts, const Bounds3& triBounds, int planeDim, float planePos, Bounds3& left, Bounds3& right);

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

void BVHAccel::computeQuality(const BVHBuildNode* node, float rootSA, float *largestOverlap) {
    if (node->nPrimitives == 0) {
        BVHBuildNode* left = node->children[0];
        BVHBuildNode* right = node->children[1];
        // compute overlap of children nodes and the ratio of its surface area vs this node's surface area
        float qa = Intersect(left->bounds, right->bounds).SurfaceArea() / rootSA;
        if (qa > *largestOverlap)
            *largestOverlap = qa;

        computeQuality(node->children[0], rootSA, largestOverlap);
        computeQuality(node->children[1], rootSA, largestOverlap);
    }
}

BVHAccel::BVHAccel(const std::vector<std::shared_ptr<Triangle>>& p, int maxPrimsInNode, SplitMethod splitMethod, float internalCost, bool reevaluateCost)
    : maxPrimsInNode(maxPrimsInNode), splitMethod(splitMethod), primitives(p), internalCost(internalCost), reevaluateCost(reevaluateCost) {
    if (primitives.empty()) return;
    // Build BVH from primitives

    // Initialize primitiveInfo array from primitives
    std::vector<BVHPrimitiveInfo> primitiveInfo(primitives.size());
    for (size_t i = 0; i < primitives.size(); ++i)
        primitiveInfo[i] = { i, primitives[i]->bounds };

    // Build BVH tree for primitives using PrimitiveInfo
    int totalNodes = 0;
    int addedSplits = 0;
    std::vector<std::shared_ptr<Triangle>> orderedPrims;
    orderedPrims.reserve(primitives.size());
    BVHBuildNode* root;
    // we only support SplitMethod::EqualCounts
    root = recursiveBuild(primitiveInfo, 0, primitives.size(), &totalNodes, &addedSplits, orderedPrims);
    primitives.swap(orderedPrims);
    primitiveInfo.resize(0);
    // compute estimated node quality
    float largestOverlap = 0.0f;
    computeQuality(root, root->bounds.SurfaceArea(), &largestOverlap);
    std::cerr << "BVH created with " << totalNodes << " nodes for " << (int)primitives.size() << std::endl;
    if (addedSplits > 0) {
        std::cerr << "  Including " << addedSplits << " additional splits" << std::endl;
    }
    std::cerr << "BVH largest overlap is " << largestOverlap << std::endl;
    std::cerr << convertedNodes << " nodes converted to leaves" << std::endl;
    std::cerr << trimmedNodes << " nodes trimmed" << std::endl;

    // Compute representation of depth-first traversal of BVH tree
    nodes.resize(totalNodes);
    int offset = 1;
    flattenBVHTree(root, 0, &offset);
}

BVHAccel::~BVHAccel() { }

struct BucketInfo {
    int count = 0;
    Bounds3 bounds;
};

float BVHAccel::sahCost(const BVHBuildNode* node, float rootSA) const {
    float nodeSARatio = node->bounds.SurfaceArea() / rootSA;
    if (node->nPrimitives > 0) {
        return node->nPrimitives * nodeSARatio;
    }

    return internalCost * nodeSARatio + 
        sahCost(node->children[0], rootSA) + 
        sahCost(node->children[1], rootSA);
}

BVHBuildNode* BVHAccel::recursiveBuild(
    std::vector<BVHPrimitiveInfo>& primitiveInfo,
    int start, int end, int* totalNodes, int* addedSplits,
    std::vector<std::shared_ptr<Triangle>>& orderedPrims) {
    BVHBuildNode* node = new BVHBuildNode();
    (*totalNodes)++;
    // save this value as we may need it later when reevaluating SAH cost
    int firstPrimOffset = orderedPrims.size();
    int curTotalNodes = *totalNodes;
    // compute bounds of all primitives in BVH node
    Bounds3 bounds;
    for (int i = start; i < end; i++)
        bounds = Union(bounds, primitiveInfo[i].bounds);
    int nPrimitives = end - start;
    if (nPrimitives == 1) {
        // create leaf BVHBuildNode
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
            for (int i = start; i < end; i++) {
                int primNum = primitiveInfo[i].primitiveNumber;
                orderedPrims.push_back(primitives[primNum]);
            }
            node->InitLeaf(firstPrimOffset, nPrimitives, bounds);
            return node;
        } else {
            // Partition primitives based on SplitMethod
            switch (splitMethod) {
            case SplitMethod::EqualCounts: {
                // Partition primitives into equal sized subsets
                mid = (start + end) / 2;
                std::nth_element(&primitiveInfo[start], &primitiveInfo[mid], &primitiveInfo[end - 1] + 1,
                    [dim](const BVHPrimitiveInfo& a, const BVHPrimitiveInfo& b) {
                        return a.centroid[dim] < b.centroid[dim];
                    });
                break;
            }
            case SplitMethod::SAH: 
            case SplitMethod::SBVH:
            default: {
                // Partition primitives using approximate SAH
                if (nPrimitives <= 2) {
                    // Partition primitives into equal sized subsets
                    mid = (start + end) / 2;
                    std::nth_element(&primitiveInfo[start], &primitiveInfo[mid], &primitiveInfo[end - 1] + 1,
                        [dim](const BVHPrimitiveInfo& a, const BVHPrimitiveInfo& b) {
                            return a.centroid[dim] < b.centroid[dim];
                        });
                } else {
                    // Allocate BucketInfo for SAH partition buckets
                    constexpr int nBuckets = 12;
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
                    for (int i = 0; i < nBuckets - 1; ++i) {
                        Bounds3 b0, b1;
                        int count0 = 0, count1 = 0;
                        for (int j = 0; j <= i; ++j) {
                            b0 = Union(b0, buckets[j].bounds);
                            count0 += buckets[j].count;
                        }
                        for (int j = i + 1; j < nBuckets; ++j) {
                            b1 = Union(b1, buckets[j].bounds);
                            count1 += buckets[j].count;
                        }
                        cost[i] = internalCost + (count0 * b0.SurfaceArea() + count1 * b1.SurfaceArea()) / bounds.SurfaceArea();
                    }

                    // Find bucket to split that minimizes SAH metric
                    float minCost = cost[0];
                    int minCostSplitBucket = 0;
                    for (int i = 1; i < nBuckets - 1; ++i) {
                        if (cost[i] < minCost) {
                            minCost = cost[i];
                            minCostSplitBucket = i;
                        }
                    }

                    // find best spatial split
                    // TODO only do this if overlap between object splits is big enough
                    bool useSpatialSplit = false;
                    if (splitMethod == SplitMethod::SBVH) {
                        // reset bucket info
                        for (int i = 0; i < nBuckets; i++) {
                            buckets[i].bounds = Bounds3();
                            buckets[i].count = 0;
                        }

                        // Initialize BucketInfo for partition buckets, split primitives if necessary
                        for (int i = start; i < end; ++i) {
                            Bounds3 curBounds = primitiveInfo[i].bounds;

                            // identify all buckets, reference overlaps
                            int firstBin = nBuckets * bounds.Offset(primitiveInfo[i].bounds.pMin)[dim];
                            if (firstBin == nBuckets) firstBin = nBuckets - 1;
                            int lastBin = nBuckets * bounds.Offset(primitiveInfo[i].bounds.pMax)[dim];
                            if (lastBin == nBuckets) lastBin = nBuckets - 1;
                            // update buckets using split references
                            float bucketSize = bounds.Diagonal()[dim] / nBuckets;
                            for (int b = firstBin; b < lastBin; b++) {
                                Bounds3 left, right;
                                split(primitives[primitiveInfo[i].primitiveNumber]->v, curBounds, dim, bounds.pMin[dim] + bucketSize * b, left, right);
                                // left split will be part of current bin
                                buckets[b].count++;
                                buckets[b].bounds = Union(buckets[b].bounds, left);
                                // update current triangle's bounds
                                curBounds = right;
                            }
                            buckets[lastBin].count++;
                            buckets[lastBin].bounds = Union(buckets[lastBin].bounds, curBounds);
                        }

                        // Compute costs for splitting after each bucket
                        for (int i = 0; i < nBuckets - 1; ++i) {
                            Bounds3 b0, b1;
                            int count0 = 0, count1 = 0;
                            for (int j = 0; j <= i; ++j) {
                                b0 = Union(b0, buckets[j].bounds);
                                count0 += buckets[j].count;
                            }
                            for (int j = i + 1; j < nBuckets; ++j) {
                                b1 = Union(b1, buckets[j].bounds);
                                count1 += buckets[j].count;
                            }
                            cost[i] = internalCost + (count0 * b0.SurfaceArea() + count1 * b1.SurfaceArea()) / bounds.SurfaceArea();
                        }

                        // Find bucket to split that minimizes SAH metric
                        // minCost already contains the best object split cost
                        float sbvhMistCost = cost[0];
                        int sbvhMinCostSplitBucket = 0;
                        for (int i = 1; i < nBuckets - 1; ++i) {
                            if (cost[i] < sbvhMistCost) {
                                sbvhMistCost = cost[i];
                                sbvhMinCostSplitBucket = i;
                            }
                        }

                        if (sbvhMistCost < minCost) {
                            minCost = sbvhMistCost;
                            minCostSplitBucket = sbvhMinCostSplitBucket;
                            useSpatialSplit = true;
                        }
                    }

                    // Either create leaf or split primitives at selected SAH bucket
                    float leafCost = nPrimitives;
                    if (nPrimitives <= maxPrimsInNode && leafCost <= minCost) {
                        // Create leaf BVHBuildNode
                        for (int i = start; i < end; i++) {
                            int primNum = primitiveInfo[i].primitiveNumber;
                            orderedPrims.push_back(primitives[primNum]);
                        }
                        node->InitLeaf(firstPrimOffset, nPrimitives, bounds);
                        return node;
                    } else if (useSpatialSplit) {
                        // split all triangles that straddle the split plane
                        // compute split plane
                        float bucketSize = bounds.Diagonal()[dim] / nBuckets;
                        float splitPos = bounds.pMin[dim] + (minCostSplitBucket + 1) * bucketSize;

                        int s = start;
                        int e = end;
                        int splits = 0;
                        while (s < e) {
                            BVHPrimitiveInfo* pInfo = &primitiveInfo[s];
                            if (pInfo->bounds.pMax[dim] < splitPos) {
                                // primitive is on the left of the split plane
                                ++s;
                            } else if (pInfo->bounds.pMin[dim] > splitPos) {
                                // primitive is on the right of the split plane
                                std::swap(primitiveInfo[s], primitiveInfo[--e]);
                            } else {
                                // primitive straddles the splitting plane
                                Bounds3 left, right;
                                split(primitives[pInfo->primitiveNumber]->v, pInfo->bounds, dim, splitPos, left, right);
                                // ignore degenerate cases: splits that generate empty bounds
                                if (left.Diagonal()[left.MaximumExtent()] < 0.00001f) {
                                    //std::cerr << "empty left split s=" << s << " e = " << e << std::endl;
                                    // no need to split the triangle, move it to the right instead
                                    std::swap(primitiveInfo[s], primitiveInfo[--e]);
                                } else if (right.Diagonal()[right.MaximumExtent()] < 0.00001f) {
                                    //std::cerr << "empty right split" << std::endl;
                                    // no need to split the triangle, move it to the left instead
                                    ++s;
                                } else {
                                    //TODO check if we should unsplit the reference instead to improve the SAH cost

                                    // replace primitiveInfo[s] with left split
                                    pInfo->bounds = left;
                                    // insert right split at e
                                    BVHPrimitiveInfo rInfo(pInfo->primitiveNumber, right);
                                    primitiveInfo.insert(primitiveInfo.begin() + e, rInfo);
                                    ++s;
                                    ++splits;
                                }
                            }
                        }

                        // splits contains the number of added splits
                        *addedSplits += splits;
                        mid = e;
                        end += splits;
                        //std::cerr << " added " << splits << " splits" << ". split bucket = " << minCostSplitBucket << std::endl;
                    } else {
                        BVHPrimitiveInfo* pmid = std::partition(&primitiveInfo[start], &primitiveInfo[end - 1] + 1,
                            [=](const BVHPrimitiveInfo& pi) {
                                int b = nBuckets * centroidBounds.Offset(pi.centroid)[dim];
                                if (b == nBuckets) b = nBuckets - 1;
                                return b <= minCostSplitBucket;
                            });
                        mid = pmid - &primitiveInfo[0];
                    }
                }
                break;
            }
            }

            int leftSplits = 0;
            BVHBuildNode* child0 = recursiveBuild(primitiveInfo, start, mid, totalNodes, &leftSplits, orderedPrims);
            int rightSplits = 0;
            BVHBuildNode* child1 = recursiveBuild(primitiveInfo, mid + leftSplits, end + leftSplits, totalNodes, &rightSplits, orderedPrims);
            node->InitInterior(dim, child0, child1);
            *addedSplits += leftSplits + rightSplits;

            // TODO update reevaluateCost to allow resetting split primitives
            if (reevaluateCost && splitMethod == SplitMethod::SAH) {
                // reevaluate node splitting cost and create a leaf instead
                float splitCost = sahCost(node, node->bounds.SurfaceArea());
                float leafCost = nPrimitives;
                if (splitCost > leafCost) {
                    convertedNodes++;
                    trimmedNodes += (*totalNodes) - curTotalNodes;
                    // remove all intermediate nodes added so far
                    (*totalNodes) = curTotalNodes;
                    // primitives have already been written to primitiveInfo and orderedPrims
                    node->InitLeaf(firstPrimOffset, nPrimitives, bounds);
                }
            }
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

void split(const vec3 *triVerts, const Bounds3 &triBounds, int planeDim, float planePos, Bounds3& left, Bounds3& right) {
    int dim = planeDim;
    float pos = planePos;

    vec3 v1 = triVerts[2];
    for (int i = 0; i < 3; i++) {
        vec3 v0 = v1;
        v1 = triVerts[i];
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
    left = Intersect(left, triBounds);
    right = Intersect(right, triBounds);

    // validate that left and right bounds are inside triangle bounds and their union produces triangle bounds
    if (!triBounds.Contains(left))
        std::cerr << "tri.bounds !contains left" << std::endl;
    if (!triBounds.Contains(right))
        std::cerr << "tri.bounds !contains right" << std::endl;
    Bounds3 bounds = Union(left, right);
    if (triBounds.pMin != bounds.pMin || triBounds.pMax != bounds.pMax)
        std::cerr << "tri.bounds != union(left, right)" << std::endl;
}
