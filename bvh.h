#pragma once

#include <memory>
#include <vector>

struct BVHBuildNode;
struct Triangle;
struct LinearBVHNode;
// BVHAccel Forward Declarations
struct BVHPrimitiveInfo;

class BVHAccel {
public:
    // BVHAccel Public Types
    enum class SplitMethod { SAH, SBVH, EqualCounts };

    // BVHAccel Public Methods
    BVHAccel(const std::vector<std::shared_ptr<Triangle>>& p, 
             int maxPrimsInNode = 1, 
             SplitMethod splitMethod = SplitMethod::EqualCounts,
             float internalCost = 0.125f, bool reevaluateCost = false);
    ~BVHAccel();

    // BVHAccel Public Data
    std::vector<std::shared_ptr<Triangle>> primitives;
    std::vector<LinearBVHNode> nodes;
private:
    // BVHAccel Private Methods
    BVHBuildNode* recursiveBuild(
        std::vector<BVHPrimitiveInfo>& primitiveInfo, 
        int start, int end, int* totalNodes, int *addedSplits,
        std::vector<std::shared_ptr<Triangle>>& orderedPrims);

    void flattenBVHTree(BVHBuildNode* node, int offset, int* firstChildOffset);
    void computeQuality(const BVHBuildNode* node, float rootSA, float* largestOverlap);
    float sahCost(const BVHBuildNode* node, float rootSA) const;

    // BVHAccel Private Data
    const int maxPrimsInNode;
    const SplitMethod splitMethod;
    const float internalCost;
    const bool reevaluateCost;
    int convertedNodes = 0;
    int trimmedNodes = 0;
};
