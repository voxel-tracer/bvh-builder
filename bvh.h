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
    enum class SplitMethod { SAH, HLBVH, Middle, EqualCounts }; // start with EqualCounts as it corresponds to our custom builder

    // BVHAccel Public Methods
    BVHAccel(const std::vector<std::shared_ptr<Triangle>>& p, 
             int maxPrimsInNode = 1, 
             SplitMethod splitMethod = SplitMethod::EqualCounts);
    ~BVHAccel();

    // BVHAccel Public Data
    std::vector<std::shared_ptr<Triangle>> primitives;
    std::vector<LinearBVHNode> nodes;
private:
    // BVHAccel Private Methods
    BVHBuildNode* recursiveBuild(
        std::vector<BVHPrimitiveInfo>& primitiveInfo, 
        int start, int end, int* totalNodes, 
        std::vector<std::shared_ptr<Triangle>>& orderedPrims);

    void flattenBVHTree(BVHBuildNode* node, int offset, int* firstChildOffset);

    // BVHAccel Private Data
    const int maxPrimsInNode;
    const SplitMethod splitMethod;
};
