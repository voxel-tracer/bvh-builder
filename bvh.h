#pragma once

#include <memory>
#include <vector>
#include <string>

struct BVHBuildNode;
struct Triangle;

// BVHAccel Forward Declarations
struct BVHPrimitiveInfo;
//struct MortonPrimitive; // only needed by HLBVH
struct LinearBVHNode;

class BVHAccel {
public:
    // BVHAccel Public Types
    enum class SplitMethod { SAH, HLBVH, Middle, EqualCounts }; // start with EqualCounts as it corresponds to our custom builder

    // BVHAccel Public Methods
    BVHAccel(const std::string out, const std::vector<std::shared_ptr<Triangle>>& p, 
             int maxPrimsInNode = 1, 
             SplitMethod splitMethod = SplitMethod::EqualCounts);
    ~BVHAccel();

private:
    // BVHAccel Private Methods
    BVHBuildNode* recursiveBuild(
        std::vector<BVHPrimitiveInfo>& primitiveInfo, 
        int start, int end, int* totalNodes, 
        std::vector<std::shared_ptr<Triangle>>& orderedPrims);
    void flattenBVHTree(BVHBuildNode* node, int offset);
    int maxOffset(BVHBuildNode* node, int offset);
    void save(std::string out, int totalNodes);

    // BVHAccel Private Data
    const int maxPrimsInNode;
    const SplitMethod splitMethod;
    std::vector<std::shared_ptr<Triangle>> primitives;
    LinearBVHNode* nodes = nullptr;
};
