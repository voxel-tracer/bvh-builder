#define __host__
#define __device__
#include "trisplit.h"

struct TriSplitInfo {
    float priority;
    Bounds3 bestNode;
    int s;
};

struct SplittingNode {
    int nodeIdx = -1;
    int level;
    int dim; // splitting dimension
    float pos; // splitting position
};

struct TriangleBounds {
    int idx;
    Bounds3 bounds;
};

bool findBestSplittingNode(const std::vector<LinearBVHNode>& nodes, int nodeIdx, int nodeLvl, const TriangleBounds &tri, SplittingNode *best) {
    const LinearBVHNode *node = &nodes[nodeIdx];

    if (!node->bounds.Contains(tri.bounds))
        return false; // node cannot be an ancestor of triangle

    bool found = false;
    if (node->nPrimitives > 0) {
        // TODO confirm that tri.idx is contained in [node->primitiveOffset, node->primitiveOffset + node->nPrimitives[
        // even though the triangle is completely inside the leaf node, it may not be its parent node
        if (tri.idx >= node->primitivesOffset && tri.idx < (node->primitivesOffset + node->nPrimitives)) {
            // it is this triangle's BVH node
            found = true;
        } else {
            return false;
        }
    } else {
        // for internal nodes, check if the triangle is part of the subtree
        found = findBestSplittingNode(nodes, node->firstChildOffset, nodeLvl + 1, tri, best) ||
            findBestSplittingNode(nodes, node->firstChildOffset + 1, nodeLvl + 1, tri, best);
    }

    if (found && node->bounds.MedianSplits(tri.bounds)) {
        // this node is the best splitting node so far (when going bottom up)
        best->nodeIdx = nodeIdx;
        best->level = nodeLvl;
        best->dim = node->bounds.MaximumExtent();
        best->pos = node->bounds.Center()[best->dim];
    } else if (node->nPrimitives > 0) {
        // this shouldn't happen, a leaf node should always split its triangle (unless it contains multiple triangles)
        std::cerr << "non splitting leaf node found!" << std::endl;
    }

    return found;
}

float computeNodeImportance(std::shared_ptr<Triangle> tri, const std::vector<LinearBVHNode>& nodes, TriSplitInfo* info, int triIdx) {
    SplittingNode best;
    TriangleBounds bounds = { triIdx, tri->bounds };
    bool found = findBestSplittingNode(nodes, 0, 1, bounds, &best);
    // found should always be true
    //if (!found) {
    //    std::cerr << "ERROR couldn't find triangle " << triIdx << " (" << tri->id << ")" << std::endl;
    //} else if (best.nodeIdx == -1) {
    //    std::cerr << "ERROR couldn't find best splitting node for triangle " << triIdx << " (" << tri->id << ")" << std::endl;
    //}

    info->bestNode = nodes[best.nodeIdx].bounds;
    //if (!info->bestNode.Contains(tri->bounds))
    //    std::cerr << "ERROR bestNode doesn't contain its triangle" << std::endl;

    return powf(2, -best.level);
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

float sumf(const std::vector<TriSplitInfo>& splitInfos) {
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
    const std::vector<std::shared_ptr<Triangle>>& primitives,
    const std::vector<LinearBVHNode>& nodes,
    std::vector<TriSplitInfo>& splitInfos,
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

void split(std::shared_ptr<Triangle> tri, Bounds3 bestNode, Bounds3& left, Bounds3& right) {
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

void split(const std::vector<std::shared_ptr<Triangle>>& primitives, 
    const std::vector<LinearBVHNode>& nodes, float Sexcess,
    std::vector<std::shared_ptr<Triangle>>& splitPrimitives) {

    // start splitting the triangles
    int Sbudget = floor(primitives.size() * Sexcess);
    std::cerr << "start splitting triangles. We will add at most " << Sbudget << " extra splits" << std::endl;

    // compute for each triangle its num splits and most important node idx
    std::vector<TriSplitInfo> splitInfos(primitives.size());
    int numSplits = computeNumSplits(primitives, nodes, splitInfos, Sbudget);
    std::cerr << "adding an additional " << numSplits << " splits" << std::endl;

    // allocate enough triangles to store all triangles including the splits
    int total = primitives.size() + numSplits;

    // actually split all triangles according to their numSplit (.s)
    int nextSplit = primitives.size();

    // copy all primitives to splitPrimitives
    splitPrimitives.reserve(total);
    std::copy(primitives.begin(), primitives.end(), std::back_inserter(splitPrimitives));
    splitPrimitives.resize(total);
    splitInfos.resize(total);

    //int skipSplits = 0;
    //int leftToSplit = numSplits;
    //bool doneSplitting = false;
    for (int i = 0; /*!doneSplitting &&*/ i < total; i++) {
        std::shared_ptr<Triangle> tri = splitPrimitives[i];
        TriSplitInfo& info = splitInfos[i];
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
            //if (left.Diagonal()[left.MaximumExtent()] < 0.00001f)
            //    std::cerr << " left bound empty idx= " << i << ", s= " << info.s << std::endl;
            //if (right.Diagonal()[right.MaximumExtent()] < 0.00001f)
            //    std::cerr << " right bound empty idx= " << i << ", s= " << info.s << std::endl;
            if (left.Diagonal()[left.MaximumExtent()] < 0.00001f ||
                right.Diagonal()[right.MaximumExtent()] < 0.00001f) {
                // left or right split has empty bounds, 
                // it means that splitting current tri was a no op and result will just duplicate it
                // so we should just stop splitting this tri
                //std::cerr << "no op split for tri " << i << " and s = " << info.s << std::endl;
                total -= info.s;

                info.s = 0;
            }
            else {
                // compute split count for left and right splits
                int sleft, sright;
                updateSplitCount(info.s, left, right, sleft, sright);

                // replace current triangle with left split
                tri->bounds = left;
                info.s = sleft;

                // add new split to the end of the array
                splitPrimitives[nextSplit] = std::make_shared<Triangle>(*tri, right);
                splitInfos[nextSplit] = { info.priority, info.bestNode, sright };
                nextSplit++;

                //if (--leftToSplit == 0)
                //    doneSplitting = true;
            }
        }

        splitPrimitives.resize(total);
    }
}