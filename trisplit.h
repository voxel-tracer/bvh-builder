#pragma once

#include <vector>

#include "geometry.h"


void split(
    const std::vector<std::shared_ptr<Triangle>>& primitives,
    const std::vector<LinearBVHNode>& nodes, float Sexcess,
    std::vector<std::shared_ptr<Triangle>>& splitPrimitives);