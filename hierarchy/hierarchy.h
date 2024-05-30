#pragma once

#include <vector>

#define GLM_ENABLE_EXPERIMENTAL

#include <glm/glm.hpp>
#include <glm/gtx/quaternion.hpp>

#include "aabb.h"
#include "geometry/geometry.h"

#define SPLIT_COUNT 4

struct Node {
    AABB aabb;

    int left_child = -1;
    int right_child = -1;
    int first_primitive_id = -1;
    int primitive_count = -1;
};

class BVH {
private:
    std::vector<std::unique_ptr<Primitive>>* _primitives;

    std::vector<Node> _nodes;
    int _root;

public:
    BVH();
    BVH(std::vector<std::unique_ptr<Primitive>>* primitives);

    int add_node(int first_primitive_id, int primitive_count);

    std::pair<float, int> get_t_for(int node_id, glm::vec3 O, glm::vec3 D,
                                    float border_t, int border_primitive_id) const;
    std::pair<float, int> get_t(glm::vec3 O, glm::vec3 D) const;
};
