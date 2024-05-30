#include <algorithm>

#include "hierarchy.h"
#include "aabb.h"

glm::vec3 another_fast_rotate(glm::quat quat, glm::vec3 vec) {
    glm::vec3 cut_quat = {quat.x, quat.y, quat.z};
    glm::vec3 t = 2.f * glm::cross(cut_quat, vec);
    return vec + quat.w * t + glm::cross(cut_quat, t);
}

void AABB::make_infinite() {
    _is_infinite = true;
}

void AABB::replace_with(AABB aabb) {
    _min = aabb._min;
    _max = aabb._max;
    _is_infinite = aabb._is_infinite;
}

void AABB::extend(glm::vec3 point) {
    _min = glm::min(_min, point);
    _max = glm::max(_max, point);
}

void AABB::extend(AABB aabb) {
    _min = glm::min(_min, aabb._min);
    _max = glm::max(_max, aabb._max);
    _is_infinite |= aabb._is_infinite;
}

void AABB::move(glm::vec3 vector) {
    _min += vector;
    _max += vector;
}

void AABB::rotate(glm::quat quat) {
    if (_is_infinite)
        return;

    glm::vec3 choices[2];
    choices[0] = _min;
    choices[1] = _max;

    AABB new_aabb;

    for (auto &chosen_x : choices)
        for (auto &chosen_y : choices)
            for (auto &chosen_z : choices)
                new_aabb.extend(another_fast_rotate(quat, glm::vec3(chosen_x.x, chosen_y.y, chosen_z.z)));

    replace_with(new_aabb);
}

std::tuple<int, float> AABB::get_cut() {
    glm::vec3 diff = _max - _min, mid = (_max + _min) / 2.f;

    if (diff.x >= diff.y and diff.x >= diff.z)
        return {0, mid.x};

    if (diff.y > diff.z)
        return {1, mid.x};

    return {2, mid.z};
}

bool AABB::divide(int cut_axis, float cut_coord) {
    if (_is_infinite)
        return false;

    glm::vec3 center = (_max + _min) / 2.f;

    if (cut_axis == 0)
        return center.x < cut_coord;
    if (cut_axis == 1)
        return center.y < cut_coord;
    return center.z < cut_coord;
}

bool AABB::is_intersected(glm::vec3 O, glm::vec3 D, float border_t) const {
    if (_is_infinite)
        return true;

    glm::vec3 position = (_max + _min) / 2.f;
    glm::vec3 box_s = glm::abs(_max - _min) / 2.f;

    float t1, t2;
    std::tie(t1, t2) = Box(box_s).get_ts(O - position, glm::normalize(D));

    return (t1 < border_t || border_t < 0) && t2 > 0;
}

BVH::BVH() {
    _primitives = nullptr;
    _root = -1;
}

BVH::BVH(std::vector<std::unique_ptr<Primitive>>* primitives) {
    _primitives = primitives;
    _root = add_node(0, (int) primitives->size());
}

int BVH::add_node(int first_primitive_id, int primitive_count) {
    int id = (int) _nodes.size();
    _nodes.emplace_back();

    for (int primitive_id = first_primitive_id; primitive_id < first_primitive_id + primitive_count; ++primitive_id) {
        auto &primitive = (*_primitives)[primitive_id];
        auto p_aabb = primitive->get_aabb();

        AABB aabb = _nodes[id].aabb;
        aabb.extend(p_aabb);
        _nodes[id].aabb = aabb;
    }

    if (primitive_count <= SPLIT_COUNT) {
        _nodes[id].first_primitive_id = first_primitive_id;
        _nodes[id].primitive_count = primitive_count;
        return id;
    }

    int cut_axis;
    float cut_coord;
    std::tie(cut_axis, cut_coord) = _nodes[id].aabb.get_cut();

    auto partition_function = [&cut_axis, &cut_coord](std::unique_ptr<Primitive> &primitive){
        return primitive->get_aabb().divide(cut_axis, cut_coord);
    };

    auto begin = _primitives->begin() + first_primitive_id;
    auto end = _primitives->begin() + first_primitive_id + primitive_count;
    auto part = std::partition(begin, end, partition_function);

    begin = _primitives->begin() + first_primitive_id;
    end = _primitives->begin() + first_primitive_id + primitive_count;

    int left_primitive_count = (int) (part - begin);
    int right_first_primitive_id = (int) (part - _primitives->begin());
    int right_primitive_count = (int) (end - part);

    if (begin == part || part == end) {
        _nodes[id].first_primitive_id = first_primitive_id;
        _nodes[id].primitive_count = primitive_count;
        return id;
    }

    int left_child = add_node(first_primitive_id, left_primitive_count);
    _nodes[id].left_child = left_child;

    int right_child = add_node(right_first_primitive_id, right_primitive_count);
    _nodes[id].right_child = right_child;

    return id;
}

std::pair<float, int> BVH::get_t_for(int node_id, glm::vec3 O, glm::vec3 D,
                                     float border_t, int border_primitive_id) const {
    const Node &node = _nodes[node_id];
    D = glm::normalize(D);

    if (not node.aabb.is_intersected(O, D, border_t))
        return {border_t, border_primitive_id};

    if (node.right_child != -1) {
        float tmp_t = border_t;
        int tmp_primitive_id = border_primitive_id;

        std::tie(tmp_t, tmp_primitive_id) = get_t_for(node.right_child, O, D, border_t, border_primitive_id);

        if (tmp_t > 0 && (tmp_t < border_t || border_t < 0)) {
            border_t = tmp_t;
            border_primitive_id = tmp_primitive_id;
        }

        std::tie(tmp_t, tmp_primitive_id) = get_t_for(node.left_child, O, D, border_t, border_primitive_id);

        if (tmp_t > 0 && (tmp_t < border_t || border_t < 0)) {
            border_t = tmp_t;
            border_primitive_id = tmp_primitive_id;
        }
    }

    else {
        for (int primitive_id = node.first_primitive_id;
             primitive_id < node.first_primitive_id + node.primitive_count;
             ++primitive_id) {
            float tmp_t = (*_primitives)[primitive_id]->get_t(O, D);

            if (tmp_t > 0 && (tmp_t < border_t || border_t < 0)) {
                border_t = tmp_t;
                border_primitive_id = primitive_id;
            }
        }
    }

    return {border_t, border_primitive_id};
}

std::pair<float, int> BVH::get_t(glm::vec3 O, glm::vec3 D) const {
    return get_t_for(_root, O, D, -F_INF, -1);
}
