#pragma once

#include <glm/glm.hpp>
#include <glm/gtx/quaternion.hpp>

#define F_INF std::numeric_limits<float>::infinity()

class AABB {
private:
    glm::vec3 _min = glm::vec3(F_INF, F_INF, F_INF);
    glm::vec3 _max = glm::vec3(-F_INF, -F_INF, -F_INF);

    bool _is_infinite = false;

public:
    AABB() = default;

    void make_infinite();

    void replace_with(AABB aabb);

    void extend(glm::vec3 point);
    void extend(AABB aabb);

    void move(glm::vec3 vector);
    void rotate(glm::quat quat);

    std::tuple<int, float> get_cut();
    bool divide(int cut_axis, float cut_coord);

    bool is_intersected(glm::vec3 O, glm::vec3 D, float border_t) const;
};