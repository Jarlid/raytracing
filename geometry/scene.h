#pragma once

#include "distributions/distributions.h"

class Scene {
private:
    int _dimension_width = 0, _dimension_height = 0;

    glm::vec3 _bg_color;

    glm::vec3 _camera_position;

    glm::vec3 _camera_right;
    glm::vec3 _camera_up;
    glm::vec3 _camera_forward;

    float _camera_fov_x = 0;

    std::vector<Primitive*> _primitives;

    int _ray_depth = 1;

    int _sample_num = 1;

    RandomEngine* _random_engine = new RandomEngine();

    Distribution* _distribution;

public:
    explicit Scene(std::ifstream* in_stream);

    std::pair<float, Primitive*> get_t(glm::vec3 O, glm::vec3 D) const;

    glm::vec3 get_color(glm::vec3 O, glm::vec3 D, int recursion_depth) const;
    glm::vec3 get_color(glm::vec3 O, glm::vec3 D) const;

    std::vector<uint8_t> render() const;

    int width() const;
    int height() const;

    RandomEngine* random_engine() const;
    Distribution* distribution() const;
};
