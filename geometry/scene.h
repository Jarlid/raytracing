#pragma once

#include <memory>

#include <rapidjson/include/rapidjson/document.h>

#include "distributions/distributions.h"

#define DEFAULT_RAY_DEPTH 6

class Scene {
private:
    int _dimension_width = 0, _dimension_height = 0;

    glm::vec3 _bg_color{0};

    glm::vec3 _camera_position{};

    glm::vec3 _camera_right{1, 0, 0};
    glm::vec3 _camera_up{0, 1, 0};
    glm::vec3 _camera_forward{0, 0, -1};

    float _camera_fov_y = 0;

    bool _camera_initialized = false;

    std::vector<std::unique_ptr<Primitive>> _primitives;

    int _ray_depth = DEFAULT_RAY_DEPTH;

    std::vector<std::unique_ptr<RandomEngine>> _random_engines;

    std::unique_ptr<Distribution> _distribution;

    std::vector<std::vector<unsigned char>> _buffers;

public:
    void initialize_node(const rapidjson::Document& document, int node_num, glm::mat4 current_transform);
    void initialize_node(const rapidjson::Document& document, int node_num);

    explicit Scene(const rapidjson::Document& document, int width, int height);

    std::pair<float, Primitive*> get_t(glm::vec3 O, glm::vec3 D) const;

    glm::vec3 get_color(glm::vec3 O, glm::vec3 D, int recursion_depth) const;
    glm::vec3 get_color(glm::vec3 O, glm::vec3 D) const;

    std::vector<uint8_t> render(int samples) const;

    int width() const;
    int height() const;

    RandomEngine* random_engine(int thread) const;
    Distribution* distribution() const;
};
