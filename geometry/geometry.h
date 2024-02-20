#pragma once

#include <istream>
#include <vector>

#define GLM_ENABLE_EXPERIMENTAL

#include <glm/glm.hpp>
#include <glm/gtx/quaternion.hpp>

//#include "../glm/glm/vec3.hpp"
//#include "../glm/glm/gtx/quaternion.hpp"

// TODO: убрать утечки памяти

class Geometry{};

struct Plane: Geometry {
    float nx, ny, nz;

    explicit Plane(std::istream* in_stream);
};

struct Ellipsoid: Geometry {
    float rx, ry, rz;

    explicit Ellipsoid(std::istream* in_stream);
};

struct Box: Geometry {
    float sx, sy, sz;

    explicit Box(std::istream* in_stream);
};

class Primitive {
private:
    Geometry* _geometry;

    glm::vec3* _color;

    glm::vec3* _position = new glm::vec3(0, 0, 0);
    glm::quat* _rotation = new glm::quat(1, 0, 0, 0);

public:
    explicit Primitive(std::istream* in_stream);
};

class Scene {
private:
    int _dimension_width = 0, _dimension_height = 0;

    glm::vec3* _bg_color;

    glm::vec3* _camera_position;

    glm::vec3* _camera_right;
    glm::vec3* _camera_up;
    glm::vec3* _camera_forward;

    float _camera_fov_x = 0;

    std::vector<Primitive*> _primitives;

public:
    explicit Scene(std::ifstream* in_stream);
};