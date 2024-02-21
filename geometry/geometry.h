#pragma once

#include <istream>
#include <vector>

#define GLM_ENABLE_EXPERIMENTAL

#include <glm/glm.hpp>
#include <glm/gtx/quaternion.hpp>

//#include "../glm/glm/vec3.hpp"
//#include "../glm/glm/gtx/quaternion.hpp"

// TODO: убрать утечки памяти

struct Geometry{
    virtual float get_t(glm::vec3 O, glm::vec3 D) = 0;
    // Формула: P = O + t * D
};

struct Plane: Geometry {
private:
    glm::vec3* _n;

public:
    explicit Plane(std::istream* in_stream);

    float get_t(glm::vec3 O, glm::vec3 D) override;
};

struct Ellipsoid: Geometry {
private:
    glm::vec3* _r;

public:
    explicit Ellipsoid(std::istream* in_stream);

    float get_t(glm::vec3 O, glm::vec3 D) override;
};

struct Box: Geometry {
private:
    glm::vec3* _s;

public:
    explicit Box(std::istream* in_stream);

    float get_t(glm::vec3 O, glm::vec3 D) override;
};

class Primitive {
private:
    Geometry* _geometry;

    glm::vec3* _color;

    glm::vec3* _position = new glm::vec3(0, 0, 0);
    glm::quat* _rotation = new glm::quat(1, 0, 0, 0);

public:
    explicit Primitive(std::istream* in_stream);

    float get_t(glm::vec3 O, glm::vec3 D);

    glm::vec3* color();
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

    std::vector<uint8_t> render() const;

    int width() const;
    int height() const;
};
