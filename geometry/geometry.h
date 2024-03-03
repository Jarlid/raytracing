#pragma once

#include <istream>
#include <vector>
#include <random>

#define GLM_ENABLE_EXPERIMENTAL

#include <glm/glm.hpp>
#include <glm/gtx/quaternion.hpp>

//#include "../glm/glm/vec3.hpp"
//#include "../glm/glm/gtx/quaternion.hpp"

// TODO: убрать утечки памяти

#define GAMMA (1 / 2.2)
#define EPSILON 0.0001f

struct Geometry{
    virtual float get_t(glm::vec3 O, glm::vec3 D) = 0;
    virtual glm::vec3* get_normal(glm::vec3 P) = 0;
    // Формула: P = O + t * D
};

struct Plane: Geometry {
private:
    glm::vec3* _n;

public:
    explicit Plane(std::istream* in_stream);

    float get_t(glm::vec3 O, glm::vec3 D) override;
    glm::vec3* get_normal(glm::vec3 P) override;
};

struct Ellipsoid: Geometry {
private:
    glm::vec3* _r;

public:
    explicit Ellipsoid(std::istream* in_stream);

    float get_t(glm::vec3 O, glm::vec3 D) override;
    glm::vec3* get_normal(glm::vec3 P) override;
};

struct Box: Geometry {
private:
    glm::vec3* _s;

public:
    explicit Box(std::istream* in_stream);

    float get_t(glm::vec3 O, glm::vec3 D) override;
    glm::vec3* get_normal(glm::vec3 P) override;
};

enum class Material {
    METALLIC,
    DIELECTRIC,
    DIFFUSER
};

class Scene;

class Primitive {
private:
    Geometry* _geometry;

    glm::vec3* _color;

    glm::vec3* _position = new glm::vec3(0, 0, 0);
    glm::quat* _rotation = new glm::quat(1, 0, 0, 0);

    Material _material = Material::DIFFUSER;
    float _ior = 1; // коэффициент преломления (имеет смысл только для диэлектриков)

public:
    explicit Primitive(std::istream* in_stream);

    float get_t(glm::vec3 O, glm::vec3 D);
    glm::vec3* get_normal(glm::vec3 P);

    glm::vec3* get_color(glm::vec3 O, glm::vec3 D, float t, const Scene& scene, int recursion_depth);
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

    int _ray_depth = 1;

    int _sample_num = 1;

    std::default_random_engine* _random_engine = new std::default_random_engine();

public:
    explicit Scene(std::ifstream* in_stream);

    std::pair<float, Primitive*> get_t(glm::vec3 O, glm::vec3 D) const;

    glm::vec3* get_color(glm::vec3 O, glm::vec3 D, int recursion_depth) const;
    glm::vec3* get_color(glm::vec3 O, glm::vec3 D) const;

    std::vector<uint8_t> render() const;

    int width() const;
    int height() const;

    std::default_random_engine* random_engine() const;
};
