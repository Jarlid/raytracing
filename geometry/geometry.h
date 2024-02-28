#pragma once

#include <istream>
#include <vector>

#define GLM_ENABLE_EXPERIMENTAL

#include <glm/glm.hpp>
#include <glm/gtx/quaternion.hpp>

//#include "../glm/glm/vec3.hpp"
//#include "../glm/glm/gtx/quaternion.hpp"

// TODO: убрать утечки памяти

#define GAMMA (1 / 2.2)
#define EPSILON 0.0001

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

    glm::vec3* get_color(glm::vec3 O, glm::vec3 D, float t, const Scene& scene);
};

enum class LightSourceType {
    DIRECTIONAL,
    POINT
};

class LightSource {
private:
    LightSourceType _light_source_type;
    glm::vec3* _light_intensity;

    glm::vec3* _light_direction;

    glm::vec3* _light_position;
    float c0 = 1, c1 = 0, c2 = 0; // light attenuation
public:
    explicit LightSource(std::istream* in_stream);

    std::pair<glm::vec3*, float> light_direction(glm::vec3 P);

    glm::vec3* diffused_light(glm::vec3 P, glm::vec3 N);
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
    std::vector<LightSource*> _light_sources;

    glm::vec3* _ambient_light = new glm::vec3(0);

public:
    explicit Scene(std::ifstream* in_stream);

    std::pair<float, Primitive*> get_t(glm::vec3 O, glm::vec3 D, float eps) const;
    std::pair<float, Primitive*> get_t(glm::vec3 O, glm::vec3 D) const;

    std::vector<uint8_t> render() const;

    std::vector<LightSource*> light_sources() const;
    glm::vec3* ambient_light() const;

    int width() const;
    int height() const;
};
