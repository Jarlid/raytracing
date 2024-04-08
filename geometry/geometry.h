#pragma once

#include <istream>
#include <vector>
#include <tuple>
#include <random>
#include <memory>

#define GLM_ENABLE_EXPERIMENTAL

#include <glm/glm.hpp>
#include <glm/gtx/quaternion.hpp>

#include "distributions/engine.h"

#define GAMMA (1 / 2.2)
#define EPSILON 0.0001f
#define F_INF std::numeric_limits<float>::infinity()

struct Geometry {
private:
    std::uniform_real_distribution<float>* _forbidden_base_distribution = nullptr;
    // TODO: убрать указатель. Возможно, стоит заменить на std::optional.

protected:
    float base_distribution(RandomEngine& random_engine, bool including_negative);

public:
    virtual ~Geometry() = default;

    virtual bool is_plane();

    virtual std::pair<float, float> get_ts(glm::vec3 O, glm::vec3 D) = 0; // Формула: P = O + t * D
    // It's so hard to not joke about TS...
    // But either way, I'm gonna go preorder TTPD.
    virtual glm::vec3 get_normal(glm::vec3 P) = 0;

    virtual glm::vec3 get_random_point(RandomEngine& random_engine) = 0;
    virtual float get_point_pdf(glm::vec3 P) = 0;
};

struct Plane: Geometry {
private:
    glm::vec3 _n{};

public:
    explicit Plane(std::istream& in_stream);

    bool is_plane() override;

    std::pair<float, float> get_ts(glm::vec3 O, glm::vec3 D) override;
    glm::vec3 get_normal(glm::vec3 P) override;

    glm::vec3 get_random_point(RandomEngine& random_engine) override;
    float get_point_pdf(glm::vec3 P) override;
};

struct Ellipsoid: Geometry {
private:
    glm::vec3 _r{};

public:
    explicit Ellipsoid(std::istream& in_stream);

    std::pair<float, float> get_ts(glm::vec3 O, glm::vec3 D) override;
    glm::vec3 get_normal(glm::vec3 P) override;

    glm::vec3 get_random_point(RandomEngine& random_engine) override;
    float get_point_pdf(glm::vec3 P) override;
};

struct Box: Geometry {
private:
    glm::vec3 _s{};

public:
    explicit Box(std::istream& in_stream);

    std::pair<float, float> get_ts(glm::vec3 O, glm::vec3 D) override;
    glm::vec3 get_normal(glm::vec3 P) override;

    glm::vec3 get_random_point(RandomEngine& random_engine) override;
    float get_point_pdf(glm::vec3 P) override;
};

struct Triangle: Geometry {
private:
    glm::vec3 _a{}, _b{}, _c{};

public:
    explicit Triangle(std::istream& in_stream);

    std::pair<float, float> get_ts(glm::vec3 O, glm::vec3 D) override;
    glm::vec3 get_normal(glm::vec3 P) override;

    glm::vec3 get_random_point(RandomEngine& random_engine) override;
    float get_point_pdf(glm::vec3 P) override;
};

enum class Material {
    METALLIC,
    DIELECTRIC,
    DIFFUSIVE
};

class Scene;

class Primitive {
private:
    std::unique_ptr<Geometry> _geometry;

    glm::vec3 _color{0};

    glm::vec3 _position{0};
    glm::quat _rotation{1, 0, 0, 0};

    Material _material = Material::DIFFUSIVE;
    float _ior = 1; // коэффициент преломления (имеет смысл только для диэлектриков)

    glm::vec3 _emission{0};

public:
    explicit Primitive(std::istream& in_stream);

    bool is_plane();
    bool has_emission();

    std::pair<float, float> get_ts(glm::vec3 O, glm::vec3 D);
    float get_t(glm::vec3 O, glm::vec3 D);
    glm::vec3 get_normal(glm::vec3 P);

    glm::vec3 get_color(glm::vec3 O, glm::vec3 D, float t, const Scene& scene, int recursion_depth);

    glm::vec3 get_random_point(RandomEngine& random_engine);
    float get_point_pdf(glm::vec3 P);
};
