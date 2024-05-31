#pragma once

#include <random>

#include <glm/glm.hpp>

#include "engine.h"
#include "geometry/geometry.h"
#include "hierarchy/hierarchy.h"

struct Distribution {
    virtual ~Distribution() = default;

    virtual glm::vec3 sample(glm::vec3 P, glm::vec3 N, RandomEngine& random_engine) = 0;
    virtual float pdf(glm::vec3 P, glm::vec3 N, glm::vec3 D) = 0;
};

struct UniformHemisphere: Distribution {
private:
    std::uniform_real_distribution<float> _base_distribution = std::uniform_real_distribution<float>(-1, 1);

public:
    glm::vec3 sample(glm::vec3 P, glm::vec3 N, RandomEngine& random_engine) override;
    float pdf(glm::vec3 P, glm::vec3 N, glm::vec3 D) override;
};

struct CosineHemisphere: Distribution {
private:
    std::uniform_real_distribution<float> _base_distribution = std::uniform_real_distribution<float>(-1, 1);

public:
    glm::vec3 sample(glm::vec3 P, glm::vec3 N, RandomEngine& random_engine) override;
    float pdf(glm::vec3 P, glm::vec3 N, glm::vec3 D) override;
};

struct LightSource: Distribution {
private:
    Primitive& _primitive;

public:
    explicit LightSource(Primitive& primitive);

    glm::vec3 sample(glm::vec3 P, glm::vec3 N, RandomEngine& random_engine) override;
    float pdf(glm::vec3 P, glm::vec3 N, glm::vec3 D) override;
};

struct Mix: Distribution {
private:
    std::vector<std::unique_ptr<Distribution>> _inner_distributions;
    std::uniform_real_distribution<float> _base_distribution = std::uniform_real_distribution<float>(0, 1);

public:
    explicit Mix(std::vector<std::unique_ptr<Distribution>> distributions);

    glm::vec3 sample(glm::vec3 P, glm::vec3 N, RandomEngine& random_engine) override;
    float pdf(glm::vec3 P, glm::vec3 N, glm::vec3 D) override;
};

struct LightSourceMix: Distribution {
private:
    std::vector<Primitive*> _primitives;
    BVH _bvh;

    std::vector<std::unique_ptr<Distribution>> _inner_distributions;
    std::uniform_real_distribution<float> _base_distribution = std::uniform_real_distribution<float>(0, 1);

public:
    explicit LightSourceMix(std::vector<Primitive*> primitives);

    glm::vec3 sample(glm::vec3 P, glm::vec3 N, RandomEngine& random_engine) override;
    float pdf(glm::vec3 P, glm::vec3 N, glm::vec3 D) override;
};
