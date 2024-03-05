#pragma once

#include <random>

#include <glm/glm.hpp>

typedef std::default_random_engine RandomEngine;

struct Distribution {
    virtual glm::vec3* sample(glm::vec3 P, glm::vec3 N, RandomEngine& random_engine) = 0;
    virtual float pdf(glm::vec3 P, glm::vec3 N, glm::vec3 D) = 0;
};

class UniformHemisphere: Distribution {
    std::uniform_real_distribution<float> _base_distribution = std::uniform_real_distribution<float>(-1, 1);

public:
    glm::vec3* sample(glm::vec3 P, glm::vec3 N, RandomEngine& random_engine) override;
    float pdf(glm::vec3 P, glm::vec3 N, glm::vec3 D) override;
};

class CosineHemisphere: Distribution {
    std::uniform_real_distribution<float> _base_distribution = std::uniform_real_distribution<float>(-1, 1);

public:
    glm::vec3* sample(glm::vec3 P, glm::vec3 N, RandomEngine& random_engine) override;
    float pdf(glm::vec3 P, glm::vec3 N, glm::vec3 D) override;
};
