#include "distributions.h"

#include <cmath>

glm::vec3* UniformHemisphere::sample(glm::vec3 P, glm::vec3 N, RandomEngine& random_engine) {
    float x, y, z;
    do {
        x = _base_distribution(random_engine);
        y = _base_distribution(random_engine);
        z = _base_distribution(random_engine);
    } while (x * x + y * y + z * z > 1 or x == 0 and y == 0 and z == 0);

    auto D = new glm::vec3(x, y, z);
    *D = glm::normalize(*D);
    if (glm::dot(*D, N) < 0)
        *D = -*D;
    return D;
}

float UniformHemisphere::pdf(glm::vec3 P, glm::vec3 N, glm::vec3 D) {
    return M_1_PI / 2;
}

glm::vec3* CosineHemisphere::sample(glm::vec3 P, glm::vec3 N, RandomEngine& random_engine) {
    float x, y, z;
    do {
        x = _base_distribution(random_engine);
        y = _base_distribution(random_engine);
        z = _base_distribution(random_engine);
    } while (x * x + y * y + z * z > 1 or x == 0 and y == 0 and z == 0);

    auto D = new glm::vec3(x, y, z);
    *D = glm::normalize(*D);
    *D = glm::normalize(*D + N);
    return D;
}

float CosineHemisphere::pdf(glm::vec3 P, glm::vec3 N, glm::vec3 D) {
    return fmaxf(0, glm::dot(N, D) / (float) M_PI);
}
