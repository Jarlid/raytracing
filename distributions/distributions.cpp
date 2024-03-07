#include "distributions.h"

#include <cmath>

glm::vec3 UniformHemisphere::sample(glm::vec3 P, glm::vec3 N, RandomEngine& random_engine) {
    float x, y, z;
    do {
        x = _base_distribution(random_engine);
        y = _base_distribution(random_engine);
        z = _base_distribution(random_engine);
    } while (x * x + y * y + z * z > 1 or x == 0 and y == 0 and z == 0);

    glm::vec3 D = glm::normalize(glm::vec3(x, y, z));
    if (glm::dot(D, N) < 0)
        D = -D;
    return D;
}

float UniformHemisphere::pdf(glm::vec3 P, glm::vec3 N, glm::vec3 D) {
    return M_1_PI / 2;
}

glm::vec3 CosineHemisphere::sample(glm::vec3 P, glm::vec3 N, RandomEngine& random_engine) {
    float x, y, z;
    do {
        x = _base_distribution(random_engine);
        y = _base_distribution(random_engine);
        z = _base_distribution(random_engine);
    } while (x * x + y * y + z * z > 1 or x == 0 and y == 0 and z == 0);

    return glm::normalize(glm::normalize(glm::vec3(x, y, z)) + N);
}

float CosineHemisphere::pdf(glm::vec3 P, glm::vec3 N, glm::vec3 D) {
    return fmaxf(0, glm::dot(N, D) / (float) M_PI);
}

LightSource::LightSource(Primitive* primitive) {
    _primitive = primitive;
}

glm::vec3 LightSource::sample(glm::vec3 P, glm::vec3 N, RandomEngine& random_engine) {
    return glm::normalize(_primitive->get_random_point(random_engine) - P);
}

float LightSource::pdf(glm::vec3 P, glm::vec3 N, glm::vec3 D) {
    float pdf = 0;
    D = glm::normalize(D);

    float t1, t2;
    std::tie(t1, t2) = _primitive->get_ts(P, D);

    std::vector<float> ts;
    if (t1 > 0)
        ts.push_back(t1);
    if (t2 > 0)
        ts.push_back(t2);

    for (auto t : ts) {
        glm::vec3 tD = t * D;
        glm::vec3 prim_P = P + tD;
        glm::vec3 prim_N = _primitive->get_normal(prim_P);

        pdf += _primitive->get_point_pdf(P + tD) * powf(glm::length(tD), 2) /
                abs(glm::dot(D, prim_N));
    }

    return pdf;
}
