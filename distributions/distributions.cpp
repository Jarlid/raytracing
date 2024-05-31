#include "distributions.h"

#include <cmath>
#include <utility>

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

LightSource::LightSource(Primitive& primitive) : _primitive(primitive) {
}

glm::vec3 LightSource::sample(glm::vec3 P, glm::vec3 N, RandomEngine& random_engine) {
    return glm::normalize(_primitive.get_random_point(random_engine) - P);
}

float LightSource::pdf(glm::vec3 P, glm::vec3 N, glm::vec3 D) {
    float pdf = 0;
    D = glm::normalize(D);

    float t1, t2;
    std::tie(t1, t2) = _primitive.get_ts(P, D);

    float ts[2];
    int ts_size = 0;

    if (t1 > 0)
        ts[ts_size++] = t1;
    if (t2 > 0)
        ts[ts_size++] = t2;

    for (int ts_i = 0; ts_i < ts_size; ++ts_i) {
        float t = ts[ts_i];

        glm::vec3 tD = t * D;
        glm::vec3 prim_P = P + tD;
        glm::vec3 prim_N = _primitive.get_normal(prim_P);

        float addon = _primitive.get_point_pdf(prim_P) * powf(glm::length(tD), 2) /
                      std::abs((float) glm::dot(D, prim_N));
        if (std::isnan(addon) or std::isinf(addon))
            return F_INF;
        pdf += addon;
    }

    return pdf;
}

Mix::Mix(std::vector<std::unique_ptr<Distribution>> distributions) {
    _inner_distributions = std::move(distributions);
}

glm::vec3 Mix::sample(glm::vec3 P, glm::vec3 N, RandomEngine &random_engine) {
    int i = (int) (_base_distribution(random_engine) * (float) _inner_distributions.size());
    return _inner_distributions[i]->sample(P, N, random_engine);
}

float Mix::pdf(glm::vec3 P, glm::vec3 N, glm::vec3 D) {
    float pdf = 0;

    for (auto & _inner_distribution : _inner_distributions)
        pdf += _inner_distribution->pdf(P, N, D);

    return pdf / (float) _inner_distributions.size();
}

LightSourceMix::LightSourceMix(std::vector<Primitive*> primitives) {
    _primitives = std::move(primitives);
    _bvh = BVH(&_primitives);

    for (auto primitive: _primitives)
        _inner_distributions.push_back(std::make_unique<LightSource>(*primitive));
}

glm::vec3 LightSourceMix::sample(glm::vec3 P, glm::vec3 N, RandomEngine &random_engine) {
    int i = (int) (_base_distribution(random_engine) * (float) _inner_distributions.size());
    return _inner_distributions[i]->sample(P, N, random_engine);
}

float LightSourceMix::pdf(glm::vec3 P, glm::vec3 N, glm::vec3 D) {
    std::vector<int> intersected_primitive_ids = _bvh.get_intersected_primitives(P, D);

    float pdf = 0;
    for (int id: intersected_primitive_ids)
        pdf += _inner_distributions[id]->pdf(P, N, D);

    return pdf / (float) _inner_distributions.size();
}
