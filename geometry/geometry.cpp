#include "geometry.h"

#include <cmath>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <tuple>
#include <random>

glm::vec3* stream_vec3(std::istream* in_stream) {
    float x, y, z = 0;
    *in_stream >> x >> y >> z;
    return new glm::vec3(x, y, z);
}

glm::quat* stream_quat(std::istream* in_stream) {
    float x, y, z, w = 0;
    *in_stream >> x >> y >> z >> w;
    return new glm::quat(w, x, y, z);
}

float pick_t(float t1, float t2) {
    if (t1 > t2)
        std::swap(t1, t2);

    if (t1 <= 0)
        return t2;
    return t1;
}

glm::vec3 saturate(glm::vec3 const & color) {
    return glm::clamp(color, glm::vec3(0.0), glm::vec3(1.0));
}

glm::vec3 aces_tonemap(glm::vec3 const & x) {
    return saturate((x * (2.51f * x + 0.03f)) / (x * (2.43f * x + 0.59f) + 0.14f));
}

uint8_t fix_color(float value) {
    auto gamma_corrected = powf(value, GAMMA);
    auto tmp = std::round(gamma_corrected * 255);
    if (tmp > 255)
        tmp = 255;
    return (uint8_t) tmp;
}

void check_color(const glm::vec3 color, const int error_code) {
    if (color.r < 0 or color.g < 0 or color.b < 0)
        exit(error_code);
}

Plane::Plane(std::istream* in_stream) {
    _n = stream_vec3(in_stream);
}

float Plane::get_t(glm::vec3 O, glm::vec3 D) {
    return -1 * glm::dot(O, *_n) / glm::dot(D, *_n);
}

glm::vec3* Plane::get_normal(glm::vec3 P) {
    return new glm::vec3(glm::normalize(*_n));
}

Ellipsoid::Ellipsoid(std::istream* in_stream) {
    _r = stream_vec3(in_stream);
}

float Ellipsoid::get_t(glm::vec3 O, glm::vec3 D) {
    float a = glm::dot(D / *_r, D / *_r);
    float b = 2 * glm::dot(O / *_r, D / *_r);
    float c = glm::dot(O / *_r, O / *_r) - 1;
    float d = b * b - 4 * a * c;

    if (d < 0)
        return -1;
    if (d == 0)
        return -0.5f * b / a;

    float t1 = (-1 * b + sqrtf(d)) / (2 * a), t2 = (-1 * b - sqrtf(d)) / (2 * a);
    return pick_t(t1, t2);
}

glm::vec3* Ellipsoid::get_normal(glm::vec3 P) {
    return new glm::vec3(glm::normalize(P / *_r / *_r));
}

Box::Box(std::istream* in_stream) {
    _s = stream_vec3(in_stream);
}

float Box::get_t(glm::vec3 O, glm::vec3 D) {
    glm::vec3 T1 = (*_s - O) / D, T2 = (- *_s - O) / D;

    float tx1 = fminf(T1.x, T2.x), tx2 = fmaxf(T1.x, T2.x);
    float ty1 = fminf(T1.y, T2.y), ty2 = fmaxf(T1.y, T2.y);
    float tz1 = fminf(T1.z, T2.z), tz2 = fmaxf(T1.z, T2.z);

    float t1 = fmaxf(tx1, fmaxf(ty1, tz1));
    float t2 = fminf(tx2, fminf(ty2, tz2));

    if (t1 > t2)
        return -1;

    return pick_t(t1, t2);
}

glm::vec3* Box::get_normal(glm::vec3 P) {
    glm::vec3 nP = P / *_s;
    glm::vec3 anP = glm::abs(P / *_s);

    if (anP.x >= anP.y and anP.x >= anP.z) {
        if (nP.x > 0)
            return new glm::vec3(1, 0, 0);
        return new glm::vec3(-1, 0, 0);
    }

    if (anP.y >= anP.x and anP.y >= anP.z){
        if (nP.y > 0)
            return new glm::vec3(0, 1, 0);
        return new glm::vec3(0, -1, 0);
    }

    if (nP.z > 0)
        return new glm::vec3(0, 0, 1);
    return new glm::vec3(0, 0, -1);
}

Primitive::Primitive(std::istream* in_stream) {
    std::string command;

    while (*in_stream >> command) {
        std::transform(command.begin(), command.end(), command.begin(), ::toupper);

        if (command == "PLANE")
            _geometry = new Plane(in_stream);
        else if (command == "ELLIPSOID")
            _geometry = new Ellipsoid(in_stream);
        else if (command == "BOX")
            _geometry = new Box(in_stream);
        else if (command == "COLOR")
            _color = stream_vec3(in_stream);
        else if (command == "POSITION")
            _position = stream_vec3(in_stream);
        else if (command == "ROTATION")
            _rotation = stream_quat(in_stream);
        else if (command == "METALLIC")
            _material = Material::METALLIC;
        else if (command == "DIELECTRIC")
            _material = Material::DIELECTRIC;
        else if (command == "IOR")
            *in_stream >> _ior;
        else if (command == "EMISSION")
            _emission = stream_vec3(in_stream);
    }
}

float Primitive::get_t(glm::vec3 O, glm::vec3 D) {
    O -= *_position;

    glm::quat reverse_rotation = glm::conjugate(*_rotation);

    O = glm::rotate(reverse_rotation, O);
    D = glm::rotate(reverse_rotation, D);

    return _geometry->get_t(O, D);
}

glm::vec3* Primitive::get_normal(glm::vec3 P) {
    glm::quat reverse_rotation = glm::conjugate(*_rotation);

    P -= *_position;
    P = glm::rotate(reverse_rotation, P);

    return new glm::vec3(glm::normalize(glm::rotate(*_rotation, *_geometry->get_normal(P))));
}

glm::vec3 *Primitive::get_color(glm::vec3 O, glm::vec3 D, float t, const Scene& scene, int recursion_depth) {
    glm::vec3 P = O + t * D;
    glm::vec3 N = *get_normal(P);
    bool inside = false;

    if (glm::dot(D, N) > 0) {
        N = -N;
        inside = true;
    }

    if (_material == Material::METALLIC) {
        glm::vec3 new_D = glm::normalize(D - 2.f * glm::dot(N, D) * N);

        glm::vec3 reflected_color = *_emission +
                *scene.get_color(P + new_D * EPSILON, new_D, recursion_depth + 1);
        reflected_color *= *_color;

        check_color(reflected_color, 348765);
        return new glm::vec3(reflected_color);
    }

    if (_material == Material::DIELECTRIC) {
        float cos_t1 = - glm::dot(N, glm::normalize(D));
        float sin_t1 = sqrtf(1 - cos_t1 * cos_t1);

        float n1 = 1, n2 = _ior;
        if (inside)
            std::swap(n1, n2);

        float sin_t2 = sin_t1 * n1 / n2;

        float R0 = powf((n1 - n2) / (n1 + n2), 2);
        float R = R0 + (1 - R0) * powf(1 - cos_t1, 5);

        std::uniform_real_distribution<float> distribution(0, 1);
        float U = distribution(*scene.random_engine());

        if (sin_t2 > 1 or U < R) {
            glm::vec3 reflected_D = glm::normalize(D - 2.f * glm::dot(N, D) * N);
            glm::vec3 reflected_color = *_emission +
                    *scene.get_color(P + reflected_D * EPSILON, reflected_D, recursion_depth + 1);

            check_color(reflected_color, 397234);
            return new glm::vec3(reflected_color);
        }

        float cos_t2 = sqrtf(1 - sin_t2 * sin_t2);

        glm::vec3 refracted_D = glm::normalize(D * n1 / n2 + N * (cos_t1 * n1 / n2 - cos_t2));
        glm::vec3 refracted_color = *_emission +
                *scene.get_color(P + refracted_D * EPSILON, refracted_D, recursion_depth + 1);

        if (not inside)
            refracted_color *= *_color;

        check_color(refracted_color, 387463);
        return new glm::vec3(refracted_color);
    }

    // if (_material == Material::DIFFUSER)
    auto distribution = UniformHemisphere();
    glm::vec3 new_D = *distribution.sample(P, N, *scene.random_engine());
    float pdf = distribution.pdf(P, N, new_D);

    glm::vec3* L_in = scene.get_color(P + new_D * EPSILON, new_D, recursion_depth + 1);
    glm::vec3 color = *_emission + *_color * *L_in * glm::dot(new_D, N) / (float) M_PI / pdf;

    check_color(color, 873464);
    return new glm::vec3(color);
}

Scene::Scene(std::ifstream* in_stream) {
    std::string command;
    std::stringstream* primitive_stream;

    bool new_primitive = false;

    while (*in_stream >> command) {
        std::transform(command.begin(), command.end(), command.begin(), ::toupper);

        bool not_primitive_or_light = false;

        if (command == "DIMENSIONS") {
            not_primitive_or_light = true;
            *in_stream >> _dimension_width >> _dimension_height;
        }
        else if (command == "BG_COLOR") {
            not_primitive_or_light = true;
            _bg_color = stream_vec3(in_stream);
        }
        else if (command == "CAMERA_POSITION") {
            not_primitive_or_light = true;
            _camera_position = stream_vec3(in_stream);
        }
        else if (command == "CAMERA_RIGHT" || command == "CAMERA_WRONG") {
            // We support not only camera rights, but also camera wrongs!
            not_primitive_or_light = true;
            _camera_right = stream_vec3(in_stream);
        }
        else if (command == "CAMERA_UP") {
            not_primitive_or_light = true;
            _camera_up = stream_vec3(in_stream);
        }
        else if (command == "CAMERA_FORWARD") {
            not_primitive_or_light = true;
            _camera_forward = stream_vec3(in_stream);
        }
        else if (command == "CAMERA_FOV_X") {
            not_primitive_or_light = true;
            *in_stream >> _camera_fov_x;
        }
        else if (command == "RAY_DEPTH") {
            not_primitive_or_light = true;
            *in_stream >> _ray_depth;
        }
        else if (command == "SAMPLES") {
            not_primitive_or_light = true;
            *in_stream >> _sample_num;
        }

        if (not_primitive_or_light || command == "NEW_PRIMITIVE") {
            if (new_primitive) {
                new_primitive = false;
                _primitives.push_back(new Primitive(primitive_stream));
            }
        }
        if (not_primitive_or_light)
            continue;

        if (command == "NEW_PRIMITIVE") {
            new_primitive = true;
            primitive_stream = new std::stringstream();
        }
        else if (new_primitive)
            *primitive_stream << command << " ";
    }

    if (new_primitive)
        _primitives.push_back(new Primitive(primitive_stream));
}

std::pair<float, Primitive*> Scene::get_t(glm::vec3 O, glm::vec3 D) const {
    float t = -1;
    Primitive* curr_primitive = nullptr;

    for (auto primitive : _primitives) {
        float new_t = primitive->get_t(O, D);
        if (t <= 0 or (new_t > 0 and new_t < t)) {
            t = new_t;
            curr_primitive = primitive;
        }
    }

    return {t, curr_primitive};
}

glm::vec3* Scene::get_color(glm::vec3 O, glm::vec3 D, int recursion_depth) const {
    if (recursion_depth > _ray_depth)
        return new glm::vec3(0);

    D = glm::normalize(D);

    float t;
    Primitive* curr_primitive;
    std::tie(t, curr_primitive) = get_t(O, D);

    glm::vec3* curr_color;
    if (t <= 0 || curr_primitive == nullptr)
        curr_color = _bg_color;
    else
        curr_color = curr_primitive->get_color(O, D, t, *this, recursion_depth);

    check_color(*curr_color, 937433);
    return new glm::vec3(*curr_color);
}

glm::vec3* Scene::get_color(glm::vec3 O, glm::vec3 D) const {
    return get_color(O, D, 1);
}

std::vector<uint8_t> Scene::render() const {
    std::vector<uint8_t> render(3 * _dimension_width * _dimension_height, 0);

    float tan_half_fov_x = tanf(_camera_fov_x / 2);
    float tan_half_fov_y = tan_half_fov_x / (float) _dimension_width * (float) _dimension_height;

    std::uniform_real_distribution<float> distribution(0, 1);

    for (int p_x = 0; p_x < _dimension_width; ++p_x) {
        for (int p_y = 0; p_y < _dimension_height; ++p_y) {
            auto color_sum = new glm::vec3(0);

            for (int curr_sample = 0; curr_sample < _sample_num; ++curr_sample) {
                float u_x = distribution(*_random_engine);
                float u_y = distribution(*_random_engine);

                float c_x = (2 * ((float) p_x + u_x) / (float) _dimension_width - 1) * tan_half_fov_x;
                float c_y = -1 * (2 * ((float) p_y + u_y) / (float) _dimension_height - 1) * tan_half_fov_y;
                float c_z = 1;

                glm::vec3 D = glm::normalize(c_x * *_camera_right + c_y * *_camera_up + c_z * *_camera_forward);
                glm::vec3 O = *_camera_position;

                *color_sum += *get_color(O, D);
            }

            *color_sum /= _sample_num;
            check_color(*color_sum, 384763);

            color_sum = new glm::vec3(aces_tonemap(*color_sum));

            int it = 3 * (p_x + p_y * _dimension_width);
            render[it + 0] = fix_color(color_sum->r);
            render[it + 1] = fix_color(color_sum->g);
            render[it + 2] = fix_color(color_sum->b);
        }
    }

    return render;
}

int Scene::width() const {
    return _dimension_width;
}

int Scene::height() const {
    return _dimension_height;
}

std::default_random_engine* Scene::random_engine() const {
    return _random_engine;
}
