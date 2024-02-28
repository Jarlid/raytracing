#include "geometry.h"

#include <cmath>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <tuple>
#include <iostream>

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

void check_color(const glm::vec3 color) {
    if (color.r < 0 or color.g < 0 or color.b < 0)
        exit(1);
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
    return new glm::vec3(glm::normalize(P / *_r));
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
        glm::vec3 new_D = D - 2.f * glm::dot(N, D) * N;

        glm::vec3* color = scene.get_color(P, new_D, EPSILON, recursion_depth + 1);
        *color *= *_color;

        check_color(*color);
        return color;
    }

    if (_material == Material::DIELECTRIC) {
        return _color; // TODO
    }

    // if (_material == Material::DIFFUSER)
    auto total_light = new glm::vec3(*scene.ambient_light());
    for (auto light_source: scene.light_sources()) {
        glm::vec3* light_direction;
        float R;
        std::tie(light_direction, R) = light_source->light_direction(P);

        t = scene.get_t(P, *light_direction, EPSILON).first;
        if (t < EPSILON or t * glm::length(*light_direction) > R)
            *total_light += *light_source->diffused_light(P, N);
    }

    *total_light *= *_color;

    check_color(*total_light);
    return total_light;
}

LightSource::LightSource(std::istream *in_stream) {
    std::string command;

    while (*in_stream >> command) {
        std::transform(command.begin(), command.end(), command.begin(), ::toupper);

        if (command == "LIGHT_INTENSITY")
            _light_intensity = stream_vec3(in_stream);
        else if (command == "LIGHT_DIRECTION") {
            _light_source_type = LightSourceType::DIRECTIONAL;
            _light_direction = stream_vec3(in_stream);
        }
        else if (command == "LIGHT_POSITION") {
            _light_source_type = LightSourceType::POINT;
            _light_position = stream_vec3(in_stream);
        }
        else if (command == "LIGHT_ATTENUATION") {
            _light_source_type = LightSourceType::POINT;
            *in_stream >> c0 >> c1 >> c2;
        }
    }
}

std::pair<glm::vec3*, float> LightSource::light_direction(glm::vec3 P) {
    if (_light_source_type == LightSourceType::DIRECTIONAL)
        return {_light_direction, std::numeric_limits<float>::infinity()};
    auto light_direction = new glm::vec3(*_light_position - P);
    return {light_direction, glm::length(*light_direction)};
}

glm::vec3 *LightSource::diffused_light(glm::vec3 P, glm::vec3 N) {
    glm::vec3* light_direction_now;
    float R;
    std::tie(light_direction_now, R) = light_direction(P);

    glm::vec3 light_intensity = _light_source_type == LightSourceType::POINT ?
                                *_light_intensity / (c0 + c1 * R + c2 * R * R) : *_light_intensity;

    float cos = glm::dot(glm::normalize(*light_direction_now), glm::normalize(N));
    if (cos < 0)
        return new glm::vec3(0);
    return new glm::vec3(cos * light_intensity);
}

Scene::Scene(std::ifstream* in_stream) {
    std::string command;
    std::stringstream* primitive_stream;
    std::stringstream* light_stream;

    bool new_primitive = false;
    bool new_light = false;

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
        else if (command == "AMBIENT_LIGHT") {
            not_primitive_or_light = true;
            _ambient_light = stream_vec3(in_stream);
        }
        else if (command == "RAY_DEPTH") {
            not_primitive_or_light = true;
            *in_stream >> _ray_depth;
        }

        if (not_primitive_or_light || command == "NEW_PRIMITIVE" || command == "NEW_LIGHT") {
            if (new_primitive) {
                new_primitive = false;
                _primitives.push_back(new Primitive(primitive_stream));
            }
            if (new_light) {
                new_light = false;
                _light_sources.push_back(new LightSource(light_stream));
            }
        }
        if (not_primitive_or_light)
            continue;

        if (command == "NEW_PRIMITIVE") {
            new_primitive = true;
            primitive_stream = new std::stringstream();
        }
        else if (command == "NEW_LIGHT") {
            new_light = true;
            light_stream = new std::stringstream();
        }
        else if (new_primitive)
            *primitive_stream << command << " ";
        else if (new_light)
            *light_stream << command << " ";
    }

    if (new_primitive)
        _primitives.push_back(new Primitive(primitive_stream));
}

std::pair<float, Primitive*> Scene::get_t(glm::vec3 O, glm::vec3 D, float eps) const {
    float t = -1;
    Primitive* curr_primitive = nullptr;

    for (auto primitive : _primitives) {
        float new_t = primitive->get_t(O, D);
        if (t <= eps or (new_t > eps and new_t < t)) {
            t = new_t;
            curr_primitive = primitive;
        }
    }

    return {t, curr_primitive};
}

glm::vec3* Scene::get_color(glm::vec3 O, glm::vec3 D, float eps, int recursion_depth) const {
    if (recursion_depth > _ray_depth)
        return new glm::vec3(0);

    float t;
    Primitive* curr_primitive;
    std::tie(t, curr_primitive) = get_t(O, D, eps);

    glm::vec3* curr_color;
    if (t <= eps || curr_primitive == nullptr)
        curr_color = _bg_color;
    else
        curr_color = curr_primitive->get_color(O, D, t, *this, recursion_depth);

    check_color(*curr_color);
    return new glm::vec3(*curr_color);
}

glm::vec3* Scene::get_color(glm::vec3 O, glm::vec3 D) const {
    return get_color(O, D, 0, 1);
}

std::vector<uint8_t> Scene::render() const {
    std::vector<uint8_t> render(3 * _dimension_width * _dimension_height, 0);

    float tan_half_fov_x = tanf(_camera_fov_x / 2);
    float tan_half_fov_y = tan_half_fov_x / (float) _dimension_width * (float) _dimension_height;

    for (int p_x = 0; p_x < _dimension_width; ++p_x) {
        for (int p_y = 0; p_y < _dimension_height; ++p_y) {
            float c_x = (2 * ((float) p_x + 0.5f) / (float) _dimension_width - 1) * tan_half_fov_x;
            float c_y = -1 * (2 * ((float) p_y + 0.5f) / (float) _dimension_height - 1) * tan_half_fov_y;
            float c_z = 1;

            glm::vec3 D = glm::normalize(c_x * *_camera_right + c_y * *_camera_up + c_z * *_camera_forward);
            glm::vec3 O = *_camera_position;

            glm::vec3* curr_color = get_color(O, D);
            curr_color = new glm::vec3(aces_tonemap(*curr_color));

            int it = 3 * (p_x + p_y * _dimension_width);
            render[it + 0] = fix_color(curr_color->r);
            render[it + 1] = fix_color(curr_color->g);
            render[it + 2] = fix_color(curr_color->b);
        }
    }

    return render;
}

std::vector<LightSource*> Scene::light_sources() const {
    return _light_sources;
}

glm::vec3 *Scene::ambient_light() const {
    return _ambient_light;
}

int Scene::width() const {
    return _dimension_width;
}

int Scene::height() const {
    return _dimension_height;
}
