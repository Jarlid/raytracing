#include "geometry.h"

#include <cmath>
#include <fstream>
#include <sstream>
#include <algorithm>

glm::vec3* vec3(std::istream* in_stream) {
    float x, y, z = 0;
    *in_stream >> x >> y >> z;
    return new glm::vec3(x, y, z);
}

glm::quat* quat(std::istream* in_stream) {
    float x, y, z, w = 0;
    *in_stream >> x >> y >> z >> w;
    return new glm::quat(w, x, y, z);
}

float scalar(const glm::vec3 a, const glm::vec3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

glm::quat* conjugate_quat(const glm::quat q) {
    return new glm::quat(q.w, -q.x, -q.y, -q.z);
}

glm::vec3* rotate(const glm::vec3 v, const glm::quat q) {
    glm::quat tmp = q * glm::quat(0, v) * *conjugate_quat(q);
    return new glm::vec3(tmp.x, tmp.y, tmp.z);
}

float pick_t(float t1, float t2) {
    if (t1 > t2)
        std::swap(t1, t2);

    if (t1 <= 0)
        return t2;
    return t1;
}

uint8_t fix_color(float value) {
    auto tmp = uint8_t(std::round(value * 255));
    if (tmp > 255)
        tmp = 255;
    return tmp;
}

Plane::Plane(std::istream* in_stream) {
    _n = vec3(in_stream);
}

float Plane::get_t(glm::vec3 O, glm::vec3 D) {
    return -1 * scalar(O, *_n) / scalar(D, *_n);
}

Ellipsoid::Ellipsoid(std::istream* in_stream) {
    _r = vec3(in_stream);
}

float Ellipsoid::get_t(glm::vec3 O, glm::vec3 D) {
    float a = scalar(D / *_r, D / *_r);
    float b = 2 * scalar(O / *_r, D / *_r);
    float c = scalar(O / *_r, O / *_r) - 1;
    float d = b * b - 4 * a * c;

    if (d < 0)
        return -1;
    if (d == 0)
        return -0.5f * b / a;

    float t1 = (-1 * b + sqrtf(d)) / (2 * a), t2 = (-1 * b - sqrtf(d)) / (2 * a);
    return pick_t(t1, t2);
}

Box::Box(std::istream* in_stream) {
    _s = vec3(in_stream);
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
            _color = vec3(in_stream);
        else if (command == "POSITION")
            _position = vec3(in_stream);
        else if (command == "ROTATION")
            _rotation = quat(in_stream);
    }
}

float Primitive::get_t(glm::vec3 O, glm::vec3 D) {
    O -= *_position;

    rotate(O, *conjugate_quat(*_rotation));
    rotate(D, *conjugate_quat(*_rotation));

    return _geometry->get_t(O, D);
}

glm::vec3* Primitive::color() {
    return _color;
}

Scene::Scene(std::ifstream* in_stream) {
    std::string command;
    std::stringstream* primitive_stream;

    bool new_primitive = false;

    while (*in_stream >> command) {
        std::transform(command.begin(), command.end(), command.begin(), ::toupper);

        bool not_primitive = false;

        if (command == "DIMENSIONS") {
            not_primitive = true;
            *in_stream >> _dimension_width >> _dimension_height;
        }
        else if (command == "BG_COLOR") {
            not_primitive = true;
            _bg_color = vec3(in_stream);
        }
        else if (command == "CAMERA_POSITION") {
            not_primitive = true;
            _camera_position = vec3(in_stream);
        }
        else if (command == "CAMERA_RIGHT" || command == "CAMERA_WRONG") {
            // We support not only camera rights, but also camera wrongs!
            not_primitive = true;
            _camera_right = vec3(in_stream);
        }
        else if (command == "CAMERA_UP") {
            not_primitive = true;
            _camera_up = vec3(in_stream);
        }
        else if (command == "CAMERA_FORWARD") {
            not_primitive = true;
            _camera_forward = vec3(in_stream);
        }
        else if (command == "CAMERA_FOV_X") {
            not_primitive = true;
            *in_stream >> _camera_fov_x;
        }

        if ((not_primitive || command == "NEW_PRIMITIVE") && new_primitive) {
            new_primitive = false;
            _primitives.push_back(new Primitive(primitive_stream));
        }
        if (not_primitive)
            continue;

        if (command == "NEW_PRIMITIVE") {
            new_primitive = true;
            primitive_stream = new std::stringstream();
        }
        else if (new_primitive) {
            *primitive_stream << command;
            *primitive_stream << " ";
        }
    }

    if (new_primitive)
        _primitives.push_back(new Primitive(primitive_stream));
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

            float t = -1;
            glm::vec3* curr_color;

            for (auto primitive : _primitives) {
                float new_t = primitive->get_t(O, D);
                if (t <= 0 or (new_t > 0 and new_t < t)) {
                    t = new_t;
                    curr_color = primitive->color();
                }
            }
            if (t <= 0) {
                curr_color = _bg_color;
            }

            int it = 3 * (p_x + p_y * _dimension_width);
            render[it + 0] = fix_color(curr_color->r);
            render[it + 1] = fix_color(curr_color->g);
            render[it + 2] = fix_color(curr_color->b);
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
