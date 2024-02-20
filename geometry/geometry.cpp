#include "geometry.h"

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

Plane::Plane(std::istream* in_stream) {
    nx = ny = nz = 0;
    *in_stream >> nx >> ny >> nz;
}

Ellipsoid::Ellipsoid(std::istream* in_stream) {
    rx = ry = rz = 0;
    *in_stream >> rx >> ry >> rz;
}

Box::Box(std::istream* in_stream) {
    sx = sy = sz = 0;
    *in_stream >> sx >> sy >> sz;
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
