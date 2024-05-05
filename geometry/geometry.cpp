#include "geometry.h"
#include "scene.h"

#include <cmath>
#include <fstream>
#include <memory>
#include <sstream>
#include <algorithm>
#include <tuple>

#include <omp.h>

glm::vec3 stream_vec3(std::istream& in_stream) {
    float x, y, z = 0;
    in_stream >> x >> y >> z;
    return {x, y, z};
}

glm::quat stream_quat(std::istream& in_stream) {
    float x, y, z, w = 0;
    in_stream >> x >> y >> z >> w;
    return {w, x, y, z};
}

glm::vec3 fast_rotate(glm::quat quat, glm::vec3 vec) {
    glm::vec3 cut_quat = {quat.x, quat.y, quat.z};
    glm::vec3 t = 2.f * glm::cross(cut_quat, vec);
    return vec + quat.w * t + glm::cross(cut_quat, t);
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

float Geometry::base_distribution(RandomEngine& random_engine, bool including_negative) {
    if (_forbidden_base_distribution == nullptr)
        _forbidden_base_distribution = new std::uniform_real_distribution<float>(-1, 1);
    if (including_negative)
        return (*_forbidden_base_distribution)(random_engine);
    return std::abs((*_forbidden_base_distribution)(random_engine));
}

bool Geometry::is_plane() {
    return false;
}

Plane::Plane(std::istream& in_stream) {
    _n = stream_vec3(in_stream);
}

bool Plane::is_plane() {
    return true;
}

std::pair<float, float> Plane::get_ts(glm::vec3 O, glm::vec3 D) {
    return {-F_INF, -1 * glm::dot(O, _n) / glm::dot(D, _n)};
}

glm::vec3 Plane::get_normal(glm::vec3 P) {
    return glm::normalize(_n);
}

glm::vec3 Plane::get_random_point(RandomEngine& random_engine) {
    exit(498723);
}

float Plane::get_point_pdf(glm::vec3 P) {
    exit(498724);
}

Ellipsoid::Ellipsoid(std::istream& in_stream) {
    _r = stream_vec3(in_stream);
}

std::pair<float, float> Ellipsoid::get_ts(glm::vec3 O, glm::vec3 D) {
    float a = glm::dot(D / _r, D / _r);
    float b = 2 * glm::dot(O / _r, D / _r);
    float c = glm::dot(O / _r, O / _r) - 1;
    float d = b * b - 4 * a * c;

    if (d < 0)
        return {-F_INF, -F_INF};
    if (d == 0)
        return {-F_INF, -0.5f * b / a};

    float t1 = (-1 * b + sqrtf(d)) / (2 * a), t2 = (-1 * b - sqrtf(d)) / (2 * a);
    if (t1 < t2)
        return {t1, t2};
    return {t2, t1};
}

glm::vec3 Ellipsoid::get_normal(glm::vec3 P) {
    return glm::normalize(P / _r / _r);
}

glm::vec3 Ellipsoid::get_random_point(RandomEngine &random_engine) {
    float x, y, z;
    do {
        x = base_distribution(random_engine, true);
        y = base_distribution(random_engine, true);
        z = base_distribution(random_engine, true);
    } while (x * x + y * y + z * z > 1 or x == 0 and y == 0 and z == 0);

    return glm::normalize(glm::vec3(x, y, z)) * _r;
}

float Ellipsoid::get_point_pdf(glm::vec3 P) {
    glm::vec3 N = glm::normalize(P / _r);
    glm::vec3 R2 = _r * _r;
    return 0.25f / (float) M_PI / sqrtf(N.x * N.x * R2.y * R2.z + R2.x * N.y * N.y * R2.z + R2.x * R2.y * N.z * N.z);
}

Box::Box(std::istream& in_stream) {
    _s = stream_vec3(in_stream);
}

std::pair<float, float> Box::get_ts(glm::vec3 O, glm::vec3 D) {
    glm::vec3 T1 = (_s - O) / D, T2 = (-_s - O) / D;

    float tx1 = fminf(T1.x, T2.x), tx2 = fmaxf(T1.x, T2.x);
    float ty1 = fminf(T1.y, T2.y), ty2 = fmaxf(T1.y, T2.y);
    float tz1 = fminf(T1.z, T2.z), tz2 = fmaxf(T1.z, T2.z);

    float t1 = fmaxf(tx1, fmaxf(ty1, tz1));
    float t2 = fminf(tx2, fminf(ty2, tz2));

    if (t1 > t2)
        return {-F_INF, -F_INF};
    return {t1, t2};
}

glm::vec3 Box::get_normal(glm::vec3 P) {
    glm::vec3 nP = P / _s;
    glm::vec3 anP = glm::abs(P / _s);

    if (anP.x >= anP.y and anP.x >= anP.z) {
        if (nP.x > 0)
            return {1, 0, 0};
        return {-1, 0, 0};
    }

    if (anP.y >= anP.x and anP.y >= anP.z){
        if (nP.y > 0)
            return {0, 1, 0};
        return {0, -1, 0};
    }

    if (nP.z > 0)
        return {0, 0, 1};
    return {0, 0, -1};
}

glm::vec3 Box::get_random_point(RandomEngine& random_engine) {
    float Wx = _s.y * _s.z, Wy = _s.x * _s.z, Wz = _s.y * _s.z;
    float U = base_distribution(random_engine, false) * (Wx + Wy + Wz);

    float U1 = base_distribution(random_engine, true);
    float U2 = base_distribution(random_engine, true);
    float coin = base_distribution(random_engine, true) < 0 ? -1 : 1;

    if (U < Wx)
        return {coin * _s.x, U1 * _s.y, U2 * _s.z};
    if (U < Wx + Wy)
        return {U1 * _s.x, coin * _s.y, U2 * _s.z};
    return {U1 * _s.x, U2 * _s.y, coin * _s.z};

}

float Box::get_point_pdf(glm::vec3 P) {
    return 1.f / 8 / (_s.y * _s.z + _s.x * _s.z + _s.y * _s.z);
}

Triangle::Triangle(std::istream &in_stream) {
    _a = stream_vec3(in_stream);
    _b = stream_vec3(in_stream);
    _c = stream_vec3(in_stream);
}

Triangle::Triangle(glm::vec3 a, glm::vec3 b, glm::vec3 c) {
    _a = a;
    _b = b;
    _c = c;
}

std::pair<float, float> Triangle::get_ts(glm::vec3 O, glm::vec3 D) {
    glm::mat3x3 linear_system(_b - _a, _c - _a, -D);
    glm::vec3 u_v_t = glm::inverse(linear_system) * (O - _a);

    float u = u_v_t.x, v = u_v_t.y, t = u_v_t.z;
    if (u < 0 or v < 0 or u + v > 1)
        return {-F_INF, -F_INF};
    return {-F_INF, t};
}

glm::vec3 Triangle::get_normal(glm::vec3 P) {
    return glm::normalize(glm::cross(_b - _a, _c - _a));
}

glm::vec3 Triangle::get_random_point(RandomEngine &random_engine) {
    float u = base_distribution(random_engine, false);
    float v = base_distribution(random_engine, false);

    if (u + v > 1) {
        u = 1 - u;
        v = 1 - v;
    }

    return _a + u * (_b - _a) + v * (_c - _a);
}

float Triangle::get_point_pdf(glm::vec3 P) {
    return 2 / glm::length(glm::cross(_b - _a, _c - _a));
}

Primitive::Primitive(std::istream& in_stream) {
    std::string command;

    while (in_stream >> command) {
        std::transform(command.begin(), command.end(), command.begin(), ::toupper);

        if (command == "PLANE")
            _geometry = std::unique_ptr<Geometry>(new Plane(in_stream));
        else if (command == "ELLIPSOID")
            _geometry = std::unique_ptr<Geometry>(new Ellipsoid(in_stream));
        else if (command == "BOX")
            _geometry = std::unique_ptr<Geometry>(new Box(in_stream));
        else if (command == "TRIANGLE")
            _geometry = std::unique_ptr<Geometry>(new Triangle(in_stream));
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
            in_stream >> _ior;
        else if (command == "EMISSION")
            _emission = stream_vec3(in_stream);
    }
}

Primitive::Primitive(glm::vec3 a, glm::vec3 b, glm::vec3 c) {
    _geometry = std::unique_ptr<Geometry>(new Triangle(a, b, c));
}

void Primitive::set_color(glm::vec3 color) {
    _color = color;
}

void Primitive::set_material(Material material) {
    _material = material;
}

void Primitive::set_emission(glm::vec3 emission) {
    _emission = emission;
}

bool Primitive::is_plane() {
    return _geometry->is_plane();
}

bool Primitive::has_emission() {
    return glm::length(_emission) > 0;
}

std::pair<float, float> Primitive::get_ts(glm::vec3 O, glm::vec3 D) {
    O -= _position;

    glm::quat reverse_rotation = glm::conjugate(_rotation);

    O = fast_rotate(reverse_rotation, O);
    D = fast_rotate(reverse_rotation, D);

    return _geometry->get_ts(O, D);
}

float Primitive::get_t(glm::vec3 O, glm::vec3 D) {
    float t1, t2;
    std::tie(t1, t2) = get_ts(O, D);
    return pick_t(t1, t2);
}

glm::vec3 Primitive::get_normal(glm::vec3 P) {
    glm::quat reverse_rotation = glm::conjugate(_rotation);

    P -= _position;
    P = fast_rotate(reverse_rotation, P);

    return glm::normalize(fast_rotate(_rotation, _geometry->get_normal(P)));
}

glm::vec3 Primitive::get_color(glm::vec3 O, glm::vec3 D, float t, const Scene& scene, int recursion_depth) {
    glm::vec3 P = O + t * D;
    glm::vec3 N = get_normal(P);
    bool inside = false;

    if (glm::dot(D, N) > 0) {
        N = -N;
        inside = true;
    }

    if (_material == Material::METALLIC) {
        glm::vec3 new_D = glm::normalize(D - 2.f * glm::dot(N, D) * N);

        glm::vec3 reflected_color = _emission +
                scene.get_color(P + new_D * EPSILON, new_D, recursion_depth + 1);
        reflected_color *= _color;

        check_color(reflected_color, 348765);
        return reflected_color;
    }

    if (_material == Material::DIELECTRIC) {
        float cos_t1 = -glm::dot(N, glm::normalize(D));
        float sin_t1 = sqrtf(1 - cos_t1 * cos_t1);

        float n1 = 1, n2 = _ior;
        if (inside)
            std::swap(n1, n2);

        float sin_t2 = sin_t1 * n1 / n2;

        float R0 = powf((n1 - n2) / (n1 + n2), 2);
        float R = R0 + (1 - R0) * powf(1 - cos_t1, 5);

        std::uniform_real_distribution<float> distribution(0, 1);
        float U = distribution(*scene.random_engine(omp_get_thread_num()));

        if (sin_t2 > 1 or U < R) {
            glm::vec3 reflected_D = glm::normalize(D - 2.f * glm::dot(N, D) * N);
            glm::vec3 reflected_color = _emission +
                    scene.get_color(P + reflected_D * EPSILON, reflected_D, recursion_depth + 1);

            check_color(reflected_color, 397234);
            return reflected_color;
        }

        float cos_t2 = sqrtf(1 - sin_t2 * sin_t2);

        glm::vec3 refracted_D = glm::normalize(D * n1 / n2 + N * (cos_t1 * n1 / n2 - cos_t2));
        glm::vec3 refracted_color = _emission +
                scene.get_color(P + refracted_D * EPSILON, refracted_D, recursion_depth + 1);

        if (not inside)
            refracted_color *= _color;

        check_color(refracted_color, 387463);
        return refracted_color;
    }

    // if (_material == Material::DIFFUSIVE)

    Distribution* distribution = scene.distribution();
    glm::vec3 new_D = distribution->sample(P, N, *scene.random_engine(omp_get_thread_num()));

    if (glm::dot(new_D, N) < 0) {
        check_color(_emission, 873462);
        return _emission;
    }

    float pdf = distribution->pdf(P + new_D * EPSILON, N, new_D);
    if (pdf == 0 or pdf == F_INF) {
        check_color(_emission, 873463);
        return _emission;
    }

    glm::vec3 L_in = scene.get_color(P + new_D * EPSILON, new_D, recursion_depth + 1);

    glm::vec3 color = _emission + _color * L_in * fmaxf(0, glm::dot(new_D, N)) / (float) M_PI / pdf;

    if (std::isnan(color.x))
        color.x = 0;
    if (std::isnan(color.y))
        color.y = 0;
    if (std::isnan(color.z))
        color.z = 0;

    check_color(color, 873464);
    return color;
}

glm::vec3 Primitive::get_random_point(RandomEngine& random_engine) {
    glm::vec3 P = _geometry->get_random_point(random_engine);

    P = fast_rotate(_rotation, P);
    P += _position;

    return P;
}

float Primitive::get_point_pdf(glm::vec3 P) {
    glm::quat reverse_rotation = glm::conjugate(_rotation);

    P -= _position;
    P = fast_rotate(reverse_rotation, P);

    return _geometry->get_point_pdf(P);
}

template<typename T>
glm::vec3 get_vec3_from_array(rapidjson::GenericArray<true, T> array) {
    if (array.Size() != 3)
        exit(983573450);

    return {array[0].GetFloat(), array[1].GetFloat(), array[2].GetFloat()};
}

template<typename T>
glm::vec4 get_vec4_from_array(rapidjson::GenericArray<true, T> array) {
    if (array.Size() != 4)
        exit(983573451);

    return {array[0].GetFloat(), array[1].GetFloat(), array[2].GetFloat(), array[3].GetFloat()};
}

template<typename T>
glm::quat get_quat_from_array(rapidjson::GenericArray<true, T> array) {
    if (array.Size() != 4)
        exit(983573452);

    return {array[3].GetFloat(), array[0].GetFloat(), array[1].GetFloat(), array[2].GetFloat()};
}

template<typename T>
glm::mat4 get_mat4_from_array(rapidjson::GenericArray<true, T> array) {
    if (array.Size() != 16)
        exit(983573453);

    return {
            array[0].GetFloat(), array[1].GetFloat(), array[2].GetFloat(), array[3].GetFloat(),
            array[4].GetFloat(), array[5].GetFloat(), array[6].GetFloat(), array[7].GetFloat(),
            array[8].GetFloat(), array[9].GetFloat(), array[10].GetFloat(), array[11].GetFloat(),
            array[12].GetFloat(), array[13].GetFloat(), array[14].GetFloat(), array[15].GetFloat(),
    };
}

template<typename T>
std::vector<unsigned int> get_numbers_from_buffer(unsigned char* buffer, int count) {
    auto* typed_buffer = reinterpret_cast<T*>(buffer);
    std::vector<unsigned int> numbers;

    for (int i = 0; i < count; ++i)
        numbers.push_back((unsigned int) typed_buffer[i]);

    return numbers;
}

std::vector<glm::vec3> get_vec3s_from_buffer(unsigned char* buffer, int count) {
    auto* typed_buffer = reinterpret_cast<float*>(buffer);
    std::vector<glm::vec3> vec3s;

    for (int i = 0; i < 3 * count; i += 3)
        vec3s.emplace_back(typed_buffer[i], typed_buffer[i + 1], typed_buffer[i + 2]);

    return vec3s;
}

void Scene::initialize_node(const rapidjson::Document& document, int node_num, glm::mat4 current_transform) {
    if (document["nodes"][node_num].HasMember("matrix"))
        current_transform *= get_mat4_from_array(document["nodes"][node_num]["matrix"].GetArray());
    else {
        glm::mat4 additional_transform(1.0f);

        if (document["nodes"][node_num].HasMember("scale")) {
            glm::mat4 scale_matrix =
                    glm::scale(glm::mat4(1.0f),
                               get_vec3_from_array(document["nodes"][node_num]["scale"].GetArray()));
            additional_transform = scale_matrix * additional_transform;
        }

        if (document["nodes"][node_num].HasMember("rotation")) {
            glm::mat4 rotation_matrix =
                    glm::mat4_cast(get_quat_from_array(document["nodes"][node_num]["rotation"].GetArray()));
            additional_transform = rotation_matrix * additional_transform;
        }

        if (document["nodes"][node_num].HasMember("translation")) {
            glm::mat4 translation_matrix =
                    glm::translate(glm::mat4(1.0f),
                                   get_vec3_from_array(document["nodes"][node_num]["translation"].GetArray()));
            additional_transform = translation_matrix * additional_transform;
        }

        current_transform *= additional_transform;
    }

    if (not _camera_initialized and document["nodes"][node_num].HasMember("camera")) {
        int camera_num = document["nodes"][node_num]["camera"].GetInt();

        if (document["cameras"][camera_num]["type"] == "perspective") {
            _camera_fov_y = document["cameras"][camera_num]["perspective"]["yfov"].GetFloat();

            glm::vec3 new_camera_position = glm::vec3(current_transform * glm::vec4(_camera_position, 1));

            _camera_right = glm::normalize(glm::vec3(current_transform * glm::vec4(_camera_right + _camera_position, 1)) - new_camera_position);
            _camera_up = glm::normalize(glm::vec3(current_transform * glm::vec4(_camera_up + _camera_position, 1)) - new_camera_position);
            _camera_forward = glm::normalize(glm::vec3(current_transform * glm::vec4(_camera_forward + _camera_position, 1)) - new_camera_position);

            _camera_position = new_camera_position;

            _camera_initialized = true;
        }
    }

    if (document["nodes"][node_num].HasMember("mesh")) {
        int mesh_num = document["nodes"][node_num]["mesh"].GetInt();
        auto mesh_primitives = document["meshes"][mesh_num]["primitives"].GetArray();

        for (int primitive_num = 0; primitive_num < mesh_primitives.Size(); ++primitive_num) {
            int indices_accessor_num = mesh_primitives[primitive_num]["indices"].GetInt();
            int indices_buffer_view_num = document["accessors"][indices_accessor_num]["bufferView"].GetInt();

            int indices_buffer_num = document["bufferViews"][indices_buffer_view_num]["buffer"].GetInt();
            int indices_buffer_offset = document["bufferViews"][indices_buffer_view_num]["byteOffset"].GetInt() +
                                        document["accessors"][indices_accessor_num]["byteOffset"].GetInt();

            unsigned char* indices_buffer = _buffers[indices_buffer_num].data() + indices_buffer_offset;

            int indices_count = document["accessors"][indices_accessor_num]["count"].GetInt();

            std::vector<unsigned int> indices;
            if (document["accessors"][indices_accessor_num]["componentType"].GetInt() == 5123)
                indices = get_numbers_from_buffer<unsigned short>(indices_buffer, indices_count);
            else if (document["accessors"][indices_accessor_num]["componentType"].GetInt() == 5125)
                indices = get_numbers_from_buffer<unsigned int>(indices_buffer, indices_count);
            else
                exit(937423984);

            int positions_accessor_num = mesh_primitives[primitive_num]["attributes"]["POSITION"].GetInt();
            int positions_buffer_view_num = document["accessors"][positions_accessor_num]["bufferView"].GetInt();

            int positions_buffer_num = document["bufferViews"][positions_buffer_view_num]["buffer"].GetInt();
            int positions_buffer_offset = document["bufferViews"][positions_buffer_view_num]["byteOffset"].GetInt() +
                                          document["accessors"][positions_accessor_num]["byteOffset"].GetInt();

            unsigned char* positions_buffer = _buffers[positions_buffer_num].data() + positions_buffer_offset;

            int positions_count = document["accessors"][positions_accessor_num]["count"].GetInt();

            if (document["accessors"][positions_accessor_num]["componentType"].GetInt() != 5126)
                exit(937423985);
            std::vector<glm::vec3> positions = get_vec3s_from_buffer(positions_buffer, positions_count);

            // Adding material properties.

            glm::vec4 color(1, 1, 1, 1);
            Material material = Material::METALLIC;
            glm::vec3 emission(0);

            if (mesh_primitives[primitive_num].HasMember("material")) {
                int material_num = mesh_primitives[primitive_num]["material"].GetInt();

                if (document["materials"][material_num]["pbrMetallicRoughness"].HasMember("baseColorFactor"))
                    color = get_vec4_from_array(
                            document["materials"][material_num]["pbrMetallicRoughness"]["baseColorFactor"].GetArray()
                    );

                if (color.a < 1)
                    material = Material::DIELECTRIC;
                else if (document["materials"][material_num]["pbrMetallicRoughness"].HasMember("metallicFactor")
                         and
                         document["materials"][material_num]["pbrMetallicRoughness"]["metallicFactor"].GetFloat() == 0)
                    material = Material::DIFFUSIVE;

                if (document["materials"][material_num].HasMember("emissiveFactor"))
                    emission = get_vec3_from_array(
                            document["materials"][material_num]["emissiveFactor"].GetArray()
                    );

                if (document["extensionsUsed"].HasMember("KHR_materials_emissive_strength") and
                    document["materials"][material_num]["extensions"]["KHR_materials_emissive_strength"]
                    .HasMember("emissiveStrength"))
                    emission *= document["materials"][material_num]
                            ["extensions"]["KHR_materials_emissive_strength"]["emissiveStrength"].GetFloat();
            }

            for (int i = 0; i < indices.size(); i += 3) {
                glm::vec3 a = glm::vec3(current_transform * glm::vec4(positions[indices[i]], 1));
                glm::vec3 b = glm::vec3(current_transform * glm::vec4(positions[indices[i + 1]], 1));
                glm::vec3 c = glm::vec3(current_transform * glm::vec4(positions[indices[i + 2]], 1));

                auto primitive = std::make_unique<Primitive>(a, b, c);
                primitive->set_color(glm::vec3(color));
                primitive->set_material(material);
                primitive->set_emission(emission);

                if (primitive->has_emission())
                    _light_sources.push_back(
                        std::move(std::unique_ptr<Distribution>(new LightSource(*primitive))));

                _primitives.push_back(std::move(primitive));
            }
        }
    }

    if (document["nodes"][node_num].HasMember("children")) {
        auto children_array = document["nodes"][node_num]["children"].GetArray();
        for (int i = 0; i < children_array.Size(); ++i)
            initialize_node(document, document["nodes"][node_num]["children"][i].GetInt(), current_transform);
    }
}

void Scene::initialize_node(const rapidjson::Document& document, int node_num) {
    initialize_node(document, node_num, glm::mat4(1.0f));
}

Scene::Scene(const rapidjson::Document& document, int width, int height) {
    _dimension_width = width;
    _dimension_height = height;

    int max_thead_num = omp_get_max_threads();
    for (int i = 0; i < max_thead_num; ++i)
        _random_engines.push_back(std::make_unique<RandomEngine>(i));

    auto buffers_array = document["buffers"].GetArray();
    for (int i = 0; i < buffers_array.Size(); ++i) {
        auto buffer_file = buffers_array[i]["uri"].GetString();

        std::ifstream in_stream(buffer_file, std::ios::binary);
        std::vector<unsigned char> buffer(std::istreambuf_iterator<char>(in_stream), {});
        _buffers.push_back(buffer);
    }

    if (document.HasMember("scene")) {
        int scene_num = document["scene"].GetInt();
        auto nodes_array = document["scenes"][scene_num]["nodes"].GetArray();

        for (int i = 0; i < nodes_array.Size(); ++i) {
            int node_num = nodes_array[i].GetInt();
            initialize_node(document, node_num);
        }
    }
    else {
        auto nodes_array = document["nodes"].GetArray();

        for (int i = 0; i < nodes_array.Size(); ++i) {
            initialize_node(document, i);
        }
    }

    if (_light_sources.empty())
        _distribution = std::unique_ptr<Distribution>(new CosineHemisphere());
    else if (_light_sources.size() == 1) {
        std::vector<std::unique_ptr<Distribution>> distributions;

        distributions.push_back(std::unique_ptr<Distribution>(new CosineHemisphere()));
        distributions.push_back(std::move(_light_sources[0]));

        _distribution = std::unique_ptr<Distribution>(new Mix(std::move(distributions)));
    }
    else {
        std::vector<std::unique_ptr<Distribution>> distributions;

        distributions.push_back(std::unique_ptr<Distribution>(new CosineHemisphere()));
        distributions.push_back(std::unique_ptr<Distribution>(new Mix(std::move(_light_sources))));

        _distribution = std::unique_ptr<Distribution>(new Mix(std::move(distributions)));
    }
}

/*
Scene::Scene(std::ifstream& in_stream) {
    int max_thead_num = omp_get_max_threads();
    for (int i = 0; i < max_thead_num; ++i)
        _random_engines.push_back(std::make_unique<RandomEngine>(i));

    std::string command;
    std::stringstream primitive_stream;

    bool new_primitive = false;

    std::vector<std::unique_ptr<Distribution>> light_sources;

    while (in_stream >> command) {
        std::transform(command.begin(), command.end(), command.begin(), ::toupper);

        bool not_primitive_or_light = false;

        if (command == "DIMENSIONS") {
            not_primitive_or_light = true;
            in_stream >> _dimension_width >> _dimension_height;
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
            in_stream >> _camera_fov_x;
        }
        else if (command == "RAY_DEPTH") {
            not_primitive_or_light = true;
            in_stream >> _ray_depth;
        }
        else if (command == "SAMPLES") {
            not_primitive_or_light = true;
            in_stream >> _sample_num;
        }

        if (not_primitive_or_light || command == "NEW_PRIMITIVE") {
            if (new_primitive) {
                new_primitive = false;

                std::unique_ptr<Primitive> primitive = std::make_unique<Primitive>(primitive_stream);

                if (not primitive->is_plane() and primitive->has_emission())
                    light_sources.push_back(
                            std::move(std::unique_ptr<Distribution>(new LightSource(*primitive))));

                _primitives.push_back(std::move(primitive));
            }
        }
        if (not_primitive_or_light)
            continue;

        if (command == "NEW_PRIMITIVE") {
            new_primitive = true;
            primitive_stream = std::stringstream();
        }
        else if (new_primitive)
            primitive_stream << command << " ";
    }

    if (new_primitive) {
        auto primitive = std::make_unique<Primitive>(primitive_stream);

        if (not primitive->is_plane() and primitive->has_emission())
            light_sources.push_back(std::move(std::unique_ptr<Distribution>(new LightSource(*primitive))));

        _primitives.push_back(std::move(primitive));
    }

    _distribution = std::unique_ptr<Distribution>(new CosineHemisphere());

    if (light_sources.empty())
        _distribution = std::unique_ptr<Distribution>(new CosineHemisphere());
    else if (light_sources.size() == 1) {
        std::vector<std::unique_ptr<Distribution>> distributions;

        distributions.push_back(std::unique_ptr<Distribution>(new CosineHemisphere()));
        distributions.push_back(std::move(light_sources[0]));

        _distribution = std::unique_ptr<Distribution>(new Mix(std::move(distributions)));
    }
    else {
        std::vector<std::unique_ptr<Distribution>> distributions;

        distributions.push_back(std::unique_ptr<Distribution>(new CosineHemisphere()));
        distributions.push_back(std::unique_ptr<Distribution>(new Mix(std::move(light_sources))));

        _distribution = std::unique_ptr<Distribution>(new Mix(std::move(distributions)));
    }
}
*/

std::pair<float, Primitive*> Scene::get_t(glm::vec3 O, glm::vec3 D) const {
    float t = -1;
    Primitive* curr_primitive = nullptr;

    for (auto& primitive : _primitives) {
        float new_t = primitive->get_t(O, D);
        if (t <= 0 or (new_t > 0 and new_t < t)) {
            t = new_t;
            curr_primitive = primitive.get();
        }
    }

    return {t, curr_primitive};
}

glm::vec3 Scene::get_color(glm::vec3 O, glm::vec3 D, int recursion_depth) const {
    if (recursion_depth > _ray_depth)
        return glm::vec3(0);

    D = glm::normalize(D);

    float t;
    Primitive* curr_primitive;
    std::tie(t, curr_primitive) = get_t(O, D);

    glm::vec3 curr_color;
    if (t <= 0 || curr_primitive == nullptr)
        curr_color = _bg_color;
    else
        curr_color = curr_primitive->get_color(O, D, t, *this, recursion_depth);

    check_color(curr_color, 937433);
    return curr_color;
}

glm::vec3 Scene::get_color(glm::vec3 O, glm::vec3 D) const {
    return get_color(O, D, 1);
}

std::vector<uint8_t> Scene::render(int samples) const {
    std::vector<uint8_t> render(3 * _dimension_width * _dimension_height, 0);

    float tan_half_fov_y = tanf(_camera_fov_y / 2);
    float tan_half_fov_x = tan_half_fov_y / (float) _dimension_height * (float) _dimension_width;

    std::uniform_real_distribution<float> distribution(0, 1);

    #pragma omp parallel for collapse(2) default(none) shared(samples, distribution, tan_half_fov_x, tan_half_fov_y, render)
    for (int p_x = 0; p_x < _dimension_width; ++p_x) {
        for (int p_y = 0; p_y < _dimension_height; ++p_y) {
            auto color_sum = glm::vec3(0);

            for (int curr_sample = 0; curr_sample < samples; ++curr_sample) {
                auto curr_random_engine = random_engine(omp_get_thread_num());
                float u_x = distribution(*curr_random_engine);
                float u_y = distribution(*curr_random_engine);

                float c_x = (2 * ((float) p_x + u_x) / (float) _dimension_width - 1) * tan_half_fov_x;
                float c_y = -1 * (2 * ((float) p_y + u_y) / (float) _dimension_height - 1) * tan_half_fov_y;
                float c_z = 1;

                glm::vec3 D = glm::normalize(c_x * _camera_right + c_y * _camera_up + c_z * _camera_forward);
                glm::vec3 O = _camera_position;

                color_sum += get_color(O, D);
            }

            color_sum /= samples;
            check_color(color_sum, 384763);

            color_sum = aces_tonemap(color_sum);

            int it = 3 * (p_x + p_y * _dimension_width);
            render[it + 0] = fix_color(color_sum.r);
            render[it + 1] = fix_color(color_sum.g);
            render[it + 2] = fix_color(color_sum.b);
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

RandomEngine* Scene::random_engine(int thread) const {
    return _random_engines[thread].get();
}

Distribution* Scene::distribution() const {
    return _distribution.get();
}
