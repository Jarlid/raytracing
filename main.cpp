#include <iostream>
#include <fstream>

#include <filesystem>

#include <rapidjson/include/rapidjson/istreamwrapper.h>
#include <rapidjson/include/rapidjson/document.h>

#include "geometry/scene.h"

int main(int argc, char* argv[]) {
    if (argc != 6)
        exit(238748323);

    char* in_file = argv[1];
    char* out_file = argv[5];

    std::cout << in_file << std::endl << out_file << std::endl;

    std::ifstream in_stream(in_file);

    rapidjson::IStreamWrapper wrapped_in_stream { in_stream };
    rapidjson::Document document;
    document.ParseStream( wrapped_in_stream );

    std::filesystem::path old_path = std::filesystem::current_path();
    std::filesystem::path new_path(in_file);
    if (new_path.has_parent_path())
        std::filesystem::current_path(new_path.parent_path());

    Scene scene(document, std::stoi(argv[2]), std::stoi(argv[3]));

    std::filesystem::current_path(old_path);

    std::ofstream out_stream(out_file, std::ios_base::binary);

    out_stream << "P6" << std::endl;
    out_stream << scene.width() << " " << scene.height() << std::endl;
    out_stream << 255 << std::endl;
    out_stream.write((char*) scene.render(std::stoi(argv[4])).data(),
                     3 * scene.width() * scene.height());

    out_stream.close();

    return 0;
}
