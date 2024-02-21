#include <iostream>
#include <fstream>

#include "geometry/geometry.h"

int main(int argc, char* argv[]) {
    if (argc != 3)
        exit(1);

    char* in_file = argv[1];
    char* out_file = argv[2];

    std::cout << in_file << std::endl << out_file << std::endl;

    std::ifstream in_stream(in_file, std::ios_base::in);
    Scene scene(&in_stream);
    in_stream.close();

    std::ofstream out_stream(out_file, std::ios_base::out);

    out_stream << "P6" << std::endl;
    out_stream << scene.width() << " " << scene.height() << std::endl;
    out_stream << 255 << std::endl;
    for (auto elem : scene.render())
        out_stream << elem;

    out_stream.close();

    return 0;
}
