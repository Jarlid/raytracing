#include <iostream>

int main(int argc, char* argv[]) {
    if (argc != 3)
        exit(1);

    char* in_file = argv[1];
    char* out_file = argv[2];

    std::cout << in_file << std::endl << out_file << std::endl;
    return 0;
}
