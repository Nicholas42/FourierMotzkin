#include <iostream>
#include "fourier_motzkin.hpp"

int main(int argc, char **argv)
{
    if (argc != 2)
    {
        std::cout << "Usage: " << argv[0] << " FILENAME\n";
        return 1;
    }

    std::cout << FourierMotzkin::fourier_motzkin(argv[1]);
    return 0;
}
