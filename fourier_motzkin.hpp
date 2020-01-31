#ifndef FOURIER_MOTZKIN_HPP
#define FOURIER_MOTZKIN_HPP

#include <vector>
#include <fstream>
#include <iterator>
#include "ineq.hpp"

namespace FourierMotzkin
{
struct certificate
{
    //Solution or certificate vector
    std::vector<value_t> vec;

    bool valid;

    friend std::ostream &operator<<(std::ostream &out, const certificate &cert);
};

certificate fourier_motzkin(const std::string &filename);
}   // END NAMESPACE FourierMotzkin
#endif   // FOURIER_MOTZKIN_HPP
