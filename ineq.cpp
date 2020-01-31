#include "ineq.hpp"
#include <iostream>

namespace FourierMotzkin
{
    /*Construct Inequality from another one without a given variable. Save the old inequality in _parents.
     * This function is needed to eliminate variables with coeffizients of 0.
    */
InequalitySystem::Inequality::Inequality(const Inequality &ineq, size_t index,
                                         size_t without) :
        _coeffs(),
        _rhs(ineq._rhs),
        _parents{index}
{
    _coeffs.reserve(ineq.num_vars() - 1);

    for (size_t i = 0; i < ineq._coeffs.size(); ++i)
    {
        if (i != without)
        {
            _coeffs.push_back(ineq._coeffs[i]);
        }
    }
}

/*Add two inequalities and save these two as parents.
 * This function is used when the absolute value of the last variable is 1 and the signs are different.
 * Therefore the last variable is eliminated.
*/
InequalitySystem::Inequality::Inequality(const Inequality &ineq1, size_t index1,
                                         const Inequality &ineq2, size_t index2,
                                         size_t without) :
        _coeffs(),
        _rhs(ineq1._rhs + ineq2._rhs),
        _parents{index1, index2}
{
    _coeffs.reserve(ineq1.num_vars() - 1);

    for (size_t i = 0; i < ineq1.num_vars(); ++i)
    {
        if (i != without)
        {
            _coeffs.push_back(ineq1._coeffs[i] + ineq2._coeffs[i]);
        }
    }
}

bool InequalitySystem::Inequality::is_valid(const std::vector<value_t> &vars) const
{
    return cmp()(evaluate_lhs(vars), _rhs);
}

value_t InequalitySystem::Inequality::evaluate_lhs(const std::vector<value_t> &vars) const
{
    assert(vars.size() == num_vars() && "Number of variables does not match inequality.");
    return std::inner_product(_coeffs.begin(), _coeffs.end(), vars.begin(), value_t{0});
}

void InequalitySystem::Inequality::scale(value_t scalar)
{
    assert(scalar > 0 && "Only positive scaling is allowed.");
    std::transform(_coeffs.begin(),
                   _coeffs.end(),
                   _coeffs.begin(),
                   [scalar](value_t coeff) { return coeff * scalar; });

    _rhs *= scalar;
    _scaling_factor /= scalar;
}

InequalitySystem::Inequality InequalitySystem::read_ineq(std::istream &in, value_t rhs)
{
    Inequality ineq{_num_vars, rhs};
    size_t i = 0;
    for (; i < ineq.num_vars() && in >> ineq._coeffs[i]; ++i)
    {
    }

    if (i < ineq.num_vars())
    {
        std::stringstream ss;
        ss << "Not enough coefficients in instream. Expected " << ineq.num_vars();
        ss << ", received " << i << ".\n";
        throw std::runtime_error(ss.str());
    }

    return ineq;
}

Sign InequalitySystem::Inequality::get_sign(size_t index) const
{
    if (_coeffs[index] == 0)
    {
        return Sign::Zero;
    }
    if (_coeffs[index] > 0)
    {
        return Sign::Positive;
    }

    return Sign::Negative;
}

void InequalitySystem::Inequality::normalize_on(size_t index)
{
    scale(value_t{1} / std::abs(_coeffs[index]));
}

//We need the signs of the given variable of all inequalities to combine them.
std::array<std::vector<size_t>, 3> InequalitySystem::partition(size_t index) const
{
    std::array<std::vector<size_t>, 3> part;

    for (size_t i = 0; i < _ineqs.size(); ++i)
    {
        part[(size_t) _ineqs[i].get_sign(index)].push_back(i);
    }

    return part;
}

//Eliminate last variable by constructing a new inequality system.
InequalitySystem InequalitySystem::reduce_on(size_t index)
{
    assert(_num_vars > 0);
    auto part = partition(index);

    //The positive and negative case need to have the same absolute value to add up to 0.
    for (auto sign: {Sign::Positive, Sign::Negative})
    {
        for (auto i: part[(size_t) sign])
        {
            _ineqs[i].normalize_on(index);
        }
    }

    const auto num_ineqs =
        part[(size_t) Sign::Zero].size() +
        part[(size_t) Sign::Positive].size() * part[(size_t) Sign::Negative].size();

    InequalitySystem ret{_num_vars - 1, num_ineqs};

    for (const auto pos: part[(size_t) Sign::Positive])
    {
        for (const auto neg: part[(size_t) Sign::Negative])
        {
            ret._ineqs.emplace_back(_ineqs[pos], pos, _ineqs[neg], neg, index);
            assert(ret._ineqs.back().num_vars() == ret.num_vars());
        }
    }

    for (auto zero: part[(size_t) Sign::Zero])
    {
        ret._ineqs.emplace_back(_ineqs[zero], zero, index);
        assert(ret._ineqs.back().num_vars() == ret.num_vars());
    }

    return ret;
}

bool InequalitySystem::is_valid(const std::vector<value_t> &vars) const
{
    return std::all_of(_ineqs.begin(), _ineqs.end(), [&vars](const auto &ineq) {
        return ineq.is_valid(vars);
    });
}


value_t InequalitySystem::get_max(const std::vector<size_t> &to_eval,
                                  const std::vector<value_t> &vars) const
{
    auto max = std::numeric_limits<value_t>::lowest();

    for (auto i: to_eval)
    {
        max = std::max(max, _ineqs[i].evaluate_lhs(vars) - _ineqs[i].rhs());
    }

    return max;
}

value_t InequalitySystem::get_min(const std::vector<size_t> &to_eval,
                                  const std::vector<value_t> &vars) const
{
    auto min = std::numeric_limits<value_t>::max();

    for (auto i: to_eval)
    {
        min = std::min(min, _ineqs[i].rhs() - _ineqs[i].evaluate_lhs(vars));
    }

    return min;
}

//Read an inequality system.
std::istream &operator>>(std::istream &in, InequalitySystem &sys)
{
    std::string line;
    std::getline(in, line);
    std::stringstream line_wise(line);

    std::vector<value_t> b;

    value_t tmp_val;
    while (line_wise >> tmp_val)
    {
        b.push_back(tmp_val);
    }

    if (b.size() != sys._num_ineqs)
    {
        std::stringstream ss;
        ss << "Vector b has wrong size in input. Expected " << sys._num_ineqs;
        ss << ", received " << b.size() << ".\n";
        throw std::runtime_error(ss.str());
    }

    size_t cur_ind = 0;
    while (std::getline(in, line))
    {
        line_wise.clear();
        line_wise.str(line);

        sys._ineqs.emplace_back(sys.read_ineq(line_wise, b[cur_ind++]));
    }

    return in;
}

value_t InequalitySystem::calc_variable(size_t index,
                                        const std::vector<value_t> &known_vars) const
{
    auto part = partition(index);

    if (part[(size_t) Sign::Positive].size() + part[(size_t) Sign::Negative].size() == 0)
    {
        // Everything is feasible, so we take 0 to avoid overflows.
        return 0;
    }

    //One of the sizes is non-zero, therefore the calculated variable is not infinity.
    auto ret = part[(size_t) Sign::Positive].size() > part[(size_t) Sign::Negative].size()
                   ? get_min(part[(size_t) Sign::Positive], known_vars)
                   : get_max(part[(size_t) Sign::Negative], known_vars);

    return ret;
}

size_t InequalitySystem::find_invalid(const std::vector<value_t> &vars) const
{
    auto it = std::find_if(_ineqs.begin(), _ineqs.end(), [&vars](const Inequality &ineq) {
        return !ineq.is_valid(vars);
    });

    assert(it != _ineqs.end());
    return it - _ineqs.begin();
}

std::vector<size_t> InequalitySystem::get_parents(size_t index) const
{
    return _ineqs[index].get_parents();
}

//This function is just for checking the result.
bool InequalitySystem::check_counterexample(
    const std::vector<value_t> &counterexample) const
{
    assert(counterexample.size() == _ineqs.size());
    value_t inner_prod = std::inner_product(
        _ineqs.begin(),
        _ineqs.end(),
        counterexample.begin(),
        value_t{0},
        std::plus<>(),
        [](const auto &ineq, const auto val) { return ineq.rhs() * val; });

    bool matrix_is_zero = true;

    for (size_t j = 0; j < _num_vars && matrix_is_zero; ++j)
    {
        matrix_is_zero &= (std::inner_product(_ineqs.begin(),
                                              _ineqs.end(),
                                              counterexample.begin(),
                                              value_t{0},
                                              std::plus<>(),
                                              [j](const auto &ineq, const auto val) {
                                                  return ineq._coeffs[j] * val;
                                              }) == 0);
    }
    return (inner_prod < 0) && matrix_is_zero;
}
}   // END NAMESPACE FourierMotzkin
