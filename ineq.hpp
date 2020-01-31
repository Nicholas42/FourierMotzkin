#ifndef INEQ_HPP
#define INEQ_HPP
#include <vector>
#include <array>
#include <functional>
#include <istream>
#include <sstream>
#include <exception>
#include <cassert>
#include <numeric>
#include <limits>
#include <algorithm>

namespace FourierMotzkin
{
using value_t = double;

enum class Sign
{
    Zero,
    Positive,
    Negative
};

class InequalitySystem
{
    class Inequality
    {
        friend InequalitySystem;
        using cmp = std::less_equal<value_t>;

      public:
        Inequality(size_t num_vars, value_t rhs) : _coeffs(num_vars, value_t{}), _rhs(rhs)
        {}

        Inequality(const Inequality &ineq, size_t index, size_t without);

        Inequality(const Inequality &ineq1, size_t index1, const Inequality &ineq2,
                   size_t index2, size_t without);

        bool is_valid(const std::vector<value_t> &vars = {}) const;

        value_t evaluate_lhs(const std::vector<value_t> &vars = {}) const;

        void scale(value_t scalar);

        Sign get_sign(size_t index = 0) const;

        //divide by the absolute value of the coefficient at the given index.
        void normalize_on(size_t index = 0);

        size_t num_vars() const
        {
            return _coeffs.size();
        }

        value_t rhs() const
        {
            return _rhs;
        }

        std::vector<size_t> get_parents() const
        {
            return _parents;
        }

      private:
        std::vector<value_t> _coeffs;
        value_t _rhs;

        /*Saves the inequalities of the last iteration
        *from which this inequality was constructed.
        */
        std::vector<size_t> _parents;

        /* The factor by which the original inequality is multiplied is needed to
        * construct a solution vector to the original problem in the backtracking process.
        */
        value_t _scaling_factor = 1;
    };

  public:
    InequalitySystem(size_t num_vars, size_t num_ineqs) :
            _num_vars(num_vars),
            _num_ineqs(num_ineqs)
    {
        _ineqs.reserve(num_ineqs);
    }

    //Construct a new inequality system by elimination of the last variable
    InequalitySystem reduce_on(size_t index = 0);

    bool is_valid(const std::vector<value_t> &vars = {}) const;

    std::array<std::vector<size_t>, 3> partition(size_t index) const;

    value_t get_max(const std::vector<size_t> &to_eval,
                    const std::vector<value_t> &vars = {}) const;

    value_t get_min(const std::vector<size_t> &to_eval,
                    const std::vector<value_t> &vars = {}) const;

    size_t num_vars() const
    {
        return _num_vars;
    }

    size_t num_ineqs() const
    {
        return _num_ineqs;
    }

    size_t find_invalid(const std::vector<value_t> &vars = {}) const;

    std::vector<size_t> get_parents(size_t index) const;

    value_t calc_variable(size_t index, const std::vector<value_t> &known_vars) const;

    bool check_counterexample(const std::vector<value_t> &counterexample) const;

    value_t get_scaling_factor(size_t index) const
    {
        return _ineqs[index]._scaling_factor;
    }

    Inequality read_ineq(std::istream &in, value_t rhs);

    friend std::istream &operator>>(std::istream &in, InequalitySystem &sys);

  private:
    std::vector<Inequality> _ineqs;
    size_t _num_vars;
    size_t _num_ineqs;
};

}   // END NAMESPACE FourierMotzkin
#endif
