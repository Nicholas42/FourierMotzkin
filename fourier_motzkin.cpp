#include "fourier_motzkin.hpp"

namespace FourierMotzkin
{
InequalitySystem read_file(const std::string &filename)
{
    std::ifstream file(filename);
    std::string line;

    size_t rows, columns;

    if (std::getline(file, line))
    {
        std::stringstream ss(line);
        ss >> rows >> columns;
    }
    else
    {
        throw std::runtime_error("Invalid file format.");
    }

    //The vector c can be ignored.
    file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    InequalitySystem ret{columns, rows};

    file >> ret;

    return ret;
}

//Eliminate all variables of a given inequality system.
std::vector<InequalitySystem> reduce_system(InequalitySystem &&sys)
{
    std::vector<InequalitySystem> ret;

    ret.reserve(sys.num_vars());
    ret.push_back(std::move(sys));

    for (size_t i = ret.front().num_vars(); i > 0; --i)
    {
        ret.push_back(ret.back().reduce_on(i - 1));
    }

    return ret;
}

/* Given is an inequality system and a linear combination to obtain a counterexample.
 * Calculate a linear combination for the inequality system with one more variable.
*/
std::vector<value_t> compute_parents(const InequalitySystem &sys,
                                     std::vector<value_t> &&parents, size_t length)
{
    std::vector<value_t> ret(length);
    for (size_t p = 0; p < parents.size(); ++p)
    {
        for (auto np: sys.get_parents(p))
        {
            //We need to divide by the scaling factor to get the linear combination for the original problem.
            ret[np] += parents[p] / sys.get_scaling_factor(p);
        }
    }

    return ret;
}

//Compute a vector to prove that a system is infeasible (Farka's Lemma).
std::vector<value_t> compute_counter_example(const std::vector<InequalitySystem> &steps)
{
    std::vector<value_t> ret(steps.front().num_vars(), value_t{0});

    std::vector<value_t> parents(steps.back().num_ineqs());

    //Find an infeasible inequality in the last step. The linear combination is just one inequality in this step.
    parents[steps.back().find_invalid()]++;

    for (auto it = steps.rbegin(); std::next(it) != steps.rend(); ++it)
    {
        parents = compute_parents(*it, std::move(parents), std::next(it)->num_ineqs());
    }

    return parents;
}

//Find a possible solution to the original problem.
std::vector<value_t> recover_variables(const std::vector<InequalitySystem> &steps)
{
    auto num_vars = steps.front().num_vars();

    std::vector<value_t> ret;
    ret.reserve(num_vars);

    for (size_t i = steps.size() - 1; i > 0; --i)
    {
        ret.push_back(0);
        //Calculate a possible last variable.
        ret.back() = steps[i - 1].calc_variable(num_vars - i, ret);
    }

    assert(steps.front().is_valid(ret));

    return ret;
}

//Return either a possible solution or a counter example using Farka's Lemma to a given instance.
certificate fourier_motzkin(const std::string &filename)
{
    std::vector<InequalitySystem> steps = reduce_system(read_file(filename));

    assert(steps.back().num_vars() == 0);

    if (!steps.back().is_valid())
    {
        return {compute_counter_example(steps), false};
    }
    else
    {
        return {recover_variables(steps), true};
    }
}

/*Print the output in the specified format:
 * no solution: "empty [certificate vector]
 * solution: [solution vector]
*/
std::ostream &operator<<(std::ostream &out, const certificate &cert)
{
    if (!cert.valid)
    {
        out << "empty ";
    }

    std::copy(cert.vec.begin(), cert.vec.end(), std::ostream_iterator<value_t>(out, " "));
    out << std::endl;

    return out;
}

}   // END NAMESPACE FourierMotzkin
