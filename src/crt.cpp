#include "crt.hpp"
#include "math.hpp"

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <utility>

namespace factoring {
namespace {

/** @brief State needed to generate the next refined x value. */
struct refinement_state {
    mpz_class x;
    mpz_class difference;
    unsigned long two_exponent;
};

/**
 * @brief Count the powers of two dividing a nonzero integer.
 * @param value Nonzero integer to inspect.
 * @return The exponent of two in the value's prime factorization.
 */
unsigned long two_adic_valuation(const mpz_class& value)
{
    if (value == 0) {
        throw std::invalid_argument(
            "The 2-adic valuation is undefined for zero.");
    }

    mpz_class magnitude = value < 0 ? -value : value;
    unsigned long exponent = 0;
    while (mpz_even_p(magnitude.get_mpz_t()) != 0) {
        mpz_tdiv_q_2exp(magnitude.get_mpz_t(), magnitude.get_mpz_t(), 1);
        ++exponent;
    }
    return exponent;
}

/**
 * @brief Determine whether a nonnegative integer is a perfect square.
 * @param value Integer to inspect.
 * @return True when value equals the square of its integer square root.
 */
bool is_perfect_square(const mpz_class& value)
{
    if (value < 0) {
        return false;
    }

    mpz_class square_root;
    mpz_sqrt(square_root.get_mpz_t(), value.get_mpz_t());
    return square_root * square_root == value;
}

/**
 * @brief Verify the active base is {2} and return its root.
 * @param number Integer N defining the root modulo 2.
 * @param factor_base Active factor base supplied by the caller.
 * @return The root of number modulo 2.
 */
mpz_class validate_singleton_base(
    const mpz_class& number,
    const std::vector<factor_base_entry>& factor_base)
{
    if (factor_base.size() != 1 || factor_base.front().prime != 2) {
        throw std::invalid_argument(
            "Relation generation requires the singleton factor base {2}.");
    }

    const mpz_class root = positive_mod(number, mpz_class(2));
    if (factor_base.front().root1 != root ||
        factor_base.front().root2 != root) {
        throw std::invalid_argument("The factor-base root modulo 2 is invalid.");
    }
    return root;
}

/**
 * @brief Compute the initial odd x and its 2-adic refinement state.
 * @param target Multiplied value g*N under the square root.
 * @return State for the first odd integer at or above sqrt(target), strictly
 *         above it when target is itself an odd square.
 */
refinement_state make_initial_state(const mpz_class& target)
{
    mpz_class x;
    mpz_sqrt(x.get_mpz_t(), target.get_mpz_t());
    if (x * x < target) {
        ++x;
    }
    if (mpz_odd_p(x.get_mpz_t()) == 0) {
        ++x;
    }
    if (x * x == target) {
        x += 2;
    }

    const mpz_class difference = x * x - target;
    const unsigned long two_exponent = two_adic_valuation(difference);
    if (two_exponent <= 2) {
        throw std::invalid_argument(
            "The refinement requires v_2(x^2 - g*N) to be greater than 2.");
    }
    return {x, difference, two_exponent};
}

/**
 * @brief Factor one relation by the selected trial primes.
 * @param state Current x and x^2-g*N value.
 * @param root Root modulo the singleton prime 2.
 * @param multiplier Odd multiplier g in the relation polynomial.
 * @param exponent_primes Trial primes used for this relation.
 * @return Relation data, including the unresolved cofactor.
 */
crt_relation make_relation(
    const refinement_state& state, const mpz_class& root,
    const mpz_class& multiplier,
    const std::vector<mpz_class>& exponent_primes)
{
    mpz_class remaining = state.difference < 0 ? -state.difference
                                                : state.difference;
    std::vector<mpz_class> relation_primes = exponent_primes;
    std::vector<unsigned long> exponents(relation_primes.size(), 0);
    for (std::size_t index = 0; index < relation_primes.size(); ++index) {
        const mpz_class& prime = relation_primes[index];
        while (remaining != 0 && remaining % prime == 0) {
            remaining /= prime;
            ++exponents[index];
        }
    }

    if (remaining > 1 &&
        mpz_probab_prime_p(remaining.get_mpz_t(), 25) != 0) {
        const auto position = std::lower_bound(
            relation_primes.begin(), relation_primes.end(), remaining);
        const std::size_t index = static_cast<std::size_t>(
            position - relation_primes.begin());
        if (position == relation_primes.end() || *position != remaining) {
            relation_primes.insert(position, remaining);
            exponents.insert(exponents.begin() +
                                 static_cast<std::ptrdiff_t>(index),
                             1);
        } else {
            ++exponents[index];
        }
    }

    std::vector<std::pair<mpz_class, mpz_class>> selected_roots{
        {mpz_class(2), root}};
    return {state.x, multiplier, mpz_class(2), std::move(relation_primes),
            std::move(exponents), remaining, std::move(selected_roots)};
}

/**
 * @brief Apply the parity-dependent t rule and update the 2-adic state.
 * @param target Integer g*N under the square.
 * @param state Current x and the exponent l.
 * @return State for x+2t with any new powers of two included in l'.
 */
refinement_state advance_refinement(const mpz_class& target,
                                    const refinement_state& state)
{
    const unsigned long t_exponent =
        state.two_exponent % 2 == 0
            ? state.two_exponent
            : state.two_exponent - 1;
    mpz_class t1 = 1;
    mpz_mul_2exp(t1.get_mpz_t(), t1.get_mpz_t(), t_exponent);
    mpz_class t = t1;
    if (t1 > state.x) {
        const mpz_class t2 = t1 - state.x;
        if (t2 < t) {
            t = t2;
        }
    }

    const mpz_class next_x = state.x + 2 * t;
    const mpz_class next_difference = next_x * next_x - target;
    mpz_class next_cofactor;
    mpz_tdiv_q_2exp(next_cofactor.get_mpz_t(), next_difference.get_mpz_t(),
                    state.two_exponent);
    unsigned long next_two_exponent = state.two_exponent;
    while (mpz_even_p(next_cofactor.get_mpz_t()) != 0) {
        mpz_tdiv_q_2exp(next_cofactor.get_mpz_t(),
                        next_cofactor.get_mpz_t(), 1);
        ++next_two_exponent;
    }
    return {next_x, next_difference, next_two_exponent};
}

} // namespace

mpz_class select_relation_multiplier(const mpz_class& number)
{
    if (number <= 1 || mpz_even_p(number.get_mpz_t()) != 0) {
        throw std::invalid_argument(
            "Multiplier selection requires an odd integer greater than 1.");
    }

    const mpz_class residue = positive_mod(number, mpz_class(8));
    for (unsigned long candidate = 1; candidate < 8; candidate += 2) {
        if (positive_mod(mpz_class(candidate) * residue, mpz_class(8)) == 1) {
            return mpz_class(candidate);
        }
    }
    throw std::logic_error("No odd relation multiplier was found modulo 8.");
}

/**
 * @brief Generate refined relations for x^2 - multiplier*number.
 * @param number Integer N in the relation polynomial x^2 - g*N.
 * @param multiplier Odd multiplier g in the relation polynomial.
 * @param factor_base Singleton factor base containing prime 2 and its root.
 * @param trial_division_primes Primes used to factor each absolute difference.
 * @param max_relations Maximum number of x values to examine.
 * @return Relations beginning at adjusted ceil(sqrt(g*N)); generation stops
 *         when the remaining cofactor is a perfect square.
 */
std::vector<crt_relation> generate_crt_relations(
    const mpz_class& number,
    const mpz_class& multiplier,
    const std::vector<factor_base_entry>& factor_base,
    const std::set<mpz_class>& trial_division_primes,
    std::size_t max_relations)
{
    if (number <= 1 || max_relations == 0) {
        return {};
    }
    if (multiplier <= 0 || mpz_even_p(multiplier.get_mpz_t()) != 0) {
        throw std::invalid_argument("The relation multiplier must be odd.");
    }

    const mpz_class target = multiplier * number;
    const mpz_class root = validate_singleton_base(target, factor_base);
    std::set<mpz_class> relation_primes = trial_division_primes;
    relation_primes.insert(2);
    const std::vector<mpz_class> exponent_primes(
        relation_primes.begin(), relation_primes.end());

    refinement_state state = make_initial_state(target);
    std::vector<crt_relation> relations;
    relations.reserve(max_relations);
    bool remaining_cofactor_is_square = false;
    for (std::size_t index = 0;
         index < max_relations && !remaining_cofactor_is_square; ++index) {
        crt_relation relation =
            make_relation(state, root, multiplier, exponent_primes);
        remaining_cofactor_is_square =
            is_perfect_square(relation.remaining_cofactor);
        relations.push_back(std::move(relation));
        if (!remaining_cofactor_is_square && index + 1 < max_relations) {
            state = advance_refinement(target, state);
        }
    }
    return relations;
}

std::vector<crt_relation> generate_crt_relations(
    const mpz_class& number,
    const std::vector<factor_base_entry>& factor_base,
    const std::set<mpz_class>& trial_division_primes,
    std::size_t max_relations)
{
    return generate_crt_relations(number, mpz_class(1), factor_base,
                                  trial_division_primes, max_relations);
}

} // namespace factoring
