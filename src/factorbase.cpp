#include "factor_base.hpp"
#include "math.hpp"

#include <stdexcept>

namespace factoring {
namespace {

/**
 * @brief Check whether N has a square root modulo a candidate prime.
 * @param number Integer being factored.
 * @param prime Candidate prime.
 * @return True when N is a quadratic residue modulo prime.
 */
bool is_eligible_prime(const mpz_class& number, const mpz_class& prime)
{
    bool eligible = prime == 2;
    if (!eligible) {
        mpz_class root;
        eligible = square_root_mod_prime(positive_mod(number, prime), prime,
                                         root);
    }
    return eligible;
}

} // namespace

/**
 * @brief Construct the minimal factor base containing only the prime 2.
 * @param number Integer N whose residue determines the root modulo 2.
 * @return A singleton factor base containing the root of N modulo 2.
 */
std::vector<factor_base_entry> make_factor_base(const mpz_class& number)
{
    if (number <= 1) {
        throw std::invalid_argument("N must be greater than 1.");
    }

    const mpz_class root = positive_mod(number, mpz_class(2));
    return {{mpz_class(2), root, root}};
}

/**
 * @brief Collect eligible primes for trial division independently of the base.
 * @param number Integer whose square residues define eligibility.
 * @param count Number of eligible primes to collect.
 * @return The smallest eligible primes in ascending order.
 */
std::set<mpz_class> make_trial_division_primes(const mpz_class& number,
                                               std::size_t count)
{
    if (number <= 1) {
        throw std::invalid_argument("N must be greater than 1.");
    }

    std::set<mpz_class> primes;
    mpz_class candidate = 1;
    while (primes.size() < count) {
        mpz_nextprime(candidate.get_mpz_t(), candidate.get_mpz_t());
        if (is_eligible_prime(number, candidate)) {
            primes.insert(candidate);
        }
    }
    return primes;
}

} // namespace factoring
