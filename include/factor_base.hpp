#pragma once

#include <gmpxx.h>

#include <set>
#include <vector>

namespace factoring {

struct factor_base_entry {
    /**
     * Prime p and the two roots r satisfying r^2 == N (mod p).
     */
    mpz_class prime;
    mpz_class root1;
    mpz_class root2;
};

/**
 * @brief Build the singleton factor base containing only prime 2.
 * @param number Integer being factored; must be greater than 1.
 * @return The prime 2 and its root modulo 2.
 */
std::vector<factor_base_entry> make_factor_base(const mpz_class& number);

/**
 * @brief Return the first requested number of eligible primes for trial division.
 * @param number Integer whose quadratic residues define eligibility.
 * @param count Number of eligible primes to collect.
 * @return The smallest eligible primes in ascending order.
 */
std::set<mpz_class> make_trial_division_primes(const mpz_class& number,
                                               std::size_t count);

} // namespace factoring
