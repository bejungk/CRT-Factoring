#pragma once

#include "factor_base.hpp"

#include <cstddef>
#include <set>
#include <utility>
#include <vector>

namespace factoring {

struct crt_relation {
    /**
     * Value x in the relation x^2 - multiplier*number.
     */
    mpz_class x;
    /** Odd multiplier g in the relation x^2 - g*N. */
    mpz_class multiplier;
    mpz_class modulus;
    /** Prime order corresponding to exponent_vector. */
    std::vector<mpz_class> exponent_primes;
    /**
     * Exponents for each listed prime. A probable-prime remaining cofactor is
     * included with exponent one, while remaining_cofactor is retained for
     * reporting.
     */
    std::vector<unsigned long> exponent_vector;
    /**
     * Portion left after trial division by the selected eligible primes.
     * If probable prime, it is also represented in exponent_primes/vector.
     */
    mpz_class remaining_cofactor;
    /**
     * (prime, selected root) pair used by this relation's factor base.
     */
    std::vector<std::pair<mpz_class, mpz_class>> prime_roots;
};

/**
 * @brief Choose a small odd multiplier making g*number congruent to 1 mod 8.
 * @param number Odd integer N used in the relation polynomial.
 * @return The unique multiplier in {1,3,5,7} satisfying g*N == 1 (mod 8).
 */
mpz_class select_relation_multiplier(const mpz_class& number);

/**
 * @brief Generate refined relations for x^2 - multiplier*number.
 * @param number Integer N in the relation polynomial.
 * @param multiplier Odd multiplier g; the polynomial is x^2 - g*N.
 * @param factor_base Singleton factor base containing prime 2 for g*N.
 * @param trial_division_primes Primes to trial-divide from each relation's
 *        absolute difference.
 * @param max_relations Maximum number of x values to examine.
 * @return Relations beginning at adjusted ceil(sqrt(g*number)); generation ends
 *         early when trial division leaves a perfect-square cofactor.
 */
std::vector<crt_relation> generate_crt_relations(
    const mpz_class& number,
    const mpz_class& multiplier,
    const std::vector<factor_base_entry>& factor_base,
    const std::set<mpz_class>& trial_division_primes,
    std::size_t max_relations);

/**
 * @brief Generate refined relations for the unscaled polynomial x^2 - number.
 * @details This compatibility overload uses multiplier 1. Use
 *          select_relation_multiplier() and the overload above to avoid the
 *          initial low 2-adic valuations for odd inputs.
 */
std::vector<crt_relation> generate_crt_relations(
    const mpz_class& number,
    const std::vector<factor_base_entry>& factor_base,
    const std::set<mpz_class>& trial_division_primes,
    std::size_t max_relations);

} // namespace factoring
