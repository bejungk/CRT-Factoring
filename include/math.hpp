#pragma once

#include <gmpxx.h>

namespace factoring {

/**
 * @brief Reduce an integer to its canonical residue modulo a positive modulus.
 * @param value Integer to reduce.
 * @param modulus Positive modulus.
 * @return The representative in [0, modulus).
 */
mpz_class positive_mod(const mpz_class& value, const mpz_class& modulus);

/**
 * @brief Check whether root is a square root of number modulo prime.
 * @param root Candidate residue.
 * @param number Integer whose residue is being square-rooted.
 * @param prime Modulus.
 * @return True if root^2 == number (mod prime).
 */
bool is_valid_root(const mpz_class& root, const mpz_class& number,
                   const mpz_class& prime);

/**
 * @brief Find a square root of a quadratic residue modulo a prime.
 * @param value Residue whose square root is required.
 * @param prime Prime modulus.
 * @param root Output parameter receiving a root when one is found.
 * @return True if a root was produced, otherwise false.
 */
bool square_root_mod_prime(const mpz_class& value, const mpz_class& prime,
                           mpz_class& root);

} // namespace factoring
