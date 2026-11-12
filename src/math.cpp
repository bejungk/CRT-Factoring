#include "math.hpp"

namespace factoring {

/**
 * @brief Reduce an integer to its canonical residue modulo a positive modulus.
 * @param value Integer to reduce.
 * @param modulus Positive modulus.
 * @return The representative in [0, modulus).
 */
mpz_class positive_mod(const mpz_class& value, const mpz_class& modulus)
{
    // GMP remainder can be negative; shift it into the canonical residue range.
    mpz_class result = value % modulus;
    if (result < 0) {
        result += modulus;
    }
    return result;
}

/**
 * @brief Check whether root is a square root of number modulo prime.
 * @param root Candidate residue.
 * @param number Integer whose residue is being square-rooted.
 * @param prime Modulus.
 * @return True if root^2 == number (mod prime).
 */
bool is_valid_root(const mpz_class& root, const mpz_class& number,
                   const mpz_class& prime)
{
    // Compare residues so the check does not depend on the chosen representative.
    return positive_mod(root * root, prime) == positive_mod(number, prime);
}

/**
 * @brief Find a square root of a quadratic residue modulo a prime.
 * @param value Residue whose square root is required.
 * @param prime Prime modulus.
 * @param root Output parameter receiving a root when one is found.
 * @return True if a root was produced, otherwise false.
 */
bool square_root_mod_prime(const mpz_class& value, const mpz_class& prime,
                           mpz_class& root)
{
    // Tonelli-Shanks for odd primes, with shortcuts for p=2 and p=3 (mod 4).
    // Reference: D. Shanks, "Five Number-Theoretic Algorithms" (1973).
    if (prime == 2) {
        root = positive_mod(value, prime);
        return true;
    }
    if (mpz_legendre(value.get_mpz_t(), prime.get_mpz_t()) != 1) {
        return false;
    }

    // Factor p-1 as q*2^s, q odd; reduce t's order until it reaches 1.
    mpz_class q = prime - 1;
    unsigned long s = 0;
    while (mpz_even_p(q.get_mpz_t())) {
        q >>= 1;
        ++s;
    }
    if (s == 1) {
        const mpz_class exponent = (prime + 1) / 4;
        mpz_powm(root.get_mpz_t(), value.get_mpz_t(), exponent.get_mpz_t(),
                 prime.get_mpz_t());
        return true;
    }

    mpz_class z = 2;
    while (mpz_legendre(z.get_mpz_t(), prime.get_mpz_t()) != -1) {
        ++z;
    }
    mpz_class c;
    mpz_powm(c.get_mpz_t(), z.get_mpz_t(), q.get_mpz_t(), prime.get_mpz_t());
    const mpz_class root_exponent = (q + 1) / 2;
    mpz_powm(root.get_mpz_t(), value.get_mpz_t(), root_exponent.get_mpz_t(),
             prime.get_mpz_t());
    mpz_class t;
    mpz_powm(t.get_mpz_t(), value.get_mpz_t(), q.get_mpz_t(),
             prime.get_mpz_t());
    unsigned long m = s;
    while (t != 1) {
        mpz_class t_squared = (t * t) % prime;
        unsigned long i = 1;
        while (i < m && t_squared != 1) {
            t_squared = (t_squared * t_squared) % prime;
            ++i;
        }
        if (i == m) {
            return false;
        }
        mpz_class exponent = 1;
        exponent <<= (m - i - 1);
        mpz_class b;
        mpz_powm(b.get_mpz_t(), c.get_mpz_t(), exponent.get_mpz_t(),
                 prime.get_mpz_t());
        const mpz_class b_squared = (b * b) % prime;
        root = (root * b) % prime;
        t = (t * b_squared) % prime;
        c = b_squared;
        m = i;
    }
    return true;
}

} // namespace factoring
