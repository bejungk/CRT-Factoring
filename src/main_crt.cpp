#include "crt.hpp"
#include "math.hpp"

#include <algorithm>
#include <exception>
#include <iostream>
#include <random>
#include <set>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

namespace {

using relation_factor_signature = std::tuple<
    bool, mpz_class, std::vector<std::pair<mpz_class, unsigned long>>,
    mpz_class>;

/**
 * @brief Build a canonical factor signature for duplicate-relation checking.
 * @param exponent_primes Primes represented in the relation's exponent vector.
 * @param exponent_vector Multiplicity for each corresponding prime.
 * @param cofactor Unresolved remainder after trial division.
 * @param multiplier g used in x^2 - g*N.
 * @param signed_difference Signed value of x^2-g*N, used to distinguish signs.
 * @return Sign, nonzero prime/exponent pairs, and unresolved composite part.
 */
relation_factor_signature make_relation_signature(
    const std::vector<mpz_class>& exponent_primes,
    const std::vector<unsigned long>& exponent_vector,
    const mpz_class& cofactor,
    const mpz_class& multiplier,
    const mpz_class& signed_difference)
{
    std::vector<std::pair<mpz_class, unsigned long>> factors;
    for (std::size_t i = 0; i < exponent_primes.size(); ++i) {
        if (exponent_vector[i] != 0) {
            factors.emplace_back(exponent_primes[i], exponent_vector[i]);
        }
    }
    std::sort(factors.begin(), factors.end(),
              [](const auto& first, const auto& second) {
                  return first.first < second.first;
              });

    // Probable-prime cofactors already appear as a factor with exponent one.
    const bool cofactor_is_recorded_as_factor =
        cofactor > 1 && mpz_probab_prime_p(cofactor.get_mpz_t(), 25) != 0;
    const mpz_class unresolved_cofactor =
        cofactor_is_recorded_as_factor ? mpz_class(1) : cofactor;
    return {signed_difference < 0, multiplier, std::move(factors),
            unresolved_cofactor};
}

class gmp_random_state {
public:
    /**
     * @brief Initialize and own a GMP Mersenne-Twister random state.
     * @param seed Seed passed to GMP's random-state initializer.
     */
    explicit gmp_random_state(unsigned long seed)
    {
        gmp_randinit_mt(state_);
        mpz_class seed_value(seed);
        gmp_randseed(state_, seed_value.get_mpz_t());
    }

    /**
     * @brief Release the owned GMP random state.
     */
    ~gmp_random_state() { gmp_randclear(state_); }

    /**
     * @brief Return the owned state for GMP random-number functions.
     * @return Mutable reference to the state.
     */
    gmp_randstate_t& get() { return state_; }

    /**
     * @brief Disable copying because the GMP state has unique ownership.
     */
    gmp_random_state(const gmp_random_state&) = delete;
    /**
     * @brief Disable copy assignment for the uniquely owned GMP state.
     */
    gmp_random_state& operator=(const gmp_random_state&) = delete;

private:
    gmp_randstate_t state_;
};

/**
 * @brief Generate a prime with exactly 256 bits.
 * @param generator GMP random state used to sample candidates.
 * @return A prime in the requested range.
 */
mpz_class random_256_bit_prime(gmp_randstate_t& generator)
{
    // Set the top bit so every candidate and its next prime are 256-bit,
    // retrying only if the next prime crosses the 256-bit upper bound.
    const mpz_class lower_bound = mpz_class(1) << 255;
    const mpz_class upper_bound = mpz_class(1) << 256;
    mpz_class prime;
    bool found_prime = false;
    while (!found_prime) {
        mpz_class candidate;
        mpz_urandomb(candidate.get_mpz_t(), generator, 255);
        candidate += lower_bound;

        mpz_nextprime(prime.get_mpz_t(), candidate.get_mpz_t());
        if (prime < upper_bound) {
            found_prime = true;
        }
    }
    return prime;
}

} // namespace

/**
 * @brief Run the random semiprime, factor-base, and relation-generation demo.
 * @return Zero on success, or one when relation generation fails.
 */
int main()
{
    constexpr std::size_t target_relation_count = 20'000;
    constexpr std::size_t trial_division_prime_count = 200;
    std::random_device entropy;
    gmp_random_state generator(static_cast<unsigned long>(entropy()));
    try {
        const mpz_class first_prime = random_256_bit_prime(generator.get());
        mpz_class second_prime = random_256_bit_prime(generator.get());
        while (second_prime == first_prime) {
            second_prime = random_256_bit_prime(generator.get());
        }
        const mpz_class number = first_prime * second_prime;
        const mpz_class multiplier =
            factoring::select_relation_multiplier(number);
        const mpz_class relation_target = multiplier * number;
        const std::vector<factoring::factor_base_entry> factor_base =
            factoring::make_factor_base(relation_target);
        std::set<mpz_class> trial_division_primes =
            factoring::make_trial_division_primes(
                relation_target, trial_division_prime_count);
        std::cout << "Using " << trial_division_primes.size()
                  << " eligible primes for trial division with multiplier g="
                  << multiplier << ".\n";
        std::set<mpz_class> seen_cofactors;
        std::set<relation_factor_signature> seen_relation_signatures;
        std::size_t duplicate_cofactor_count = 0;
        std::size_t duplicate_relation_count = 0;
        const std::vector<factoring::crt_relation> relations =
            factoring::generate_crt_relations(
                number, multiplier, factor_base, trial_division_primes,
                target_relation_count);
        for (std::size_t i = 0; i < relations.size(); ++i) {
            const factoring::crt_relation& relation = relations[i];
            const bool duplicate_cofactor =
                !seen_cofactors.insert(relation.remaining_cofactor).second;
            const mpz_class signed_difference =
                relation.x * relation.x - relation.multiplier * number;
            const bool duplicate_relation =
                !seen_relation_signatures
                     .insert(make_relation_signature(
                         relation.exponent_primes, relation.exponent_vector,
                         relation.remaining_cofactor, relation.multiplier,
                         signed_difference))
                     .second;
            duplicate_cofactor_count += duplicate_cofactor ? 1 : 0;
            duplicate_relation_count += duplicate_relation ? 1 : 0;
        }

        std::cout << "Generated " << relations.size() << '/'
                  << target_relation_count << " relation(s) from the singleton "
                  << "factor base {2}.\n"
                  << "Duplicate cofactors: " << duplicate_cofactor_count << '\n'
                  << "Duplicate relations: " << duplicate_relation_count << '\n';
    } catch (const std::exception& error) {
        std::cerr << "Test failed: " << error.what() << '\n';
        return 1;
    }

    return 0;
}
