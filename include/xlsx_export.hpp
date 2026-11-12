#pragma once

#include "crt.hpp"

#include <cstdint>
#include <string>
#include <vector>

namespace factoring {

/** @brief One cofactor observed in a semiprime experiment. */
struct cofactor_experiment_result {
    std::size_t experiment;
    std::size_t combination;
    mpz_class first_prime;
    mpz_class second_prime;
    mpz_class number;
    mpz_class modulus;
    mpz_class x;
    mpz_class cofactor;
    bool duplicate_cofactor;
    bool duplicate_relation;
    std::vector<mpz_class> exponent_primes;
    std::vector<unsigned long> exponent_vector;
};

/**
 * @brief Write the summary, factor base, and CRT relation data to an XLSX workbook.
 * @param path Destination file path.
 * @param number Integer N represented by the workbook.
 * @param bound Reported factor-base bound value.
 * @param factor_base Selected primes and modular roots.
 * @param relations CRT relations to include.
 */
void export_results_xlsx(
    const std::string& path, const mpz_class& number, std::uint32_t bound,
    const std::vector<factor_base_entry>& factor_base,
    const std::vector<crt_relation>& relations);

/**
 * @brief Write collected semiprime experiment cofactors to a one-sheet workbook.
 * @param path Destination XLSX file.
 * @param results Cofactor observations from all experiments.
 */
void export_cofactors_xlsx(
    const std::string& path,
    const std::vector<cofactor_experiment_result>& results);

} // namespace factoring
