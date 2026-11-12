#include "crt.hpp"

#include <iostream>
#include <stdexcept>
#include <vector>

namespace {

/**
 * @brief Check the refinement sequence and perfect-square stop condition.
 * @return True when all expected values match.
 */
bool check_refinement_vector()
{
    const mpz_class number = 5 * 109;
    const std::vector<mpz_class> expected_x{25, 57};
    const std::vector<mpz_class> expected_remaining_cofactors{5, 169};
    const std::vector<factoring::factor_base_entry> factor_base =
        factoring::make_factor_base(number);
    const std::vector<factoring::crt_relation> relations =
        factoring::generate_crt_relations(number, factor_base, {}, 20);

    bool matches = relations.size() == expected_x.size();
    for (std::size_t index = 0;
         index < relations.size() && index < expected_x.size(); ++index) {
        matches = matches && relations[index].x == expected_x[index] &&
                  relations[index].remaining_cofactor ==
                      expected_remaining_cofactors[index];
    }
    return matches;
}

/**
 * @brief Check the odd-l shift rule on N=11*43.
 * @return True when both expected relations match.
 */
bool check_odd_exponent_vector()
{
    const mpz_class number = 11 * 43;
    const std::vector<mpz_class> expected_x{23, 31};
    const std::vector<mpz_class> expected_remaining_cofactors{7, 61};
    const std::vector<factoring::factor_base_entry> factor_base =
        factoring::make_factor_base(number);
    const std::vector<factoring::crt_relation> relations =
        factoring::generate_crt_relations(number, factor_base, {}, 2);

    bool matches = relations.size() == expected_x.size();
    for (std::size_t index = 0;
         index < relations.size() && index < expected_x.size(); ++index) {
        matches = matches && relations[index].x == expected_x[index] &&
                  relations[index].remaining_cofactor ==
                      expected_remaining_cofactors[index];
    }
    return matches;
}

/**
 * @brief Check that the 13*91 example is rejected because its initial l is 1.
 * @return True when relation generation reports the unsupported case.
 */
bool check_l_equals_one_vector()
{
    const mpz_class number = 13 * 91;
    const mpz_class x = 35;
    const mpz_class difference = x * x - number;
    const mpz_class odd_part = difference / 2;
    const std::vector<factoring::factor_base_entry> factor_base =
        factoring::make_factor_base(number);

    bool rejected = false;
    try {
        factoring::generate_crt_relations(number, factor_base, {}, 1);
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    return x * x == 1225 && difference == 42 && odd_part == 21 && rejected;
}

/**
 * @brief Check that the 13*89 example is rejected because its initial l is 2.
 * @return True when relation generation reports the non-refinable case.
 */
bool check_l_equals_two_vector()
{
    const mpz_class number = 13 * 89;
    const mpz_class x = 35;
    const mpz_class difference = x * x - number;
    const mpz_class odd_part = difference / 4;
    const std::vector<factoring::factor_base_entry> factor_base =
        factoring::make_factor_base(number);

    bool rejected = false;
    try {
        factoring::generate_crt_relations(number, factor_base, {}, 1);
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    return x * x == 1225 && difference == 68 && odd_part == 17 && rejected;
}

/**
 * @brief Check immediate square-cofactor termination for N=73*89.
 * @return True when x=81 yields the single expected relation.
 */
bool check_initial_square_cofactor_vector()
{
    const mpz_class number = 73 * 89;
    const std::vector<factoring::factor_base_entry> factor_base =
        factoring::make_factor_base(number);
    const std::vector<factoring::crt_relation> relations =
        factoring::generate_crt_relations(number, factor_base, {}, 20);

    return relations.size() == 1 && relations.front().x == 81 &&
           relations.front().x * relations.front().x - number == 64 &&
           relations.front().exponent_vector.size() == 1 &&
           relations.front().exponent_vector.front() == 6 &&
           relations.front().remaining_cofactor == 1;
}

} // namespace

/**
 * @brief Run the small-number relation-refinement test vectors.
 * @return Zero when all vectors pass, otherwise one.
 */
int main()
{
    const bool refinement_matches = check_refinement_vector();
    const bool odd_exponent_matches = check_odd_exponent_vector();
    const bool l_one_is_rejected = check_l_equals_one_vector();
    const bool l_two_is_rejected = check_l_equals_two_vector();
    const bool initial_square_stops = check_initial_square_cofactor_vector();
    if (!refinement_matches || !odd_exponent_matches || !l_one_is_rejected ||
        !l_two_is_rejected || !initial_square_stops) {
        std::cerr << "CRT refinement test vector mismatch.\n";
        return 1;
    }

    std::cout << "CRT refinement test vectors passed.\n";
    return 0;
}
