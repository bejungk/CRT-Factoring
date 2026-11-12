#include "math.hpp"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace {
using factoring::positive_mod;

// This build links the GMP C DLL, not an ABI-compatible gmpxx stream library.
std::ostream& operator<<(std::ostream& out, const mpz_class& value)
{
    return out << value.get_str();
}

struct options {
    mpz_class number = 1022117; // 1009 * 1013; factors are never search inputs.
    std::size_t base_size = 5;
    unsigned long levels = 3;
    unsigned long radius = 32;
    std::size_t max_candidates = 250000;
    mpz_class max_multiplier = 1000000000;
    std::size_t target = 0; // Defaults to |P|.
    std::string csv;
};

struct root_set {
    mpz_class modulus;
    std::vector<mpz_class> roots;
};

struct relation {
    mpz_class x, g, modulus, b, u, k, cofactor;
    std::vector<unsigned long> forced, exponents;
};

struct search_result {
    mpz_class g0, cofactor = 0, direct_factor = 0;
    std::vector<mpz_class> primes;
    std::vector<relation> relations;
    std::size_t candidates = 0, duplicates = 0, profiles = 0;
    std::size_t unit_cofactors = 0, discovered_at = 0;
    std::string stopped;
};

mpz_class power(const mpz_class& p, unsigned long e)
{
    mpz_class result;
    mpz_pow_ui(result.get_mpz_t(), p.get_mpz_t(), e);
    return result;
}

mpz_class floor_div(const mpz_class& a, const mpz_class& b)
{
    mpz_class result;
    mpz_fdiv_q(result.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t());
    return result;
}

mpz_class gcd(const mpz_class& a, const mpz_class& b)
{
    mpz_class result;
    mpz_gcd(result.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t());
    return result;
}

mpz_class inverse(const mpz_class& a, const mpz_class& modulus)
{
    mpz_class result;
    if (!mpz_invert(result.get_mpz_t(), a.get_mpz_t(), modulus.get_mpz_t())) {
        throw std::logic_error("Noninvertible CRT/Hensel coefficient.");
    }
    return result;
}

// Odd primes use nonsingular Hensel lifting. At 2, explicitly retain every
// surviving child, including all four roots when the exponent is at least 3.
root_set lift_roots(const mpz_class& d, const mpz_class& p, unsigned long e)
{
    if (e == 0 || gcd(d, p) != 1) {
        throw std::invalid_argument("Lifting requires e > 0 and gcd(D,p)=1.");
    }
    root_set result{p, {}};
    mpz_class r;
    if (!factoring::square_root_mod_prime(positive_mod(d, p), p, r)) {
        return {power(p, e), {}};
    }
    result.roots.push_back(r);
    if (p != 2) {
        result.roots.push_back(p - r);
    }
    for (unsigned long j = 1; j < e; ++j) {
        const mpz_class next_modulus = result.modulus * p;
        std::vector<mpz_class> next;
        for (const auto& old : result.roots) {
            if (p == 2) {
                for (unsigned long digit = 0; digit < 2; ++digit) {
                    const mpz_class child = old + digit * result.modulus;
                    if (positive_mod(child * child - d, next_modulus) == 0) {
                        next.push_back(child);
                    }
                }
            } else {
                const mpz_class quotient = (old * old - d) / result.modulus;
                const mpz_class digit = positive_mod(-quotient * inverse(2 * old, p), p);
                next.push_back(old + result.modulus * digit);
            }
        }
        result = {next_modulus, std::move(next)};
    }
    std::sort(result.roots.begin(), result.roots.end());
    return result;
}

root_set combine(const root_set& a, const root_set& b)
{
    root_set result{a.modulus * b.modulus, {}};
    const mpz_class inv = inverse(a.modulus, b.modulus);
    for (const auto& x : a.roots) {
        for (const auto& y : b.roots) {
            const mpz_class t = positive_mod((y - x) * inv, b.modulus);
            result.roots.push_back(x + a.modulus * t);
        }
    }
    std::sort(result.roots.begin(), result.roots.end());
    return result;
}

void extract_cofactor(relation& row, const mpz_class& n,
                      const std::vector<mpz_class>& primes)
{
    row.cofactor = row.x * row.x - row.g * n;
    if (row.cofactor <= 0) {
        throw std::logic_error("Only positive relations are collected.");
    }
    row.exponents.assign(primes.size(), 0);
    for (std::size_t i = 0; i < primes.size(); ++i) {
        while (mpz_divisible_p(row.cofactor.get_mpz_t(), primes[i].get_mpz_t())) {
            row.cofactor /= primes[i];
            ++row.exponents[i];
        }
    }
}

bool valid_relation(const relation& row, const mpz_class& n,
                    const std::vector<mpz_class>& primes)
{
    if (row.cofactor <= 0 || row.exponents.size() != primes.size()) return false;
    mpz_class reconstructed = row.cofactor;
    for (std::size_t i = 0; i < primes.size(); ++i) {
        if (row.cofactor % primes[i] == 0) return false;
        reconstructed *= power(primes[i], row.exponents[i]);
    }
    return row.x * row.x - row.g * n == reconstructed;
}

search_result search(const options& opt)
{
    search_result out;
    // Every odd N is its own inverse modulo 8.
    out.g0 = positive_mod(opt.number, mpz_class(8));
    if (opt.max_multiplier < out.g0) {
        throw std::invalid_argument("max-multiplier is smaller than the initial g.");
    }
    const mpz_class d = out.g0 * opt.number;
    mpz_class prime = 1;
    while (out.primes.size() < opt.base_size) {
        mpz_nextprime(prime.get_mpz_t(), prime.get_mpz_t());
        const mpz_class divisor = gcd(prime, opt.number);
        if (divisor > 1) {
            if (divisor < opt.number) out.direct_factor = divisor;
            continue;
        }
        if (gcd(prime, d) != 1) continue;
        mpz_class root;
        if (factoring::square_root_mod_prime(positive_mod(d, prime), prime, root)) {
            out.primes.push_back(prime);
        }
    }

    // Cache lifts for each prime. The profile varies each exponent independently.
    std::vector<std::vector<root_set>> lifts;
    for (const auto& p : out.primes) {
        std::vector<root_set> levels;
        for (unsigned long j = 0; j < opt.levels; ++j) {
            levels.push_back(lift_roots(d, p, (p == 2 ? 3UL : 1UL) + j));
        }
        lifts.push_back(std::move(levels));
    }
    using signature = std::tuple<mpz_class, mpz_class, std::vector<unsigned long>>;
    std::set<signature> seen;
    std::map<mpz_class, std::vector<relation>> buckets;
    std::vector<unsigned long> profile(out.primes.size(), 0);
    mpz_class sqrt_d;
    mpz_sqrt(sqrt_d.get_mpz_t(), d.get_mpz_t());

    bool profiles_left = true;
    while (profiles_left && out.candidates < opt.max_candidates &&
           out.relations.size() < opt.target) {
        ++out.profiles;
        root_set combined{1, {mpz_class(0)}};
        std::vector<unsigned long> forced;
        for (std::size_t i = 0; i < profile.size(); ++i) {
            combined = combine(combined, lifts[i][profile[i]]);
            forced.push_back(profile[i] + (out.primes[i] == 2 ? 3UL : 1UL));
        }
        const mpz_class& m = combined.modulus;
        for (const auto& b : combined.roots) {
            const mpz_class center = floor_div(sqrt_d - b, m);
            for (unsigned long step = 0; step <= 2 * opt.radius; ++step) {
                // Visit the center, then +1,-1,+2,-2,... (both sides of sqrt).
                const long offset = step == 0 ? 0 :
                    (step % 2 ? static_cast<long>((step + 1) / 2)
                              : -static_cast<long>(step / 2));
                const mpz_class u = center + offset;
                const mpz_class x = b + m * u;
                if (x <= 0) continue;
                const mpz_class difference = x * x - d;
                if (difference <= 0) continue;
                const mpz_class f = difference / m;
                // Largest k yielding c_raw > 0, capped by the g bound.
                const mpz_class max_k = floor_div(opt.max_multiplier - out.g0, m);
                const mpz_class best_k = std::min(floor_div(f - 1, opt.number), max_k);
                std::vector<mpz_class> choices{0};
                if (best_k > 0) choices.push_back(best_k);
                for (const auto& k : choices) {
                    if (out.candidates >= opt.max_candidates ||
                        out.relations.size() >= opt.target) break;
                    ++out.candidates;
                    relation row{x, out.g0 + m * k, m, b, u, k, 0, forced, {}};
                    extract_cofactor(row, opt.number, out.primes);
                    if (!valid_relation(row, opt.number, out.primes)) {
                        throw std::logic_error("Relation reconstruction failed.");
                    }
                    if (row.cofactor == 1) {
                        ++out.unit_cofactors;
                        continue; // Study a genuine common external cofactor.
                    }
                    const mpz_class x_mod = positive_mod(row.x, opt.number);
                    const mpz_class opposite = opt.number - x_mod;
                    const mpz_class canonical_root = std::min(x_mod, opposite);
                    if (!seen.emplace(row.cofactor, canonical_root, row.exponents).second) {
                        ++out.duplicates;
                        continue;
                    }
                    if (out.cofactor == 0) {
                        auto& bucket = buckets[row.cofactor];
                        bucket.push_back(row);
                        if (bucket.size() == 2) {
                            out.cofactor = row.cofactor;
                            out.relations = std::move(bucket);
                            out.discovered_at = out.candidates;
                            buckets.clear();
                        }
                    } else if (row.cofactor == out.cofactor) {
                        out.relations.push_back(std::move(row));
                    }
                }
                if (out.candidates >= opt.max_candidates || out.relations.size() >= opt.target) break;
            }
            if (out.candidates >= opt.max_candidates || out.relations.size() >= opt.target) break;
        }
        profiles_left = false;
        for (std::size_t i = 0; i < profile.size(); ++i) {
            if (++profile[i] < opt.levels) { profiles_left = true; break; }
            profile[i] = 0;
        }
    }
    out.stopped = out.relations.size() >= opt.target ? "target reached" :
        out.candidates >= opt.max_candidates ? "candidate budget exhausted" : "lift profiles exhausted";
    return out;
}

struct dependency_report {
    std::size_t rank = 0, nullity = 0, tested = 0;
    mpz_class factor = 0, subset = 0, x = 0, y = 0;
};

dependency_report analyze(const mpz_class& n, const std::vector<mpz_class>& primes,
                          const std::vector<relation>& rows)
{
    dependency_report report;
    if (rows.empty()) return report;
    const mpz_class c = rows.front().cofactor;
    mpz_class sqrt_c;
    mpz_sqrt(sqrt_c.get_mpz_t(), c.get_mpz_t());
    const bool c_square = sqrt_c * sqrt_c == c;
    std::map<unsigned, std::pair<std::uint64_t, mpz_class>> pivots;
    std::vector<mpz_class> kernel;
    for (std::size_t j = 0; j < rows.size(); ++j) {
        if (rows[j].cofactor != c || !valid_relation(rows[j], n, primes)) {
            throw std::logic_error("Invalid common-cofactor group.");
        }
        std::uint64_t bits = c_square ? 0 : (std::uint64_t(1) << primes.size());
        for (std::size_t i = 0; i < primes.size(); ++i) {
            bits |= std::uint64_t(rows[j].exponents[i] % 2) << i;
        }
        mpz_class subset = 0;
        mpz_setbit(subset.get_mpz_t(), static_cast<mp_bitcnt_t>(j));
        bool inserted = false;
        for (unsigned i = 0; i <= primes.size(); ++i) {
            if (!(bits & (std::uint64_t(1) << i))) continue;
            const auto found = pivots.find(i);
            if (found == pivots.end()) {
                pivots.emplace(i, std::make_pair(bits, subset));
                inserted = true;
                break;
            }
            bits ^= found->second.first;
            subset ^= found->second.second;
        }
        if (!inserted) kernel.push_back(subset);
    }
    report.rank = pivots.size();
    report.nullity = kernel.size();
    auto check = [&](const mpz_class& subset) {
        ++report.tested;
        mpz_class x = 1, y = 1;
        unsigned long count = 0;
        std::vector<unsigned long> sums(primes.size(), 0);
        for (std::size_t j = 0; j < rows.size(); ++j) {
            if (!mpz_tstbit(subset.get_mpz_t(), static_cast<mp_bitcnt_t>(j))) continue;
            ++count;
            x = positive_mod(x * rows[j].x, n);
            for (std::size_t i = 0; i < primes.size(); ++i) sums[i] += rows[j].exponents[i];
        }
        for (std::size_t i = 0; i < primes.size(); ++i) {
            if (sums[i] % 2) throw std::logic_error("Odd prime exponent in dependency.");
            mpz_class term;
            mpz_powm_ui(term.get_mpz_t(), primes[i].get_mpz_t(), sums[i] / 2, n.get_mpz_t());
            y = positive_mod(y * term, n);
        }
        mpz_class base = c, term;
        if (c_square) mpz_sqrt(base.get_mpz_t(), c.get_mpz_t());
        else if (count % 2) throw std::logic_error("Odd cofactor exponent in dependency.");
        mpz_powm_ui(term.get_mpz_t(), base.get_mpz_t(), c_square ? count : count / 2, n.get_mpz_t());
        y = positive_mod(y * term, n);
        if (positive_mod(x * x - y * y, n) != 0) {
            throw std::logic_error("Dependency failed square-congruence check.");
        }
        for (const mpz_class difference : {mpz_class(x - y), mpz_class(x + y)}) {
            const mpz_class divisor = gcd(difference, n);
            if (divisor > 1 && divisor < n && report.factor == 0) {
                report.factor = divisor; report.subset = subset; report.x = x; report.y = y;
            }
        }
    };
    // Check every basis dependency, then combinations of the first twelve.
    for (const auto& subset : kernel) check(subset);
    const std::size_t limit = std::min<std::size_t>(12, kernel.size());
    for (std::size_t mask = 1; mask < (std::size_t(1) << limit) && report.factor == 0; ++mask) {
        if ((mask & (mask - 1)) == 0) continue;
        mpz_class subset = 0;
        for (std::size_t i = 0; i < limit; ++i) if (mask & (std::size_t(1) << i)) subset ^= kernel[i];
        check(subset);
    }
    return report;
}

void require(bool condition, const char* message)
{
    if (!condition) throw std::logic_error(message);
}

void self_test()
{
    // Compare all lifted/combined roots with exhaustive enumeration on small moduli.
    for (const auto& spec : std::vector<std::pair<mpz_class, unsigned long>>{{2, 5}, {3, 4}, {13, 2}, {19, 2}}) {
        const auto roots = lift_roots(73, spec.first, spec.second);
        require(roots.modulus == power(spec.first, spec.second), "Wrong lift modulus.");
        std::vector<mpz_class> expected;
        for (mpz_class r = 0; r < roots.modulus; ++r) {
            if (positive_mod(r * r - 73, roots.modulus) == 0) expected.push_back(r);
        }
        require(roots.roots == expected, "Hensel roots differ from exhaustive roots.");
    }
    const auto crt = combine(lift_roots(73, 2, 3), lift_roots(73, 3, 2));
    std::vector<mpz_class> expected;
    for (mpz_class r = 0; r < 72; ++r) if (positive_mod(r * r - 73, 72) == 0) expected.push_back(r);
    require(crt.modulus == 72 && crt.roots == expected && expected.size() == 8,
            "CRT must retain all eight roots modulo 72.");
    require(lift_roots(5, 2, 3).roots.empty(), "A nonsquare modulo 8 must not lift.");
    const std::vector<mpz_class> primes{2};
    relation first{10, 1}, second{13, 1};
    extract_cofactor(first, 77, primes); extract_cofactor(second, 77, primes);
    require(first.cofactor == 23 && second.cofactor == 23, "Cofactor extraction failed.");
    require(analyze(77, primes, {first}).nullity == 0,
            "Nonsquare c requires its own parity coordinate.");
    const auto useful = analyze(77, primes, {first, second});
    require(useful.rank == 1 && useful.nullity == 1 && useful.factor > 1,
            "Known common-cofactor pair should factor 77.");
    relation scaled{20, 4}; extract_cofactor(scaled, 77, primes);
    require(analyze(77, primes, {first, scaled}).factor == 0,
            "Scaled seed must give only a trivial dependency.");
    relation square{9, 1}; extract_cofactor(square, 77, primes);
    require(analyze(77, primes, {square}).factor > 1, "Square cofactor dependency failed.");
    relation nonunit_square{52, 35}; extract_cofactor(nonunit_square, 77, primes);
    require(nonunit_square.cofactor == 9 && analyze(77, primes, {nonunit_square}).factor > 1,
            "An odd number of square-cofactor rows must include sqrt(c) in Y.");
    // A complete bounded search verifies the discovery/locking contract and certificates.
    options opt; opt.target = opt.base_size;
    const auto result = search(opt);
    require(result.cofactor > 1 && result.relations.size() >= opt.target,
            "Default discovery search did not reach its target.");
    for (const auto& row : result.relations) {
        require(valid_relation(row, opt.number, result.primes), "Bad search relation.");
        require(row.cofactor == result.cofactor && row.g == result.g0 + row.modulus * row.k &&
                row.x == row.b + row.modulus * row.u, "Search certificate mismatch.");
        mpz_class m = 1;
        for (std::size_t i = 0; i < result.primes.size(); ++i) {
            require(row.exponents[i] >= row.forced[i], "Forced prime power was lost.");
            m *= power(result.primes[i], row.forced[i]);
        }
        require(m == row.modulus, "Forced modulus mismatch.");
    }
    opt.max_candidates = 1;
    const auto bounded = search(opt);
    require(bounded.candidates == 1 && bounded.relations.empty(), "Candidate budget ignored.");
    std::cout << "Fixed-cofactor invariant tests passed.\n";
}

template<class T> std::string joined(const std::vector<T>& values)
{
    std::ostringstream out;
    for (std::size_t i = 0; i < values.size(); ++i) { if (i) out << ';'; out << values[i]; }
    return out.str();
}

unsigned long small_option(const std::string& value)
{
    mpz_class parsed;
    if (parsed.set_str(value, 10) != 0 || parsed < 0 || !mpz_fits_ulong_p(parsed.get_mpz_t())) {
        throw std::invalid_argument("Invalid nonnegative integer: " + value);
    }
    return parsed.get_ui();
}
} // namespace

int main(int argc, char** argv)
{
    try {
        options opt;
        for (int i = 1; i < argc; ++i) {
            const std::string flag = argv[i];
            if (flag == "--self-test") { self_test(); return 0; }
            if (flag == "--help") {
                std::cout << "cofactor_experiment [--number N] [--base-size 5] [--relations |P|]\n"
                    "  [--lift-levels 3] [--radius 32] [--max-candidates 250000]\n"
                    "  [--max-multiplier 1000000000] [--csv path] [--self-test]\n"
                    "Discovers the first repeated c > 1, then locks c. Positive relations only.\n";
                return 0;
            }
            if (++i == argc) throw std::invalid_argument("Missing value for " + flag);
            const std::string value = argv[i];
            if (flag == "--number" || flag == "--max-multiplier") {
                mpz_class& destination = flag == "--number" ? opt.number : opt.max_multiplier;
                if (destination.set_str(value, 10) != 0) throw std::invalid_argument("Invalid integer.");
            } else if (flag == "--base-size") opt.base_size = small_option(value);
            else if (flag == "--relations") opt.target = small_option(value);
            else if (flag == "--lift-levels") opt.levels = small_option(value);
            else if (flag == "--radius") opt.radius = small_option(value);
            else if (flag == "--max-candidates") opt.max_candidates = small_option(value);
            else if (flag == "--csv") opt.csv = value;
            else throw std::invalid_argument("Unknown option: " + flag);
        }
        if (opt.number <= 1 || mpz_even_p(opt.number.get_mpz_t()) ||
            opt.base_size < 2 || opt.base_size > 12 || opt.levels < 1 || opt.levels > 8 ||
            opt.radius > 10000 || opt.max_candidates == 0 || opt.max_multiplier < 1) {
            throw std::invalid_argument("Require odd N > 1, 2..12 primes, 1..8 lift levels, "
                                        "radius <= 10000, and positive search bounds.");
        }
        if (!opt.target) opt.target = opt.base_size;
        if (opt.target < 2 || opt.target > 1024) throw std::invalid_argument("Require 2..1024 relations.");
        const auto start = std::chrono::steady_clock::now();
        const auto result = search(opt);
        const auto report = analyze(opt.number, result.primes, result.relations);
        const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
        std::cout << "N=" << opt.number << " g0=" << result.g0 << " P={" << joined(result.primes) << "}\n"
                  << "c=" << result.cofactor << " (0 means not discovered), relations=" << result.relations.size()
                  << '/' << opt.target << ", stop=" << result.stopped << '\n'
                  << "candidates=" << result.candidates << " locked_at=" << result.discovered_at
                  << " profiles=" << result.profiles << " duplicates=" << result.duplicates
                  << " skipped_c1=" << result.unit_cofactors << " seconds=" << seconds << '\n';
        for (std::size_t i = 0; i < result.relations.size(); ++i) {
            const auto& row = result.relations[i];
            std::cout << i + 1 << ": A=" << row.x << " g=" << row.g << " exponents=["
                      << joined(row.exponents) << "] c=" << row.cofactor << " M=" << row.modulus << '\n';
        }
        std::cout << "parity_rank=" << report.rank << " nullity=" << report.nullity
                  << " dependencies_tested=" << report.tested << '\n';
        if (report.factor > 1) {
            std::cout << "Dependency factor: " << report.factor << " * " << mpz_class(opt.number / report.factor)
                      << " = " << opt.number << "; X=" << report.x << " Y=" << report.y << "; rows=";
            for (std::size_t j = 0; j < result.relations.size(); ++j) {
                if (mpz_tstbit(report.subset.get_mpz_t(), static_cast<mp_bitcnt_t>(j))) std::cout << j + 1 << ' ';
            }
            std::cout << '\n';
        } else std::cout << "No nontrivial factor from the tested dependencies.\n";
        if (result.direct_factor > 1) std::cout << "Factor noticed while selecting P: " << result.direct_factor << '\n';
        if (result.cofactor > 1) {
            const mpz_class divisor = gcd(result.cofactor, opt.number);
            if (divisor > 1 && divisor < opt.number) std::cout << "Direct gcd(c,N) factor: " << divisor << '\n';
        }
        if (!opt.csv.empty()) {
            std::ofstream csv(opt.csv);
            if (!csv) throw std::runtime_error("Cannot open CSV output.");
            csv << "N,g0,primes,A,g,M,b,u,k,forced_exponents,actual_exponents,c\n";
            for (const auto& row : result.relations) {
                csv << opt.number << ',' << result.g0 << ',' << joined(result.primes) << ',' << row.x << ','
                    << row.g << ',' << row.modulus << ',' << row.b << ',' << row.u << ',' << row.k << ','
                    << joined(row.forced) << ',' << joined(row.exponents) << ',' << row.cofactor << '\n';
            }
            if (!csv) throw std::runtime_error("Failed to write CSV output.");
        }
        return result.relations.size() >= opt.target ? 0 : 2;
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }
}
