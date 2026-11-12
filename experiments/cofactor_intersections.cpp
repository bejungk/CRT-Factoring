#include "math.hpp"

#include <algorithm>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <optional>
#include <queue>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>
#include <omp.h>

namespace {
using factoring::positive_mod;

// Only the GMP C DLL is linked in this Windows build.
std::ostream& operator<<(std::ostream& out, const mpz_class& value)
{
    return out << value.get_str();
}

struct options {
    mpz_class n = 77, g1 = 0, g2 = 0;
    unsigned long ell_min = 3, ell_max = 12;
    unsigned long max_candidates = 100000;
    bool keep_searching = false;
    std::string output;
};

struct branch {
    mpz_class a, modulus, b, u, c;
    unsigned long ell;
};

struct match {
    unsigned long ell;
    mpz_class a1, a2, b1, b2, u1, u2, c, gcd_minus, gcd_plus;
    bool useful() const { return gcd_minus > 1 || gcd_plus > 1; }
};

struct result {
    unsigned long candidates = 0, exact_branches = 0, trivial = 0;
    mpz_class last_a = 0;
    std::vector<match> matches;
    bool found_factor = false;
};

mpz_class power2(unsigned long ell) { return mpz_class(1) << ell; }

mpz_class floor_div(const mpz_class& a, const mpz_class& b)
{
    mpz_class q;
    mpz_fdiv_q(q.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t());
    return q;
}

mpz_class integer_sqrt(const mpz_class& value)
{
    mpz_class root;
    mpz_sqrt(root.get_mpz_t(), value.get_mpz_t());
    return root;
}

mpz_class gcd(const mpz_class& a, const mpz_class& b)
{
    mpz_class divisor;
    mpz_gcd(divisor.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t());
    return divisor;
}

std::vector<mpz_class> lift_two_roots(const mpz_class& d, unsigned long ell)
{
    if (ell < 3 || positive_mod(d, mpz_class(8)) != 1) {
        throw std::invalid_argument("Require ell >= 3 and g*N = 1 (mod 8).");
    }
    mpz_class modulus = 8;
    std::vector<mpz_class> roots{1, 3, 5, 7};
    for (unsigned long e = 3; e < ell; ++e) {
        const mpz_class next_modulus = 2 * modulus;
        std::vector<mpz_class> next;
        for (const auto& r : roots) {
            for (unsigned long digit = 0; digit < 2; ++digit) {
                const mpz_class child = r + digit * modulus;
                if (positive_mod(child * child - d, next_modulus) == 0) next.push_back(child);
            }
        }
        roots = std::move(next);
        modulus = next_modulus;
    }
    std::sort(roots.begin(), roots.end());
    return roots;
}

// Check one C1(u1) directly: equality requires A2^2=g2*N+2^(ell+2)*c.
// This avoids a Cartesian search or a hash table of the second polynomial.
std::optional<match> intersection(const options& opt, const branch& left)
{
    const mpz_class m2 = 4 * left.modulus;
    const mpz_class radicand = opt.g2 * opt.n + m2 * left.c;
    const mpz_class a2 = integer_sqrt(radicand);
    if (a2 * a2 != radicand) return std::nullopt;
    const mpz_class b2 = positive_mod(a2, m2);
    const mpz_class u2 = (a2 - b2) / m2;
    // Independently check the polynomial value and the difference-of-squares
    // identity before recording a match.
    const mpz_class c2 = m2 * u2 * u2 + 2 * b2 * u2 + (b2 * b2 - opt.g2 * opt.n) / m2;
    if (left.c <= 0 || mpz_even_p(left.c.get_mpz_t()) || c2 != left.c ||
        left.a * left.a - opt.g1 * opt.n != left.modulus * left.c ||
        a2 * a2 - opt.g2 * opt.n != m2 * left.c ||
        (a2 - 2 * left.a) * (a2 + 2 * left.a) != (opt.g2 - 4 * opt.g1) * opt.n) {
        throw std::logic_error("Intersection certificate failed.");
    }
    mpz_class minus = gcd(a2 - 2 * left.a, opt.n);
    mpz_class plus = gcd(a2 + 2 * left.a, opt.n);
    // Zero marks a trivial gcd, so useful() never counts gcd=N as a factor.
    if (minus == opt.n || minus == 1) minus = 0;
    if (plus == opt.n || plus == 1) plus = 0;
    return match{left.ell, left.a, a2, left.b, b2, left.u, u2, left.c, minus, plus};
}

result search(const options& opt)
{
    const mpz_class d1 = opt.g1 * opt.n;
    const mpz_class sqrt_d1 = integer_sqrt(d1);
    auto later = [](const branch& a, const branch& b) {
        return a.a > b.a; // Smallest positive A is visited first across all ell.
    };
    std::priority_queue<branch, std::vector<branch>, decltype(later)> pending(later);
    result out;
    for (unsigned long ell = opt.ell_min; ell <= opt.ell_max; ++ell) {
        const mpz_class m = power2(ell);
        for (const auto& b : lift_two_roots(d1, ell)) {
            const mpz_class u = floor_div(sqrt_d1 - b, m) + 1;
            const mpz_class a = b + m * u;
            const mpz_class c = (a * a - d1) / m;
            // At p=2 this parity is constant along the entire branch:
            // C(u+1)-C(u)=2*A+M is even. Two of the four branches qualify.
            if (mpz_even_p(c.get_mpz_t())) continue;
            pending.push({a, m, b, u, c, ell});
            ++out.exact_branches;
        }
    }
    while (!pending.empty() && out.candidates < opt.max_candidates) {
        branch left = pending.top();
        pending.pop();
        ++out.candidates;
        out.last_a = left.a;
        if (const auto found = intersection(opt, left)) {
            out.matches.push_back(*found);
            if (found->useful()) {
                out.found_factor = true;
                if (!opt.keep_searching) break;
            } else ++out.trivial;
        }
        left.c += 2 * left.a + left.modulus;
        left.a += left.modulus;
        ++left.u;
        pending.push(std::move(left));
    }
    return out;
}

void require(bool ok, const char* message)
{
    if (!ok) throw std::logic_error(message);
}

using match_key = std::tuple<unsigned long, mpz_class, mpz_class, mpz_class>;
match_key key(const match& value) { return {value.ell, value.a1, value.a2, value.c}; }

void self_test()
{
    for (unsigned long ell = 3; ell <= 8; ++ell) {
        const mpz_class m = power2(ell);
        std::vector<mpz_class> brute;
        for (mpz_class r = 0; r < m; ++r) if (positive_mod(r * r - 385, m) == 0) brute.push_back(r);
        require(lift_two_roots(385, ell) == brute && brute.size() == 4, "Incorrect 2-adic roots.");
    }
    options opt;
    opt.g1 = opt.g2 = 5;
    opt.ell_min = opt.ell_max = 3;
    const auto small = search(opt);
    require(small.found_factor && small.candidates == 3 && small.exact_branches == 2,
            "Did not find the earliest useful intersection across both exact branches.");
    const auto& pair = small.matches.front();
    require(pair.a1 == 29 && pair.a2 == 47 && pair.c == 57 &&
            pair.b1 == 5 && pair.b2 == 15 && pair.u1 == 3 && pair.u2 == 1 &&
            pair.gcd_minus == 11 && pair.gcd_plus == 7,
            "Known polynomial intersection differs.");
    const auto earlier_example = intersection(opt, branch{43, 8, 3, 5, 183, 3});
    require(earlier_example && earlier_example->a2 == 79 &&
            earlier_example->gcd_minus == 7 && earlier_example->gcd_plus == 11,
            "The original c=183 example failed.");
    branch trivial{59, 8, 3, 7, 387, 3};
    const auto trivial_match = intersection(opt, trivial);
    require(trivial_match && !trivial_match->useful(), "Trivial square congruence counted as a factor.");
    opt.max_candidates = 1;
    require(search(opt).matches.empty(), "Search ignored candidate budget.");

    // Compare the heap enumeration to a plain exhaustive odd-A scan. Include
    // different g1,g2 to exercise the general equation, not only equal g.
    for (unsigned long g2 : {5UL, 13UL}) {
        opt.g2 = g2; opt.ell_max = 6; opt.max_candidates = 200; opt.keep_searching = true;
        const auto scanned = search(opt);
        std::vector<match_key> actual, expected;
        for (const auto& value : scanned.matches) actual.push_back(key(value));
        for (mpz_class a = 21; a <= scanned.last_a; a += 2) {
            mpz_class c = a * a - opt.g1 * opt.n;
            unsigned long ell = 0;
            while (mpz_even_p(c.get_mpz_t())) { c /= 2; ++ell; }
            if (ell < opt.ell_min || ell > opt.ell_max) continue;
            const mpz_class target = 4 * a * a + (opt.g2 - 4 * opt.g1) * opt.n;
            const mpz_class a2 = integer_sqrt(target);
            if (a2 * a2 == target) expected.emplace_back(ell, a, a2, c);
        }
        require(actual == expected, "Intersection search differs from exhaustive scan.");
    }
    std::cout << "Cofactor-intersection invariant tests passed.\n";
}

unsigned long unsigned_option(const std::string& text)
{
    mpz_class value;
    if (value.set_str(text, 10) || value < 0 || !mpz_fits_ulong_p(value.get_mpz_t())) {
        throw std::invalid_argument("Invalid nonnegative integer: " + text);
    }
    return value.get_ui();
}

void print_result(std::ostream& out, const options& opt, const result& found, double seconds)
{
    out << "N=" << opt.n << " P={2} g1=" << opt.g1 << " g2=" << opt.g2
        << " ell1=" << opt.ell_min << ".." << opt.ell_max << " ell2=ell1+2\n"
        << "candidate_limit=" << opt.max_candidates << " exact_root_branches=" << found.exact_branches
        << " candidates=" << found.candidates << " last_A1=" << found.last_a
        << " matches=" << found.matches.size() << " trivial_matches=" << found.trivial
        << " seconds=" << seconds << '\n';
    for (const auto& row : found.matches) {
        out << "A1=" << row.a1 << " ell1=" << row.ell << " b1=" << row.b1 << " u1=" << row.u1
            << " A2=" << row.a2 << " ell2=" << row.ell + 2 << " b2=" << row.b2 << " u2=" << row.u2
            << " c=" << row.c << '\n';
        out << "  " << row.a1 << "^2 - " << opt.g1 << "*" << opt.n << " = 2^" << row.ell << "*" << row.c << '\n'
            << "  " << row.a2 << "^2 - " << opt.g2 << "*" << opt.n << " = 2^" << row.ell + 2 << "*" << row.c << '\n';
        if (row.useful()) {
            const mpz_class factor = row.gcd_minus > 1 ? row.gcd_minus : row.gcd_plus;
            out << "  nontrivial factor: " << factor << " * " << mpz_class(opt.n / factor)
                << "; gcd_minus=" << row.gcd_minus << " gcd_plus=" << row.gcd_plus << '\n';
        } else out << "  trivial square congruence\n";
    }
    out << "status=" << (found.found_factor ? "nontrivial factor found" : "no nontrivial factor within candidate budget") << '\n';
}

// Exact trial-division primality is inexpensive for this deliberately bounded
// test set (all factors are below 1.1e9). It avoids labeling probable primes as
// proven semiprimes. These routines construct inputs, not search candidates.
bool is_small_prime(unsigned long n)
{
    if (n < 2) return false;
    if (n % 2 == 0) return n == 2;
    for (unsigned long d = 3; d <= n / d; d += 2) if (n % d == 0) return false;
    return true;
}

unsigned long next_small_prime(unsigned long n)
{
    while (!is_small_prime(n)) ++n;
    return n;
}

struct batch_case {
    unsigned long p, q, target_ratio;
    options input;
    result found;
    double seconds = 0;
    int worker = -1;
    std::string error;
};

std::vector<batch_case> separated_cases(const options& bounds)
{
    const unsigned long scales[]{1000, 10000, 100000, 1000000, 10000000};
    const unsigned long ratios[]{3, 10, 30, 100};
    std::vector<batch_case> cases;
    for (unsigned long s = 0; s < 5; ++s) {
        for (unsigned long r = 0; r < 4; ++r) {
            const unsigned long p = next_small_prime(scales[s] + 37 * (r + 1) + 101 * s);
            const unsigned long q = next_small_prime(ratios[r] * p + 17 + 2 * s);
            options input = bounds;
            input.n = mpz_class(p) * q;
            input.g1 = input.g2 = positive_mod(input.n, mpz_class(8));
            input.output.clear();
            cases.push_back({p, q, ratios[r], std::move(input), {}});
        }
    }
    return cases;
}

void parallel_self_test()
{
    options bounds;
    bounds.max_candidates = 2000;
    auto cases = separated_cases(bounds);
    std::vector<result> serial(4), parallel(4);
    for (int i = 0; i < 4; ++i) serial[i] = search(cases[i].input);
    std::vector<std::string> errors(4);
    omp_set_dynamic(0);
#pragma omp parallel for num_threads(4) schedule(dynamic, 1)
    for (int i = 0; i < 4; ++i) {
        try { parallel[i] = search(cases[i].input); }
        catch (const std::exception& error) { errors[i] = error.what(); }
    }
    for (int i = 0; i < 4; ++i) {
        require(errors[i].empty(), "Parallel search threw an exception.");
        const auto& a = serial[i]; const auto& b = parallel[i];
        require(a.candidates == b.candidates && a.last_a == b.last_a &&
                a.found_factor == b.found_factor && a.matches.size() == b.matches.size(),
                "OpenMP search differs from serial search.");
        for (std::size_t j = 0; j < a.matches.size(); ++j) {
            require(key(a.matches[j]) == key(b.matches[j]) &&
                    a.matches[j].gcd_minus == b.matches[j].gcd_minus &&
                    a.matches[j].gcd_plus == b.matches[j].gcd_plus,
                    "OpenMP produced different intersection certificates.");
        }
    }
    for (const auto& item : cases) {
        require(is_small_prime(item.p) && is_small_prime(item.q) && item.q > item.p &&
                item.input.n == mpz_class(item.p) * item.q, "Invalid separated-factor input.");
    }
    std::cout << "OpenMP/serial reproducibility and 20 semiprime input checks passed.\n";
}

int run_separated_batch(const options& bounds, unsigned long requested_threads,
                        const std::string& directory)
{
    auto cases = separated_cases(bounds);
    const int threads = requested_threads ? static_cast<int>(requested_threads) :
        std::min(20, std::min(omp_get_max_threads(), omp_get_num_procs()));
    const std::filesystem::path destination(directory);
    std::filesystem::create_directories(destination);
    omp_set_dynamic(0);
    int actual_threads = 0;
    std::cout << "Local OpenMP batch: 20 experiments; processors=" << omp_get_num_procs()
              << "; requested_workers=" << threads << "; candidates_per_case="
              << bounds.max_candidates << '\n' << std::flush;
    const auto start = std::chrono::steady_clock::now();
#pragma omp parallel num_threads(threads)
    {
#pragma omp single
        actual_threads = omp_get_num_threads();
#pragma omp for schedule(dynamic, 1)
        for (int i = 0; i < 20; ++i) {
            auto& item = cases[i];
            item.worker = omp_get_thread_num();
            const auto case_start = std::chrono::steady_clock::now();
            try {
                // Only the ordinary search options, including N, are passed in.
                // Neither p nor q is accessible to the search.
                item.found = search(item.input);
                for (const auto& pair : item.found.matches) {
                    for (const auto& divisor : {pair.gcd_minus, pair.gcd_plus}) {
                        if (divisor > 1 && divisor != item.p && divisor != item.q) {
                            throw std::logic_error("Reported factor disagrees with verified input primes.");
                        }
                    }
                }
            } catch (const std::exception& error) { item.error = error.what(); }
            item.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - case_start).count();
#pragma omp critical(cofactor_batch_progress)
            {
                std::cout << "case=" << i + 1 << " worker=" << item.worker
                          << " candidates=" << item.found.candidates
                          << " status=" << (!item.error.empty() ? "error" :
                              item.found.found_factor ? "factored" : "budget exhausted")
                          << '\n' << std::flush;
            }
        }
    }
    const double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    unsigned long successes = 0, failures = 0;
    unsigned long long total_candidates = 0;
    double sum_case_seconds = 0;
    for (const auto& item : cases) {
        successes += item.found.found_factor ? 1 : 0;
        failures += item.error.empty() ? 0 : 1;
        total_candidates += item.found.candidates;
        sum_case_seconds += item.seconds;
    }
    std::ofstream summary(destination / "summary.md");
    if (!summary) throw std::runtime_error("Cannot write batch summary.");
    summary << "# Twenty separated-factor semiprimes: local OpenMP run\n\n"
            << "Executed locally with " << actual_threads << " OpenMP workers (runtime reports "
            << omp_get_num_procs() << " processors). OpenMP macro: " << _OPENMP << ".\n\n"
            << "Budget: " << bounds.max_candidates << " candidates per case; ell1=" << bounds.ell_min
            << ".." << bounds.ell_max << "; ell2=ell1+2; g1=g2 is the small multiplier making g*N=1 mod 8.\n\n"
            << "Results: " << successes << "/20 factored; " << 20 - successes - failures
            << " exhausted the budget; " << failures << " errors. Total candidates: " << total_candidates << ".\n\n"
            << std::fixed << std::setprecision(6)
            << "Batch wall time: " << elapsed << " seconds. Sum of per-case elapsed times: "
            << sum_case_seconds << " seconds. Per-case times include contention during this parallel run; "
            << "the sum is not a measured serial baseline or CPU time.\n\n"
            << "Inputs were fixed before searching: smaller-prime seed = scale + 37*(ratio_index+1) + 101*scale_index, "
            << "for scales 10^3..10^7 and ratios 3,10,30,100; p is the first prime at or above that seed. "
            << "q is the first prime at or above ratio*p + 17 + 2*scale_index. "
            << "Indices start at zero. All p and q are proven prime by exhaustive trial division "
            << "through their square roots (all factors <1.1e9). Search receives N, not its factors.\n\n"
            << "| Case | p | q | q/p | N | bits | g | candidates | outcome | seconds | worker |\n"
            << "|---:|---:|---:|---:|---:|---:|---:|---:|---|---:|---:|\n";
    for (std::size_t i = 0; i < cases.size(); ++i) {
        const auto& item = cases[i];
        const std::string name = "case_" + std::to_string(i + 1) + ".txt";
        std::ofstream detail(destination / name);
        if (!detail) throw std::runtime_error("Cannot write batch case report.");
        print_result(detail, item.input, item.found, item.seconds);
        detail << "expected_p=" << item.p << " expected_q=" << item.q
               << " target_ratio=" << item.target_ratio << " worker=" << item.worker << '\n';
        if (!item.error.empty()) detail << "error=" << item.error << '\n';
        if (!detail) throw std::runtime_error("Failed to write case report.");
        summary << "| [" << i + 1 << "](" << name << ") | " << item.p << " | " << item.q
                << " | " << static_cast<double>(item.q) / item.p << " | " << item.input.n
                << " | " << mpz_sizeinbase(item.input.n.get_mpz_t(), 2) << " | " << item.input.g1
                << " | " << item.found.candidates << " | " << (!item.error.empty() ? "ERROR" :
                    item.found.found_factor ? "factored" : "budget exhausted")
                << " | " << item.seconds << " | " << item.worker << " |\n";
    }
    summary << "\nThis is a structured exploratory sample, not a representative random benchmark. "
            << "Search-budget exhaustion is not a proof of absence of intersections. "
            << "The algorithm is still a constrained difference-of-squares scan.\n";
    if (!summary) throw std::runtime_error("Failed to write summary.");
    std::cout << "Batch finished: " << successes << "/20 factored, " << failures << " errors, "
              << total_candidates << " candidates, " << actual_threads << " actual workers, "
              << elapsed << " seconds wall time. Reports: " << destination.string() << '\n';
    return failures ? 1 : 0; // An exhausted search budget is a completed experiment.
}
} // namespace

int main(int argc, char** argv)
{
    try {
        options opt;
        bool batch = false, custom_number_or_multiplier = false;
        unsigned long threads = 0;
        std::string batch_directory = "experiments/results/separated_openmp";
        for (int i = 1; i < argc; ++i) {
            const std::string flag = argv[i];
            if (flag == "--self-test") { self_test(); return 0; }
            if (flag == "--parallel-self-test") { parallel_self_test(); return 0; }
            if (flag == "--batch-separated") { batch = true; continue; }
            if (flag == "--keep-searching") { opt.keep_searching = true; continue; }
            if (flag == "--help") {
                std::cout << "cofactor_intersections [--number 77] [--g1 G] [--g2 G]\n"
                    "  [--ell-min 3] [--ell-max 12] [--max-candidates 100000]\n"
                    "  [--keep-searching] [--output path] [--self-test]\n"
                    "  [--batch-separated --threads 0 --batch-directory path]\n"
                    "Batch mode runs 20 fixed separated-factor inputs using OpenMP; threads=0 selects automatically.\n"
                    "Defaults g1=g2 to the small multiplier making g*N=1 mod 8.\n"
                    "Searches positive odd equal cofactors with ell2=ell1+2.\n";
                return 0;
            }
            if (++i == argc) throw std::invalid_argument("Missing value for " + flag);
            const std::string value = argv[i];
            if (flag == "--number" || flag == "--g1" || flag == "--g2") {
                custom_number_or_multiplier = true;
                mpz_class& destination = flag == "--number" ? opt.n : flag == "--g1" ? opt.g1 : opt.g2;
                if (destination.set_str(value, 10) || destination <= 0) throw std::invalid_argument("Expected a positive integer.");
            } else if (flag == "--ell-min") opt.ell_min = unsigned_option(value);
            else if (flag == "--ell-max") opt.ell_max = unsigned_option(value);
            else if (flag == "--max-candidates") opt.max_candidates = unsigned_option(value);
            else if (flag == "--output") opt.output = value;
            else if (flag == "--threads") threads = unsigned_option(value);
            else if (flag == "--batch-directory") batch_directory = value;
            else throw std::invalid_argument("Unknown option: " + flag);
        }
        if (opt.n <= 1 || mpz_even_p(opt.n.get_mpz_t()) || opt.ell_min < 3 ||
            opt.ell_max < opt.ell_min || opt.ell_max > 256 || opt.max_candidates == 0) {
            throw std::invalid_argument("Require odd N>1, 3<=ell-min<=ell-max<=256 and a positive candidate budget.");
        }
        if (threads > 20) throw std::invalid_argument("Use 0..20 workers for this 20-case batch.");
        if (batch) {
            if (custom_number_or_multiplier || !opt.output.empty() || opt.keep_searching) {
                throw std::invalid_argument("Batch fixes inputs and stops at the first useful pair; use --batch-directory for reports.");
            }
            return run_separated_batch(opt, threads, batch_directory);
        }
        const mpz_class default_g = positive_mod(opt.n, mpz_class(8));
        if (opt.g1 == 0) opt.g1 = default_g;
        if (opt.g2 == 0) opt.g2 = default_g;
        if (positive_mod(opt.g1 * opt.n, mpz_class(8)) != 1 ||
            positive_mod(opt.g2 * opt.n, mpz_class(8)) != 1) {
            throw std::invalid_argument("Both multipliers must satisfy g_i*N=1 (mod 8).");
        }
        const auto start = std::chrono::steady_clock::now();
        const auto found = search(opt);
        const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
        print_result(std::cout, opt, found, seconds);
        if (!opt.output.empty()) {
            std::ofstream output(opt.output);
            if (!output) throw std::runtime_error("Cannot open output file.");
            print_result(output, opt, found, seconds);
            if (!output) throw std::runtime_error("Could not write output file.");
        }
        return found.found_factor ? 0 : 2;
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }
}
