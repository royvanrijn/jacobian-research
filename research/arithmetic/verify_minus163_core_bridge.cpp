#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>

using U = std::uint64_t;
using Wide = __uint128_t;
constexpr int p = 163;
static_assert((U{1} << 43) < 341550071728321ULL);
constexpr int profile_count = 729;
constexpr std::array<U, 6> tracked{3, 5, 7, 11, 13, 17};
using Cutoffs = std::array<std::array<U, 2>, profile_count>;


U modular_power(U a, U exponent, U modulus) {
    U result = 1;
    while (exponent) {
        if (exponent & 1)
            result = Wide(result) * a % modulus;
        exponent >>= 1;
        a = Wide(a) * a % modulus;
    }
    return result;
}

bool is_large_prime(U n) {
    // Jaeschke's bound for the first seven prime bases exceeds 2^43.
    constexpr std::array<U, 7> bases{2, 3, 5, 7, 11, 13, 17};
    if (n <= bases.back() || !(n & 1) || n >= (U{1} << 43))
        throw std::logic_error("large primality test outside its domain");
    const int twos = std::countr_zero(n - 1);
    const U odd = (n - 1) >> twos;
    for (U base : bases) {
        U power = modular_power(base, odd, n);
        if (power == 1 || power == n - 1)
            continue;
        for (int step = 1; step < twos && power != n - 1; ++step)
            power = Wide(power) * power % n;
        if (power != n - 1)
            return false;
    }
    return true;
}

struct WitnessBank {
    const std::array<bool, p>& squares;
    const std::vector<bool>& prime_flags;
    std::vector<unsigned> small_primes;
    std::vector<U> summands;

    explicit WitnessBank(const std::array<bool, p>& squares,
                         const std::vector<bool>& prime_flags)
        : squares(squares), prime_flags(prime_flags) {
        for (unsigned q = 3; q <= 163; q += 2)
            if (prime_flags[q])
                small_primes.push_back(q);
        constexpr U search_bound = 100000;
        std::vector<bool> split(search_bound + 1, true);
        for (U q = 2; q <= search_bound; ++q)
            if (prime_flags[q] && !squares[q % p])
                for (U multiple = q; multiple <= search_bound; multiple += q)
                    split[multiple] = false;
        for (U a = 1; a <= search_bound; a += 2)
            if (split[a])
                summands.push_back(a);
    }

    bool certify_core(U b) const {
        if (!squares[b % p])
            return false;
        for (U q : small_primes) {
            if (b % q)
                continue;
            if (!squares[q % p])
                return false;
            do {
                b /= q;
            } while (b % q == 0);
        }
        if (b == 1)
            return true;
        if (b < prime_flags.size())
            return prime_flags[b];
        // Division by split primes has preserved the initial square residue.
        return is_large_prime(b);
    }

    U witness(U n) const {
        for (U a : summands) {
            if (a >= n)
                break;
            if (certify_core(n - a))
                return a;
        }
        return 0;
    }
};

struct Runner {
    const std::vector<unsigned>& primes;
    const Cutoffs& cutoffs;
    const WitnessBank& bank;
    U scale, lower, upper;
    int t;
    U leaves = 0;
    U calls = 0;
    U largest_witness = 0;
    U first_failure = 0;
    U witness_sum = 0;

    void odd_factors(U m, U tau, U sigma, bool square, int code, std::size_t next) {
        ++calls;
        const U n = m << t;
        const Wide score = Wide(t + 3) * tau * m * scale;
        const auto row = cutoffs[code];
        if (n >= lower && score >= Wide(row[square]) * sigma) {
            ++leaves;
            const U a = bank.witness(n);
            if (!a && (!first_failure || n < first_failure))
                first_failure = n;
            largest_witness = std::max(largest_witness, a);
            witness_sum += a;
        }
        const Wide threshold = Wide(std::min(row[0],row[1])) * sigma;
        const auto need = static_cast<U>((threshold + score - 1) / score);
        const int degree = std::max(1, static_cast<int>(std::bit_width(need - 1)));
        const U ceiling = upper / n;
        if (next == primes.size() || ceiling < primes[next])
            return;
        for (std::size_t j = next; j < primes.size(); ++j) {
            const U q = primes[j];
            // tau(v) <= 2^Omega(v). Any continuation needing at least this
            // many prime factors has size at least q^degree.
            U quotient = ceiling;
            for (int e = 0; e < degree && quotient; ++e)
                quotient /= q;
            if (!quotient)
                break;
            U power = 1;
            U series = 1;
            for (U a = 1; power <= ceiling / q; ++a) {
                power *= q;
                series += power;
                odd_factors(m * power, tau * (a+1), sigma * series,
                            square != bool(a & 1U), code, j+1);
            }
        }
    }

    void small_factors(std::size_t i, U m, U tau, U sigma, bool square, int code) {
        if (i == tracked.size()) {
            odd_factors(m, tau, sigma, square, code, 0);
            return;
        }
        U power = 1;
        U series = 1;
        const U q = tracked[i];
        const U ceiling = upper / (m << t);
        for (U a = 0; ; ++a) {
            small_factors(i+1, m * power, tau * (a+1), sigma * series,
                          square != bool(a & 1U), 3*code + std::min(a,U{2}));
            if (power > ceiling / q)
                break;
            power *= q;
            series += power;
        }
    }
};

int main(int argc, char** argv) {
    if (argc != 2)
        throw std::invalid_argument("provide the arithmetic threshold table");
    std::ifstream input(argv[1]);
    int first, last, profiles;
    U scale, cap;
    input >> first >> last >> scale >> cap >> profiles;
    if (!input || first != 34 || last != 42 || profiles != profile_count ||
        cap != (1UL << 25) || scale != (1UL << 20))
        throw std::runtime_error("invalid threshold header");
    std::vector<Cutoffs> all_cuts(last - first + 1);
    for (auto& bin : all_cuts)
        for (auto& thresholds : bin) {
            input >> thresholds[0] >> thresholds[1];
            if (!input || !std::min(thresholds[0], thresholds[1]))
                throw std::runtime_error("invalid threshold row");
        }

    std::vector<bool> is_prime(cap+1, true);
    is_prime[0] = is_prime[1] = false;
    for (U q = 2; q*q <= cap; ++q)
        if (is_prime[q])
            for (U multiple = q*q; multiple <= cap; multiple += q)
                is_prime[multiple] = false;
    std::array<bool,p> square{};
    for (int a = 1; a < p; ++a)
        square[a*a % p] = true;
    std::vector<unsigned> primes;
    for (U q = tracked.back()+2; q <= cap; q+=2)
        if (is_prime[q] && q != p && !square[q % p])
            primes.push_back(q);

    const WitnessBank bank(square, is_prime);
    if (is_large_prime(3215031751) || !is_large_prime(4294967291))
        throw std::logic_error("primality controls failed");

    struct Job {
        int k, t;
        U leaves = 0, calls = 0, largest_witness = 0, first_failure = 0, witness_sum = 0;
    };
    std::vector<Job> jobs;
    for (int k = first; k <= last; ++k)
        for (int t = 1; t <= k; ++t)
            jobs.push_back({k, t});
    #pragma omp parallel for schedule(dynamic)
    for (std::size_t j = 0; j < jobs.size(); ++j) {
        auto& job = jobs[j];
        const int k = job.k;
        Runner runner{primes, all_cuts[k-first], bank, scale, U{1} << k, (U{2} << k)-1, job.t};
        runner.small_factors(0, 1, 1, 1, job.t % 2 == 0, 0);
        job.leaves = runner.leaves;
        job.calls = runner.calls;
        job.largest_witness = runner.largest_witness;
        job.first_failure = runner.first_failure;
        job.witness_sum = runner.witness_sum;
    }
    bool ok = true;
    for (int k = first; k <= last; ++k) {
        Job total{k, 0};
        for (auto job : jobs) {
            if (job.k != k)
                continue;
            total.leaves += job.leaves;
            total.calls += job.calls;
            total.largest_witness = std::max(total.largest_witness, job.largest_witness);
            total.witness_sum += job.witness_sum;
            if (job.first_failure && (!total.first_failure || job.first_failure < total.first_failure))
                total.first_failure = job.first_failure;
        }
        std::cout << "{\"k\":" << k << ",\"targets\":" << total.leaves
                  << ",\"visited_prefixes\":" << total.calls
                  << ",\"largest_smaller_summand\":" << total.largest_witness
                  << ",\"summand_checksum\":" << total.witness_sum
                  << ",\"first_failure\":" << total.first_failure << "}\n";
        ok &= total.first_failure == 0;
    }
    return ok ? 0 : 1;
}
