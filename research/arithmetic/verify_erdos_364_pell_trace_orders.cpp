#include <array>
#include <bit>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <vector>

#include <gmpxx.h>

using Integer = std::uint64_t;
using Matrix = std::array<mpz_class, 4>;
static_assert(sizeof(unsigned long) >= sizeof(Integer));

std::vector<Integer> continued_fraction(Integer radicand) {
    Integer root = 0;
    while (root + 1 <= radicand / (root + 1)) {
        ++root;
    }
    if (root * root == radicand) {
        throw std::runtime_error("square radicand");
    }
    Integer m = 0, denominator = 1, term = root;
    std::vector<Integer> coefficients{root};
    do {
        m = denominator * term - m;
        denominator = (radicand - m * m) / denominator;
        term = (root + m) / denominator;
        coefficients.push_back(term);
    } while (term != 2 * root);
    coefficients.pop_back();
    if (coefficients.size() % 2 != 0) {
        throw std::runtime_error("expected an even Pell period");
    }
    return coefficients;
}

Matrix convergents(const std::vector<Integer>& coefficients, std::size_t lo, std::size_t hi) {
    if (hi - lo <= 24) {
        Matrix matrix{1, 0, 0, 1};
        for (std::size_t i = lo; i < hi; ++i) {
            Integer term = coefficients[i];
            Matrix next{matrix[0] * term + matrix[1], matrix[0],
                        matrix[2] * term + matrix[3], matrix[2]};
            matrix.swap(next);
        }
        return matrix;
    }
    auto mid = lo + (hi - lo) / 2;
    auto left = convergents(coefficients, lo, mid);
    auto right = convergents(coefficients, mid, hi);
    return {left[0] * right[0] + left[1] * right[2],
            left[0] * right[1] + left[1] * right[3],
            left[2] * right[0] + left[3] * right[2],
            left[2] * right[1] + left[3] * right[3]};
}

std::vector<unsigned> least_prime_factors(unsigned limit) {
    std::vector<unsigned> factors(limit + 1);
    for (unsigned p = 2; p <= limit; ++p) {
        if (factors[p] != 0) {
            continue;
        }
        for (unsigned n = p; n <= limit; n += p) {
            if (factors[n] == 0) {
                factors[n] = p;
            }
        }
    }
    return factors;
}

std::vector<Integer> squarefree_factors(Integer radicand) {
    std::vector<Integer> factors;
    for (Integer p = 2; p <= radicand / p; ++p) {
        if (radicand % p == 0) {
            factors.push_back(p);
            radicand /= p;
            if (radicand % p == 0) {
                throw std::runtime_error("radicand is not squarefree");
            }
        }
    }
    if (radicand > 1) {
        factors.push_back(radicand);
    }
    return factors;
}

Integer multiply_mod(Integer a, Integer b, Integer m) {
    return static_cast<__uint128_t>(a) * b % m;
}

Integer power_mod(Integer a, Integer e, Integer m) {
    Integer result = 1;
    while (e != 0) {
        if (e % 2 != 0) {
            result = multiply_mod(result, a, m);
        }
        e /= 2;
        a = multiply_mod(a, a, m);
    }
    return result;
}

Integer trace_mod(Integer t, Integer n, Integer m) {
    if (n == 1) {
        return t;
    }
    Integer a = 1, b = t;
    for (auto bit = std::bit_floor(n); bit != 0; bit >>= 1) {
        Integer even = (2 * multiply_mod(a, a, m) + m - 1) % m;
        Integer odd = (2 * multiply_mod(a, b, m) + m - t) % m;
        if ((n & bit) == 0) {
            a = even;
            b = odd;
        } else {
            a = odd;
            b = (2 * multiply_mod(b, b, m) + m - 1) % m;
        }
    }
    return a;
}

struct RemainderTree {
    std::vector<mpz_class> product;
    std::size_t leaves;
    explicit RemainderTree(const std::vector<unsigned>& primes)
        : leaves(1) {
        while (leaves < primes.size()) {
            leaves *= 2;
        }
        product.resize(2 * leaves, 1);
        for (std::size_t i = 0; i < primes.size(); ++i) {
            product[leaves + i] = Integer(primes[i]) * primes[i];
        }
        for (std::size_t i = leaves - 1; i != 0; --i) {
            product[i] = product[2 * i] * product[2 * i + 1];
        }
    }

    void fill(std::size_t node, const mpz_class& value, std::vector<Integer>& residues) const {
        if (node >= leaves) {
            if (node - leaves < residues.size()) {
                residues[node - leaves] = mpz_fdiv_ui(value.get_mpz_t(), product[node].get_ui());
            }
            return;
        }
        const mpz_class remainder = value % product[node];
        fill(2 * node, remainder, residues);
        fill(2 * node + 1, remainder, residues);
    }
};

bool in_cubic_chain(Integer p) {
    return p % 8 == 3 || p % 8 == 5;
}

void run(Integer radicand, unsigned cutoff, bool list) {
    auto factors = squarefree_factors(radicand);
    auto coefficients = continued_fraction(radicand);
    auto matrix = convergents(coefficients, 0, coefficients.size());
    const auto& t = matrix[0];
    const auto& u = matrix[2];
    if (t*t - radicand*u*u != 1 || mpz_fdiv_ui(t.get_mpz_t(), 4) != 0 || mpz_even_p(u.get_mpz_t())) {
        throw std::runtime_error("unexpected Pell coordinates");
    }
    auto spf = least_prime_factors(cutoff + 1);
    std::vector<bool> forced(cutoff + 2);
    mpz_class period = 1;
    unsigned chain_roots = 0;
    for (auto p : factors) {
        if (mpz_fdiv_ui(u.get_mpz_t(), p) != 0) {
            period *= p;
            chain_roots += in_cubic_chain(p);
            if (p <= cutoff + 1) {
                forced[p] = true;
            }
        }
    }
    const auto initial = period;
    std::vector<unsigned> primes;
    std::vector<unsigned> odd_exponents;
    for (unsigned p = 3; p <= cutoff; p += 2) {
        if (spf[p] != p || radicand % p == 0) {
            continue;
        }
        bool split = power_mod(radicand % p, (p - 1) / 2, p) == 1;
        unsigned group_size = split ? p - 1 : p + 1;
        if (group_size % 4 == 0) {
            primes.push_back(p);
            while (group_size % 2 == 0) {
                group_size /= 2;
            }
            odd_exponents.push_back(group_size);
        }
    }
    const auto tree = RemainderTree(primes);
    std::vector<Integer> trace_residues(primes.size());
    tree.fill(1, t, trace_residues);
    unsigned count = 0, last = 0;
    for (std::size_t i = 0; i < primes.size(); ++i) {
        unsigned p = primes[i];
        unsigned remaining = odd_exponents[i], exponent = 1;
        while (remaining > 1) {
            unsigned q = spf[remaining];
            if (forced[q]) {
                exponent *= q;
            }
            do {
                remaining /= q;
            } while (remaining % q == 0);
        }
        Integer residue = trace_mod(trace_residues[i], exponent, Integer(p) * p);
        if (residue % p == 0 && residue != 0 && !forced[p]) {
            forced[p] = true;
            period *= p;
            ++count;
            last = p;
            unsigned rank = exponent;
            for (unsigned rest = exponent; rest > 1; rest /= spf[rest]) {
                unsigned q = spf[rest];
                if (trace_mod(trace_residues[i] % p, rank / q, p) == 0) {
                    rank /= q;
                }
            }
            bool can_have_predecessor = rank > 1 && spf[rank] == rank
                && in_cubic_chain(rank) && mpz_fdiv_ui(t.get_mpz_t(), rank) != 0;
            if (in_cubic_chain(p) && rank > 1 && !can_have_predecessor) {
                ++chain_roots;
            }
            if (list) {
                std::cout << p << ' ' << rank << ' '
                          << trace_mod(trace_residues[i], rank, Integer(p) * p) << '\n';
            }
        }
    }
    auto decimal = period.get_str();
    std::cout << radicand << ' ' << coefficients.size() << ' ' << initial << ' '
              << count << ' ' << last << ' ' << decimal.size() << ' '
              << chain_roots << '\n';
}

int main(int argc, char** argv) {
    if (argc < 3 || argc > 4) {
        std::cerr << "usage: pell-sweep RADICAND PRIME_CUTOFF [list]\n";
        return 2;
    }
    auto radicand = std::stoull(argv[1]);
    auto cutoff = std::stoul(argv[2]);
    if (radicand < 7 || radicand % 8 != 7 || radicand > 100000000000000ULL ||
        cutoff < 3 || cutoff > 100000000) {
        return 2;
    }
    run(radicand, cutoff, argc == 4);
}

