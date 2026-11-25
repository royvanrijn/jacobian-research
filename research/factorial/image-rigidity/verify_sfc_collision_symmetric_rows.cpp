#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

// Search ordinary orders 2, 3 and 4 on three-column supports with
// U/V swap symmetry and variable copies of the split and merged rows.
// Run without repo dependencies: g++ -std=c++20 -O3 this_file.cpp -o lattice.
using U = std::uint32_t;
constexpr int count_moments = 19;
using Moment = std::array<U, count_moments>;
using Row = std::array<int, 3>;
using Params = std::array<int, 7>;

constexpr std::array<Row, count_moments> moment_exponents = {{
    {2,0,0}, {1,1,0}, {1,0,1}, {0,0,2},
    {3,0,0}, {2,1,0}, {2,0,1}, {1,1,1}, {1,0,2}, {0,0,3},
    {4,0,0}, {3,1,0}, {2,2,0}, {3,0,1}, {2,1,1}, {2,0,2},
    {1,1,2}, {1,0,3}, {0,0,4}
}};

constexpr bool is_prime(U n) {
    if (n < 2) return false;
    for (U d = 2; d*d <= n; ++d) {
        if (n % d == 0) return false;
    }
    return true;
}

template <U prime>
struct Sieve {
    static_assert(is_prime(prime));
    using Polynomial = std::array<U, 5>;
    std::vector<U> factorial;
    explicit Sieve(int bound) : factorial(4 * bound + 1, 1) {
        for (U i = 1; i < factorial.size(); ++i) {
            factorial[i] = multiply(factorial[i - 1], i);
        }
    }

    static U multiply(U left, U right) {
        return (static_cast<std::uint64_t>(left) * right) % prime;
    }

    static U residue(std::int64_t x) {
        x %= prime;
        return x < 0 ? x + prime : x;
    }

    static U subtract(U left, U right) {
        return left >= right ? left - right : left + prime - right;
    }

    Moment row_moments(const Row& row) const {
        Moment result;
        for (int i = 0; i < count_moments; ++i) {
            const Row& exponents = moment_exponents[i];
            result[i] = factorial[row[0]*exponents[0] + row[1]*exponents[1] + row[2]*exponents[2]];
        }
        return result;
    }

    static Moment moment_product(const Moment& left, const Moment& right) {
        Moment result;
        for (int i = 0; i < count_moments; ++i) {
            result[i] = multiply(left[i], right[i]);
        }
        return result;
    }

    Moment pair_moments(int a, int b, int c) const {
        return moment_product(row_moments({a,b,c}), row_moments({b,a,c}));
    }

    static std::vector<Moment> moment_powers(const Moment& base, int bound) {
        std::vector<Moment> powers(bound + 1);
        powers[0].fill(1);
        for (int i = 1; i <= bound; ++i) {
            powers[i] = moment_product(powers[i-1], base);
        }
        return powers;
    }

    static bool projectively_coprime(Polynomial left, Polynomial right, int left_degree, int right_degree) {
        // Keep the homogeneous degrees when testing the point at infinity.
        if (left[left_degree] == 0 && right[right_degree] == 0) {
            return false;
        }
        while (left_degree >= 0 && left[left_degree] == 0) --left_degree;
        while (right_degree >= 0 && right[right_degree] == 0) --right_degree;
        while (right_degree >= 0) {
            if (right_degree == 0) return true;
            while (left_degree >= right_degree) {
                const U left_high = left[left_degree], right_high = right[right_degree];
                const int shift = left_degree - right_degree;
                for (int i = 0; i < left_degree; ++i) {
                    const U correction = i < shift ? 0 : multiply(left_high, right[i-shift]);
                    left[i] = subtract(multiply(right_high, left[i]), correction);
                }
                left[left_degree] = 0;
                while (left_degree >= 0 && left[left_degree] == 0) --left_degree;
            }
            left.swap(right);
            std::swap(left_degree, right_degree);
        }
        return left_degree == 0;
    }

    static bool excludes(const Moment& m) {
        const auto [m200,m110,m101,m002,
                    m300,m210,m201,m111,m102,m003,
                    m400,m310,m220,m301,m211,m202,m112,m103,m004] = [&]() {
            std::array<std::int64_t, count_moments> result;
            for (int i = 0; i < count_moments; ++i) result[i] = static_cast<std::int64_t>(m[i]);
            return result;
        }();
        const U v2 = residue(2*m200 - 2*m110);
        if (v2 == 0) return false;
        const U v4 = residue(2*m400 - 8*m310 + 6*m220);
        const std::array<U, 3> m2{residue(2*m200 + 2*m110), residue(4*m101), residue(m002)};
        const std::array<U, 4> m3{residue(2*m300 + 6*m210), residue(6*m201 + 6*m111), residue(6*m102), residue(m003)};
        const std::array<U, 2> cross3{residue(2*m300 - 2*m210), residue(2*m201 - 2*m111)};
        const std::array<U, 5> m4{
            residue(2*m400 + 8*m310 + 6*m220), residue(8*m301 + 24*m211),
            residue(12*m202 + 12*m112), residue(8*m103), residue(m004)
        };
        const std::array<U, 3> cross4{
            residue(2*m400 - 2*m220), residue(4*m301 - 4*m211),
            residue(2*m202 - 2*m112)
        };
        // These polynomials eliminate z² from the second moment in the
        // coefficient chart (r+z,r-z,t), and retain homogeneous degrees 3,4.
        Polynomial cubic{};
        for (int j = 0; j <= 3; ++j) cubic[j] = multiply(v2, m3[j]);
        for (int j = 0; j <= 2; ++j) {
            for (int k = 0; k <= 1; ++k) {
                cubic[j+k] = subtract(cubic[j+k], multiply(multiply(3, m2[j]), cross3[k]));
            }
        }
        Polynomial quartic;
        const U v2_squared = multiply(v2, v2);
        for (int j = 0; j <= 4; ++j) quartic[j] = multiply(v2_squared, m4[j]);
        for (int j = 0; j <= 2; ++j) {
            const U negative = multiply(multiply(6, v2), m2[j]);
            const U positive = multiply(v4, m2[j]);
            for (int k = 0; k <= 2; ++k) {
                quartic[j+k] = residue(static_cast<std::int64_t>(quartic[j+k])
                    - multiply(negative, cross4[k]) + multiply(positive, m2[k]));
            }
        }
        return projectively_coprime(cubic, quartic, 3, 4);
    }

    bool excludes_parameters(const Params& p) const {
        const auto [u,v,a,b,c,d,e] = p;
        const auto split = moment_powers(pair_moments(5,2,3), u);
        const auto merged = moment_powers(row_moments({7,7,6}), v);
        return excludes(moment_product(moment_product(split[u], merged[v]),
                                        moment_product(pair_moments(a,b,c), row_moments({d,d,e}))));
    }
};


template <U prime>
bool self_test() {
    using SmallSieve = Sieve<prime>;
    using Poly = typename SmallSieve::Polynomial;
    // The first pair crosses at t=1; the second pair crosses at infinity.
    if (SmallSieve::projectively_coprime(Poly{prime-1,1}, Poly{prime-1,0,1}, 1, 2)
        || SmallSieve::projectively_coprime(Poly{1,1}, Poly{1,0,1}, 3, 4)
        || !SmallSieve::projectively_coprime(Poly{1,1}, Poly{1,0,1}, 1, 2)) {
        return false;
    }

    // Independent real R,B with weights (1,4,1) at (-1,0,1) give
    // E[(R+iB)^n]=0 for n=2,3,4. Represent this as a three-column
    // swap-symmetric law with U=4+R+B, V=4+R-B, W=8.
    Moment collision{};
    constexpr std::array<int, 3> weights{1,4,1};
    for (int r = -1; r <= 1; ++r) {
        for (int b = -1; b <= 1; ++b) {
            const Row values{4+r+b, 4+r-b, 8};
            for (int i = 0; i < count_moments; ++i) {
                U weight = weights[r+1]*weights[b+1];
                for (int j = 0; j < 3; ++j) {
                    for (int k = 0; k < moment_exponents[i][j]; ++k) {
                        weight = SmallSieve::multiply(weight, values[j]);
                    }
                }
                collision[i] = SmallSieve::residue(collision[i] + weight);
            }
        }
    }
    return !SmallSieve::excludes(collision);
}

template <U prime>
void filter_residue_classes(const Params& p, int common_period,
                            std::vector<std::array<int,2>>& unresolved) {
    if (unresolved.empty()) return;
    const auto [u,v,a,b,c,d,e] = p;
    const Sieve<prime> sieve(std::max({7, a, b, c, d, e}));
    const auto split = sieve.moment_powers(sieve.pair_moments(5,2,3), prime-2);
    const auto merged = sieve.moment_powers(sieve.row_moments({7,7,6}), prime-2);
    const auto base = sieve.moment_product(sieve.pair_moments(a,b,c), sieve.row_moments({d,d,e}));
    // Nonzero row factorials are units. It suffices to check every lift
    // of the common residue pair to one full pair of periods mod prime-1.
    const auto needs_another_prime = [&](const auto& congruences) {
        for (U i = congruences[0]; i < prime-1; i += common_period) {
            for (U j = congruences[1]; j < prime-1; j += common_period) {
                const auto moment = sieve.moment_product(base, sieve.moment_product(split[i], merged[j]));
                if (!sieve.excludes(moment)) return true;
            }
        }
        return false;
    };
    std::erase_if(unresolved, [&](const auto& c) { return !needs_another_prime(c); });
}

bool infinite_multiplicities(const Params& p, bool verbose) {
    constexpr int period = 36;
    std::vector<std::array<int,2>> unresolved;
    for (int u = 0; u < period; ++u) {
        for (int v = 0; v < period; ++v) {
            unresolved.push_back({u,v});
        }
    }
    const auto apply = [&]<U prime>() {
        static_assert((prime - 1) % period == 0);
        filter_residue_classes<prime>(p, period, unresolved);
        if (verbose) std::cout << "after " << prime << ": " << unresolved.size() << "\n";
    };
    apply.template operator()<37>();
    apply.template operator()<73>();
    apply.template operator()<109>();
    apply.template operator()<181>();
    apply.template operator()<397>();
    apply.template operator()<433>();
    apply.template operator()<577>();
    apply.template operator()<613>();
    apply.template operator()<757>();
    apply.template operator()<829>();
    apply.template operator()<937>();
    apply.template operator()<1009>();
    apply.template operator()<1117>();
    apply.template operator()<1153>();
    apply.template operator()<1297>();
    apply.template operator()<1549>();
    if (verbose && !unresolved.empty()) {
        std::cout << "unresolved classes\n";
        for (const auto& [u,v] : unresolved) std::cout << u << " " << v << "\n";
    }
    return unresolved.empty();
}

int main(int argc, char** argv) {
    constexpr U prime = 65521;

    if (!self_test<prime>()) {
        std::cerr << "self-test failed\n";
        return 1;
    }
    if (argc == 2 && std::string(argv[1]) == "--infinite") {
        std::cout << "base pair (5,4,4),(4,5,4), self row (5,5,4)\n";
        return infinite_multiplicities({0,0,5,4,4,5,4}, true) ? 0 : 1;
    }
    if (argc == 2 && std::string(argv[1]) == "--infinite-box") {
        int tested = 0;
        int closed = 0;
        for (int a = 3; a <= 9; ++a)
        for (int b = 1; b < a; ++b)
        for (int c = 1; c <= 8; ++c)
        for (int d = 1; d <= 9; ++d)
        for (int e = 1; e <= 9; ++e) {
            ++tested;
            if (infinite_multiplicities({0,0,a,b,c,d,e}, false)) {
                ++closed;
            } else {
                std::cout << "unresolved base " << a << " " << b << " " << c
                          << " " << d << " " << e << "\n";
            }
        }
        std::cout << "infinite families tested " << tested << "\n";
        std::cout << "infinite families excluded " << closed << "\n";
        return 0;
    }
    if (argc != 4) {
        std::cerr << "usage: lattice MAX_COPIES MAX_PAIR_ENTRY MAX_SELF_ENTRY"
                     " or lattice --infinite or lattice --infinite-box\n";
        return 2;
    }
    int maximum_copies = std::stoi(argv[1]);
    int maximum_pair = std::stoi(argv[2]);
    int maximum_self = std::stoi(argv[3]);
    int bound = std::max({7, maximum_pair, maximum_self});
    if (maximum_pair < 3 || maximum_self < 1 || maximum_copies < 0 ||
        bound > (65449 - 1) / 4) return 2;
    Sieve<prime> first(bound);
    Sieve<65519> second(bound);
    Sieve<65497> third(bound);
    Sieve<65479> fourth(bound);
    Sieve<65449> fifth(bound);
    const std::vector<Moment> split = first.moment_powers(first.pair_moments(5,2,3), maximum_copies);
    const std::vector<Moment> merged = first.moment_powers(first.row_moments({7,7,6}), maximum_copies);
    std::vector<Moment> products;
    for (int u = 0; u <= maximum_copies; ++u) {
        for (int v = 0; v <= maximum_copies; ++v) {
            products.push_back(first.moment_product(split[u], merged[v]));
        }
    }
    std::uint64_t tested = 0;
    std::array<std::uint64_t, 5> excluded{};
    for (int a = 3; a <= maximum_pair; ++a)
    for (int b = 1; b < a; ++b)
    for (int c = 1; c <= maximum_pair - 1; ++c) {
        const Moment pair = first.pair_moments(a,b,c);
        for (int d = 1; d <= maximum_self; ++d)
        for (int e = 1; e <= maximum_self; ++e) {
            const Moment base = first.moment_product(pair, first.row_moments({d,d,e}));
            size_t index = 0;
            for (int u = 0; u <= maximum_copies; ++u)
            for (int v = 0; v <= maximum_copies; ++v) {
                const Moment moment = first.moment_product(base, products[index++]);
                ++tested;
                const Params p{u,v,a,b,c,d,e};
                if (first.excludes(moment)) ++excluded[0];
                else if (second.excludes_parameters(p)) ++excluded[1];
                else if (third.excludes_parameters(p)) ++excluded[2];
                else if (fourth.excludes_parameters(p)) ++excluded[3];
                else if (fifth.excludes_parameters(p)) ++excluded[4];
                else {
                    std::cout << "survivor";
                    for (int coordinate : p) std::cout << " " << coordinate;
                    std::cout << "\n";
                }
            }
        }
    }
    std::cout << "tested " << tested << "\n";
    std::cout << "excluded";
    for (std::uint64_t n : excluded) std::cout << " " << n;
    std::cout << "\n";
}
