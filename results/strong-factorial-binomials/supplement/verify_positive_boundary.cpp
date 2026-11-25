// Exact certificate for the finite positive-support boundary ranges.
#include <algorithm>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <vector>

namespace {

constexpr int first_terminal = 25;
constexpr int modulus = 1000000007;

struct Support {
    int terminal, exponent, gap;
};

std::vector<Support> support_range(int group) {
    std::vector<Support> result;
    auto add = [&](int exponent, int gap, int cutoff) {
        if (gap * gap >= 21 * (exponent + gap)) return;
        for (int terminal = first_terminal; terminal < cutoff; ++terminal) {
            result.push_back({terminal, exponent, gap});
        }
    };
    if (group == 0) {
        for (int exponent = 2; exponent <= 8; ++exponent) {
            for (int gap = 3 * exponent; gap <= 24; ++gap) {
                add(exponent, gap, std::max(4000, 128 * (exponent + gap)));
            }
        }
    } else if (group == 1) {
        for (int gap = 10; gap <= 19; ++gap) add(1, gap, 5000);
    } else {
        add(1, 6, 10000);
    }
    return result;
}

int floor_sqrt(int value) {
    int root = 0;
    while (std::int64_t(root + 1) * (root + 1) <= value) ++root;
    return root;
}

class PrimeFactors {
public:
    explicit PrimeFactors(int limit) : least_factor(limit + 1) {
        for (int p = 2; p <= limit; ++p) {
            if (least_factor[p]) continue;
            least_factor[p] = p;
            if (std::int64_t(p) * p > limit) continue;
            for (int j = p * p; j <= limit; j += p) {
                if (!least_factor[j]) least_factor[j] = p;
            }
        }
    }

    bool is_prime(int number) const {
        return number >= 2 && least_factor[number] == number;
    }

    bool block_has_factor_above(int start, int width, int bound) const {
        for (int offset = 1; offset <= width; ++offset) {
            int number = start + offset;
            while (number > bound) {
                const int p = least_factor[number];
                if (p > bound) return true;
                number /= p;
            }
        }
        return false;
    }

private:
    std::vector<int> least_factor;
};

struct Edge {
    int rise, run, width;
};

std::vector<Edge> sparse_polygon(int order, int exponent, int gap, int p) {
    std::vector<std::pair<int, int>> hull;
    auto add = [&](int index, int ordinate) {
        if (!hull.empty() && hull.back().first == index) return;
        while (hull.size() >= 2) {
            const auto [i, y] = hull[hull.size() - 2];
            const auto [j, z] = hull.back();
            if (std::int64_t(z - y) * (index - j)
                < std::int64_t(ordinate - z) * (j - i)) break;
            hull.pop_back();
        }
        hull.emplace_back(index, ordinate);
    };

    // Since p > aN, p > N and p^2 > (a+d)N, the ordinates
    // are floor((a*order + d*k)/p). Keep the right end of each plateau.
    add(0, 0);
    const int jumps = (exponent + gap) * order / p;
    for (int j = 1; j <= jumps; ++j) {
        add((j * p - exponent * order - 1) / gap, j - 1);
    }
    add(order, jumps);

    std::vector<Edge> result;
    for (std::size_t j = 1; j < hull.size(); ++j) {
        const int width = hull[j].first - hull[j - 1].first;
        const int height = hull[j].second - hull[j - 1].second;
        const int common = std::gcd(width, height);
        result.push_back({height / common, width / common, width});
    }
    return result;
}

void retain_possible_degrees(std::vector<int>& degrees, const Support& support, int p) {
    const auto [terminal, exponent, gap] = support;
    const auto lower = sparse_polygon(terminal - 1, exponent, gap, p);
    const auto upper = sparse_polygon(terminal, exponent, gap, p);
    const int limit = degrees.back();
    std::vector<char> reachable(limit + 1), next(limit + 1);
    reachable[0] = true;
    for (const Edge& left : lower) {
        for (const Edge& right : upper) {
            if (left.rise != right.rise || left.run != right.run) continue;
            const int step = left.run;
            const int choices = std::min(left.width, right.width) / step + 1;
            if (choices == 1 || step > limit) continue;
            // For each residue class, OR over the allowed translates
            // by a sliding count on the old vector.
            for (int residue = 0; residue < step; ++residue) {
                int count = 0;
                for (int k = residue; k <= limit; k += step) {
                    count += reachable[k];
                    if (k - choices * step >= 0) count -= reachable[k - choices * step];
                    next[k] = count > 0;
                }
            }
            reachable.swap(next);
        }
    }
    std::erase_if(degrees, [&](int h) { return !reachable[h]; });
}

void retain_prime_support_degrees(std::vector<int>& degrees,
                                  const Support& support,
                                  const PrimeFactors& factors) {
    const auto [terminal, exponent, gap] = support;
    const int x = exponent * terminal;
    const int maximum = (exponent + gap) * terminal;
    const int shape_bound = exponent + 2 * gap - 1;
    const int bound = std::max(shape_bound, floor_sqrt(maximum));
    std::erase_if(degrees, [&](int degree) {
        const int width = gap * degree;
        // The uniform starting-point extension assumes d < 3a.
        // Only the direct starting-point inequality is needed here.
        const bool starting_point = x + 1 > width + 13
                                     && 262 * (x + 1) > 279 * width;
        if (starting_point && 2 * width >= shape_bound
            && std::int64_t(2 * width) * (2 * width) >= maximum) return true;
        const int prime_bound = std::min(gap * terminal,
            starting_point ? bound : std::max(2 * width, bound));
        return factors.block_has_factor_above(x, width, prime_bound)
            || factors.block_has_factor_above(x - exponent, width, prime_bound);
    });
}

int product_mod(int a, int b) {
    return std::int64_t(a) * b % modulus;
}

int inverse_mod(int value) {
    int result = 1;
    for (int k = modulus - 2; k > 0; k /= 2) {
        if (k % 2 != 0) result = product_mod(result, value);
        value = product_mod(value, value);
    }
    return result;
}

std::vector<int> reciprocal_mod(int order, int exponent, int gap,
                                const std::vector<int>& inverses) {
    std::vector<int> result(order + 1);
    int coefficient = 1;
    int argument = exponent * order;
    result[order] = 1;
    for (int k = 1; k <= order; ++k) {
        coefficient = product_mod(product_mod(coefficient, order - k + 1), inverses[k]);
        for (int j = 0; j < gap; ++j) coefficient = product_mod(coefficient, ++argument);
        result[order - k] = coefficient;
    }
    return result;
}

bool coprime_mod(std::vector<int> dividend, std::vector<int> divisor) {
    while (divisor.size() > 1) {
        const int inverse_lead = inverse_mod(divisor.back());
        const int degree = divisor.size() - 1;
        for (int i = int(dividend.size()) - 1; i >= degree; --i) {
            if (!dividend[i]) continue;
            const int scale = product_mod(dividend[i], inverse_lead);
            const int shift = i - degree;
            for (int j = 0; j < degree; ++j) {
                int value = dividend[shift + j] - product_mod(scale, divisor[j]);
                if (value < 0) value += modulus;
                dividend[shift + j] = value;
            }
            dividend[i] = 0;
        }
        while (!dividend.empty() && dividend.back() == 0) dividend.pop_back();
        if (dividend.empty()) return false;
        dividend.swap(divisor);
    }
    return divisor.size() == 1;
}

bool check_fixtures(const std::vector<int>& inverses) {
    if (reciprocal_mod(2, 1, 2, inverses) != std::vector<int>{360, 24, 1}) return false;
    if (coprime_mod({modulus - 1, 0, 1}, {modulus - 1, 1})) return false;
    if (!coprime_mod({1, 0, 1}, {modulus - 1, 1})) return false;
    // The two hulls have common slopes 2/3 and 1 of minimum widths
    // 3 and 1, respectively. The only positive even size they allow is 4.
    std::vector<int> degrees;
    for (int h = 2; h <= 24; h += 2) degrees.push_back(h);
    retain_possible_degrees(degrees, {25, 1, 19}, 29);
    return degrees == std::vector<int>{4};
}

} // namespace

int main() {
    for (int p = 2; std::int64_t(p) * p <= modulus; ++p) {
        if (modulus % p == 0) return 1;
    }

    std::vector<std::vector<Support>> groups;
    int maximum_argument = 0;
    int maximum_terminal = 0;
    for (int group = 0; group < 3; ++group) {
        groups.push_back(support_range(group));
        for (const auto& [terminal, exponent, gap] : groups.back()) {
            const int b = exponent + gap;
            // These checks justify using p>aN in the sparse hull formula.
            if (gap >= terminal || b > exponent * exponent * terminal
                || exponent >= gap * (terminal - 2)) return 1;
            maximum_argument = std::max(maximum_argument, b * terminal);
            maximum_terminal = std::max(maximum_terminal, terminal);
        }
    }
    if (maximum_argument >= modulus || maximum_terminal < 2) return 1;
    const PrimeFactors factors(maximum_argument);
    std::vector<int> inverses(maximum_terminal + 1);
    for (int k = 1; k <= maximum_terminal; ++k) inverses[k] = inverse_mod(k);
    if (!check_fixtures(inverses)) return 1;

    std::cout << "Prime factor table upper endpoint: " << maximum_argument << '\n';
    const char* labels[] = {"relative gap", "lower one, gaps 10 to 19", "lower one, gap six"};
    bool passed = true;
    for (std::size_t group = 0; group < groups.size(); ++group) {
        const auto& supports = groups[group];
        std::int64_t block = 0, polygon = 0, modular = 0, unresolved = 0;
        #pragma omp parallel for schedule(dynamic) reduction(+:block,polygon,modular,unresolved)
        for (std::size_t index = 0; index < supports.size(); ++index) {
            const auto& support = supports[index];
            const auto [terminal, exponent, gap] = support;
            const int x = exponent * terminal;
            const int maximum = (exponent + gap) * terminal;
            std::vector<int> degrees;
            for (int h = 2; h <= terminal - 1; h += 2) degrees.push_back(h);
            bool used_block_test = false;
            for (int p = (x + 1) | 1; p <= maximum && !degrees.empty(); p += 2) {
                if (!factors.is_prime(p)) continue;
                retain_possible_degrees(degrees, support, p);
                if (!degrees.empty() && !used_block_test) {
                    retain_prime_support_degrees(degrees, support, factors);
                    used_block_test = true;
                    if (degrees.empty()) ++block;
                } else if (degrees.empty()) {
                    ++polygon;
                }
            }
            if (degrees.empty()) continue;
            if (!used_block_test) {
                retain_prime_support_degrees(degrees, support, factors);
                if (degrees.empty()) {
                    ++block;
                    continue;
                }
            }
            if (coprime_mod(reciprocal_mod(terminal, exponent, gap, inverses),
                            reciprocal_mod(terminal - 1, exponent, gap, inverses))) {
                ++modular;
            } else {
                #pragma omp critical
                std::cerr << "Unresolved (N,a,d): " << terminal << ' ' << exponent << ' ' << gap << '\n';
                ++unresolved;
            }
        }
        std::cout << labels[group] << ": " << supports.size() << " supports, "
                  << block << " prime-block exclusions, " << polygon << " Newton exclusions, "
                  << modular << " modular exclusions, " << unresolved << " unresolved\n";
        passed = passed && unresolved == 0
                 && block + polygon + modular == std::int64_t(supports.size());
    }
    return passed ? 0 : 1;
}

