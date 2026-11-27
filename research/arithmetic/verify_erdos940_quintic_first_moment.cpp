#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {

constexpr int prime_cutoff = 20'000;
constexpr std::uint64_t scale = 1'000'000'000'000;
constexpr std::uint64_t root_scale = 1'000'000;
constexpr int max_orbits = 7;
using Wide = unsigned __int128;
using Counts = std::array<std::uint64_t, max_orbits>;
using LocalFactors = std::array<std::uint64_t, 4>;

Wide power(std::uint64_t value, int exponent) {
    Wide result = 1;
    for (int i = 0; i < exponent; ++i) {
        result *= value;
    }
    return result;
}

int fifth_power_modulo(int x, int modulus) {
    const std::int64_t square = static_cast<std::int64_t>(x) * x % modulus;
    return static_cast<int>((square * square % modulus) * x % modulus);
}

// Round downward throughout. The binary search lower-bounds p^(-1/5).
std::uint64_t prime_weight_lower(int p) {
    std::uint64_t low = 0;
    std::uint64_t high = root_scale;
    const Wide limit = power(root_scale, 5);
    while (high - low > 1) {
        const std::uint64_t midpoint = (low + high) / 2;
        if (power(midpoint, 5) * p <= limit) {
            low = midpoint;
        } else {
            high = midpoint;
        }
    }
    std::uint64_t result = 0;
    for (int j = 1; j < 5; ++j) {
        result += static_cast<std::uint64_t>(
            scale * power(low, j) / (p * power(root_scale, j)));
    }
    return result;
}

std::uint64_t scaled_product(std::uint64_t left, std::uint64_t right) {
    const Wide result = static_cast<Wide>(left) * right / scale;
    if (result > std::numeric_limits<std::uint64_t>::max()) {
        throw std::overflow_error("Fixed-point multiplication overflow");
    }
    return static_cast<std::uint64_t>(result);
}

struct PowerOrbits {
    int modulus;
    bool prime_power;
    std::vector<int> orbits;
    std::vector<int> representatives;
    std::vector<int> unit_orbits;
    int ways[max_orbits][max_orbits][max_orbits]{};

    explicit PowerOrbits(int p)
        : modulus(p == 5 ? 25 : p), prime_power(p == 5), orbits(modulus, -1) {
        std::vector<int> subgroup;
        std::vector<bool> seen(modulus);
        const int residue_prime = p;
        for (int x = 1; x < modulus; ++x) {
            if (x % residue_prime == 0) {
                continue;
            }
            const int value = fifth_power_modulo(x, modulus);
            if (!seen[value]) {
                seen[value] = true;
                subgroup.push_back(value);
            }
        }
        for (int x = 0; x < modulus; ++x) {
            if (orbits[x] >= 0) {
                continue;
            }
            const int index = static_cast<int>(representatives.size());
            representatives.push_back(x);
            if (x % residue_prime != 0) {
                unit_orbits.push_back(index);
            }
            for (int h : subgroup) {
                orbits[static_cast<std::int64_t>(x) * h % modulus] = index;
            }
        }
        if (representatives.size() > max_orbits || unit_orbits.size() != 5) {
            throw std::runtime_error("Unexpected fifth-power orbit count");
        }
        for (int z = 0; z < static_cast<int>(representatives.size()); ++z) {
            for (int x = 0; x < modulus; ++x) {
                const int y = (representatives[z] - x + modulus) % modulus;
                ++ways[z][orbits[x]][orbits[y]];
            }
        }
    }

    Counts add_unit_coefficient(Counts previous, int coefficient_orbit) const {
        Counts result{};
        const int count = static_cast<int>(representatives.size());
        const int zero_weight = prime_power ? 5 : 1;
        for (int z = 0; z < count; ++z) {
            for (int x = 0; x < count; ++x) {
                const int completions =
                    zero_weight * ways[z][x][0] + 5 * ways[z][x][coefficient_orbit];
                result[z] += previous[x] * completions;
            }
        }
        return result;
    }

    void minimize(int active, int start, Counts counts, Counts& minima) const {
        if (active >= 4) {
            const std::uint64_t singular =
                prime_power ? static_cast<std::uint64_t>(power(5, active)) : 1;
            const auto end = counts.begin() + representatives.size();
            minima[active] = std::min({minima[active], counts[0] - singular,
                *std::min_element(counts.begin() + 1, end)});
        }
        if (active < 5) {
            for (int next = start; next < static_cast<int>(unit_orbits.size()); ++next) {
                minimize(active + 1, next, add_unit_coefficient(counts, unit_orbits[next]), minima);
            }
        }
    }

    LocalFactors bounds() const {
        Counts minima;
        minima.fill(std::numeric_limits<std::uint64_t>::max());
        minimize(0, 0, Counts{1}, minima);
        LocalFactors density{};
        density[0] = static_cast<std::uint64_t>(
            static_cast<Wide>(scale) * minima[5] / power(modulus, 4));
        density[1] = static_cast<std::uint64_t>(
            static_cast<Wide>(scale) * minima[4] / power(modulus, 3));
        if (prime_power && (minima[5] != 203'125 || minima[4] != 5'625)) {
            throw std::runtime_error("Unexpected minimum modulo 25");
        }
        if (modulus == 11 && (minima[5] != 3'125 || minima[4] != 0)) {
            throw std::runtime_error("Unexpected minimum modulo 11");
        }
        return density;
    }
};

LocalFactors local_densities(int p) {
    if (p == 5 || p % 5 == 1) {
        return PowerOrbits(p).bounds();
    }
    LocalFactors density{};
    for (int k = 0; k <= 3; ++k) {
        const Wide denominator = power(p, 4 - k);
        density[k] = static_cast<std::uint64_t>(scale * (denominator - 1) / denominator);
    }
    return density;
}

}  // namespace

int main() {
    std::vector<bool> prime(prime_cutoff + 1, true);
    for (int p = 2; p * p <= prime_cutoff; ++p) {
        if (prime[p]) {
            for (int multiple = p * p; multiple <= prime_cutoff; multiple += p) {
                prime[multiple] = false;
            }
        }
    }

    std::uint64_t product = scale;
    constexpr std::array<int, 4> choices{1, 5, 10, 10};
    for (int p = 2; p <= prime_cutoff; ++p) {
        if (!prime[p]) {
            continue;
        }
        const LocalFactors density = local_densities(p);
        const std::uint64_t weight = prime_weight_lower(p);
        std::uint64_t factor = 0;
        std::uint64_t weight_power = scale;
        for (int k = 0; k <= 3; ++k) {
            factor += choices[k] * scaled_product(weight_power, density[k]);
            weight_power = scaled_product(weight_power, weight);
        }
        product = scaled_product(product, factor);
    }

    std::cout << "Lower numerator for the product through " << prime_cutoff
              << " at scale " << scale << ": " << product << '\n';
    if (product < 3'000 * scale) {
        std::cerr << "The product lower bound did not reach 3000\n";
        return 1;
    }
}

