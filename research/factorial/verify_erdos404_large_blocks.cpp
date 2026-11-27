
#include <algorithm>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <string_view>
#include <vector>

namespace {

using Word = std::uint64_t;
using Wide = unsigned __int128;
using Bits = std::vector<Word>;

Word power_of_five(unsigned exponent) {
    Word result = 1;
    for (unsigned i = 0; i < exponent; ++i) {
        result *= 5;
    }
    return result;
}

unsigned factorial_valuation(Word index) {
    unsigned valuation = 0;
    while ((index /= 5) != 0) {
        valuation += index;
    }
    return valuation;
}

Word multiply_mod(Word a, Word b, Word modulus) {
    return static_cast<Wide>(a) * b % modulus;
}

Word inverse_unit(Word unit, Word modulus) {
    // Euler's theorem for the positive power-of-five modulus.
    Word exponent = 4 * (modulus / 5) - 1;
    Word inverse = 1;
    while (exponent != 0) {
        if (exponent % 2 != 0) {
            inverse = multiply_mod(inverse, unit, modulus);
        }
        unit = multiply_mod(unit, unit, modulus);
        exponent /= 2;
    }
    return inverse;
}

Word shifted_start(Word offset, Word shift, Word modulus) {
    if (offset < shift) {
        offset += modulus;
    }
    return offset - shift;
}

Word last_mask(Word modulus) {
    const unsigned tail = modulus % 64;
    return tail == 0 ? ~Word{0} : (Word{1} << tail) - 1;
}

Word cyclic_word(const Bits& bits, Word start, Word modulus) {
    const Word index = start / 64;
    const unsigned offset = start % 64;
    Word result = bits[index] >> offset;
    if (offset != 0 && index + 1 < bits.size()) {
        result |= bits[index + 1] << (64 - offset);
    }
    const Word tail = modulus - start;
    if (tail < 64) {
        result = (result & ((Word{1} << tail) - 1)) | (bits[0] << tail);
    }
    return result;
}

Bits suffix_subsums(Word split, Word boundary, Word modulus) {
    Bits reachable((modulus + 63) / 64);
    Bits next(reachable.size());
    Word weight = 1;
    reachable[0] = 1;
    for (Word index = split; index < boundary; ++index) {
        if (index != split) {
            weight = multiply_mod(weight, index, modulus);
        }
        if (weight == 0) {
            continue;
        }
        #pragma omp parallel for schedule(static) if(reachable.size() > 100000)
        for (Word word = 0; word < reachable.size(); ++word) {
            const Word start = shifted_start(64 * word, weight, modulus);
            next[word] = reachable[word] | cyclic_word(reachable, start, modulus);
        }
        next.back() &= last_mask(modulus);
        reachable.swap(next);
    }
    return reachable;
}

bool covers(Word boundary, Word width, Word suffix_width) {
    const Word start = boundary - width;
    const Word split = boundary - suffix_width;
    const Word low = power_of_five(factorial_valuation(split) - factorial_valuation(start));
    const Word high = power_of_five(factorial_valuation(boundary) - factorial_valuation(split));
    const Word modulus = low * high;

    Word unit = 1;
    for (Word index = start + 1; index <= split; ++index) {
        Word factor = index;
        while (factor % 5 == 0) {
            factor /= 5;
        }
        unit = multiply_mod(unit, factor, high);
    }
    const Word inverse = inverse_unit(unit, high);

    // Prefix values use start!, whereas suffix values use split!.
    std::vector<Word> prefix{0};
    prefix.reserve(Word{1} << (split - start));
    Word weight = 1;
    for (Word index = start; index < split; ++index) {
        if (index != start) {
            weight = multiply_mod(weight, index, modulus);
        }
        const auto length = prefix.size();
        for (Word i = 0; i < length; ++i) {
            prefix.push_back((prefix[i] + weight) % modulus);
        }
    }
    const Word total_low = prefix.back() % low;
    std::vector<std::vector<Word>> shifts(low);
    for (Word value : prefix) {
        shifts[value % low].push_back(multiply_mod(value / low, inverse, high));
    }
    for (auto& group : shifts) {
        if (group.empty()) {
            return false;
        }
        std::sort(group.begin(), group.end());
        group.erase(std::unique(group.begin(), group.end()), group.end());
        if (group.size() < 4) {
            continue;
        }
        // Visit well separated shifts first. The coprime step visits each once.
        const Word length = group.size();
        Word step = length / 2;
        while (std::gcd(step, length) != 1) {
            ++step;
        }
        std::vector<Word> reordered;
        reordered.reserve(length);
        for (Word i = 0, index = 0; i < length; ++i, index = (index + step) % length) {
            reordered.push_back(group[index]);
        }
        group.swap(reordered);
    }

    Bits missing = suffix_subsums(split, boundary, high);
    for (auto& word : missing) {
        word = ~word;
    }
    missing.back() &= last_mask(high);
    if (std::all_of(missing.begin(), missing.end(), [](Word word) { return word == 0; })) {
        return true;
    }

    int failed = 0;
    #pragma omp parallel for schedule(dynamic, 1) reduction(|:failed)
    for (Word residue = 0; residue < low; ++residue) {
        // Taking the complementary subset interchanges these two low cosets.
        if (residue > shifted_start(total_low, residue, low)) {
            continue;
        }
        for (Word word = 0; word < missing.size(); ++word) {
            Word candidates = word + 1 == missing.size() ? last_mask(high) : ~Word{0};
            for (Word shift : shifts[residue]) {
                candidates &= cyclic_word(missing, shifted_start(64 * word, shift, high), high);
                if (candidates == 0) {
                    break;
                }
            }
            if (candidates != 0) {
                failed = 1;
                break;
            }
        }
    }
    return failed == 0;
}

bool check_jumps(unsigned jump, Word limit, Word width, Word suffix_width) {
    const Word step = power_of_five(jump);
    unsigned checked = 0;
    for (Word boundary = step; boundary < limit; boundary += step) {
        if (boundary / step % 5 == 0) {
            continue;
        }
        if (!covers(boundary, width, suffix_width)) {
            std::cerr << "Incomplete coverage at " << boundary << " with width " << width << '\n';
            return false;
        }
        ++checked;
    }
    std::cout << "Jump " << jump << ", width " << width
              << ", verified boundaries " << checked << std::endl;
    return true;
}

} // namespace

int main(int argc, char* argv[]) {
    if (argc > 2 || (argc == 2 && std::string_view(argv[1]) != "9")) {
        std::cerr << "Run without arguments for the certificate below 5^8, or pass 9\n";
        return 2;
    }
    if (covers(15625, 35, 25) || covers(87500, 35, 35)) {
        std::cerr << "A negative control unexpectedly covered every residue\n";
        return 1;
    }
    const Word limit = power_of_five(argc == 2 ? 9 : 8);
    if (!check_jumps(5, limit, 40, 35) ||
        !check_jumps(6, limit, 40, 25) ||
        !check_jumps(7, limit, 45, 25) ||
        (argc == 2 && !check_jumps(8, limit, 50, 35))) {
        return 1;
    }
}
