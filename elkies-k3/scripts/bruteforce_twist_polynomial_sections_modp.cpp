#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

namespace {

int mod(std::int64_t value, int prime) {
    value %= prime;
    if (value < 0) value += prime;
    return static_cast<int>(value);
}

int inverse(int value, int prime) {
    for (int candidate = 1; candidate < prime; ++candidate) {
        if (mod(value * candidate, prime) == 1) return candidate;
    }
    throw std::runtime_error("noninvertible field element");
}

int evaluate(const std::vector<int>& coefficients, int value, int prime) {
    int result = 0;
    for (auto it = coefficients.rbegin(); it != coefficients.rend(); ++it) {
        result = mod(static_cast<std::int64_t>(result) * value + *it, prime);
    }
    return result;
}

std::vector<int> convolution(
    const std::vector<int>& left,
    const std::vector<int>& right,
    int prime
) {
    std::vector<int> result(left.size() + right.size() - 1, 0);
    for (std::size_t i = 0; i < left.size(); ++i) {
        for (std::size_t j = 0; j < right.size(); ++j) {
            result[i + j] = mod(
                result[i + j] + static_cast<std::int64_t>(left[i]) * right[j],
                prime
            );
        }
    }
    return result;
}

std::vector<std::vector<int>> invert_matrix(
    std::vector<std::vector<int>> matrix,
    int prime
) {
    const int size = static_cast<int>(matrix.size());
    std::vector<std::vector<int>> augmented(size, std::vector<int>(2 * size, 0));
    for (int row = 0; row < size; ++row) {
        for (int column = 0; column < size; ++column) {
            augmented[row][column] = matrix[row][column];
        }
        augmented[row][size + row] = 1;
    }
    for (int column = 0; column < size; ++column) {
        int pivot = column;
        while (pivot < size && !augmented[pivot][column]) ++pivot;
        if (pivot == size) throw std::runtime_error("singular interpolation matrix");
        std::swap(augmented[pivot], augmented[column]);
        const int scale = inverse(augmented[column][column], prime);
        for (int entry = 0; entry < 2 * size; ++entry) {
            augmented[column][entry] = mod(
                static_cast<std::int64_t>(augmented[column][entry]) * scale,
                prime
            );
        }
        for (int row = 0; row < size; ++row) {
            if (row == column) continue;
            const int multiplier = augmented[row][column];
            for (int entry = 0; entry < 2 * size; ++entry) {
                augmented[row][entry] = mod(
                    augmented[row][entry]
                        - static_cast<std::int64_t>(multiplier)
                            * augmented[column][entry],
                    prime
                );
            }
        }
    }
    std::vector<std::vector<int>> result(size, std::vector<int>(size));
    for (int row = 0; row < size; ++row) {
        for (int column = 0; column < size; ++column) {
            result[row][column] = augmented[row][size + column];
        }
    }
    return result;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "usage: bruteforce_twist_polynomial_sections_modp INPUT\n";
        return 2;
    }
    std::ifstream input(argv[1]);
    if (!input) throw std::runtime_error("cannot open input");

    int prime = 0;
    int chi = 0;
    int block_count = 0;
    input >> prime >> chi >> block_count;
    if (prime < 3 || chi <= 0 || block_count <= 0) {
        throw std::runtime_error("invalid header");
    }
    int a_count = 0;
    input >> a_count;
    std::vector<int> A(a_count);
    for (int& value : A) input >> value;
    int b_count = 0;
    input >> b_count;
    std::vector<int> B(b_count);
    for (int& value : B) input >> value;

    std::map<int, int> representative_y;
    for (int index = 0; index < block_count; ++index) {
        int block_index = 0;
        int leading_x = 0;
        int leading_y = 0;
        input >> block_index >> leading_x >> leading_y;
        representative_y.emplace(leading_x, leading_y);
    }
    if (!input) throw std::runtime_error("truncated input");

    const int x_degree = 2 * chi;
    const int y_degree = 3 * chi;
    if (x_degree != 8 || y_degree != 12) {
        throw std::runtime_error("this bounded enumerator currently requires chi=4");
    }

    std::vector<bool> square(prime, false);
    for (int value = 0; value < prime; ++value) {
        square[mod(static_cast<std::int64_t>(value) * value, prime)] = true;
    }
    std::vector<int> a_values(prime), b_values(prime);
    for (int value = 0; value < prime; ++value) {
        a_values[value] = evaluate(A, value, prime);
        b_values[value] = evaluate(B, value, prime);
    }

    std::uint64_t tested = 0;
    std::uint64_t passed_value_sieve = 0;
    std::uint64_t solutions = 0;
    const bool interpolate_values = prime > x_degree;
    for (const auto& [leading_x, leading_y] : representative_y) {
        auto test_polynomial = [&](const std::vector<int>& X, int first_value) {
            ++tested;
            bool passes = true;
            for (int value = first_value; value < prime && passes; ++value) {
                const int x_value = evaluate(X, value, prime);
                const int rhs_value = mod(
                    static_cast<std::int64_t>(x_value) * x_value % prime * x_value
                        + static_cast<std::int64_t>(a_values[value]) * x_value
                        + b_values[value],
                    prime
                );
                passes = square[rhs_value];
            }
            if (passes) {
                ++passed_value_sieve;
                const std::vector<int> X2 = convolution(X, X, prime);
                std::vector<int> rhs = convolution(X2, X, prime);
                const std::vector<int> AX = convolution(A, X, prime);
                rhs.resize(2 * y_degree + 1, 0);
                for (std::size_t degree = 0; degree < AX.size(); ++degree) {
                    rhs[degree] = mod(rhs[degree] + AX[degree], prime);
                }
                for (std::size_t degree = 0; degree < B.size(); ++degree) {
                    rhs[degree] = mod(rhs[degree] + B[degree], prime);
                }

                std::vector<int> Y(y_degree + 1, 0);
                Y[y_degree] = leading_y;
                const int denominator_inverse = inverse(2 * leading_y, prime);
                for (int degree = 2 * y_degree - 1; degree >= y_degree; --degree) {
                    const int index = degree - y_degree;
                    int known = 0;
                    for (int left = index + 1; left <= y_degree; ++left) {
                        const int right = degree - left;
                        if (right > index && right <= y_degree) {
                            known = mod(
                                known + static_cast<std::int64_t>(Y[left]) * Y[right],
                                prime
                            );
                        }
                    }
                    Y[index] = mod(
                        static_cast<std::int64_t>(rhs[degree] - known)
                            * denominator_inverse,
                        prime
                    );
                }
                const std::vector<int> Y2 = convolution(Y, Y, prime);
                if (Y2 == rhs) {
                    ++solutions;
                    std::cout << "SOLUTION " << leading_x << ' ' << leading_y;
                    for (int value : X) std::cout << ' ' << value;
                    for (int value : Y) std::cout << ' ' << value;
                    std::cout << '\n';
                }
            }
        };

        if (!interpolate_values) {
            std::vector<int> X(x_degree + 1, 0);
            X[x_degree] = leading_x;
            bool finished = false;
            while (!finished) {
                test_polynomial(X, 0);
                int position = 0;
                while (position < x_degree) {
                    if (++X[position] < prime) break;
                    X[position] = 0;
                    ++position;
                }
                finished = position == x_degree;
            }
        } else {
            std::vector<std::vector<int>> vandermonde(
                x_degree, std::vector<int>(x_degree, 1)
            );
            for (int row = 0; row < x_degree; ++row) {
                for (int column = 1; column < x_degree; ++column) {
                    vandermonde[row][column] = mod(
                        static_cast<std::int64_t>(vandermonde[row][column - 1]) * row,
                        prime
                    );
                }
            }
            const auto inverse_vandermonde = invert_matrix(vandermonde, prime);
            std::vector<std::vector<int>> allowed(x_degree);
            for (int value = 0; value < x_degree; ++value) {
                for (int x_value = 0; x_value < prime; ++x_value) {
                    const int rhs_value = mod(
                        static_cast<std::int64_t>(x_value) * x_value % prime * x_value
                            + static_cast<std::int64_t>(a_values[value]) * x_value
                            + b_values[value],
                        prime
                    );
                    if (square[rhs_value]) allowed[value].push_back(x_value);
                }
                if (allowed[value].empty()) {
                    throw std::runtime_error("empty interpolation value set");
                }
            }
            std::vector<std::vector<std::vector<int>>> contribution(x_degree);
            for (int value = 0; value < x_degree; ++value) {
                contribution[value].resize(
                    allowed[value].size(), std::vector<int>(x_degree, 0)
                );
                int leading_term = 1;
                for (int exponent = 0; exponent < x_degree; ++exponent) {
                    leading_term = mod(
                        static_cast<std::int64_t>(leading_term) * value, prime
                    );
                }
                for (std::size_t option = 0; option < allowed[value].size(); ++option) {
                    const int adjusted = mod(
                        allowed[value][option]
                            - static_cast<std::int64_t>(leading_x) * leading_term,
                        prime
                    );
                    for (int coefficient = 0; coefficient < x_degree; ++coefficient) {
                        contribution[value][option][coefficient] = mod(
                            static_cast<std::int64_t>(
                                inverse_vandermonde[coefficient][value]
                            ) * adjusted,
                            prime
                        );
                    }
                }
            }
            std::vector<std::size_t> digits(x_degree, 0);
            std::vector<int> X(x_degree + 1, 0);
            X[x_degree] = leading_x;
            for (int value = 0; value < x_degree; ++value) {
                for (int coefficient = 0; coefficient < x_degree; ++coefficient) {
                    X[coefficient] = mod(
                        X[coefficient] + contribution[value][0][coefficient], prime
                    );
                }
            }
            bool finished = false;
            while (!finished) {
                test_polynomial(X, x_degree);
                int position = 0;
                while (position < x_degree) {
                    const std::size_t old_digit = digits[position];
                    const std::size_t new_digit = old_digit + 1;
                    for (int coefficient = 0; coefficient < x_degree; ++coefficient) {
                        X[coefficient] = mod(
                            X[coefficient]
                                - contribution[position][old_digit][coefficient],
                            prime
                        );
                    }
                    if (new_digit < allowed[position].size()) {
                        digits[position] = new_digit;
                        for (int coefficient = 0; coefficient < x_degree; ++coefficient) {
                            X[coefficient] = mod(
                                X[coefficient]
                                    + contribution[position][new_digit][coefficient],
                                prime
                            );
                        }
                        break;
                    }
                    digits[position] = 0;
                    for (int coefficient = 0; coefficient < x_degree; ++coefficient) {
                        X[coefficient] = mod(
                            X[coefficient] + contribution[position][0][coefficient],
                            prime
                        );
                    }
                    ++position;
                }
                finished = position == x_degree;
            }
        }
        std::cerr << "completed leading_x=" << leading_x << " tested=" << tested
                  << " solutions=" << solutions << '\n';
    }
    std::cout << "SUMMARY " << tested << ' ' << passed_value_sieve << ' '
              << solutions << ' '
              << (interpolate_values ? "value_interpolation" : "coefficient_odometer")
              << '\n';
    return 0;
}
