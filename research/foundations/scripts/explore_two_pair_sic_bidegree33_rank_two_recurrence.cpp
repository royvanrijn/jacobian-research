// Modular recurrence probe for a generic rank-two bidegree-(3,3) point.
//
// This is an exploratory helper, not a recurrence certificate.  It computes
// the factorial-normalized moments
//
//   nu_m = mu_m / (3m+1)!
//
// from a factorization C=U*W and searches for a scalar P-recurrence of
// prescribed order and coefficient degree over a prime field.

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using i64 = std::int64_t;
using Polynomial = std::vector<i64>;

static i64 prime;

static i64 mod(i64 value) {
    value %= prime;
    return value < 0 ? value + prime : value;
}

static i64 power_mod(i64 base, i64 exponent) {
    i64 answer = 1;
    base = mod(base);
    while (exponent) {
        if (exponent & 1) {
            answer = answer * base % prime;
        }
        base = base * base % prime;
        exponent >>= 1;
    }
    return answer;
}

static Polynomial multiply_cubic(
    const Polynomial& input,
    const Polynomial& cubic
) {
    Polynomial output(input.size() + 3, 0);
    for (std::size_t index = 0; index < input.size(); ++index) {
        for (int shift = 0; shift <= 3; ++shift) {
            output[index + shift] = mod(
                output[index + shift] + input[index] * cubic[shift]
            );
        }
    }
    return output;
}

static std::vector<i64> moments(int maximum, int point) {
    // A deterministic exact-rank-two factor point.  The small entries are
    // interpreted modulo prime only after forming C=U*W.
    const i64 u_points[2][4][2] = {
        {
            {1, 0},
            {0, 1},
            {2, 3},
            {5, 7},
        },
        {
            {1, 2},
            {3, 5},
            {7, 11},
            {13, 17},
        },
    };
    const i64 w_points[2][2][4] = {
        {
            {11, 13, 17, 19},
            {23, 29, 31, 37},
        },
        {
            {19, 23, 29, 31},
            {37, 41, 43, 47},
        },
    };
    if (point < 0) {
        throw std::runtime_error("point must be nonnegative");
    }
    Polynomial a1(4), a2(4), p1(4), p2(4);
    for (int index = 0; index < 4; ++index) {
        if (point < 2) {
            a1[index] = mod(u_points[point][index][0]);
            a2[index] = mod(u_points[point][index][1]);
            p1[index] = mod(w_points[point][0][index]);
            p2[index] = mod(w_points[point][1][index]);
        } else if (point < 10000) {
            const i64 parameter = point - 2;
            a1[index] = mod(
                u_points[0][index][0]
                + parameter * u_points[1][index][0]
            );
            a2[index] = mod(
                u_points[0][index][1]
                + parameter * u_points[1][index][1]
            );
            p1[index] = mod(
                w_points[0][0][index]
                + parameter * w_points[1][0][index]
            );
            p2[index] = mod(
                w_points[0][1][index]
                + parameter * w_points[1][1][index]
            );
        } else {
            const i64 scale = point - 9999;
            a1[index] = mod(scale * u_points[0][index][0]);
            a2[index] = mod(scale * u_points[0][index][1]);
            p1[index] = mod(w_points[0][0][index]);
            p2[index] = mod(w_points[0][1][index]);
        }
    }

    std::vector<i64> factorial(3 * maximum + 2, 1);
    for (int index = 1; index < static_cast<int>(factorial.size()); ++index) {
        factorial[index] = factorial[index - 1] * index % prime;
    }

    std::vector<Polynomial> dual{{1}};
    std::vector<Polynomial> coordinate{{1}};
    std::vector<i64> sequence;
    sequence.reserve(maximum + 1);

    for (int m = 0; m <= maximum; ++m) {
        std::vector<i64> binomial(m + 1, 1);
        for (int k = 1; k < m; ++k) {
            binomial[k] = mod(
                binomial[k - 1] * (m - k + 1) % prime
                * power_mod(k, prime - 2)
            );
        }

        i64 mu = 0;
        for (int k = 0; k <= m; ++k) {
            i64 pairing = 0;
            for (int degree = 0; degree <= 3 * m; ++degree) {
                pairing = mod(
                    pairing
                    + dual[k][degree]
                    * coordinate[k][degree] % prime
                    * factorial[degree] % prime
                    * factorial[3 * m - degree]
                );
            }
            mu = mod(mu + binomial[k] * pairing);
        }
        sequence.push_back(
            mu * power_mod(factorial[3 * m + 1], prime - 2) % prime
        );

        if (m == maximum) {
            break;
        }
        std::vector<Polynomial> next_dual(m + 2);
        std::vector<Polynomial> next_coordinate(m + 2);
        for (int k = 0; k <= m; ++k) {
            next_dual[k] = multiply_cubic(dual[k], a2);
            next_coordinate[k] = multiply_cubic(coordinate[k], p2);
        }
        next_dual[m + 1] = multiply_cubic(dual[m], a1);
        next_coordinate[m + 1] = multiply_cubic(coordinate[m], p1);
        dual = std::move(next_dual);
        coordinate = std::move(next_coordinate);
    }
    return sequence;
}

static std::vector<i64> null_vector(std::vector<std::vector<i64>> matrix) {
    const int rows = static_cast<int>(matrix.size());
    const int columns = static_cast<int>(matrix[0].size());
    std::vector<int> pivots;
    int row = 0;
    for (int column = 0; column < columns && row < rows; ++column) {
        int pivot = row;
        while (pivot < rows && matrix[pivot][column] == 0) {
            ++pivot;
        }
        if (pivot == rows) {
            continue;
        }
        std::swap(matrix[pivot], matrix[row]);
        i64 inverse = power_mod(matrix[row][column], prime - 2);
        for (int j = column; j < columns; ++j) {
            matrix[row][j] = matrix[row][j] * inverse % prime;
        }
        for (int i = row + 1; i < rows; ++i) {
            if (matrix[i][column] == 0) {
                continue;
            }
            i64 scale = matrix[i][column];
            for (int j = column; j < columns; ++j) {
                matrix[i][j] = mod(
                    matrix[i][j] - scale * matrix[row][j]
                );
            }
        }
        pivots.push_back(column);
        ++row;
    }
    if (row >= columns) {
        return {};
    }
    std::vector<bool> is_pivot(columns, false);
    for (int pivot : pivots) {
        is_pivot[pivot] = true;
    }
    int free_column = columns - 1;
    while (free_column >= 0 && is_pivot[free_column]) {
        --free_column;
    }
    if (free_column < 0) {
        return {};
    }
    std::vector<i64> answer(columns, 0);
    answer[free_column] = 1;
    for (int i = static_cast<int>(pivots.size()) - 1; i >= 0; --i) {
        i64 value = 0;
        for (int j = pivots[i] + 1; j < columns; ++j) {
            value = mod(value + matrix[i][j] * answer[j]);
        }
        answer[pivots[i]] = mod(-value);
    }
    return answer;
}

static bool verify_recurrence(
    const std::vector<i64>& sequence,
    const std::vector<i64>& coefficients,
    int order,
    int degree
) {
    for (int m = 0; m + order < static_cast<int>(sequence.size()); ++m) {
        i64 total = 0;
        for (int shift = 0; shift <= order; ++shift) {
            i64 m_power = 1;
            for (int exponent = 0; exponent <= degree; ++exponent) {
                total = mod(
                    total
                    + coefficients[shift * (degree + 1) + exponent]
                    * m_power % prime
                    * sequence[m + shift]
                );
                m_power = m_power * m % prime;
            }
        }
        if (total != 0) {
            return false;
        }
    }
    return true;
}

int main(int argc, char** argv) {
    if (argc >= 2 && std::string(argv[1]) == "--moments") {
        if (argc != 4 && argc != 5) {
            std::cerr
                << "usage: recurrence-probe --moments PRIME MAXIMUM "
                << "[POINT]\n";
            return 2;
        }
        prime = std::atoll(argv[2]);
        const int maximum = std::atoi(argv[3]);
        const int point = argc == 5 ? std::atoi(argv[4]) : 0;
        if (3 * maximum + 1 >= prime) {
            throw std::runtime_error(
                "prime must exceed 3*maximum+1"
            );
        }
        const std::vector<i64> sequence = moments(maximum, point);
        std::cout
            << "MOMENTS maximum=" << maximum
            << " point=" << point << "\n";
        for (int index = 0; index <= maximum; ++index) {
            std::cout << index << " " << sequence[index] << "\n";
        }
        return 0;
    }
    if (argc != 5 && argc != 6) {
        std::cerr
            << "usage: recurrence-probe PRIME MAXIMUM ORDER DEGREE [POINT]\n";
        return 2;
    }
    prime = std::atoll(argv[1]);
    const int maximum = std::atoi(argv[2]);
    const int order = std::atoi(argv[3]);
    const int degree = std::atoi(argv[4]);
    const int point = argc == 6 ? std::atoi(argv[5]) : 0;
    if (3 * maximum + 1 >= prime) {
        throw std::runtime_error("prime must exceed 3*maximum+1");
    }
    const int columns = (order + 1) * (degree + 1);
    if (maximum - order + 1 < columns + 8) {
        throw std::runtime_error("not enough moments for fit plus holdout");
    }

    std::vector<i64> sequence = moments(maximum, point);
    std::vector<std::vector<i64>> matrix(
        columns - 1,
        std::vector<i64>(columns, 0)
    );
    for (int m = 0; m < columns - 1; ++m) {
        for (int shift = 0; shift <= order; ++shift) {
            i64 m_power = 1;
            for (int exponent = 0; exponent <= degree; ++exponent) {
                matrix[m][shift * (degree + 1) + exponent] =
                    sequence[m + shift] * m_power % prime;
                m_power = m_power * m % prime;
            }
        }
    }
    std::vector<i64> coefficients = null_vector(std::move(matrix));
    if (
        coefficients.empty()
        || !verify_recurrence(sequence, coefficients, order, degree)
    ) {
        std::cout << "NONE\n";
        return 1;
    }
    std::cout << "FOUND order=" << order << " degree=" << degree << "\n";
    for (int shift = 0; shift <= order; ++shift) {
        std::cout << shift;
        for (int exponent = 0; exponent <= degree; ++exponent) {
            std::cout << " "
                      << coefficients[shift * (degree + 1) + exponent];
        }
        std::cout << "\n";
    }
    return 0;
}
