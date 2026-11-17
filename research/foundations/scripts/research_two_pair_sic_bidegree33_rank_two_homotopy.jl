#!/usr/bin/env julia

"""
Numerical algebraic-geometry scout for the generic rank-two cubic chart.

The chart is

    F(x,y) = (1+x)B(y) + x^2(lambda+x)D(y),

with B(0)=1 and the coefficient of y^3 in D eliminated by the first
moment.  The seven equations mu_2,...,mu_8 are square in the remaining
seven variables.  This script constructs them without an intermediate
CAS export and reports their mixed volume.  With `--solve`, it also tracks
the polyhedral homotopy and prints numerical solutions together with
later-moment residuals.

Everything produced here is numerical reconnaissance.  It is not an
exact component classification or a characteristic-zero certificate.
"""

using HomotopyContinuation
using LinearAlgebra
using Printf

@var lambda a1 a2 a3 b0 b1 b2
const VARIABLES = [lambda, a1, a2, a3, b0, b1, b2]


function convolve2(left, right)
    answer = Matrix{Any}(
        undef,
        size(left, 1) + size(right, 1) - 1,
        size(left, 2) + size(right, 2) - 1,
    )
    fill!(answer, 0)
    for left_j in axes(left, 2), left_i in axes(left, 1)
        iszero(left[left_i, left_j]) && continue
        for right_j in axes(right, 2), right_i in axes(right, 1)
            iszero(right[right_i, right_j]) && continue
            answer[left_i + right_i - 1, left_j + right_j - 1] +=
                left[left_i, left_j] * right[right_i, right_j]
        end
    end
    return answer
end


function coefficient_matrix()
    d3 = -1 - a1 / 3 - lambda * b2 / 3
    B = [1, a1, a2, a3]
    D = [b0, b1, b2, d3]
    coefficient = Matrix{Any}(undef, 4, 4)
    fill!(coefficient, 0)
    for y_degree in 0:3
        coefficient[1, y_degree + 1] += B[y_degree + 1]
        coefficient[2, y_degree + 1] += B[y_degree + 1]
        coefficient[3, y_degree + 1] += lambda * D[y_degree + 1]
        coefficient[4, y_degree + 1] += D[y_degree + 1]
    end
    return coefficient
end


function moments_through(maximum::Int)
    coefficient = coefficient_matrix()
    power = Matrix{Any}(undef, 1, 1)
    power[1, 1] = 1
    moments = Vector{Any}(undef, maximum)
    for order in 1:maximum
        power = convolve2(power, coefficient)
        value = 0
        for diagonal in 0:(3 * order)
            value += (
                factorial(big(diagonal))
                * factorial(big(3 * order - diagonal))
                * power[diagonal + 1, diagonal + 1]
            )
        end
        moments[order] = value
        println("generated mu_", order)
        flush(stdout)
    end
    return moments
end


function main(arguments)
    do_solve = "--solve" in arguments
    maximum = do_solve ? 14 : 8
    moments = moments_through(maximum)
    equations = moments[2:8]
    system = System(equations; variables=VARIABLES)
    volume = mixed_volume(system)
    println("MIXED_VOLUME ", volume)
    do_solve || return

    result = solve(
        system;
        start_system=:polyhedral,
        show_progress=true,
        threading=true,
    )
    candidates = solutions(result; only_nonsingular=false)
    @printf(
        "SOLVE paths=%d finite=%d nonsingular=%d solutions=%d\n",
        npaths(result),
        nfinite(result),
        nnonsingular(result),
        length(candidates),
    )
    for (index, candidate) in enumerate(candidates)
        later = ComplexF64[
            evaluate(moment, VARIABLES => candidate)
            for moment in moments[9:14]
        ]
        @printf(
            "CANDIDATE %d norm=%.9e max_mu9_14=%.9e\n",
            index,
            norm(candidate),
            maximum(abs.(later)),
        )
        println(join(candidate, " "))
        println(join(later, " "))
    end
end


main(ARGS)
