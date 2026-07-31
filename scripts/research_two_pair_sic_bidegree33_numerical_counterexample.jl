#!/usr/bin/env julia

"""
Numerical counterexample scout on a normalized non-null quadratic chart.

The Clebsch--Gordan coordinates and basis polynomials are the same as in
`explore_two_pair_sic_bidegree33_full_anchor.py`.  One higher-weight
coordinate is fixed to one, the quadratic component is fixed to `2*X*T`,
and complex Newton iteration solves eleven selected normalized moments in
the remaining eleven variables.  Every reported point is only a numerical
lead; the script evaluates additional moments but proves no zero identity.
"""

using LinearAlgebra
using Printf
using Random

const PARAMETER_NAMES = [
    "s0", "s1", "s2", "s3", "s4", "s5", "s6",
    "t0", "t1", "t2", "t3", "t4",
]

const WEIGHTS = [3, 2, 1, 0, -1, -2, -3, 2, 1, 0, -1, -2]

const Q_POLYNOMIALS = [
    [1],
    [-3, 3],
    [3, -9, 3],
    [-1, 9, -9, 1],
    [-3, 9, -3],
    [-3, 3],
    [-1],
    [1, 1],
    [-2, 0, 2],
    [1, -3, -3, 1],
    [2, 0, -2],
    [1, 1],
]

const QUADRATIC_Q = [-1, -1, 1, 1]


function basis_matrix(weight::Int, q_polynomial::Vector{Int})
    result = zeros(ComplexF64, 4, 4)
    x_shift = max(weight, 0)
    y_shift = max(-weight, 0)
    for (offset, coefficient) in enumerate(q_polynomial)
        q_degree = offset - 1
        result[x_shift + q_degree + 1, y_shift + q_degree + 1] += coefficient
    end
    return result
end


const PARAMETER_BASES = [
    basis_matrix(weight, polynomial)
    for (weight, polynomial) in zip(WEIGHTS, Q_POLYNOMIALS)
]

const QUADRATIC_BASE = basis_matrix(0, QUADRATIC_Q)


function convolve2(left::Matrix{ComplexF64}, right::Matrix{ComplexF64})
    answer = zeros(
        ComplexF64,
        size(left, 1) + size(right, 1) - 1,
        size(left, 2) + size(right, 2) - 1,
    )
    for left_j in axes(left, 2), left_i in axes(left, 1)
        left_value = left[left_i, left_j]
        iszero(left_value) && continue
        for right_j in axes(right, 2), right_i in axes(right, 1)
            right_value = right[right_i, right_j]
            iszero(right_value) && continue
            answer[left_i + right_i - 1, left_j + right_j - 1] += (
                left_value * right_value
            )
        end
    end
    return answer
end


function beta_weights(order::Int)
    degree = 3 * order
    weights = zeros(Float64, degree + 1)
    weights[1] = 1 / (degree + 1)
    for diagonal in 0:(degree - 1)
        weights[diagonal + 2] = (
            weights[diagonal + 1]
            * (diagonal + 1)
            / (degree - diagonal)
        )
    end
    return weights
end


function full_parameters(
    variables::Vector{ComplexF64},
    chart_index::Int,
)
    parameters = zeros(ComplexF64, length(PARAMETER_NAMES))
    parameters[chart_index] = 1
    source_index = 1
    for target_index in eachindex(parameters)
        target_index == chart_index && continue
        parameters[target_index] = variables[source_index]
        source_index += 1
    end
    return parameters
end


function coefficient_matrix(
    variables::Vector{ComplexF64},
    chart_index::Int,
)
    parameters = full_parameters(variables, chart_index)
    coefficient = copy(QUADRATIC_BASE)
    for index in eachindex(parameters)
        coefficient .+= parameters[index] .* PARAMETER_BASES[index]
    end
    return coefficient
end


function powers_through(
    coefficient::Matrix{ComplexF64},
    maximum::Int,
)
    powers = Vector{Matrix{ComplexF64}}(undef, maximum + 1)
    powers[1] = ones(ComplexF64, 1, 1)
    for order in 1:maximum
        powers[order + 1] = convolve2(powers[order], coefficient)
    end
    return powers
end


function moment_from_power(power::Matrix{ComplexF64}, order::Int)
    weights = beta_weights(order)
    return sum(
        weights[diagonal + 1] * power[diagonal + 1, diagonal + 1]
        for diagonal in 0:(3 * order)
    )
end


function moment_derivative(
    previous_power::Matrix{ComplexF64},
    basis::Matrix{ComplexF64},
    order::Int,
)
    weights = beta_weights(order)
    answer = 0.0 + 0.0im
    for basis_j in axes(basis, 2), basis_i in axes(basis, 1)
        basis_value = basis[basis_i, basis_j]
        iszero(basis_value) && continue
        x_degree = basis_i - 1
        y_degree = basis_j - 1
        minimum_diagonal = max(x_degree, y_degree)
        maximum_diagonal = min(
            x_degree + size(previous_power, 1) - 1,
            y_degree + size(previous_power, 2) - 1,
        )
        for diagonal in minimum_diagonal:maximum_diagonal
            answer += (
                order
                * weights[diagonal + 1]
                * basis_value
                * previous_power[
                    diagonal - x_degree + 1,
                    diagonal - y_degree + 1,
                ]
            )
        end
    end
    return answer
end


function residual_and_jacobian(
    variables::Vector{ComplexF64},
    chart_index::Int,
    orders::Vector{Int};
    with_jacobian::Bool=true,
)
    coefficient = coefficient_matrix(variables, chart_index)
    powers = powers_through(coefficient, maximum(orders))
    residual = ComplexF64[
        moment_from_power(powers[order + 1], order)
        for order in orders
    ]
    if !with_jacobian
        return residual, zeros(ComplexF64, 0, 0)
    end
    free_indices = [
        index
        for index in eachindex(PARAMETER_NAMES)
        if index != chart_index
    ]
    jacobian = zeros(ComplexF64, length(orders), length(free_indices))
    for (row, order) in enumerate(orders)
        previous_power = powers[order]
        for (column, parameter_index) in enumerate(free_indices)
            jacobian[row, column] = moment_derivative(
                previous_power,
                PARAMETER_BASES[parameter_index],
                order,
            )
        end
    end
    return residual, jacobian
end


function newton_solve(
    initial::Vector{ComplexF64},
    chart_index::Int,
    orders::Vector{Int};
    tolerance::Float64=1e-11,
    maximum_iterations::Int=80,
)
    variables = copy(initial)
    best_variables = copy(variables)
    best_norm = Inf
    for iteration in 1:maximum_iterations
        residual, jacobian = residual_and_jacobian(
            variables,
            chart_index,
            orders,
        )
        residual_norm = norm(residual)
        if residual_norm < best_norm
            best_norm = residual_norm
            best_variables .= variables
        end
        if residual_norm < tolerance
            return variables, residual_norm, iteration, true
        end
        if !isfinite(residual_norm) || rank(jacobian; rtol=1e-12) < length(variables)
            break
        end
        step = -(jacobian \ residual)
        step_norm = norm(step)
        if step_norm > 20 * (1 + norm(variables))
            step .*= 20 * (1 + norm(variables)) / step_norm
        end
        accepted = false
        scale = 1.0
        for _ in 1:14
            candidate = variables .+ scale .* step
            candidate_residual, _ = residual_and_jacobian(
                candidate,
                chart_index,
                orders;
                with_jacobian=false,
            )
            if norm(candidate_residual) < residual_norm
                variables .= candidate
                accepted = true
                break
            end
            scale /= 2
        end
        accepted || break
    end
    return best_variables, best_norm, maximum_iterations, false
end


function least_squares_refine(
    initial::Vector{ComplexF64},
    chart_index::Int,
    orders::Vector{Int};
    tolerance::Float64=1e-12,
    maximum_iterations::Int=160,
)
    variables = copy(initial)
    best_variables = copy(variables)
    best_norm = Inf
    damping = 1e-8
    for iteration in 1:maximum_iterations
        residual, jacobian = residual_and_jacobian(
            variables,
            chart_index,
            orders,
        )
        residual_norm = norm(residual)
        if residual_norm < best_norm
            best_norm = residual_norm
            best_variables .= variables
        end
        if residual_norm < tolerance
            return variables, residual_norm, iteration, true
        end
        normal_matrix = jacobian' * jacobian
        diagonal_scale = max(opnorm(normal_matrix, Inf), 1e-20)
        step = -(
            normal_matrix
            + damping * diagonal_scale * I
        ) \ (jacobian' * residual)
        if !all(isfinite, step)
            damping *= 10
            continue
        end
        accepted = false
        scale = 1.0
        for _ in 1:16
            candidate = variables .+ scale .* step
            candidate_residual, _ = residual_and_jacobian(
                candidate,
                chart_index,
                orders;
                with_jacobian=false,
            )
            if norm(candidate_residual) < residual_norm
                variables .= candidate
                damping = max(damping / 3, 1e-16)
                accepted = true
                break
            end
            scale /= 2
        end
        if !accepted
            damping *= 10
            damping > 1e16 && break
        end
    end
    return best_variables, best_norm, maximum_iterations, false
end


function parse_orders(specification::String)
    if occursin(":", specification)
        endpoints = parse.(Int, split(specification, ":"))
        length(endpoints) == 2 || error("orders must be A:B or comma-separated")
        return collect(endpoints[1]:endpoints[2])
    end
    return parse.(Int, split(specification, ","))
end


function parse_arguments(arguments::Vector{String})
    options = Dict{String,String}(
        "chart" => "s0",
        "orders" => "2:12",
        "check-through" => "30",
        "trials" => "100",
        "seed" => "33001",
        "scale" => "0.75",
        "refine-through" => "0",
    )
    for argument in arguments
        startswith(argument, "--") || error("unknown argument $argument")
        pieces = split(argument[3:end], "="; limit=2)
        length(pieces) == 2 || error("arguments must use --name=value")
        options[pieces[1]] = pieces[2]
    end
    chart_index = findfirst(==(options["chart"]), PARAMETER_NAMES)
    isnothing(chart_index) && error("unknown chart $(options["chart"])")
    return (
        chart_index=chart_index,
        orders=parse_orders(options["orders"]),
        check_through=parse(Int, options["check-through"]),
        trials=parse(Int, options["trials"]),
        seed=parse(Int, options["seed"]),
        scale=parse(Float64, options["scale"]),
        refine_through=parse(Int, options["refine-through"]),
    )
end


function print_candidate(
    variables::Vector{ComplexF64},
    chart_index::Int,
    residual_norm::Float64,
    iteration::Int,
    check_through::Int,
)
    parameters = full_parameters(variables, chart_index)
    check_orders = collect(1:check_through)
    moments, _ = residual_and_jacobian(
        variables,
        chart_index,
        check_orders;
        with_jacobian=false,
    )
    @printf(
        "CANDIDATE residual=%.6e iteration=%d max_check=%.6e\n",
        residual_norm,
        iteration,
        maximum(abs.(moments)),
    )
    for (name, value) in zip(PARAMETER_NAMES, parameters)
        @printf("  %s = %.17g %+.17gim\n", name, real(value), imag(value))
    end
    for (order, value) in zip(check_orders, moments)
        @printf(
            "  nu_%d = %.9e %+.9eim abs=%.9e\n",
            order,
            real(value),
            imag(value),
            abs(value),
        )
    end
end


function main(arguments::Vector{String})
    options = parse_arguments(arguments)
    length(options.orders) == length(PARAMETER_NAMES) - 1 || error(
        "the Newton system must contain exactly eleven moments",
    )
    Random.seed!(options.seed)
    best_norm = Inf
    best_variables = zeros(ComplexF64, length(PARAMETER_NAMES) - 1)
    successes = 0
    best_refined_norm = Inf
    best_refined_variables = copy(best_variables)
    for trial in 1:options.trials
        initial = options.scale .* (
            randn(ComplexF64, length(PARAMETER_NAMES) - 1)
        )
        variables, residual_norm, iteration, converged = newton_solve(
            initial,
            options.chart_index,
            options.orders,
        )
        if residual_norm < best_norm
            best_norm = residual_norm
            best_variables .= variables
            @printf(
                "BEST trial=%d residual=%.6e converged=%s\n",
                trial,
                residual_norm,
                converged,
            )
            flush(stdout)
        end
        if converged
            successes += 1
            print_candidate(
                variables,
                options.chart_index,
                residual_norm,
                iteration,
                options.check_through,
            )
            flush(stdout)
            if options.refine_through > 0
                refine_orders = collect(2:options.refine_through)
                (
                    refined_variables,
                    refined_norm,
                    refined_iteration,
                    refined_converged,
                ) = least_squares_refine(
                    variables,
                    options.chart_index,
                    refine_orders,
                )
                if refined_norm < best_refined_norm
                    best_refined_norm = refined_norm
                    best_refined_variables .= refined_variables
                    @printf(
                        "BEST_REFINED trial=%d residual=%.6e converged=%s\n",
                        trial,
                        refined_norm,
                        refined_converged,
                    )
                    print_candidate(
                        refined_variables,
                        options.chart_index,
                        refined_norm,
                        refined_iteration,
                        options.check_through,
                    )
                    flush(stdout)
                end
            end
        end
    end
    @printf(
        "SUMMARY trials=%d successes=%d best_residual=%.6e\n",
        options.trials,
        successes,
        best_norm,
    )
    if options.refine_through > 0
        @printf(
            "REFINE_SUMMARY through=%d best_residual=%.6e\n",
            options.refine_through,
            best_refined_norm,
        )
    end
    if successes == 0
        print_candidate(
            best_variables,
            options.chart_index,
            best_norm,
            0,
            options.check_through,
        )
    end
end


main(ARGS)
