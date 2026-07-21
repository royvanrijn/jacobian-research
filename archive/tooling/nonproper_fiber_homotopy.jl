# Numerical adversarial benchmark for HomotopyContinuation.jl.
# Run after: import Pkg; Pkg.add("HomotopyContinuation")
using HomotopyContinuation
using LinearAlgebra
using Printf

@var x y z s
u = 1 + x*y
F = System([
    u^3*z + y^2*u*(4 + 3*x*y) - (-1//4 + s),
    y + 3*x*u^2*z + 3*x*y^2*(4 + 3*x*y),
    2*x - 3*x^2*y - x^3*z,
]; variables=[x, y, z], parameters=[s])

function main()
    # Exact start fiber at s=0.  Keeping paths separate preserves their labels.
    paths = [
        ComplexF64[-1, 3/2, 13/2],
        ComplexF64[ 1, -3/2, 13/2],
        ComplexF64[ 0, 0, -1/4],
    ]
    active = trues(length(paths))
    current_s = 0.0
    deltas = 10.0 .^ (-(1:10))

    println("delta,path,status,observed_norm,predicted_norm,relative_error,projective_distance,jacobian_condition")
    for delta in deltas
    target_s = 0.25 - delta
    for i in eachindex(paths)
        active[i] || continue
        result = solve(F, [paths[i]];
            start_parameters=[current_s], target_parameters=[target_s])
        sols = solutions(result; only_nonsingular=false)
        if isempty(sols)
            active[i] = false
            @printf("%.1e,%d,LOST,NaN,NaN,NaN,NaN,NaN\n", delta, i)
            continue
        end
        paths[i] = ComplexF64.(sols[1])

        if i <= 2
            sign_t = i == 1 ? 1.0 : -1.0
            t = sign_t * sqrt(delta)
            predicted = ComplexF64[-1/(2t), 3t, 26delta]
            relerr = norm(paths[i] - predicted) / norm(predicted)
            # Normalize [x:y:z:1] in the affine chart X=-1.  The exact
            # projective representative there is [-1:6t^2:52t^3:2t].
            observed_projective = ComplexF64[-1, -paths[i][2]/paths[i][1],
                                              -paths[i][3]/paths[i][1],
                                              -1/paths[i][1]]
            predicted_projective = ComplexF64[-1, 6t^2, 52t^3, 2t]
            pdist = norm(observed_projective - predicted_projective)
        else
            predicted = ComplexF64[0, 0, -delta]
            relerr = norm(paths[i] - predicted) / max(norm(predicted), eps())
            pdist = 0.0
        end
        xv, yv, zv = paths[i]
        J = ComplexF64[
            3*xv^2*yv^3*zv + 6*xv*yv^4 + 6*xv*yv^2*zv + 7*yv^3 + 3*yv*zv  3*xv^3*yv^2*zv + 12*xv^2*yv^3 + 6*xv^2*yv*zv + 21*xv*yv^2 + 3*xv*zv + 8*yv  (1+xv*yv)^3;
            9*xv^2*yv^2*zv + 18*xv*yv^3 + 12*xv*yv*zv + 12*yv^2 + 3*zv  6*xv^3*yv*zv + 27*xv^2*yv^2 + 6*xv^2*zv + 24*xv*yv + 1  3*xv*(1+xv*yv)^2;
            2 - 6*xv*yv - 3*xv^2*zv  -3*xv^2  -xv^3
        ]
        @printf("%.1e,%d,OK,%.16e,%.16e,%.3e,%.3e,%.3e\n",
                delta, i, norm(paths[i]), norm(predicted), relerr, pdist, cond(J))
    end
        current_s = target_s
    end

    println("Expected: paths 1-2 have norm ~ 1/(2sqrt(delta)); path 3 has norm=delta.")
    println("Expected projective limit for paths 1-2: [-1:0:0:0].")
    println("Expected condition-number growth on paths 1-2: order delta^-3.")
    println("Final affine tracker state: ", active,
            " (LOST is an observed finite-precision divergence threshold, not a singular endpoint).")
    active[3] || error("the bounded path was unexpectedly lost")
end

main()
