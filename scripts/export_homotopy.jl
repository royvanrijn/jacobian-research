# HomotopyContinuation.jl adapter for the normalized 3D chart residue.
# Run in Julia after: import Pkg; Pkg.add("HomotopyContinuation")
using HomotopyContinuation
@var a20 a11 a30 b10 b01 b21
F = System([
    2*a20*b01 - a11*b10,
    3*a30*b01 - 2*a11*b21,
    a11*b01 - 1//2,
    b01 - 1,
    b10 - 4,
    b21 + 3
])
result = solve(F)
println("finite solutions: ", length(solutions(result)))
for solution in solutions(result)
    println(solution)
end
