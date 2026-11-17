#!/usr/bin/env julia

# Compute the closed-cycle Picard--Fuchs operator of the fixed rank-two fiber
# with Brochet--Chyzak--Lairez reduction.  The open interval boundary must be
# audited separately before this can certify the moment sequence.

using MultivariateCreativeTelescoping

const Q = """
216-648*t+648*t^2-216*t^3
+91*u-5*u*t-263*u*t^2+177*u*t^3
+23*u^2+44*u^2*t+145*u^2*t^2-212*u^2*t^3
+11*u^3-4*u^3*t+102*u^3*t^2+245*u^3*t^3
+13*u^4*t+5*u^4*t^2+131*u^4*t^3
+17*u^5*t^2+20*u^5*t^3
+19*u^6*t^3
"""

# In the original beta coordinates, CT_u(1/(1-z*Q/u^3)) is the
# u-contour integral of this rational differential form, up to 2*pi*i.
const ORIGINAL_INTEGRAND =
    "u^2/(u^3-z*(" * replace(Q, '\n' => "") * "))"

# The birational beta compression
#
#   x=u,  y=u*t/(1-t)
#
# gives P=Phi(x,y)/(x+y)^3 and (du/u)dt=dx*dy/(x+y)^2.  Therefore the
# same generating form is (x+y)dxdy/((x+y)^3-z*Phi).  Its denominator has
# total degree six instead of nine, which is materially smaller for the
# package's projective homogenization.
const PHI = """
19*x^3*y^3+17*x^3*y^2+13*x^3*y+11*x^3
+37*x^2*y^3+31*x^2*y^2+29*x^2*y+23*x^2
+149*x*y^3+127*x*y^2+113*x*y+91*x
+354*y^3+302*y^2+268*y+216
"""
const COMPACT_INTEGRAND =
    "(x+y)/((x+y)^3-z*(" * replace(PHI, '\n' => "") * "))"

method = isempty(ARGS) ? "crt" : ARGS[1]
rho = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 1
coordinates = length(ARGS) >= 3 ? ARGS[3] : "compact"
output_path = length(ARGS) >= 4 ? ARGS[4] : nothing
integrand = if coordinates == "compact"
    COMPACT_INTEGRAND
elseif coordinates == "original"
    ORIGINAL_INTEGRAND
else
    error("coordinates must be compact or original")
end
debug_param = pf_debug_param(debug=Val(true))

result, algebra = if method == "direct"
    picard_fuchs(
        integrand;
        parameters=["z"],
        rho=rho,
        use_trivial_syzygies=true,
        debug_param=debug_param,
    )
else
    picard_fuchs_crt(
        integrand;
        parameters=["z"],
        rho=rho,
        use_trivial_syzygies=true,
        tracer=true,
        crt_parallel=false,
        debug_param=debug_param,
    )
end

output = output_path === nothing ? stdout : open(output_path, "w")
redirect_stdout(output) do
    println("PICARD_FUCHS_BEGIN")
    for operator in result
        prettyprint(operator, algebra)
    end
    println("PICARD_FUCHS_END")
end
output_path === nothing || close(output)

println("PASS coordinates=" * coordinates)
println("PASS operator_count=" * string(length(result)))
for (index, operator) in enumerate(result)
    println(
        "PASS operator_" * string(index) *
        "_ore_terms=" * string(length(operator)),
    )
end
output_path === nothing || println("PASS wrote " * output_path)
