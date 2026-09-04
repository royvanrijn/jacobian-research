#!/usr/bin/env julia

# Build the complete global norm-square S-squareclass envelope for ICARM 356.
#
# The expensive class and unit computations are both requested with GRH=false.
# The subsequent p-Selmer call reuses those proved contexts in the same Hecke
# session.  Output is a line-oriented interchange format: a separate local
# condition worker can consume the explicit power-basis representatives using
# only an nfinit, without importing Hecke objects or requiring a PARI BNF.

using Hecke
using Random

const PROTOCOL = "ELKIES356HECKEGLOBAL2-v1"
const RANDOM_SEED = 20260904

function stage(name, status; fields...)
  suffix = join(("|$(key)=$(value)" for (key, value) in fields), "")
  println("$(PROTOCOL)|stage=$(name)|status=$(status)$(suffix)")
  flush(stdout)
end

Random.seed!(RANDOM_SEED)

Qx, x = polynomial_ring(QQ, :x)
A = ZZ(24391876744717707263532695900840552395172973498186560300)
B = ZZ(46943906433780620456844832699051340439698711588743845207309557656274241785479710000)
f = x^3 - x^2 - A*x - B

field_discriminant = -ZZ(752277724884429603556499627496577219114608773920079457931321069663132744982183297505569516372870314046056607683593379327741881715)
ramified_primes = map(ZZ, [
  5,
  29,
  751,
  28960331,
  1204882855601765528877267647500895974865482613,
  197980272243427555346397293722916980361535459279712115031762027678304939,
])
selmer_rational_primes = map(ZZ, [
  2,
  3,
  5,
  13,
  23,
  29,
  37,
  41,
  139,
  751,
  28960331,
  1204882855601765528877267647500895974865482613,
  197980272243427555346397293722916980361535459279712115031762027678304939,
])

@assert all(is_prime, ramified_primes)
@assert all(is_prime, selmer_rational_primes)
@assert prod(ramified_primes) == abs(field_discriminant)

stage("input", "PASS";
  hecke_version=Base.pkgversion(Hecke),
  random_seed=RANDOM_SEED,
  curve=356,
  grh=false,
)

K, theta = number_field(f, :theta)
equation = equation_order(K)
equation_discriminant = discriminant(equation)
index_square = divexact(equation_discriminant, field_discriminant)
is_index_square, defining_order_index = is_square_with_sqrt(index_square)
@assert is_index_square

stage("maximal_order", "start")
O = maximal_order(
  equation;
  discriminant=field_discriminant,
  ramified_primes=ramified_primes,
)
@assert discriminant(O) == field_discriminant
stage("maximal_order", "PASS";
  discriminant=field_discriminant,
  defining_order_index=defining_order_index,
)

stage("class_group", "start"; grh=false)
C, mC = class_group(O; GRH=false)
class_invariants = elementary_divisors(C)
class_2rank = count(is_even, class_invariants)
stage("class_group", "PASS";
  grh=false,
  invariants=join(class_invariants, ","),
  two_rank=class_2rank,
)

stage("unit_group", "start"; grh=false)
U, mU = unit_group(O; GRH=false)
unit_invariants = elementary_divisors(U)
stage("unit_group", "PASS";
  grh=false,
  invariants=join(unit_invariants, ","),
  rank=unit_group_rank(O),
)

# Both exact contexts above remain attached to O.  Hecke's p-Selmer routine
# asks for them again internally, and therefore reuses the proved class and
# unit maps rather than starting its default GRH-dependent computation.
ideal_type = typeof(ideal(O, 2))
S = ideal_type[]
for p in selmer_rational_primes
  append!(S, collect(keys(factor(p*O))))
end
unique!(S)
stage("number_field_pselmer", "start"; prime_ideal_count=length(S))
Sel, mSel = pselmer_group(2, S; algo=:raw)
selmer_invariants = elementary_divisors(Sel)
@assert all(==(2), selmer_invariants)
stage("number_field_pselmer", "PASS";
  dimension=length(selmer_invariants),
  prime_ideal_count=length(S),
)

q_support = vcat(ZZ(-1), selmer_rational_primes)
QSel, mQSel = pselmer_group(2, q_support)
norm_images = [preimage(mQSel, norm(mSel(g))) for g in gens(Sel)]
norm_hom = hom(Sel, QSel, norm_images)
norm_kernel, into_selmer = kernel(norm_hom)
norm_snf, into_kernel = snf(norm_kernel)
norm_invariants = elementary_divisors(norm_snf)
@assert all(==(2), norm_invariants)

stage("norm_kernel", "PASS";
  dimension=length(norm_invariants),
  rational_target_dimension=length(elementary_divisors(QSel)),
)

for (index, generator) in enumerate(gens(norm_snf))
  alpha = mSel(into_selmer(into_kernel(generator)))
  alpha_norm = norm(alpha)
  norm_is_square, _ = is_square_with_sqrt(alpha_norm)
  @assert norm_is_square
  coefficients = [coeff(alpha, i) for i in 0:2]
  println(
    "$(PROTOCOL)|stage=norm_basis|status=PASS|index=$(index)" *
    "|coefficients=$(join(coefficients, ","))|norm=$(alpha_norm)"
  )
end
flush(stdout)

stage("complete", "PASS";
  global_s_squareclass_dimension=length(selmer_invariants),
  norm_square_envelope_dimension=length(norm_invariants),
  class_group_two_rank=class_2rank,
)
