# Global Euler images, Gaussian pullback lattices and phase-power pole packets

16 November2026. Internal all-index proofs; no novelty or independent,
formal or external assurance is claimed. The
[time-zero theorem](TIME_ZERO_SATURATION_AND_PRINCIPAL_PART_RECONSTRUCTION.md)
restores full Gaussian sources and actual images by one principal-part map.
Here finite point-module tests close the global-lattice questions for two
entire previously retained architectures and give an exact phase-power
image transformation. The [checker](verify_gaussian_origin_and_global_lattices.py)
and [source audit](gaussian_origin_global_lattice_source_audit.json) retain
finite identities, not a large integration/kernel campaign.

## 1. A scalar test for the global intermediate extension

Let D=K[t]<partial_t>, characteristic zero. For a scalar polynomial
operator L whose leading derivative coefficient has no zeros except0,
consider D/DL. If L acts injectively on the LEFT point module delta_0,
there is no quotient delta_0. If right multiplication by L is injective
on the RIGHT point module D/tD, then the left ideal DL is t-saturated:

    tP=QL => Q in tD => P in DL.

The implication uses reduction modulo tD followed by cancellation of t
in the Weyl domain. Thus D/DL has no submodule supported at0 either.
At other finite points it is an ordinary connection, so there are no
point submodules or quotients there. The
[retained finite-length criterion](RATIONAL_CORES_AND_BOUNDARY_SUBSPACES.md)
identifies it as the intermediate extension across0 of its generic
connection. No holonomic-dual computation is needed for these tests.

If a selected formal series f has minimal order equal to ord L, the
surjection D/DL -> Df has only a possible punctual kernel. The saturation
test excludes it. This gives the EXACT global annihilator Ann_D(f)=DL,
not just a minimal rational differential equation.

## 2. Complete scalar Euler-unit classification

Take ANY polynomial b(z) in K[z] of degree r>=1 and put

    theta=t partial_t, A_b=1-t b(theta), L_b=partial_t A_b,
    h_b=sum_(m>=0) a_m t^m,
    a_0=1, a_(m+1)=b(m)a_m.

Then A_b h_b=1 and L_b h_b=0. The weighted proof in the
[squared-spectator note](SQUARED_SPECTATOR_INVARIANCE_AND_UNBOUNDED_GAUSSIAN_IMAGES.md)
works unchanged: with weight r*ord_t(coeff)-theta-order, the initial
symbol is the prime binomial 1-b_r T xi^r. Therefore A_b is irreducible
over K((t)) and over K(t), with slope1/r and rank r.

On delta_0 with basis v_j=partial_t^j v_0, direct normal ordering gives

    A_b v_j=v_j+j b(-j-1)v_(j-1),
    L_b v_j=v_(j+1)+j b(-j-1)v_j.                          (1)

On the right module D/tD the corresponding formulas are

    partial_t^j A_b = partial_t^j-j b(j-1)partial_t^(j-1),
    partial_t^j L_b = partial_t^(j+1)-(j+1)b(j)partial_t^j.
                                                               (2)

Missing negative indices are0. A_b is unipotent triangular in both
tests; L_b has nonzero highest-degree term in both. Hence

    S_b=D/DA_b is a GLOBAL SIMPLE module,
    N_b=D/DL_b is the intermediate extension of its connection.

Simplicity of S_b follows from generic irreducibility and absence of all
punctual submodules/quotients. S_b is already t-invertible. Indeed its
cyclic vector u satisfies

    u=t^k product_(j=0..k-1)b(theta+j)u, for every k>=1.

This makes t surjective on every derivative of u, and saturation gives
injectivity. Its Laurent companion basis is u,theta u,...,theta^(r-1)u.

There is an EXACT global sequence

    0 -> K[t] -> N_b -> S_b ->0,  1 |-> [A_b].             (3)

Right cancellation of A_b proves the first injection; taking the quotient
by [A_b] gives S_b. The formal map [1]->h_b sends [A_b] to1.
It has the following complete actual-image alternative.

* If b has a nonnegative integer root n, let n be the first. Then h_b
  is a polynomial of degree n. The sequence(3) splits, and its actual
  formal image is K[t], with kernel one copy of S_b. Explicitly, if
  q=h_b and e=[1], the vector e-q[A_b] is killed by A_b and supplies
  the splitting.
* If b has NO nonnegative integer root, h_b has exact Gevrey order r
  and radius0 over C. Indeed |b(m)| is bounded above and below by
  nonzero constants times m^r for large m. In particular h_b is not
  rational. Rationality over any characteristic-zero extension descends
  by finite-rank linear recurrence/Hankel linear algebra. There is no
  rational particular solution A_b q=1, since its Laurent expansion
  would be the UNIQUE h_b. Thus(3) is nonsplit even generically.
  Since S_b is simple, its unit submodule is essential. The actual
  period kernel intersects that unit trivially and must therefore be0:

        D h_b = D/D(partial_t A_b), EXACTLY.               (4)

This distinguishes the genuine selected series from the larger source
in the terminating case. It does not infer actual rank from a nonminimal
equation. It is a complete global classification of this scalar Euler
family, not a claim every Gaussian phase is of this form.

## 3. Full global structure of the squared pole-packet images

Let b have no nonnegative integer roots and let rho_1,...,rho_s be
distinct nonzero scalars. Set A_i=1-rho_i*t*b(theta) and let h_i be
its unique formal solution A_i h_i=1. Distinct leading prime symbols
give pairwise nonisomorphic simple S_i=D/DA_i. The EXACT presentation
of M=sum D h_i has generators e_0,e_1,...,e_s and relations

    partial_t e_0=0, A_i e_i=e_0.                          (5)

To prove faithfulness, its source is an extension of K[t] by the direct
sum of the S_i. Each single summand extension is nonsplit by section2.
A nonzero kernel disjoint from the unit would contain one of these
simple summands and split its pullback extension, a contradiction.
This proves the actual presentation without computing a source-image
kernel. All unit and scalar operator identities are retained explicitly.

Consequently the earlier generic statement strengthens globally to

    0 -> K[t] -> M -> direct_sum_i S_i ->0.                 (6)

The quotient is now proved GLOBAL semisimple; the unit is essential,
the extension nonsplit, and M has Jordan--Holder length1+s. It is its
own intermediate extension across0. Since each S_i is t-invertible,
localizing(6) gives

    M[t^-1]/M = K[t,t^-1]/K[t] = delta_0.                 (7)

Thus the full fixed-connection lattice problem on the punctured affine
line has ONLY TWO embedded extensions between its intermediate and
meromorphic extensions: M and M[t^-1]. Only M lies in the specified
formal power-series embedding. The source boundary has one unit
direction regardless of rank. The characteristic cycle is

    (1+r*s)[zero section] + (r+1)*s[conormal at0].          (8)

For S_i, this follows from the principal leading symbol
-rho_i*b_r*t^(r+1)*xi^r; additivity and the unit give(8).
This conormal multiplicity includes irregularity and is NOT a count
of boundary quotients M/j_!*M, which are0 here.

Taking b=b_ell from the squared-spectator theorem closes its formerly
unclaimed global t=0 lattice and global quotient structure for EVERY
ell and EVERY finite pole packet. The old sealed note is preserved;
its cautious generic-only boundary is strengthened here.

## 4. The earlier ternary even-rank tower selects a hyperplane glue

The prior
[rational Gaussian pullback construction](RATIONAL_GAUSSIAN_PULLBACKS_AND_IRREGULAR_RANK.md)
ALREADY proves unbounded actual ranks2nu in Gaussian dimension3. It is
not a new discovery of this note or the squared-spectator construction.
The latter's additional invariant is exact ordinary-scheme equality
with changing rank. Here the new result is the old tower's GLOBAL lattice.

Retain its parameters nu>=1, beta,r0 nonzero, kappa=beta*nu/r0,
delta=a(0)/beta and polynomial h(Z,T)=a(Z)+b(Z)T, with
b=beta*(1-Z/r0)^(nu+1) and the stated polynomiality condition on a.
Write

    f_j=Y_(j/nu,delta), Q_j=(1-Z/r0)^j, 1<=j<=nu.

The old all-order actual identity gives Phi[Q_j]=f_j, with f_j(0)=1.
Every f_j has exact minimal order2 and operator

    L_alpha=kappa^2*t^3*D_t^2
      +(2(alpha+1)kappa^2*t^2-kappa*delta*t-1)D_t
      +alpha(alpha+1)kappa^2*t-alpha*kappa*delta,
    alpha=j/nu.                                         (9)

The old generic decomposition has rank2nu and separates these nu
blocks by their distinct regular-singular exponent classes at infinity.
Its full connection is ordinary on the punctured affine line.

Let N_j=D f_j. For EVERY alpha, left action on delta_0 is

    L_alpha v_m=-v_(m+1)+kappa*delta*(m+1-alpha)v_m
      -kappa^2*m*(m-alpha)*(m+1-alpha)v_(m-1).             (10)

Right action modulo tD is

    D_t^m L_alpha=-D_t^(m+1)-kappa*delta*(m+alpha)D_t^m
      +kappa^2*m*(m+alpha)*(m+alpha-1)D_t^(m-1).          (11)

Both are injective by their top term; the right cokernel is one-dimensional.
Section1 and the already proved minimal order give N_j=D/DL_alpha
EXACTLY. Each is an intermediate extension with N_j/tN_j=K.
Therefore

    J=direct_sum_j N_j

is the intermediate extension of the full actual generic connection,
and L=J[t^-1]=M[t^-1]. Here M denotes the full Gaussian image.
Each N_j[t^-1]/N_j is delta_0 by the principal-part argument:
t-injectivity and N_j/tN_j=K give successive quotients t^aN_j/t^(a+1)N_j=K;
the actual constant-one period makes the maps to K[t]/t^k isomorphisms
for every k. Thus the proof applies to these individual cyclic modules,
without assuming that they are themselves full Gaussian sources.
Thus L/J=delta_0 tensor K^nu. The ACTUAL scalar principal-part map on
these summands is SUM of coordinates. The universal origin theorem gives

    M/J=delta_0 tensor ker(sum:K^nu->K)=delta_0^(nu-1).
                                                               (12)

The full image has the explicit GLOBAL generating set

    f_1,...,f_nu, (f_j-f_nu)/t for 1<=j<nu.              (13)

The divided differences in(13) are original POLYNOMIAL marker periods,
not new formal amplitudes. Set R_j=(Q_j-Q_nu)/Z, a polynomial, and use
the exact W connection of the old phase P=W*h+G:

    nabla_W(-R_j)=Q_j-Q_nu-t*h*R_j.

Hence (f_j-f_nu)/t=Phi[h*R_j]. These original markers realize every
boundary direction in(12). Point-module derivatives supply the full
tower, so(13) generates every original amplitude image. This uses the
retained full generic upper bound, not only the unit or selected periods.

For nu>1 the unit is NOT essential in this genuine pure Gaussian image.
It belongs to N_nu by the retained inhomogeneous alpha=1 identity.
For j<nu, N_j has nonintegral infinity exponent class j/nu and cannot
contain a nonzero rational function; in particular N_j intersect K[t]=0.
Nevertheless the full M is indecomposable and has a unique normalized
formal period. For indecomposability, M/tM is one-dimensional, whereas
every nonzero submodule of K[[t]] has a least occurring t-order and hence
a nonzero quotient by t. A direct-sum decomposition would contradict
that dimension. The boundary
glue(12) joins the generically separated blocks. It is nonsplit because
M has no t-torsion. Thus neither indecomposability nor unique formal Hom
implies essentiality of the unit in all pure Gaussian images.

This closes the entire existing pullback tower's global origin boundary,
not all repeated-root Gaussian phases. Its other-family hypotheses and
the smaller ternary quartic/cubic endpoints are unchanged.

## 5. Exact monomial recomposition of any rational pole packet

Let B be ANY pure-null Gaussian phase with full image
K[s,prod_(lambda in S)(1-lambda*s)^-1]. Fix k>=2 and put P=B^k.
The persistence of pure failures under ALL polynomial recompositions is
[already closed](POLYNOMIAL_RECOMPOSITION_REFLECTION_AND_UNIVERSAL_FAILURE.md);
it is not a new conclusion here. The new result is the full actual image.

Set S^[k]={lambda^k:lambda in S}, with duplicates removed, and

    a_m=(k*m)!/m!,
    b_k(z)=k*product_(j=1..k-1)(k*z+j),
    H_rho(t)=sum_m a_m*rho^m*t^m,
    A_rho=1-rho*t*b_k(theta).

Then the full image is EXACTLY

    M_(B^k)=sum_(rho in S^[k]) D H_rho.                   (14)

If S is empty it is K[t]. Otherwise (5)-(8) apply with r=k-1:

    rank_actual=1+(k-1)|S^[k]|,
    M/K[t]=direct_sum_(rho in S^[k])D/DA_rho,
    M[t^-1]/M=delta_0.                                   (15)

The quotient is global semisimple, the unit essential, and M is the
intermediate extension at0. Pole packets COLLAPSE under lambda->lambda^k;
their cardinality must not be multiplied blindly.

Here is the full-amplitude proof. Every original marked B-period has a
polynomial part plus terms p_lambda(n)lambda^n in its ordinary Taylor
coefficient, with polynomial p_lambda. Replacing B by B^k takes its
n=k*m coefficient and multiplies by (k*m)!/m!. Thus every marked
P-period is polynomial plus a sum of polynomial Euler operators on H_rho.
This gives inclusion in the right side of(14), including all amplitudes.

To show each rho survives, choose a marked B-period with a pole in that
class. Derivative markers Q*B^j, 0<=j<k, give class coefficient polynomials

    product_(i=1..j)(k*m+i)
       *sum_(lambda^k=rho)lambda^j*p_lambda(k*m+j).       (16)

At the highest polynomial degree the inner leading coefficients form
the Vandermonde system (lambda^j). There are at most k distinct lambda
in the class; its first that-many rows are invertible. Thus some j
gives a nonzero polynomial. A nonzero polynomial in theta cannot kill
H_rho modulo K[t]: its weighted initial monomial is not divisible by
the prime symbol of A_rho. In the global semisimple quotient of the
right side, every simple rho-direction is therefore attained. Its unit
is attained by purity. This proves equality(14), not merely a generic
rank bound or a claim about a single selected subsequence.

The classical hypergeometric label is

    H_rho = _kF_0(1,1/k,...,(k-1)/k;;k^k*rho*t).

Its defining series and Euler equation are cross-checked with
[NIST DLMF16.2.1](https://dlmf.nist.gov/16.2#E1) and
[16.8.3](https://dlmf.nist.gov/16.8#E3). These are formal divergent
series for k>=2, not convergent Wick integrals. Rank minimality,
faithfulness and global lattice are proved above, not imported from
a displayed generalized hypergeometric equation.

For the retained ternary quartic B, S={1}. Every k>=2 therefore gives
a pure-null phase B^k in THREE Gaussian variables with EXACT actual rank k,
global cyclic annihilator D(partial_t A_1), and a linear bad marker Z.
The generating h equals F_(1+Z); its coefficient is (k*m)!/m!.
Thus every integer rank>=2 is realized in dimension3, alongside the
already known even-rank tower. No minimal-degree claim or degree sweep
is made; B^k has degree4k.

## 6. Remaining scope and credits

The universal principal-part hyperplane and these two complete global
lattice descriptions are general structural progress. They do not prove
all actual Gaussian images are intermediate extensions: the pullback
tower itself has the nonzero boundary(12). Nor is the unit universally
essential. Regular nonrational realizations outside the retained classes,
exhaustive safe normal forms and general actual generic connections remain.

The older ternary even-rank tower, factorial recomposition failure and
general boundary-subspace theorem are explicitly credited as prior work.
The new content is origin reconstruction, exact global presentations,
the selected hyperplane glue, and the full phase-power pole-packet image.
No prior assurance, endpoint or capped job is changed. Oaku identities
and torsion/specialization data remain sealed. Sabbah/Maisonobe retain
convergent rational-slope credit; arbitrary nonconvergent formal/all-ring
extension priority remains separately unassessed. The broad goal stays active.
