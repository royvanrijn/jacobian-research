# Rank-one actual images and nilpotent cutoff obstructions

15 November 2026. Internal written deductions; independent review remains
outstanding. No novelty claim is made. General integration-module and
rationality algorithms are **not implemented or executed**. The
[checker](verify_rank_one_images.py) runs exact scalar, pole-isolation,
Bernstein identity and finite Gaussian transcription controls. It does
not certify the universal proofs.

## 1. A structural trichotomy for the actual image

Let H be a holonomic module over D=K[t]<partial_t>, and let a D-linear
formal period map have image

    M subset K[[t]].

Assume 1 belongs to M. This is automatic for the retained full Gaussian,
full SIC and Gaussian--Laurent GVC period maps when the original positive
pure moments vanish. In the polynomial-valued cases, keep the original
output x independent and use coefficient field K(x). All rank statements
below use that field, rather than evaluation at x=0.

There is an effective theoretical trichotomy over exactly presented
effective characteristic-zero fields:

1. M=K[t]: every full-amplitude image is polynomial and the original
   all-multiplier eventual property holds.
2. M=K[t,d(t)^(-1)] for a nonconstant squarefree d with d(0)!=0:
   the actual generic rank is one, but finite poles cause failure.
3. M contains a nonrational series: its actual generic rank is at least
   two, and the all-multiplier property fails.

Compute a finite **polynomial Weyl-module** generating set of H, not a
rational basis after inverting t-polynomials. Test each generator image
for rationality by section 3. If all are rational, their exact pole
denominators determine d; d=1 is case 1. Otherwise case 3 holds.
For the original eventual property, the equivalence between polynomiality
of the full image and all original multipliers is the retained
[Gaussian criterion](FINITE_MULTIPLIER_MIXED_DECISION.md) or
[operator amplitude-transfer theorem](CONSTANT_COEFFICIENT_OPERATOR_DECISION.md).

For GVC this trichotomy concerns the **full localized/integrated Laurent
period image**. Rationality of just all original diagonal H_Q is not
substituted for rationality of every full module generator: Laurent
index shifts can introduce integrations in time. Polynomial eventuality
transfers by the earlier proof, while rationality requires its own test.

Thus actual rank one is completely classified. This does not determine
the full kernel of H, compute arbitrary higher actual rank, or classify
the remaining nonrational images. Nonrational need not mean irregular:
log(1-t) is a regular-singular nonrational control.

## 2. Every rank-one image containing 1 is a finite localization

The embedding in K[[t]] makes M torsion-free over K[t]. Its localization
therefore injects into K((t)). If dim_(K(t)) M[t-polynomials^(-1)]=1,
its element 1 spans that localization. Consequently every element of
M is rational in t. Conversely, a nonzero image contained in K(t) has
generic rank one.

Since H is finitely generated over D, so is its image. Choose finitely
many rational image generators. Their denominators have finitely many
poles; let S be their union over an algebraic closure of K.
No image is singular at t=0 because M subset K[[t]], so 0 is not in S.

Work first over an algebraically closed field. If f in M has a pole
of order k at a, multiply f by

    (t-a)^(k-1) product_(b!=a)(t-b)^(pole_order_b(f)).

The result has exactly one simple pole, at a, plus a polynomial.
Since K[t] subset M, subtract its polynomial part and divide by its
nonzero scalar residue to obtain 1/(t-a) in M. Its derivatives give
every higher pole order. Partial fractions now prove

    M=K[t, product_(a in S)(t-a)^(-1)].                      (1)

Finite extension and descent give the same statement over the original
effective K, with a squarefree polynomial d in K[t] whose roots are S.
Equivalently one can perform the isolation with coprime irreducible
denominators and Bezout identities. The scalar-extension argument is
valid because all spaces are K-vector spaces and extension is faithful.

The quotient is precisely the direct sum of principal-part modules:

    0 -> K[t] -> M -> direct_sum_(a in S) delta_a -> 0.       (2)

Here delta_a=D/D(t-a); send its generator to the class of 1/(t-a).
The derivative iterates map to nonzero multiples of all negative powers.
Every nonzero principal part generates this simple point module.
Different poles have disjoint supports, giving the direct sum.

No nonzero summand of (2) splits back into M. A nonzero point-module
generator in M would be killed by t-a, impossible in the rational
function field. More strongly, K[t] is the unique simple submodule and
is essential in M: every nonzero D-submodule contains a nonzero rational
function, clearing its denominator gives a nonzero polynomial, and
derivatives then give 1.

Over an algebraically closed field the D-module length is 1+|S| and,
when S is nonempty, the socle filtration has two steps, with socle K[t]
and semisimple quotient direct_sum delta_a. The characteristic cycle
has the zero section with multiplicity one and each pole conormal with
multiplicity one. This is a regular holonomic localization of the
trivial connection. There is no irregular rank-one pure image branch.

The generic connection has trivial monodromy, even when S is nonempty.
Generic rank, rationality and connection monodromy therefore all miss
the finite-point quotient in (2). Actual pole data decide safety.
The hypothesis 1 in M is indispensable: K[t]exp(t) has generic rank one
and an irregular singularity at infinity, but contains no nonzero
polynomial and is not of form (1).

## 3. A finite rationality certificate for a specified scalar series

Let f in K[[t]] be computable and let

    L=sum_(j=0)^r p_j(t) partial_t^j

be a supplied nonzero annihilator with p_j in K[t]. If r=0 then f=0.
Otherwise every finite pole of a rational solution lies among the zeros
of p_r: at an ordinary point, the highest derivative of a pole has
strictly largest pole order and cannot cancel.

For each such algebraic point a!=0 put

    mu_a=min_(p_j!=0)(ord_a(p_j)-j),
    B_a(z)=sum_(ord_a(p_j)-j=mu_a) lc_a(p_j) (z)_j.          (3)

The falling factorial (z)_j is evaluated at negative integers as well.
B_a is nonzero. A rational solution with a pole of order k at a has
leading term c(t-a)^(-k); the lowest-order term of Lf forces

    B_a(-k)=0.

Let k_a be the largest positive k allowed by this finite negative-integer
root test, or zero if there are none. Omit a=0 because the selected f is
a formal power series there. Define

    D(t)=product_(a!=0)(t-a)^k_a.                            (4)

Use irreducible-factor root fields to compute (3); conjugate points have
the same bound, so D may be taken in K[t], up to a nonzero scalar, and
D(0)!=0. This works over the effective finitely generated fields already
used in the retained endpoint tests.

Every possible rational f has Df polynomial. Conversely Df polynomial
implies f rational. Gauge the scalar operator by D to obtain a nonzero
polynomial-coefficient annihilator of g=Df, and apply the earlier
[scalar polynomiality test](FINITE_MULTIPLIER_MIXED_DECISION.md).
It uses the infinity degree bound plus the origin indicial zero test,
so it checks the entire series, including a delayed first nonzero term.

If g is polynomial, output f=g/D in reduced form and its actual poles.
If g is nonpolynomial, the finite denominator bounds and endpoint
certificate prove f nonrational. Apparent singularities of L can give
unneeded candidate factors, but reduction of g/D removes them; they
are not actual image poles.

This is the scalar endpoint version of classical rational-solution
methods. [Oaku--Takayama--Tsai, JPAA164(2001),199--220](https://www.lab.twcu.ac.jp/oaku/OakuTakayamaTsai2001.pdf),
section 3, bounds denominator exponents by local b-functions and reduces
rational solutions to polynomial solutions. Our one-variable proof (3)
is included directly and does not require running their multivariate
producer implementation. The general effective-field/module algorithm
here remains unimplemented.

## 4. One ordinary b-polynomial does not record the number of pole defects

If d is nonconstant and squarefree in one variable, its ordinary global
Bernstein-Sato polynomial is always

    b_d(s)=s+1,                                             (5)

irrespective of its roots or degree. Indeed choose A,B in K[t] with

    A d+B d'=1.

Then the finite differential operator

    P(s)=A(t)(s+1)+B(t) partial_t

satisfies P(s)d^(s+1)=(s+1)d^s. The factor s+1 is necessary:
evaluating an alleged equation at s=-1 would otherwise express a
nonzero multiple of d^(-1) as a polynomial differential operator on 1,
which has no finite pole.

Consequently (5) cannot distinguish a one-pole image from a many-pole
image, although their quotients in (2) and module lengths differ.
One must retain the divisor and the actual image, not only a scalar
b-polynomial. No assertion identifying the socle filtration above with
an unnormalized Hodge/twistor filtration is needed for this proof.

## 5. A genuine Gaussian family with nonstabilizing ordinary cutoff ideals

Use E(ZW)=E(T^2)=1 and parameter e. The retained affine-W theorem supplies

    P_e=W(1+eZ)-e(1+3eZ/2+e^2 Z^2/2)T^2,
    E exp(tP_e)=1,
    E[Z exp(tP_e)]=t/(1-et).                               (6)

These hold over Q[e] coefficientwise, including e=0, where P_0=W.
For e!=0 this is the retained quartic failure; its full actual image is
K[t,(1-et)^(-1)]. At e=0 the image is K[t]. The generic rank remains one,
while the pole t=1/e escapes to infinity and the module length drops
from two to one.

Let I_N subset Q[e] be the **full ordinary** ideal of all marked moments

    E(Z^a W^b T^c P_e^m),       m>N(a+b+c),

with all a,b,c>=0. Thus I_N is exactly the specified-cutoff ideal from the
[joint-source theorem](HOLONOMIC_SOURCES_AND_CUTOFF_LOCI.md).
Its complete all-N answer is

    I_0=(1),       I_N=(e^N) for every integer N>=1.          (7)

Here is the proof, without a finite search over markers. Give Z-minus-W
charge a-b. Every monomial of P_e has e-degree equal to its charge plus
one. Equivalently, over Q(e),

    P_e=e P_1(eZ,e^(-1)W,T).                                (8)

Wick contraction with Z^aW^bT^c therefore forces every surviving term
of a mixed m-th moment to have the **same** e-exponent

    m+b-a.                                                 (9)

If this exponent is negative, the moment is zero since P_e is polynomial
in e. A nonzero constant marker has only pure moments, already zero.
For a nonconstant marker and N>=1, the forbidden inequality gives

    m+b-a >= (N-1)a+(N+1)b+Nc+1 >= N.                       (10)

For N=1 the last inequality is still valid. Thus I_N subset (e^N).
The marker Z with m=N+1 attains the reverse inclusion by (6):

    E(Z P_e^(N+1))=(N+1)! e^N.

For N=0 the same marker at m=1 has value 1. This proves (7).

All closed cutoff loci for N>=1 are the same reduced point e=0.
Their ordinary ideals nevertheless decrease strictly forever, and
Q[e]/I_N has length and nilpotence index N. This is a genuine Gaussian
example of scheme nonstabilization, even though the safe point set is
closed and has common sufficient integer coefficient one.
It is not a counterexample to uniform bounds over complex field points.

There is a sharper scheme obstruction. Over the Artin coefficient algebra

    A_L=Q[e]/(e^L),       L>=1,

P_e is pure-null and every multiplier is eventually zero. Equations
(7)--(10) give the uniform bound m>L deg Q over A_L. It is sharp:
the marker Z at m=L has nonzero value L! e^(L-1).
Thus the least real or integer coefficient is exactly L, for every L.

The polynomial degree stays four (at most four for small L), and the
Gaussian dimension stays three. Hence there is **no fixed-dimension,
fixed-degree cutoff coefficient valid across all nonreduced coefficient
algebras**. Arbitrarily thick infinitesimal parameter directions require
arbitrarily large slopes. The open uniformity problem over complex
coefficients is distinct from this now-excluded scheme-uniform extension.
No new field-valued GMC counterexample is claimed.

## 6. The November revision of the filtration paper and the next boundary

The latest listed [Lőrincz--Yang revision,2609.05215v2](https://arxiv.org/html/2609.05215v2)
is **11 November 2026**, superseding the September v1. It removes the
regularity assumption in its algebraic simple-localization theorem.
Its weight-level algorithm 5.9 is stated for O_X localized at f; the
general theorem does not assert every actual period image has that form.
The revision also supplies a counterexample to the unrestricted
Torrelli equivalence, so the old positive phrasing cannot be imported.
Selected sections and exact boundaries are recorded in the source audit.

For the rank-one branch, (1) proves precisely the localization structure
needed to compare with that framework, with S=O on a curve.
For arbitrary higher-rank Gaussian images, neither simplicity of an
appropriate source nor equality with a simple localization has been
proved. No weight/Hodge producer package was run, no scalar b-root test
replaces the actual kernel/pole map, and no global coefficient-safe-locus
uniformity theorem follows from the paper.

The rational rank-one image branch and the nonreduced uniformity boundary
are now closed internally. Remaining image geometry starts with
nonrational actual rank at least two, including regular logarithmic
extensions as well as possible irregular pieces. The normalized
constant-b higher-rank rational-quotient problem is unchanged.

There is additional genuine-family structure to investigate: the universal
polynomial Gaussian and operator periods are jointly holonomic in their
coefficient parameters as well as time and sources, by the same polynomial
or localized rational exponential connection. The control
1+s*((1+t)^a-1) is not jointly D-finite in a: its successive a derivatives
contain independent powers of log(1+t), which is transcendental over
the rational function field. Thus that control omits this extra structure.
Joint parameter holonomicity has not here been proved to supply the
required specialization-valid degree bounds or global safe-locus
stabilization; actual period kernels at exceptional parameters still
require control.

Global safe-locus stabilization over C, valid operator compression,
concrete lower-power frontiers and independent proof review remain open.

[Source and assurance audit](rank_one_image_source_audit.json).

