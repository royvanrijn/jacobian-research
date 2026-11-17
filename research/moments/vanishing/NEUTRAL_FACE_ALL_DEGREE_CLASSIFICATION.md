# All-degree classification of an affine neutral Gaussian face

15 November 2026. Internal written theorem depending on the independently
unreviewed [logarithmic curve criterion](LOGARITHMIC_CURVE_PERIOD_CRITERION.md)
and the retained one-radial Gaussian nullcone theorem. The
[checker](verify_neutral_face_classification.py) checks exact residue,
normalization and coefficient identities; it does not certify the analytic
contour argument. No independent/external review, formal verification or
novelty claim is made.

## 1. A full weakly one-sided architecture, in every degree

Use E(ZW)=E(T^2)=1 and all other Wick pairings zero, and assign circular
weight a-b to Z^a W^b T^c. Let P be any complex polynomial whose monomials
all have nonnegative weights, or all have nonpositive weights. Higher powers
of both Z and W are allowed. Assume only that its weight-zero part is
affine in U=ZW:

    P_0=U*H(T)+r(T),       H,r arbitrary complex polynomials.       (1)

**Theorem.** The following are equivalent:

1. E(P^m)=0 for all sufficiently large m.
2. E(P^m)=0 for every positive m.
3. H=r=0, equivalently P_0=0.

Under these conditions, for every polynomial Q,

    E(QP^m)=0 whenever m>deg Q.                                 (2)

The constant one in this uniform mixed cutoff is sharp within the class.
For P=Z^2 W and Q=W^d, the equality-index value is

    E(QP^d)=(2d)! !=0       for every d>=1.

There is no degree bound on P,H,r or Q. This is a complete classification
of the stated weakly one-sided architecture, including eventual pure
vanishing. It is not a classification of polynomials containing both
positive and negative weights, or of arbitrary higher-radial-degree
neutral faces.

## 2. Project to the neutral face before integrating

Suppose all weights are nonnegative. Weight is additive, so the
weight-zero part of P^m is exactly P_0^m. Gaussian expectation kills every
other weight. The same statement holds when all weights are nonpositive.
Thus

    E(P^m)=E(P_0^m) for every m>=0.                             (3)

This is an all-index identity; no finite moment inference is involved.
Since U=ZW has E(U^j)=j! and is independent of T, expanding (1) gives
the formal pure exponential generating function

    I(t)=E exp(tP)=E_T[exp(t*r(T))/(1-t*H(T))].                  (4)

For example the coefficient indexed by j U-factors is
binomial(m,j)*j!/m!=1/(m-j)!, exactly the coefficient in (4).
The denominator is a t-adic unit. All coefficients are polynomials in T,
so coefficientwise Gaussian integration is legitimate without convergence.

## 3. Resolvent rigidity: the series (4) is polynomial only trivially

**Lemma.** For arbitrary complex polynomials H,r, the formal series

    M_(H,r)(t)=E_T[exp(t*r(T))/(1-t*H(T))]

belongs to C[t] if and only if H=r=0. In that case it is identically one.

### Nonconstant H: every moving residue decays

Let d=deg H>=1. Over C(t), the polynomial 1-tH(x) is squarefree: a common
root with its x derivative would satisfy H'(x)=0 and H(x)=1/t, which
is impossible for transcendental t. Remove the finitely many exceptional
specialized base values. Its roots alpha_i(t) are algebraic and approach
the finite roots of H as t tends to infinity.

If a is a root of H of multiplicity k, write

    H(x)=c*(x-a)^k*(1+O(x-a)),       c!=0.

After a finite Puiseux cover,

    alpha_i-a=C_i*t^(-1/k)*(1+o(1)),       C_i!=0,
    Res_(x=alpha_i) dx/(1-tH(x))
        =-1/(t H'(alpha_i))=O(|t|^(-1/k)).                      (5)

As k<=d, every residue is O(|t|^(-1/d)). The differential has at most a
simple pole at infinity because d>=1. Consequently set

    rho=t^(-1/d),
    omega=t^(1/d)*dx/(1-tH(x)),
    Phi=t*r(x)-x^2/2.                                          (6)

Then omega is logarithmic at every puncture, including infinity, with
bounded algebraic residues. The normalization rho decays. The criterion
therefore bounds every zero-exponential part of the original period by

    O(|t|^(-1/d)*(1+log|t|)).                                   (7)

Its proof excludes every nonzero polynomial period, as well as constants:
a nonzero polynomial has a nondecaying zero-exponential asymptotic, and
the separated nonzero-action blocks cannot supply that polynomial. This
is a direct consequence of the growth estimate, not a new convergence
assumption about (4).

### The full formal connection and its moving-form derivative

Put R=C(t)[x,1/(1-tH(x))] and

    V=R dx/(d_x+d_x Phi)(R).

This is finite-dimensional twisted de Rham cohomology. Define the formal
horizontal functional

    J_0(f dx)=E_T[f(t,T)*exp(t*r(T))].                           (8)

Every localized denominator has constant t coefficient one, so (8) is
well-defined coefficientwise in C((t)). Gaussian integration by parts
annihilates every exact form. The actual tangent connection is

    nabla_0[f dx]=[(f_t+r(x)*f) dx].

In particular it differentiates 1/(1-tH(x)) as well as the phase. For
(6), pass to a finite algebraic base extension, take

    nabla=nabla_0+(rho'/rho) id,
    J=rho*J_0,          v=[omega].                              (9)

Then J is horizontal and J(v)=M_(H,r). The scalar gauge does not discard
the moving-pole terms: omega itself depends on t and must be differentiated
as a form. All formal series may be read in a finite Puiseux Laurent field.

If M_(H,r) were a polynomial f(t), its values on the cyclic derivatives
of v would be f,f',f'',... . These all lie in the coefficient field, so
J restricts to a rational horizontal functional on the cyclic submodule.
At an ordinary point extend that functional to the full fibre and solve
the dual analytic connection locally. The retained Hien--Roucairol
comparison in the curve criterion gives a flat cycle with period f(t).
Here Phi is nonconstant in x over C(t), since its x^2 coefficient includes
the fixed nonzero Gaussian contribution -1/2; projection and phase are
algebraically independent. Equation (7) and action separation prohibit a
nonzero polynomial f. But M_(H,r)(0)=1, so f cannot be zero. This proves
the lemma whenever H is nonconstant.

### Constant H: use univariate Gaussian nullity

For H=c, polynomiality of M implies

    E exp(t*r(T))=(1-ct)*M(t) in C[t].                          (10)

The retained one-radial nullcone theorem FWT1 implies that no nonzero
univariate Gaussian polynomial s(T) has all positive moments zero.
Indeed add an unused independent real Gaussian S and use circular
coordinates A=(T+iS)/sqrt(2), B=(T-iS)/sqrt(2). If deg s=e>=1, the
highest homogeneous terms of s((A+B)/sqrt(2)) contain both A^e and B^e
with nonzero coefficients, so their support is not strictly one-sided.
A nonzero constant also has a nonzero moment. The rank-one specialization
of FWT1 therefore forces s=0 whenever every positive moment vanishes.
This is also the immediate consequence of Wilson's retained Theorem 1.

If the left side of (10) is polynomial of degree N, choose an integer
k>N. Every positive moment of s=r^k is zero, because its indices are
km>N. Hence r=0. Now M=1/(1-ct), which is polynomial only when c=0.
This proves the resolvent lemma in every case.

## 4. Complete proof and the entire trapped W-linear branch

Eventual pure vanishing of P is equivalent to polynomiality of I(t).
Sections 2--3 force H=r=0, proving 1=>3 in the theorem. The remaining
implications and (2) follow from strict weight: every monomial of P^m
has absolute weight at least m with the fixed common sign, while every
monomial of Q has absolute weight at most deg Q. The displayed sharp
example follows from the circular Wick identity.

An immediate corollary closes a broad formal branch left by the earlier
[squarefree linear-T theorem](SQUAREFREE_LINEAR_T_ALL_DEGREE.md).
Let h(Z,T),G(Z,T) be completely arbitrary polynomials with h(0,T)=0, and

    P=W*h(Z,T)+G(Z,T).

Its weights are all nonnegative. Its neutral face is

    P_0=ZW*H(T)+r(T),
    H(T)=h_Z(0,T),       r(T)=G(0,T).

Thus eventual pure vanishing, all positive pure vanishing, and

    h divisible by Z^2,       G divisible by Z                  (11)

are equivalent. Every such member satisfies the sharp uniform mixed
cutoff m>deg Q. This covers arbitrary degrees in both Z and T, including
all constant and repeated-root linear-T branches with a(0)=b(0)=0.
The formal fixed point is g=0; (4) is also exactly GL1 on that branch,
but the neutral-face projection proves it without localizing away g.

## 5. Literature, assurance and remaining scope

The only one-variable algebraic input in the constant-H case is the
already retained [Long theorem FWT1](../../foundations/extended-geometry/FACTORIALLY_WEIGHTED_MULTITORUS_THEOREM.md).
The selected primary [Wilson v1 body](https://arxiv.org/html/2607.23887v1),
Theorem 1 and its weight lemma, was revisited for this corollary. Its
record still identifies v1 as the available arXiv version. A newly
rechecked [Liu record](https://zenodo.org/records/21674142) points through
its latest-version link to the [quadratic v3 record](https://zenodo.org/records/21688251).
Only those record descriptions were read in this pass; no PDF theorem,
producer package or finite-threshold calculation was imported. The
quadratic classification and bounds already have separate retained notes.
These selected checks do not claim exhaustive current literature coverage.

The nonconstant-H proof imports the internally written moving-puncture
logarithmic criterion, including its inverse-cut action separation and
uniform growth argument. Its independent analytic-review obligation remains
open. Exact residue/normalization checks do not remove that obligation.

The completed scope now includes an affine-radial neutral face and the
entire h(0,T)=0 W-linear branch in arbitrary degrees. Do not recreate
their low-degree subcases. Remaining structural questions concern neutral
radial degree at least two, genuinely mixed-sign angular supports, and
the constant/repeated-root linear-T branches with h(0,T) nonzero. Their
classification may include both safe and failing members. The unrelated
least cubic prefix, expanded radical cofactors and earlier completed
witnesses are unchanged. [Source audit](neutral_face_source_audit.json).
