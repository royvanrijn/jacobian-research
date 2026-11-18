# Uniform denominator restrictions for degree-twenty-six trace parameters

A positive-natural-density set of primes is forbidden in the
denominators of the canonical trace parameters s,t of every rational
point on the degree-twenty-six accessory curve. At each of those
primes, every affine Qp point already has both parameters integral.
Prime43 gives an exact joint root-free certificate for this set.

This is a uniform arithmetic restriction, not a rational-point
obstruction. The curve has52 affine points modulo43, including47
smooth ones, so it has infinitely many smooth good Q43 points.
No integrality of the reconstructed numerator coefficients is inferred.

## 1. A general Newton denominator criterion

Let
\[
 F(s,t)=\sum_{i=0}^d f_i(t)s^i\in\mathbb Z[s,t],\qquad
 \deg f_i\le D-i,
\]
and define its leading diagonal polynomial
\[
 A(x)=\sum_{i=0}^d[t^{D-i}]f_i(t)\,x^i.
\]
Assume A(0) and its degree-d coefficient are nonzero. Write
L(t)=f_d(t). For a prime p not dividing those two integers, its
degree remains d and L has degree D-d modulo p.

**Theorem 1.** At every affine Qp point of F=0, precisely the
following denominator alternatives are possible:

1. Both s,t are p-integral.
2. v_p(s)=v_p(t)<0, and the unit s/t reduces to a root of A modulo p.
3. v_p(s)<0<=v_p(t), and t reduces to a root of L modulo p.

In particular, if A and L have no roots over Fp, every affine Qp
point has s,t in Zp.

**Proof.** Suppose first v_p(t)<0. If v_p(s)>=v_p(t), put x=s/t
and v=1/t. The row-degree hypothesis makes
\[
 G(x,v)=v^D F(x/v,1/v)\in\mathbb Z_p[x,v],\qquad G(x,0)=A(x).
\]
Here x is integral and v lies in pZp, so A(xbar)=0. Since A(0)
is a unit, xbar is nonzero; hence v_p(s)=v_p(t).

If instead v_p(s)<v_p(t)<0, every monomial s^i*t^j has valuation
at least
\[
      Dv_p(t)+i\bigl(v_p(s)-v_p(t)\bigr).
\]
This bound is strictly smallest for i=d. Within that row the
leading monomial s^d*t^(D-d) is uniquely smallest, and its
coefficient is a p-unit. It cannot cancel. This case is impossible.

Finally suppose t is integral and s is not. If L(t) were a unit,
the term L(t)*s^d would have strictly smaller valuation than every
lower s-degree term, again impossible. Thus L(tbar)=0. The cases
include zero coordinates with v_p(0)=infinity. This proves all
alternatives and the root-free consequence. ∎

No assertion that the permitted alternatives actually occur is made.
The proof uses only a discrete valuation. It works over any discretely
valued field when the coefficients are integral, the two diagonal
constants are units, and the residue-polynomial hypotheses are taken
over its residue field. Completeness and a finite residue field are
unnecessary.

## 2. The exact degree-twenty-six leading data

Let F be the frozen primitive integer polynomial of the
[high-genus accessory theorem](Q26_HIGH_GENUS_AFFINE_TRACE_CURVE.md).
Its [complete coefficient input](q26_high_genus_affine_trace_curve_inputs.json)
has
\[
                d=12,\qquad D=144,\qquad \deg f_i=144-i.
\]
The following exact coefficient identity holds:
\[
        L(t)=25\,B(t)^2\,(t^2-t-1)^{10},\qquad \deg B=56.
 \tag{1}
\]
Every coefficient of the integer polynomial B is retained in the
[new input](q26_trace_parameter_denominator_inputs.json), and the
checker verifies the full product(1). Irreducibility of B is neither
assumed nor needed. The two leading-diagonal constants are
\[
 A(0)=2^3 3^4 13^2 17^2,
 \qquad [x^{12}]A=3^6 5^2 7^2 13^2 17^2.                 \tag{2}
\]

Thus Theorem1 applies outside E={2,3,5,7,13,17}. For every rational
point, a prime p outside E dividing the denominator of t divides
the denominator of s to exactly the same order, and A has a
nonzero root modulo p. If only s has a p-denominator, the integral
residue of t is a root of B(t)*(t^2-t-1) modulo p.

## 3. An exact joint43 certificate

**Proposition 2.** A, B and G=t^2-t-1 all retain their degrees
modulo43, are separable there, and have no F43 roots. Their leading
coefficients modulo43 are36,16,1, while A(0)=22 modulo43.
Consequently every affine Q43 point of F=0 has s,t in Z43.

**Proof.** For each of the three polynomials, the exact checker
verifies both
\[
       \gcd(\overline f,\overline f')=1,
       \qquad \gcd(\overline f,X^{43}-X)=1.
\]
It independently evaluates every one of the43 residue arguments
and retains all nonzero values. The full leading coefficient
factorization(1) makes L root-free as well. Theorem1 gives the
conclusion. ∎

A complete direct census of the43^2 affine residue pairs gives52
points,47 with nonzero gradient. For example(s,t)=(32,1) has
gradient(25,4) modulo43. Thus the integrality theorem is compatible
with infinitely many local points; it does not assert Q43
nonexistence. Hensel lifting and the finite bad-complement argument
of the [local-branch theorem](Q26_SMALL_PRIME_LOCAL_BRANCHES.md)
give infinitely many good Q43 models. The residue point itself
is not assumed to have good critical-pair reduction.

## 4. Positive density of forbidden parameter-denominator primes

**Theorem 3.** There is a set of primes of positive natural density
such that s and t are p-integral for every rational point F(s,t)=0.
More strongly, every affine Qp point is integral at each such prime.

**Proof.** Let M be the joint splitting field over Q of A, B and G.
The unit and separability checks at43 imply that each splitting
field is unramified there, hence their compositum is unramified.
For the nonmonic factors, an integral monic root rescaling is
legitimate because their leading coefficients are43-units.
The standard reduction/Frobenius correspondence in
[Milne ANT3.00, Corollary8.22 and Theorem8.23,p138](sources/milne-ant300.pdf)
shows that the Frobenius conjugacy class at43 fixes no root in
any of the three actions.

The same joint conjugacy class at any other unramified prime,
away from finitely many leading/discriminant primes, has the
same no-fixed-root property simultaneously for A,B,G.
No independence between the three splitting fields is assumed.
By the natural-density form of
[Milne's Chebotarev theorem, Theorem8.31,p140](sources/milne-ant300.pdf),
this class occurs with positive density. Removing finitely many
primes does not change that density. At every remaining prime
Theorem1 and(1) force both parameters integral. ∎

No explicit density or full Galois group of M is computed. This
does not prove integrality away from a finite set: primes where
the leading polynomials have residue roots can still support
denominators. It is a restriction on s,t, not automatically on
the reconstructed q or N coefficients. Their rational formulas
can have poles near recurrence-exception residues even when
s and t are integral.

## 5. No rational normalization point lies outside this affine chart

Let C be the smooth projective normalization of F=0. Its affine
coordinate boundary here means the poles of s or of t. It does
not mean poles of the ratio s/t at affine points with t=0.

**Corollary 4.** C has no Q43 point at a pole of either coordinate,
and hence no Q-rational point on that boundary.

**Proof.** An integral polynomial with unit leading coefficient has
only integral Q43 roots: a root of negative valuation would make
its leading term uniquely smallest. Proposition2 therefore implies
that neither A nor L has a Q43 root.

At a Q43-rational point P of the base-changed smooth curve, the
local ring is a discrete valuation ring with residue field Q43.
Every nonzero constant is a unit. Apply the general valuation
criterion to its function field and the valuation at P. The
residue polynomials are A and L over Q43, both root-free, so
s and t have nonnegative order at P. Thus neither has a pole
there. Every Q-rational point embeds into C(Q43), giving the
final claim. ∎

The singular rational point(1:0:0) of the naive degree144 plane
closure is not a counterexample: it has no Q43 branch on the
normalization. A rational plane point at infinity does not itself
give a rational normalization point.

## Evidence and remaining question

The [checker](verify_q26_trace_parameter_denominators.py),
[inputs](q26_trace_parameter_denominator_inputs.json),
[receipt](q26_trace_parameter_denominator_verification.json) and
[source ledger](q26_trace_parameter_denominator_source_ledger.json)
retain all leading identities, exact root-free/separability checks,
every residue value and the complete52-point census. The general
valuation proof and Chebotarev import supply the uniform and
positive-density statements; a sample of primes does not prove them.

The rational-point set remains finite and undetermined. No rational
point, nonexistence theorem, low-genus quotient, Brauer-Manin
obstruction, integrality of all model coefficients, formal
verification, independent external replay, external human review
or novelty claim is established.
