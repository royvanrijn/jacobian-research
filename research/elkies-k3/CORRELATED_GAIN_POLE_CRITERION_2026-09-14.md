# Correlated gains: branch values, poles and the next height gate

Exploratory structural derivation, with elementary proofs below. This is not
a new MW17 construction, an independently reviewed theorem certificate, or
a change to [mathematical status](../MATH_STATUS.json). The
[positive objective](CORRELATED_QUADRATIC_GAINS_2026-09-12.md) remains open.
The existing Q80 polynomial genus-one exclusion is an inherited dependency,
including its written-proof verification boundary.

## What must coincide

Put K=Q(t), L=K(w), w^2=d. The usual character decomposition gives

```
E(L) tensor Q = E(K) tensor Q  direct_sum  E^d(K) tensor Q.
```

Thus the target requires rank at least two in one twist. In a biquadratic
extension, directions in the two different singleton characters do not
descend to its third quadratic subfield: the involution fixing that field
negates both directions. Taking their sum does not make it invariant.
This explains precisely why the exhibited pair in
[Elkies, section 3](https://arxiv.org/html/2608.25406v1) does not already
compress to the requested cover. This is the standard character argument,
also retained in [Theorem F4](RANK_MUTATION_AND_LIFT_THEOREMS.md#theorem-f4-multiquadratic-character-decomposition-and-base-genus).

There is no loss of rational rank in working with actual anti-invariant
points: P maps to P-sigma(P), whose kernel after tensoring with Q is the
inherited subgroup. This operation can increase height. It does not assert
that a chosen nonzero trace can be halved over K, or that a small-height
anti-invariant representative always exists.

At a smooth branch fibre an anti-invariant point specializes into E[2],
because sigma fixes the base point and P=-sigma(P). Crucially, E[2]
includes O. A polynomial abscissa forces a nonzero 2-torsion value;
a pole permits O. Two different nonzero values force full cubic splitting
over that branch residue field. A value O does not impose that splitting.
The polynomial exclusions therefore cannot be transferred to all rational
abscissas without treating this extra branch behavior.
This behavior is already used in the retained
[single-branch reciprocity proof, section 3](Q80_SINGLE_BRANCH_RECIPROCITY_2026-09-14.md#3-one-disagreement-produces-the-forbidden-divisor).
The present contribution is its explicit polynomial pole presentation and
height stratification, not discovery of the branch-pole phenomenon.

## A complete pole presentation in the smooth-branch setting

Assume a short nonisotrivial K3 model

```
y^2=x^3+A(t)*x+B(t),  deg A<=8, deg B<=12,
Delta=4*A^3+27*B^2 squarefree of degree24.
```

Infinity is smooth. Let d be squarefree of degree b=2 or4, coprime to
Delta, with its literal rational scalar retained. Infinity is unbranched.
For any finite set of proposed points a rational change of base coordinate
can arrange that infinity is neither a branch value, singular fibre, nor
zero-section intersection; transport the K3 model and d together.

Every nonzero anti-invariant point then has a presentation

```
d=e*f,   e a monic squarefree divisor of d,
x=N/(e*H^2),          y=w*R/(e^2*H^3),
gcd(N,e*H)=1,         gcd(H,f)=1,
N^3+A*e^2*N*H^4+B*e^3*H^6 = f*R^2.                  (1)
```

All data are in Q[t], with H,R nonzero. N may be zero, in which case
coprimality requires e*H to be constant. The scalar of d lies in f.
There is no requirement that H be coprime to e. Write a=deg e, m=deg H.
The bounds for no zero intersection at infinity are

```
deg N <= 4+a+2*m,
deg R <= 6+2*a+3*m-b/2.                              (2)
```

Conversely, (1), (2) and the stated smoothness and coprimality hypotheses
produce an anti-invariant section on this same cover. Two presentations
may have different e and H, but must use the same literal d.

Proof of completeness: x is invariant and y/w belongs to K. At an
unbranched pole the minimal Weierstrass pole orders are 2k,3k, so the
denominator of x has even order. At a branch use a local parameter s
with sigma(s)=-s. The formal parameter z=-x/y is anti-invariant, so its
vanishing order, if positive, is odd, say 2k+1. Then x has pole order
2k+1 on the original t-line and y/w has pole order 3k+2. Put one copy
of that branch factor in e and k copies in H. At an unbranched pole put
k copies in H. Clearing denominators gives (1). A branch not in e has
no pole, hence gcd(H,f)=1. Constants in the denominator can be absorbed
in N and R without changing d. The infinity trivializations are
x_infinity=t^(-4)*x and y_infinity=t^(-6)*y, which give (2).
Direct substitution proves the converse. In particular the often tempting
ansatz x=N/H^2 misses branch poles of odd order on the t-line.

At a root of e, (1) imposes N^3=f*R^2 in the residue field, with N,f,R
nonzero. Equivalently f*N is a square there. This is a condition on the
leading pole coefficients, replacing the finite cubic-root condition at
that branch. Local compatibility remains insufficient for the global
coefficient identity.

## Exact heights and independence

The degree-two pullback has chi=4 and only irreducible fibres: the24 nodal
parent fibres are unramified and the branch fibres are smooth. Counting
zero intersections in (1) gives

```
P.O = a+2*m,       h(P)=8+2*a+4*m.                   (3)
```

An unbranched pole of H contributes twice its multiplicity. At a branch
pole the contribution is 1+2*ord(H). Infinity contributes zero by (2).
Formula (3) now follows from the
[Shioda height formula, section 11.8](https://arxiv.org/html/0907.0298v3#S11.SS8).
It also shows there is no nonzero torsion on this pullback.

For two sections compute h(P1+P2), including all poles, and set

```
c=(h(P1+P2)-h(P1)-h(P2))/2,
det = h(P1)*h(P2)-c^2.                               (4)
```

Positive determinant certifies independence modulo the entire inherited
group. Equal heights and unequal abscissas already suffice: dependence
would force P1=+/-P2 because torsion is absent. Unequal heights require
the determinant or an equivalent exact independence argument. A point
and its multiple must not be counted as two gains.

The conditional construction criterion is therefore explicit: two
solutions of (1)--(2) with common d and positive determinant, on an actual
arithmetic MW17 parent, give rank at least19 on the single quadratic
cover. To obtain infinitely many rational specializations, the conic
must have a rational point, or the genus-one base must have a rational
point and a certified nontorsion Jacobian point. The nonconstant parent
then supplies the usual specialization implication. These are separate
obligations; neither follows from coefficient identities alone.

## What the Q80 closure implies for the next height layer

On the hypotheses above, height8 is exactly a=m=0, the polynomial chart.
The [Q80 genus-one closure](Q80_GENUS_ONE_POLYNOMIAL_CLOSURE_2026-09-14.md)
therefore says that its anti-invariant lattice cannot contain two
independent vectors of height8. If that lattice has rank at least two,
its second successive minimum is at least10. This is conditional on the
retained closure; it does not exclude a lattice with no height8 vectors.

For b=4 the first layers are:

| Height | (a,m) | Necessary pole geometry |
|---:|---|---|
|8|(0,0)|No zero intersection; excluded for a Q80 pair|
|10|(1,0)|One rational branch pole; d must have a rational linear factor|
|12|(0,1) or (2,0)|One rational unbranched pole, or a degree-two branch divisor|
|14|(1,1) or (3,0)|A linear branch divisor plus H of degree one, or a cubic branch divisor|
|16|(0,2), (2,1), or (4,0)|The corresponding allowed branch and pole divisors|

This table lists pole budgets before the inherited arithmetic tests.
The single-branch reciprocity theorem additionally excludes a=3 on Q80
for b=4, at every m: such a point has exactly one nonzero branch value.
It also excludes a=1 for b=2, and requires any two different branch-value
patterns to disagree at least twice. In particular a prospective (8,10)
pair must disagree at a second branch besides the height10 point's pole.
These are retained exclusions, not new computations.

All divisors and degrees here are on the original projective line. In
particular, an irreducible branch quartic has a=0 or4. Below height16
only heights8 and12 are possible, and height12 requires a rational
unbranched pole. Its second successive minimum on Q80 is consequently
at least12. These are necessary conditions, not existence claims.

A first genuinely new Q80 coefficient gate is (a,m)=(1,0), with
deg N<=5, deg R<=6, common quartic d=e*f and deg f=3. A pair may have
heights (8,10) or (10,10); both merit distinct independence checks.
For an irreducible quartic this entire height10 layer is empty before
any calculation. Branch divisors meeting Delta require a separate
analysis with fibre corrections and are outside this note.

## An exact control for the branch-pole term

Reuse the existing genus-one control, which does not have an MW17 parent:

```
d=t^4+t+1, s=t^4+1, A=d*(2*s+1)-1, B=d*s^2,
P=(0,s*w), h(P)=8.
```

Its double has the presentation

```
e=d, f=1, H=2*s, N=A^2,
R=-8*d^2*s^4-A^3.
```

The addition formula gives x(2P)=A^2/(4*d*s^2), and substitution yields
exactly (1). Here (a,m)=(4,4), so (3) gives32, agreeing with h(2P)=4*h(P).
Thus odd branch poles occur even on a retained example. This control
is a multiple, not a new independent direction. The control parent,
smoothness and original height are inherited from the
[common-quartic control note](COMMON_QUARTIC_SINGULARITY_LOCUS_2026-09-14.md).

Reproduce the elementary control identity with the repository's
`.venv/bin/python` (checked with SymPy 1.14.0). The system `python3`
initially failed to import SymPy; using the existing repository environment
required no dependency installation. All four assertions below passed:

```python
import sympy as S
t = S.symbols('t')
d, s = t**4+t+1, t**4+1
A, B = d*(2*s+1)-1, d*s**2
e, H, N, R = d, 2*s, A**2, -8*d**2*s**4-A**3
assert S.expand(N**3+A*e**2*N*H**4+B*e**3*H**6-R**2) == 0
assert S.gcd(N,e*H) == 1
assert (S.degree(N,t),S.degree(R,t)) == (16,24)
assert 8+2*S.degree(e,t)+4*S.degree(H,t) == 32
```

No height10 or height12 point on an MW17 parent was constructed. The
advance is a complete pole presentation in this smooth-branch setting
and a precise next height gate. A parent-side geometric condition that
forces two global solutions of (1) with an infinite rational base remains
unknown; calling the twist-rank reformulation itself that mechanism would
overstate the result. No coefficient search or new exclusion campaign ran.
