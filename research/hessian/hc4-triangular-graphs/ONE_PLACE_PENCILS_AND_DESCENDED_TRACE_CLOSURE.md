# One-place curve reconstruction closes every descended rank-one trace

15 November 2026. Claims **HC4MYIZ0/1/2**.
This closes the multibranch descended-trace cases left open in
[IY2](CONTRACTIBLE_RADICAL_CURVES_AND_PENCIL_DEFECT.md).
All helper degrees remain unrestricted. The unrestricted HC4/JC2
problem is not resolved.

## IZ0: an intrinsic formula for the generic genus of a one-place pencil

Let b in C[x,y] be irreducible and nonconstant. Let A=V(b), and assume
the smooth projective completion of its affine normalization adds exactly
one point. This is the one-place-at-infinity condition. Let g0 be the
genus of that smooth projective normalization, and define

    delta(A)=dim_C H0(A, normalization_* O_normalization/O_A),

the total finite normalization deficit. Then the geometric generic
fibre of b has exactly one puncture and genus

    G(b)=g0+delta(A).                                      (1)

In particular these two generic-pencil invariants are the same for
*any two isomorphic such affine plane curves*, even if their embeddings
are not equivalent under an ambient automorphism.

Here are the exact newly imported external foundations.

1. [Vistoli, arXiv:math/9904035v3, Theorem 1 and Corollary 1](https://arxiv.org/html/math/9904035v3):
   for a primitive complex plane polynomial with one fibre having one
   branch at infinity, its proper minimal normal compactification over
   A1 has integral fibres, a single boundary section, and fibres smooth
   at that section. The introduction states the equivalent simultaneous
   resolution-at-infinity consequence.
2. [Stacks 0B9T](https://stacks.math.columbia.edu/tag/0B9T):
   coherent Euler characteristic is locally constant for a finitely
   presented flat sheaf in a proper finitely presented family, and
   commutes with base change.
3. [Stacks 0C1R](https://stacks.math.columbia.edu/tag/0C1R):
   normalization of a reduced algebraic curve is finite, with quotient
   supported at its singular points and lengths equal to its local
   delta invariants.

These are external inputs, not symbolic-checker conclusions.
The deductions, including the hypotheses of the compactification and
flatness statements, are as follows.

An irreducible b is noncomposite: a composition b=P(h) with deg P>1
factors over C into nonconstant h-root factors. The retained IL0
noncomposite argument therefore makes its geometric generic fibre
integral. Thus Vistoli's primitivity hypothesis holds.

The one-place condition is intrinsic to the affine curve. In any plane
embedding the coordinates extend to the smooth projective
normalization. Their poles can occur only at the one removed point;
at least one has a pole there, since otherwise all are constant on
the projective curve. Hence its projective plane closure has just
one normalized point above the line at infinity.

Use the proper normal integral surface X over A1 supplied by Vistoli.
It is of finite presentation. It is flat over A1: on an affine chart
its ring is a torsion-free C[t]-module, because the surface is integral
and dominates A1. Every finitely generated submodule is free over that
PID, and the module is their filtered union, hence flat. This verifies
flatness rather than assuming it from properness or smooth general fibres.

Let Y0 be the proper zero fibre on X. It is integral, agrees with A on
the affine part, and adds one smooth point. Its normalization is the
smooth projective normalization of A; the finite singularity deficit
is exactly delta(A), with no extra term at infinity. The normalization
sequence gives

    chi(O_Y0)=1-g0-delta(A).                               (2)

By proper flat Euler-characteristic constancy, (2) is also the
Euler characteristic of the geometric generic proper fibre. That
fibre is smooth: the affine generic fibre is smooth in characteristic
zero and Vistoli supplies smoothness at its unique boundary point.
It has Euler characteristic 1-G(b). This proves (1) and the one-puncture
assertion.

If A and C are isomorphic, their normalizations and smooth projective
completions are isomorphic, and the finite quotient in the normalization
sequence has the same length. Thus both g0 and delta agree.
No ambient extension of the curve isomorphism, plane degree equality,
Puiseux-degree bound or finite computation is required.

### Polynomial-normalization cross-check and the genus distinction

When the affine normalization is A1, g0=0, so G(b)=delta(A).
The singular zero fibre can have normalization genus zero even though
a smooth pencil member has positive genus. A node and a cusp both
have delta=1 and generic genus one, despite their different topology.

For A=Spec R with R contained in its normalization C[t], the semigroup
of degrees of nonzero elements of R is also intrinsic: a curve
isomorphism extends to t->a t+c on the normalization, preserving degrees.
This reproduces the mechanism of
[Assi-Garcia-Sanchez, arXiv:1407.0490v1, Proposition 19](https://arxiv.org/html/1407.0490v1).
Its Section 5 final genus paragraph gives a consistent one-place
semigroup-genus description. The proof above uses the compactification
and normalization sequence, and does not import an algorithm, degree
search, equality of embedded degrees, or an ambient equivalence theorem.

## IZ1: reconstruction over the complete inverse image of one one-place curve

Let K:A2_z->A2_T be complex plane Keller. Let A be an irreducible
closed target curve with one place at infinity. Suppose the *entire
reduced* inverse image C=K^(-1)(A) maps isomorphically onto A.
Then K is a polynomial automorphism.

Let b define A. The isomorphism makes C irreducible, with one place
at infinity as well. Let g define C. HW0 supplies square-free Keller
pullbacks. Complete prime support therefore gives

    b(K)=c g,       c in C*.                               (3)

There is no omitted source component or extra multiplicity in (3).
Both polynomials are irreducible and hence have geometrically integral
generic fibres by the retained IL0 argument.

IZ0 gives those two generic fibres equal genus G and exactly one
puncture. Apply the exact IY0 boundary identity to (3):

    2G-1=d(2G-1)+E,

where d is the Keller field degree. The nonconstant finite morphism of
smooth projective completions has a point over the target puncture.
It must be the only source puncture, so E=0. Since 2G-1 is a nonzero
odd integer, d=1. The retained coprime-pullback denominator argument
then gives the actual polynomial Keller inverse.

This extends IY1 beyond contractible curves. Internal nodes, multiple
cusps, other finite singularities and positive normalization genus
are all allowed, provided the full curve isomorphism and one-place
condition hold. The source and target plane embeddings can differ.

The proof does not say that an affine etale curve map is automatically
finite. The unique-puncture hypotheses, and the nonconstant projective
extension, prove that no source puncture has finite image in this
particular argument. The nonfinite Gm->A1 controls retained in IY
remain valid outside it.

## IZ2: the entire target-descended rank-one trace family closes

Let K=(K1,K2) be complex plane Keller. Suppose an actual closed embedding

    j(z,w)=(w1,w2,X(z,w),Y1(z,w),Y2(z,w))

satisfies

    X^2+w1Y1+w2Y2=D(z)+w dot K(z).                        (4)

Put f=X(z,0), so D=f^2. If for even one positive integer m

    D^m belongs to C(K1,K2),                              (5)

then K, and the original cotangent gradient, are polynomial
automorphisms. This includes every target-polynomial or target-rational
trace and every positive-power descent. All helper degrees and
higher free jets are unrestricted.

If f is zero or constant, the earlier critical/noncritical constant-trace
theorems apply. Suppose f is nonconstant. Since f^(2m) belongs to the
target field k, k(f)/k is finite cyclic Galois: all roots of its separable
power equation are root-of-unity multiples of f and lie in k(f).
HW0 excludes a nontrivial abelian Galois intermediate field. Hence
f belongs to k. The retained denominator argument R intersect k=C[K]
then gives

    f=a(K),       D=a(K)^2                                (6)

for a nonconstant target polynomial a.

IX1 uses the actual closed zero-section and first free jet to make
K a closed immersion on V(f). IX2 proves that it identifies the
*entire reduced* curve a(K)=0 with the *entire reduced* target curve
a=0. Thus each irreducible target component has its whole inverse
image isomorphic to it.

Suppose K were nonautomorphic. By the retained IX2/Jelonek-Lason
argument, every component of a=0 has affine normalization A1:
singleton fibres on that target curve place it in the nonproperness
set, whose components are polynomially parametric. Such a component
has exactly one place at infinity. Apply IZ1 to any one component.
It forces K invertible, a contradiction.

This closes every descended case in (5), including all internally
multibranch components previously left open in IY2. There is no
residual nodal or other polynomially parametric zero-curve exception
within this family.

For the original potential h=D(z)+w dot K(z), its gradient inverse is
explicit once K is inverted: its free outputs recover z through K,
and the other outputs recover w from
eta=grad D(z)+DK(z)^T w. No nonlinear source Hessian invariance is used.

### A componentwise criterion even without global trace descent

For any irreducible component C of the reduced f=0 curve, the actual
radical-curve closed immersion identifies C with its closed image A.
If C has one place at infinity and is the entire reduced K inverse
image of A, IZ1 forces the inverse. It suffices that the irreducible
equation b of A have just one source prime in b(K), namely C;
square-freeness is automatic by HW0.

Thus under hypothetical nonautomorphy a one-place radical component
must share its target curve with an additional source prime outside
the radical curve. It cannot share it with a second component of f=0,
because K is injective on that entire curve. Source components outside
f=0 are not silently discarded.

### Constant affine frames and polynomial offsets

The same conclusion holds when a centred ambient linear free pair pulls
back to A w+v(z), A in GL2(C), and the restricted helper quadratic has
rank one. In the actual source coordinates (z,b) the potential is

    D'(z)+b dot K'(z),
    D'=D-v^T A^(-T)K,       K'=A^(-T)K.

If any positive power of D' belongs to C(K')=C(K), the rank-one theorem
applies. This is a literal cotangent identity and a retained actual
closed embedding; it is not a claim of Hessian invariance under an
arbitrary nonlinear source change.

## Assurance and the remaining frontier

See [source context](sources/one-place-trace-context.json) for the
newly activated compactification input, the two explicit foundational
Stacks inputs and all retained byte pins. They are not independently
reproved or formally verified. The intrinsic genus deduction, exact
full-preimage argument and HC4 consequences are written above.

Run

    .venv/bin/python research/hessian/hc4-triangular-graphs/verify_one_place_trace.py

for [exact controls](receipt-one-place-trace.json), and see
[scoped validation](validation-one-place-trace.json).
Checks distinguish the zero-fibre normalization genus from the generic
genus, verify normalization-ring identities and genuine nonlinear
nodal-trace embeddings, and retain full-preimage and boundary-loss
hypotheses. They are not a finite proof of the new external results.

A hypothetical remaining rank-one retained-free embedding therefore
has a genuinely non-descended normal trace: no positive power of f
or D belongs to the target field. Any one-place radical component
must be unsaturated as described above. Higher helper ranks, unaligned
intrinsic projections, actual endpoint representatives and unrestricted
HC4/JC2 remain open. Earlier IX/IY descended multibranch frontiers
are historical; they are closed by IZ2, with their source bytes retained.
