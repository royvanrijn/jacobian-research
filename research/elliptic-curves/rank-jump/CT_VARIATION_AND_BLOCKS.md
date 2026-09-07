# Simultaneous Cassels–Tate blocks, 5 September 2026

The fixed-cubic control now has an exact **simultaneous 10+2 decomposition of
two restricted Cassels–Tate forms**. The ten-dimensional component cannot be
split into smaller common orthogonal nondegenerate components. This is stronger
than separately putting each alternating matrix into symplectic normal form.
It is a block of changing **solubility obstructions**, not a construction of ten
rational directions. No prospective rank selector follows yet.

The [certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_ct_variation_v1.json)
and [portable formula inputs](../../artifacts/generated-results/elliptic-curves/rank_jump_ct_variation_inputs_v1.json)
are replayed by `python3 elliptic-curves/rank-jump/ct_variation.py check`.
The [protocol](CT_VARIATION_PROTOCOL.json) explicitly records retrospective
selection. No new descent, curve search, or parameter sweep was run.

## Common classes, different obstructions

Use the labelled cubic algebra, anchor subgroup and seven retained local
subspaces from [the preceding note](LOCAL_COLLISION_AND_RECIPROCITY.md).
The anchor is the historic Fermigier rank-at-least-20 curve. Its twenty
independent rational Kummer classes span (W). In the deformation

\[
\alpha=\theta+u\theta^2,
\]

the generic arithmetic rank is zero; these twenty anchor classes are **not**
twenty generic sections of the deformation. (W_u\subset W) is the complete
locally admissible intersection inside this fixed twenty-dimensional space.
The retained CT matrix (B_u) is complete on (W_u), not asserted complete on
the entire Selmer group.

For each of the 21 pairs set (C=W_u\cap W_v) and
\(D_{u,v}=B_u|_C+B_v|_C\), with values in \(\mathbf F_2\).
The common cubic-algebra labels identify the *classes*, making this difference
well-defined. Comparing the untransported matrix arrays would not do so.

Across the fifteen pairs of distinct nonzero deformations, the rank of (D)
is 6 through 12. Thus matching 2-torsion representations and restricting to
classes that pass both sets of local conditions does not preserve their CT
pairings. Three examples are:

| (u,v) | \(\dim C\) | \(\operatorname{rk}B_u|_C\) | \(\operatorname{rk}B_v|_C\) | \(\operatorname{rk}D\) |
|---|---:|---:|---:|---:|
| (-3,-1) | 15 | 12 | 14 | 12 |
| (-2,1) | 7 | 4 | 6 | 6 |
| (-1,1) | 12 | 10 | 10 | 12 |

Two precise limitations on explanatory models follow by linear algebra.

1. If a proposed **linear class feature map** (L:C\to\mathbf F_2^k) explains
   (D(x,y)=b(Lx,Ly)), then \(\ker L\subset\operatorname{rad}D\) and
   (k\ge\operatorname{rk}D). This rejects only this specified factorization
   model. Three generators of a varying number field are not automatically
   three such class features; higher cochain constructions remain possible.
2. An expression of (D) as a sum of elementary alternating rank-two forms
   needs at least \(\operatorname{rk}D/2\) summands. The already computed local
   admissibility functionals vanish on (C), so bilinear expressions made
   solely from those functionals cannot give these nonzero differences.

These are **solubility-obstruction** results. They are not visibility results,
nor rank bounds on the whole curves.

## A ten-dimensional component that cannot be separated

The pair ((-1,1)) is the unique nonzero pair in this panel for which (D) is
nondegenerate on all of (C). Put (A=B_{-1}|_C), (B=B_1|_C), (D=A+B).
Define (T) by (A(x,y)=D(x,Ty)). It is self-adjoint for (D), and

\[
\mu_T(t)=t^5(t+1),\qquad
\chi_T(t)=t^{10}(t+1)^2.
\]

The primary decomposition is

\[
C=N\perp U,\qquad N=\ker T^5,\quad U=\ker(T+I),\qquad
\dim N=10,\quad\dim U=2.
\]

It is orthogonal for **both** CT forms and (D). On (U), (T=I), (A=D)
and (B=0). On (N), (T) consists of two length-five nilpotent chains.
The dimensions of \(\ker(T^j|_N)\), for (j=0,\ldots,5), are
\(0,2,4,6,8,10\). Hence (A|_N) has rank 8 and (B|_N) rank 10.

Here are explicit chains, as packed masks in the original twenty anchor
classes, with index starting at zero:

\[
\begin{aligned}
(e_0,\ldots,e_4)&=(593,201198,186130,317529,631775),\\
(f_0,\ldots,f_4)&=(491700,992370,111568,244478,206147).
\end{aligned}
\]

The replay proves (Te_i=e_{i+1}), (Tf_i=f_{i+1}), with the last vector
mapping to zero. It also verifies the full pairing table:

\[
D(e_i,f_j)=[i+j=4],\qquad A(e_i,f_j)=[i+j=3],\qquad B=A+D.
\]

Each chain pairs trivially with itself for all three forms. The two masks

\[
154245,\quad845062
\]

span (U). These are representatives of cohomology classes; no assertion is
made that their anchor point representatives specialize to rational points
on either deformed curve.

**Indecomposability proof.** Any common orthogonal decomposition into
nondegenerate subspaces for (D) is (T)-stable. If (N=V\perp V') were
such a proper decomposition, at least one summand would contain a length-five
cyclic chain. Every cyclic (T)-subspace is (D)-isotropic: for even (k),

\[
D(x,T^kx)=D(T^{k/2}x,T^{k/2}x)=0,
\]

and for odd (k=2j+1), this equals (A(T^jx,T^jx)=0).
Self-adjointness then makes all pairs within the cyclic chain vanish. A
nondegenerate alternating space containing an isotropic five-dimensional
space has dimension at least ten. That summand must therefore be all of
(N), a contradiction. The component is indecomposable in this precise
simultaneous orthogonal sense.

This is **not** a decomposition of the full arithmetic Selmer object, a
canonical grouping of Mordell–Weil points, or evidence that a low-degree curve
generates the chains. It is an intrinsic component of this restricted pair of
bilinear forms under the retained class identification.

The common-space radicals are especially easy to overinterpret here.
At (u=-1), \(\operatorname{rad}A\) is spanned by (e_4,f_4), but neither
nonzero direction in their span survives pairing with the full (W_{-1}):
the preceding rectangular audit gives \(C\cap\operatorname{rad}B_{-1}=0\).
At (u=1), the common-space radical is (U), but only one of its two
dimensions survives pairing with full (W_1). Even that surviving class is
only a necessary candidate for rational solubility.

There is also an obstruction to the simplest arithmetic interpretation of
the chains. The shared (E[2]) has full (S_3=\mathrm{GL}_2(\mathbf F_2))
action, whose commuting endomorphisms are only 0 and (I). Consequently (T)
cannot be the restriction of a map on cohomology induced by an endomorphism
of this Galois module: those maps are zero or identity, whereas (T) has
minimal polynomial (t^5(t+1)). Likewise the non-Galois cubic field has no
nonidentity rational field automorphism. A possible arithmetic explanation
must use more data than a fixed automorphism of the shared descent algebra
or torsion module. This does not exclude higher correspondences or
pair-dependent norm constructions.

## Individual negative local terms have no intrinsic support

Fisher's pairing formula uses products of Hilbert symbols

\[
\prod_p(a,\gamma(x_p))_p,
\]

with (a) the leading coefficient of the second quartic. Its theorem and
finite-support argument are in [Fisher, Theorem 3.1 and Remark 3.3](https://antsmath.org/ANTSXV/papers/ANTS-XV_fisher.pdf).
Replacing \(\gamma\) by (c\gamma), (c\in\mathbf Q^\times), changes each
term by \((a,c)_p\), but changes the global product by 1, by rational Hilbert
reciprocity. This is a change of representative of the formula, not a change
of the curve or its CT pairing.

For the first three retained pairs in each of six uniformly formatted
transcripts, the audit tried (c=-1,2,3,5). It recomputed the original exact
quartic and gamma evaluations and all local Hilbert symbols. The retained
support contains (2,3,5,\infty) and all numerator/denominator primes of
(a); outside it every change symbol is 1. No new factorization was needed.

Of 72 rescalings, **57 change the set of negative local terms; none changes
the global pairing**. At (u=-2), for masks (111,282), the original unique
negative term is at (7{,}819{,}109). Scaling gamma by 2 makes the unique
negative term occur at (301{,}565{,}794{,}147) instead. The pairing remains 1.
At the anchor, all twelve tested rescalings leave the local terms unchanged
because the selected leading coefficients are rational squares.

A fixed normalization can still be used as a computational convention.
Counts and locations of its negative terms need a proof of invariance before
being interpreted as arithmetic support of a soluble or obstructed block.
This does not invalidate intrinsic local Kummer conditions or a correctly
normalized global governing-field construction.

## Mechanisms and the next mathematical gap

The current ranking is:

1. **Higher descent with simultaneous cancellation of CT obstructions**
   (solubility): strongest remaining structural mechanism. The present
   certificate shows a large inseparable component of *variation*, while the
   anchor has twenty rationally soluble classes. It does not explain the
   condition that makes a large component soluble.
2. **Correlated local incidence at root collisions** (incidence): proved in
   the preceding note, including the reciprocal relation among new-prime
   cuts. It explains which inherited classes remain admissible and stops
   before the CT obstruction.
3. **A common auxiliary curve or cover** (incidence and solubility if proved):
   remains plausible, but the frozen small bisection dictionary failed to
   supply a shared quadratic carrier. No new positive construction here.

Weak or rejected explanations are preservation of the cubic field alone,
bilinear use of local cuts that vanish on the common space, individual
negative-term prime counts, and chart visibility treated as a rank predictor.

A relevant positive theorem exists for **quadratic twists with full rational
2-torsion**: [Smith, Theorem 3.2 and §3.1](https://arxiv.org/pdf/1607.07860)
expresses appropriate CT differences through splitting in pair-specific
governing extensions constructed using norm equations. Our (S_3) pencil is
not that twist family. Passing to its splitting field supplies full
2-torsion but does not make the varying curves quadratic twists. The theorem
therefore cannot simply be applied to our 10-dimensional component.

The next bounded experiment should seek an **arithmetic lift of one certified
chain relation**, rather than fit an arbitrary matrix. Fix the masks of
(e_3,f_0,e_4) above before any further computation. At (u=-1,1), their
common-space pairings satisfy

\[
A(e_3,f_0)=B(e_3,f_0)=1,\qquad
A(e_4,f_0)=0,\quad B(e_4,f_0)=1.
\]

The useful target is to derive these two entries from the fixed cubic
classes, the factors (1-u\theta_i), and explicit norm/cochain data with
reciprocity removing normalization dependence. Four retained-entry
reconstructions, no parameter enumeration and no point search, are sufficient
for the first test. Merely recomputing Fisher values passes an implementation
check but **fails the explanatory endpoint**. A formula must explain why the
same operator (T) links the masks; otherwise the pencil remains a linear
algebra observation without its missing arithmetic implication.

Agent 1 can eventually use a cheap, point-blind version of such an invariant
as a **solubility filter**, combined with independent incidence evidence.
There is no justified new selection feature to hand over now. The missing
chain remains

\[
\text{point-blind specialization condition}
\Longrightarrow\text{several independent globally soluble classes}
\Longrightarrow\text{large Mordell–Weil jump}.
\]

Eliminating CT obstructions would still leave the possible higher-divisible
Sha obstruction. Rational witnesses or a theorem closing that implication
are required.
