# The exceptional divisor of the all-rank affine completion

This computes an actual geometric feature of KMFI104's Keller map.
For rank>=2 its exceptional target hyperplane is fully in the image:
a general point has exactly two affine preimages and a codimension-two
subspace has exactly one. In the remaining PSL2(11) rank-four case,
the normalization has eleven split sheets over that hyperplane.
Nine are omitted unramified divisors.

These statements concern the specified higher-weight representative
with weights (1,r+1,...,2r) and generic branch centering.
Finite normalization depends on the target affine ring.
The numerical bounds here are not absolute bounds for every Keller
realization or every birational target model.

## 1. The two actual source components over C0

Let S=X(p)+gamma(1,p), det DS=gamma^r, rank DX=r, r>=2, and suppose
its stabilizer normally generates its group. Choose the center p*
where D=X(A^r) is smooth, X is etale onto D, h in
d(S)=gamma^(r+1)h is nonzero on gamma0, and no other branch component
of the normal closure passes through X(p*).
The branch degree may be greater than one: X is a local
isomorphism at the chosen point.

Center p*,X(p*) at0 and choose a regular Jordan frame for
P(0)=[-p|I]DX(0). It is regular nilpotent: DX has rank r,
its image contains Delta, so P has rank r-1.
Write J for the upper Jordan block with superdiagonal ones.
Use KMFI104 with a_i=r+i-1, e0=1, e_i=r+i and its finite rectifier.
Denote the completed target by (C,Z0,...,Zr), so
\[
       Y_i=C^{e_i}Z_i,\qquad C=x\gamma,\qquad
                         Z_i=E_i(\Theta(v))/x^{e_i}.        \tag{1}
\]

**KMFI106.** The actual source preimage of C0 is the disjoint union
\[
       \{x=0\}\sqcup\{\gamma=0\}
                    \cong A^{r+1}\sqcup(G_m\times A^r).    \tag{2}
\]
The first component maps isomorphically onto all C0.
The second maps isomorphically onto C0 minus {Z0=0}.
Thus C0 has two affine preimages over Z0!=0 and one over Z0=0.

Here is the polynomial proof, including its possible nonlinear term.
After p_i=gamma^(a_i)u_i, the stage at gamma0 is
\[
             E(0,u)=(1,u1+u2,\ldots,u_{r-1}+ur,ur).         \tag{3}
\]
The only possible nonlinear contribution is a multiple of u1^2
in the last component. Differentiate tr(P^r)=0 in p1 at0:
since J^(r-1) has its only nonzero entry at(1,r),
\[
       \partial_{p1}P_{r1}(0)=0
                          =\partial_{p1}^2X_{r+1}(0).
\]
The derivative of -p_r*partial_(p1)X1 vanishes because r>=2.
Thus this quadratic is zero. All other nonlinear source terms
have weight>2r; the linear and gamma p terms give (3).

The finite-jet rectifier gives on x0
\[
 Z0=y0,\quad Z_i=y_i\ (1<=i<r),\quad
                              Zr=y_r+\beta y0^{2r}.        \tag{4}
\]
Its pure-t errors start at degree2r, its lower-index linear
errors start there too, and terms with two high variables exceed
every required weight. The first output has linear term t.
Formula(4) is an automorphism.

The polynomial gamma on x0 is the nonzero constant gamma0,
so (x,gamma) is the unit ideal. On gamma0, x is invertible,
and the inverse of Theta identifies this component with
(x,u) in G_m times A^r. Formula(3) gives
\[
 Z0=1/x,\qquad Z_{\rm high}=
             \operatorname{diag}(x^{-e_1},\ldots,x^{-e_r})(I+J)u .
\]
The inverse is x=1/Z0 and
u=(I+J)^-1 diag(x^(e_i)) Z_high.
This proves (2) on the whole hyperplane, not just generically.
The rank-deficient S3 control in KMFI104 is not used to
justify the regular-Jordan hypothesis.

## 2. Full normalization inertia at the new hyperplane

Let sigma generate the full normal-closure inertia over D,
let m=ord(sigma), and put e=r+1. The visible cycle gives e|m.

**KMFI107.** At the generic point of the new hyperplane C0,
the full decomposition and inertia groups are
\[
                   D_C=I_C=\langle\sigma^e\rangle.          \tag{5}
\]
Every normalization prime there has residue degree1.
Its ramification indices are the cycle lengths of sigma^e.
If sigma^e=1, the whole normal closure splits there.

At the chosen smooth branch center, a tame finite cover with
one regular branch divisor is etale-locally Kummer
([Stacks, Tag0EYG](https://stacks.math.columbia.edu/tag/0EYG)).
Use the actual finite normal closure, without assuming S finite.
A connected local Galois component is t^m=d(Y), up to a unit
with constant residue. The etale coefficient extensions split
at the algebraically closed center and lift to the henselian
C-adic ring k(Z0,...,Zr)[[C]].

The lowest weighted part of d(Y), for weights e_i, is
\[
                             aY_0^e,\qquad a\ne0.          \tag{6}
\]
If its weight M were less than e, substituting
S_i=gamma^(e_i)E_i into d(S)=gamma^e h would make that initial
part vanish on E(0,u)=(1,(I+J)u), impossible:
restriction to Y0=1 is injective on a fixed weighted degree.
If M>e the gamma order would exceed e, since h(0,0)!=0.
Thus M=e. Its only monomials are Y0^e and Y1.
The Y1 coefficient is zero: the normal to D at0 is the last
coordinate normal, and r>=2. This proves (6).

After (1),
\[
                 d(C^{e_i}Z_i)=C^e(aZ0^e+C R).
\]
In the C-adic henselian field, the final unit is aZ0^e times
a principal unit with an mth root. Since a has all roots in k,
the equation reduces to t^m=(CZ0)^e.
It factors into e extensions t^(m/e)=zeta CZ0, each totally
ramified with unchanged residue field.
Its image in the original inertia is generated by sigma^e.
This proves (5), including the residue assertion.

No residue S3 or extra branch inertia at C0 is inferred from
missing sheets: the regular target substitution pulls back the
cyclic local cover at the original smooth branch center.

The number of omitted unramified primes over C0 is exactly
|Fix(sigma^e)|-2, by KMFI106. It is at least r: the visible
e-cycle becomes e fixed sheets and the actual unramified tail
gives another fixed sheet; only two are retained.

## 3. Exact class-group growth

Let Z_S be the finite normalization for the original sweep,
and Z_F that for the constructed Keller map.
Let q be the number of cycles of sigma^e, including fixed points.

**KMFI108.** Both class groups are free of finite rank, and
\[
 \operatorname{rank}Cl(Z_F)=\operatorname{rank}Cl(Z_S)+q-1.  \tag{7}
\]
The basis of Cl(Z_S) consists of its prime divisors omitted
from the polynomial sweep source.

The birational morphism A^(r+1)->Z_S exists because integral
elements are polynomials on the normal source.
It contracts no source divisor: outside gamma0, S is etale;
gamma0 dominates D. Its quasi-finite locus contains every
source codimension-one point. Zariski Main identifies this
locus with an open of Z_S; its source complement has
codimension>=2, hence its class group is0 and its units are k*.
The divisor restriction sequence makes Cl(Z_S) freely generated
by its omitted prime divisors, using the principal-divisor
argument in KMFI87. No open immersion of the entire
non-quasi-finite source is asserted.

Over C!=0, (1) gives (Z_F)_C=Z_S times G_m.
The polynomial/Laurent normalization is normal and finite
with the correct fraction field, proving this equality.
Its class group is Cl(Z_S): horizontal prime divisors become
principal over k(Z_S)[C,C^-1], leaving vertical classes;
a rational function with vertical divisor only is a*C^j,
so vertical relations are exactly those on Z_S.

Units of Z_S are k*, since they pull back to polynomial source
units. Thus units of (Z_F)_C are precisely k* C^Z.
Removing the q primes over C0 gives
\[
 0\longrightarrow
 \mathbf Z^q/\mathbf Z(e'_1,\ldots,e'_q)
 \longrightarrow Cl(Z_F)\longrightarrow Cl(Z_S)\longrightarrow0,
\]
where e'_j are the cycle lengths of sigma^e.
At least two e'_j are1, so the relation is primitive.
The kernel is free of rank q-1, proving (7).
Equivalently, the actual boundary over C!=0 has the original
omitted divisors and the now-omitted visible gamma0 divisor;
over C0 it has q-2 more primes. These are distinct divisors.

## 4. The remaining PSL2(11) field

For a hypothetical KMFI105 input, r4, e5 and sigma has type5^2*1.
Thus sigma^5=1 and KMFI107 gives:

- C0 is unramified, with eleven residue-degree-one normalization primes;
- exactly two are retained and nine are omitted unramified;
- Cl(Z_F) has rank Cl(Z_S)+10, hence at least11, because Z_S
  already omits one index5 prime over D.

The two order5 primes are both omitted from the actual Keller
source over C!=0 and have canonical coefficients4.
The nine new unramified primes have coefficient zero.
Additional boundary elsewhere is allowed.
This does not assert that Z_S itself has class rank11.

All of C0 is in the image: fibers have size2 off the linear
subspace H={C=Z0=0}, and size1 on H.
Off H the strict local group is trivial, but only two of its
eleven fixed sheets are retained.
On H the local group is C5: the regular substitution maps to
the original smooth branch center with cyclic group C5;
the transformed main branch passes through H and keeps its
nontrivial inertia. Its one fixed sheet is the actual x0 section.
Thus KMFI94's bound is strict off H because unramified omissions
are mandatory. None of these constraints is an action exclusion.

## 5. A nonalternating positive control

Apply the higher-weight algorithm to KMFI60's finite binary
degree8 sweep, centered at (1,0), using a=(2,3), e=(1,3,4).
The marked first-stage point is
\[
              \theta0=(4s/3,3(1-s)/16,3/8),\qquad s^2=-3.
\]
The retained model specifies E, the inverse Jacobian, the curve
jet and unipotent matrix: there are no unresolved coefficients.

Its original full inertia is3^2*1^2, whose cube is the identity.
All eight primes over C0 split; exactly two are retained and
six are omitted unramified.
The original finite sweep has ordinary A3 source, so Cl(Z_S)=0.
Equation(7) proves **Cl(Z_F)=Z^7**.
There is one ramified normalization prime outside C0, the
original visible index3 prime. The canonical class is twice
one primitive boundary basis element; its other six coefficients
are zero. The action is the already-established binary eight,
not PSL2(11).

Only the ordinary degree4 jet is needed for this rectifier:
higher ordinary terms remain above every output weight after
the triangular substitution, which preserves the maximal ideal.
The incomplete unsimplified symbolic probe is retained as
exploration, with no result inferred from it.

    .venv/bin/python research/jacobian/monodromy-theory/verify_completion_boundary.py

[Helper](completion_boundary.py) · [model](completion_boundary_model.json) ·
[receipt](completion_boundary_receipt.json) ·
[primary pins](sources/completion-boundary-manifest.json).

The fiber, Kummer and class-group statements have written universal
proofs. Symbolic identities, complete group profiles and the fully
specified positive input are supplemental controls.
Formal verification, independent whole-proof replay and external
review are not claimed. Exceptional PSL2(11) affine realization
and the complete action spectrum remain open.
