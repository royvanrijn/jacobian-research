# One exceptional fibre closes both positive-p base charts

14 November 2026. Let `Psi=L*A^2+M*A+N*B`, with `L*N != 0`, be the
scaled Meng--Yang potential. Consider either genuine polynomial graph

\[
 y=R(x,p)\quad\hbox{or}\quad x=R(y,p),\qquad \deg_p R>0,
\]

and then any actual polynomial automorphism of the three retained duals
over the free-base polynomial ring. Degrees and word lengths are arbitrary.

**HC4MYGZ2. In both charts, a nonzero constant four-variable Hessian
determinant forces a polynomial gradient inverse.** The opposite chart
is now completely closed, including every finite-base escape and every
number of distinct critical roots. The earlier six-root frontier is retired.
If the boundary has no simple critical line, unit Hessian is impossible.
The entirely repeated/simple-boundary case keeps the distinction between
invertibility and nonexistence of unit Hessian.

The proof compares the **same exceptional base-zero value fibre** in
the native graph and in the necessary Tang form. It needs no classification
of escapes, Newton faces, finite critical projections or covering degrees.
It also gives a shorter independent proof of the preceding positive-p
y-chart theorem. It retains Tang and the already imported boundary-line
theorem as inputs; it does not resolve the two plane Keller endpoints.

## 1. A general exceptional-fibre obstruction (HC4MYGZ0)

Work over C, with compactly supported Euler characteristic. For a
polynomial `f(s,z1,z2,z3)`, let

\[
 C_f=V(f_{z_1},f_{z_2},f_{z_3}).
\]

Suppose the following three conditions hold, for a fixed scalar s0:

1. For generic s and generic lambda, `chi(f(s,-)=lambda) != 1`.
2. Every irreducible component of C_f projects nonconstantly to s,
   and C_f has no point over s0. The empty scheme is allowed.
3. For generic lambda, `chi(f(s0,-)=lambda) <= 1`.

Then **f cannot have nonzero constant Hessian determinant**. These
conditions and conclusion are invariant under actual polynomial retained-dual
automorphisms over C[s]. Relative critical schemes and whole value fibres
are preserved by [HC4MYGK0](RELATIVE_CRITICAL_LOCUS.md); no nonlinear
invariance of the Hessian away from that scheme is assumed.

Here is the proof. Translate s0 to zero. Under unit Hessian, C_f is smooth
of pure dimension one when nonempty, since its three Jacobian rows are
rows of the invertible Hessian. If grad f were an automorphism, C_f would
be the inverse image of a target coordinate line, hence an affine line.
Its base projection would be a nonconstant polynomial on A1 avoiding
zero, impossible over C. An empty C_f is also impossible for an automorphism.

For a nonautomorphic gradient, the retained [HCTANG2 structural
alternative](../tang-hc4-audit/STRUCTURAL_REVIEW.md) supplies a constant
jointly affine source two-plane. Its projection to s has rank zero or one.
The following adaptations use only constant affine changes preserving
the free-base fibration, so they preserve the fibres being compared.

In rank zero the form is

\[
 f=D(s,U)+P F(s,U)+Q G(s,U),\qquad J(F,G)\in C^\times.
\]

The common zero set of F,G is finite: any common curve would force
their differentials dependent at a smooth point. For generic fixed s,
the polynomials F(s,U),G(s,U) are therefore coprime. Complete that row
over C[U] by Bezout. The resulting three-variable fibre coordinate has
every value fibre A2, contradicting condition 1. This is a fibre-variety
isomorphism, not a Hessian-preserving nonlinear change.

In rank one the form is

\[
 f=E(U,Z)+s A(U,Z)+V B(U,Z),\qquad J(A,B)\in C^\times.       \tag{1}
\]

The relative critical equations are B=0 and
`grad E+s grad A+V grad B=0`. Their last two equations solve uniquely
for s,V as polynomials in U,Z, since their determinant is J(A,B).
They identify the **whole** C_f with the smooth, nonempty plane curve
C=B=0, carrying the base projection to a polynomial function sigma on C.
On this curve

\[
                    dE=-\sigma\,dA.                         \tag{2}
\]

By condition 2, sigma is a nonconstant unit on every component. The
Keller condition makes dA nowhere zero on C. Hence E is nonconstant
on every component, and the generic number d0 of points E=lambda on C
is a positive integer (at least the number of components).

Also chi(C)<=0. Indeed, on the smooth projective completion of any
component, the nonconstant function sigma has both a zero and a pole.
Neither lies on the affine component, since sigma is a unit there.
At least two points are omitted, so its Euler characteristic
`2-2g-number_of_omitted_points` is nonpositive. Sum over components.
This does not require sigma to be finite or etale, nor the curve to
be rational.

Now fix exactly s=0 in (1). Off B=0, the equation f=lambda solves V,
giving A2 minus C. On B=0, each point E=lambda supplies an affine V
line. For generic lambda,

\[
 \chi(f(0,-)=\lambda)=1-\chi(C)+d_0\ \ge 2.                  \tag{3}
\]

This contradicts condition 3. All Tang alternatives and the automorphic
case have been exhausted. In particular the comparison is at the
exceptional base value, not an extrapolation from generic base fibres.

## 2. Both native exceptional fibres (HC4MYGZ1)

Use s for the retained base, put `R0(p)=R(0,p)` and `H(p)=L*p^2+M*p`.
Assume R0 is not identically zero, and let k be its number of distinct
complex zeros. At s=0 the two graph potentials are exactly

\[
 \begin{array}{ll}
 y=R(s,p):& H(p)+4NpR_0(p)^2+N R_0(p)q,\\
 x=R(s,p):& H(T)+2N R_0(p)r,
       \quad T=p+3R_0(p)q-R_0(p)^3r.
 \end{array}                                                   \tag{4}
\]

For the y chart and p off the k roots, the value equation solves q and
r is free. For the x chart, use (p,T,r) on that same open p set; this
is a regular isomorphism there, with inverse
`q=(T-p+R0^3*r)/(3*R0)`. The value equation solves r and T is free.
In both cases this stratum is `(A1 minus k points) times A1`.
At a root rho of R0, the whole q,r plane has the constant value H(rho).
Consequently the exact Euler function is, for either chart,

\[
 \chi(f(0,-)=\lambda)=1-k+
             \#\{\rho:R_0(\rho)=0,\ H(\rho)=\lambda\}.       \tag{5}
\]

Repeated roots are counted once, and distinct roots with the same H
value each contribute one A2. In particular the generic value is1-k<=1.
The case of a nonzero constant R0 has k=0 and generic fibre A2.
The identically zero R0 is excluded here and is handled separately below.

The same formulas determine the boundary relative critical set. In the
y chart a zero rho of R0 gives

\[
 f_p=H'(\rho)+NqR_0'(\rho),\quad f_q=f_r=0.
\]

In the x chart it gives

\[
 f_p=H'(\rho)(1+3R_0'(\rho)q)+2N R_0'(\rho)r,
                       \quad f_q=f_r=0.                      \tag{6}
\]

There are no critical points off the roots: f_q=N*R0 in the y chart;
in the x chart `f_r+(R0^2/3)*f_q=2N*R0`. Thus if R0 has no simple
root and `R0(-M/(2L)) != 0`, there is **no** relative critical point
at s=0, including nonreduced roots. A simple root instead supplies
the retained literal A1 component. A multiple root with H'(rho)=0
supplies a critical two-plane and forbids unit Hessian.

## 3. Apply the obstruction in all degrees (HC4MYGZ2)

Suppose deg_p R>0 and an actual retained-dual pullback has unit Hessian.
If R0=0, [HC4MYGU0](RELATIVE_COORDINATE_ENDPOINT.md)'s composite-trace
obstruction forbids it. If R0 has a simple root, the previous boundary-line
theorems [HC4MYGR2](CRITICAL_COMPONENTS.md) and
[HC4MYGX0](OPPOSITE_NONLINEAR_FRONTIER.md) give a polynomial gradient
inverse. Multiple roots at H'=0 are excluded by their critical two-planes.
It remains to rule out

\[
 R_0\ne0,\quad R_0\text{ has no simple root},\quad
                         R_0(-M/(2L))\ne0.                    \tag{7}
\]

Here condition 3 of HC4MYGZ0 follows directly from (5), and C_f has
no base-zero point by (6). For every nonzero base its entire critical
set is explicitly classified in
[HC4MYGK1](RELATIVE_CRITICAL_LOCUS.md) and HC4MYGX0:

\[
                3+2sR(s,p)=0,\qquad R_p\ne0,                 \tag{8}
\]

with q,r uniquely recovered and the stated square-coordinate T fixed.
In the opposite chart the additional R=0 locus has no relative critical
point for s nonzero. Since `F_p=2sR_p != 0` in (8), every nonempty
component projects nonconstantly to s. Vertical content factors have
R_p=0 and yield no critical component. This verifies condition 2 without
extending a rational computation chart across any missing point.

For generic fixed nonzero s, the generic three-dual value fibre has Euler
`1+r>=2`, where r is the number of distinct generic p roots of F=3+2sR.
The y-chart stratum proof is [HC4MYGR2](CRITICAL_COMPONENTS.md).
The opposite-chart proof in HC4MYGX0 includes the additional R-zero
strata; their contributions cancel, giving the same1+r. Positive p
degree guarantees r>=1. Actual retained-dual automorphisms preserve
both this generic fibre and the exceptional fibre (5).

All three conditions of HC4MYGZ0 hold. Thus (7) has no unit-Hessian
pullback, irrespective of escapes, singularities in a discarded plane
closure, repeated factors, or degree. This proves HC4MYGZ2. In particular
the former opposite r>=6 finite-escape case is closed, rather than moved
to a larger distinct-root threshold. The full normalization and pole-budget
results remain valid narrower evidence and regression controls.

## 4. The two entire p-only base families (HC4MYGZ3)

Allow p-independent R as well. In **either** entire family y=R(x,p) or
x=R(y,p), existence of a nonautomorphic unit-Hessian retained-dual
member is exactly existence of a unit-ideal plane Keller pair.
Positive p degree contributes no such member by HC4MYGZ2. The base-only
members are exactly [HC4MYGU1/2/3](RELATIVE_COORDINATE_ENDPOINT.md)'s
relative-coordinate endpoint; both orientations realize that same endpoint.
The fixed genuine graph x=1 already realizes it in the opposite family.

This concerns nonautomorphic members. It does not assert that every
positive-p graph is unit-Hessian-free. It neither constructs a unit-ideal
plane Keller pair nor proves plane Keller surjectivity. Over C, absence
of that endpoint is exactly plane Keller surjectivity, as previously
proved; no unrestricted JC2 equivalence is newly claimed.

The arguments extend to every characteristic-zero coefficient field:
embed the finitely generated coefficient field of the graph and actual
automorphism/inverse into C, and descend any polynomial gradient inverse
by uniqueness. Boundary conditions are interpreted over the algebraic
closure. q/r-dependent omitted bases, actual free-base/dual mixing,
uncovered embeddings and the distinct nonlinear-p singleton-zero
properness endpoint remain open. So do unrestricted HC4 and JC2.

[verify_exceptional_fibre.py](verify_exceptional_fibre.py) and its
[receipt](receipt-exceptional-fibre.json) check the two complete boundary
potentials, open-stratum inverse maps, critical equations, Euler identities,
Tang projection adaptations and whole rank-one reconstruction. The
unrestricted curve, topology and exhaustion arguments above are written
proofs. Retained Tang and the classical boundary-line inputs are explicit;
no new classification theorem, formal verification, independent external
review or novelty priority is claimed.
