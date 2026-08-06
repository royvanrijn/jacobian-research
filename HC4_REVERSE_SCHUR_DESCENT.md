# Reverse Schur descent for `HC4`

## Status

This note implements the scalar reverse-Schur programme and closes its
**cone-pencil stratum**.  It does not prove `HC4`.

The proved result is:

> **Theorem `HC4RSD1` (scalar cone-pencil obstruction).**  Let \(K\) have
> characteristic zero, let \(x=(x_1,\ldots,x_4)\), and put
> 
> \[
>  \Phi(t,x)=\frac{\lambda}{2}t^2+tA(x)+B(x).
> \]
> 
> Suppose that
> 
> \[
>  \det\operatorname{Hess}_{t,x}\Phi=c\in K^\times,
>  \qquad
>  \det\operatorname{Hess}_x(B+sA)=0,                 \tag{0.1}
> \]
> 
> and that the generic Hessian pencil has an \(x\)-constant kernel line:
> there is a nonzero \(v(s)\in K(s)^4\) with
> 
> \[
>  \operatorname{Hess}_x(B+sA)v(s)=0.                \tag{0.2}
> \]
> 
> Then \(A\) and \(B\) have a common constant Hessian-kernel direction
> \(v_0\in K^4\) with \(D_{v_0}A\ne0\).  Every scalar Schur descendant
> 
> \[
>  \psi_{\kappa,\mu}
>  =B+\frac{\kappa}{2}A^2+\mu A,
>  \qquad \kappa\ne0,                                \tag{0.3}
> \]
> 
> has injective gradient.  In particular it has no collision and cannot
> meet an affine-degree-two or affine-degree-three projective-polar packet.

The hypothesis (0.2) is automatic when \(A\) and \(B\) are homogeneous
four-variable forms of the same degree, by the low-dimensional
Gordan--Noether theorem used elsewhere in this repository.  It is also
directly checkable for a given polynomial pencil.  The theorem does **not**
assert that an arbitrary
nonhomogeneous singular-Hessian polynomial in four variables is a cone.
Classifying pencils with an \(x\)-moving kernel line remains the scalar
frontier.

The exact equation builders are in
[`jcsearch/reverse_schur_descent.py`](jcsearch/reverse_schur_descent.py).
They support both scalar and simultaneous matrix pivots.  The verifier
[`scripts/verify_hc4_reverse_schur_descent.py`](scripts/verify_hc4_reverse_schur_descent.py)
writes the generated ledger
[`artifacts/generated-results/hc4_reverse_schur_descent.json`](artifacts/generated-results/hc4_reverse_schur_descent.json).

## 1. Exact scalar equations

Put

\[
 g=\nabla A,
 \qquad M(s,x)=\operatorname{Hess}_x(B+sA),
 \qquad D(s,x)=\det M(s,x).
\]

The parent Hessian is

\[
 \operatorname{Hess}\Phi
 =\begin{pmatrix}
   \lambda&g^{\mathsf T}\\
   g&M(t,x)
  \end{pmatrix},
\]

so

\[
 \det\operatorname{Hess}\Phi
 =\lambda D(t,x)-g^{\mathsf T}\operatorname{adj}(M(t,x))g. \tag{1.1}
\]

Thus a singular pencil removes the first term, independently of the parent
quadratic coefficient.  The exact bordered-unit equation is

\[
 \boxed{-g^{\mathsf T}\operatorname{adj}(M(s,x))g=c.}       \tag{1.2}
\]

For \(\lambda\ne0\), fixing the pivot-gradient value \(y\) gives the critical
point

\[
 t^*=\frac{y-A}{\lambda}
\]

and the reduced potential

\[
 \psi_y=B-\frac{(A-y)^2}{2\lambda}.                         \tag{1.3}
\]

Exact block elimination gives

\[
 \boxed{
 \det\operatorname{Hess}\Phi(t^*,x)
 =\lambda\det\operatorname{Hess}\psi_y.}                   \tag{1.4}
\]

Equation (1.3) is (0.3) with

\[
 \kappa=-\lambda^{-1},\qquad \mu=y/\lambda.
\]

The same family (0.3) is obtained by quadratically repairing an affine
pivot.  Under (0.1), the exact determinant is

\[
 \det\operatorname{Hess}\psi_{\kappa,\mu}=-\kappa c.       \tag{1.5}
\]

The implementation keeps (0.1), (1.2), and the collision equations as
separate coefficient systems.  In particular, a bounded search can report
which gate failed instead of treating a zero-dimensional Gröbner basis as a
classification.

## 2. The corank gate

The generic rank of \(M\) is forced before any classification.

If \(\operatorname{rank}M\le2\), adjoining one border row and column raises
rank by at most two, so the \(5\)-by-\(5\) parent Hessian has rank at most
four.  This contradicts \(c\ne0\).  Combined with \(D=0\), this gives

\[
 \boxed{\operatorname{rank}_{K(s,x)}M=3.}                    \tag{2.1}
\]

Hence the adjugate has rank one.  Under (0.2), choose the kernel generator
\(v(s)\) and write

\[
 \operatorname{adj}(M)=q(s,x)v(s)v(s)^{\mathsf T}.           \tag{2.2}
\]

Over the UFD \(K(s)[x]\), equations (1.2) and (2.2) give

\[
 -q(s,x)\bigl(v(s)^{\mathsf T}\nabla A\bigr)^2=c.            \tag{2.3}
\]

Both factors on the left must be units.  Therefore

\[
 v(s)^{\mathsf T}\nabla A\in K(s)^\times,
 \qquad \operatorname{Hess}(A)v(s)=0.                       \tag{2.4}
\]

Substituting (2.4) into \(M(s,x)v(s)=0\) also gives

\[
 \operatorname{Hess}(B)v(s)=0.                              \tag{2.5}
\]

Let

\[
 V=\{v\in K^4:\operatorname{Hess}(A)v
                 =\operatorname{Hess}(B)v=0\}.
\]

Equations (2.4)--(2.5) put \(v(s)\) in \(V\otimes_K K(s)\).
Moreover, the constant linear functional
\(v\mapsto D_vA\) is nonzero on that scalar extension.  It is therefore
nonzero on \(V\) itself.  Choose \(v_0\in V\) with

\[
 \alpha=D_{v_0}A\ne0.                                       \tag{2.6}
\]

This is the promised common constant direction.  Notice where the cone
hypothesis was used: it made \(v(s)\) independent of \(x\), so differentiating
\(v(s)^{\mathsf T}\nabla A\) gives exactly
\(\operatorname{Hess}(A)v(s)\).  This step is unavailable for a moving
polynomial vector field \(v(s,x)\).

## 3. Integrability and collision exclusion

Make a linear change of the four \(x\)-coordinates so that
\(v_0=\partial_z\) and write \(u=(u_1,u_2,u_3)\).  Integrating (2.6) and
its \(B\)-analogue gives

\[
 A=\alpha z+a(u),\qquad B=\beta z+b(u),\qquad \alpha\ne0.   \tag{3.1}
\]

The pencil Hessian is

\[
 M(s,x)=
 \begin{pmatrix}
  \operatorname{Hess}_u(b+sa)&0\\
  0&0
 \end{pmatrix}.
\]

Expanding (1.2) along the last row and column yields

\[
 \boxed{
 \det\operatorname{Hess}_u(b+sa)=-c/\alpha^2\in K^\times
 }                                                           \tag{3.2}
\]

for every \(s\).

Now suppose two points \(p,q\) have equal gradients under (0.3).  Its last
gradient coordinate is

\[
 \partial_z\psi_{\kappa,\mu}
 =\beta+\alpha(\kappa A+\mu).                               \tag{3.3}
\]

Since \(\alpha\kappa\ne0\), equality in (3.3) gives \(A(p)=A(q)\).  Put

\[
 s_0=\kappa A(p)+\mu=\kappa A(q)+\mu.
\]

The other three equal-gradient equations become exactly

\[
 \nabla_u(b+s_0a)(u_p)=\nabla_u(b+s_0a)(u_q).              \tag{3.4}
\]

By (3.2), this is a three-variable constant-Hessian gradient map.  The known
truth of `HC3` makes it injective, so \(u_p=u_q\).  Equation \(A(p)=A(q)\)
and (3.1) then give \(z_p=z_q\).  Thus \(p=q\), proving the theorem.

This is stronger than excluding one normalized collision: every geometric
fiber contains at most one point. After base change to an algebraic closure,
Ax--Grothendieck makes the injective gradient a polynomial automorphism.
Consequently every fiber algebra is reduced of length one.

## 4. Projective-polar and compactified-gradient intersection

The current degree-five `HC4` atlas retains 318 affine-degree-two and 306
affine-degree-three numerical signatures after the smooth rank-three vertex
obstruction.  The cone-pencil theorem reaches none of them:

| reverse-Schur stratum | affine degree 2 | affine degree 3 |
|---|---:|---:|
| scalar cone pencil | 0 | 0 |

This is not a numerical deletion from the atlas. It is a structural
preimage statement: every classified reverse-Schur descendant is injective,
whereas the packets encode collision-normalized candidates. Their
intersection is therefore empty before lower Rees or torsion data are
computed; no degree bound on the descendant is being asserted.

Likewise, no rational reconstruction is attempted.  There is no surviving
scalar cone component from which to reconstruct a point over
\(\mathbb Q\).

## 5. Matrix-pivot implementation

For \(r\) pivot variables, let \(A=(A_1,\ldots,A_r)^{\mathsf T}\), let
\(\Lambda\) be an invertible symmetric matrix, and put

\[
 \Phi(t,x)=\frac12t^{\mathsf T}\Lambda t+t^{\mathsf T}A(x)+B(x).
\]

With

\[
 J=D_xA,
 \qquad M(s,x)=\operatorname{Hess}_x(B+s^{\mathsf T}A),
\]

the parent Hessian is

\[
 \begin{pmatrix}
  \Lambda&J\\
  J^{\mathsf T}&M(t,x)
 \end{pmatrix}.                                             \tag{5.1}
\]

Fixing the pivot-gradient value \(y\) gives

\[
 t^*=\Lambda^{-1}(y-A)
\]

and

\[
 \psi_y
 =B-\frac12(A-y)^{\mathsf T}\Lambda^{-1}(A-y).              \tag{5.2}
\]

The implemented exact identity is

\[
 \boxed{
 \det\operatorname{Hess}\Phi(t^*,x)
 =\det(\Lambda)\det\operatorname{Hess}\psi_y.
 }                                                           \tag{5.3}
\]

The equation builder generates every minor imposing

\[
 \operatorname{rank}M(s,x)\le4-r,                            \tag{5.4}
\]

then extracts its coefficients in \(s,x\).  It also exposes the symmetry and
third-derivative equations needed when a search begins from a matrix pencil
rather than from potentials \(A_i,B\).

No classification of the moving kernel planes in (5.4) is claimed.  This is
the main structured extension left by the implementation.

## 6. Reproduction and remaining frontier

Run:

```bash
.venv/bin/python scripts/verify_hc4_reverse_schur_descent.py
```

The command verifies:

1. the scalar and matrix block-Schur identities;
2. the generic scalar corank-one gate;
3. the rank-one adjugate and bordered-unit factorization on an exact
   common-kernel calibration;
4. rejection of the integrable moving-kernel pencil
   \(B=x_1x_3+z^2/2, A=x_2x_3\), whose kernel is
   \((-s,1,0,0)\) but whose bordered determinant is \(x_3^2\);
5. the empty intersection of the scalar cone stratum with the live
   affine-degree-two/three atlas rows; and
6. the simultaneous two-pivot equation schema on an exact corank-two
   calibration.

The next exact classification problems are now separated:

- `HC4RSD2` closes every affine-in-\(x\) kernel line whose normalized line
  is fixed in the pencil parameter, and `HC4RSD3` excludes the complete
  parameter-moving affine branch;
- `HC4RSD4` closes the first degree-unbounded nonlinear family
  \(v=(P(z,w),1,0,0)\): the bordered unit forces its transverse dependence
  onto one linear form. `HC4RSD5` extends this to every fixed primitive
  two-component generator in a constant support plane: Piola eliminates
  active-variable dependence, then a constant-determinant polynomial frame
  forces the component pair onto an affine line and reduces it to
  `HC4RSD4`;
- fixed kernels with three or four nonlinear components, and
  parameter-moving nonlinear kernel lines;
- scalar pencils with \(D(s,x)\ne0\), where the two terms in (1.1) cancel;
- matrix pencils with jointly moving kernel planes; and
- only after one of those survives, its lower compactified-gradient algebra
  and rational reconstruction.

`HC4RSD6` in
[`HC4_AFFINE_PIVOT_COVERAGE_GATE.md`](HC4_AFFINE_PIVOT_COVERAGE_GATE.md)
now supplies the first converse coverage gate. An affine scalar singular
pivot exists exactly when a constant covector has nonzero constant
inverse-Hessian norm. On the open rank-three quintic packet this requires
the degree-eight Schur vector to have coefficient rank at most two. The
next finite calculation is the intersection of that rank locus with the
nonsquarefree Hessian-discriminant and Schur-divisibility equations, rather
than an unrestricted nonlinear-kernel classification.

For inherited collisions, `HC4RSD7` closes even that residual affine
coverage locus. The metric numerator is the ternary Hessian determinant on
each affine pivot fiber; `HC3` makes the tangential gradient injective there.
Since parent-gradient equality at a common pivot value forces equal pivot
values, no nontrivial collision transfers through any affine zero-corner
scalar pivot, whether the reduced Hessian is singular or is repaired in the
nonsingular exact-remainder branch. Nonlinear, mixed/coisotropic, matrix,
and nonlinear exact-remainder pivots remain.

`HC4RSD8` in
[`HC4_QUADRATIC_PIVOT_RANK_OBSTRUCTION.md`](HC4_QUADRATIC_PIVOT_RANK_OBSTRUCTION.md)
starts the nonlinear scalar branch. If \(A\) is quadratic and the pencil is
identically singular, then

\[
\operatorname{rank}\operatorname{Hess}A\leq2.
\]

Rank four is incompatible with pencil singularity. In rank three, splitting
off the null direction and clearing the active \(3\)-by-\(3\) block gives
\(P^2=-c\det(sQ_3+H)\), contradicting parity of the degree in \(s\). This
argument allows an \(x\)-moving kernel line. The bordered unit additionally
forces \(dA\) to be nonzero on the kernel of its constant Hessian, so every
survivor has affine normal form
\(A=w+u^{\mathsf T}Q_ru/2\). The remaining quadratic target is therefore the
rank-one/rank-two moving-kernel locus; the common-kernel part is already
covered by `HC4RSD1`--`HC4RSD5`. In rank two, the passive binary Hessian of
\(B\) is singular. If it vanishes, the pencil determinant is the square of
the active/passive cross determinant and the surviving kernel lies in a
fixed two-plane, so `HC4RSD5` applies. Thus the only new rank-two stratum has
passive Hessian rank one.

`HC4RSD9` in
[`HC4_QUADRATIC_RANK_TWO_PIVOTS.md`](HC4_QUADRATIC_RANK_TWO_PIVOTS.md)
closes that last rank-two stratum. In hyperbolic coordinates \(A=xy+w\),
the remaining determinant faces give the complete form

\[
B=xz+\frac{\rho}{2}(y+h(x)A)^2
  +\beta(x)y+\gamma(x)A+\delta(x).
\]

Every descendant with nonzero quadratic repair coefficient has a triangular
gradient inverse.

`HC4RSD10` in
[`HC4_QUADRATIC_RANK_ONE_PIVOTS.md`](HC4_QUADRATIC_RANK_ONE_PIVOTS.md)
closes the final rank-one stratum. Normalize \(A=x^2/2+w\), and let \(E\)
be the passive \(3\)-by-\(3\) Hessian of \(B\). The leading pencil and parent
faces give

\[
\det E=0,
\qquad
a^{\mathsf T}\operatorname{adj}(E)a=0,
\quad a=(0,0,1)^{\mathsf T}.
\]

Passive ranks zero and two contradict the generic corank-one bordered unit,
so \(E\) has rank one. The rank-one polynomial-Hessian normal form and the
exact parent identity

\[
\det\operatorname{Hess}\Phi
=\rho\det(a,d,\ell)^2
\]

make the passive frame a polynomial unit. Its Wronskian vanishes, fixing the
projective direction of the first two entries of \(\ell\), and exact
integration gives

\[
B=xz+\frac{\rho}{2}(y+h(x)w)^2
  +\alpha(x)y+\gamma(x)w+\delta(x).
\]

Every descendant again has a triangular polynomial gradient inverse. Thus
all pivot-Hessian ranks zero through four are closed for quadratic scalar
pivots with an identically singular reduced pencil. What remains is the
higher-degree nonlinear scalar branch, nonsingular pencils with exact
determinant-term cancellation, and genuinely mixed/coisotropic or moving
matrix pivots.

HC4RSD11--HC4RSD16 in
[`HC4_SCALAR_CANCELLATION_DICHOTOMY.md`](HC4_SCALAR_CANCELLATION_DICHOTOMY.md)
resolve most of that exact-cancellation qualification. If the scalar pivot
corner is nonzero, completing the square is an exact equivalence with a
four-variable constant-Hessian pencil \(\psi+sA\), including its gradient
collisions. If the corner is zero and \(D_vA\) is a nonzero constant, the
parent gradient factors fiberwise through a ternary constant-Hessian
gradient and is a polynomial automorphism. The bordered unit puts every
quadratic zero-corner pivot in this graph-coordinate class.

For a quadratic nonzero-corner pencil, pivot-Hessian rank four is
impossible. Rank three has the complete triangular form

\[
A=xy+\frac12z^2,\qquad
\psi=wx+y(\alpha z+\beta(x))+G(x,z).
\]

In rank two, the passive-rank-one component reduces to HC2, while passive
rank zero is exactly the cotangent lift

\[
\psi=z\,b(x,y)+w\,c(x,y)+D(x,y),
\qquad \det J(b,c)\in K^\times,
\]

so invertibility is equivalent to JC2. In rank one, the two ternary
singular-Hessian normal forms reduce to HC2 or the same cotangent JC2
packet. The bordered unit freezes a rationally \(x\)-moving kernel in the
constant-kernel type. In the exceptional type its transverse-coordinate
coefficient likewise forces the distinguished covector to be projectively
constant, unless the chart has already collapsed to the constant-kernel
type. Thus every quadratic scalar pencil reduces to HC2 or exactly JC2.

HC4RSD17--HC4RSD40 in
[`HC4_HIGHER_DEGREE_PENCIL_OBSTRUCTIONS.md`](HC4_HIGHER_DEGREE_PENCIL_OBSTRUCTIONS.md)
continue the same cancellation analysis without assuming a quadratic
direction. The global pencil is exactly a polynomial nilpotent deformation

\[
N=(\operatorname{Hess}\psi)^{-1}\operatorname{Hess}A,
\qquad N^{\mathsf T}\operatorname{Hess}\psi
=\operatorname{Hess}\psi\,N,
\]

with both metric products Hessian-integrable. Every direction of generic
Hessian rank one reduces, in every degree, to HC2 or the JC2 cotangent
packet. For a cubic direction whose leading Hessian has rank three, the
sole moving-kernel form is

\[
A=wz+y b(z)+G(x,z),
\qquad v=(0,1,0,-b'(z)).
\]

The (s^3)-face makes (v) null for the base Hessian. Exact integration
and the next two faces force a pair of incompatible polynomial units when
(b''\ne0). Hence the moving cubic packet is empty. In the remaining
constant-kernel component, the full ternary bordered determinant restricts
to a binary tangent pencil. Its two singular forms are proportional and
share one ruling. A fixed ruling makes the border coefficient a cylinder,
which reduces to HC2 or exactly the JC2 cotangent lift. The universal-field
identity fixes all homogeneous rulings in every degree. The complete cubic
and quartic correction calculations fix every arbitrary border coefficient
through degree four; in particular the only apparent non-pure quartic Schur
solution dies by the coefficient (-3b^4x^2y^2z^2/4). The quintic
simple-root square and repeated-root Schur charts then close every non-pure
leading quintic, including the two apparent exceptions at their next faces.
The pure-fifth passive-Hessian chart is closed by HC4BL5: its curved
corrections have unit ideals and every surviving lower tail aligns with one
constant passive form. Thus every degree-five border coefficient closes.
The simple-root square also closes squarefree leading binary forms in every
degree. Its double-root valuation continuation closes the generic
discriminant stratum as well. The arbitrary-multiplicity root-order formula
then exposes every resonance. In degree six, exact weighted faces close all
non-pure leading sextics. For the remaining pure-sixth top, the first three
collision faces stabilize one constant passive direction through the
quintic correction. If that correction has nonzero passive curvature, the
factorized lower-direction equation and its repeated-root charts force the
whole potential to remain binary. Completing the square in the remaining
passive-affine quintic term then closes every curved completed quartic. The
last two-linear-form tower through degrees five and four is also closed:
its independent packet has incompatible finite, infinity, and base faces,
while every dependent-rank endpoint is a fixed cylinder by the exact
(D D''-2(D')^2) and Wronskian-square obstructions. Thus every scalar
degree-six leading direction closes. The valuation sieve and seventeen
complete discriminant-open weighted faces then close every non-pure septic,
including all same-weight cubic transverse tails. The pure-seventh residue
has first faces (49x^{12}\det\operatorname{Hess}_{y,z}(c_6)) and
((49/2)x^{12}(H_6)_{yy}(2(c_5)_{zz}-(8/7)k^2x^3)), giving the normal form
(c_6=H_6(x,y)+kx^5z) and the forced curved-chart quintic square. The complete
degree-eighteen face is a Hessian factor minus
((7x^2(P_4)_y-4k(H_6)_y)^2); it kills (z^4), fixes the quartic (z^3)-tail to
(k^3xz^3/49), and leaves one binary root-divisibility identity. For
(k\ne0), the immutable low-(x) coefficients force ((H_6)_{yy}) to contain
(x^2); two terminal product--square factorizations leave exactly six
packets, all killed in the next two faces. At (k=0), retaining the cubic
(z^3)-tail gives a coupled Wronskian whose complete nonzero solution is the
ordered two-linear-form packet
((H_6)_{yy}=LE^3,(P_4)_y=LE^2,D_2=LE); all five projective charts close.
The zero-((P_4)_y) recurrence also closes by four normalized coefficients
and the exact nonlinear-coordinate identities. On the remaining
passive-affine sextic boundary, a shifted passive-Hessian face closes every
curved quintic correction. Its two-linear-form tower is now reduced to
eight degree-fifteen quartic direction packets. Their common quartic-polar
rank-one system has only a zero stratum and a square-Hessian resonance. The
resonances all have immutable lower coefficients, while the only nonempty
zero strata satisfy a global affine-transverse cylinder identity. Thus all
eight packets close and every scalar degree-seven leading direction is
fixed. The next degree-based scalar frontier is the repeated-root locus in
degree at least eight.
