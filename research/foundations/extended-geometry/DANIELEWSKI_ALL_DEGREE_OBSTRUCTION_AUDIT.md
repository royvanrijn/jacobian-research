# Audit of the all-degree Danielewski contraction obstruction

Consider

\[
X=\mathbb A^1_a\times D,\qquad
D=\{cw=x(x^2+1)\},
\]

with

\[
B=k[c,x,w]/(cw-x^3-x)
\]

and residue bracket

\[
\{c,x\}=c,\qquad
\{c,w\}=3x^2+1,\qquad
\{x,w\}=w.
\]

The proposed all-degree statement is that there are no
\(A,F,G\in B[a]\) satisfying

\[
J=
A_a\{F,G\}-F_a\{A,G\}+G_a\{A,F\}=1.              \tag{1}
\]

The all-degree statement is still open.  The original quadratic
normalization omitted three linear \(x\)-coefficients; restoring them gives
a corrected 30-variable ideal.  Exact modular reconstruction over
\(\mathbb Q\) now gives the unit Gröbner basis, closing the full quadratic
space.  This audit records that repair, sets up the first cubic test, and
gives a cubic family showing why the current two structural gates are
insufficient by themselves.

## 1. This is a restricted Jacobian-conjecture problem

There is an open immersion

\[
\Phi\colon\mathbb A^3_{a,b,c}\longrightarrow X,
\qquad
(a,b,c)\longmapsto
(a,c,bc,b((bc)^2+1)).                              \tag{2}
\]

Over \(\mathbb C\), its complement is the disjoint union

\[
E_+\sqcup E_-,
\qquad
E_\pm=\{c=0,\ x=\pm i\}\simeq\mathbb A^2_{a,w}.    \tag{3}
\]

Moreover,

\[
\Phi^*(da\wedge\omega)=-da\wedge db\wedge dc.
\]

Consequently, a solution of (1) gives a polynomial Keller map

\[
K=(A,F,G)\circ\Phi\colon\mathbb A^3\to\mathbb A^3
\]

with ordinary determinant \(-1\).

This Keller map cannot be an automorphism.  Indeed, if it were, then

\[
s=\Phi\circ K^{-1}\colon\mathbb A^3\to X
\]

would be a section of the étale morphism
\((A,F,G)\colon X\to\mathbb A^3\).  A section of a separated étale
morphism is both open and closed.  Since \(X\) is connected, its image
would be all of \(X\), contradicting
\(\operatorname{im}(s)=\operatorname{im}(\Phi)=X\setminus(E_+\sqcup E_-)\).

Thus:

\[
\boxed{\text{Every solution of (1) produces a counterexample to
the three-dimensional Jacobian conjecture.}}       \tag{4}
\]

Equivalently, (1) asks for a Keller triple in the proper subring

\[
k[a,c,bc,b((bc)^2+1)]\subset k[a,b,c]
\]

which extends étale across both missing boundary planes.  This does not
make the proposed theorem inaccessible, but it explains why a short
cohomological proof would be unusually strong.

The connection is historically natural.  Danielewski completions of
affine open subsets were already used to formulate potential plane
Jacobian counterexamples.  The known literature also contains nonproper
étale endomorphisms of related affine pseudo-planes, so nonproper
étaleness alone is not a contradiction.

## 2. What the two proved structural identities say

Write the top \(a\)-coefficients as

\[
A=A_ra^r+\cdots,\quad
F=F_sa^s+\cdots,\quad
G=G_ta^t+\cdots .
\]

The coefficient of \(a^{r+s+t-1}\) in (1) is

\[
T_{r,s,t}
=rA_r\{F_s,G_t\}
-sF_s\{A_r,G_t\}
+tG_t\{A_r,F_s\}.                                 \tag{5}
\]

For \(r>0\),

\[
T_{r,s,t}
=\frac{A_rF_sG_t}{r}
\frac{
\left\{F_s^r/A_r^s,\ G_t^r/A_r^t\right\}
}{
(F_s^r/A_r^s)(G_t^r/A_r^t)
}.                                                  \tag{6}
\]

Hence the two weighted ratios Poisson-commute when the exponent in (5)
is positive.  They lie in a common rational pencil.

Independently,

\[
\rho\!\left(A\{F,G\}\right)=(a+\lambda,\mu),       \tag{7}
\]

where

\[
\rho\colon B\longrightarrow
H^2_{\mathrm{dR}}(D)
\simeq k[x]/\partial_x((x^3+x)k[x])
\simeq k[1]\oplus k[x].
\]

The recurrence is

\[
[x^{m+2}]=-\frac{m+1}{m+3}[x^m].                  \tag{8}
\]

Two corrections to the informal all-degree heuristic are important.

First, (7) requires one fixed nonzero cohomology class, \([1]\), in a
two-dimensional cohomology group.  It does not require two independent
classes at the next face.  Its two vanishing-cycle periods are nonzero and
opposite, but they are the periods of one class.

Second, (6) classifies the transcendence degree of the leading pencil, not
the pencil itself.  A rational parameter can have poles on \(E_+\) or
\(E_-\), and arbitrary rational pencils on a rational surface are not
reduced to finitely many forms by Lüroth's theorem.

## 3. Quadratic normalization gap

The full ambient-quadratic space is

\[
\operatorname{span}
\{1,a,c,x,w,a^2,c^2,x^2,w^2,ac,ax,aw,cx,cw,xw\}.
\]

Only \(a,c,w\) have nonzero source-linear terms after (2), so a target
linear change normalizes their coefficient matrix.  It does not remove
the coefficient of \(x\), because \(x=bc\) has zero source-linear jet.
The correctly normalized target therefore has the form

\[
(A,F,G)=(a,c,w)
+v_xx+\sum_{j=1}^{10}v_jN_j,                      \tag{9}
\]

with \(v_x,v_j\in k^3\) and \(N_j\) the ten quadratic ambient
monomials.  The existing exact calculation sets \(v_x=0\).

Restoring \(v_x\) gives 33 coefficients before, and 30 variables after,
the same three linear eliminations.  There are again 89 distinct reduced
equations.  The ideal is the unit ideal modulo

\[
32003,\qquad 32009,\qquad 32027,                  \tag{10}
\]

and Singular's exact modular-over-\(\mathbb Q\)
`modGB("slimgb",I,1)` reconstruction gives

\[
\operatorname{GB}(I_{\mathbb Q})=\{1\}.           \tag{11}
\]

An independent `msolve` calculation over \(\mathbb Q\) also reports the
empty variety.  Thus the normalization gap is closed and the full
ambient-quadratic no-go is restored.

## 4. A normalized cubic compatibility witness

The failure of a purely leading-face/de Rham contradiction already appears
in ambient degree three.  Let \(K,L\ne0\) and put

\[
\begin{aligned}
A&=a+w+K a^2c,\\
F&=-w+\frac{15}{4}a x^2+L a^2c,\\
G&=c.
\end{aligned}                                      \tag{12}
\]

The linear Jacobian of \((A,F,G)\) in the local coordinates
\((a,c,w)\) at the origin is one.  The \(a\)-degree profile is
\((r,s,t)=(2,2,0)\), and the top coefficients are

\[
A_2=Kc,\qquad F_2=Lc,\qquad G_0=c.
\]

Thus (5) vanishes identically and the leading ratios collapse to a single
trivial pencil.

Nevertheless, the exact de Rham compiler gives

\[
\rho(A\{F,G\})=(a,0),\qquad
\rho(J)=(1,0).                                     \tag{13}
\]

So (12) satisfies every leading-face and compiled-flux condition.  But
direct expansion gives

\[
\begin{aligned}
J={}&1+\frac{45}{4}x^4+\frac{27}{4}x^2
-15Ka^2c^2x\\
&+(6K+6L)acx^2+(2K+2L)ac
-\frac{15}{2}acx,
\end{aligned}                                      \tag{14}
\]

which is not one.

This is not a counterexample to the proposed theorem.  It is a
counterexample to the missing implication

\[
\text{leading pencil collapse + correct de Rham flux}
\Longrightarrow J=1\text{ is impossible}.          \tag{15}
\]

Any successful descent must retain a pointwise coefficient invariant that
detects the exact error in (14).

## 5. The full normalized cubic test

The space of ambient polynomials of degree at most three in
\((a,c,x,w)\) has 35 monomials.  Its intersection with the hypersurface
ideal in this degree is generated by

\[
cw-x^3-x,
\]

so its image in \(B[a]\) has dimension 34.

After target translation and normalization of the invertible source-linear
jet to \((a,c,w)\), the linear \(x\)-coefficient again remains free.
Each target row has

\[
1+10+20-1=30
\]

residual coefficients.  The universal cubic ansatz therefore begins with

\[
\boxed{90\text{ coefficient variables}.}           \tag{16}
\]

Applying the de Rham compiler before the pointwise determinant expansion
gives

\[
\deg_a\rho_0(A\{F,G\})\le5,\qquad
\deg_a\rho_1(A\{F,G\})\le6.
\]

After discarding the two free constant terms, equation (7) gives

\[
\boxed{11\text{ nonconstant flux equations}:\
5\text{ in }[1]\text{ and }6\text{ in }[x].}        \tag{17}
\]

These equations are necessary but, by (12)--(14), not sufficient.  The
next computation should combine them with the lowest pointwise
coefficients of \(J-1\), then branch by \(a\)-degree profile.  A raw
90-variable Gröbner elimination is not the appropriate first step.

## 6. Boundary-root and pole audit

At a missing plane \(E_\alpha\), \(\alpha=\pm i\), the bracket restricts
to

\[
\{H_1,H_2\}|_{E_\alpha}
=P'(\alpha)
(H_{1,c}H_{2,w}-H_{1,w}H_{2,c})
+w(H_{1,x}H_{2,w}-H_{1,w}H_{2,x}).                \tag{18}
\]

But

\[
P'(i)=P'(-i)=-2.                                   \tag{19}
\]

Therefore the two local normal ledgers do not differ by a sign in their
leading \(P'\)-term.  They can only be separated by evaluation of the
odd/even \(x\)-parts of the coefficients or by global period data.

For an affine-linear leading coefficient

\[
H=u c+v x+q w+d,
\]

the condition that \(H\) vanish generically on \(E_\alpha\) is

\[
q=0,\qquad d=-v\alpha.                             \tag{20}
\]

Thus an \(a^2\)-leading coefficient can already put a pole of the ratios
in (6) on one boundary component over \(\mathbb C\).  Over a
\(\mathbb Q\)-defined chart, conjugation forces the corresponding
conditions to occur together unless the coefficients are extended.
For \(a\)-degree one the leading coefficient is quadratic on \(D\), and
there are still more boundary-vanishing types.

The pole scenario in the proposed lead is therefore genuine.  It is not
removed by Lüroth reduction.

## 7. Gauge caution

Constant target \(\mathrm{SL}_3\) changes, target translations, and affine
changes \(a\mapsto\lambda a+\mu\) preserve the ambient cubic filtration.
General Danielewski automorphisms do not.  For example, exponentiating the
standard locally nilpotent derivation sends

\[
x\mapsto x+\tau c,\qquad
w\mapsto \frac{P(x+\tau c)}{c},
\]

and substitution into a cubic target can raise its ambient degree.

Hence a degree-three exhaustion may quotient only by the
filtration-preserving stabilizer, or it must track the degree change
explicitly.  Quotienting by the full Danielewski automorphism group before
enumeration can remove or import cubic candidates.

## 8. Verdict and next lemma

The proposed all-degree no-go theorem is plausible but unproved.  The
current evidence supports the following narrower conclusion:

\[
\boxed{
\begin{gathered}
\text{the corrected degree-2 ideal is exactly unit over }\mathbb Q;\\
\text{degree }3\text{ has 90 normalized variables and 11 flux gates};\\
\text{leading collapse and flux alone are compatible.}
\end{gathered}}
\]

The decisive missing statement must be pointwise.  A useful target is:

> For each cubic \(a\)-degree profile and each boundary-pole type of the
> leading pencil, the coefficient ideal generated by (5), (7), and a
> bounded set of lowest pointwise coefficients of \(J-1\) is the unit
> ideal.

This is finite, respects the potential one-boundary pole failure, and can
be attacked before the full determinant ledger.  If it fails, its
components provide structured cubic counterexample ansätze rather than an
uninterpreted 90-variable search.

## 9. Reproduction

Run

```bash
.venv/bin/python scripts/verify_danielewski_cubic_frontier.py
.venv/bin/python scripts/verify_danielewski_full_quadratic_obstruction.py
```

The first checker verifies the open-immersion Jacobian identity, the 90-variable
cubic count, the 11 compiled flux equations, the normalized compatibility
family (12), and the equal local \(P'\)-values versus independent global
periods at \(x=\pm i\).  The second restores the three omitted \(x\)
coefficients, verifies the 30-variable/89-equation ledger, and reproduces
the exact characteristic-zero unit basis.

## 10. Literature context

- A. Dubouloz, B. Kunyavskii, and A. Regeta,
  *Bracket Width of Simple Lie Algebras*, Documenta Mathematica 26
  (2021), 1601--1627.  Its Danielewski section records the residue Poisson
  algebra, the nonzero top de Rham class, and a Jacobian-type hypothesis
  for rational maps induced from pairs of Danielewski functions.
- A. Dubouloz and K. Palka,
  *The Jacobian Conjecture fails for pseudo-planes*, Advances in
  Mathematics 339 (2018), 248--284.  It provides nonproper étale
  endomorphisms on related affine surfaces and shows that nonproper
  étaleness is not, by itself, contradictory.
- A. Dubouloz, *Surfaces de Danielewski et Conjecture Jacobienne*,
  workshop abstract, Basel (2009).  It explicitly formulates the
  Danielewski-completion viewpoint on potential Jacobian counterexamples
  and warns that naive differential cohomology does not decide the
  étale-morphism question.
