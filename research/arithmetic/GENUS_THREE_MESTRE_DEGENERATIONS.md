# An obstruction to genus-three degenerations of the balanced Mestre construction

The standard six affine sections in Mestre's construction do not provide
a genus-three hyperelliptic curve merely by specializing their constants.
The obstruction below treats the whole balanced affine ansatz, rather than
a finite set of Mestre parameters. It does not exclude more general marked
norm equations or settle the geometrically simple torsion construction.

The primary construction is described in
[Levin–Yan–Wiljanen, Section 5](sources/levin-yan-wiljanen-large-three-rank.pdf),
printed pp. 10–11, following
[Mestre](sources/mestre-icm-1994.pdf).
There one forms a sextic $p$, a quadratic $g$, and the cubic $r=g^3-p$.
The six roots of $p$ give points on $Y^3=r(X)$. The displayed affine
Mestre family has three positive and three negative slopes, equal sums
of squares after centering, and constant $r_3=1$. The theorem below
uses only polynomial identities and the complete rational points on
$y^2=x^3+1$; no torsion-independence assertion is imported from that
construction.

## The complete rational-point input

**Lemma HG15a.** The rational points on $E:y^2=x^3+1$ are
\[
 E(\mathbb Q)=\{O,(-1,0),(0,\pm1),(2,\pm3)\}.
\]

The retained checker uses SageMath 10.9 with rational rank databases
disabled. The exact mwrank upper bound and proved rank are zero. The
rational torsion subgroup has invariant factor [6], and its six points
are exactly those displayed. This is an internal computational proof
input, not an independently replayed or formally verified result.

## The balanced affine ansatz

Let $\alpha_1,\alpha_2,\alpha_3,\beta_1,\beta_2,\beta_3\in\mathbb Q$
be such that each triple consists of distinct numbers and
\[
 \sum_i\alpha_i=\sum_i\beta_i=0,\qquad
 \prod_i(Z-\alpha_i)=Z^3+PZ-Q,\qquad
 \prod_i(Z-\beta_i)=Z^3+PZ-R.
\]
In particular, $P=-\tfrac12\sum_i\alpha_i^2<0$.
Suppose $S=Q+R=c^3$ for a nonzero $c\in\mathbb Q$.
With an indeterminate $t$, set
\[
 \begin{split}
 p(X,t)&=\prod_i(X-\alpha_i-t)\prod_i(X-\beta_i+t),\\
 g(X,t)&=X^2-t^2+\frac{2P}{3},\\
 r(X,t)&=g(X,t)^3-p(X,t),\\
 \mathcal D(t)&=\operatorname{disc}_X r(X,t).
 \end{split}
\]

**Theorem HG15.** The polynomial $\mathcal D$ has degree ten,
leading coefficient $1024P^4$, and no multiple rational root.
For every $a\in\mathbb Q^\times$, the normalization of
$y^2=a\mathcal D(t)$ is therefore never a genus-three curve.

Here the normalization is taken componentwise if the equation is
geometrically reducible. The assertion does not say that every
$\mathcal D$ is squarefree: multiple roots of degree at least two
are not excluded. Nor does it exclude changing the six root functions,
dropping the cube condition on $S$, or using a different torsion
construction.

### Coefficient identities

Put $\delta=Q-R$. Direct expansion gives
\[
 \begin{split}
 r_3&=S,\\
 r_2&=-4Pt^2+3\delta t+\frac{P^2}{3},\\
 r_1&=(3t^2+P)S,\\
 r_0&=4Pt^4-\frac{P^2t^2}{3}
       +\delta(t^3+Pt)+\frac{8P^3}{27}
       -\frac{S^2-\delta^2}{4}.
 \end{split}                                                    \tag{15}
\]
Only $-4r_2^3r_0$ in the cubic discriminant can have degree ten.
Its leading coefficient is $1024P^4\ne0$.

### A multiple discriminant root forces a collision at $g=0$

Suppose $t_0\in\mathbb Q$ is a multiple root of $\mathcal D$.
The cubic $r(X,t_0)$ cannot have a triple root $\xi$.
Such a root would be rational, and
$r(X,t_0)=c^3(X-\xi)^3$. At each of the three distinct rational
numbers $x_i=\alpha_i+t_0$ we would have
$g(x_i,t_0)^3=c^3(x_i-\xi)^3$. Injectivity of cubing on $\mathbb Q$
then gives three distinct zeros of the monic quadratic
$g(X,t_0)-c(X-\xi)$, which is impossible.

Thus $r(X,t_0)$ has a double root $\xi$ and a distinct simple root.
The double root is rational. For a cubic with one double and one
simple root, the first variation of its discriminant is a nonzero
scalar times the variation of its value at the double root.
Consequently $\mathcal D'(t_0)=0$ implies
\[
 r(\xi,t_0)=r_X(\xi,t_0)=r_t(\xi,t_0)=0.                         \tag{16}
\]
The case $\xi=0$ is impossible, since $r_X(0,t_0)=PS\ne0$.
Scale $X,t,\alpha_i,\beta_i,c$ by $\xi$ and $P$ by $\xi^2$.
All hypotheses persist, so it suffices to treat $\xi=1$.
Below write $t$ for this normalized $t_0$.

Let $A=3t^2+P+3$ and
\[
 H=A^2-36t^2
   =\bigl(P+3(t-1)^2\bigr)\bigl(P+3(t+1)^2\bigr).
\]
The two derivative equations in (16) are linear in $S,\delta$, with
coefficient matrix
\[
 \begin{pmatrix}A&6t\\6t&A\end{pmatrix}.
\]
If $H\ne0$, solving them gives
\[
 S=\frac{2P(P+12t^2)(15t^2-P-3)}{3H},
\]
and substituting both solutions into $r(1,t)$ gives
\[
 r(1,t)=-\frac{P(P+12t^2)(3t^2-2P-3)^3}{27H}.                \tag{17}
\]
Since $P,S\ne0$, equations (16)–(17) force
$3t^2-2P-3=0$. This is precisely $g(1,t)=0$.

The singular linear system must also be considered. If
$P=-3(t-1)^2$, the consistency condition is
\[
 r_t(1,t)-r_X(1,t)
       =-6(t+1)(t-1)^2(3t-1)^2=0.
\]
At $t=1$, $P=0$, which is excluded. At $t=1/3$, $g(1,t)=0$.
At $t=-1$, $P=-12$ and the derivative equations force
$S+\delta=32$, hence $Q=16$; but
$Z^3-12Z-16=(Z-4)(Z+2)^2$, contrary to distinctness of the
$\alpha_i$. The other case $P=-3(t+1)^2$ is symmetric:
the consistency factor is $-6(t-1)(t+1)^2(3t+1)^2$;
its exceptional nonzero-$P$ solution gives $R=16$ and a repeated
$\beta_i$. Thus in every admissible case,
\[
                              g(\xi,t_0)=0.                 \tag{18}
\]

Now $p=g^3-r$ and (16)–(18) imply $p=p_X=p_t=0$.
Since the factors within each parallel triple are distinct, exactly
one positive-slope and one negative-slope root meet at $\xi$.

### The collision contradicts the Fermat cubic

Let $a,b$ be the other two positive-slope roots at $t_0$, and
$d,e$ the other two negative-slope roots. Then
\[
 a+b=3t_0-\xi,\quad d+e=-3t_0-\xi,\qquad
 a^2+b^2=d^2+e^2=2\xi^2.                                   \tag{19}
\]
Indeed each triple has square sum $-2P+3t_0^2=3\xi^2$,
by (18). Cauchy's inequality applied to both pairs gives
$|3t_0/\xi|\le1$.
Newton's identities, or a direct sum of cubes using (19), give
\[
 \frac{S}{\xi^3}=9\left(\frac{t_0}{\xi}\right)^2-1.
\]
As $S=c^3$, the pair
\[
 (x,y)=\left(\frac c\xi,\frac{3t_0}{\xi}\right)
\]
is a rational point on $y^2=x^3+1$ with $x\ne0$ and $|y|\le1$.
Lemma HG15a leaves only $(x,y)=(-1,0)$. Thus $t_0=0$,
and (19) implies $(a-b)^2=3\xi^2$, impossible for nonzero rational
$\xi,a,b$. This excludes every multiple rational root of
$\mathcal D$.

Finally write $\mathcal D=h^2f$ over $\mathbb Q[t]$, where $f$ is
squarefree. A genus-three normalization would require
$\deg f=8$, since $\deg\mathcal D=10$ is even. This would give
$\deg h=1$ and hence a multiple rational root, already excluded.
A nonsquare constant twist $a$ changes none of these genus calculations.
If $f$ is constant, the equation supplies no geometrically integral
hyperelliptic curve of genus three.

## Consequences for the construction campaign

The six affine roots in Levin–Yan–Wiljanen Section 5 have equal
square sums because their cubic leading coefficient $r_3$ is constant.
A translation of $t$ makes the two triples separately centered; a
rational common affine change of $t$ handles common opposite slopes.
Whenever both constant triples remain distinct, their nonzero
cube-leading-coefficient case is therefore covered by HG15.
Specializations with duplicate constant sections are excluded from the
theorem and may lose the required torsion independence.

The alternative four-constant-root construction with the pair
$(\sqrt T,-\sqrt T)$ is outside HG15. Removing a double discriminant
factor there does produce genus-three curves; it does not by itself
prove that three rational torsion directions survive. In particular,
point counts divisible by 27 and absolutely simple finite-field
reduction do not establish rational three-rank three.
The remaining construction targets are kept in the
[marked-torsion note](GENUS_THREE_MARKED_TORSION_DEFORMATIONS.md).

## Exact verification and its scope

Run from the repository root:

    sage -python research/arithmetic/verify_genus_three_mestre_degenerations.py

The checker expands the generic coefficient identities, the leading
discriminant term, both nonsingular-system substitutions, every
singular-system consistency factor, the collision sum-of-cubes identity,
and the complete rank-zero rational-point input. These are exact
symbolic and elliptic arithmetic checks. The written reduction from a
multiple discriminant root to the listed systems, its completeness,
and the genus argument remain proof steps for mathematical review.
The proof and replay are internal; historical priority, independent
replay, formal verification and external review are unestablished.
