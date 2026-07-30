# The Dvorsky--Long five-variable GVC and SIC counterexample

## 1. Statement and scope

Over a field of characteristic zero, put

\[
 P=(t+c)(ad+bt),\qquad
 \Lambda=\partial_t(\partial_a\partial_d-\partial_b\partial_c),
 \qquad Q=-c.                                                 \tag{1.1}
\]

Then

\[
 \boxed{\Lambda^m(P^m)=0\quad(m\geq1)}                       \tag{1.2}
\]

but

\[
 \boxed{\Lambda^m(QP^m)
 =(-1)^{m+1}(m+1)!(m!)^2t\ne0\quad(m\geq2).}                 \tag{1.3}
\]

Consequently the unrestricted constant-coefficient Generalized Vanishing
Conjecture fails in five variables.  The same formulas give a counterexample
to the Special Image Conjecture in five contraction pairs.

This does **not** lower the ordinary-Laplacian or homogeneous quartic
Hessian-nilpotent witness dimensions.  The displayed operator has order
three, whereas the ordinary Laplacian has order two.

## 2. Binomial-coefficient proof

Write

\[
 D=\partial_a\partial_d-\partial_b\partial_c,
 \qquad \Lambda^m=\partial_t^mD^m.                            \tag{2.1}
\]

The two factors in (2.1) commute.  Expand

\[
 P^m=\sum_{i,k=0}^m
 \binom mi\binom mk
 c^it^{m-i+k}a^{m-k}d^{m-k}b^k                              \tag{2.2}
\]

and

\[
 D^m=\sum_{j=0}^m(-1)^j\binom mj
 (\partial_a\partial_d)^{m-j}(\partial_b\partial_c)^j.       \tag{2.3}
\]

In (2.2)--(2.3), the \(a,d\) derivatives require \(k\leq j\), while
the \(b\) derivatives require \(k\geq j\).  Thus only \(k=j\) survives.
The factorials simplify to

\[
 \binom mj^2((m-j)!)^2j!\,(i)_j
 =(m!)^2\binom ij.
\]

It follows that

\[
\begin{aligned}
 D^m(P^m)
 &=(m!)^2\sum_{i=0}^m\binom mi t^{m-i}
       \sum_{j=0}^i(-1)^j\binom ijc^{i-j}t^j\\
 &=(m!)^2\sum_{i=0}^m\binom mi t^{m-i}(c-t)^i\\
 &=(m!)^2c^m.                                                \tag{2.4}
\end{aligned}
\]

Applying \(\partial_t^m\) proves (1.2).

For \(QP^m=-cP^m\), the same calculation replaces
\(\binom ijc^{i-j}\) by \(\binom{i+1}{j}c^{i+1-j}\).  The full binomial sum
would include \(j=i+1\).  The operator sum (2.3) omits only the boundary term
\((i,j)=(m,m+1)\).  Therefore

\[
 D^m(-cP^m)
 =(m!)^2\left((-1)^{m+1}t^{m+1}-c^m(c-t)\right).             \tag{2.5}
\]

For \(m\geq2\), \(\partial_t^m\) kills the second term of (2.5), while

\[
 \partial_t^m(t^{m+1})=(m+1)!\,t.
\]

This proves (1.3).  At \(m=1\), the value is \(c+2t\), which is also
nonzero.

The quantifiers now match the standard GVC definition exactly: the
hypothesis holds for every \(m\), but for the single fixed multiplier
\(Q=-c\), the asserted eventual vanishing fails at every \(m\geq2\).

## 3. The five-pair SIC consequence

Write

\[
 z=(t,a,b,c,d),\qquad
 \zeta=(\tau,\alpha,\beta,\gamma,\delta),
\]

and consider exactly the Image-Conjecture subspace

\[
 \mathcal M_5
 =\sum_{x\in\{t,a,b,c,d\}}(\partial_x-\zeta_x)
   \mathbb C[\zeta,z].                                      \tag{3.1}
\]

The following is a compact explicit Mathieu-subspace counterexample:

\[
\boxed{
\begin{aligned}
 f&=\tau(\alpha\delta-\beta\gamma)(t+c)(ad+bt),\\
 g&=-c.
\end{aligned}}                                               \tag{3.2}
\]

In particular, \(f\) has degree six and eight terms when expanded, all with
coefficients \(\pm1\), while \(g\) is a linear monomial.

To verify the claim, let

\[
 \mathcal E(\zeta^\mu q(z))=\partial_z^\mu q(z).
\]

Zhao's image-kernel identity gives

\[
 \mathcal M_5=\ker\mathcal E,                                \tag{3.3}
\]

since replacing every generator \(\zeta_x-\partial_x\) by its negative does
not change its image.  Equivalently, if

\[
 h=\sum_{\mu,\nu}c_{\mu,\nu}\zeta^\mu z^\nu,
\]

then

\[
 h\in\mathcal M_5
 \quad\Longleftrightarrow\quad
 \sum_{\mu}(\mu+\rho)!\,c_{\mu,\mu+\rho}=0
 \quad\text{for every }\rho\in\mathbb N^5.                  \tag{3.4}
\]

Thus membership is a finite list of exact coefficient identities for each
fixed polynomial.

Now let

\[
 \lambda(\zeta)=\tau(\alpha\delta-\beta\gamma),\qquad
 P(z)=(t+c)(ad+bt).
\]

Then, term by term,

\[
 \mathcal E(f^m)=\Lambda^m(P^m),\qquad
 \mathcal E(gf^m)=\Lambda^m(gP^m).                           \tag{3.5}
\]

Thus every positive power of \(f\) lies in \(\ker\mathcal E\), while
\(gf^m\notin\ker\mathcal E\) for every \(m\geq1\).  Explicitly, the value
for \(m=1\) is \(c+2t\), and for every \(m\geq2\) it is the nonzero polynomial
in (1.3).  Hence
\(\ker\mathcal E\) is not a Mathieu--Zhao space already for five contraction
pairs:

\[
 \boxed{
 f^m\in\mathcal M_5\ (m\geq1),\qquad
 gf^m\notin\mathcal M_5\ (m\geq1).
 }                                                           \tag{3.6}
\]

The ambient SIC polynomial ring has ten variables, five \(\zeta\)'s and five
\(z\)'s.  This pair count must not be confused with the five-variable count
in GVC(5).

## 4. Provenance

Christopher D. Long's
[*Counterexamples to the xz-Conjecture and the Mathieu Conjecture for
SU(2)*](https://arxiv.org/abs/2607.19012) gives the dehomogenized \(SU(2)\)
seed

\[
 F=(1+c)(ad+b),\qquad G=-c.
\]

Alexander Dvorsky posted (1.1) and the resulting GVC(5) claim on 23 July 2026
in the [Secret Blogging Seminar
discussion](https://sbseminar.wordpress.com/2026/07/20/the-new-counterexample-to-the-jacobian-conjecture/),
and subsequently recorded the SIC lift.  The appropriate attribution is:

> Long's \(SU(2)\) seed; Dvorsky's homogenization and five-variable GVC/SIC
> lift.

Dvorsky also clarified in that discussion that Long's \(SU(2)\) example was
motivated by, but does not follow algebraically from, the announced
three-dimensional Jacobian counterexample.  No such implication is claimed
here.

## 5. Reproduction and scoreboard effect

Run

```bash
python3 scripts/audit_dvorsky_gvc5_counterexample.py
```

The dependency-free checker uses exact sparse integer polynomials and
verifies (2.4), (2.5), (1.2), and (1.3) through \(m=8\).  The displayed
binomial calculation, not the bounded replay, proves the all-order result.

The resulting ambient-dimension ledger is

\[
 2\leq n_{\rm GVC}\leq5,\qquad
 r_{\rm SIC}=2,\qquad
 2\leq n_{\Delta{\rm GVC}}\leq40,\qquad
 6\leq n_{\rm HN,4}\leq40.                                  \tag{5.1}
\]

The sharp SIC endpoint now comes from the separate
[`TWO_PAIR_IMAGE_MATHIEU_COUNTEREXAMPLE.md`](TWO_PAIR_IMAGE_MATHIEU_COUNTEREXAMPLE.md),
together with SIC(1).  The GVC lower endpoint uses GVC(1).  The last two
upper endpoints come respectively from the repository's
ordinary-Laplacian witness and the externally certified MacFarlane
homogeneous-cotangent quartic HN witness.

## 6. A bounded four-variable descent search

The first finite descent attempted here keeps the architecture of Long's
\(SU(2)\) seed and Dvorsky's homogenization visible.  Identify \(c\) with a
linear form so that, up to an irrelevant scalar sign,

\[
 P=(rt+ua+vb+wd)(ad+bt),                                   \tag{6.1}
\]

and replace the quadratic part of the symbol by the general ternary
quadratic form

\[
\begin{aligned}
 \Lambda&=\partial_tR(\partial_a,\partial_b,\partial_d),\\
 R&=A\partial_a^2+B\partial_a\partial_b+C\partial_a\partial_d
   +D\partial_b^2+E\partial_b\partial_d+F\partial_d^2.       \tag{6.2}
\end{aligned}
\]

Thus (6.2), rather than a selected diagonal or binomial subfamily, is the
full order-three symbol in the symmetry-preserving slice
\(\partial_t\operatorname{Sym}^2\langle\partial_a,\partial_b,\partial_d\rangle\).
The exact lattice search uses

\[
 (r,u,v,w)\in\{-1,0,1\}^4,\qquad
 (A,B,C,D,E,F)\in\{-2,-1,0,1,2\}^6,                        \tag{6.3}
\]

with sign normalization on both vectors and primitive normalization on the
second.  It checks all \(297{,}920\) resulting pairs.  There are \(7{,}152\)
solutions of the initial scalar scheme

\[
 \Lambda^m(P^m)=0,\qquad 1\leq m\leq4.
\]

Only \(7{,}120\) survive through \(m=12\): all \(32\) delayed failures occur
first at \(m=6\).  For example,

\[
\begin{aligned}
 P&=(t-a-b+d)(ad+bt),\\
 \Lambda&=\partial_t(
 \partial_a\partial_b-\partial_b^2
 -\partial_b\partial_d-\partial_d^2)
\end{aligned}
\]

has zero pure contractions for \(1\leq m\leq5\), but

\[
 \Lambda^6(P^6)=-22{,}394{,}880{,}000.
\]

This is a concrete warning that the requested \(m\leq4\) scheme has
resonant false positives.

For every one of the \(7{,}120\) order-twelve survivors, the search also
computes the full linear map

\[
 Q\longmapsto\Lambda^m(QP^m),\qquad
 Q\in\langle t,a,b,d\rangle,
\]

for \(5\leq m\leq12\).  No survivor admits a fixed linear multiplier that
is nonzero at every order in that window.  This uses the whole multiplier
space: over characteristic zero a fixed generic \(Q\) avoids finitely many
nonzero kernels, so it is enough to test whether the map itself is nonzero
at each order.

Run

```bash
python3 scripts/search_dvorsky_gvc4_bounded.py
```

The generated certificate is
[`dvorsky_gvc4_bounded_search.json`](../artifacts/generated-results/dvorsky_gvc4_bounded_search.json).
This is an exhaustive negative result only in the declared lattice slice.
It does not exclude rational points outside the coefficient boxes, symbols
outside (6.2), nonlinear multipliers, or later behavior of the cutoff
survivors.  In particular it does **not** establish GVC(4) or SIC(4), and
the unrestricted GVC frontier remains the five-variable Long--Dvorsky
witness.  The overall SIC frontier has independently dropped to three pairs.
