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

Let

\[
 \lambda(w)=w_t(w_aw_d-w_bw_c),\qquad
 f(w,z)=\lambda(w)P(z),
\]

where \(z=(t,a,b,c,d)\), and let

\[
 \mathcal E(w^\alpha q(z))=\partial_z^\alpha q(z).
\]

Then, term by term,

\[
 \mathcal E(f^m)=\Lambda^m(P^m),\qquad
 \mathcal E(Qf^m)=\Lambda^m(QP^m).                           \tag{3.1}
\]

Thus every positive power of \(f\) lies in \(\ker\mathcal E\), while
\(Qf^m\notin\ker\mathcal E\) for every \(m\geq2\).  Hence
\(\ker\mathcal E\) is not a Mathieu--Zhao space already for five contraction
pairs:

\[
 \boxed{\neg\operatorname{SIC}(5).}                          \tag{3.2}
\]

The ambient SIC polynomial ring has ten variables, five \(w\)'s and five
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
 2\leq r_{\rm SIC}\leq5,\qquad
 2\leq n_{\Delta{\rm GVC}}\leq40,\qquad
 6\leq n_{\rm HN,4}\leq42.                                  \tag{5.1}
\]

The lower endpoints use GVC(1) and SIC(1).  The last two upper endpoints
remain the repository's ordinary-Laplacian and homogeneous quartic HN
witnesses.
