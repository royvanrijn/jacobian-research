# Factorial-moment translation and sharp finite witnesses

## 1. Result and scope

Let

\[
 \mathcal L(U^\alpha)=\alpha!
 =\prod_{i=1}^n\alpha_i!
 \qquad
 \left(\mathcal L:\mathbb C[U_1,\ldots,U_n]\longrightarrow\mathbb C\right).
                                                               \tag{1.1}
\]

The Factorial Conjecture asks whether

\[
 \mathcal L(f^m)=0\quad\text{for every }m\geq1
 \qquad\Longrightarrow\qquad f=0.                              \tag{1.2}
\]

Van den Essen--Wright--Zhao introduced this conjecture while deriving it
from the Image Conjecture.  They proved the one-variable case and several
support-restricted cases in
[*On the Image Conjecture*](https://arxiv.org/abs/1008.3962).
Edo--van den Essen later proposed the
[*Strong Factorial Conjecture*](https://arxiv.org/abs/1304.3956): a nonzero
polynomial should not have \(\mathcal N(f)\) consecutive zero factorial
moments, where \(\mathcal N(f)\) is its number of nonzero monomials.

The experiment below gives three exact conclusions.

1. The Dvorsky--Long GVC(5) witness has an all-order factorial-functional
   translation, but the translated polynomial depends on \(m\).  Diagonal
   extraction is not multiplicative, already at \(m=2\), so this does not
   give a polynomial \(f\) satisfying (1.2).
2. For every \(r\geq2\), an \(r\)-term homogeneous linear form has exactly
   its first \(r-1\) factorial moments zero.  This is the longest possible
   initial zero string among homogeneous linear forms with \(r\) nonzero
   terms.
3. For every odd \(r\geq3\), tensoring that cyclotomic filter with the
   Dvorsky diagonal shadow gives a \(2r\)-term quartic with exactly its
   first \(2r-1\) factorial moments zero.

Thus there is **no counterexample to the Factorial Conjecture here**:
every displayed fixed polynomial has a later nonzero moment.  There are,
however, two sharp finite negative results:

- no moment cutoff depending only on total degree, uniformly in the ambient
  dimension, can imply \(f=0\), even in degree one;
- the monomial-count threshold in the Strong Factorial Conjecture cannot be
  decreased from \(\mathcal N(f)\) consecutive moments to
  \(\mathcal N(f)-1\).

The calculations are reproduced by
[`verify_factorial_moment_witnesses.py`](../scripts/verify_factorial_moment_witnesses.py).

## 2. The factorial functional as a moment

The elementary integral representation

\[
 \mathcal L(f)=
 \int_{\mathbb R_{\geq0}^n}
 f(u_1,\ldots,u_n)e^{-(u_1+\cdots+u_n)}\,du_1\cdots du_n       \tag{2.1}
\]

identifies \(\mathcal L\) with expectation under independent
mean-one exponential variables.  It also identifies the factorial sector
inside circular Gaussian moments: if
\(\mathbb E(W_i^{a_i}Z_i^{b_i})=\delta_{a_i,b_i}a_i!\), then

\[
 \mathcal L(f(U_1,\ldots,U_n))
 =\mathbb E\bigl[f(W_1Z_1,\ldots,W_nZ_n)\bigr].                \tag{2.2}
\]

This explains both the promise and the limitation of importing Gaussian
witnesses.  A torus-invariant Gaussian polynomial is already a factorial
polynomial after \(U_i=W_iZ_i\).  Long's five-term witness is not
torus-invariant and also uses an unpaired real Gaussian source.  Angular
constant-term extraction must therefore be applied separately to every
power.  Since constant-term extraction is not multiplicative, the resulting
sequence need not consist of powers of one factorial polynomial.

The same obstruction can be computed exactly, without the unpaired source,
from the Dvorsky--Long GVC(5) witness.

## 3. Exact diagonal translation of the GVC(5) witness

Use the variable order \((t,a,b,c,d)\) and put

\[
\begin{aligned}
 \lambda(w)&=w_t(w_aw_d-w_bw_c),\\
 P(z)&=(t+c)(ad+bt),\\
 p(w,z)&=\lambda(w)P(z).
\end{aligned}                                                  \tag{3.1}
\]

The all-order identity in
[`DVORSKY_GVC5_COUNTEREXAMPLE.md`](DVORSKY_GVC5_COUNTEREXAMPLE.md) is

\[
 \lambda(\partial_z)^mP(z)^m=0\qquad(m\geq1).                  \tag{3.2}
\]

Define torus-diagonal extraction by

\[
 \operatorname{Diag}(w^\alpha z^\beta)
 =\begin{cases}
 U^\alpha,&\alpha=\beta,\\
 0,&\alpha\ne\beta.
 \end{cases}                                                   \tag{3.3}
\]

Both \(\lambda\) and \(P\) are homogeneous of degree three.  Hence
componentwise derivative survival in (3.2), together with equality of total
degrees, forces \(\alpha=\beta\).  Termwise contraction gives

\[
 \lambda(\partial_z)^mP(z)^m
 =\mathcal L\!\left(\operatorname{Diag}(p^m)\right).           \tag{3.4}
\]

Write \(X=U_aU_d\) and \(Y=U_bU_c\).  Direct coefficient matching gives

\[
 h_m:=\operatorname{Diag}(p^m)
 =U_t^m\sum_{k=0}^m(-1)^k\binom{m}{k}^3X^{m-k}Y^k.             \tag{3.5}
\]

Indeed, choosing \(k\) negative symbol terms, \(k\) copies of \(bct\), and
the matching terms from the two binomial factors gives the three binomial
coefficients.  Applying (1.1) yields

\[
\begin{aligned}
 \mathcal L(h_m)
 &=m!\sum_{k=0}^m(-1)^k\binom{m}{k}^3
       ((m-k)!)^2(k!)^2\\
 &=(m!)^3\sum_{k=0}^m(-1)^k\binom{m}{k}=0.                    \tag{3.6}
\end{aligned}
\]

Thus (3.5) is an all-order sequence of nonzero factorial-null polynomials.
It is not the power sequence of one polynomial.  In fact

\[
\begin{aligned}
 h_1&=U_t(X-Y),\\
 h_2&=U_t^2(X^2-8XY+Y^2),\\
 h_1^2&=U_t^2(X^2-2XY+Y^2),
\end{aligned}                                                  \tag{3.7}
\]

and

\[
 \mathcal L(h_2)=0,\qquad \mathcal L(h_1^2)=12.                \tag{3.8}
\]

Equation (3.8) is the first exact obstruction to promoting the Image/GVC
witness to a Factorial Conjecture witness: the operation producing \(h_m\)
does not commute with taking powers.

## 4. Smallest homogeneous linear witnesses for a prescribed prefix

Let

\[
 g=c_1X_1+\cdots+c_sX_s,\qquad c_1\cdots c_s\ne0.              \tag{4.1}
\]

Expansion of \(g^m\) cancels the multinomial denominators against the
factorials in (1.1):

\[
 \frac{\mathcal L(g^m)}{m!}
 =\sum_{\alpha_1+\cdots+\alpha_s=m}c_1^{\alpha_1}\cdots
 c_s^{\alpha_s}
 =h_m(c_1,\ldots,c_s),                                        \tag{4.2}
\]

where \(h_m\) is the complete homogeneous symmetric polynomial.  Therefore

\[
 \sum_{m\geq0}\frac{\mathcal L(g^m)}{m!}T^m
 =\prod_{j=1}^s(1-c_jT)^{-1}.                                 \tag{4.3}
\]

Suppose the first \(R\) moments vanish.  If \(R\geq s\), comparison in

\[
 \left(\sum_{m\geq0}h_mT^m\right)
 \left(\sum_{j=0}^s(-1)^je_jT^j\right)=1                      \tag{4.4}
\]

successively forces \(e_1=\cdots=e_s=0\).  But
\(e_s=c_1\cdots c_s\ne0\), a contradiction.  Consequently

\[
 R\leq s-1.                                                    \tag{4.5}
\]

This bound is attained.  Let \(\zeta_r\) be a primitive \(r\)-th root of
unity and set

\[
 G_r=\sum_{j=0}^{r-1}\zeta_r^jX_j.                             \tag{4.6}
\]

Since the coefficient multiset is the complete set of \(r\)-th roots,

\[
 \prod_{j=0}^{r-1}(1-\zeta_r^jT)=1-T^r.
\]

Equations (4.3) and (4.6) give the all-order formula

\[
 \boxed{
 \mathcal L(G_r^m)=
 \begin{cases}
 m!,&r\mid m,\\
 0,&r\nmid m.
 \end{cases}}                                                  \tag{4.7}
\]

For a requested initial zero range \(1\leq m\leq R\), the least possible
number of nonzero terms among homogeneous linear forms is therefore \(R+1\),
and \(G_{R+1}\) attains it.  The first nontrivial case is

\[
 G_3=X_0+\omega X_1+\omega^2X_2,\qquad
 \mathcal L(G_3)=\mathcal L(G_3^2)=0,\quad
 \mathcal L(G_3^3)=6.                                         \tag{4.8}
\]

Because \(r\) is arbitrary while every \(G_r\) has degree one, (4.7)
excludes a finite cutoff depending only on polynomial degree and not on
ambient dimension.

## 5. A witness-derived quartic saturating the strong threshold

The first Dvorsky diagonal shadow in (3.7) has moments

\[
\begin{aligned}
 A&=U_t(U_aU_d-U_bU_c),\\
 \mathcal L(A^m)
 &=(m!)^3\sum_{k=0}^m\frac{(-1)^k}{\binom{m}{k}}\\
 &=(m!)^3(1+(-1)^m)\frac{m+1}{m+2}.                           \tag{5.1}
\end{aligned}
\]

The last equality is the elementary alternating reciprocal-binomial
identity.  In particular, precisely the odd moments vanish.

Take the variables of \(G_r\) disjoint from those of \(A\), let \(r\) be
odd, and define

\[
 F_r=A\,G_r.                                                   \tag{5.2}
\]

The factorial functional factors across disjoint variable sets.  Combining
(4.7) and (5.1), a moment of \(F_r\) can be nonzero only when its order is
both even and divisible by \(r\).  Since \(r\) is odd, the first such order
is \(2r\):

\[
\boxed{
 \begin{aligned}
  \mathcal L(F_r^m)&=0 &&(1\leq m<2r),\\
  \mathcal L(F_r^{2r})
    &=((2r)!)^4\frac{2r+1}{r+1}\ne0.
 \end{aligned}}                                                \tag{5.3}
\]

The polynomial \(F_r\) has degree four and exactly \(2r\) monomials.  Thus
its first \(\mathcal N(F_r)-1\) moments vanish.  The case \(r=3\) is a
six-term quartic in eight variables with

\[
 \mathcal L(F_3^m)=0\ (1\leq m\leq5),\qquad
 \mathcal L(F_3^6)=470\,292\,480\,000.                         \tag{5.4}
\]

This does not violate the Strong Factorial Conjecture; it proves that its
proposed run length is sharp even inside a family directly descended from
the diagonal shadow of the GVC(5) witness.

## 6. Reproduction and remaining search

Run

```bash
python3 scripts/verify_factorial_moment_witnesses.py
```

The dependency-free checker:

- expands \(p^m\), extracts its diagonal, and verifies (3.5)--(3.6) through
  \(m=8\);
- verifies the nonmultiplicativity certificate (3.7)--(3.8);
- reconstructs the cyclotomic moments in exact prime cyclotomic rings for
  \(r=2,3,5,7,11\), through order \(2r\);
- verifies the first nonzero moments of the quartics \(F_r\) for
  \(r=3,5,7,11\);
- writes
  [`factorial_moment_witnesses.json`](../artifacts/generated-results/factorial_moment_witnesses.json).

The all-order proofs are the displayed binomial and generating-series
identities, not the finite replay.

The remaining genuine target is a fixed polynomial with all positive
factorial moments zero.  The present translation identifies the precise
gate: replace the nonmultiplicative sequence
\(\operatorname{Diag}(p^m)\) by powers of one polynomial without destroying
the cancellations.  Finite prefixes alone cannot provide that conclusion,
and (4.7)--(5.3) show why degree-only bounded searches across growing
dimension can be arbitrarily misleading.
