# The factorially weighted multitorus theorem

## 1. Status and provenance

This note imports, credits, and audits the main theorem of Christopher D.
Long,
[*A Factorially Weighted Constant-Term Theorem on Algebraic
Tori*](https://github.com/octonion/mathematics/blob/main/gmc/gmc2_stronger_arbitrary_torus.tex),
dated 22 July 2026.  The
[public repository announcement](https://github.com/octonion/mathematics/commit/2f1ba8ae5d3d282141aca2285326a2d31ce20cd2)
predates the first checked-in GMC(2) formalization in this repository.
This records public provenance; it does not adjudicate private discovery
chronology.

Long's manuscript states that it was produced with AI assistance and had not
yet been human verified or formalized.  The argument below is a local
mathematical audit and reformulation.  It is not peer review and does not
transfer authorship.  The external repository is
[MIT licensed](https://github.com/octonion/mathematics/blob/main/LICENSE).

The rank-one Gaussian specialization overlaps the independently maintained
[GMC(2) lower-face theorem](TWO_REAL_GMC_LOWER_FACE_THEOREM.md) and its
[Lean package](../formal/gmc2/README.md).  The genuinely new imported scope
for this repository is:

1. arbitrary angular torus rank with a **single** factorial radial variable;
2. classification from eventual, rather than all-order, pure vanishing;
3. nonvanishing at \(kp\) for every sufficiently large prime \(p\);
4. the resulting Mathieu--Zhao theorem for the whole kernel.

The same geometry also combines with the repository's
[prime-separating radial theorem](PRIME_SEPARATING_RADIAL_MOMENTS.md).
Section 5 records the resulting arbitrary-torus extension.

## 2. Long's theorem

For \(r\geq1\), put

\[
 \mathcal A_r
 =\mathbb C[u,z_1^{\pm1},\ldots,z_r^{\pm1}]
\]

and define

\[
 \Gamma_r(u^nz^\mu)
 =n!\,\mathbf 1_{\mu=0},
 \qquad n\in\mathbb N,\quad \mu\in\mathbb Z^r.
 \tag{2.1}
\]

If

\[
 F=\sum_{n,\mu}c_{n,\mu}u^nz^\mu,
\]

write

\[
 \Omega(F)
 =\{\mu:c_{n,\mu}\neq0\text{ for some }n\}.
 \tag{2.2}
\]

> **Theorem 2.1 (Long's factorially weighted multitorus theorem).**
> For \(F\in\mathcal A_r\), the following are equivalent:
>
> 1. \(\Gamma_r(F^m)=0\) for every sufficiently large \(m\);
> 2. \(\Gamma_r(F^m)=0\) for every \(m\geq1\);
> 3. \(0\notin\operatorname{conv}\Omega(F)\).
>
> Under these conditions, for every fixed \(G\in\mathcal A_r\),
> \[
>  \Gamma_r(GF^m)=0\qquad(m\gg0).
>  \tag{2.3}
> \]
> Hence \(\ker\Gamma_r\) is a Mathieu--Zhao subspace.

The quantitative form is stronger than the failure of eventual vanishing.

> **Theorem 2.2 (Long's prime-index nonvanishing).**
> If \(F\neq0\) and
> \(0\in\operatorname{conv}\Omega(F)\), then there is \(k\geq1\) such that
> \[
>  \Gamma_r(F^{kp})\neq0
>  \tag{2.4}
> \]
> for every sufficiently large rational prime \(p\).

These are one-radial theorems.  The number of angular weights is arbitrary,
but the functional sees only the single radial exponent \(n\).

## 3. Exposed-face proof

Write the support of \(F\) as a finite subset

\[
 S\subset\mathbb N\times\mathbb Z^r.
\]

Assume \(0\in\operatorname{conv}\Omega(F)\).  Minimize the average radial
degree

\[
 \rho
 =\min\left\{
   \sum_{(n,\mu)\in S}t_{n,\mu}n:
   t_{n,\mu}\geq0,
   \sum t_{n,\mu}=1,
   \sum t_{n,\mu}\mu=0
 \right\}.
 \tag{3.1}
\]

Linear-programming duality gives \(\nu\in\mathbb R^r\) such that

\[
 n+\langle\nu,\mu\rangle\geq\rho
 \qquad((n,\mu)\in S).
 \tag{3.2}
\]

Let \(S_0\) be the equality set.  Complementary slackness shows that the
origin lies in the convex hull of the weights occurring in \(S_0\).
Projection \(S_0\to\mathbb Z^r\) is injective: if two face terms have the
same \(\mu\), equation (3.2) forces their radial exponents to agree.
Consequently the face polynomial

\[
 h(z)=\sum_{(n,\mu)\in S_0}c_{n,\mu}z^\mu
 \tag{3.3}
\]

has exactly the projected support and
\(0\in\operatorname{conv}\operatorname{supp}(h)\).

By the multivariable theorem of Duistermaat and van der Kallen, some
\(k\geq1\) satisfies

\[
 C:=\operatorname{CT}_z(h^k)\neq0.
 \tag{3.4}
\]

Every weight-zero product of \(k\) face terms has radial exponent

\[
 d=k\rho\in\mathbb N.
\]

Terms using a factor outside the face have strictly larger
\((n,\mu)\mapsto n+\langle\nu,\mu\rangle\)-degree.  Thus, for
\(R=F^k\),

\[
 [u^dz^0]R=C\neq0,
 \tag{3.5}
\]

and every monomial \(u^Nz^\eta\) of \(R\) satisfies

\[
 N+\langle\nu,\eta\rangle\geq d.
 \tag{3.6}
\]

In particular every weight-zero monomial of \(R^p\) has radial degree at
least \(pd\).  Normalize its moment in characteristic zero:

\[
 \frac{\Gamma_r(R^p)}{(pd)!}
 =\sum_{j\geq pd}[u^jz^0]R^p\,\frac{j!}{(pd)!}.
 \tag{3.7}
\]

Adjoin \(C^{-1}\) to the finitely generated coefficient domain.  Finite-type
specialization gives a characteristic-\(p\) residue field in which \(C\)
remains nonzero for every sufficiently large prime \(p\).  Frobenius gives

\[
 \overline{R^p}
 =\sum_{n,\mu}\bar r_{n,\mu}^{\,p}u^{pn}z^{p\mu}.
 \tag{3.8}
\]

After reducing (3.7), only original weight-zero terms remain.  The
contribution at \(n=d\) is \(\bar C^{\,p}\).  For \(n>d\),

\[
 \frac{(pn)!}{(pd)!}
\]

contains \(p(d+1)\) and vanishes modulo \(p\).  Hence

\[
 \frac{\Gamma_r(R^p)}{(pd)!}
 \equiv\bar C^{\,p}\neq0\pmod p.
 \tag{3.9}
\]

This proves Theorem 2.2.  It contradicts eventual pure vanishing and proves
the nontrivial direction of Theorem 2.1.

Conversely, if the origin is outside the convex hull, strict separation gives
a linear functional positive on every weight of \(F\).  Therefore \(F^m\)
has no weight-zero term, and the bounded weights of a fixed \(G\) cannot
cancel the linearly escaping weights of \(F^m\).  This proves both all-order
pure vanishing and eventual mixed vanishing.

## 4. Gaussian and nullcone specialization

For independent standard real Gaussians \(X,Y\), put

\[
 Z=\frac{X+iY}{\sqrt2},
 \qquad
 W=\frac{X-iY}{\sqrt2}.
\]

Then

\[
 \mathbb E(Z^aW^b)=a!\,\mathbf1_{a=b}.
\]

The injective algebra map

\[
 \iota:\mathbb C[Z,W]\hookrightarrow\mathcal A_1,
 \qquad
 \iota(Z)=z,\quad
 \iota(W)=uz^{-1},
 \tag{4.1}
\]

intertwines Gaussian expectation and \(\Gamma_1\).  Theorem 2.1 therefore
recovers GMC(2), while Theorem 2.2 strengthens the local statement:
two-sided angular support is detected at every sufficiently large prime
dilation of one fixed power.

If all angular weights are positive, the one-parameter subgroup

\[
 (Z,W)\longmapsto(tZ,t^{-1}W)
\]

drives the polynomial to zero; if they are negative, use its inverse.  In
the original \(X,Y\) coordinates this is a subgroup of
\(\operatorname{SO}_2(\mathbb C)\).  Thus pure Gaussian moment-nullness in
two variables is equivalent to membership in the orthogonal nullcone.
This is a corollary of the already proved one-sided-support classification,
not a new Gaussian dimension bound.

## 5. Repository synthesis: prime-separating multitorus functionals

Let \(a=(a_n)_{n\geq0}\) be a sequence of nonzero rational numbers and define

\[
 \Phi_{a,r}(u^nz^\mu)=a_n\,\mathbf1_{\mu=0}.
 \tag{5.1}
\]

Recall the
[normalized prime-separation conditions](PRIME_SEPARATING_RADIAL_MOMENTS.md#3-the-correct-arithmetic-hypothesis):
for every \(0\leq n<M\), arbitrarily large primes \(p\) satisfy

\[
 \frac{a_d}{a_{np}}\in\mathbb Z_{(p)}
 \quad(np\leq d\leq Mp)
 \tag{5.2}
\]

and

\[
 \frac{a_{jp}}{a_{np}}\in p\mathbb Z_{(p)}
 \quad(n<j\leq M).
 \tag{5.3}
\]

> **Theorem 5.1 (prime-separating multitorus theorem).**
> If \(a\) is prime-separating and
> \[
>  \Phi_{a,r}(F^m)=0\qquad(m\gg0),
> \]
> then \(0\notin\operatorname{conv}\Omega(F)\).  Consequently
> \(\ker\Phi_{a,r}\) is Mathieu--Zhao for every torus rank \(r\).

The geometric extraction in Section 3 is unchanged.  For the resulting
\(R=F^k\), choose \(M>\deg_uR\), a prime that is simultaneously good for
the coefficient ring and satisfies (5.2)--(5.3), and normalize
\(\Phi_{a,r}(R^p)\) by \(a_{pd}\).  Condition (5.2) makes every term
integral before reduction.  Frobenius kills coefficients whose radial
degree is not divisible by \(p\).  Condition (5.3) kills every surviving
degree \(pj\) with \(j>d\), leaving \(C^p\neq0\).  This contradicts eventual
vanishing.  Strict torus separation again gives the mixed conclusion.

For \(a_n=n!\), the admissible primes are all sufficiently large primes and
the argument recovers Long's stronger Theorem 2.2.  Theorem 5.1 is the
combination of Long's multitorus exposed-face geometry with the repository's
prime-separating radial arithmetic.

## 6. Search and implication-chain consequences

The theorem is a complete sieve for one-radial torus models.

| Candidate model | Consequence |
|---|---|
| Exact functional-intertwining embedding into some \((\mathcal A_r,\Gamma_r)\) | Cannot yield a Mathieu--Zhao counterexample |
| \(0\notin\operatorname{conv}\Omega(F)\) | Pure moments vanish for support reasons and every fixed mixed moment vanishes eventually |
| \(0\in\operatorname{conv}\Omega(F)\) | Pure moments do not vanish eventually; one fixed power has nonzero moments at all sufficiently large prime dilations |
| Only a long finite zero prefix is known | Do not promote it; the prime-index theorem supplies an infinite obstruction subsequence |

Accordingly, a search should test for an exact one-radial functional model
before computing moment prefixes.  If such a model exists, no coefficient
search inside it is needed.  A surviving counterexample architecture must
break at least one hypothesis: it must retain several genuinely independent
radial directions, use a non-prime-separating radial functional, or fail to
intertwine multiplication and the moment functional.

This gives no new arrow from a multitorus theorem to higher-dimensional GMC,
SIC, or the Jacobian Conjecture.  In particular:

- ordinary \(2q\)-real Gaussian expectation has \(q\) independent radial
  variables and weight
  \(\prod_i n_i!\), not one \(n!\);
- the two-pair Image-Mathieu contraction functional is multifactorial;
- interval, beta, and compact-Haar functionals need not be
  prime-separating and already support counterexamples.

The known higher-dimensional witnesses exploit precisely these excluded
multi-radial or non-factorial mechanisms.  Torus rank is therefore not the
obstruction to the lower-face method; independent radial rank is.

## 7. Reproduction boundary

The dependency-free checker
[`audit_factorially_weighted_multitorus.py`](../scripts/audit_factorially_weighted_multitorus.py)
replays finite rank-two examples, the normalized prime congruence, strict
support separation, and the Gaussian embedding.  These are bounded exact
regressions, not the proof.  The proof is Sections 3 and 5, with two named
external inputs:

1. the multivariable Duistermaat--van der Kallen constant-term theorem;
2. good finite-characteristic specialization of finite-type integral
   domains.

The checked-in Lean package formalizes the rank-one Gaussian specialization
only.  No arbitrary-rank formalization is claimed.
