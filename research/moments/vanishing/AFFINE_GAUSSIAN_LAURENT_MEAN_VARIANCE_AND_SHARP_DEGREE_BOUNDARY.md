# Affine Gaussian--Laurent mean--variance classification and sharp degree boundary

Let K be any characteristic-zero field, let n,r>=0, and put
R=K[t_1^{+-1},...,t_r^{+-1}] (R=K when r=0),
A=R[X_1,...,X_n]. Let gamma_n be the native Gaussian polynomial functional, and let Lambda=CT_t composed with gamma_n, with CT the identity if r=0. A phase is pure if every positive pure moment vanishes; it is safe if every fixed A-marker has eventually vanishing mixed moments.

Write an affine Gaussian phase as

\[
 P=a(t)+\sum_{i=1}^n b_i(t)X_i,\qquad
 V(t)=\sum_{i=1}^n b_i(t)^2,
 \qquad
 \Sigma(P)=\operatorname{supp}_t(a)\cup\operatorname{supp}_t(V).
 \tag{1}
\]

Supports mean the actual Laurent supports after collecting coefficients. In particular cancellations in V are respected. Put supp(0)=empty and conv(empty)=empty. The theorem does not require V to have a square root, to be nonzero, or to have a positivity interpretation.

**Theorem.** For this given affine phase over K, the following are equivalent:

1. Lambda(P^m)=0 for every m>=1.
2. Lambda(P^m)=0 for all sufficiently large m.
3. For every fixed h in A, Lambda(hP^m)=0 for all sufficiently large m.
4. Zero does not belong to conv Sigma(P).

If zero belongs to conv Sigma(P), then one integer e>=1 satisfies

\[
 \Lambda(P^{ep})\ne0
 \quad\hbox{for EVERY sufficiently large rational prime }p.
 \tag{2}
\]

Thus the failure of the separation condition has an explicit fixed-power prime-index witness; it is not an inference from a finite prefix of moments. Both equivalences and (2) hold over every characteristic-zero coefficient field and for r=0 as well as r>=1.

## 1. Exact angular hull of the paired reduction

Use one formal circular pair Z,W with C(Z^u W^v)=u! if u=v and zero otherwise, and set

\[
 L=a(t)+Z+\frac{V(t)}2W.
 \tag{3}
\]

The R-valued formal Gaussian identities give

\[
 \gamma_n(e^{sP})=\exp(sa+s^2V/2)=C(e^{sL}).
 \tag{4}
\]

Consequently Lambda(P^m)=(C tensor CT_t)(L^m) for every m>=0. This includes V=0, including nonzero isotropic coefficient vectors over fields permitting them.

The embedding Z->z_0, W->u z_0^{-1}, t_i->z_i sends L to

\[
 a(z)+z_0+\frac{V(z)}2u z_0^{-1}
 \quad\hbox{in }K[u,z_0^{+-1},...,z_r^{+-1}].
 \tag{5}
\]

It intertwines the circular/Laurent functional with the one-radial factorial functional Gamma_{r+1}. There is exactly one radial factorial variable, even when r grows. Its angular support is precisely

\[
 \Omega=\{(1,0)\}
 \cup\{(0,\mu):\mu\in\operatorname{supp}a\}
 \cup\{(-1,\mu):\mu\in\operatorname{supp}V\}.
 \tag{6}
\]

The three first-coordinate weights are distinct, so no cancellation between the three groups occurs.

We prove the exact geometric equivalence

\[
 0\in\operatorname{conv}\Omega
 \quad\Longleftrightarrow\quad
 0\in\operatorname{conv}\Sigma(P).
 \tag{7}
\]

For the forward implication, write a convex combination for zero with weight alpha on (1,0), weights beta_mu on the a-points and weights eta_mu on the V-points. First-coordinate balance says alpha=sum eta_mu. The sum of the beta and eta weights is positive: otherwise their nonnegativity would force both alpha=0 and total weight zero. Their last r coordinates sum to zero, so normalizing these weights gives a convex combination for zero in Sigma(P).

Conversely, take nonnegative weights beta_mu,eta_mu on tagged copies of the two supports, with total one and sum beta_mu mu + sum eta_mu mu=0. If an exponent belongs to both supports, assign its weight to either tagged copy. Set alpha=sum eta_mu. The weights alpha,beta_mu,eta_mu have total 1+alpha>0, and their weighted sum on (6) is zero in every coordinate. Dividing by 1+alpha proves the reverse implication. If Sigma(P) is empty, Omega is the singleton {(1,0)}; both origin-inclusion assertions in (7) are false. This proves (7) including every empty-support case.

## 2. All marker directions and prime nonvanishing

Long's one-radial factorial theorem, with arbitrary angular rank, says that all positive pure vanishing, eventual pure vanishing, angular Newton separation and all-marker safety are equivalent. Its prime-index refinement gives a fixed integer e>=1 with nonzero moments at every sufficiently large prime index ep when the angular hull contains zero. The source states the theorem over C; its given-phase support and exact-moment conclusions hold over any characteristic-zero field by embedding the finitely generated coefficient field into C. Strict separation then gives safety for every marker over the original field. No embedding of the whole field K is assumed.

Apply that theorem to (5), and use (4) and (7). This proves the equivalence of assertions 1, 2 and 4, and gives (2). It also shows that when these conditions hold, L is safe against every fixed paired-polynomial/Laurent marker. To transfer that safety to ALL original Gaussian markers, the following step is necessary; equality of pure moments alone would not suffice.

Expand any fixed h as a finite Hermite sum

\[
 h=\sum_\alpha q_\alpha(t)\operatorname{He}_\alpha(X),
 \qquad q_\alpha\in R.
\]

The multivariate Hermite generating identity gives, with d=|alpha| and m>=d,

\[
 \gamma_n(\operatorname{He}_\alpha P^m)
   =(m)_d\,b^\alpha\gamma_n(P^{m-d}).
 \tag{8}
\]

Here (m)_d=m!/(m-d)! is the FALLING factorial for m>=d, with (m)_0=1. For m<d the marked expectation is zero. This follows by extracting the v^alpha and s^m coefficients from

\[
 \gamma_n\left(e^{v\cdot X-\sum v_i^2/2}e^{sP}\right)
  =\exp(sa+s^2V/2+s\,b\cdot v).
\]

Combining (4) and (8), the alpha term of Lambda(hP^m) is

\[
 (m)_d\,(C\otimes\operatorname{CT}_t)
       \bigl(q_\alpha(t)b^\alpha L^{m-d}\bigr).
 \tag{9}
\]

Its marker q_alpha b^alpha is fixed, independent of m. Safety of L kills (9) at all sufficiently large shifted orders. The finite set of alpha gives one common threshold. Hence assertion 4 implies assertion 3. Conversely assertion 3 implies assertion 2 by the marker h=1. This closes every stated direction.

## 3. Laurent-rank-zero specialization and the sharp Gaussian-degree boundary

For r=0, every nonempty Laurent support is the singleton zero vector. The theorem therefore specializes to the exact all-field affine criterion

\[
 P\hbox{ pure}\quad\Longleftrightarrow\quad
 P\hbox{ eventually pure}\quad\Longleftrightarrow\quad
 P\hbox{ safe}\quad\Longleftrightarrow\quad
 a=0\ \hbox{and}\ \sum_i b_i^2=0.
 \tag{10}
\]

This is an affine-phase statement only. It is not an arbitrary-degree Gaussian classification over K. For example, a nonzero isotropic linear form may satisfy (10) over a non-formally-real field.

**Corollary (sharp dimension/Gaussian-degree boundary).** Fix a characteristic-zero field K and integers n,d>=0 and r>=1. The assertion that EVERY pure phase of Gaussian degree at most d is safe against ALL polynomial/Laurent markers holds if and only if n<=2 or d<=1. These bounded-degree phase classes are not asserted to be algebras.

For n<=2, the full native Gaussian--Laurent Mathieu theorem proves the assertion for every d. For d<=1, the mean--variance theorem above proves it for every n. Conversely, when n>=3 and d>=2, the assertion fails: the native rational quadratic phase

\[
 \frac{1+t_1}{t_1}
 [X_1^2+X_2^2-t_1(t_1+2)X_3^2]
\]

has every positive pure moment zero and, with the SAME fixed unit marker t_1, has marked value 2^m m! for EVERY m>=1. It ignores extra coordinates and belongs to each class d>=2 over the original K. Thus failure first occurs at Gaussian degree two when n>=3 and r>=1. Only Gaussian degree is bounded; Laurent exponents are unrestricted. No minimal ordinary total degree, priority claim, or r=0 arbitrary-degree boundary is asserted.

## Exact source scope

The sole external theorem input is Christopher D. Long's one-radial factorial theorem and its fixed-power prime-index nonvanishing, in the pinned manuscript *A Factorially Weighted Constant-Term Theorem on Algebraic Tori and the Gaussian Moments Conjecture in Dimension Two*, commit 2f1ba8ae5d3d282141aca2285326a2d31ce20cd2, thm:main and thm:prime-nonvanishing. The source is available at
https://raw.githubusercontent.com/octonion/mathematics/2f1ba8ae5d3d282141aca2285326a2d31ce20cd2/gmc/gmc2_stronger_arbitrary_torus.tex .
The paired algebra embedding, finite field-of-coefficients transfer, n<=2 positive theorem and native quadratic counterpair are given in [The sharp Gaussian--Laurent boundary and circular homogenization](SHARP_GAUSSIAN_LAURENT_BOUNDARY_AND_CIRCULAR_HOMOGENIZATION.md). The exact existing Long input is recorded in [Laurent and binary Gaussian vanishing equivalences](LAURENT_AND_BINARY_GAUSSIAN_VANISHING_EQUIVALENCES.md). Equations (4), (7) and the full Hermite marker transfer are proved above. These are explicit algebraic corollaries with the stated phase, marker, degree and field hypotheses; no finite moment computation or priority assertion is used.
