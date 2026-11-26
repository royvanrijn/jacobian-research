# Independent affine powers: a Schur consequence for Strong Factorial windows

**Theorem 1.** Let \(K\) be any characteristic-zero field and
\(\Lambda_m(X^\alpha)=\alpha!\). Suppose
\[
 F=C\prod_{i=1}^m(a_iX_i+b_i)^{d_i}\ne0,\qquad d_i\ge1,\ C\ne0. \tag{1}
\]
Let \(r\) be the number of indices with \(a_ib_i\ne0\).
Then at most \(r\) positive integers \(n\) satisfy \(\Lambda_m(F^n)=0\).
In particular every window of \(r+1\) consecutive positive moments contains a nonzero one, and \(F\) satisfies the Strong Factorial Conjecture for its actual support. This includes structured families in three variables,
rather than arbitrary SFC(3).

## 1. Exact primary input

For every \(N\ge1\),
\[
 E_N(Z)=\sum_{j=0}^N Z^j/j!
\]
is irreducible in \(\mathbb Q[Z]\). The checked primary research source is
Michael Filaseta, *A generalization of an irreducibility theorem of I.
Schur*, [author PDF](https://people.math.sc.edu/filaseta/papers/schurpaper.pdf),
Theorem 1 on PDF page 1. The exact [retained author PDF](sources/filaseta-schur-generalization-author.pdf) supplies the imported statement. Its stronger formulation permits integer
coefficients \(a_j\) with \(|a_0|=|a_N|=1\).
A complete proof of this precise theorem is supplied within the same
paper: the Dumas Newton-polygon theorem, Lemma 1, the
Sylvester--Schur consecutive-integer prime-divisor result (Theorem 6),
and the deduction following Lemma 3, on PDF pages 4--7.
These statements and the complete supplied deduction were read.
The original attribution is Schur, *Einige Satze uber Primzahlen mit
Anwendungen auf Irreduzibilitatsfragen I*, Sitzungsberichte der
Preussischen Akademie der Wissenschaften, Physikalisch-Mathematische
Klasse (1929), 125--136, reference 19. No claim is made here to have
read that original 1929 scan.

The only arithmetic input needed below is irreducibility of \(E_N\).
Since distinct \(E_N,E_M\) have different degrees, they are coprime in
\(\mathbb Q[Z]\). There is a rational polynomial Bezout identity
\(U E_N+V E_M=1\) for each \(N\ne M\). The identity remains valid in
\(K[Z]\), so a fixed element of any extension field of \(K\) is a root
of at most one \(E_N\).

## 2. Factorial moments and the exact support

If \(a_ib_i\ne0\), direct binomial expansion gives
\[
 \Lambda_1((a_iX_i+b_i)^N)
   =a_i^N N! E_N(b_i/a_i). \tag{2}
\]
For a fixed \(i\), the positive integers \(d_in\) are distinct.
The preceding coprimeness therefore permits a zero in (2), with
\(N=d_in\), at at most one positive integer \(n\).
If \(a_i=0\), the nonzero factor is constant and every moment is nonzero.
If \(b_i=0\), it is a monomial and every moment is nonzero.
The zero polynomial factor is excluded by \(F\ne0\).

Independence of distinct variables gives
\[
 \Lambda_m(F^n)=C^n
    \prod_i\Lambda_1((a_iX_i+b_i)^{d_in}). \tag{3}
\]
Its zero index set is contained in the union of \(r\) sets of size at
most one. This proves the stated bound.

The actual number of monomials is
\[
 N(F)=\prod_{i:a_ib_i\ne0}(d_i+1). \tag{4}
\]
Every binomial coefficient in an active factor is nonzero in
characteristic zero, and distinct-variable products cannot collide.
If \(r=0\), (3) never vanishes. If \(r\ge1\), then
\(N(F)\ge2^r\ge r+1\). Thus any window of \(N(F)\) consecutive positive
orders contains an order outside the at-most-\(r\) zero set.
This is exactly the Strong Factorial window for (1).

The same argument allows extra monomial factors in additional disjoint
variables. It does not permit an extra monomial multiplier in an active
affine variable: that would change (2) to a shifted factorial polynomial.
It also does not replace a single-variable factor by an arbitrary
linear form in a block of variables.

## 3. Scope comparison

Rocks' primary author preprint, Section 3.3, constructs new SFC examples
by adding a new independent linear variable to a known example.
The product theorem above uses multiplication and the stronger
at-most-one-zero property of each factor. A bounded repository search
did not identify this exact product statement; that is not a novelty
or completeness assertion.

For \(m=2,d_1=d_2=1\), (1) is a separated affine-product quadratic.
A separate elementary proof requires no Schur input: each affine
moment sequence has no zeros at adjacent orders or at orders two apart,
and hence their product cannot have three consecutive zero moments.
The Schur argument above strengthens the conclusion to at most two
zero orders in total, and supplies the indicated arbitrary-dimensional
extension. It does not settle a general polynomial of total degree two
in two or three variables.


## 4. Sharpness and the dimension-two target

The bound \(r+1\) is sharp over \(\mathbb C\). For each \(1\le i\le r\),
choose any root \(\zeta_i\) of \(E_i\). Such a root exists and is nonzero,
because \(E_i\) has positive degree and constant coefficient one.
Pairwise coprimality shows that
\(\Lambda_1((X+\zeta_i)^n)=n!E_n(\zeta_i)\) vanishes exactly at \(n=i\).
Consequently
\[
 F=\prod_{i=1}^r(X_i+\zeta_i)
\]
has zero moment indices exactly \(1,\ldots,r\), and its moment at \(r+1\)
is nonzero. Its support has exactly \(2^r\) terms. Thus neither the bound
on the number of zero indices nor the uniform window \(r+1\) can be reduced
for this class. This sharpness assertion concerns complex coefficients.

For binary quadratics, the active case with two affine factors on independent
variables is \((ax+b)(cy+d)\). The theorem gives a nonzero moment in every
length-three window, although its support has four terms. A product of oblique
affine forms in the same two variables is a larger class and is not covered
by the factorization of moments in (3).

The [homogeneous quadratic theorem](BINARY_HOMOGENEOUS_QUADRATIC_STRONG_FACTORIAL.md)
is a separate result. Together these statements do not close arbitrary
inhomogeneous binary quadratics or full \(SFC(2)\). Their structured
three-variable specializations do not close full \(SFC(3)\).
Schur's irreducibility is an indispensable classical input here.
The supplied proof and its application were checked internally; no external
review, formal verification, or priority claim is made.
