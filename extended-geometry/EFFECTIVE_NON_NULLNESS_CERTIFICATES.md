# Effective certificates of non-nullness

## 1. Executive conclusion

The lower-face proof in
[`TWO_REAL_GMC_LOWER_FACE_THEOREM.md`](TWO_REAL_GMC_LOWER_FACE_THEOREM.md)
is algorithmic for exact rational, integer, or explicitly presented algebraic
coefficients.

For rational input, a certificate can consist of

\[
(\rho,\theta,S_0,r,n,c,p),\qquad
n=\rho r,\qquad c=\operatorname{CT}_T(P_0^r)\ne0,
\]

where \(S_0\) is the exposed balanced lower face and \(p\) is a prime at
which the input and \(c\) have good reduction.  The certificate returns

\[
m=rp
\]

and proves directly that

\[
\mathbb E(P^m)\ne0.
\]

The steps involving the lower face and the prime are effective with modest
complexity.  The main complexity problem is finding, or bounding, the first
\(r\) for which \(\operatorname{CT}(P_0^r)\ne0\).

There is a direct route to a **singly-exponential degree bound** for \(r\):
the generating series

\[
H_f(t)=\sum_{j\geq0}\operatorname{CT}(f^j)t^j
\]

of a one-variable Laurent polynomial \(f=P_0\) is an algebraic power series,
being a sum of residues of a bivariate rational function.  Specializing the
explicit residue and sum-of-roots estimates of Bostan--Dumont--Salvy gives
the following bound.  If

\[
f(z)=\sum_{j=-a}^{b}u_jz^j,\qquad u_{-a}u_b\ne0,\qquad a,b>0,
\]

then some \(r\) with \(\operatorname{CT}(f^r)\ne0\) satisfies

\[
r\leq
2(a+b)\binom{a+b-1}{a-1}
\leq 2(a+b)2^{a+b-1}.                            \tag{1.1}
\]

Section 5 gives the derivation.  It should be promoted to a named lemma, with
the specialization independently checked in a CAS and then formalized, before
being used as a black-box repository theorem.

Thus:

* algorithmic termination is unconditional, by Duistermaat--van der Kallen;
* existing effective residue theory gives the explicit singly-exponential
  bound (1.1);
* a polynomial bound for \(r\) remains the sharper and genuinely interesting
  open question;
* once \(r\) is known, \(p\), and hence \(m=rp\), has an elementary polynomial
  bound in \(r\) and the coefficient height.

## 2. Exact input model

An effective statement needs an exact coefficient model.  The cleanest first
version is

\[
P=\sum_{k\in S}T^kB_k(U)\in\mathbb Q[T^{\pm1},U].
\]

Record:

* \(s\): the number of nonzero monomials of \(P\);
* \(D_T=\max |k|\): the angular degree;
* \(D_U=\max\deg B_k\): the radial degree;
* \(H\): a bit-height bound for rational coefficients.

The same method works over a number field given by a defining polynomial and
coordinates in a fixed basis.  It is not an algorithm for arbitrary complex
numbers supplied only as black boxes.  For a general finite-type
\(\mathbb Z\)-domain, existence of good reductions is enough for the abstract
proof, but effectivity also requires an explicit presentation and effective
prime-ideal arithmetic.

## 3. The certificate algorithm

Write

\[
\nu_k=\operatorname{ord}_U B_k,\qquad
b_k=[U^{\nu_k}]B_k.
\]

Assume that the angular support is balanced:

\[
0\in\operatorname{conv}(S).
\]

### Step A: expose the lower balanced face

Solve

\[
\rho=\min\left\{
\sum_k\lambda_k\nu_k:
\lambda_k\geq0,\ 
\sum_k\lambda_k=1,\ 
\sum_k k\lambda_k=0
\right\}.
\]

Produce a rational supporting line

\[
\nu_k\geq\rho+\theta k
\]

and set

\[
S_0=\{k:\nu_k=\rho+\theta k\},\qquad
P_0(T)=\sum_{k\in S_0}b_kT^k.
\]

Because this is the lower hull of finitely many lattice points
\((k,\nu_k)\subset\mathbb Z^2\), it can be computed by rational convex-hull
arithmetic.  In fact, only the lower hull above angular coordinate zero is
needed.  This stage is polynomial time in the ordinary bit model.

The face certificate is checked using only rational inequalities, equality
on \(S_0\), and the presence of weights on both weak sides of zero.

### Step B: find a nonzero angular moment

Search over \(r=1,2,\ldots\), computing

\[
c_r=\operatorname{CT}_T(P_0^r).
\]

Stop at the first \(r\) with \(c_r\ne0\), and put \(c=c_r\).
Duistermaat--van der Kallen proves termination because
\(0\in\operatorname{conv}(S_0)\).

For a fixed proposed \(r\), \(c_r\) can be computed by dynamic programming on
angular exponent, sparse multiplication, or coefficient extraction from a
binary powering circuit.  This is polynomial in \(r\), the angular width, and
the coefficient bit-size, although it is only pseudo-polynomial if angular
exponents are allowed to be enormous binary integers unrelated to the degree
parameter.

Every weight-zero monomial of \(P^r\) supported on the face has the same
radial degree

\[
n=\rho r\in\mathbb Z_{\geq0},
\]

and

\[
F_r(U):=\operatorname{CT}_T(P^r)
=cU^n+\text{terms of degree \(>n\)}.
\]

### Step C: choose a good prime

Clear denominators once, or work in the localization of \(\mathbb Z\) at the
input denominators.  Choose any rational prime \(p\) which divides neither
the relevant denominators nor the numerator of \(c\).

No factorization is required: enumerate primes, test the finitely many
integers modulo \(p\), and stop at the first good one.

### Step D: return \(m=rp\)

Let

\[
F_{rp}(U)=\operatorname{CT}_T(P^{rp}).
\]

Every term of \(F_{rp}\) has radial degree at least \(np\).  In characteristic
\(p\),

\[
F_{rp}(U)\equiv F_r(U)^p.
\]

Normalize the Gaussian moment by \((np)!\):

\[
\frac{\mathbb E(P^{rp})}{(np)!}
=\sum_{d\geq np}[U^d]F_{rp}(U)\frac{d!}{(np)!}.
\]

Modulo \(p\):

* the coefficient at degree \(np\) is \(c^p\);
* a degree not divisible by \(p\) has zero coefficient by Frobenius;
* for \(d=jp\) with \(j>n\), the quotient \(d!/(np)!\) is divisible by \(p\).

Therefore

\[
\frac{\mathbb E(P^{rp})}{(np)!}\equiv c^p\not\equiv0\pmod p,
\]

so \(\mathbb E(P^{rp})\ne0\) in characteristic zero.

This is a direct nonvanishing proof, not merely a contradiction under the
assumption that all moments vanish.

## 4. Complexity of the prime

Suppose first that the face coefficients are integers and
\(|b_k|\leq 2^H\), with \(q=|S_0|\).  Since the sum of the multinomial
coefficients is \(q^r\),

\[
|c|\leq q^r 2^{rH},
\qquad
\log |c|=O(r(H+\log q)).
\]

If every prime up to \(x\) divided a nonzero integer \(C\), their product
would divide \(C\).  Standard elementary bounds for the primorial therefore
give a prime not dividing \(C\) of numerical size

\[
p=O(\log |C|).
\]

After clearing rational denominators, the same statement applies to the
product of the denominator data and the numerator of \(c\).  With a naive
common denominator this gives

\[
p=\operatorname{poly}(r,s,H).
\]

Consequently, if

\[
r\leq R(s,D_T,H),
\]

then one obtains, conservatively,

\[
m=rp\leq
\operatorname{poly}(R,s,H).
\]

In particular, a singly-exponential bound for \(r\) gives a
singly-exponential bound for \(m\).  The prime is not the asymptotic
bottleneck.

Over a number field, replace \(c\) by a cleared integral representative and
avoid the rational primes dividing its norm, the input denominator ideal,
and any required discriminant data.  The same primorial argument bounds a
good rational prime in terms of the logarithmic norm/height, although the
residue field may be \(\mathbb F_{p^f}\).

## 5. An explicit singly-exponential bound for \(r\)

Let

\[
f(z)=P_0(z)=\sum_{j=-a}^{b}u_jz^j,
\qquad u_{-a}u_b\ne0,
\]

where \(a,b>0\), and put \(d=a+b\).  Then

\[
H_f(t)=\sum_{r\geq0}\operatorname{CT}(f^r)t^r
=\operatorname{CT}_z\frac{1}{1-tf(z)}.
\]

After multiplication by \(z^a\), the denominator is

\[
z^a-tg(z),\qquad \deg g=d.
\]

More explicitly,

\[
H_f(t)
=\frac{1}{2\pi i}\int
\frac{z^{a-1}}{z^a-tg(z)}\,dz.                   \tag{5.1}
\]

For small \(t\), the denominator

\[
Q(t,z)=z^a-tg(z)
\]

has exactly \(a\) small branches tending to zero.  Equation (5.1) is the sum
of the residues of \(z^{a-1}/Q(t,z)\) at those branches.

The polynomial \(Q\) is square-free over \(K(t)\).  Indeed, a common
positive-\(z\)-degree factor of \(Q\) and
\(\partial Q/\partial z\) would divide
\(aQ-z\,\partial Q/\partial z=t(zg'(z)-ag(z))\).
It would therefore be represented by a factor in \(K[z]\); comparison of the
constant and \(t\)-coefficients of \(z^a-tg(z)\) would make it divide both
\(z^a\) and \(g(z)\), impossible because \(g(0)\ne0\).

Apply Theorem 8 of Bostan--Dumont--Salvy to

\[
\frac{z^{a-1}}{z^a-tg(z)}.
\]

The denominator has bidegree \((1,d)\) in \((t,z)\).  Their residue
polynomial \(R(t,y)\), whose roots include all \(d\) residues, satisfies

\[
\deg_y R\leq d,\qquad \deg_t R\leq 2d.            \tag{5.2}
\]

Next apply their Theorem 12 to the sum of \(a\) roots of \(R\).  It gives an
annihilating polynomial \(Q_H(t,y)\) for \(H_f(t)\) with

\[
\deg_y Q_H\leq\binom da,\qquad
\deg_t Q_H\leq
2d\binom{d-1}{a-1}.                              \tag{5.3}
\]

This construction gives an annihilator, which need not be irreducible after
a special coefficient substitution.  Factor it over \(K(t)\) and select the
irreducible factor of the branch \(H_f\); its \(t\)-degree cannot increase.

Let \(Q(t,y)\in K[t,y]\) be the irreducible polynomial of the branch
\(H_f(t)\).  If \(H_f\ne1\), then \(Q(t,1)\ne0\).  If

\[
H_f(t)-1=\gamma t^r+O(t^{r+1}),\qquad\gamma\ne0,
\]

then

\[
Q(t,1)-Q(t,H_f(t))
\]

is divisible in \(K[[t]]\) by \(1-H_f(t)\).  Since
\(Q(t,H_f(t))=0\),

\[
r\leq\operatorname{ord}_{t=0}Q(t,1)
\leq\deg_t Q.
\]

Combining this observation with (5.3) proves

\[
r\leq2d\binom{d-1}{a-1}
\leq2d\,2^{d-1}.                                 \tag{5.4}
\]

If \(f\) has a nonzero constant term, \(r=1\).  Otherwise a balanced support
has \(a,b>0\), so (5.4) covers every case.

The coefficient height does not enter the exponent bound.  It enters the
cost of exact arithmetic and the size of the good prime.  The bound is
polynomial when one of \(a,b\) is fixed and singly exponential in the worst
case \(a\asymp b\).

This route is preferable to invoking a general effective Nullstellensatz:
it exploits the fact that the lower face is a Laurent polynomial in one
angular variable and targets the special moment sequence itself.

## 6. Why torus nullcone algorithms do not by themselves bound \(r\)

For fixed weights \(k_1,\ldots,k_q\), the coefficient vector
\((b_1,\ldots,b_q)\) carries a one-dimensional torus action.  Nullcone
membership is combinatorial:

\[
(b_i)\text{ is outside the nullcone}
\iff
0\in\operatorname{conv}\{k_i:b_i\ne0\}.
\]

This is decidable by linear programming, and much more general orbit
problems for torus actions admit polynomial-time algorithms.

But the invariant ring contains all zero-weight monomials, whereas
\(\operatorname{CT}(f^r)\) is one distinguished linear combination of the
invariant monomials of degree \(r\).  A low-degree separating invariant need
not be a nonzero constant-term moment.  Duistermaat--van der Kallen says that
the whole moment sequence detects exactly the same nullcone, but its original
proof is analytic/geometric and does not state a cutoff.

Therefore a degree bound for torus invariants cannot simply be substituted
for a bound on \(r\).

## 7. Moment ideals: what follows and what does not

Let

\[
M_r=\operatorname{CT}\left(\sum_i x_iT^{k_i}\right)^r
\in K[x_1,\ldots,x_q].
\]

For a fixed support, Duistermaat--van der Kallen implies

\[
V(M_1,M_2,\ldots)=\mathcal N_T,
\]

the torus nullcone.  Noetherianity gives some finite cutoff set-theoretically,
and a uniform first-nonzero bound \(R\) gives the concrete equality

\[
V(M_1,\ldots,M_R)=\mathcal N_T.
\]

Equivalently, the nullcone ideal is contained in

\[
\sqrt{(M_1,\ldots,M_R)}.
\]

This is the correct immediate “finite generation” consequence.

It does **not** follow that the full moment ideal

\[
(M_1,M_2,\ldots)
\]

is generated by \(M_1,\ldots,M_R\).  Radical generation and literal ideal
generation are different; the repository's cubic nullcone theorem already
exhibits this distinction.  Bounds for literal moment-ideal generation need
additional recurrence or ideal-membership arguments.

If angular degree and radial degree are bounded, there are only finitely many
possible lower-face supports and radial-order patterns.  A uniform \(r\)
bound over those patterns converts the pointwise certificate into a finite
radical-nullcone test for the whole bounded-degree family.

## 8. Automated nullcone testing

There are two different algorithms.

1. **Combinatorial torus test.** Inspect the nonzero angular weights.  This
   decides whether the coefficient point is in the torus nullcone in
   polynomial time.
2. **Moment witness.** If it is outside the nullcone, compute the lower face,
   find \(r\), choose \(p\), and return \(m=rp\) together with a proof that
   the actual Gaussian moment is nonzero.

The first is the faster decision procedure.  The second produces a witness
in the language of Gaussian moments and is valuable when a downstream
system trusts moment computations but not an external convex-geometric
classification theorem.

## 9. Lean-checkable certificate format

The existing package [`formalization/gmc2`](../formalization/gmc2) already
formalizes the modular core.  Its current global proof imports three axioms:
Duistermaat--van der Kallen, finite-type specialization, and supporting-face
extraction.

A per-input certificate can avoid all three:

```text
input:
  exact rational coefficients of P

face:
  rho, theta, face indices
  proofs/checks of every supporting inequality
  a zero-sum exponent vector alpha of total size r

constant term:
  r
  c = CT(P0^r), represented exactly or by an arithmetic circuit
  proof/check that c != 0

prime:
  p
  primality certificate
  denominator residues are nonzero mod p
  c mod p != 0

output:
  n = rho*r
  m = r*p
```

The checker verifies:

* rational supporting inequalities;
* exact or modular constant-term extraction;
* Frobenius/constant-term commutation;
* the factorial-quotient divisibility lemma;
* nonzero normalized moment modulo \(p\).

For compactness, \(P_0^r\) and the relevant coefficient can be represented by
a repeated-squaring arithmetic circuit rather than an expanded Laurent
polynomial.  A Pratt certificate or a deterministic verified primality
procedure can certify \(p\).  The mathematical kernel is already represented
by `FactorialQuotient.lean`, `ConstantTerm.lean`, and
`PrimeIsolation.lean`.

## 10. Recommended research program

### A. Audit and package the explicit residue bound

Independently reproduce the specialization of Theorems 8 and 12 to

\[
\operatorname{CT}_z(1-tz^{-a}g(z))^{-1}.
\]

Then promote (5.4) to a checked lemma:

\[
0\in\operatorname{conv}(\operatorname{supp}f),\quad
f\ne0
\Longrightarrow
\operatorname{CT}(f^r)\ne0
\quad\text{for some }
1\leq r\leq2(a+b)\binom{a+b-1}{a-1}.
\]

The audit should explicitly check the orientation/sign in the residue sum,
the count of \(a\) small branches, square-freeness of \(z^a-tg(z)\), and
nongeneric coefficient specializations.

### B. Search for lower-bound families

For supports in \([-D,D]\), solve

\[
\operatorname{CT}(f^1)=\cdots=\operatorname{CT}(f^N)=0
\]

on coefficient charts while requiring both signs to remain present.  Record
the largest achievable \(N\).  This distinguishes plausible polynomial
bounds from genuinely exponential behavior.

### C. Implement certificate emission

Add a script which:

1. parses \(P\in\mathbb Q[T^{\pm1},U]\);
2. computes the lower hull and balanced face;
3. streams constant terms until it finds \(r\);
4. finds the first good prime;
5. emits JSON containing the certificate and \(m=rp\);
6. independently verifies the normalized moment modulo \(p\).

### D. Separate three bounds

Track independently:

* the decision complexity of non-nullness;
* the least witness exponent \(r\);
* the size and verification time of a proof certificate.

A large \(r\) does not force a large certificate if the coefficient
extraction is supplied as a succinct arithmetic circuit, but verification
time still depends on the chosen proof representation.

## 11. References

* J. J. Duistermaat and W. van der Kallen,
  [*Constant terms in powers of a Laurent polynomial*](https://webspace.science.uu.nl/~kalle101/powers.pdf),
  Indagationes Mathematicae 9 (1998), 221--231.  The one-variable proof uses
  the nontrivial singularity of the constant-term generating function; the
  paper does not state an effective first-index bound.
* A. Bostan, L. Dumont, and B. Salvy,
  [*Algebraic Diagonals and Walks: Algorithms, Bounds, Complexity*](https://arxiv.org/abs/1510.04526),
  Journal of Symbolic Computation 83 (2017), 68--92.  It gives algorithms and
  explicit size bounds for annihilating polynomials of bivariate rational
  diagonals, exponential in input degree.
* P. Bürgisser, M. L. Doğan, V. Makam, M. Walter, and A. Wigderson,
  [*Polynomial Time Algorithms in Invariant Theory for Torus Actions*](https://drops.dagstuhl.de/storage/00lipics/lipics-vol200-ccc2021/LIPIcs.CCC.2021.32/LIPIcs.CCC.2021.32.pdf),
  CCC 2021.  This supplies the relevant complexity context for torus
  nullcones and orbit problems, but not a cutoff for the special sequence of
  constant-term moments.
* J. Kollár,
  [*Sharp Effective Nullstellensatz*](https://www.math.ucdavis.edu/~deloera/MISC/LA-BIBLIO/trunk/Kollar/kollarnullstellen.pdf),
  Journal of the AMS 1 (1988), 963--975.  This is a fallback route for
  radical/ideal certificates once a finite polynomial system has been
  specified, but it does not by itself make the infinite moment cutoff
  effective.
