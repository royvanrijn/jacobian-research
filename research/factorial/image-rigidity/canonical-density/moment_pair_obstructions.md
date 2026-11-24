# Counterexamples to positivity, degree drops, and genuine moment collisions

Internal research note, 22 November 2026. The two previously constructed
bad modifier pairs have coprime moment polynomials. A separate explicit
one-parameter family can be tuned to a genuine common-root pair of moments,
but its marginal distributions are not the linked initial factorial laws.
No SFC counterexample, or unrestricted exclusion theorem, is obtained.

## 1. What would constitute a factorial counterexample?

Let `L(X^a Y^b)=a!b!`, with nonnegative integer exponents. Write

\[
 J_m(z;r,u)=\sum_{j=0}^m(-1)^j\binom mj
   \frac{(mr+1)_j}{(m(r+u)+2)_j}z^{m-j}.
\]

For nonnegative integers `r,u`, set

\[
 f=X^rY^u((z-1)X+zY).
\]

The exact identity is

\[
 \mathcal L(f^m)=
 \frac{(m(r+u+1)+1)!(mr)!(mu)!}{(m(r+u)+1)!}
 J_m(z;r,u).                                             \tag{1}
\]

Indeed, in the exponential integral for the factorial functional put
`X=rho*t`, `Y=rho*(1-t)`. The radial integral gives the first factorial;
normalizing the remaining beta integral gives (1). Alternatively, expansion
of both sides proves it as a rational polynomial identity. It therefore
holds over every characteristic-zero field, without an ordering of that
field.

Consequently a common root of the **actual** `J_n,J_(n+1)` would give a
two-support counterexample in this binary monomial-linear family over its
algebraic field of definition. Neither endpoint `z=0,1` is a root, by the
positive rising-factorial endpoint formulas. The two coefficients of `f`
would thus be nonzero. The construction requires one common fixed integer
profile `r,u`, and consecutive powers of that same `f`.

Arbitrary positive beta mixtures do not satisfy this linked-power
requirement. Even an exact common root for two independently chosen
moment laws is not, by itself, a factorial counterexample.

## 2. The exact resultant identity, and its limitation

Suppose the monic Euclidean chain is regular at every step:

\[
 P_{k+1}=(z-\alpha_k)P_k+\beta_kP_{k-1},
 \qquad \deg P_k=k,\quad P_0=1,\quad \beta_k\ne0.
\]

Then

\[
 \boxed{\operatorname{Res}(P_{n+1},P_n)
                  =\prod_{k=1}^n\beta_k^{\,k}.}           \tag{2}
\]

To prove it, evaluate the first polynomial at the roots of `P_k`. Since
`k(k+1)` is even, reversing the order of the two resultants introduces
no sign. Thus

\[
 \operatorname{Res}(P_{k+1},P_k)
    =\beta_k^k\operatorname{Res}(P_k,P_{k-1}),
\]

and iteration proves (2). A negative multiplier is compatible with a
nonzero resultant. For a nondegenerate real probability law, its even-degree moment
polynomial is strictly positive on the real axis. Evaluating the
resultant at real roots and conjugate root pairs shows that a coprime
consecutive pair has strictly positive resultant. In a regular chain,
negative odd-index multipliers must therefore occur an even number of
times; negative even-index multipliers do not change the product sign.

The regularity hypothesis is essential. The nominal multiplier
`beta_k` is the coefficient of degree `k-1` in the remainder. If it
vanishes, the remainder may merely have smaller degree. For monic
polynomials `A,B` of degrees `m>n`, a nonzero remainder `R=rho*Rhat`
of degree `d<n` instead gives

\[
 \operatorname{Res}(A,B)=(-1)^{mn}\rho^n
                          \operatorname{Res}(B,\widehat R).
                                                               \tag{3}
\]

A nonzero constant remainder proves coprimality. Only a complete Euclidean
gcd computation, or an equivalent resultant certificate, resolves the
common-root question. One cannot set a nominal `beta_k` to zero in (2)
while retaining nonexistent normalized later polynomials.

## 3. Both retained bad pairs are exactly coprime

For a probability law `lambda`, write

\[
 M_k[\lambda](z)=\int_0^1(z-t)^k\,d\lambda(t).
\]

The following are the full inputs for the degree-eleven/degree-ten pairs.

### 3.1 Strict modifier/transport example

With `epsilon=1/10000`, put

\[
 \nu=(1-\epsilon)\left[
  \frac{701}{1344}\operatorname{Beta}(111,143)
 +\frac{643}{1344}\operatorname{Beta}(132,122)\right]
 +\epsilon\operatorname{Beta}(111,122),
\]

\[
 \mu=(1-\epsilon)\operatorname{Beta}(122,134)
                  +\epsilon\operatorname{Beta}(112,123).
\]

### 3.2 Log-concavity example

Let `E_5(v)=sum_(j=0)^5 v^j/j!`. The law `nu` is proportional to

\[
 \operatorname{Beta}(56,67)(dt)\,
 E_5\!\left(\frac{55}{2}
                \left(2t-\frac{112}{123}\right)^2\right),
\]

and `mu=Beta(62,74)`.

For **each** pair the certificate proves

\[
 \boxed{\gcd_{\mathbb Q[z]}(M_{11}[\nu],M_{10}[\mu])=1.}  \tag{4}
\]

It records the complete rational coefficients and explicit polynomial
Bezout vectors modulo the prime `65537`. All rational denominators are
units modulo that prime, and both polynomials remain monic of their
original degrees. The recorded vectors satisfy `A*upper+B*lower=1`
exactly modulo the prime. The Sylvester determinant consequently has
nonzero reduction and is nonzero over `Q`. Its nonvanishing persists
over every characteristic-zero field extension.

The independent verifier reconstructs every input moment, using factorial
quotients for beta moments. For the tilt it expands the finite sum directly
by binomial coefficients, independently of the generator's repeated
polynomial multiplication. It also reconstructs the first two multipliers,
which agree with the previous positive and negative values.

Thus the negative following multipliers in these two examples are not
common-zero obstructions, even before imposing canonical reachability.
Changing their coefficients would constitute a different example.
This note does not claim that arbitrary finite-moment embeddings are
excluded merely because their chosen density representations are
noncanonical.

## 4. Tuning a positive mixture: degree drop and collision are different

Keep `mu=Beta(2,2)`, take represented degree `k=3`, and define

\[
 \nu_\theta=\theta\operatorname{Beta}(4,4)
                   +(1-\theta)\operatorname{Beta}(5,5),
 \qquad 0\le\theta\le1.
\]

Put `y=z-1/2` and

\[
 A=\frac3{20},\quad
 b_\theta=\frac{20\theta-9}{660},\quad
 c_\theta=\frac3{2288}+\frac\theta{1716}.
\]

Exact central beta moments give

\[
 P_3=y(y^2+A),\qquad
 Q_4=y^4+(A+b_\theta)y^2+c_\theta,
\]

so the remainder after subtracting `yP_3` is
`b_theta*y^2+c_theta`. Evaluating `Q_4` at `y=0` and the two roots of
`y^2+A` proves, uniformly in `theta`,

\[
 \boxed{\operatorname{Res}(Q_4,P_3)
           =c_\theta(c_\theta-A b_\theta)^2.}              \tag{5}
\]

At `theta=9/20` the nominal quadratic pivot is zero, yet

\[
 Q_4-yP_3=\frac9{5720}\ne0.
\]

Hence the two polynomials are coprime and their resultant is
`(9/5720)^3>0`. The pivot changes sign across this parameter; tuning it
through zero does not produce a common root there.

At the different parameter `theta=72/85`, however,

\[
 P_3=y\left(y^2+\frac3{20}\right),\qquad
 Q_4=\left(y^2+\frac3{20}\right)
                 \left(y^2+\frac9{748}\right).
\]

The gcd is exactly `y^2+3/20`, and the common roots are
`z=1/2 +/- i*sqrt(3/20)`. Both distributions remain strictly positive
inside `(0,1)`. This is a genuine common-root counterexample for
**independently selected positive moment laws**.

It is not an SFC counterexample. The construction does not provide a
single integer profile and a single polynomial whose consecutive powers
produce these two laws. A proof of nonreachability as canonical densities
is a separate statement; it must not be silently strengthened to exclude
all possible finite-moment representations.

## 5. Reproduction and scope

`moment_pair_certificates.json` contains the complete finite
inputs. Run `python verify_moment_pairs.py`; the verifier uses
only Python's standard library and no sibling imports. It checks the
moment reconstructions, all modular Bezout coefficients, the first two
pivots, the centered family coefficients and both special parameters.
Four exact rational Sylvester determinants at distinct parameters verify
(5) independently: both sides have degree at most three in `theta`.
Two further determinants check the degree-drop and common-root points.

A separate generator used SymPy for rational gcd and symbolic resultant
cross-checks. The retained verifier and certificate do not depend on that
generator or on SymPy.

The positive result is exact separation of three questions: failure of a
positivity method, failure of a regular Euclidean degree pattern, and an
actual simultaneous moment zero. The two old obstructions supply neither
a coprimality failure nor an SFC counterexample; the tuned mixture does fail coprimality in the
larger independent-law class. Uniform exclusion on the actual factorial
orbit remains open.
