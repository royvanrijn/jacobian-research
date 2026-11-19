# Safety in the three-height quintic Image family

Work over a field \(K\) of characteristic zero. On \(t^2+2xy=1\) use
\[
 \Lambda(F)=\operatorname{CT}_x\frac12\int_{-1}^1
             F\bigl(x,(1-t^2)/(2x),t\bigr)\,dt .
\]
The integral denotes the rational polynomial rule. A phase is **pure** if
all positive moments vanish and **safe** if every fixed polynomial marker
has eventual mixed vanishing. Least degree is the sphere filtration degree.
Use the odd phase
\[
 P=\frac{t}{x}(C(x)+t^2D(x)+t^4E(x)),
 \qquad
 \deg C\le5,\quad \deg D\le3,\quad \deg E\le1,\quad
 C(0)+D(0)+E(0)=0.
\]

**Theorem ThreeHeightQuinticSafety.** The following are equivalent.

1. \(P\) is pure.
2. The even moments of orders \(2,4,\ldots,18,546,630,650\) vanish.
3. \(P\) has only strictly positive weights or only strictly negative
   weights for \(\mathrm{wt}(x,y,t)=(1,-1,0)\). The zero phase is included.

In coefficient terms the two possibilities are
\(x^2\mid C,D,E\), or all three coefficients are constant. Every pure
phase in the box is safe for sphere markers and for original polynomial
markers of the balanced two-pair lift. The theorem covers both vanishing
and nonvanishing leading height coefficients.

## 1. The moment equations and reduction to a prime field

Put \(z=t^2\). There are unique coefficients \(a,b,c,d,e,f,g,h,i,j,k\)
such that \(P=t\,U(x,z)/x\), where
\[
 U=(1-z)(a+bz)+x(c+dz+ez^2)+x^2(f+gz)+x^3(h+iz)+x^4j+x^5k .
 \tag{1}
\]
All odd moments vanish by parity. At an even order \(m\), the equation is
\[
 \mu_m=\sum_{\ell=0}^{2m}
       \frac{[x^mz^\ell]U(x,z)^m}{m+2\ell+1}=0 .
 \tag{2}
\]
The two coefficient actions are overall scalar multiplication of \(P\)
and
\[
 (U_0,U_1,\ldots,U_5)\longmapsto
       (U_0,\rho U_1,\rho^2U_2,\ldots,\rho^5U_5),
 \qquad U_r=[x^r]U .
\]
The second action is \(P(x,t)\mapsto\rho P(\rho x,t)\), and multiplies
\(\mu_m\) by \(\rho^m\). One may use both actions over extensions.

For the characteristic-zero exclusion it suffices to exclude the
normalized cases over \(\overline{\mathbb F}_{2731}\), using the reductions
of appropriate integer-scaled equations from (2). Indeed, if a bad case
has a point over any characteristic-zero field, the polynomial equations
and its specified nonzero coefficients give a point over
\(\overline{\mathbb Q}\) by the Nullstellensatz. Take a valuation above
2731 on a number field containing that point. If \((a,b)\ne(0,0)\),
first multiply all coefficients to give \(a,b\) minimum valuation zero.
If a coefficient of some \(U_r\), \(r\ge1\), is nonzero, choose the
second scaling with
\[
 v(\rho)=-\min_{r\ge1,\ u\ne0,\ u\text{ a coefficient of }U_r}v(u)/r .
\]
A finite extension supplies such a \(\rho\). All coefficients are
integral, and the reduction retains a nonzero center and at least
one nonzero higher coefficient. If \(a=b=0\), the moments depend only on
\(c,d,e\): in (2), every factor has \(x\)-degree at least one, so only
grade one contributes. Thus if \((c,d,e)\ne0\), it suffices in this
case to normalize and reduce this triple. These arguments retain a
nonzero coordinate in each group whose nonvanishing is required.

## 2. Three high moments remove coefficients before elimination

Let \(q=2731\), a prime. At each of the orders \(m=546,630,650\)
the only denominator divisible by \(q\) in (2) is \(q\) itself. Indeed
the denominators are the odd numbers from \(m+1\) to \(5m+1\),
which contain \(q\) and lie below \(3q\). Consequently, in the reduction
of \(q\mu_m=0\), only the coefficient of \(x^mt^{2730}\) in
\(\bigl(tU(x,t^2)\bigr)^m\) survives. Here and below all these
coefficient operations take place in the residue field.

At \(m=546\), height \(2730=5m\) requires height five in every
factor. The only such terms in \(tU\) are \(-bt^5\) and \(ex t^5\).
The \(x\)-degree requirement then forces \(e^m=0\), so \(e=0\).

For \(b\ne0\), two further coefficients are isolated by the
upper supporting lines of the set of exponent pairs
\((x\text{-degree},t\text{-degree})\).

- With \(e=0\), every pair lies below \(h=5-2r/3\), with equality
  exactly for the possible terms \(-bt^5\) and \(ix^3t^3\).
  At \(m=630\), the sum of heights at total \(x\)-degree \(m\)
  reaches \(2730=(13/3)m\) only on that line. The coefficient is
  \(\binom{630}{210}(-b)^{420}i^{210}\), forcing \(i=0\).
- With \(e=i=0\), use \(h=5-4r/5\). The equality terms are
  \(-bt^5\) and \(kx^5t\). At \(m=650\) the surviving coefficient is
  \(\binom{650}{130}(-b)^{520}k^{130}\), forcing \(k=0\).

Both binomial coefficients are units since their upper indices are
less than \(q\). The case \(b=0\) uses only the first of these
three reductions.

## 3. Exact unit ideals for the low moments

The remaining equations have denominators prime to \(q\).
Cover the cases by the following charts over
\(\overline{\mathbb F}_q\). Normalize the center using overall scaling
and the first nonzero coefficient in the displayed list using the
second scaling. In the zero-center case use only overall scaling of
\(c,d,e\). Earlier entries in each list are zero on the corresponding
chart. Higher coefficients are absent from the moments when the
center is zero.

| Center and prior equations | Leaders, in order | Largest even order used for each leader |
|---|---|---|
| \(a=b=0\) | \(c,d,e\) | \(6,4,2\) |
| \(a=1,\ b=e=0\) | \(c,d,f,g,h,i,j,k\) | \(16,14,12,2,12,6,4,10\) |
| \(b=1,\ e=i=k=0\) | \(c,d,f,g,h,j\) | \(18,12,10,8,12,8\) |

In each chart the equations through the stated order generate the unit
ideal over \(\mathbb F_q\), so they have no point over its algebraic
closure. These are exact finite calculations. The
[checker](verify_three_height_quintic_image_safety.py) constructs all
17 systems from (1), independently compares evaluations with literal
two-variable powers, and checks the reduced Gröbner basis \([1]\) or
a nonzero constant. It uses deterministic finite-field linear algebra
in msolve 0.6.5 (\(-g\,2\), \(-l\,2\)). The
[receipt](three_height_quintic_image_safety_verification.json) records
the chart orders and input hashes.

For an explicit formula for the chart polynomials, in a term of
\([x^m]U^m\) let \(n_0\) be its number of grade-zero factors. If the
height after taking out \(t^m(1-t^2)^{n_0}\) is \(t^{2s}\), multiply
its multinomial coefficient by
\[
 I(m,s,n_0)=\int_0^1t^{m+2s}(1-t^2)^{n_0}\,dt
 =\frac{2^{n_0}n_0!}{\prod_{j=0}^{n_0}(m+2s+2j+1)} .
\]
Integration by parts gives the product. All divisions in these
finite calculations use nonzero denominators. The largest chart, with
\(b=c=1\), has variables \(a,d,f,g,h,j\) and orders \(2,4,\ldots,18\).
Thus the modular calculation includes the free center ratio.

Combining Sections 1–3 shows that if the specified moments vanish, a
zero center has \(c=d=e=0\), and a nonzero center has all higher
coefficients zero. These are precisely the claimed one-sided weight
possibilities. A fixed sphere marker has a finite weight range, so
after multiplying by sufficiently high powers of a one-sided phase
no term of weight zero survives.

For the original markers use the bilinears \(R,Z,W,T\) of the
[contraction identity](TWO_PAIR_IMAGE_MATHIEU_COUNTEREXAMPLE.md#2-contraction-as-circular-gaussian-expectation).
The weights of \(Z,W,R,T\) are \(1,-1,0,0\), from scaling
\(\xi_1,z_1\) inversely. Homogenizing in the radial variable \(R\)
preserves weights. The contraction is equivariant, with output in
\(K[z_1,z_2]\). For a fixed component of a marker of bidegree \((u,v)\),
the output either vanishes by degree or has degree \(v-u\), with
weights in \([-(v-u),0]\). Multiplication by large powers of the
phase lift moves the entire finite weight range outside this interval.
Each marker has finitely many bidegree components, proving safety.
Taking the marker to be one also proves purity for both branches.

## 4. A quintic with nineteen zero initial moments

**Proposition ThreeHeightQuinticNineteenMomentExample.** Over
\(\overline{\mathbb Q}\) the box contains a phase of least degree five
whose first nineteen moments vanish and whose twentieth moment is
nonzero.

Take
\[
 P_\beta=\frac{t(1-t^2)(1+\beta t^2)}x+x^4t .
\]
The \(x\)-weights are \(-1\) and \(4\). Balance requires that the
order be divisible by five, and parity further requires divisibility
by ten. At order \(10r\) the moment is
\[
 \binom{10r}{2r} B_r(\beta),
 \qquad
 B_r(\beta)=\sum_{s=0}^{8r}\binom{8r}{s}
                 I(10r,s,8r)\beta^s .
\]
The polynomials \(B_1\) and \(B_2\) have degrees 8 and 16. The checker
finds that their degrees are retained modulo \(2731\) and that their
modular gcd is one. Therefore they are coprime over \(\mathbb Q\)
by the resultant criterion. Choose any root of \(B_1\). Then
\(\mu_{10}=0\) and \(\mu_{20}\ne0\). The term \(x^4t\) survives
reduction to the unique \(t\)-degree-at-most-one normal form and
has degree five. The other term is \(2ty(1+\beta t^2)\) of degree
at most four. This checks the exact least degree. The example
rules out a universal purity test using only the first nineteen
moments even inside this box.

## 5. Exactly which odd primes the hyperbolic family omits

**Corollary HyperbolicPrimeOmissions.** For an odd prime \(p\),
the odd-height hyperbolic family omits degree \(p\) exactly when
\(2p-1\) is prime or \(2p-1=qr\) with primes \(q,r\equiv3\bmod4\),
allowing \(q=r\).

The [hyperbolic factor criterion](HYPERBOLIC_ODD_HEIGHT_IMAGE_WITNESSES.md)
expresses \(2p-1\) as the product of two integers greater than one,
both \(1\bmod4\). A proper prime divisor \(1\bmod4\) gives exactly
such a factorization, since \(2p-1\equiv1\bmod4\). If there is no such
prime divisor and \(2p-1\) is composite, all its prime factors are
\(3\bmod4\) and their total multiplicity is even. If there are at least
four, group any two and use the remaining factors for the other group.
If there are exactly two, the only proper factors are the
\(3\bmod4\) primes, so the required split is impossible. A prime
\(2p-1\) also has no split. Thus degree five is omitted via \(9=3^2\)
and degree seven via the prime \(13\).

Only the specified quintic coefficient box is classified here. General
quintics and the corresponding degree-seven problem remain open.
