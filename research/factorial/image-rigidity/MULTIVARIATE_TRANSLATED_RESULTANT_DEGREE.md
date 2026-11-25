# Exact resultant degree along every translated factorial ray

Fix \(N=n+1\ge2\). Let a nonempty list of unequal exponent
rows have smaller entries \(A_i\ge0\), gaps \(d_i\ge1\),
and fixed orientations. At every translation \(T\ge0\),
add \(T\) to both entries of every row. Equal rows contribute
scalars and can be discarded. Let \(D_+,D_-\) be the sums
of the increasing and decreasing gaps, respectively, and put
\(D=D_++D_-\), \(\Delta=D_+-D_-\), and \(C=\sum_i d_i^2\).

For \(j=n,N\), divide the factorial moment by the product of the
row factors \((j(T+A_i))!\). The resulting polynomial is
\[
 F_j(T,z)=\sum_{k=0}^j\binom jk z^k
     \prod_i\prod_{h=1}^{d_i\ell_{j,i}(k)}(j(T+A_i)+h),\qquad
 \ell_{j,i}(k)=
 \begin{cases}
   k&\text{in an increasing row},\\
   j-k&\text{in a decreasing row}.
 \end{cases}
\]
Write \(R(T)=\operatorname{Res}_z(F_n(T,z),F_N(T,z))\) in
\(\mathbb Z[T]\).

**Theorem.** The resultant is nonzero and its degree is exactly
\[
 \deg_T R=
 \begin{cases}
    nND,&\Delta\ne0,\\
    nN(D-\tfrac12),&\Delta=0.
 \end{cases}
\]
Thus the adjacent moments are coprime at every nonnegative integer
translation except possibly at a set of translations of
cardinality at most the displayed degree. The statement applies
at every terminal, with the exceptional set depending on the
terminal and the untranslated profile.

## Leading term off the balanced case

Let \(u\) be a formal variable and put
\[
 G_j(u,z)=u^{jD_-}F_j(u^{-1},u^\Delta z).
\]
Indeed \(\deg_T[z^k]F_j=jD_-+k\Delta\), so \(G_j\) has polynomial
coefficients in \(u\). The resultant scaling law gives
\[
       \operatorname{Res}_z(G_n,G_N)=u^{nND}R(u^{-1}).       \tag{1}
\]
The scalar prefactors account for \(2nND_-\) in the exponent,
and the common variable rescaling accounts for \(nN\Delta\).
At \(u=0\) one has
\[
                 G_j(0,z)=j^{jD_-}(1+j^\Delta z)^j .
\]
When \(\Delta\ne0\), these two polynomials have distinct roots,
so the left-hand side of (1) has order zero in \(u\). This gives
the asserted degree in that case.

## The exact Hermite jet in the balanced case

Assume \(\Delta=0\). Write \(j^{-jD_-}[z^k]G_j(u,z)=
\binom jk P_{j,k}(u)\), so that
\[
        P_{j,k}(u)=\prod_i\prod_{h=1}^{d_i\ell_{j,i}(k)}
                             \bigl(1+u(A_i+h/j)\bigr).
\]
For \(r\ge1\), the coefficient of \(u^r\) in its formal logarithm
is a polynomial in \(k\) of degree at most \(r+1\).
In fact it is \((-1)^{r+1}/r\) times a sum of polynomials of
degree \(r\) over a linearly varying range of \(h\). The
coefficient at \(u\) has leading \(k^2\)-coefficient \(C/(2j)\).
Exponentiating, the coefficient of \(u^p\) in \(P_{j,k}\)
has \(k\)-degree at most \(2p\) and \(k^{2p}\)-coefficient
\((C/(2j))^p/p!\). For \(r>1\), \(r+1<2r\), so only the
linear term in the logarithm contributes to that coefficient.

Now put \(u=\lambda^2\) and \(z=-1+\lambda x\). For the
\(z^k\) contribution, the coefficient of \(\lambda^\ell\)
in \((1-\lambda x)^k\) is \((-x)^\ell\binom{k}{\ell}\),
a polynomial in \(k\) of degree at most \(\ell\). The alternating
sum with weights \((-1)^k\binom jk\) kills every polynomial
of degree less than \(j\), and sends \(k^j\) to \((-1)^j j!\).
Consequently
\[
 H_j(\lambda,x)=j^{-jD_-}\lambda^{-j}
                      G_j(\lambda^2,-1+\lambda x)
\]
is polynomial, and its exact constant term in \(\lambda\) is
\[
 H_j(0,x)=j!\sum_{p=0}^{\lfloor j/2\rfloor}
                    \frac{(C/(2j))^p}{p!(j-2p)!}\,x^{j-2p}.
 \tag{2}
\]
Here only terms with \(2p+\ell=j\) contribute, and the signs
\((-1)^j\) and \((-1)^\ell\) cancel.

The monic real polynomial \(i^{-j}H_j(0,ix)\) is the
probabilists' Hermite polynomial with variance \(C/j>0\).
The relative squared scale of the \(N\)-anchor to the
\(n\)-anchor is \(n/N\), within the strict
[scaled Hermite interlacing interval](MULTIVARIATE_SCALED_HERMITE_INTERLACING.md).
Thus the two polynomials in (2) have no common zero.

The same-degree resultant transformation under
\(z=-1+\lambda x\) multiplies by \(\lambda^{nN}\).
The two normalizations dividing by \(\lambda^j\) contribute
\(\lambda^{-2nN}\). Therefore
\[
 \operatorname{Res}_x(H_n,H_N)
 = (nN)^{-nND_-}\lambda^{-nN}
                   \operatorname{Res}_z(G_n(\lambda^2,z),G_N(\lambda^2,z)).
\]
Its left-hand side has a nonzero constant term by (2).
It follows that the resultant on the left of (1) has order
\(nN/2\) in \(u\), an integer since the orders are consecutive.
This gives the balanced degree formula.

At an ordinary nonnegative translation the original normalized
polynomials have both endpoint coefficients nonzero. The resultant
criterion and rational Bézout then establish the characteristic-zero
claim. A finite bound on the number of exceptional translations
does not decide whether any such translation actually occurs.
The [divisor-prime criterion](MULTIVARIATE_TRANSLATED_DIVISOR_PRIME.md)
gives an independent certificate for individual translations on
rays admitting a narrow threshold.
