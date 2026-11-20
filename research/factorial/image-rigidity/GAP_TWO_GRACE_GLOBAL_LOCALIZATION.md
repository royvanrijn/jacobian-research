# Global gap-two factorial root localization by two circular domains

## Exact scope and the actual moment polynomial

**Theorem OrdinaryGapTwoGlobalHermiteRectangle.**

Let a>=40000 and N>=2 be integers. Every zero x of the ACTUAL reversed
ordinary factorial polynomial
\[
 A_{N,a}(x)=\sum_{k=0}^{N}\binom Nk(aN+1)_{2k}x^{N-k}
\]
satisfies, on the physical common-N and Hermite scales
\[
 w=x/N^2,\qquad Z=(w+a^2)/(2a^{3/2}),
 \qquad
 |\Re Z|<1/20,\quad |\Im Z|<4.                         \tag{G1}
\]
This is a bound for EVERY complex zero at EVERY stated finite order.
It does not assert that consecutive polynomials are coprime, that the
imaginary range can be replaced by the bulk interval, or that an Airy
cycle has been identified.

For f=X^a(A+BX^2), B!=0, its Nth ordinary factorial moment is
B^N(aN)! A_{N,a}(A/B). Thus x=A/B is the actual coefficient ratio;
the scale in(G1) is exactly the one used by the retained Hermite bulk
and coalescing-saddle notes. Both coefficients nonzero and lower/upper
orders retain their own physical shapes. No fixed-shape replacement
is made.

## 1. A self-contained closed-disc convolution lemma

For degree-n polynomials in the binomial basis put
\[
 p(z)=\sum_{k=0}^{n}\binom nk a_kz^k,\quad
 q(z)=\sum_{k=0}^{n}\binom nk b_kz^k,\quad
 (p\diamond q)(z)=\sum_{k=0}^{n}\binom nk a_kb_kz^k.
\]

**Lemma ClosedDiscGraceComposition.**

Suppose every zero of p lies in a closed disc D(c,R), with R>0,
and q(0)!=0. For every zero zeta of p diamond q there are alpha in that
disc and a zero beta of q such that zeta=-alpha*beta.

Here is a direct proof of the classical disc form. For a polynomial p of
exact degree n and a point alpha outside its closed root disc, its polar
derivative is D_alpha p(z)=n p(z)+(alpha-z)p'(z). It has exact degree n-1:
its leading coefficient is n times the old leading coefficient times
(alpha-mean(root)), and the mean of the roots is inside the convex disc.

Normalize the disc to the unit disc. If z is outside it, every number
1/(z-r_j) belongs to
\[
 D\left(\frac{\overline z}{|z|^2-1},\frac1{|z|^2-1}\right).             \tag{G2}
\]
This follows by substituting r=z-1/t into |r|<=1 and completing a square.
The disc in(G2) is convex and excludes zero. Its average
u=n^{-1}sum_j1/(z-r_j) also belongs to it. If D_alpha p(z)=0, then
alpha=z-1/u would belong to the unit disc, a contradiction. Therefore
D_alpha p has no zero outside the old disc. Affine covariance gives
the assertion for D(c,R).

Let p-hat be the symmetric multiaffine polarization,
p-hat(z_1,...,z_n)=sum_k a_k e_k(z_1,...,z_n).
The identity p-hat(alpha,z,...,z)=D_alpha p(z)/n follows by
comparing elementary symmetric coefficients. Repeatedly apply normalized
polar derivatives with alpha_1,...,alpha_n
outside the disc. Each step lowers the exact degree by one and keeps
all remaining zeros in the disc. The final constant
p-hat(alpha_1,...,alpha_n) is consequently nonzero.

Write q(z)=q_n product_j(z-beta_j), where every beta_j!=0. Elementary
symmetric functions give the exact identity
\[
 \widehat p(-\zeta/\beta_1,\ldots,-\zeta/\beta_n)
       =(p\diamond q)(\zeta)/q(0).                                   \tag{G3}
\]
At a convolution zero, at least one of -zeta/beta_j must therefore
belong to the disc. This is the required conclusion. In particular
the MINUS sign is essential, as the degree-one case already shows.

This proves the imported classical mechanism directly. Its primary
lineage is Grace--Szego; Hormander's published polar-form argument,
[Math.Scand.2(1954),55-64](https://doi.org/10.7146/math.scand.a-10395),
Lemma2 and Theorem1 on printed57-58, and the
precise closed-disc formulation in Alvarez-Nodarse--Castillo
[arXiv2510.01703v2](https://arxiv.org/html/2510.01703v2),
Theorem3.2, were inspected. The actual1954 PDF credits
LARS HORMANDER; the publisher webpage's different author metadata
is not used. Neither orthogonality on the unit circle nor a
real-rooted free-convolution theorem is a hypothesis of this lemma.
No priority claim is made for the classical composition argument.

## 2. Two discs for each Gamma-Appell factor

For gamma>0 and integer n>=1 define the monic polynomial
\[
 G_{n,\gamma}(x)=\sum_{k=0}^{n}\binom nk(\gamma)_kx^{n-k},
 \qquad h=n-1.
\]

**Lemma GammaAppellTwoDiscs.**

Every zero beta of G satisfies
\[
 -\gamma-h\le\Re\beta\le-\gamma+h,\qquad
 |\Im\beta|\le2\sqrt{h(\gamma+h)},                                    \tag{G4}
\]
and hence belongs to BOTH closed discs
\[
 D(0,\gamma+3h),\qquad D(-2\gamma,\gamma+3h).                          \tag{G5}
\]

We prove the matrix representation instead of invoking orthogonality at
a negative Laguerre parameter. For fixed n,gamma put alpha=-gamma-n.
The n-by-n complex symmetric tridiagonal matrix J has diagonal
2j-1-gamma-n for j=1,...,n, and adjacent off-diagonals
i sqrt[j(gamma+n-j)] for j=1,...,n-1.

Its successive leading characteristic polynomials q_k satisfy
\[
 q_0=1,\quad q_1=x+\gamma+n-1,\quad
 q_k=(x+\gamma+n-2k+1)q_{k-1}
             +(k-1)(\gamma+n-k+1)q_{k-2}.
\]
Coefficient comparison gives
\[
 q_k=\sum_{j=0}^{k}\binom kj(\gamma+n-k)_{k-j}x^j.                    \tag{G6}
\]
For completeness set b=gamma+n-k. After the common factor
(b+1)_{k-j-1} is removed, the coefficient check uses Pascal's
identity and
(k-j)binom(k-1,j-1)=(k-1)binom(k-2,j-1).
The leading coefficient is immediate. At k=n this is G_{n,gamma}.
Equivalently G=(-1)^n n! L_n^{(-gamma-n)}, but no positivity claim
for that Laguerre weight is used.

Write J=B+iC with B,C real symmetric. For an eigenvector v,
beta=(v*Bv+i v*Cv)/(v*v). The diagonal of B lies in
[-gamma-h,-gamma+h]. Each adjacent entry of C has modulus at most
sqrt[h(gamma+h)]; applying 2|v_j v_(j+1)|<=|v_j|^2+|v_(j+1)|^2
bounds |v*Cv|/(v*v) by twice that number. This proves(G4) for this
possibly nonnormal matrix.

Finally both |Re beta| and |Re(beta+2gamma)| are at most gamma+h.
Consequently their squared moduli are at most
(gamma+h)^2+4h(gamma+h)<=(gamma+3h)^2. This proves BOTH discs(G5),
including n=1. There is no condition gamma>=n hidden in the proof.

## 3. Exact physical factorization and transported bounds

Fix n=N and put
\[
 \gamma_1=(aN+1)/2,\quad\gamma_2=(aN+2)/2,\quad h=N-1,\quad
 M=\gamma_1\gamma_2.
\]
The duplication identity gives
\[
 A_{N,a}(x)=4^N\,[G_{N,\gamma_1}\diamond G_{N,\gamma_2}](x/4).          \tag{G7}
\]
The convolution is the UNSIGNED binomial composition from Section1,
also when written in descending monic orientation.

Let zeta=x/4 be any zero of this composition. Apply the closed-disc
lemma first with the origin disc in(G5), and then with the shifted
disc. In the second application zeta lies in a disc with center
2gamma_1 beta and radius(gamma_1+3h)|beta| for some root beta of the
second factor. Using(G4)-(G5), we obtain
\[
 |\zeta|\le M+E,\quad -M-E\le\Re\zeta\le-M+U,                        \tag{G8}
\]
where
\[
 E=3h(\gamma_1+\gamma_2)+9h^2,\quad
 U=(5\gamma_1+3\gamma_2)h+9h^2.
\]
The two uses of the lemma may select DIFFERENT beta. Only the uniform
bounds on their real parts and moduli were used; no common beta or
factorwise pairing of all roots is assumed.

Multiplication by 4/N^2 in(G8), with x=N^2w, gives the useful strict
real bounds
\[
 -12a-36<\Re(w+a^2)<16a+36.                                        \tag{G9}
\]
For example the upper expression before simplification is
16a-19a/N+36-50/N+12/N^2, and the corresponding lower expression is
-12a+9a/N-36+54/N-20/N^2. Both inequalities hold for every N>=2.

We must check U<M BEFORE subtracting squared real bounds. Put t=1/N,
so 0<t<=1/2. Direct expansion gives
\[
\begin{aligned}
 U/N^2&=4a+9-(4a+25/2)t+(7/2)t^2<4a+9,\\
 (E+U)/N^2&=7a+18-(7a+26)t+8t^2<7a+18,\\
 (2M+E)/N^2&=a^2/2+3a+9-(3a/2+27/2)t+(11/2)t^2
                         <a^2/2+3a+9 .
\end{aligned}                                                     \tag{G10}
\]
Since M/N^2>=a^2/4>4a+9 for a>=40000, U<M. Thus(G8) implies
\[
 |\Im\zeta|^2
 \le (M+E)^2-(M-U)^2
 \le (E+U)(2M+E).
\]
After the physical and Hermite rescaling,
\[
 |\Im Z|^2
 <14+\frac{120}{a}+\frac{468}{a^2}+\frac{648}{a^3}<16.                \tag{G11}
\]
The last function decreases with a and its value at40000 is less than16.
Equation(G9) similarly gives
\[
 |\Re Z|<8/\sqrt a+18/a^{3/2}
 \le8/200+18/(40000*200)<1/20.
\]
This proves the theorem on the ENTIRE actual zero set, uniformly in N.

## 4. Consequence and remaining proof obligations

For upper-moment zeros at ordinary gap two and a>=40000, global root
localization is now reduced to the finite Hermite rectangle(G1).
The retained local bulk theorem controls an inner subrectangle;
it does not yet control all of(G1). In particular the intervals between
the bulk boundary and the turning points, the two Airy endpoint zones,
and the zero-free regions beyond the endpoints still require a complete
argument. This localization theorem does NOT place every zero in the
inner bulk tube, and it does NOT close gap-two Strong Factorial.

The exact coalescing geometry and the actual-order action drop provide
candidate bridges for these remaining zones. Arbitrary positive measures,
generic Gamma convolutions, negative-parameter Laguerre orthogonality,
and real-rooted convolution interlacing supply no automatic replacement
for the stated missing estimates.

## 5. Exact equivalence for the remaining consecutive-order question

**Corollary UpperRootRectangleCoprimalityEquivalence.**

Fix integers a>=40000 and n>=1, and put N=n+1. Here the gap two is
the EXPONENT gap in X^a(A+BX^2); the two moment orders are consecutive.
Define the invertible complex affine coordinate
\[
 s_N(Z)=N^2[-a^2+2a^{3/2}Z],\quad
 U(Z)=A_{N,a}(s_N(Z)),\quad L(Z)=A_{n,a}(s_N(Z)).
\]
The same common-N argument s_N is used for BOTH polynomials; the
lower polynomial still has its own degree n and coefficient parameter an.
Let R be the open rectangle |Re Z|<1/20, |Im Z|<4. Then
\[
 \gcd_{\mathbb Q[x]}(A_{N,a},A_{n,a})=1
 \quad\Longleftrightarrow\quad
 \text{there is no }Z\in R\text{ with }U(Z)=L(Z)=0.                 \tag{G12}
\]

Indeed, coprimality over Q is equivalent to the absence of a common
complex root. The affine substitution is bijective. Every upper root
is inside R by the theorem, so every possible common root is inside R.
These observations give both implications, with no condition on the
location of lower roots.

The universal assertions obtained by quantifying (G12) over ALL
integers a>=40000,n>=1 are therefore equivalent. Both the global
consecutive-coprimality assertion and its rectangle noncollision
formulation remain OPEN. The theorem is valid over C, or equivalently
as a common-root exclusion over EVERY characteristic-zero extension
field. A failed rational gcd need not have a root in a specified
non-algebraically-closed field K, so no such fixed-K converse is asserted.
Neither (G12) nor the localization proves either open assertion.

## Exact reproduction and context

Run from the repository root:

    python3 research/factorial/image-rigidity/verify_gap_two_grace_global_localization.py \
      --output /tmp/gap_two_grace_global_localization_review.json

The standard-library verifier checks every principal-minor coefficient
identity in its declared exact fixtures, the actual duplication vectors,
the full-degree polar identity, closed inverse-disc controls and all
cleared rational inequalities/negative controls. The universal arguments
are the written proofs above. The mandatory ledger pins the proof,
checker and exact fixture input. The cached
[Lars Hormander primary scan](sources/hormander_grace_1954.pdf)
is retained only for the explicitly inspected classical context.

Updated 20 November 2026.
