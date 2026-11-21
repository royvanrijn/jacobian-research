# Gamma ordered resultants: a closed positive law class and an exact obstruction

This note proves positive polynomial-Gamma tilt closure under translated zero bias and an exact consecutive Appell collision within that enlarged class. The pure Gamma-power/product-Gamma ordered-resultant problem remains unresolved.

Write

    P_m(x)=E(x+Y)^m,  mu=EY,  V=Var(Y)>0.

For every real moment law with finite moments of every order and V>0, at every integer m>=2 coefficient subtraction gives

    P_m-(x+mu)P_(m-1)=(m-1)V B_(m-2),
    B_l(x)=sum_j binom(l,j) nu_j x^(l-j),
    nu_j=[E Y^(j+2)-mu E Y^(j+1)]/[(j+1)V].             (Z1)

## 1. The probability transform is classical zero bias

Goldstein--Reinert, *Stein's method and the zero bias transformation with application to simple random sampling*, Ann. Appl. Probab. 7 (1997), 935--952, DOI 10.1214/aoap/1043862419, Definition 1.1 and Lemma 2.1, establish the centered zero-bias law and its derivative identity. The actually read author PDF is dated 28 November 2001, printed pp.2--4: [retained author manuscript](sources/goldstein_reinert_zero_bias_author2001.pdf) ([official author URL](https://dornsife.usc.edu/larry-goldstein/wp-content/uploads/sites/221/2023/06/zer.pdf)). Applying it to Y-mu and translating back gives

    E[(Y-mu)f(Y)]=V E[f'(Y*)].                          (Z2)

Taking f(y)=y^(j+1) identifies the moments of Y* with nu_j in (Z1). Thus the law-level interpretation of the first Euclidean remainder is zero bias, not a newly introduced probability transformation. Its positivity for general laws says nothing by itself about complex-root separation.

## 2. Positive polynomial Gamma tilts are closed at the law level

**Theorem PositivePolynomialGammaTiltZeroBiasClosure.** The translated zero-bias law of X^d under every allowed positive polynomial Gamma tilt is again X^d under a positive polynomial Gamma tilt. The rescaled coefficients in (Z5) are coefficientwise nonnegative polynomials in the shape and original tilt coefficients. The following construction proves these statements.

Fix gamma>0, a positive integer d, and a nonzero polynomial

    w(z)=sum_(r=0)^K w_r z^r,   w_r>=0, w_K>0.

Let X have rate-one Gamma(gamma) law before tilting, and let the tilted expectation be E_w f=E[w(X)f(X)]/A, with A=Ew(X). Put Y=X^d under this tilted law, mu=E_w Y, V=Var_w(Y)>0. All moments exist, and the density is positive on (0,infinity).

Set a_r=w_r(gamma)_r, R_r=(gamma+r)_d, A=sum a_r, B=sum a_r R_r, so mu=B/A. For an integer t define A_t=sum_(0<=r<=min(t,K)) a_r and B_t=sum_(0<=r<=min(t,K)) a_r R_r, with both sums zero if t<0. Define

    N(z)=sum_(j=0)^(K+d-1) n_j z^j,
    n_j=[mu A_j-B_(j-d)]/(gamma)_(j+1).                (Z3)

Adjacent coefficients give the exact polynomial identity

    z N'(z)+(gamma-z)N(z)=(mu-z^d)w(z).                (Z4)

The omitted coefficients vanish: for j>=K+d the numerator in (Z3) is mu A-B=0. The top coefficient is n_(K+d-1)=w_K, as also follows from (Z4).

Every coefficient n_j is nonnegative. Indeed R_r strictly increases in r, so the initial-segment weighted mean B_u/A_u is at most the full mean B/A whenever A_u>0. Because A_j>=A_(j-d),

    mu A_j-B_(j-d)>=mu A_(j-d)-B_(j-d)>=0.

There is also a coefficientwise positive version of this calculation. With u=j-d,

    B A_j-A B_u
       =B(A_j-A_u)+sum_(s<=u<r) a_r a_s (R_r-R_s).     (Z5)

The sum ranges over 0<=r,s<=K. Each R_r-R_s has nonnegative coefficients as a polynomial in gamma, since it is the positive integer translate difference of a rising factorial. Therefore the coefficients of

    Nhat=(gamma)_(K+d) A N

are polynomials with nonnegative coefficients in gamma,w_0,...,w_K. In (Z3), the remaining multiplier is (gamma+j+1)_(K+d-j-1), also coefficientwise nonnegative. This rescaling changes no normalized law.

Integration by parts in (Z4), with vanishing endpoints because gamma>0 and exponential decay, gives

    Cov_w(X^d,f(X))=E[X N(X)f'(X)]/A.                 (Z6)

For f(X)=X^(d(j+1)), this yields

    nu_j=d E[X^(d(j+1))N(X)]/(A V).

The case j=0 gives E[X^d N(X)]=A V/d. Hence Y* is exactly X^d under the new Gamma tilt

    w_new(X)=X^d N(X),                                (Z7)

or, equivalently, X^d Nhat(X). This identifies the LAW, not just its moments:
in (Z6) take f(X)=F(X^d), with F a continuously differentiable compactly
supported function. The chain rule and the displayed normalization give
E_w[(Y-mu)F(Y)]=V E_new[F'(X^d)]. The boundary at zero vanishes because
gamma>0 and N is finite, and the boundary at infinity vanishes by exponential
decay. Thus the classical zero-bias uniqueness statement applies even when
the moment sequence of X^d is indeterminate.

This is an enlarged class closed under the first remainder-law transform for every positive integer d. For w=1 the polynomial is N(z)=sum_(j=0)^(d-1) (gamma+j+1)_(d-j-1) z^j. The pure finite product-Gamma moment class itself already fails closure for d=2. Indeed N(z)=z+gamma+1 and the new tilt is z^2(z+gamma+1). Its moments are

    nu_j=(gamma+2)_(2j) (2gamma+2j+3)/(2gamma+3).

Consequently

    nu_(j+1)/nu_j
      =(gamma+2j+2)(gamma+2j+3)(2gamma+2j+5)/(2gamma+2j+3).

At the denominator zero j=-gamma-3/2, the three numerator factors are
-gamma-1, -gamma and2, all nonzero for gamma>0. This rational function is
not a polynomial. A finite product of independent positive-shape Gamma
variables, with any fixed positive overall scale, has a POLYNOMIAL
successive moment ratio c*product_i(b_i+j). Equality at all integers j>=0
would force equality of rational functions, which is impossible here.
This proves the class nonclosure without a separate scratch premise.

This law closure is not a closed Euclidean recurrence for the desired pair. After one step the pair is P_(m-1)(w) and P_(m-2)(w_new), with different laws; (Z1) cannot simply be reapplied as if they were consecutive degrees of the same law.

## 3. A finite positive polynomial tilt already has a consecutive collision

**Theorem PositivePolynomialGammaTiltConsecutiveCollision.** The enlarged positive polynomial-tilt class has a consecutive Appell collision already in degrees4 and5, with the literal rate-one law below. The underlying two-Gamma collision is the RETAINED example in [GAMMA_CONVOLUTION_APPELL_OBSTRUCTION.tex](GAMMA_CONVOLUTION_APPELL_OBSTRUCTION.tex), Theorem `thm:finite-obstruction`, scaled by2; it is not rediscovered here. The additional statement is its exact realization by a finite positive polynomial tilt. Let gamma=5014/361, and take a mixture of rate-one Gamma laws with shapes gamma+n. The following are its six positive rational mixing weights:

| n | pi_n |
|---|---|
| 0 | 23959174167747319/24009255465548716 |
| 1 | 787871138391581/404852698201373433 |
| 24 | 337387137375/4512464765741536 |
| 25 | 578110771745/9270172181795112 |
| 88 | 14673828944125/11828739703970562912 |
| 89 | 70917583625/48018510931097432 |

They sum to one. Relative to rate-one Gamma(gamma), this is the literal degree-89 polynomial tilt

    w(z)=sum_n pi_n z^n/(gamma)_n,

with nonnegative coefficients and Ew=1. Its first six moments are

    1, 264/19, 74750/361, 22622000/6859,
    7321626000/130321, 2579379400000/2476099.           (Z8)

They are exactly the moments through degree five of Z=G_1+19G_2, with independent rate-one Gamma shapes 125/9 and 1/3249. This is twice the retained Gamma-convolution counterexample. No distributional equality beyond the displayed moments is needed or claimed.

Put u=x+264/19. Exact coefficient multiplication gives

    P_4(x)=(u^2-2u+76)(u^2+2u+12),
    P_5(x)=(u+4)(u^2-2u+76)^2.                        (Z9)

The shared quadratic has discriminant -300, and V=14. Even the first Euclidean remainder visibly retains it:

    P_5-u P_4=56 (u^2-2u+76)(u+38/7).                 (Z10)

Thus a closed positive Gamma-tilt law class does not imply nonzero consecutive Appell resultants. This is not a pure Gamma-power or product-Gamma counterexample. It rules out an induction that forgets the special coupled ancestry of the two laws. The coefficientwise positivity of the tilt transform in (Z5) remains valid.

The class distinction can also be checked intrinsically. This finite polynomial tilt has a finite exponential moment for every real argument less than one. A positive product of two or more independent Gamma variables, or X^d with d>=2, has no positive exponential moment. A single scaled Gamma law would satisfy 3 kappa_3^2=2 kappa_2 kappa_4, whereas here (kappa_2,kappa_3,kappa_4)=(14,32,324), which fails that equality. Thus the witness is genuinely outside all those pure classes.

The rational weights were found by one exact linear feasibility calculation matching the five factorial moments of the relevant negative-binomial mixing index. Discovery used SymPy; verification uses only Fraction arithmetic. The standalone replay `verify_gamma_polynomial_tilt_zero_bias.py` rechecks every weight, (Z8)--(Z10), all 91 coefficients of the degree-89 instance of (Z4), its zero-bias normalization and four remainder moments. These finite checks support the exact witness; the universal assertions (Z3)--(Z7) are proved above, not inferred from those checks.

## 4. A second exact barrier: no fixed-band Appell recurrence

For a pure product-Gamma law, c_n=product_i(b_i)_n. The identity (b)_n=sum_(sigma in S_n) b^(number of cycles of sigma) expresses c_n as the weight of d-tuples of permutations. Decompose such a tuple into joint orbits of the generated permutation group. The formal labeled-component identity shows that its n-th cumulant kappa_n is the sum of the positive monomial weights of transitive tuples. This is a formal-series argument, valid even when the moment-generating series has zero analytic radius.

For every n>=1 there is a transitive tuple: take an n-cycle in the first coordinate and identity permutations in all other coordinates. Consequently kappa_n>0 for all positive parameters, in every order. Differentiating the formal Appell generating series gives

    x P_n=P_(n+1)-sum_(j=0)^n binom(n,j)kappa_(j+1)P_(n-j).  (Z11)

In particular the unique coefficient of P_0 is -kappa_(n+1), never zero. There is no recurrence of fixed finite width in this Appell basis. This excludes using a fixed finite-band d-orthogonal recurrence for these same fixed-law polynomials. It does not exclude a multi-index family, a degree-dependent operator, a structured subresultant recurrence, or another basis.

## 5. Exact remaining target

The all-order ordered Gamma-power resultant remains open here. A useful next sufficient lemma must retain more than positivity of each moment law or of the zero-bias tilt: it must control the particular cross-resultants of the coupled laws generated from w=1, or produce a different signed subresultant/combinatorial expansion. A positive orthogonal-polynomial production matrix for an associated moment functional is not such a proof, since its characteristic polynomials differ from these Appell polynomials.

The only external primary statement used for probability-law identification
is Goldstein--Reinert Definition1.1/Lemma2.1, in the retained author manuscript
specified in Section1. Its definition, existence and uniqueness are classical.
The closure formulas, the coefficientwise positivity, the exact degree89
positive-tilt realization and the formal no-fixed-band consequence are proved
above. Bounded exact verification is separate from those universal proofs.
No all-order product-Gamma separation result, source-absence priority claim,
external review or complete formal verification is asserted.
