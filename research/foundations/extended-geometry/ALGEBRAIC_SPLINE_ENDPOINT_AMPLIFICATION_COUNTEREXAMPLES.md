# Rational endpoint amplification for algebraic spline laws

## 1. Conditional theorem and actual functions

Fix an integer D>=1. Let rho be a nonnegative compactly supported probability density on R, polynomial of degree at most D-1 on the intervals of a finite partition with distinct REAL ALGEBRAIC knots a_j. Its zero extension is required to be globally C^(D-2) for D>=2. For D=1 no continuity is required. Support always means the support of the associated measure; changing a density at isolated points is immaterial. Zero and inactive partition knots may be inserted.

Let A be a commutative Q-algebra with a normalized Q-linear functional lambda:A->Q and an ACTUAL element z in A satisfying

    lambda(P(z)) = integral_R rho(u)P(u)du    for every P in Q[u].       (1.1)

On A[t,t^-1] let Lambda be lambda followed by the constant term in t. A phase F is pure if Lambda(F^m)=0 for every integer m>=1. A fixed marker H survives every order if Lambda(HF^m)!=0 for every m>=1.

**Theorem 1.** Under these hypotheses there are an EVEN polynomial f in Q[u] with f(0)=1 and an integer L>=max(D,4) such that the ACTUAL Laurent functions

    T=1+t,    p=T*t^-2*f(Tz)^2,
    F=p^L(1+t^L*p),    H=t^-D                                  (1.2)

satisfy Lambda(F^m)=0 and Lambda(HF^m)>0 in Q for EVERY integer m>=1. In particular the SAME fixed UNIT marker H survives every order. One may choose

    f(u)=f0(u)(1+u^2)^M,    M a nonnegative integer,             (1.3)

where f0(u)=v0(u^2) and v0 is the normalized squarefree product of the distinct rational minimal polynomials of all nonzero squared partition knots a_j^2. Thus f0 and f are EVEN and every nonzero real partition knot is a SIMPLE zero of f. Given exact algebraic descriptions of the knots and the first finitely many rational moments specified in Section5, finite algebraic tests choose M and a sufficient L.

Whenever actual scalar extensions have compatible functionals lambda_K=lambda tensor_Q K, the same functions give the identical zeros and nonzero rational marked values over EVERY characteristic-zero field K. Positivity here refers to the rational numbers before scalar extension; no ordering of K is assumed.

For the additional geometric conclusions, assume X=Spec A is a smooth geometrically integral affine Q-variety and the ACTUAL fiber z=0 has a geometric point. Then F is absolutely not a unit times a proper power: over every algebraically closed extension it has no expression U*V^e with U a unit and e>=2. Its actual Laurent support contains adjacent exponents and is radially primitive, meaning that the gcd of all support differences is1. These extra conclusions require the stated actual geometry.

The exponent L depends on the law and the amplifier. No sharp minimal L, arbitrary model law, or global Mathieu ascent is asserted.

## 2. Finite integration by parts and the exact endpoint coefficient

Let

    c_j=rho^(D-1)(a_j+)-rho^(D-1)(a_j-),

using exterior derivatives zero. For D=1 these are jumps of the piecewise constant density. Continuity of every derivative through D-2, and degree at most D-1 on each piece, give the distributional identity

    rho^(D)=sum_j c_j delta_(a_j).

Equivalently, finite repeated integration by parts gives for every polynomial P and every polynomial G with G^(D)=P

    integral rho(u)P(u)du = (-1)^D sum_j c_j G(a_j).            (2.1)

Indeed lower-derivative boundary terms cancel at each knot, the exterior terms vanish by compact support, and the remaining terms are exactly the jumps displayed above. This also proves (2.1) for polynomials without compact support, by doing the finite integrations over a compact interval containing the support.

For k>=1 let G_k be the polynomial D-fold antiderivative of f(u)^(2k) with G_k^(i)(0)=0 for 0<=i<D. Applying (2.1) coefficientwise to f(Tu)^(2k) gives the polynomial identity

    lambda(p^k)=(-1)^D t^-2k T^(k-D) sum_j c_j G_k(Ta_j).       (2.2)

The apparent negative power of T when k<D cancels: every G_k(Ta_j) is divisible by T^D. All subsequent coefficient arguments use k>=D.

At a nonzero knot with f(a_j)=0, the Taylor formula at a_j gives

    G_k((1+t)a_j)=Q_(k,j)(t)
      + a_j^(2k+D) f'(a_j)^(2k) (2k)!/(2k+D)! * t^(2k+D)
      + terms of degree at least 2k+D+1,                      (2.3)

where Q_(k,j) has degree at most D-1. Derivatives numbered D through 2k+D-1 vanish there, while derivative 2k+D is (2k)! f'(a_j)^(2k). The formula remains valid for a multiple zero, with displayed coefficient zero. At a_j=0 the term G_k(Ta_j)=G_k(0)=0 and f needs no zero.

For k>=D, after multiplication by T^(k-D) the baseline in (2.3) has degree at most k-1 and the remainder starts in degree2k+D. Consequently

    Lambda(t^-D p^k)=(2k)!/(2k+D)! sum_j B_j beta_j^k,
    B_j=(-1)^D c_j a_j^D,    beta_j=(a_j f'(a_j))^2.           (2.4)

This identity is exact for every k>=D. Zero knots and knots with c_j=0 contribute nothing.

## 3. Outer support endpoints have positive weights

Let J={j:c_j!=0 and a_j!=0}. It is nonempty: a nonzero compact probability density has an interval of support, and its extreme endpoints have nonzero highest-derivative jumps, with at least one endpoint nonzero. Put

    R=max_(j in J)|a_j|>0.

Every active knot lies in the measure support. An active knot equal to R, if present, is its rightmost endpoint; an active knot equal to -R, if present, is its leftmost endpoint. The farthest active knots cannot be interior to the support, because the corresponding further extreme support endpoint would itself be active and have larger absolute value.

On the final nonzero interval before a positive right endpoint R, continuity through D-2 and matching to exterior zero force

    rho(u)=C(R-u)^(D-1),    C>0.

Thus c_R=(-1)^D C(D-1)! and B_R=C(D-1)!R^D>0. At a negative left endpoint -R the first piece is C(u+R)^(D-1), so c_-R=C(D-1)! and again B_-R=C(D-1)!R^D>0. This includes D=1. There may be one or two outer active knots; EVERY one has positive B. Inactive partition points outside the support are ignored in these comparisons.

## 4. Positive rational amplification and eventual strict positivity

Let v0(s) be the product of the DISTINCT monic minimal polynomials over Q of all nonzero squared partition knots a_j^2, divided by its nonzero constant coefficient, and put f0(u)=v0(u^2). Characteristic zero makes v0 squarefree, v0(0)=1, and v0'(a_j^2)!=0. Hence f0 is EVEN, f0(0)=1, and f0'(a_j)=2a_j v0'(a_j^2)!=0 at every nonzero knot. Additional negative or Galois roots cause no problem; a pair of actual knots +a,-a still has SIMPLE zeros. No symmetry of rho is assumed.

For each integer M>=0 set

    f_M(u)=f0(u)(1+u^2)^M,
    b_j=(a_j f0'(a_j))^2>0,
    beta_j(M)=b_j(1+a_j^2)^(2M),    j in J.                  (4.1)

The amplifier has rational coefficients, value1 at0 and is positive at EVERY real number. Hence f_M retains simple zeros at all nonzero actual knots, without assigning different multiplicities to Galois conjugates.

Choose an outer index j0 with maximal b_j among |a_j|=R. For every interior active knot,

    beta_j(M)/beta_j0(M)
      =(b_j/b_j0)*[(1+a_j^2)/(1+R^2)]^(2M) ->0.              (4.2)

A sufficiently large M makes every one of these finitely many ratios less than1. With that M fixed, put beta_max=max_(j in J)beta_j(M). Only outer endpoints attain this maximum, so

    A_max=sum_(j:beta_j(M)=beta_max) B_j>0.                  (4.3)

The maximum may have two indices; their weights add positively. Therefore

    sum_j B_j beta_j(M)^k / beta_max^k
      =A_max+sum_(beta_j<beta_max)B_j(beta_j/beta_max)^k
      ->A_max>0.                                          (4.4)

There is an integer K>=1 for which the signed endpoint sum is positive for EVERY k>=K. This conclusion follows from the geometric-ratio bound in (4.4), with no finite-prefix inference. Fix f=f_M for the rest of the proof.

## 5. Algebraicity of the jumps and finite choices of the parameters

The real jump coefficients need not be assumed rational or algebraic separately. Native rationality in (1.1) forces their algebraicity. List the distinct NONZERO partition knots as a_1,...,a_r, including inactive knots if any. Necessarily r>=1. Apply (2.1) to P(u)=u^n and G(u)=n! u^(n+D)/(n+D)!, for n=0,...,r-1. The zero knot contributes zero. Thus

    sum_(j=1)^r (c_j a_j^D) a_j^n
      =(-1)^D (n+D)!/n! * lambda(z^n),    0<=n<r.           (5.1)

The Vandermonde matrix (a_j^n) is invertible over Qbar. The right side is rational, so every c_j a_j^D and every nonzero-knot c_j is algebraic. If0 is a partition knot, its coefficient is algebraic as well, because

    c_0=-sum_(a_j!=0)c_j.                                 (5.2)

Here sum_j c_j=0 follows by evaluating the compactly supported derivative rho^(D) against a function identically1 near its support, or directly by telescoping jumps. Thus every B_j in (2.4) is real algebraic, although it need not be rational.

Given exact algebraic knot descriptions and lambda(1),...,lambda(z^(r-1)), solve (5.1), identify the active indices and perform the finite comparisons in (4.2). Increasing M until all required strict algebraic inequalities hold terminates. For this M all beta_j and A_max are algebraic. Increasing K>=1 until

    sum_(beta_j<beta_max)|B_j| (beta_j/beta_max)^K < A_max    (5.3)

also terminates. Every ratio is in (0,1); the empty sum is allowed. The bound is monotone in K and proves positivity for all later k. These are exact finite algebraic tests using the proven unrestricted endpoint identity, not tests of finitely many phase moments. No algorithm for extracting a density from an arbitrary proposed geometric model is assumed or supplied.

## 6. One fixed unit marker at every order

Choose an integer L>=max(D,K,4). Expand the ACTUAL phase (1.2):

    F^m=sum_(h=0)^m binom(m,h)t^(Lh)p^(Lm+h).               (6.1)

For each term put k=Lm+h. Before the factor t^-2k in (2.2), its baseline has degree at most k-1 and its remainder begins in degree2k+D.

For the unmarked term the requested degree is2k-Lh. It exceeds k-1 because

    k-Lh=L(m-h)+h>=1,

and lies below2k+D. Every unmarked term is zero. For the marked term the requested degree is2k-Lh+D, again greater than k-1. If h>=1 it is at most2k, since L>=D, and hence its remainder contributes zero. Only h=0 survives. Formula (2.4) yields for EVERY m>=1

    Lambda(F^m)=0,
    Lambda(HF^m)=(2Lm)!/(2Lm+D)! sum_j B_j beta_j^(Lm)>0.    (6.2)

The complete marked value is rational because F,H and the actual functional are rational. The individual algebraic jumps, knots or endpoint summands need not be rational. Compatible characteristic-zero scalar extensions preserve these exact zeros and nonzeros.

## 7. Actual simple divisor and primitive Laurent support

Assume the extra geometry from Section1. Choose an ACTUAL point x0 over Qbar with z(x0)=0; such a point exists because the nonempty fiber is of finite type over Q. Then f(0)=1 gives p(x0,t)=T/t^2. The modifying factor specializes to

    Q_L(t)=1+t^(L-2)(1+t).                                (7.1)

It has no root0 or-1. Its derivative is t^(L-3)[(L-2)+(L-1)t]. At its only possible nonzero repeated root t0=-(L-2)/(L-1),

    Q_L(t0)=1+(-1)^(L-2)(L-2)^(L-2)/(L-1)^(L-1)!=0,        (7.2)

because the added term has absolute value less than1. Hence Q_L is squarefree and has a root alpha in Qbar with alpha!=0,-1 and Q_L'(alpha)!=0.

At (x0,alpha), p is a unit and 1+t^L p has nonzero t-derivative. Smoothness gives an actual geometric prime divisor through that point where its order is1 and p has order0. Thus ord(F)=1 there. A unit has valuation0, whereas the valuation of a proper e-th power is divisible by e>=2. This proves absolute exclusion even of a UNIT times a proper power. The divisor and derivative argument remains valid over any algebraically closed extension.

The same actual point gives

    F(x0,t)=t^-2L(1+t)^L+t^(-L-2)(1+t)^(L+1).              (7.3)

For L>=4 the coefficients at exponents -2L and -2L+1 are1 and L. Actual coefficient functions are therefore nonzero at adjacent exponents. Their support differences have gcd1, which excludes pullback through t->t^e, e>1, even after multiplication by a Laurent monomial unit.

## 8. Mathematical scope

The construction uses the supplied actual rational functional, algebraic spline law and compatible scalar extensions. Smooth geometric integrality and an actual common zero supply the stronger divisor and support conclusions separately. A law for an abstract measure does not by itself supply any coordinate algebra, orbit or rational model. A geometric application must establish those hypotheses independently, including the GLOBAL regularity of the zero-extended density.

The proof consists of finite repeated integration by parts, rational minimal polynomials, Vandermonde linear algebra and a geometric-ratio bound, followed by direct Laurent coefficient and divisor arguments. It needs no localization formula, Hamiltonian circle hypothesis or external stationary-phase theorem. In particular this conditional result imports no global Section4 assertion from Duistermaat--Heckman. It gives a sufficient finite choice of L, with no claim of a sharp smallest exponent or a converse purity-under-roots principle. The unrestricted all-order conclusions are proved by (2.2)--(6.2).
