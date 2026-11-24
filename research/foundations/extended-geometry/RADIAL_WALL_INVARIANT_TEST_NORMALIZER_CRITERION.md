# Radial walls and invariant-test normalizer counterexamples

Let A be a commutative Q-algebra with a normalized Q-linear functional lambda. Extend it to A[t,t^-1] by

    Lambda(sum_n a_n t^n)=lambda(a_0).

A phase F is pure if Lambda(F^m)=0 for every m>=1. A fixed marker is unsafe if Lambda(H F^m) is nonzero for infinitely many m. The theorem below gives one Laurent unit marker nonzero at every positive order. Its geometric conclusions require a separate actual evaluation hypothesis.

## 1. Conditional theorem and exact hypotheses

Let r>=1 and b>=0 be integers. Let Omega_s be finitely many pairwise disjoint, full-dimensional bounded open convex polyhedral cells in R^r. On them let rho_s be nonnegative real polynomial densities of degree at most b, with total integral 1. Every supporting affine hyperplane has real algebraic coefficients. Supporting hyperplanes through zero are allowed. Assume that at least one cell whose closure contains zero has a nonzero density polynomial.

Let W be a finite Galois-stable geometric group of linear matrices on Qbar^r preserving a positive definite rational quadratic form q. The group may be trivial. It need not preserve the individual cells or density. Write

    R=Q[u_1,...,u_r]^W.

Here invariance means invariance after Qbar extension under the stated geometric matrix set. Assume an actual unital Q-algebra homomorphism theta:R->A satisfying, for EVERY P in R,

    lambda(theta(P))=sum_s integral_Omega_s rho_s(u)P(u)du.       (1)

Individual coordinates u_i need not have images in A. Assume compatible scalar extension of A and lambda to every characteristic-zero extension of Q under consideration.

Write rho_s=sum_j rho_sj into homogeneous components about zero. Call a cell Origin if zero belongs to its closure, and Away otherwise. Define

    j0=min{j: rho_sj is nonzero for some Origin cell}.            (2)

There exist a rational W-invariant polynomial f with f(0)=1 and an integer L>=r+b such that the following actual functions work. Put T=1+t and extend theta coefficientwise to R[T]:

    a(T)=theta(f(Tu)),
    p=T*t^-1*a(T)^2,
    F=p^L*(1+p),                 H=t^(r+j0).                    (3)

Then for every m>=1,

    Lambda(F^m)=0,               Lambda(H F^m)>0.               (4)

The full marked values are positive rational numbers. Consequently the same Q-defined phase and fixed unit have zero pure moments and nonzero marked moments after every compatible characteristic-zero extension. Every integer monomial marker t^ell with ell<r is eventually zero against F^m.

If every homogeneous piece of degree below j0 vanishes on EVERY cell, including Away cells, then every ell<r+j0 is eventually zero, and r+j0 is the least unsafe integer monomial exponent. Without this additional global vanishing hypothesis no such least-exponent assertion is made. The theorem does not classify every exponent above r+j0 and does not give a smallest L.

For the geometric strengthening assume A=O(X), where X is a smooth geometrically integral affine Q-variety. Require an actual geometric point xi of X such that

    theta(P)(xi)=P(0) for every P in R.                          (5)

Then F on X times G_m is not a unit times a proper power after any characteristic-zero extension, and its actual Laurent support in t is primitive. These are additional conclusions from (5), not consequences of a formal radial coordinate alone.

## 2. A rational invariant boundary polynomial

Ignore supporting hyperplanes through zero. Normalize the linear equation of each remaining hyperplane to have constant term 1. Its coefficients are algebraic. Take the finite Galois orbit of this finite collection and multiply their normalized equations. The result q_base is rational and has value 1 at zero. Multiplying q_base(wu) over all geometric w in W gives a rational W-invariant polynomial f_base with f_base(0)=1, vanishing on every nonradial supporting hyperplane. The empty product is 1.

Finiteness of the combined orbit follows because W is itself a finite Galois-stable set of algebraic matrices. This is a finite Q-group action condition; it is not a claim for an arbitrary finite algebraic matrix group whose Galois conjugates might generate an infinite group.

Any additional rational W-invariant factor with value 1 at zero preserves these boundary properties. In Section 5 such a factor suppresses the Away cells. No factor is required on a radial face, because its Euler boundary flux is identically zero.

## 3. Cellwise Euler identity and its exact coefficient domain

Fix any rational W-invariant f divisible by f_base with f(0)=1. For k>=1 define real polynomial functions

    J_skj(T)=T^(r+j)*integral_Omega_s rho_sj(u)f(Tu)^(2k)du,
    nu_skj=J_skj(1)=integral_Omega_s rho_sj(u)f(u)^(2k)du.         (6)

Since rho_sj is homogeneous of degree j, differentiation and the ordinary divergence theorem on a polytope give

    J_skj'(T)=T^(r+j-1)*integral_boundaryOmega_s
                 (u dot n)rho_sj(u)f(Tu)^(2k)dS.                (7)

Indeed the derivative of T^(r+j) times the integral equals T^(r+j-1) times the integral of div_u[u rho_sj(u)f(Tu)^(2k)]. All functions are polynomials and all faces are compact; face intersections have surface measure zero.

If a face lies in a hyperplane through zero, u dot n is zero identically on that face for every T. On any other face f(u)=0. The polynomial f(Tu) is then divisible by T-1 on the face, so its 2k-th power is divisible by (T-1)^(2k). Therefore the exact polynomial divisibility statement is

    J_skj(1+t)=nu_skj+O(t^(2k+1)).                              (8)

The O notation here means divisibility by the displayed power in a finite polynomial, not an asymptotic approximation.

Since W acts linearly, each homogeneous component f_d of f is W-invariant. Thus f(Tu)=sum_d T^d f_d(u) belongs to R[T], and a(T)=sum_d T^d theta(f_d) is an actual polynomial in A[T], with a(0)=1. Every coefficient of f(Tu)^(2k) is also invariant. Applying (1) coefficientwise and taking the constant term gives

    Lambda(t^ell p^k)
      =[t^(k-ell)] T^k sum_s integral_Omega_s rho_s(u)f(Tu)^(2k)du
      =sum_(s,j) nu_skj * binom(k-r-j,k-ell),                  (9)

for every integer k>=1 and ell>=-k. To verify the last equality, replace the summand by T^(k-r-j)J_skj(T) and use (8). The requested coefficient has order k-ell<=2k, so the remainder cannot contribute. Negative powers of T are expanded formally at T=1. For an integer upper argument, binom(n,d) denotes [t^d](1+t)^n and is zero when d<0.

This derivation uses only invariant-test moments (1). Individual moments of u_i and individual homogeneous density pieces are never required to be native algebraic inputs.

## 4. Pure tail and lower marker bounds

For k>=r+b every upper argument k-r-j in (9) is a nonnegative integer strictly below k. Taking ell=0 therefore proves

    Lambda(p^k)=0 for every k>=r+b.                            (10)

For any fixed integer ell<r, the lower argument k-ell exceeds every upper argument k-r-j once k is sufficiently large for the domain ell>=-k. Thus

    Lambda(t^ell p^k)=0 for all sufficiently large k.           (11)

If all rho_sj with j<j0 vanish on ALL cells, the same reasoning improves ell<r to ell<r+j0. These are marker statements for the given phase; they do not assert global Mathieu properties of A or of the full Laurent algebra.

## 5. Positive leading term and suppression of Away cells

For an Origin cell and u in that cell, convexity and zero in its closure imply epsilon*u is inside the cell for every 0<epsilon<=1. Divide rho_s(epsilon*u)>=0 by epsilon^j0 and let epsilon tend to zero. All lower pieces on Origin cells vanish by definition, so

    rho_s,j0(u)>=0 on every Origin cell.                      (12)

At least one such polynomial is nonzero, hence strictly positive on a nonempty open subset.

All Away cell closures have positive distance from zero. Choose rational c>0 small enough that 0<1-c*q<=1 on all cell closures. If there are Away cells, put

    B=max_allclosures |f_base|,
    gamma=max_Awayclosures(1-c*q)<1.

Choose an integer M large enough that B*gamma^M<1/4. If there are no Away cells use M=0. Set

    f=f_base*(1-c*q)^M.                                      (13)

It is rational and W-invariant, retains all boundary zeros and has |f|<1/4 on Away cells. Only AFTER fixing M choose a sufficiently small closed ball inside a positive Origin cell and sufficiently near zero that |f|>1/2 and rho_s,j0>0 throughout that ball. Such a ball exists along a cone of directions on which the leading homogeneous polynomial is positive. If A_plus is its positive rho_s,j0 integral, the Origin contribution at degree j0 is at least A_plus*4^-k for every k>=1.

For k>=r+b and ell=r+j0, formula (9) reduces to

    Lambda(t^(r+j0)p^k)
      =sum_Origin nu_s,k,j0
       +sum_(Away,j<=j0) nu_skj*binom(k-r-j,j0-j).             (14)

Origin terms of smaller degree vanish. Higher-degree terms have lower binomial argument exceeding the upper argument and vanish. Polynomial pieces on the bounded Away cells are integrable. Consequently the second sum has absolute value at most

    C*(1+k)^j0*16^-k                                         (15)

for a finite C independent of k. The complete marked value is therefore at least

    4^-k*[A_plus-C*(1+k)^j0*4^-k]>0                           (16)

for every k>=K for some finite integer K. This is an all-later-orders argument from exponential decay, not a finite moment-prefix inference. With arbitrary real density coefficients no effective algorithm for choosing K is asserted.

## 6. One fixed unit at every order and the exact scope of sharpness

Choose L>=max(r+b,K). Expand

    F^m=sum_(h=0)^m binom(m,h)p^(Lm+h).                       (17)

Every exponent is in both tails (10) and (16). Thus the pure moment is zero and the complete marked moment with the SAME unit t^(r+j0) is positive for every m>=1. Since the phase and marker are native Q-functions, their full marked moments are rational; nonzero rational values remain nonzero after every compatible characteristic-zero extension.

For any fixed ell<r, all exponents in (17) eventually satisfy (11), proving eventual vanishing for that marker. Under the additional global homogeneous vanishing condition, the same holds for every ell<r+j0, while (4) proves r+j0 unsafe at every order. Thus that extra hypothesis gives the stated least unsafe integer exponent.

Away-cell pieces of smaller degree can contribute to smaller markers when they do not vanish globally. Positivity of the first Origin term alone does not exclude those small nonzero contributions. No stronger sharpness conclusion is used.

## 7. Actual simple divisor and primitive Laurent support

Assume (5) and the geometric hypotheses of Section 1. Evaluation at xi gives a(T)=f(0)=1 for every T, hence

    p(xi,t)=1+t^-1,
    F(xi,t)=(1+t^-1)^L*(2+t^-1).                             (18)

At t=-1/2 the modifier 1+p vanishes, its t-derivative is -4, and p=-1 is a unit. Smoothness makes this a transverse local parameter at the actual point of X_Qbar times G_m. Some geometric prime divisor through it has valuation one in F. A unit times a proper power has valuation divisible by that power at every prime, so such a representation is impossible. Replaying the same actual smooth point after algebraically closed extension proves the absolute statement over every characteristic-zero extension.

The coefficients at exponents 0 and -1 in (18) are 2 and 2L+1. They are nonzero, so the corresponding actual coefficients of F cannot vanish. Its Laurent support contains adjacent integers and therefore has differences generating Z. This is literal support primitivity, not an assertion about an arbitrary grading or an arbitrary quotient phase.

## 8. Method credit and boundary

The interval-Laurent normalizer circuit has Long's antecedents; the polynomial-density and cellwise-Euler variants are recorded in the standard compact-product and polyhedral-density proofs. The argument here writes the required flux, coefficient domain, leading-term positivity and actual invariant-test descent in full. It imports no all-order moment checker, finite-prefix certificate, compact-group model classification or nonabelian moment law.

A group application must separately supply the full-dimensional radial law, a nonzero Origin cell, its actual rational invariant-test homomorphism, scalar-compatible native trace and actual origin evaluation. Root walls through zero are harmless for (7), but lower-dimensional radial support, wall mass or absence of an actual zero fiber is not silently removed. The theorem does not claim these hypotheses for every compact homogeneous space.
