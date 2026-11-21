# Noncentral constant affine cubic shapes: purity, safety and the exact degree spectrum

Let K be any characteristic-zero field. Let lambda in K be CONSTANT and

    lambda !=1/2.                                      (H1)

Here "noncentral" means only that the third root of the normalized cubic
is not the midpoint of the two selected interval-endpoint roots0 and1.
It does NOT mean that the underlying cubic lacks a symmetry: lambda=-1
and lambda=2 are included, though their depressed cubic is symmetric.

Put

    Phi_lambda(v)=v(v-1)(v-lambda),
    kappa_m(lambda)=integral_0^1 Phi_lambda(v)^m dv.

All integrals below are rational polynomial functionals. Suppose
E in K[x], E(0)!=0, and s,h in K[[x]] satisfy

    s(0)=-1, h(0)=2,
    P(x,t)=E(x)h(x)^3 Phi_lambda((t-s(x))/h(x))
          =C(x)+tB(x)+t^2D_t(x)+t^3E(x) in K[x,t].      (H2)

On A=K[x,y,t]/(t^2+2xy-1), set f=P/x. This is a regular element:
P(0,t)=E(0)(t+1)(t-1)(t+1-2lambda). Use

    Lambda(g)=CT_x (1/2)integral_-1^1
                        g(x,(1-t^2)/(2x),t)dt.

Purity means Lambda(f^m)=0 at every positive integer m. Safety means
that every fixed g in A has Lambda(g f^m)=0 for all sufficiently large m.

**Theorem NoncentralConstantAffineCubicClassification.**
Under(H1)--(H2), put a=E h^3 and ell=h/2. Then

    f is pure iff ell=1-x a'/a.                        (H3)

Every pure member has rational s,h and a unique expression

    h=2/Q, Q in K[x], Q(0)=1, Q squarefree.              (H4)

It is safe iff Q is constant. If Q=1, then E is constant and the
phase is a polynomial Hopf shear of a strictly negative-angular-weight
phase. Its possible least degrees are exactly2,5,8,... .

If q=degQ>=1, then degE=3q+1. EVERY nonconstant polynomial psi in K[x]
has Lambda(psi f^m)!=0 at infinitely many positive orders. The possible
least degrees in this unsafe branch are exactly6,9,12,... .
The theorem does not assert nonvanishing at EVERY order.

All assertions hold after every coefficient-field extension. There is
no hypothesis that all kappa_m are nonzero. The proof allows isolated
or sparsely occurring scalar moment zeros. No claim is made that
purity forces a cubic to have constant affine shape; no unrestricted
mixed-cubic or global original Image degree5/7 classification follows.

## 1. Regularity and the exact ordinary degree

The two regularity conditions are C0+D_t0=0 and B0+E0=0.
A literal polynomial representative of f is

    2C0*y+2B0*y*t+(C-C0)/x+t(B-B0)/x
        +t^2(D_t+C0)/x+t^3(E+B0)/x.

It has no xy monomial. Under a degree-compatible order with leading
relation monomial xy it is reduced. Its nonzero top homogeneous part
cannot be divisible by t^2+2xy, whose leading monomial is xy.
Consequently its least polynomial degree in A is exactly

    d_f=max(degC-1,degB,degD_t+1,degE+2),               (H5)

with deg0=-infinity. The full regular cubic height box is used here,
not an odd-only subclass.

## 2. An exact scalar integration-by-parts recurrence

Translate t=v-(1+lambda)/3, so Phi_lambda=t^3+p t+r, where

    p=-(lambda^2-lambda+1)/3,
    r=-(lambda+1)(2lambda-1)(lambda-2)/27,
    Delta=-4p^3-27r^2=lambda^2(1-lambda)^2.

The translated endpoints are distinct roots of Phi_lambda, even when
the third root is repeated. For an integer m>=1 set

    S_m(t)=(9m+6)t^4+(15m+12)p t^2-3r t+4(m+1)p^2.

A direct universal polynomial identity is

    S_m' Phi+m S_m Phi'
      =3(3m+2)(3m+4)Phi^2
         -27(m+1)(2m+1)r Phi-m(m+1)Delta.              (H6)

Multiplication by Phi^(m-1) makes the left side the derivative of
S_m Phi^m. Its endpoint values vanish for m>=1. Integration gives

    3(3m+2)(3m+4)kappa_(m+1)
      =27(m+1)(2m+1)r kappa_m
         +m(m+1)Delta kappa_(m-1).                    (H7)

The m=0 relation is NOT asserted: the integration-by-parts boundary
term need not vanish there. We use kappa_0=1 and
kappa_1=(2lambda-1)/12 directly.

Write

    a_m=9(m+1)(2m+1)/[(3m+2)(3m+4)],
    b_m=m(m+1)/[3(3m+2)(3m+4)].

Then kappa_(m+1)=a_m r kappa_m+b_m Delta kappa_(m-1)
for m>=1.

## 3. Fixed-gap finiteness and arbitrarily late nonzero blocks

First assume r Delta!=0. In a constant-field extension choose omega
with omega^2=-Delta and set zeta=r/omega!=0. Choose c_0=1 and
c_(m+1)=a_m c_m. The normalized VALUES
R_m=kappa_m/(c_m omega^m) satisfy

    R_(m+1)=zeta R_m-t_m R_(m-1),   m>=1,
    t_m=b_m/(a_m a_(m-1))
       =(m^2-1/9)/[27(4m^2-1)]
       =(1/108)[1+5/(36(m^2-1/4))].                   (H8)

No claim R_1=zeta R_0 is made. The actual moments need not start the
monic-polynomial recurrence at index0.

There are no adjacent zero moments, since(H7) would propagate two
consecutive zeros backward to kappa_0=0. If kappa_m=kappa_(m+r_0)=0
with m>=1 and r_0>=2, then R_(m+1)!=0. Define local continuants by

    D_(m,0)=1, D_(m,1)=zeta,
    D_(m,j)=zeta D_(m,j-1)-t_(m+j)D_(m,j-2).

The value recurrence gives
R_(m+1+j)=R_(m+1)D_(m,j), hence

    D_(m,r_0-1)(zeta)=0.                              (H9)

The h-by-h path determinant D_(m,h)(zeta) has squared off-diagonal
entries t_(m+2),...,t_(m+h). It is a rational function of m with
coefficients in Q[zeta]. For fixed h let T=1/108 and define

    D_0(z)=1, D_1(z)=z, D_h(z)=zD_(h-1)(z)-T D_(h-2)(z).

These constant continuants are squarefree in every characteristic-zero
field. An entirely algebraic proof uses the identity

    (z^2-4T)D_h'(z)
       =h z D_h(z)-2T(h+1)D_(h-1)(z),                 (H10)

verified by induction from the recurrence. A common zero of D_h and
D_h' would force D_(h-1)=0, and then successive backward recurrences
would force D_0=0, impossible.

For each fixed integer j the formal Laurent expansion in1/m is

    t_(m+j)=T[1+5/(36m^2)]+O(m^-3).

All finitely many path entries have the same first correction.
Writing sigma=1+5/(72m^2), simultaneous scaling of the path
off-diagonal entries gives sigma^h D_h(z/sigma). Since determinants
are polynomial in their squared entries, this proves formally

    D_(m,h)(z)=D_h(z)
       +(5/(72m^2))[hD_h(z)-zD_h'(z)]+O(m^-3).         (H11)

This is a rational-function Laurent expansion, not an analytic
asymptotic theorem or a limit inferred from samples.
At a fixed zeta!=0, either D_h(zeta)!=0 or the m^-2 coefficient is
-(5/72)zeta D_h'(zeta)!=0 by squarefreeness. Thus D_(m,h)(zeta)
is a nonzero rational function of m. It has only finitely many zeros
at positive integer m; its displayed denominators are nonzero there.

It follows from(H9) that, for each fixed gap r_0>=1, there are only
finitely many pairs of zero indices m,m+r_0. A finite union for
1<=r_0<=L shows that successive zero gaps tend to infinity if the
zero set is infinite. If it is finite, the late tail has no zeros.
In either case, arbitrarily late blocks of any prescribed length
consist entirely of NONZERO kappa values.

The exceptional scalar cases are exact:
- Delta=0 means lambda=0 or1. Here r and kappa_1 are nonzero;
  (H7) propagates nonvanishing at every order.
- r=0 means lambda=-1,1/2 or2. For lambda=-1 or2,
  kappa_1=-1/4 or1/4 and both parity sequences stay nonzero.
- At lambda=1/2, precisely the odd moments vanish.

Thus EVERY lambda!=1/2 has arbitrarily late nonzero blocks of every
finite length, over every characteristic-zero coefficient field.
The construction in the constant extension containing omega does
not change which original kappa values are zero.

## 4. Algebraic series cannot have arbitrarily long late zero blocks unless polynomial

Let F(u)=sum f_n u^n be algebraic over K(u).
The extension K(u)(F) is finite and separable. Differentiating its
minimal polynomial expresses F' in that field; the derivation extends
there, so all derivatives lie in a finite-dimensional K(u)-space.
A nonzero linear relation among finitely many derivatives, with
denominators cleared, gives

    sum_j p_j(u)F^(j)(u)=0, p_j in K[u].

If p_j=sum_d c_(j,d)u^d, coefficient extraction gives, for all
sufficiently large n,

    sum_(j,d) c_(j,d)(n-d+j)_(falling j) f_(n-d+j)=0.

Grouping by shift k=j-d produces a finite polynomial recurrence
sum_k Q_k(n) f_(n+k)=0. For a fixed shift, the nonzero terms have
distinct degrees j in n, so the largest cannot cancel identically.
In particular a nonzero leading coefficient Q_(k_max) exists.

Past its finitely many integer zeros, any block of k_max-k_min
consecutive zero coefficients propagates the entire tail zero.
If k_max=k_min, the tail is already zero. Therefore a nonpolynomial
algebraic series has a bounded length of late zero runs.

Combining this with Section3 gives the useful support conclusion:

    If F is algebraic and kappa_m f_m=0 eventually,
    then F is polynomial.                            (H12)

This proof works directly over any characteristic-zero field.
It needs no analytic gap theorem or whole-field embedding into C.

## 5. Purity as an algebraic diagonal series

Set a=E h^3 and ell=h/2. Their initial values are nonzero and ell0=1.
The substitution v=(t-s)/h takes the interval endpoints to0,1 modulo x.
The derivative of an antiderivative of Phi_lambda(v)^m vanishes to
order at least m at both endpoints. Moving either endpoint by O(x)
therefore changes the integral by O(x^(m+1)). Hence

    (1/2)integral_-1^1 P(x,t)^m dt
       =kappa_m a^m ell modulo x^(m+1).               (H13)

This includes lambda=0 or1, where an endpoint root is repeated.
Thus purity is exactly kappa_m [x^m]a^m ell=0 for all m>=1.

The depressed coefficients of P are rational functions:

    B/E-D_t^2/(3E^2)=p h^2,
    C/E-BD_t/(3E^2)+2D_t^3/(27E^3)=r h^3.             (H14)

The constants p,r cannot both vanish.

If r=0, the allowed cases are lambda=-1 or2. Section3 says that
every kappa_m is nonzero. For the GIVEN unit a, the coefficient
conditions [x^m]a^m ell=0 with ell0=1 determine each new coefficient
of ell uniquely; its coefficient in equation m is a0^m.
The series 1-xa'/a satisfies them, proving(H3) in this case.
No independence between the actual given a and ell is presumed.

Now assume r!=0. Equation(H14) makes h^3 and a rational, with a a
unit at zero. Put

    u=x/a(x), J_a=1-xa'/a=xu'/u,

and let x=chi(u) be the formal inverse. The inverse is algebraic,
as is h. Therefore the formal series

    Bdiag(u)=(ell/J_a)(chi(u))

is algebraic, with constant coefficient1. A formal residue
substitution using du/u=J_a dx/x gives exactly

    [x^m]a^m ell=[u^m]Bdiag(u).

Purity and(H12) make Bdiag polynomial. Also its linear coefficient
vanishes because kappa_1=(2lambda-1)/12!=0.

## 6. The omitted-value argument forces Bdiag=1

If a polynomial Bdiag with Bdiag0=1 and linear coefficient0 were
nonconstant, it would have at least TWO distinct nonzero roots over
an algebraic closure. Indeed a polynomial with just one distinct
root b!=0 and constant coefficient1 is (1-u/b)^d, whose linear
coefficient -d/b is nonzero in characteristic zero.

From the definition of Bdiag,

    h=2 J_a(x) Bdiag(u(x)).

This makes h rational. Let b!=0 be a root of Bdiag.
If a finite alpha satisfied u(alpha)=b, then alpha!=0 because u(0)=0.
The function a=x/u is regular and nonzero at alpha. The functions u'
and J_a=xu'/u are regular there and may vanish; they have no pole.
Thus h vanishes there,
making E=a/h^3 have a pole, contrary to E being polynomial.

Consequently u omits every root of Bdiag on the finite affine line.
A nonconstant rational map can omit at most one finite value:
write u=p_0/q_0 in coprime form over the algebraic closure.
If b is omitted, p_0-bq_0 is a nonzero constant. Two omitted values
force q_0 constant by subtraction and then p_0 constant too.
This contradicts u'(0)=1/a0!=0.

Thus Bdiag=1, proving ell=J_a. This completes the necessity in(H3).
Conversely, if ell=1-xa'/a, then

    [x^m]a^m ell=[x^m](a^m-x(a^m)'/m)=0,

which proves purity at every order regardless of individual kappa
zeros. The exclusion lambda=1/2 is essential to the argument:
it supplies the zero linear coefficient of Bdiag.

## 7. Rationality, squarefree denominator and degree at infinity

Since one of p h^2 or r h^3 is rational, h'/h is rational.
Equation(H3) writes

    h/2=1-x(E'/E+3h'/h),

so h is rational in all cases, including lambda=-1 or2.
The rational depressed center
s+(1+lambda)h/3=-D_t/(3E) then makes s rational.

Over an algebraic closure, a finite nonzero zero of h would give
a=E h^3 positive order. The right side of(H3) would have a nonzero
simple pole there, while h/2 would vanish. At zero h=2.
Thus h has no finite zero, so h=2/Q with Q polynomial and Q0=1.
A pole of h can only be simple because a logarithmic derivative
has at most simple poles; Q is squarefree.

If Q=1, (H3) gives E'/E=0, so E is constant.
If q=degQ>=1, write e=degE. At infinity, h/2 tends to zero and
x a'/a tends to e-3q. Thus

    e=3q+1,       deg_at_infinity a=1.                 (H15)

All assertions are rational-function valuations, not analytic limits
depending on an embedding of K into C.

## 8. Marked tails: the second algebraic-support argument

For every psi in K[x], equation(H13) gives all needed constant-term
coefficients: its error times psi(x)x^(-m) has positive x order.
Under(H3) formal residue substitution yields

    Lambda(psi(x) f^m)=kappa_m [u^m]psi(chi(u)),
    u=x/a(x).                                         (H16)

If q>=1 and this marked sequence vanished eventually for a
nonconstant psi, the algebraic series psi(chi(u)) would be polynomial
by(H12). Write it as H(u). Substitution then gives the rational identity

    psi(x)=H(u(x)).

By(H15), u=x/a is regular at infinity, so the right side is regular
there while the nonconstant polynomial psi has a pole.
This contradiction proves infinitely many nonzero marked moments.
In particular x itself proves unsafety.

It is important to reuse(H12) here: nonpolynomiality alone would
not suffice because some individual kappa_m may vanish.
Constants give zero at every positive order by purity.

## 9. The constant-Q branch is safe

If Q=1, then h=2,E is constant and s is polynomial. Write
s=-1+xH, H in K[x]. The phase is the polynomial Hopf shear
t->t-xH, x fixed, of

    E(t+1)(t-1)(t+1-2lambda)/x
       =-2E y(t+1-2lambda).                           (H17)

The base has strictly negative angular weight. Every fixed marker
has bounded Laurent x support, so its product with sufficiently
large powers has no constant x coefficient.

The needed shear invariance is direct. Put J=x partial_t-t partial_y.
For any H in K[x],g in A, integration on the Laurent chart gives

    Lambda(HJg)=CT_x xH(x)
          [g(x,0,1)-g(x,0,-1)]/2=0.

The bracket is polynomial and xH has positive x order.
The derivation HJ is locally nilpotent, and its finite exponential
preserves Lambda; its inverse transports every fixed marker.
Applying exp(-HJ) to(H17) proves safety.
No all-automorphism theorem or automorphism classification is needed.

## 10. Exact degree spectra and realizations

If q>=1, e=3q+1. If s is bounded at infinity, C,B,D_t have degree
at most e while E has degree e. Formula(H5) gives d_f=e+2=3q+3.

If s has a pole at infinity of integer order j>=1, h decays there.
The leading powers of s dominate all lower-shape h terms, giving

    degC=e+3j, degB=e+2j, degD_t=e+j.

Thus

    d_f=max(e+3j-1,e+2j,e+j+1,e+2)=3(q+j).              (H18)

Every unsafe degree is therefore divisible by3 and at least6.

If Q=1,E is constant and s is polynomial. Constant s gives degree2.
For degs=j>=1 the same leading-term comparison gives

    d_f=max(3j-1,2j,j+1,2)=3j-1.                       (H19)

All safe degrees are exactly2 modulo3, starting at2.

Both lists are realized for EACH fixed lambda!=1/2:
- For each j>=1, set Q=1,E=1,h=2,s=-1+x^j.
  This gives a regular pure safe phase of least degree3j-1.
- For each j>=1, set Q=1+x,E=Q^4,h=2/Q,s=-1+x^j.
  Here a=8Q and ell=1/Q=1-xa'/a, so the phase is pure,
  of least degree3j+3, and unsafe. In this particular family

      Lambda(x f_j^m)=8^m kappa_m(lambda).

  This may vanish at some orders, but Section3 guarantees infinitely
  many nonzero values. It does NOT become an EVERY-order formula
  of nonvanishing for arbitrary noncentral lambda.

These are phase statements in the specified affine-shape branch.
An original-ring primitive lift or absolute nonpower claim would
require its own literal geometric verification and is not asserted
for every parameter family in this theorem.

## 11. Scope and inspected classical context

The main hypothesis is ONLY lambda!=1/2 together with the regular
constant-shape presentation(H2), E0!=0 and the endpoint normalization.
There is no all-kappa-nonvanishing assumption.
Lambda=-1 and2 are INCLUDED as off-center selections of a symmetric
underlying cubic. Lambda=0 and1, with repeated cubic roots, are also
included. The midpoint selection lambda=1/2 remains outside this
classification; no claim is made that it exhausts all earlier odd
classes or that purity forces a constant shape.

The algebraic=>D-finite and coefficient-recurrence facts are classical.
The actual author-hosted scan of R.P.Stanley, *Differentiably finite
power series*, European Journal of Combinatorics1(1980),175--188,
was visually inspected at Theorem2.1 printedp178 and
Definition1.3/Theorem1.5 printedp176:
https://math.mit.edu/~rstan/pubs/pubfiles/45.pdf
(DOI10.1016/S0195-6698(80)80051-5).
Its displayed statements use complex coefficients; Section4 gives
the needed characteristic-zero proof directly. The citation is
credit, not an uninspected all-field theorem import.

The local continuants also have a classical associated-Jacobi
identification. Jet Wimp, *Explicit formulas for the associated Jacobi
polynomials and some applications*, Canadian Journal of Mathematics
39(4)(1987),983--1000, equation(12), printedp986, uses
gamma=alpha+beta+1 and association shift c. At alpha=1/3,beta=-1/3,
its monic squared recurrence coefficient is
(N^2-1/9)/(4N^2-1), N=n+c. Thus, if P_hat denotes the MONIC associated
Jacobi polynomial in that convention,

    D_(m,h)(z)=27^(-h/2) P_hat_h^(1/3,-1/3)
                               (sqrt(27)*z;c=m+1).

The first squared entry is t_(m+2), as in Section3.
See DOI10.4153/CJM-1987-050-4:
https://doi.org/10.4153/CJM-1987-050-4.
This identifies the LOCAL continuants only; it does not assert that
the actual scalar kappa values have ordinary Jacobi initial data.
The source is context, not an imported fixed-association-shift
zero theorem, spectral result or all-field claim. Section3 proves
the needed rational-function and zero-gap statements directly.

The scalar integration-by-parts identity, local determinant argument,
omitted-value proof, endpoint identities, marked-tail argument and
degree comparison are all written here. No FPYZ interval rigidity,
spectral zero-location theorem, analytic asymptotic theorem,
all-Aut classification or finite-prefix proof is a premise.
No historical-priority, external-review or formal-verification claim
is made.
