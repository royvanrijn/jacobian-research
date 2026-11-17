# Rational Gaussian pullbacks and unbounded actual irregular rank

15 November 2026. Internal formal/algebraic construction; independent review
outstanding. All-order identities and the rank argument are written below.
The [checker](verify_rational_gaussian_pullbacks.py) replays small exact
transcriptions, not a universal proof or a general period-kernel solver.

This gives genuine ternary failures with nonconstant repeated-root b,
in contrast to the empty squarefree locus and the safe constant-b normal
forms. Their exact actual generic ranks are every positive even integer.
Neither a contour estimate nor convergence of a Wick integral is assumed.

## 1. Formal Gaussian transport, with its Jacobian retained

Write E_X(X^(2k))=(2k-1)!! and E_X(X^(2k+1))=0. If

    Y=X+O(t) in C[X][[t]],

then for every coefficientwise polynomial f,

    E_X[f(Y)*Y_X*exp((X^2-Y^2)/2)] = E_X[f(X)].               (1)

Everything is coefficientwise finite. To prove this set Y_s=X+s(Y-X).
After removing the common exponential, differentiation with respect to s
shows

    d_s(f(Y_s)*Y_(s,X)*exp((X^2-Y_s^2)/2))
      =(d_X-X)(f(Y_s)*d_s(Y_s)*exp((X^2-Y_s^2)/2)).

Gaussian integration by parts kills the right side. Set s=0 and s=1.
This proof uses formal polynomial coefficients, not a real change of
contour or an integrability assumption.

The formal-change conventions in
[Alexandrov et al.,2309.12176v3,section6.1](https://arxiv.org/html/2309.12176v3)
were read as a scope check. The latest listed version is11May2025.
Its KP/topological-recursion theorems are not used. Equation(1) is proved
directly here, so no omitted contour choice from that paper is imported.

## 2. Which dz/b have rational primitives?

For a nonzero polynomial b in C[z], a rational primitive R with R'=1/b
exists iff b is constant, or

    b(z)=beta*(1-z/r)^(nu+1),    nu>=1, beta!=0, r!=0,        (2)

or the same single-root-power form with its root at zero.
The nonzero-root normalization in (2) is chosen for the construction below.

Here is a complete elementary proof. A simple root of b would give a
nonzero simple-pole residue, impossible for a rational derivative. Suppose
the distinct roots have multiplicities m_i>=2, their number is s, and
D=sum m_i. A primitive can have only these finite poles, of orders m_i-1.
Its reduced denominator therefore has degree D-s. At infinity it tends
to a constant R_infinity, and

    R-R_infinity has order D-1 at infinity,

because R'=1/b has order D. A nonzero rational function with denominator
degree D-s cannot have order greater than D-s there. Thus s<=1.
Conversely a single power of multiplicity at least2 has the stated
rational primitive; a constant b has an affine primitive.
This classifies this rational-primitive mechanism, not all repeated-root
Gaussian phases.

Passing to a finite algebraic cover does not enlarge the exact-primitive
class. If an algebraic R has derivative1/b, its field trace divided by
the extension degree is a rational primitive, since differentiation
commutes with trace. Over C the difference from that rational primitive
has derivative zero and is constant, so R itself is rational. Thus neither
an algebraic inverse nor a higher-genus finite cover bypasses this obstruction.

For (2) set

    alpha=1/nu,   kappa=beta*nu/r,
    R(z)=((1-z/r)^(-nu)-1)/kappa.                            (3)

Then R(0)=0, R'=1/b, and the inverse branch at zero is

    R^(-1)(s)=r*(1-(1+kappa*s)^(-alpha)).                    (4)

A root at zero does not admit this formal branch and is not silently
localized. Its retained zero-component classification is separate.

## 3. A polynomial family from an exactly Gaussian rational phase

Choose any polynomial a satisfying

    (z-r)^ceil(2nu/3) divides a(z),                          (5)

and put delta=a(0)/beta and h(z,T)=a(z)+b(z)*T. Define

    G(z,T)=h/(2z)*(T^2-(R(z)*h/z-delta)^2),
    P=W*h(Z,T)+G(Z,T).                                     (6)

Despite its rational expression, G is polynomial, and condition(5) is
exactly its polynomiality condition.

At z=0, R/z=1/beta+O(z), so R*h/z-delta=T+O(z). The bracket in (6) is
divisible by z. At z=r, R has a pole of order nu. Write v=ord_r(a).
The constant-T coefficient has its potentially leading pole from
R^2*a^3, of order2nu-3v. If v<ceil(2nu/3), it dominates the possible
delta*R*a^2 term because v<nu, so it cannot cancel.
If3v>=2nu, every coefficient is regular: the other possible pole order
2nu-(2v+nu+1) is also nonpositive, as are nu-2v and all remaining ones.
There are no other finite poles. This proves necessity and sufficiency.

Let g=t*h(g,T) be the formal fixed point, g(0,T)=0, and

    D=1-t*(a'(g)+b'(g)*T),      Y=R(g)/t-delta.

Because b(0)=beta!=0, g=t*(a(0)+beta*T)+O(t^2) and Y=T+O(t), with
polynomial T-coefficients. Exact differentiation and (6) give

    Y_T=1/D,       t*G(g,T)=(T^2-Y^2)/2.                    (7)

The retained Gaussian-Lagrange identity therefore becomes

    F_Q(t)=E[Q(Z)*exp(tP)]
          =E_T[Q(g)*exp(tG(g,T))/D]
          =E_X[Q(R^(-1)(t*(X+delta)))].                    (8)

Here Q can be a polynomial or a rational function regular at0. Its
substitution has polynomial X-coefficients at each t order. Equation(1)
proves the last equality coefficientwise.

In particular F_1=1: **every P in (6) is pure-null in all orders**.
There is no analytic comparison or unreviewed action-filtration step in
this identity.

For constant b=beta, take R=z/beta and delta=a(0)/beta in the same formula.
It gives exactly the already retained Gaussian-shear normal form with
u=(a-a(0))/(beta*z). Thus the construction recovers that safe boundary,
but its nonconstant power-root members below all fail.

## 4. A fixed linear multiplier detects every member's failure

For alpha>0 rational define the formal series

    Y_(alpha,delta)(t)=E_X[(1+kappa*t*(X+delta))^(-alpha)].

Equations(4),(8) imply

    F_Z=r*(1-Y_(1/nu,delta)),
    E[Z*P^n]=-r*(-1)^n*(1/nu)_n*kappa^n*E[(X+delta)^n],
                     n>=1.                                (9)

The shifted Gaussian moments are nonzero infinitely often: their
exponential generating series exp(delta*s+s^2/2) is not a polynomial.
Every Pochhammer factor is nonzero. Hence the one fixed multiplier Z
has infinitely many nonzero moments, proving GMC failure for every
member of (6), not just for generic a.

On the particularly simple subfamily a=0, delta=0,

    E[Z*P^(2k)]=-r*(1/nu)_(2k)*kappa^(2k)*(2k-1)!! !=0,
    E[Z*P^(2k+1)]=0.                                       (10)

The parameters beta,r are any nonzero complex numbers. This is a fixed
three-variable architecture in every nu, with no squarefree-b hypothesis.

For a=0 an explicitly polynomial presentation avoids rational expressions.
Let s=1-Z/r and H_nu(s)=1+s+...+s^(nu-1). Then

    P=beta*s^(nu+1)*W*T
      +beta/(2*r*nu^2)*s^(nu+1)
          *[(nu^2-s^2*H_nu(s)^2)/(1-s)]*T^3.                (11)

The quotient is polynomial since its numerator vanishes at s=1.
Its degree is2nu-1 and its leading coefficient is1, so deg P=3nu+3.
For nu=1,beta=r=1 this is the sextic

    P=(1-Z)^2*(W*T+(2-Z)*T^3/2),
    E[Z*P^2]=-2,        E[Z*P^4]=-72.

These examples do not change the existing minimal ternary failure degree4.

## 5. The full actual generic image has exact rank 2nu

This section bounds the **full image**, not merely the pure seed or the
single marked series.

### 5.1 A 2nu-dimensional source for every multiplier

The Gaussian phase before integration is

    Psi=W*(t*h-Z)+t*G(Z,T)-T^2/2.

Localizing at b(Z) is valid for the chosen formal Wick branch: b(0)!=0.
For an amplitude rational in Z and regular at0, each t coefficient has
bounded W-degree; circular contraction uses only finitely many Taylor
coefficients at Z=0. Thus Wick expectation and integration by parts extend
to this localization. No real integral of a rational amplitude is asserted.

Set S=t*h-Z and T=(S+Z-t*a)/(t*b). The coordinate Jacobian is1/(t*b).
Writing the remaining phase as H(Z,S), subtract H(Z,0) and shift W by the
polynomial (H(Z,S)-H(Z,0))/S. The phase becomes

    S*V + Phi(t,Z),
    Phi=-1/2*(R(Z)/t-delta)^2.

The normal pair S,V has one-dimensional exponential de Rham cohomology.
Its polynomial amplitudes reduce by the signed circular contractions
E(S^i V^j)=delta_ij*(-1)^i*i!, or directly by its two IBP relations.
Every original polynomial multiplier therefore gives a Laurent-rational
amplitude on the remaining curve, with the full volume dz/(t*b).
The W-dependent directions have not been omitted.

Use u=(1-Z/r)^(-1). The curve ring is C(t)[u,u^(-1)], and

    Phi=-(u^nu-1-kappa*delta*t)^2/(2*kappa^2*t^2),
    dz/(t*b)=nu/(kappa*t)*u^(nu-1) du.                       (12)

The polynomial phase has degree2nu in u with nonzero leading coefficient.
Its twisted Laurent de Rham quotient has basis

    u^(-1), 1, u, ..., u^(2nu-2).                           (13)

To check this, reduce all powers <=-2 using the lowest derivative term of
(d_u+Phi_u)u^j, and all powers >=2nu-1 by its leading phase term.
For independence, a negative lowest power of a Laurent primitive would
produce a power <=-2 that cannot cancel; a nonnegative primitive would
produce a highest power >=2nu-1. Neither can give a nonzero combination
supported in (13). Thus every original period lies in a rank2nu image.

### 5.2 Actual marked periods attain the whole bound

For j=1,...,nu, the polynomial multiplier (Z-r)^j gives, up to its
nonzero constant factor, Y_(j/nu,delta). Each is annihilated by

    kappa^2*t^3*Y''
      +(2*(alpha+1)*kappa^2*t^2-kappa*delta*t-1)*Y'
      +(alpha*(alpha+1)*kappa^2*t-alpha*kappa*delta)*Y=0.     (14)

One derivation uses w=1/(kappa*t)+delta and
f(w)=E(X+w)^(-alpha). Gaussian IBP gives f''+w*f'+alpha*f=0;
Y=(kappa*t)^(-alpha)*f(w) yields (14).
Alternatively its explicit coefficients verify (14) directly.

Every selected series has minimal homogeneous order exactly2.
For delta=0 its coefficient ratio along even indices is

    c_(2k+2)/c_(2k)
       =kappa^2*(alpha+2k)*(alpha+2k+1)/(2*(k+1)),

whose modulus tends to infinity. Its radius of convergence is zero.
For general delta use the invertible formal rational substitution

    Y_(alpha,delta)(t)
      =(1+kappa*delta*t)^(-alpha)
          *Y_(alpha,0)(t/(1+kappa*delta*t)).                 (15)

Convergence would imply convergence of Y_(alpha,0), a contradiction.
A nonzero formal power series satisfying a first-order homogeneous
equation has rational logarithmic derivative regular at the origin,
after removing its finite order of vanishing; hence it is convergent.
Thus (14)'s order cannot drop to1 or0.

At infinity (14) is regular singular with exponents alpha and alpha+1.
For the distinct j=1,...,nu their exponent classes modulo integers are
distinct. Their rank-two cyclic images therefore have a direct sum:
after formal completion at infinity, generalized formal-monodromy
eigenspaces separate the j's. Any kernel of the sum map splits by these
eigenvalues, and its intersection with each individual image is zero.
Faithfulness of completion gives zero kernel over C(t).
This is the standard regular-singular formal classification, for example
[van der Put--Saito,0902.1702v1,section1.2](https://arxiv.org/html/0902.1702v1);
no hypergeometric irreducibility theorem is needed.

The nu cyclic images have total rank2nu and are contained in the full image.
Together with the upper bound this proves

    rank(actual full period image)=2nu.                    (16)

The polynomial seed1 is already in the alpha=1 cyclic image, by the exact
inhomogeneous identity

    kappa^2*t^3*Y_(1,delta)'
       +(kappa^2*t^2-kappa*delta*t-1)*Y_(1,delta)+1=0.

Thus no extra trivial rank is being added to (16).

## 6. Genuine irregular failure with no rational-pole core

Minimality in section5 makes the Fuchs test on (14) an actual-image test.
Its normalized Y' coefficient has a pole of order3 at0, exceeding the
allowed order1 for a regular-singular second-order equation. Thus this
actual cyclic image, and hence the full connection, is irregular at its
formal origin. By (14),(16), the full connection is regular singular at infinity.

All coefficients of the curve connection in basis(13), including the
normalization derivative -1/t, are Laurent polynomials in t. The phase
leading coefficient is a nonzero multiple of t^(-2), and the Laurent
reductions invert only that coefficient and nonzero integers.
It is therefore an ordinary connection at every finite t!=0.
The original amplitudes reduce into this connection there; a full-rank
D-submodule of an ordinary connection equals it locally, since such a
connection has no punctual quotient.

Any rational element of the actual image consequently has no pole at
t!=0; its formal power series has no pole at0 either. Therefore

    M intersect C(t)=C[t].                                  (17)

The entire punctual safety submodule from the
[rational-core theorem](RATIONAL_CORES_AND_BOUNDARY_SUBSPACES.md) is zero.
The obstruction is torsion-free of generic rank2nu-1 and includes genuine
irregular directions. This explicitly realizes the alternative to the
known quartic rank-one finite-pole failure.

These statements concern this constructed family. They do not classify
all Gaussian higher-rank connections, all possible ranks, or all lattices
over them.

## 7. Polynomial specializations with arbitrarily large rank drops

In (11) put r=e and beta=e^(3nu+1). All coefficients become polynomial in e.
The W coefficient is e^(2nu)*(e-Z)^(nu+1).
Since the polynomial quotient in (11) has leading coefficient1, at e=0,

    P_(nu,0)=(-1)^nu/(2nu^2)*Z^(3nu)*T^3.

This is safe by positive circular charge, with sharp cutoff
m>deg Q/(3nu). The equality witness m=2,Q=W^(6nu) is nonzero.
Its full actual image is C[t], of rank1. For every e!=0, (16) gives rank2nu
and the member fails. All positive pure moments vanish identically as
polynomial functions of e, by (8) on e!=0 and coefficientwise descent.

Thus even at fixed Gaussian dimension3, actual period rank can fall
from an arbitrarily large even value to1 at an exceptional coefficient
point. This is an explicit family of kernel/irregularity specialization,
not a counterexample to complex-field safe-locus constructibility or
uniform bounds. No ordinary cutoff-ideal chain for this new family is
claimed to have been computed.

## 8. Structural consequences and assurance

Rational primitives of dz/b are now completely classified. Their nonzero
single-root-power branch supports an explicit polynomial Gaussian-pullback
family with a complete polynomiality condition, exact all-order pure and
marked identities, exact full even-rank tower, and polynomial rank-drop
specializations. General repeated-root b with several distinct roots,
other pure-null members for the same b, b=0 beyond retained cases, several
spectators, global uniformity and valid operator compression remain open.

The primary-source check uses the2025 formal-Gaussian conventions and the
classical regular-singular classification only. No KP/Painleve theorem,
contour comparison, saddle estimate, general minimizer or kernel producer
is imported into the all-order identity proof. No novelty or external
assurance is claimed. The small replay covers identities and low-degree
quotient transcriptions; the all-nu proofs remain independently unreviewed.

[Source and assurance audit](rational_gaussian_pullback_source_audit.json).

