# Complete classification with one affine circular variable

14 November 2026. This follow-up uses the retained [formal Gaussian--Lagrange
identity GL1](../../foundations/extended-geometry/FORMAL_GAUSSIAN_LAGRANGE_LEMMA.md).
The [sparse Pareto theorem](../../foundations/extended-geometry/THREE_REAL_GMC_MINIMAL_FAILURES.md)
already identifies the Gaussian triangular shear as a safe survivor.
The result here classifies an entire polynomial architecture, including
members with more terms. It does not claim novelty or a classification of
arbitrary ternary polynomials. All assertions below have written all-order
proofs; finite regressions check transcription only.

Let Z,W,T have Wick pairings E(ZW)=E(T^2)=1 and all other pairings zero.
Consider arbitrary complex polynomials h,v,u,r in one variable and

    P=W*h(Z)+v(Z)*T^2+u(Z)*T+r(Z).                    (1)

## Complete pure-moment classification

Every positive pure moment of P is zero if and only if P belongs to exactly
one of these cases, distinguished by h(0) and deg h:

1. **The zero branch:** h is divisible by Z^2, and v,u,r are divisible by Z.
   Every monomial has strictly positive Z-minus-W weight. For every Q,
   E(QP^m)=0 for m>deg Q.
2. **The constant branch:** h=a is a nonzero constant, v=0, u is arbitrary,
   and r=-Z*u^2/(2a). This is a polynomial Gaussian triangular shear of aW,
   so it satisfies eventual mixed vanishing.
3. **The linear branch:** h=a+bZ, with a*b nonzero,

       v=(h-h^3/a^2)/(2Z),
       u=h^2*k,       r=-(a^2/2)*Z*h*k^2,            (2)

   for an arbitrary polynomial k. Here every positive pure moment vanishes,
   while the fixed multiplier Z satisfies

       E(Z*P^m)=m!*a*b^(m-1) != 0  (m>=1).           (3)

   Thus every member is a GMC failure. Each has degree at least four,
   because v has a nonzero Z^2 coefficient multiplying T^2.

The zero polynomial is included in case 1. Exchanging Z and W gives the
reflected architecture. A complex orthogonal change preserves the Wick
functional, but this theorem does not assert that every ternary cubic can
be moved into (1).

## The formal generating function

There is a unique g(t) in t*C[[t]] with g=t*h(g). Put

    J=1-t*h'(g),       K=1-2t*v(g).

Applying GL1 at fixed T and then the one-variable Gaussian quadratic
identity gives, coefficientwise,

    M(t):=E(exp(tP))
      = exp(t*r(g)+t^2*u(g)^2/(2K)) / (J*sqrt(K)).   (4)

The square root has constant term one. The substitutions are t-adic;
no analytic convergence or positivity of complex coefficients is assumed.
GL1 allows the t-dependent numerator in this application: its proof is
coefficientwise and every fixed t-degree uses finitely many terms.

We use an elementary rigidity lemma. If S and R are rational functions
regular at zero, S(0)=0, R(0)=1 and exp(S)=R as formal series, then S=0
and R=1. Indeed S'=R'/R. A rational derivative has zero residue at every
finite pole, whereas the latter logarithmic derivative has the nonzero
integer order of each zero or pole of R as its residue. Thus R has no
finite zero or pole, so is constant; the identity then forces S'=0.

## Proof on h(0)!=0

Write z=g(t), so t=z/h(z) is an invertible formal change. Define

    A=h(z), B=h(z)-2z*v(z), C=h(z)-z*h'(z),
    S=z*r(z)/A+z^2*u(z)^2/(2*A*B).

Then (4) is exp(S)/(C/A*sqrt(B/A)). Hence M=1 and the rigidity lemma
applied to exp(2S)=C^2*B/A^3 are equivalent to

    C^2*B=A^3,               2*B*r+z*u^2=0.         (5)

Both equations are polynomial identities, with A(0)=B(0)=C(0)!=0.

The first identity forces deg h<=1, in every polynomial degree. To prove
this, let alpha_1,...,alpha_q be the distinct roots of h, with multiplicities
m_i. None is zero. After canceling the common factors of C and h,

    C/h = 1-z*sum_i m_i/(z-alpha_i).

If n=deg h>=2, the reduced numerator has degree q, with leading coefficient
1-n. At z=alpha_i it is -alpha_i*m_i*product_(j!=i)(alpha_i-alpha_j),
which is nonzero. Therefore it has a root away from all roots of h. This
contradicts C^2 dividing h^3 in (5). So n<=1.

If h=a, (5) gives v=0 and r=-z*u^2/(2a), exactly case 2.
If h=a+bz with b!=0, C=a and B=h^3/a^2. Thus v has the form (2), while
r=-a^2*z*u^2/(2h^3). Since h is coprime to z, polynomiality forces
ord_h(u)>=2. Writing u=h^2*k gives precisely the stated r. Conversely
these identities make (4) identically one. With a multiplier Z, GL1
replaces the numerator by g(t), giving

    E(Z*exp(tP))=g(t)=a*t/(1-b*t).

Coefficient extraction proves (3) for every positive m, including all
polynomial choices of k. This extends the retained rank-one quartic
classification, rather than reusing a finite-prefix survivor as evidence.

## Proof on h(0)=0

The unique formal fixed point is g=0. Let b=h'(0), c=v(0), d=u(0),
e=r(0). Equation (4) becomes

    exp(t*e+t^2*d^2/(2*(1-2ct))) /
           ((1-bt)*sqrt(1-2ct)).

Its being one implies, by rational exponential rigidity,

    (1-bt)^2*(1-2ct)=1,
    t*e+t^2*d^2/(2*(1-2ct))=0.

A product of polynomials equals one only when both factors are constant,
so b=c=0 and then e=d=0. These are exactly the divisibilities in case 1.
Their converse and all mixed conclusions follow from positive isotropic
weight. This completes the full pure-moment classification.

## Sharp mixed cutoff in the constant branch

Let h=a!=0 and u have degree j>=1 with leading coefficient c. The full
source identity, with scalar source variables alpha,beta,gamma, is

    E(exp(tP+alpha*Z+beta*W+gamma*T))
      = exp(alpha*(at+beta)+gamma^2/2
             +gamma*t*u(at+beta)
             -beta*t*u(at+beta)^2/(2a)).             (6)

Integrate W first by the circular translation identity, then integrate T;
the pure exponent cancels as displayed. Assign source weights 1, 2j+1,
j+1 to alpha,beta,gamma. Each term in the exponent has t-degree at most
its source weight. Thus

    E(Z^i*W^l*T^k*P^m)=0
        for m>i+(2j+1)l+(j+1)k,

and the uniform bound is m>(2j+1)*deg Q. The leading beta*t^(2j+1)
term in (6) also gives, for every d>=1,

    E(W^d*P^((2j+1)d))
        =((2j+1)d)!*(-c^2*a^(2j-1)/2)^d != 0.       (7)

So the uniform coefficient 2j+1 is sharp for this family. If u is constant,
P is isotropic linear and the bound is m>deg Q; u=0 is included.
A second proof of safety uses the Gaussian-preserving polynomial shear
exp(-(u(Z)/a)*(Z*d/dT-T*d/dW)) applied to aW. Its coefficient is invariant
under the locally nilpotent derivation, so the exponential is polynomial;
it fixes 2ZW+T^2 and preserves the Gaussian functional by integration by
parts. This agrees with the prior sparse-shear argument.

## Complete cubic subarchitecture and a finite gate

Impose exactly the total-degree-three bounds

    deg h<=2, deg v<=1, deg u<=2, deg r<=3.          (8)

This is a twelve-coefficient family in the twenty-dimensional space of
ternary cubics (including lower terms). Case 3 is impossible because its
T^2 coefficient has degree two. In case 2 polynomiality and deg r<=3
force deg u<=1. Thus the complete moment-null cubic locus in this family
consists of the positive-weight branch and

    P=aW+(c+dZ)T-Z*(c+dZ)^2/(2a),   a!=0.          (9)

Every member is GMC-safe, with the uniform sharp mixed cutoff m>3 deg Q.
When d!=0, a complex orthogonal shear takes (9) to the already retained
three-term Gaussian shear; the closure here is the exhaustive coefficient
classification of (8). Formula (7) gives its sharp equality cases.

More generally, in total degree at most D>=2, the first 5D-3 pure moments
suffice throughout architecture (1). This specializes to twelve for cubics
and seventeen for quartics. Neither number is claimed optimal. On A(0)!=0 put

    F=2z*r*B+z^2*u^2,      S=F/(2AB).

The derivative of log M(z/A) has denominator 2*A^2*B^2*C and numerator

    C*(F'*A*B-F*(A'*B+A*B'))
       -2*A^2*B^2*C' +3*A*B^2*C*A' -A^2*B*C*B'.  (10)

In total degree D, deg A,B,C<=D-1 and deg F<=2D. Hence (10) has degree
at most 5D-4. The first 5D-3 zero moments make M=1+O(t^(5D-2)); the
invertible change t=z/A forces its logarithmic-derivative numerator to
vanish to order at least 5D-3, above that degree. Thus it is identically
zero and M=1. On h(0)=0 the explicit formula above has a logarithmic
derivative with denominator (1-bt)*(1-2ct)^2 and numerator of degree at
most three. Four zero moments already suffice on that chart, in every
degree. The bound 5D-3 is at least four for D>=2. The proof works for
arbitrary complex coefficients and retains coefficient specializations.

In degree at most four, the same classification closes the whole
sixteen-coefficient family deg h<=3, deg v<=2, deg u<=3, deg r<=4.
Its failure locus is exactly (2) with k=k0+k1*Z and a*b!=0. Polynomial
coefficient degrees force deg k<=1 from both u=h^2*k and
r=-(a^2/2)*Z*h*k^2, and these bounds are sufficient. Every safe constant
branch still has deg u<=1, so the sharp uniform mixed cutoff across the
safe quartics is m>3 deg Q. The one-sided branch has m>deg Q. The
seventeen-moment gate certifies pure vanishing on this entire quartic
architecture; the three-way classification then decides GMC safety.

The height-two architecture here omits W^2,W^3,W^2Z,W^2T,WT,WZT,WT^2,T^3.
The [subsequent formal Airy--Riccati theorem](CUBIC_HEIGHT_AIRY_OBSTRUCTION.md)
now closes the additional T^3 coefficient in the thirteen-coefficient cubic
family, with a sufficient 35-moment gate. The remaining seven monomials all
involve higher W powers or W/T coupling. No frame-removal theorem is supplied,
so the global OP-GMC3-CUBIC remains open. Use the smaller twelve-moment gate
in this note when the T^3 coefficient is already zero.
