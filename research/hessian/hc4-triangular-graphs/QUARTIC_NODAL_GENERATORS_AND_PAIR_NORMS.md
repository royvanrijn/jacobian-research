# Quartic nodal generators and the unordered-pair norm

Let A_0={f in C[t]:f(1)=f(-1)}. If A_0=C[p,q] and p has degree
FOUR, then a polynomial target shear q->q-H(p) makes q quadratic.
Consequently

    integral_(-1)^1 p*dq != 0.                              (1)

The other generator has unrestricted degree. The proof treats separately
the irreducible six-pair Laurent cover and the split even-quartic cover.
It uses no classification of all plane embeddings of the nodal curve.

Together with the preceding degree-two/three result, nodal action is now
nonzero when either generator has degree at most FOUR. This strengthens
the actual degree-two-conductor exclusion to directional output degree
at most four. A separate index argument excludes every prime generator
degree at least five. The unrestricted action is proved separately in
[the all-degree nodal-action theorem](UNRESTRICTED_NODAL_ACTION_AND_DEGREE_TWO_CONDUCTOR_EXCLUSION.md),
using the audited classical classification.

## 1. The index and the quartic pair curve

The index argument from Section1 of
[the low-degree nodal-action proof](LOW_DEGREE_NODAL_ACTION_AND_CONDUCTOR_EXCLUSION.md)
gives the following exact setup. Subtract the common node value of p,
and scale p to make it monic. Then

    p(t)=(t^2-1)*(t^2+b*t+c), R_0=C[t], D=C[p].

Since dim_C(R_0/A_0)=1, the inclusion determinant of the D-basis
1,q,q^2,q^3 into 1,t,t^2,t^3 is

    Delta(U)=eta*U, eta in C*.                              (2)

Here U represents p. In particular q is primitive over C(p).

Depress the quartic by r=t+b/4, writing

    p=r^4+A*r^2+B*r+C.

The letters A,B,C in this paragraph are COMPLEX COEFFICIENTS, not rings.
For a pair of distinct generic roots r_i,r_j of p(r)=U put s=r_i+r_j.
Comparing the two quadratic factors of the quartic gives

    v=r_i*r_j=(s^2+A+B/s)/2,
    U=C-(s^2+A)^2/4+B^2/(4*s^2).                           (3)

When B!=0, the unordered-pair function field is C(s), of degree SIX
over C(U). Its monic equation is

    H_6(s)=s^6+2*A*s^4+(A^2-4*C+4*U)*s^2-B^2=0.             (4)

Irreducibility follows directly from the primitive equation linear in U:
the coefficient 4*s^2 is relatively prime to its constant part because
B!=0. Thus its six conjugates are exactly the six unordered pair sums.
The coordinate ring C[s,s^-1] is finite over C[U]; no classification of
quartic monodromy is being used.

The divided difference

    g(s)=(q(r_i)-q(r_j))/(r_i-r_j)

is a Laurent polynomial in s. Indeed a symmetric polynomial in r_i,r_j
is polynomial in their sum s and product v from(3). Vandermonde
determinants give the exact equality

    Norm_(C(s)/C(U))(g)=Delta(U).                            (5)

This norm is the basis INDEX, not its square. The change-of-basis matrix
evaluated at the four roots turns the Vandermonde determinant in r_i
into the one in q(r_i).

## 2. The B!=0 case: one zero and two asymptotic orders

The norm of s is the nonzero constant -B^2. For alpha!=0,

    Norm(s-alpha)=H_6(alpha)
                 =4*alpha^2*(U-U(alpha)).                  (6)

Factor a nonzero Laurent polynomial into a power of s and linear
factors s-alpha. Equations(2),(5),(6) show that g has exactly ONE
nonzero zero, simple. The existing node pair t=1,-1 has depressed
sum s_0=b/2 and gives that zero. Here s_0!=0: if b=0 the original
quartic is even, and its depressed odd coefficient B is zero.
Consequently

    g(s)=lambda*s^k*(s-s_0), lambda in C*, k an integer.     (7)

Subtract successive polynomials in p from q until its degree n is
NOT divisible by four. These target shears preserve its algebra and
action. The resulting q is nonconstant and has n>=2, because a
linear polynomial cannot have equal values at the two distinct node
preimages.

For the monomial r^j, let h_j(s) be its divided difference. Put h_1=1,
h_2=s. The universal recurrence and formula are

    h_j=s*h_(j-1)-v*h_(j-2),
    h_j=sum_l (-1)^l*binom(j-1-l,l)*v^l*s^(j-1-2l),          (8)

where 0<=l<=floor((j-1)/2). At infinity the leading coefficient of
h_j at s^(j-1) is

    ((1+i)^j-(1-i)^j)/(2^j*i),

nonzero exactly when j is not divisible by four. Every other exponent
drops by 2a+3d with a,d>=0: replacing a leading s^2 part of v by A
drops two and replacing it by B/s drops three. In particular h_j has
NO term s^(j-2).

It follows that g has infinity order n-1. Equation(7) therefore gives
k=n-2>=0. Thus g has NO pole at s=0.

On the other hand, the highest-l term in(8), with B!=0, gives

    h_(2m+1)=(-B/2)^m*s^(-m)+higher powers,
    h_(2m)=m*(-B/2)^(m-1)*s^(2-m)+higher powers.              (9)

If n>=3 is odd, the leading monomial of q contributes a unique
largest pole at zero. No lower monomial can cancel it, a contradiction.

If n is even, it is 4m+2 because n is not divisible by four. Suppose
n>=6. The coefficient of r^(n-1) would contribute a pole of order
(n-2)/2, larger than the poles from every other monomial of q,
including r^n. Hence that coefficient is zero. The s^(n-2) coefficient
of g is now zero: h_n has no such term, the r^(n-1) term is absent,
and lower monomials have infinity order at most n-3.

But(7) gives its s^(n-2) coefficient as -lambda*s_0!=0. This final
contradiction proves n=2. Thus q modulo C[p] is quadratic.

## 3. The B=0 case: both pair components must be retained

Now p=r^4+A*r^2+C is even in r. The unordered-pair cover splits,
so the degree-six Laurent norm argument must NOT be applied.
Write the unique D-remainder

    q=H(U)+beta(U)*r+gamma(U)*r^2+delta(U)*r^3.

For the two antipodal pairs let w=r^2. Their cover is

    U=w^2+A*w+C,
    O(w)=beta(U(w))+w*delta(U(w)).

For the other four pairs, the cover and divided difference are

    U=C-(s^2+A)^2/4,
    G(s)=beta(U(s))+s*gamma(U(s))+(s^2-A)*delta(U(s))/2.

The full Vandermonde product is EXACTLY

    Delta=Norm_(C(w)/C(U))(O)*Norm_(C(s)/C(U))(G).             (10)

Both O and G are nonzero polynomials; otherwise Delta=0.
The norm-degree formula for a polynomial map gives

    1=deg_w O+deg_s G.                                     (11)

If deg_w O=0, its two possible nonconstant summands have distinct
even and odd degrees, so delta=0 and beta is constant. Then
G=beta+s*gamma(U(s)); degree one forces gamma to be a nonzero constant.

If deg_w O=1, the same parity argument makes delta a nonzero constant
and beta constant. But G cannot be constant: its nonzero
delta*s^2/2 term cannot cancel s*gamma(U(s)), whose degree is either
one or at least five. This contradicts deg_s G=0.

Only the first alternative remains. Therefore q modulo C[p] is
quadratic in r in this case as well.

For an exact independent check of(10), with d=C-U, its two factors are
the norm of beta+delta*w under w^2+A*w+d and the norm of
beta+gamma*s+delta*(s^2-A)/2 under
s^4+2*A*s^2+A^2-4*d. Their product is the determinant of the quartic
basis inclusion with B=0. The accompanying checker verifies this
universal polynomial identity, not merely a list of quartics.

## 4. Return to the actual node and compute the action

In both cases q=H(p)+a_2*r^2+a_1*r+a_0 with a_2!=0.
The node preimages in the r coordinate are b/4+1 and b/4-1.
Their equal q values force a_1=-a_2*b/2. Returning to t and absorbing
a constant into H gives

    q=H(p)+a_2*t^2.

If b=0, both p and q would be even polynomials in t, contradicting
Frac(C[p,q])=C(t). Hence b!=0. Direct integration now yields

    integral_(-1)^1 p*dq=-8*b*a_2/15 != 0.                  (12)

Rescaling back from the monic p only multiplies this nonzero value
by a nonzero constant.

This also gives a concrete target reduction: after removing H(p),
subtract the even polynomial in t^2 from p. The remaining generator is
b*(t^3-t). Thus the quartic-generator class reduces by polynomial
target automorphisms to the already proved quadratic/cubic class.
No such reduction is asserted for arbitrary higher-degree generators.

## 5. Coprime degrees and the prime-degree filter

There is a separate infinite-degree restriction that does not depend
on the quartic calculation. Let A_0=C[p,q], with degrees m,n.
Both m,n>=2, because A_0 contains no nonconstant linear polynomial
and cannot be generated by one element.

The element q is integral over D=C[p]. Since q generates C(t) over
C(p), its monic minimal polynomial has degree m. Its coefficients
belong to D: they lie in C(p), are integral over D, and D is normal.
Thus 1,q,...,q^(m-1) is a D-basis of A_0. As before, the inclusion
determinant into the basis 1,t,...,t^(m-1) of C[t] has degree ONE:

    deg_U Delta=dim_C(C[t]/A_0)=1.                         (13)

Suppose gcd(m,n)=1. Let a and b be the nonzero leading coefficients
of p and q. In formal Puiseux series at U=infinity, the m roots of
p(t)=U have leading terms

    t_i=a^(-1/m)*zeta_i*U^(1/m),

where zeta_i are the distinct m-th roots of unity. Indeed, substituting
t=U^(1/m)*w reduces the equation at infinity to a*w^m=1, whose
roots are simple and hence lift to formal series. The Vandermonde
ratio for Delta is

    Delta(U)=product_(i<j) (q(t_j)-q(t_i))/(t_j-t_i).         (14)

The numbers zeta_i^n are distinct because gcd(m,n)=1. Hence every
numerator has its stated nonzero leading term: there is NO leading
cancellation in(14). Counting its powers of U gives

    deg_U Delta=(m-1)*(n-1)/2.                             (15)

Consequently coprime generator degrees can only be (2,3) or (3,2).

For completeness, the same leading coefficient is visible without
choosing roots of unity. In the monomial model p=t^m,q=t^n, the j-th
basis column, 0<=j<m, has its only nonzero entry in row n*j mod m,
equal to U^floor(n*j/m). For gcd(m,n)=1 those rows are a permutation.
Its determinant is a signed monomial and

    sum_(j=0)^(m-1) floor(n*j/m)=(m-1)*(n-1)/2.              (16)

This is also an exact control of the nonzero asymptotic coefficient
used in(15).

Now suppose m is PRIME and m>=5. Successively subtract scalar
multiples of powers of p from q whenever its degree is divisible by
m. The process strictly decreases the degree and ends with a
generator q_0 of degree n_0 not divisible by m. It cannot end with
a constant: A_0=C[p] would then be normal. It cannot end at degree
one because q_0 belongs to A_0. Thus n_0>=2 and gcd(m,n_0)=1,
contradicting(13)-(15).

Therefore NO polynomial generator of A_0 has prime degree at least
five. This excludes infinitely many degrees with NO bound on the
other generator. It does not exclude composite degrees sharing a
nontrivial divisor after all leading target shears.

## 6. Consequence for actual degree-two Keller conductors

In the global conductor-degree-two case, every transverse source
fibre has the EXACT nodal algebra generated by the two target
restrictions. If either output has directional source degree at
most FOUR, every fibre has a generator of degree two, three or four.
The preceding proof and the earlier low-degree theorem make its
action nonzero on EVERY fibre.

For an actual constant-Jacobian pair the same endpoint identities
give I'(y)=-2*Jac, so that polynomial action has a complex zero.
This is impossible. Accordingly no noninvertible plane Keller map
in the degree-two conductor class has such an output, even after
a polynomial target automorphism.

The prime-degree filter adds the following restrictions. On EVERY
transverse fibre, neither output restriction can have prime degree
at least five, and their degrees cannot be coprime except for the
pair (2,3) or (3,2). In particular, after ANY polynomial target automorphism,
the global directional degree of either output is at least SIX:
degrees at most four were just excluded, and degree five would give
a forbidden prime-degree generator on generic fibres.

The affine action I has exactly one zero. At that fibre, for every
polynomial target automorphism, both generator degrees must be
COMPOSITE and at least six, and they must have a common divisor
greater than one. Indeed action scales by the nonzero target
Jacobian, while degrees two, three and four have nonzero action
and all higher prime degrees are impossible.

These statements allow the directional degree to drop on special
fibres, but a drop to a prime degree at least five is forbidden.
The argument has no bound on the degree of the other output.
The composite common-divisor cases are settled by the separate
all-degree nodal-action theorem. This elementary argument does not
prove automatic trace integrality, HC4 or JC2.

## 7. Verification scope

The proof uses the nodal index of dimension one, Vandermonde
determinants, rational polynomial norms, and the written recurrence
and asymptotic arguments. There is no imported nodal-embedding
classification or conditional plane-curve theorem.

Exact controls verify the universal quartic inclusion/norm identity,
the separate split-pair factorization, the Laurent norm factors,
the monomial recurrence and its pole coefficients, and quartic
generator actions with arbitrarily higher added polynomials in p.
They also check exact monomial inclusion determinants, root-of-unity
resultants, and the coprime degree formula in finite ranges.
Finite controls support the formulas; the all-degree exclusion
uses the written highest-pole and coefficient arguments.
Independent mathematical review, formal verification and external
review remain absent.
