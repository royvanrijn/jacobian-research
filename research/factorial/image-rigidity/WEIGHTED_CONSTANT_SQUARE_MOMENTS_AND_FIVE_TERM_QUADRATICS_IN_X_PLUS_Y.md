# Weighted constant-square moments and five-term quadratics in x+y

Put Lambda(X^j)=j! and define, for n>=0,

 H_n(a)=Lambda((1+aX^2)^n)=sum_(i=0)^n binom(n,i)(2i)! a^i,
 K_n(a)=H_n(a)+aH_n'(a)
       =sum_(i=0)^n binom(n,i)(i+1)(2i)! a^i.                (1)

In particular K_0=1, K_1=1+4a, K_2=1+8a+72a^2.

**Theorem 1.** For every m>=1, K_m and K_(m+1) are coprime over Q.
Every K_n, n>=1, is squarefree. Consequently neither assertion acquires
a common or repeated root over any characteristic-zero extension.

**Corollary 2 (actual five-term quadratic in x+y).** For every characteristic-zero
field and all nonzero a,d,k, the phase

 f(x,y)=a(x+y)^2+d x+k

has a nonzero factorial moment among every five consecutive positive indices.
It has exactly five monomials. Interchanging x and y gives the same statement
with only the y-linear term. This does not prove arbitrary rank-one orientation,
general binary quadratic SFC, full SFC(2), or SFC(3).

The proof is written below. No finite gcd experiment, irreducibility conjecture,
or previously proved all-unit factor theorem is a premise. The even-coordinate
divided-power factorization needed here is proved explicitly.

## 1. Exact polynomial identities

For n>=1, binomial differentiation gives

 aH_n'=n(H_n-H_(n-1)),
 K_n=(n+1)H_n-nH_(n-1),
 aK_n'=n(K_n-K_(n-1)).                                    (2)

The native factorial recurrence is

 H_(n+1)=2a(n+1)(2n+1)H_n-4a n(n+1)H_(n-1)+1.              (3)

It follows either by direct coefficients, or from the one-variable factorial
derivative identity with linear coefficient0 and constant1.

Here is another direct coefficient identity. If H_n=sum c_i a^i, then
c_(i+1)=2(n-i)(2i+1)c_i. With theta=a(d/da), this gives
H_n-1=2a(n-theta)(2theta+1)H_n. Expanding yields

 4a^3 H_n''+(6-4n)a^2 H_n'+(1-2na)H_n=1.                  (4)

Using K_n=H_n+aH_n' in(4) gives the exact identity

 4a^2 K_n'=1+(4n+2)aK_n-[1+2(n+1)a]H_n.

Also (2)--(3) give, without a coefficient-dependent division,

 (n+1)[1+2(n+2)a]H_n
   =(n+2)+4a(n+1)(n+2)K_n-K_(n+1).

Eliminating H_n produces the universal polynomial identity

 4a^2(n+1)[1+2(n+2)a]K_n'
   =[1+2(n+1)a]K_(n+1)
     -2a(n+1)[3+2(n+2)a]K_n-1.                            (5)

In particular at a common zero of K_n,K_(n+1), both a and
1+2(n+2)a are nonzero and

 K_n'=-1/[4a^2(n+1)(1+2(n+2)a)].                           (6)

The nonzero assertion for a also follows from K_n(0)=1.
No root is discarded by dividing a potentially zero pivot: (5) first rules
out that possibility.

## 2. The even-index reciprocal and its exact divided-power reduction

Fix an EVEN integer n>=2. Write s_2 for binary digit sum and put

 L_n=2^{s_2(n)}(n+1)(2n)!/4^n.

Legendre's formula shows L_n is an odd rational integer. Therefore

 P_n(V)=2^{s_2(n)}K_n(V/4)/L_n

is monic. Take its reciprocal F_n(Z)=Z^nP_n(1/Z), with F_n(0)=1.
For eta_j(Z)=2^j Z^j/j!, its exact coefficient formula is

 F_n=sum_(j=0)^n b_j eta_j,
 b_j=(n-j+1)/[(n+1)D_(n,j)],
 D_(n,j)=prod_(r=0)^(j-1)(2n-2r-1).                       (7)

This follows immediately by taking the coefficient of V^(n-j) in P_n.
Every denominator in(7) is odd. Since n is even,

 b_j in Z_2,  b_j=1 modulo2 if j even,
                 b_j=0 modulo2 if j odd.                 (8)

This is NOT the all-coordinate-unit hypothesis of an earlier theorem.
Its finite-degree lifting argument is supplied next.

## 3. Fixed-degree Hensel lifting for the even-coordinate pattern

Let e=2^rho>=2 be the highest bit of the current even degree, and write
that degree as e+S with S even, 0<=S<e. Work first modulo2 in the
finite divided-power algebra with degrees below e.
Multiplication is eta_i eta_j=binom(i+j,i)eta_(i+j).
A low product crossing degree e vanishes modulo2 by its binary carry;
eta_e eta_j=eta_(e+j) modulo2 for j<e.

For a polynomial F having reduction(8), decompose

 Fbar=F_lo+eta_e F_hi.

Here F_hi is the polynomial sum of all even eta_j through degree S.
It has constant1 and unit top coordinate. Choose

 Ubar=F_hi,   Gbar=eta_e+F_lo/Ubar,

where low division is in the truncated degree-<e algebra.
The positive-degree ideal is nilpotent, so the constant-one Ubar is
invertible there. Both low polynomials have support only in even indices.
The constant of Gbar is1 and its top divided-power coordinate is1.

The fixed-degree derivative of multiplication is invertible over F_2.
Indeed, constrain the correction D to G to have degree<e, and let the
correction V to U have degree<=S. Its high component is exactly V:
the product eta_e V has odd multiplication coefficients, while all products
of two low terms crossing e have even coefficients. Its low component is

 D Ubar+G_lo V.

Given an error E_lo+eta_e E_hi, solve explicitly

 V=E_hi,  D=(E_lo-G_lo V)/Ubar.

Hence the finite coefficient correction map is invertible.
Lift successively modulo2,4,8,..., keeping the degrees bounded by e and S
and the top divided-power coordinate of G fixed1. Corrections solving
each residual exist by this same invertible reduction map. The coefficient
vectors converge in their fixed finite Z_2 lattices. The lifted factors
have constant units and reductions Gbar,Ubar. The residual factor U has
the SAME even-coordinate pattern through degree S, so the argument recurses.
If S=0 it is only a unit scalar, and the recursion stops.

After normalizing each factor to constant1, we obtain a factor Q_e of
each binary degree e>=2 occurring in n. There is no degree-one packet
because n is even. Its ordinary coefficients satisfy

 Q_e(0)=1, v([Z^j]Q_e)>=s_2(j), 1<=j<e,
 v([Z^e]Q_e)=1.

The last equality follows from its fixed unit top divided-power coordinate
and v(2^e/e!)=1. Constant-unit normalization preserves these bounds.
Taking monic reciprocals gives

 P_n=prod_e g_e,  g_e(V)=V^e Q_e(1/V),
 g_e monic in Z_2[V],  v(g_e(0))=1,
 v([V^j]g_e)>=s_2(e-j), 0<=j<e.                            (9)

The ordinary Newton polygon is the single edge from(0,1) to(e,0).
Thus every root of g_e has valuation1/e. This also shows g_e irreducible
over Q_2: a nonconstant monic Q_2 factor of degree d would have constant
valuation d/e, an integer, so e divides d. No proper degree is possible.
Characteristic zero makes it separable. Distinct packets have distinct
root valuations and share no roots. This proves that every EVEN K_n is
squarefree, with all root packets and coefficient bounds stated in(9).

## 4. An even upper index in an adjacent pair

If K_(n-1)(a)=K_n(a)=0 with n even, then a!=0 and(2) gives K_n'(a)=0.
This contradicts the even-index squarefreeness just proved.
Thus an adjacent pair cannot have its EVEN index as the upper index.

## 5. An even lower index: the forced derivative contradiction

Suppose now K_n(a)=K_(n+1)(a)=0 with n EVEN. Work over Q_2bar.
Put beta=4a. By(9), beta belongs to exactly one packet g_e and

 v(beta)=1/e,  v(a)=-2+1/e.

Since n+2 is even,

 v(2(n+2)a)=v_2(n+2)-1+1/e>0.

Therefore 1+2(n+2)a is a unit; n+1 is an odd unit as well.
Formula(6) forces

 v(K_n'(a))=2-2/e,
 v(P_n'(beta))=s_2(n)-2/e.                               (10)

Write n=D+e+S, where D is the sum of higher binary degrees, S=n mod e is
even, h=s_2(D), and ell=s_2(S). At beta, each lower-degree packet g_f
has value of valuation f/e, and each higher packet has value of valuation1.
This follows from the distinct root valuations: every difference beta-alpha
has valuation min(1/e,1/f). Thus the total valuation of all other factors is
S/e+h. Since s_2(n)=h+1+ell, (10) gives

 v(g_e'(beta))=ell+1-(S+2)/e.                              (11)

All coefficients of g_e lie in Q_2 and have integer valuations.
In g_e'(beta)=sum_(j=1)^e j c_j beta^(j-1), distinct nonzero summands have
distinct fractional parts (j-1)/e modulo1. Hence the least summand is
unique and cannot cancel.

Because 0<=S<=e-2 is EVEN, the only index with the fractional part in(11)
is

 j_0=e-S-1,    1<=j_0<=e-1,   j_0 odd.

It has v_2(j_0)=0. Equality with(11) would therefore force
v(c_(j_0))=ell. But(9) requires

 v(c_(j_0))>=s_2(e-j_0)=s_2(S+1)=ell+1,

because adding1 to the even S creates a new final binary bit.
This contradiction includes e=2,S=0 and the endpoint S=e-2.
No derivative cancellation over a ramified coefficient field is being
ignored: the packet coefficients themselves lie in Q_2.

Sections4--5 exclude BOTH orientations of every adjacent pair.
If an adjacent common root existed in any characteristic-zero field,
it would be algebraic over Q, since K_m has nonzero constant leading
coefficient (m+1)(2m)!. Its number field embeds into Q_2bar, where the
argument just given excludes it. This proves adjacent coprimeness in
Theorem1 over every characteristic-zero extension.

Finally, a multiple root of ANY K_n is nonzero and, by(2), would also be a
root of K_(n-1). This has just been excluded (K_0=1 handles n=1).
Thus all K_n are squarefree, including the odd indices. No irreducibility
claim for the odd K_n is made.

## 6. The actual five-term factorial-moment consequence

Let f=a(x+y)^2+d x+k with a,d,k nonzero, and suppose five consecutive
moments vanish starting at N>=1. Put

 G_l=Lambda_1((aX^2+dX+k)^l),
 H_l=Lambda_1((aX^2+k)^l).

The polynomial boundary identity is

 (r+1)d Lambda_2(f^r)=G_(r+1)-H_(r+1).

It follows from the difference of the x and y Stein identities, just as
in the six-term rank-one argument, with the second linear coefficient0.
Thus G_l=H_l=R_l for l=N+1,...,N+5.

Subtracting their native one-variable recurrences at index l+2 gives

 (l+1)d R_l=-k^(l+1),    l=N+1,N+2,N+3.

Scale f by k^(-1), and write its normalized coefficients again as a,d.
Now k=1 and C=-1/d is nonzero. The constant-square sequence in(1) satisfies

 H_(m-1)=C/m, H_m=C/(m+1), H_(m+1)=C/(m+2),
                         m=N+2>=3.

Using(2), these values give K_m=K_(m+1)=0, contradicting Theorem1.
This proves Corollary2 for all positive starts over all characteristic-zero
fields. The exact support is five: x^2,xy,y^2,x,1 have nonzero coefficients.
No shorter window is inferred by deleting another coefficient.

The proof is self-contained. The polynomial and differential identities,
and the even-coordinate lifting proved above, are its mathematical inputs.
Finite native gcd or factorization checks are supplementary only.
The conclusion is restricted to the stated weighted moments and the actual
five-term family; it does not settle general binary quadratic SFC, SFC(2),
or SFC(3).
