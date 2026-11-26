# A uniform conjugate-part congruence for the odd native high block

This proof establishes a native coefficient congruence. It does not prove the
formal-quotient or actual Newton-factor conjugate-part estimate still needed to
exclude the half-valuation three-forcing branch.

Let M=1 or2 modulo4, M>=1, and let z be a positive-valuation root of

 2(z+M)^3+3(z+M)^2+M=0,

with v_2(z)=1/2. Write K=Q_2(z), let sigma be its nontrivial automorphism,
and Anti(x)=(x-sigma(x))/2. Select p=M if M is odd and p=M-1 if M is even,
so p=1 modulo4. Write p=e+S, where e=2^rho>=4 is its highest binary bit;
then 0<S<e and S=1 modulo4.

Let T_n(r,z)=Lambda((rX^2+X+z)^n)+n Lambda((rX^2+X+z)^(n-1)). Define the
ordinary normalized reciprocal

 Q_n(Z)=(4Z/z)^n T_n(z/(4Z),z)/(2n)!
       =sum_(k=0)^n b_k^(n) eta_k, eta_k=2^k Z^k/k!.

Its constant coefficient is1. Put

 c=[prod_(a=S+1)^p (2a-1)]^(-1),
 E_j=(2^j/j!)[b_(e+j)^(p)-c b_j^(S)], 0<=j<=S.

**Theorem 1 (native high-block congruence).** Every E_j satisfies

                 v_2(Anti(E_j)) >= rho+3/2.                         (1)

The coefficient polynomials use the same z in both degrees p and S. The
rational unit c is indispensable in comparing these two divided-power blocks.
The proof is uniform in M,e,S,j and uses no finite-prefix inference.

## 1. The ramified field and powers of z

The cubic has a width-one negative Newton edge and a width-two positive edge.
Equivalently, after setting w=2z its reduction has the simple unit root1.
Hensel lifting gives the negative Q_2 root and leaves a ramified quadratic
factor. Thus K is quadratic. Put t=Tr(z), D=-Norm(z), so z^2=t z+D,
and put pi=z-t/2. Then pi^2=D+t^2/4 and v(pi)=1/2.

Let W=M(M+1)/2, an odd integer, and u_0=-W/3. Directly from the cubic,
with u=z^2/2,

 (2M+1)(3u+W)+2z(u+3W)=0,
 u-u_0=-16Wz/[3(3(2M+1)+2z)].                                      (2)

Its valuation is9/2. In u-u_0=(D/2-u_0)+(t/2)z the first summand has integer
valuation and the second half-integer valuation. Therefore v(t)=5. Comparing
coefficients of the cubic and its quadratic factor gives

 D=2u_0/[1+2t/(3(2M+1))].                                         (3)

In particular D/(2u_0)-1 has valuation>=6 and pi^2/(2u_0)-1 has valuation>=6.

Write z^(-1)=a+b pi, with a=-t/(2D), b=1/D. Thus v(a)=3 and v(b)=-1.
Expanding its powers proves

 l odd:  v(Anti(z^(-l)))=-l/2,
          Anti(z^(-l))=pi/(2u_0)^((l+1)/2) * (1+epsilon_l),
          v(epsilon_l)>=6;
 l even: v(Anti(z^(-l)))>=-l/2+9/2.                                (4)

Indeed, for odd l the term (b pi)^l is uniquely least; every term involving
at least two a's has a relative gain>=7. Equations(3) supply the stated
relative error6. For even l, an odd number h of a's is required for an
anti term. At h=1 its gain over -l/2 is at least 7/2+v_2(l)>=9/2;
at h>=3 the gain is at least21/2. These estimates hold for all l, however large.

## 2. Native coefficient collapse and the binomial difference

For 0<=k<=n set i=n-k. Literal multinomial expansion followed by combining the
T addition gives

 b_k^(n)= [1+sum_(l=1)^k 2(i+l)binom(k,l)(2i+1)_(l-1)z^(-l)]/D_(n,k),
 D_(n,k)=prod_(a=0)^(k-1)(2n-2a-1).                                (5)

All D_(n,k) are rational odd units. The notation (x)_h denotes a rising
factorial. The formula also holds for k=0 with empty sum.
For k=e+j in degree p and k=j in degree S, the value i=S-j is identical and
D_(p,e+j)=c^(-1)D_(S,j). Thus, up to a common rational unit, E_j is

 (2^j/j!) sum_(l>=1) 2(i+l)(2i+1)_(l-1) Delta_l z^(-l),
 Delta_l=binom(e+j,l)-binom(j,l).                                  (6)

Binomial coefficients outside their range are zero. Vandermonde and
v_2 binom(e,a)=rho-v_2(a), 1<=a<=e, imply

 v_2(Delta_l)>=rho-floor(log_2 l).                                 (7)

When j=0 and l is odd, the sharper value v_2(Delta_l)=rho applies whenever
that term is nonzero. The ordinary scale 2^j/j! contributes s_2(j).
An interval of l-1 integers has factorial valuation at least v_2((l-1)!),
so before the extra anti gain in(4) a summand of(6) has valuation at least

 rho+s_2(j)+l/2-s_2(l-1)-floor(log_2 l).                            (8)

## 3. All large shifts and the harmless small shifts

For every even l>=2,

 l/2-s_2(l-1)-floor(log_2 l)>=-2.

This follows by direct inspection at l=2,4,6,8,10,12,14; for l>=16 use
s_2(l-1)<=1+floor(log_2 l). The extra9/2 in(4) makes every even-l anti
summand in(6) have valuation at least rho+5/2, stronger than(1).

For every odd l>=9,

 l/2-s_2(l-1)-floor(log_2 l)>=1/2.

Check l=9,11,13,15 directly. If r=floor(log_2 l)>=4, l-1 is even,
s_2(l-1)<=r, and l/2>=2^(r-1)+1/2>=2r+1/2.
For j>0, s_2(j)>=1, so(8) proves(1). For j=0 use the sharper odd-l
binomial valuation following(7), which also proves(1). Thus only l=1,3,5,7
need closer attention.

For l=1, Delta_1=e and the valuation is
rho+1/2+v_2(i+1)+s_2(j)>=rho+3/2:
if j=0 then i=S is odd, and otherwise s_2(j)>=1.

For the remaining odd shifts, using(4), their native anti gains before
Delta_l and the ordinary scale are respectively

 l=3: 1/2+v_2(i+3)+v_2(i+1),
 l=5: 1/2+v_2(i+5)+v_2(i+1)+v_2(i+2),
 l=7: 1/2+v_2(i+7)+v_2(i+1)+v_2(i+2)+v_2(i+3).                    (9)

For l=3 and i odd the gain is>=7/2. For i even, j is odd; if s_2(j)>=2,
(7) and(9) prove(1). Hence l=3 can fail the separate-term bound only at j=1.
For l=5 and i odd the gain is>=5/2; j is then even, and if j=0 use its
sharper binomial bound, otherwise s_2(j)>=1. For i even the gain is>=3/2;
again only j=1 is potentially bad. For l=7 and i odd the gain is>=9/2.
For i even, the gain is>=3/2 and only j=1 might be bad using(7), but

 Delta_7=binom(e+1,7)=binom(e,6)+binom(e,7)

has valuation rho-1 whenever present (e>=8), proving(1) in that case too.

## 4. The only two low terms cancel at j=1

It remains to combine l=3 and l=5 for j=1. Here i=S-1=0 modulo4. By(4),
modulo errors at least six valuation units higher, their anti sum divided
by e*pi is A+B, where

 A=(i+3)(i+1)(2i+1)/u_0^2 * (e^2-1)/6,
 B=(i+5)(i+1)(i+2)(2i+1)(2i+3)/u_0^3
      * (e+1)(e-1)(e-2)(e-3)/120.                                (10)

Both have valuation -1, but 2A and 2B are odd 2-adic integers. Therefore
v_2(A+B)>=0. The ordinary coefficient scale at j=1 is2, so the combined
anti term has valuation at least rho+3/2. The relative errors>=6 in(4)
are individually far above that threshold and cannot spoil it. This treats
the only failure of the separate-term bounds. Summing all terms proves(1).

## 5. Scope

This theorem supplies an all-range version of the high-block diagnostic.
It does not itself show that the formal quotient Q_p/Q_S has the desired
anti congruence, nor that this congruence survives degree-bounded Newton
factorization at the necessary weighted precision. Those implications need
separate proofs; no three-zero forcing or new Strong Factorial theorem is
claimed here. The exact local branch exclusions and finite exceptional-pair certificates
are independent of this native estimate.
