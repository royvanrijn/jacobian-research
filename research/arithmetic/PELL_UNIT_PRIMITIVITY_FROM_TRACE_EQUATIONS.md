# Fundamental Pell units from shifted power traces

An integral unit supplied by a norm identity can be a proper power in
the specialized field. For the two Pell constructions below, the trace
equation forces every fixed prime divisor of a proper-power exponent
into a finite Diophantine set. Counting the large exponents then gives
almost injective specialization in the parameter.

## A trace criterion for almost all units to be fundamental

**Theorem PF1.** Fix \(\sigma=\pm1\), a nonconstant polynomial
\(R\in\mathbb Q[T]\), an integer \(r\ge2\), and
\(\lambda,\kappa\in\mathbb Q\) with \(\lambda\ne0\).
Suppose that \(r\ge3\) when \(\sigma=1\) and that
\(\kappa\ne2,-2\) in that case. For each parameter \(j\) in any subset of
the positive integers, suppose there is a real quadratic field \(K_j\)
and a unit \(U_j>1\) of norm \(\sigma\) with
\[
       \operatorname{Tr}(U_j)=\lambda R(j)^r+\kappa.
       \tag{1}
\]
If \(\sigma=1\), let the reference unit be the least totally positive
unit greater than one. If \(\sigma=-1\), let it be the least positive
unit greater than one. Then for every \(\eta>0\), the number of \(j\le T\)
for which \(U_j\) differs from the reference unit is
\(O_{\lambda,\kappa,R,r,\eta}(T^\eta)\), with an effective bound.
On an unbounded parameter set, \(K_j\) is therefore distinct on the
complement of these exceptions after finitely many initial parameters
are omitted.

Write \(U_j=\epsilon_j^{h_j}\) with \(h_j\ge1\) and with
\(\epsilon_j\) the reference unit. In the negative-norm case
\(\epsilon_j\) has norm \(-1\) and \(h_j\) is odd. If a prime \(p\)
divides \(h_j\), let \(V=\epsilon_j^{h_j/p}\) and \(z=\operatorname{Tr}(V)\).
Then \(z\in\mathbb Z\), \(N(V)=\sigma\), and
\[
              D_p(z,\sigma)-\kappa=\lambda R(j)^r.               \tag{2}
\]
The Dickson trace polynomial is characterized by
\(D_p(u+\sigma/u,\sigma)=u^p+(\sigma/u)^p\).
Its finite critical values belong to
\(\{2(\sqrt{\sigma})^p,-2(\sqrt{\sigma})^p\}\).
Indeed the zeros of the derivative of \(u^p+(\sigma/u)^p\)
have \(u^{2p}=\sigma^p\). At the branch points \(u^2=\sigma\) of
\(u+\sigma/u\), the derivative quotient extends with nonzero value
\(p^2u^{p-1}\) by cancellation. Thus \(D_p(z,\sigma)-\kappa\) is
squarefree: for positive norm this uses \(\kappa\ne\pm2\), and for
negative norm it uses odd \(p\) and rational \(\kappa\).

For a fixed \(p\), (2) is an \(S\)-integral superelliptic equation
in \(z,R(j)\), over \(\mathbb Q\), with \(S\) consisting of the
infinite place and the finitely many denominator primes. For \(r\ge3\)
it has degree \(p\ge2\), and for \(r=2\) the hypotheses give \(p\ge3\).
[Bérczes–Evertse–Győry, Theorems 2.1 and 2.2](sources/berczes-evertse-gyory-hyper-superelliptic-excerpt.pdf)
therefore bound both heights effectively. There are only finitely
many corresponding \(j\) because \(R\) is nonconstant.

For the tail put \(d=r\deg R\). Up to \(j=T\), we have
\(U_j\le C(2T)^d\), from (1) and \(U_j+\sigma/U_j=\operatorname{Tr}(U_j)\).
Every \(\epsilon_j\) is at least \((1+\sqrt5)/2\) by its positive integral
trace, so \(h_j=O(\log(2T))\). For \(p>B\) dividing \(h_j\), the positive
integer \(z\) satisfies \(z\ll T^{d/B}\). Each pair \(p,z\) can give
at most \(d\) roots \(j\) of (1) with its prescribed trace. There are
\(O(T^{d/B}\log(2T))\) pairs. Take \(B>2d/\eta\) and include the
finite sets for \(p\le B\). The two bounds give the claimed estimate.
Finally, (1) is eventually strictly monotone. Two reference units in
the same field are equal, so their traces and large parameters are equal.
This argument measures exceptions by parameter height. The bound also
applies after intersection with any prime progression or a fixed
positive-density Frobenius condition.

## Quartic units and their quadratic and elliptic counts

**Corollary PF2.** Let \(1\le j\le T\) run through the parameter of
[CP1](QUARTIC_PELL_TWO_CLASS_DIRECTIONS.md). For every \(\eta>0\),
all but \(O_\eta(T^\eta)\) parameters give the fundamental positive
unit of \(K_j\), of norm \(-1\).
The same number of parameters suffices to delete all repetitions
and initial exceptions. In particular CP1 produces
\(\gg X^{1/8}\) distinct real quadratic fields of discriminant at most
\(X\) with its two ordinary class directions and its wild unit.
The [PM2](PELL_NORM_SECTIONS_MORDELL_RANK.md) Mordell curves give
\(\gg X^{1/16}\) distinct rational isogeny classes of conductor at
most \(X\), each of rank at least three. On the almost injective subset,
the logarithm of the fundamental unit is \(4\log j+O(1)\).

In CP1 use \(t=(N_0+Lj)/128\). If \(P=a t^4+b t^2+c\) then
\[
      2P=2a(t^2+b/(2a))^2+2c-b^2/(2a).
\]
This is (1) with \(\sigma=-1\), \(r=2\) and a quadratic \(R\).
The unit is positive and greater than one for sufficiently large \(j\).
Its trace has growth of order \(j^4\) and
\(\operatorname{Disc}(K_j)\le4(P(t)^2+1)\ll j^8\).
The PM2 conductor bound is \(O(j^{16})\), and its positive
\(j\)-invariant-zero models from distinct quadratic fields are
nonisogenous over \(\mathbb Q\).
The stronger parameter count requires only these height upper bounds.
It leaves the distribution of the square divisors of \(P(t)^2+1\)
and the multiplicities of equal conductors to separate estimates.

## Moving prime classes on an almost injective prime slice

**Corollary PF3.** Fix \(S,(a_p),B,c,m\) as in
[PI3](PELL_PERIOD_RING_CLASS_INERTIA.md) or
[PM3](PELL_NORM_SECTIONS_MORDELL_RANK.md).
On the respective prime progression one may omit the auxiliary
Chebyshev nonpower congruences and retain all their conclusions
after deleting \(O_\eta(T^\eta)\) exceptions among primes \(q\le T\).
In [AW3](PELL_PRIME_PATH_WILD_AFFINE_FIELDS.md) and PM3 the same
deletion works on each prescribed translation Frobenius slice.
On the retained set the displayed norm-one unit is the fundamental
totally positive unit, and the fields for distinct sufficiently
large \(q\) are distinct.

Consequently PI3, AW3 and PM3 give
\(\gg X^{1/(2m)}/\log X\) distinct real quadratic fields up to
quadratic discriminant \(X\). The degree-\(m\) affine fields on
each Frobenius slice give
\(\gg X^{1/(m(m-1))}/\log X\) distinct fields up to their absolute
discriminant \(X\).
The PM3 Mordell pair gives \(\gg X^{1/(4m)}/\log X\) distinct
rational isogeny classes of conductor at most \(X\).
On every retained slice the logarithm of the fundamental totally
positive unit is \(m\log q+O(1)\).

In the notation of those constructions the displayed unit has trace
\[
    2(1+B^2n)
       =\frac{B^2}{Bc-1}q^m
        +2-\frac{B^2c^2}{Bc-1}.
\]
Here \(m\ge3\). The distance of the constant term from \(2\) is
\(-(Bc)^2/(Bc-1)\), and its distance from \(-2\) is
\(-(Bc-2)^2/(Bc-1)\). Both are nonzero. Thus PF1 applies with
\(\sigma=1\), \(r=m\) and \(R(q)=q\). On the retained slice, comparing
the radical coefficients in the given unit and the fundamental unit
is now an equality. In particular the required valuation of the
fundamental unit is exact, and all the ideal-class and conductor
arguments from PI3 and PM3 apply. In PM3 the fundamental totally
positive unit is noncubic even in \(K^\times\): any cube root in
the field is integral with integral inverse, positive, and has norm one.

For AW3 and PM3 keep \(q\equiv1\bmod m\) and all the prescribed
local congruences other than the nonpower selectors. They form one
reduced residue class modulo a fixed \(M\) divisible by \(m\).
The Kummer element \(C=c(Bc-2)/(Bc-1)\) and its pole at \(s\nmid m\)
are unaffected. The maximal abelian subfield of the AW3 auxiliary
normal closure remains \(\mathbb Q(\zeta_m)\), so the cyclotomic
progression and the translation order are compatible. In fact the
conditional proportion of order \(e\) inside that progression is
\(\varphi(e)/m\): at the fixed cyclotomic Frobenius the \(m\)
translations are all available, and \(\varphi(e)\) have order \(e\).
Chebotarev counts unions of conjugacy classes. Hence the retained
count up to \(T\) is
\(\bigl(\varphi(e)/(m\varphi(M))+o(1)\bigr)\operatorname{Li}(T)\).
For PI3 without the Frobenius condition use the prime number theorem
for its fixed progression. The deleted set has smaller order in
both counts by PF1 with, for example, \(\eta=1/2\).

The underlying quadratic discriminant is \(O(q^{2m})\).
The affine exponent of that discriminant is \((m-1)/2\),
and all the extra wild factors are constant in the fixed slice.
The quadratic field is recovered from the degree-\(m\) affine field
by AW3's ramification criterion. The PM3 conductor is \(O(q^{4m})\).
These give the three height thresholds. The estimates keep all
parameters of the slice fixed except \(q\). They give no uniform
prime count when \(S\), the exponents or the auxiliary progression
vary with the height.

## Critical trace shifts and the remaining powers of two

**Corollary PF4.** In the positive-norm case of PF1 allow
\(\kappa=\pm2\), still with \(r\ge3\). If \(\kappa=2\), the
conclusion of PF1 continues to hold. If \(\kappa=-2\), then
outside \(O_\eta(T^\eta)\) exceptions the exponent of the
fundamental totally positive unit belongs to \(\{1,2,4\}\),
and exponent four can be removed from this list unless \(r=4\).
The constants are effective. There are infinite examples of
nonfundamental specialized units on the critical chart, including
fourth powers when \(r=4\).

For an odd prime \(p\), the factorization at either critical value
has one simple factor and all other factors double:
\[
       D_p(z,1)-2s=(z-2s)Q_{p,s}(z)^2,\qquad s\in\{1,-1\}.
       \tag{3}
\]
Here \(Q_{p,s}\in\mathbb Q[z]\) has distinct roots and
\(Q_{p,s}(2s)\ne0\). In the variable \(u\) with \(z=u+u^{-1}\),
the left side is \((u^p-s)^2/u^p\). Of the roots of \(u^p=s\),
only \(u=s\) is a branch point of \(u+u^{-1}\), and the other
roots pair by inversion. This gives the simple root at \(z=2s\)
and the \((p-1)/2\) double roots.

Enlarge the finite set of denominator primes to contain the
numerator primes of \(\lambda\) and the primes of the resultant in
(3). Outside that set the two factors cannot both have positive
valuation. With \(k=r/\gcd(r,2)\), (2) then implies
\[
       z-2s=a u_0^r,\qquad Q_{p,s}(z)=b v_0^k .
\]
The constants \(a,b\) range over finite effective lists of products
of primes in the set, using exponents reduced modulo \(r,k\)
and allowing the two signs. The variables \(u_0,v_0\) are
\(S\)-integers. Zero factors in (3) only give \(R(j)=0\),
a finite set. Substituting the first equation into the second gives
\(Q_{p,s}(2s+a u_0^r)=b v_0^k\), a polynomial of degree
\(r(p-1)/2\ge3\) with no repeated roots, since
\(Q_{p,s}(2s)\ne0\) and \(Q_{p,s}\) has simple roots.
The same height theorem as in PF1 applies, now with exponent \(k\).
All fixed odd prime divisors are thereby controlled, and the
large-prime tail bound of PF1 still applies. For \(\kappa=2,p=2\)
the polynomial \(D_2(z,1)-2=z^2-4\) also meets Theorem 2.1.

For \(\kappa=-2\), an exponent divisible by eight would give
\[
        (D_4(z,1))^2=\lambda R(j)^r,
        \qquad z=\operatorname{Tr}(U_j^{1/8}),
\]
where the eighth root is taken in the unit group.
Valuations again give \(D_4(z,1)=b v_0^k\) with \(b\) from a
finite list and \(v_0\) an \(S\)-integer. This is an effective
height problem since \(D_4\) has four distinct roots and \(k\ge2\).
With exponent divisible by four, the corresponding polynomial
is \(D_2(z,1)=z^2-2\), which has two distinct roots.
Here Theorem 2.1 suffices when \(k\ge3\), and the only case
\(k=2\) allowed by \(r\ge3\) is \(r=4\).

The exponent-two case is genuine. Take \(R(j)=j^2,\lambda=1\)
and let \(V_j>1\) be the norm-one unit of integral trace \(j^r\).
Then \(U_j=V_j^2\) has trace \(j^{2r}-2\).
For the fourth-power case use \(r=4,R(j)=j,\lambda=4\) and
take positive solutions \(z^2-2j^2=2\), for instance from
\(z+j\sqrt2=(2+\sqrt2)(3+2\sqrt2)^a\) for \(a\ge1\).
The norm-one unit with trace \(z\) has fourth power of trace
\((z^2-2)^2-2=4j^4-2\).
Both constructions give real quadratic fields for large parameters.
The fourth-power example only uses the displayed Pell sequence of
parameters, so the general exceptional-set bound gives no density
statement within that sequence.

## Verification

The [exact checker](verify_pell_unit_primitivity.py) verifies the
trace reparametrizations and the Dickson derivative identity on explicit
instances. The written argument proves the derivative assertion and the critical
factorization for all indices and uses the primary height theorem for
the fixed-index sets. It uses no enumeration of those integral points.
