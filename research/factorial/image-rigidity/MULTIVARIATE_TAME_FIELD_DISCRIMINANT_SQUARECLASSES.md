# Tame discriminant squareclasses from stepped factorial faces

Let \(P_n,P_N\) be adjacent ordinary binary factorial moments.
Use the [common-field discriminant ledger](MULTIVARIATE_COMMON_FACTOR_FIELD_DISCRIMINANT.md)
for their primitive integral squarefree parts. For each prime,
a shared branch is described by a root valuation \(a/b\)
in lowest terms and an irreducible stepped residual factor \(g\)
of degree \(f\). A local part of a common root field in this
branch has ramification index \(e=bu\) and residue degree
\(F=fv\). The prior valuation ledger bounds \(uv\) by the
residual multiplicity.

**Theorem.** In any such local part tame at an odd prime \(p\),
write \(\chi_p\) for the Legendre symbol on units. Its discriminant
valuation is \(F(e-1)\). After removing this power of \(p\), its
unit symbol is
\[
 (-1)^{e(F-1)}
 \chi_p\bigl((-1)^{e(e-1)/2+e-1}e^e\bigr)^{F}\ U,
 \qquad
 U=
 \begin{cases}
  1,& e\text{ odd},\\
  \chi_p((-1)^f g(0))^{v},& b\text{ even},\\
  \text{one of }1,-1,& e\text{ even and } b\text{ odd}.
 \end{cases}                                                    \tag{1}
\]
At \(p=2\) a tame part has odd \(e\), and the discriminant
unit squareclass is represented modulo eight by
\[
                 5^{F-1}\bigl((-1)^{(e-1)/2}e\bigr)^F .          \tag{2}
\]
These formulas augment the finite ledger at any chosen
primes. They can strictly improve the test at a ramified
prime. They give unramified splitting constraints also at primes
not dividing the gcd of the two array discriminants.

No assumptions on the count, heights, signs or gaps of the
original rows are needed. A separate corollary gives a
stronger denominator growth condition near a diffuse contact
when the possible quadratic root fields have bounded
discriminant.

## The residual unit and the relative discriminant

Let \(E=\mathbb Q_p(\alpha)\) be the candidate local field
of one root \(\alpha\). Use \(U_0\) for its maximal unramified
subextension and let \(\xi_0\) be an arbitrary uniformizer
of \(E\). Its residue field has size \(p^F\), the same as
that of \(U_0\). The residue \(\overline \xi\) of
\(\xi_0^e/p\) is nonzero. Lift it to its Teichmüller unit
\(\xi\in U_0\). The unit \(\xi p/\xi_0^e\) is congruent
to one in the residue field. Since \(p\nmid e\), Hensel's
lemma gives an \(e\)-th root of this unit congruent to one
in \(E\). Multiplying \(\xi_0\) by this root gives
\[
                              \omega^e=p\xi .
\]
This is the usual tame uniformizer construction.
[Serre, Local Fields, Chapters II and III](https://www.sas.rochester.edu/mth/sites/doug-ravenel/otherpapers/Serre-Local.pdf)

The equation is Eisenstein over \(U_0\), so \(\omega\) generates
\(E\) over \(U_0\). Its powers through \(e-1\) form an integral
basis over the integers of \(U_0\). Indeed a linear combination
of these powers has terms of pairwise distinct valuations modulo
integers, and is integral only if each coefficient is integral.
The relative discriminant up to a unit square therefore equals
the polynomial discriminant
\[
    C_e(p\xi)^{e-1},\qquad C_e=(-1)^{e(e-1)/2+e-1}e^e .
                                                               \tag{3}
\]
In fact the resultant with the derivative is
\(e^e(\prod_i\omega_i)^{e-1}\), and the product of roots is
\((-1)^{e-1}p\xi\).

In the branch of \(\alpha\), the residue \(z_0\) of
\(\alpha^b/p^a\) is a root of \(g\) in \(\mathbb F_{p^F}\).
Since \(v_p(\alpha)=a/b\) and \(e=bu\), write
\(\alpha=\omega^{au}\varepsilon\) for a unit \(\varepsilon\)
of \(E\). This yields the exact residue identity
\[
                         z_0=\overline \xi^{\,a}
                                         \overline\varepsilon^{\,b}.
                                                               \tag{4}
\]
For odd \(e\) the power \(\xi^{e-1}\) in (3) is a square.
If \(b\) is even, \(a\) is odd and (4) determines the
squareclass of \(\overline \xi\) as that of \(z_0\).
In the latter case its norm's symbol is
\(\chi_p((-1)^fg(0))^{v}\), because the norm of \(z_0\)
from \(\mathbb F_{p^f}\) is the signed constant coefficient
of \(g\). If \(e\) is even and \(b\) is odd, (4) by itself
allows either symbol.

The discriminant tower identity takes the \(e\)-th power
of the unramified discriminant and the norm of (3) from
\(U_0\). The unramified discriminant is a unit. Its
square root can be taken as an embedding determinant
in \(U_0\), so Frobenius acts on it by \((-1)^{F-1}\).
For odd \(p\) the unit's symbol is thus \((-1)^{F-1}\).
For \(p=2\) it is a square if \(F\) is odd and has class
\(5\) if \(F\) is even. Indeed \(\mathbb Q_2(\sqrt 5)\)
is the unramified quadratic extension, by reduction of
\(x^2+x-1\), and an embedding determinant changed in sign
by Frobenius generates that extension. Taking the norm of
(3) now gives (1). For \(p=2\), \(e\) is odd and the
squareclasses of \(C_e\) and \(\xi^{e-1}\) are respectively
\((-1)^{(e-1)/2}e\) and one, giving (2).

## The full tame completion test

The same residue identity determines the possible tame
completions up to isomorphism. Fix a permitted \(e=bu\) and
\(F=fv\), set \(K=\mathbb F_{p^F}\), and take \(z_0\in K\) to
be a root of \(g\). A tame completion of these indices can
be represented by a class \([\lambda]\in K^\times/(K^\times)^e\).
Over the unramified field \(U_0\) of degree \(F\) it is generated
by \(\omega^e=p\,\operatorname{Teich}(\lambda)\).
Identify these classes under Frobenius to obtain the
\(\mathbb Q_p\)-isomorphism types. To see this, an isomorphism
restricts to a power of Frobenius on \(U_0\), and the ratio
of the image of one uniformizer to the other is a unit.
Taking its \(e\)-th power and reducing gives exactly the
stated class relation. Conversely, such a class relation
lifts through Teichmüller units and gives an isomorphism
after the corresponding Frobenius action.

The branch allows precisely the classes whose images satisfy
\[
                    [\lambda]^a=[z_0]
                  \quad\text{in }K^\times/(K^\times)^b,           \tag{5}
\]
with the choice of conjugate \(z_0\) taken up to Frobenius.
Necessity is (4). For sufficiency as a local condition,
choose \(\bar\varepsilon\) with
\(z_0=\lambda^a\bar\varepsilon^b\) and form
\(\alpha_0=\omega^{au}\operatorname{Teich}(\bar\varepsilon)\).
One can add an element of sufficiently large valuation to
\(\alpha_0\) without changing its valuation or the residue
of \(\alpha_0^b/p^a\). Choose the sum outside all proper
intermediate fields, which are finitely many proper
\(\mathbb Q_p\)-linear subspaces with empty interior.
The sum then generates the completion and satisfies the
prescribed local conditions. For \(e=b,F=f\) the class is
unique up to Frobenius, since \(a\) is invertible modulo \(b\)
when \(b>1\). The case \(b=1\) is unramified.
In particular, a common residual factor of multiplicity
one fixes the completion if that component is tame.

The test retains information beyond the discriminant unit.
At \(p=7\) the Eisenstein cubics \(X^3-7\) and \(X^3-14\)
define nonisomorphic completions because \(1\) and \(2\)
lie in distinct cube classes modulo seven. Their discriminants
differ by a factor of four, so their discriminant squareclasses
are the same. When candidate global fields have been enumerated,
their completions can be compared directly against (5).

## The augmented ledger

For each candidate total field degree, enumerate local
factor degrees and ramification indices as in the cited
valuation ledger. At a tame part include both its exponent
and the unit symbol from (1) or (2). One may use both
available symbols in the ambiguous odd-prime case. For a
wild part keep the permitted valuations from the older
ledger and permit any unit symbol. The total discriminant
valuation is the sum of the local exponents and the unit
squareclass is their product. This follows by decomposing
the local trace form into the direct sum over the local
fields of the one global root field. The local integral
bases differ from a localized global integral basis by a
unit determinant, which does not change the squareclass.

Thus a global candidate \(\Delta\) from the older ledger
must satisfy the unit test on
\(\Delta/p^{v_p(\Delta)}\) at each chosen prime as well.
This applies also when \(p\) is outside the gcd of the
two array discriminants. Such a prime has zero allowed
valuation, so all its local parts must be unramified.
The same test there is the product of the signs of the
residue Frobenius cycles. This condition can be computed
using the actual shared stepped residual factors, including
the allocation of the candidate degree. The resulting
sieve retains only necessary conditions for a physical
common root.

## A strict ordinary fixture

Take the two fully shared rows \((17,14),(2,6)\) at orders
two and three. The primitive arrays in increasing
coefficient order are
\[
        (336226,31465,6930),\qquad
        (27657385,3346776,774774,222768).
\]
Their signed discriminants are \(-8330138495\) and
\(-815488706270069368607712096\), whose gcd in absolute
value is seven. Even positivity restricts a common
factor to degree two. The only possible negative
field discriminant by the older ledger is \(-7\).
Indeed the \(7\)-adic valuations of the two arrays
are \((0,1,1)\) and \((1,0,1,1)\). Their only shared
root valuation is \(-1/2\), with width two and common
stepped residual factor \(X+3\) of multiplicity one.
This permits \(e=b=2,\ F=f=1\), with discriminant
valuation one, so the older ledger permits \(-7\).

But the residue \(z_0=-3=4\bmod7\) is a square.
The tame quadratic discriminant unit in (1) is
therefore a square. The required unit \((-7)/7=-1\)
is a nonsquare modulo seven. Equation (1) excludes the
common factor. This example tests the added information
at a discriminant prime itself.

## Quadratic contact denominator growth with bounded fields

Suppose a sequence of actual collisions for arbitrary
ordinary profiles gives irreducible quadratic common
factors whose root fields have bounded discriminants.
Let \(T\) be the reciprocal trace of the two roots.
Assume \(T\) approaches an algebraic reference value
\(T_*\in(-2,2)\) with error \(O(\eta^\gamma)\), in the
notation of the
[moving trace comparison](MULTIVARIATE_MOVING_PRIME_TRACE_RESULTANT.md).
Its rational denominator is the projective denominator
\(D_H\) in degree two. Pass to a subsequence with one
fixed field discriminant \(\Delta=-D<0\).

The ratio of the two roots is a root of \(x^2-Tx+1\)
in their quadratic field. For \(T\in(-2,2)\) this ratio
is nonreal, so \(T^2-4=\Delta w^2\) for rational \(w\ne0\).
The positive coordinate \(z=|w/(T+2)|\) satisfies
\[
 T=\frac{2(1+\Delta z^2)}{1-\Delta z^2}.
\]
For \(z=a/b\) in lowest terms, the denominator of \(T\)
is \((b^2+D a^2)/g\) with \(g\mid4D\). Indeed any common
divisor of \(b^2+D a^2\) and \(2(b^2-D a^2)\) divides
\(4D a^2\), and it is coprime to \(a\).
In particular \(D_H\ge b^2/(4D)\).

The displayed transformation has a local analytic
inverse at the positive reference coordinate \(z_*\).
Thus \(|z-z_*|=O(\eta^\gamma)\), and \(z_*\) is algebraic.
If \(z_*\) is irrational, the
[theorem of Roth](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/EFB89E2873019B246F004FAA06E05A7F/S0025579300000644a.pdf/rational-approximations-to-algebraic-numbers.pdf)
gives \(|z-z_*|\ge c_\epsilon b^{-2-\epsilon}\)
for every \(\epsilon>0\), and hence
\[
                          D_H\gg \eta^{-2\gamma/(2+\epsilon)} .
\]
If \(z_*\) is rational, either the trace equals its
reference value or the elementary rational difference
bound instead gives \(D_H\gg\eta^{-2\gamma}\).
Only finitely many discriminants are under consideration,
so the estimates hold uniformly after taking the weakest
constant among the corresponding possibilities. The
irrational estimate from Roth is qualitative. It
strengthens direct rational approximation to the trace
by using the fixed field of the ratio.

The field tests constrain arbitrary ordinary moving
profiles. Uniform control of all field discriminants or
all trace denominators at the remaining contacts is still
needed for an unrestricted exclusion.
