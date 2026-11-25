# Finite-place embedding determinants in primitive power norms

Fix integers \(q,k\ge2\) and a sign \(s\in\{\pm1\}\). Consider
primitive, positive inputs \(U=Ax^q,V=By^q\le X\) with
\(T=|U+sV|=CR^k\le X\), \(T\ne0\), and \(\gcd(U,V)=1\).
Fix \(A,B,C\). Let \(K\) be the degree-\(q\) étale algebra
\(\mathbb Q[Z]/(Z^q+sA^{q-1}B)\), with \(\theta\) the class of \(Z\).
All norms, ideal groups and class groups below are products over its
field factors. The discriminant \(D_K\) is that of the product of
their rings of integers.
We may assume \(\gcd(R,qkABC)=1\). To relax the restriction
\(\gcd(R,qk)=1\), record the complete \(qk\)-parts of \(T\) and move
them into \(C\) at subpower cost, with no increase in the cap on
\(R\). Primitivity supplies \(\gcd(T,AB)=1\). The disjointness
\(\gcd(R,C)=1\) holds when the varying packet lies on original full
valuation positions, as in the application below.

Fix divisor locks \(F\mid U,\ G\mid V\) with
\(\gcd(FG,qkAB)=1\), and set
\[
 M=\prod_{p\mid FG}p^{\lceil v_p(FG)/q\rceil},\qquad
 C_0=\prod_{\substack{p\mid C\\p\nmid qkAB}}p^{v_p(C)}.
\]
If \(R\le Z\), a sufficient condition for the fixed-record count
\[
 O_{q,k,\varepsilon}\bigl(X^\varepsilon
                     |\mathrm{Cl}(K)[k]|\bigr)               \tag{VN1}
\]
is
\[
 M^{\binom q2} C_0^{\binom{q-1}2}>q!2^{q/k} Z.                 \tag{VN2}
\]
The refinement below replaces the left side by
\((FG)^{(q-1)/2}C_0^{\binom{q-1}2}\) with a constant depending on
\(q,k\), and allows arbitrary lock primes.
One may replace the right side by \(q!2^{q/k}Z/\sqrt{|D_K|}\).
For full valuation locks \(M=(FG)^{1/q}\).
The bounds include distinct target occupancy with all three letters
varying under the specified records.

## Recovery of the variable ideal packet

Set \(\alpha=Ax-\theta y\), so \(|N\alpha|=A^{q-1}T\). At every
good prime \(p\nmid qkAB\) dividing \(T\), the selected root of the
generator polynomial is \(Ax/y\bmod p\), a root distinct from
the others by the discriminant. All \(x,y,A,\theta\) used here are
local units. Exactly one degree-one prime of \(K\) can divide
\((\alpha)\) over \(p\), and its exponent is \(v_p(T)\).
Fix the choice of this prime at each good \(p\mid C\). At the
remaining fixed primes dividing \(qkABC\), the norm exponents are
known, and the allocations of them among primes of \(K\) cost
\(X^{o_X(1)}\), uniformly: the choices are bounded by a fixed
power of a divisor function of \(qkA^{q-1}C\). The good-arm
assignments cost at most \(q^{\omega(C)}=X^{o_X(1)}\).
This gives a fixed nonzero integral ideal \(\mathfrak a\) with
\[
            (\alpha)=\mathfrak a\mathfrak b^k,\qquad
      N\mathfrak b=R,\quad \gcd(R,qkAB C)=1 .
\]
Each occurring class of \(\mathfrak b\) lies in a coset of
\(\mathrm{Cl}(K)[k]\). In each occurring class choose an actual
reference \((\alpha_*,\mathfrak b_*)\). Then
\[
   \mathfrak b=(\gamma)\mathfrak b_*,\qquad
   \alpha/\alpha_*=u\gamma^k,\qquad
   \gamma\in\mathfrak b_*^{-1}.
\]
Record one of \(O_{q,k}(1)\) representatives for \(u\) modulo
\(k\)-th powers of units, using the fixed rank and torsion bounds
for the unit group of the factors of the degree-\(q\) algebra.

## Archimedean and nonarchimedean determinant bounds

Write \(\theta_i=\sigma_i(\theta)\) for all embeddings. There are
\(O_q((\log X)^q)\) possible dyadic boxes on the \(q\) magnitudes
\(|\sigma_i\alpha|\): each is at most \(2X\) and at least
\((2X)^{-(q-1)}\), the latter by the nonzero integral norm.
Work in one box. For any \(q\) candidates \(\gamma_j\), let
\[
 D=\det_{i,j}(\sigma_i\gamma_j)
       = \pm \frac{\sqrt{D_K}}{N\mathfrak b_*}\,I,\qquad I\in\mathbb Z ,
\]
with a fixed choice of square root. This is valid also in the
product algebra by taking the concatenated integral bases of the
fractional ideals. It is the determinant of
an integral-coordinate matrix relative to an integral basis of
the fractional ideal \(\mathfrak b_*^{-1}\), using its usual
covolume formula. For each \(i\), the values
\(|\sigma_i\gamma_j|\) differ by at most \(2^{1/k}\) as \(j\)
varies, since the unit and reference element are fixed.
As \(|N\gamma_j|=R_j/N\mathfrak b_*\), the Leibniz formula gives
\[
                  |D|\le q!2^{q/k}Z/N\mathfrak b_* .           \tag{VN3}
\]

At a prime \(p\mid FG C_0\), choose an embedding of the normal
closure of \(K\) into the algebraic \(p\)-adic closure. The roots
\(\theta_i\) have distinct reductions, since \(p\nmid qAB\).
We work in the completion of a maximal unramified extension.
All quantities called units below are indeed units there.
All \(k\)-th roots of units can be taken in this extension since
\(p\nmid k\), by first choosing the residue root, then Hensel.
Branches of unit roots in \(\mu_k\) cost at most \(k^q\) choices
per prime, hence \(X^{o_X(1)}\) across all the present primes.
All local valuations of \(\sigma_i\gamma_j\) are zero, because
the ideals \(\mathfrak b,\mathfrak b_*\) have norms coprime to
\(p\): for \(p\mid FG\) use primitivity, and for \(p\mid C_0\)
use \((R,C)=1\).

For \(p\mid F\), put \(t_j=x_j/y_j\in p^m\mathbb Z_p\), with
\(m=\lceil v_p(F)/q\rceil\). The equality for \(\alpha\) gives
\[
  \sigma_i(\gamma_j)^k
     = y_j\frac{-\theta_i}{\sigma_i(\alpha_*u)}
                \bigl(1-(A/\theta_i)t_j\bigr).
\]
Every row prefactor is a unit, as is each \(y_j\). Choose a
\(k\)-th root of \(y_j\) as column factor, and fix root branches
of the row prefactors. After recording the \(k\)-th-root-of-unity
choices for each point, remove the unit row and column factors.
Each column is then a vector of fixed integral power series
in \(t_j\). Integrality of the binomial series follows from
\(1/k\in\mathbb Z_p\), by approximating it with nonnegative
integers in each binomial coefficient. The determinant has
valuation at least \(m\binom q2\), since expanding in the
Taylor-coordinate vectors gives zero for repeated powers and
the least sum of \(q\) distinct nonnegative powers is \(\binom q2\).
This argument can use finite truncations to the required precision.
For \(p\mid G\) use \(y_j/x_j\) and \(x_j\) as the unit scale;
the row constants are \(A/\sigma_i(\alpha_*u)\).

For \(p\mid C_0\), put \(e=v_p(C)\).
Exactly one row \(i_0\), fixed by the degree-one prime record,
satisfies \(v_p(Ax_j-\theta_{i_0}y_j)=e\) for every \(j\).
Thus \(t_j=x_j/y_j\) is congruent to the fixed reference
\(t_*\bmod p^e\) and \(y_j\) is a unit. All the other rows have
\(A t_*-\theta_i\) a unit, since the reductions are distinct.
Remove the column roots of \(y_j\) and make the same branch
records and Taylor expansion in \(t_j-t_*\) for these \(q-1\)
rows. The exceptional row still consists of units after removing
the column factors. Expand along that row. Each minor of the
other \(q-1\) rows has valuation at least \(e\binom{q-1}2\),
so the full determinant does as well. The \(F,G,C_0\) supports
are disjoint by primitivity.

At each of these primes, both \(D_K\) and \(N\mathfrak b_*\)
are units. Consequently the integer \(I\) is divisible by
\(M^{\binom q2} C_0^{\binom{q-1}2}\). By (VN2) and (VN3) every
such \(q\)-tuple determinant is zero.

## A bounded count from the pullback line

The projective points \([\sigma_i\gamma]_{i=1}^q\) lie in the
pullback of a line under the coordinatewise \(k\)-th power map
of \(\mathbb P^{q-1}\). On the line the coordinates are
the linear forms \((A x-\theta_i y)/\sigma_i(\alpha_*u)\),
with distinct zeroes on its \([x:y]\)-parameter. This pullback
is a geometrically irreducible curve of degree at most
\(k^{q-2}\) and spans the ambient space. In fact the \((q-1)\)
ratios to the last linear form are independent modulo \(k\)-th
powers in the rational function field over an algebraic
closure: take valuations at their distinct simple zeroes.
Kummer theory gives full \((\mathbb Z/k)^{q-1}\) degree and
distinct characters on their \(k\)-th-root ratios, proving
irreducibility and linear independence of the \(q\) projective
coordinates. The power-map pullback is cut out by \(q-2\)
degree-\(k\) equations. The finite map and the dimension theorem
for the common zeroes show that every component dominates the
line, so generic irreducibility accounts for the whole curve.
Bézout supplies the degree bound.

Vanishing of all \(q\)-tuple determinants puts the sample in a
hyperplane, which contains at most \(k^{q-2}\) points of the
curve. Different primitive positive pairs give different
projective points, because a point determines the line parameter
\([x:y]\) by taking \(k\)-th powers, and \(\gcd(x,y)=1\).
This proves (VN1) after summing local, dyadic, unit and class
records.

## Uniform class cost and a strictly saving sector

If \(q=p^a,k=p^b\), the existing cyclic \(p\)-tower proof in
ABC_PRIMITIVE_NORM_RAY_PACKING.md gives
\[
          |\mathrm{Cl}(K)[k]|
                \le \exp(O_{q,k}(1+\omega(AB)))=X^{o_X(1)} .
\]
Let \(K_j\) run over the field factors. The splitting field
\(E=\mathbb Q(\mu_q,\theta_1)\) of the polynomial is a cyclic
extension of \(p\)-power degree over the fixed field
\(\mathbb Q(\mu_q)\), because \(q=p^a\). Its cyclic degree-\(p\)
steps ramify only over \(pAB\). In such a step the dimension of
class \(p\)-torsion is at most \(p\) times its invariant dimension.
For invariant classes the norm obstruction is a quotient of
the group of base units modulo norms, of bounded \(p\)-rank.
Hilbert 90 gives invariant ideal representatives after removing
this obstruction. The invariant ideals are products of
extended ideals and powers of ramified prime ideals. Iterating
the resulting rank bound from the fixed base yields
\(\log|\mathrm{Cl}(E)[p]|=O_q(1+\omega(AB))\).
[Prime-power class-group descent](ABC_PRIMITIVE_NORM_RAY_PACKING.md#prime-power-towers-remove-the-class-number-cost)

Each \(K_j\) is embedded in \(E\). The kernel of extension of
ideal classes from \(K_j\) to \(E\) injects into
\(H^1(\operatorname{Gal}(E/K_j),\mathcal O_E^\times)\).
Indeed an ideal becoming principal gives the cocycle of a
principal generator divided by its conjugates, up to inverting
the cocycle convention. A coboundary changes that generator
by a unit into a \(K_j\)-element, which generates the original
ideal by injectivity of extension on ideals. This cohomology
group is killed by the group order and needs a bounded number
of generators: the cocycles are a subgroup of a bounded
product of finitely generated unit groups of bounded rank and
torsion. The kernel therefore has bounded size. This controls
\(\mathrm{Cl}(K_j)[p]\) by the tower bound. Pass to \(p^b\)
by filtration and take the product over \(j\).

As a concrete uniform \(abc\) sector take \(q=2,k=4,T=c=U+V\),
with squarefree kernels \(A,B\) satisfying
\(\operatorname{rad}(AB)\le X^{\lambda/300}\). Suppose the actual
target splits into coprime \(C R^4\) with
\[
  \operatorname{rad}(C)\le X^{(1/4+1/300)\lambda},
             \qquad R\le X^{7/80},
\]
and suppose there are actual full locks on \(U,V\) off \(2AB\)
with \(\operatorname{rad}(FG)\le X^{\lambda/15}\) and
\(FG\ge X^{9/50}\). The binomial is irreducible since
\(Z^2+AB\) has no rational root and \(A,B>0\).
A power of \(2\) in \(R\) can be moved into \(C\) with a
subpower overhead, since the support \(\{2\}\) is fixed.
The determinant threshold follows from
\(M=\sqrt{FG}\ge X^{9/100}>X^{7/80}\). Enumerate the actual
\(A,B,C,F,G\) at their joint radical costs using the fixed-support
count and divisor-many assignments. This gives
\[
        \ll_{\lambda,\varepsilon} X^{97\lambda/300+\varepsilon}
              = X^{(1/3-1/100)\lambda+\varepsilon}.
\]
The estimate holds at every fixed budget, with all three height
and radical caps retained. It does not require \(U,V\) or \(T\)
to be short.

For a model of the covered profile, give \(U,V\) radical mass
\(1/3\) each on dense positions \(e\equiv2\pmod4\) close to scaled
degree \(3\). Give the target defect \(C\) mass \(1/4\) on odd
positions close to scaled degree \(8/3\), and the target fourth-power
packet mass \(1/12\) at degrees near \(4\). All three numerical
heights are one in the centered model. With pair eligible degree
at least \(2.9\) and individual normalized prime sizes at most
\(1/240\), a greedy full pair selection in \([1/16,1/15]\)
supplies \(FG\ge X^{9/50}\). The new hypothesis also survives small
perturbations of the profile; a fixed bounded packet of positions
need not carry positive mass. The Euclidean \((2,4,4)\) direct
grid cannot assign a fourth-power slot to either whole \(U\) or
whole \(V\), and the fixed \((2,2,4)\) signature is spherical.
For example \((x+iy)=(u+iv)^4\) with coprime \(u,v\) of opposite
parity gives \(x^2+y^2=(u^2+v^2)^4\) and \(\gcd(x,y)=1\).
There are \(\gg X^{1/4-\varepsilon}\) distinct targets \(c\le X\)
so obtained, by taking \(u,v\) in a fixed positive-angle box of
size \(X^{1/8}\) (in particular avoid the finitely many zero rays)
and using the divisor bound on representations by two squares.
Thus counting the fixed coefficient spherical equation alone
cannot replace the new paid locks.

## Refinement: numerical weight on coefficient supports

The coprimality restriction of \(FG\) to \(ABqk\) can be
removed. At \(p\mid F\) put \(z_j=A x_j/y_j\). If \(p\mid A\),
factor \(-\theta_i\) out of \(z_j-\theta_i\); the small quantity
\(z_j/\theta_i\) has normalized valuation
\(v_p(Ax_j^q)/q\ge v_p(F)/q\), since all \(y_j\) are units and
\(v_p(\theta_i)=(q-1)v_p(A)/q\). Consequently the row constants
\(-\theta_i/\sigma_i(\alpha_*u)\) remain units. The same
factorization works at \(p\nmid A\), where \(p\mid F\) implies
\(p\nmid B\) and the \(\theta_i\) are units. The case \(p\mid G\)
uses \(y_j/x_j\), with the dominant \(A x_j\) term and the same
small relative valuation \(v_p(G)/q\). Work in the completion
of an algebraic \(p\)-adic closure, allowing fractional
normalized valuations and ramified roots. All the \(\gamma_j\)
are still units there because their ideal supports avoid \(FG\).
At \(p\nmid k\), the integral binomial argument now gives
\[
       v_p(D)\ge \frac{q-1}{2}v_p(FG)
\]
without requiring distinct reductions of the \(\theta_i\).
Those reductions were needed only in the good-\(C\) step,
which is unchanged.

For \(p\mid k\) set \(b_{p,k}=v_p(k)+1/(p-1)\). The \(n\)-th
coefficient of \((1+z)^{1/k}\) has valuation at least
\(-n b_{p,k}\), as follows directly from the product formula
for that coefficient and \(v_p(n!)\le n/(p-1)\).
If \(v_p(FG)/q>b_{p,k}\), the series converges, and the proof
loses at most \(\binom q2 b_{p,k}\) in the determinant valuation.
If \(v_p(FG)/q\le b_{p,k}\), the trivial unit-entry determinant
bound suffices for the same relaxed inequality. The total
multiplicative loss in the determinant lower bound is thus the
fixed constant \(\prod_{p\mid k}p^{\binom q2 b_{p,k}}\).

For any nonzero \(D=\sqrt{D_K}\, I/N\mathfrak b_*\) the real
absolute value satisfies
\[
  |D|\ge\frac{1}{N\mathfrak b_*}\prod_{p\mid FG C_0}p^{v_p(D)},
  \qquad
  v_p(D)=v_p(I)+v_p(D_K)/2
\]
on the indicated primes, since \(p\nmid N\mathfrak b_*\).
This is just the prime-factor decomposition of
\(|I|\sqrt{|D_K|}\), with all other factors discarded. In
particular one need not subtract a discriminant loss at a
ramified lock. Thus the sufficient condition becomes
\[
       (FG)^{(q-1)/2} C_0^{\binom{q-1}2}>c_{q,k}Z .
\]
This form includes all full factors supported on \(AB\), at
subpower additional cost after \(A,B\) were charged. It also
includes arbitrarily large powers of the fixed primes dividing
\(k\).

For instance in the degree bands with \(q=2,k=4\) where
the three rates are close to \(3\) and the target packet has
mass \(s>0\), the bound on \(R\) is close to \(X^{3s/4}\).
If the joint \(U,V\) parity kernel has normalized mass \(f\),
the full locks on its support are already paid and contribute
weight close to \(3f/2\) to \(\sqrt{FG}\). Suppose the remaining
eligible primes on the two summands supply enough radical
mass for the extra selection, with individual normalized masses
tending to zero, or bounded by a small constant compared to the
saving margin. The extra selector then costs at most
\((s/2-f)_+\) up to degree-band and rounding errors. Its
hypotheses concern the actual prime powers, including the
coprime split \(T=C R^4\). The total record exponent is
\[
                  1/3-s+\max(f,s/2)+o(1).
\]
This gives a strict saving whenever \(s\) and \(s-f\) are
bounded below and the rates sufficiently close to \(3\) by
comparison with those margins. This allows a non-negligible
kernel cost even when the two complete letters and the target
all have full numerical height. In a centered equal-rate
profile the one-third charged parity selector on all three
letters retains no rate advantage, whereas the fourth-power
target packet gives the gain here.

For the local calculation with fractional valuations at \(p\mid F\),
one can explicitly use \(z_j=(A/\theta_1)x_j/y_j\). In row \(i\),
the multiplier of \(z_j\) is \(\theta_1/\theta_i\), a unit.
The analogous \(G\) chart uses \(z_j=(\theta_1/A)y_j/x_j\).
This yields the common small parameter needed for the Taylor-column
argument.

## Other binary norms

For a different separable degree-\(q\) binary form \(P\), the main
good-prime proof also applies to its norm algebra. Assume the form
coefficients and input magnitudes are at most polynomial in \(X\).
Using the suitably rescaled root of \(P(t,1)\), write
\(|N\alpha|=a^{q-1}|P(x,y)|\), with \(a\) its nonzero leading
coefficient. At the primes of the lock require a common projective
residue \([x:y]\) to precision \(m_p\), with \(P(x,y)\) a unit
and the norm form having good reduction. Choose the unit
homogeneous coordinate and Taylor-expand about that common
residue. These primes contribute \(\prod p^{m_p\binom q2}\).
Good fixed primes of the coefficient \(C\) contribute the same
\(C_0^{\binom{q-1}2}\) as before. The class cost is the product of
\(k\)-torsion sizes of the norm fields, with no claim that it is
subpower uniformly in arbitrary growing forms. This extends the
mechanism to exact power packets in primitive binary form values.

## Selecting the complete divisors

The general selector uses the same accounting. Let \(f\) be the
joint coefficient radical cap for \(A,B\), and \(g\) the radical
cap for \(C\), in units of \(L\). The full parts of \(U,V\) on
the support of \(AB\) can be fixed with a subpower number of
choices after the coefficient record is paid. Write \(w_0\) for
their numerical weight, \(w_C\) for the numerical weight of
\(C_0\), and \(W=\log_X(T/C)\) in the prescribed box. If an
extra pair selector of normalized radical cost \(t\) has
numerical weight \(w\) and
\[
              \frac{q-1}{2}(w_0+w)
                 +\binom{q-1}{2}w_C>\frac Wk
\]
by a fixed margin, the exponent at matching prime-power orders
is at most \((f+g+t)L+\epsilon\) with the usual box allowances.
Coefficients supported on a fixed union of primes are counted
by the fixed-support estimate. If a prime of \(C\) is also
in a proposed \(R\), freeze its full part of \(R^k\) on that
support and include it in \(C\). This incurs subpower cost in
the height limit, decreases the required bound on \(R\), and
can only improve the good-\(C\) term.

Darmon and Granville use the ideal-class and Kummer reduction
for weighted power equations and describe their spherical covers.
Heath-Brown's local determinant gives the general predecessor for
the divisibility step. The result here is a finite-place counterpart
of the full-embedding real determinant for binomial norms. The primitive ideal and
class records are shared, while the determinant here obtains
precision from the paid complete divisors and the fixed
target coefficient. The closest unrestricted fixed-coefficient
spherical count has a positive height exponent
\(1/q+1/q+1/k-1\) whenever this is positive and a primitive
point exists. The congruences above are additional conditions
on that curve. In comparison with the high exact-order
valuation grids, the spherical example uses the fixed orders
\(2,2,4\) and only a packet of the target.
[Darmon–Granville, §§2.1 and 7.1](https://www.math.mcgill.ca/darmon/pub/Articles/Research/12.Granville/pub12.pdf)
· [Heath-Brown, §3](https://arxiv.org/abs/math/0405392)
· [Full embedding real determinant](ABC_PRIMITIVE_NORM_RAY_PACKING.md#full-embedding-determinant-for-binomial-rays)
· [Arango-Piñeros, Theorem 1.3, version 2](https://arxiv.org/html/2508.13093v2)
· [Two high valuation grids](ABC_TWO_VALUATION_GRIDS.md)

These estimates leave diffuse prime valuation profiles open.
For example, on formal profiles with equal mass spread over
increasingly many prime orders on each letter, any common
\(q\ge2\) recovers at most one valuation position per letter.
Then paying the two binomial coefficients can already cost
\(2/3-o(1)\) of the total radical budget. This formal
obstruction concerns the use of a common exact order.
