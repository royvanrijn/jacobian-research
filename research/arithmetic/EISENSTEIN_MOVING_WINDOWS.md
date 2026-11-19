# Pointwise Eisenstein counts in moving summand windows

Let \(A\) consist of \(1\) and the positive integers whose prime factors
are all \(1\bmod 3\). In an ordered pair write \(s_1,s_2\) for the
**actual summands**. For \(n\equiv2\bmod6\) require \(s_1,s_2\in A\).
For \(n\equiv4\bmod6\) fix either \((s_1,s_2)\in A\times(3A)\) or
\((s_1,s_2)\in(3A)\times A\). Let \(r(n;x,H)\) count the distinct
ordered value pairs in the specified family with
\[
 s_1+s_2=n,\qquad |s_1-x|\le H/2 .
\]
In particular, when \(s_1=3a\) it is \(3a\) that belongs to the
specified interval. Put
\[
 J(n)=\prod_{\substack{p\mid n,\ p\ge5\\p\equiv2\bmod3}}
                 \frac{p-1}{p-2}.
\]

**Theorem EisLocationWindows.** For every fixed \(0<\varepsilon<1\)
there are positive \(c_\varepsilon,C_\varepsilon,N_\varepsilon\) such
that for every admissible integer \(n\ge N_\varepsilon\), uniformly in
real \(x,H\) satisfying
\[
 0<x<n,\qquad M=\min(x,n-x),\qquad 1\le H\le M,\qquad
 H^3\ge n^{1+\varepsilon}M ,
 \tag{1}
\]
one has
\[
 c_\varepsilon\frac{H J(n)}{\log n}\ \le\
 r(n;x,H)\ \le\ C_\varepsilon\frac{H J(n)}{\log n}.
 \tag{2}
\]
Each prescribed mixed orientation satisfies (2) separately.

At \(x=n/2\) condition (1) covers every fixed power \(H=n^\theta\)
with \(2/3<\theta<1\). When the prescribed summand has size comparable
to \(H\), it covers every fixed \(1/2<\theta<1\). The transition comes
from the cost of localizing the angular test on a quaternion shell. That
cost is \(O(\sqrt{M/H})\), and remains bounded for \(M\asymp H\).

## Projection for directional tests

Write a quaternion in the discriminant-\(3\) order as \(\alpha=z+yj\)
with \(z,y\in\mathbf Z[\omega]\), \(\omega^2+\omega+1=0\). Its reduced
norm is \(m=N(z)+N(y)\). Write \(\mathcal Q_m\) for this shell,
\(u_\alpha=\alpha/\sqrt m\), and \(U=\mathcal O^\times\) for its group
of twelve norm-one units. Let \(G\) be the product of the finite
projective reduction groups used by the local tests. Always include the
fixed factor
\[
 G_3=(\mathcal O/3\mathcal O)^\times/\mathbf F_3^\times.
\]
All targets at which this reduction is used are prime to \(3\).
Only the two units \(\pm1\) reduce to scalars at this factor, so the
kernel of \(U\to\Gamma\), where \(\Gamma\) is the image in \(G\), is
exactly \(\{\pm1\}\). This includes the case with no variable level.
[Ramified reduction and its six cosets](../../results/eisenstein-restricted-prime-sums/sections/mixed-local.tex).

**Lemma (directional projection).** Let \(\rho=\operatorname{Sym}^{2k}\)
of \(\mathrm{SU}(2)\), \(k\ge0\). This representation is trivial on the central signs and
therefore defines a representation of \(\Gamma\) by lifting its elements
to \(U\). Let \(f:G\to\mathbf C\) and let \(v,w\) be unit vectors of
\(\rho\). Define
\[
 F(g)=\frac1{|\Gamma|}\sum_{\gamma\in\Gamma}
                 \rho(\gamma)v\, f(\gamma^{-1}g).
\]
Then \(F(\gamma g)=\rho(\gamma)F(g)\), and
\[
 \|F\|_{\Gamma\backslash G}\le
        \frac{\|f\|_{\ell^2(G)}}{\sqrt{|\Gamma|}} .
\]
Here the equivariant norm sums squared vector norms once per left coset.
For the all-ideal Hecke operator \(T_m\), at \(m\) coprime to the finite
level, the exact identity is
\[
 |U|\langle w,(T_mF)(1)\rangle
  =\sum_{\alpha\in\mathcal Q_m} f(\alpha_{\rm fin})
                  \langle w,\rho(u_\alpha)^{-1}v\rangle .
 \tag{3}
\]
The inner product is linear in its second argument and \(\alpha_{\rm fin}\)
is the projective reduction.

Reindexing proves equivariance. On vector-valued functions on \(G\),
the average defining \(F\) is the orthogonal projection of the function
\(g\mapsto f(g)v\) under the unitary action by
\(F_0(g)\mapsto\rho(\gamma)F_0(\gamma^{-1}g)\). For an equivariant
function its squared full norm is \(|\Gamma|\) times its squared coset
norm. For (3), use
\[
 (T_mF)(1)=\frac1{|U|}\sum_{\alpha\in\mathcal Q_m}
                        \rho(u_\alpha)^{-1}F(\alpha_{\rm fin}).
\]
For each \(\gamma\) choose a unit lift \(\gamma_0\) and change variables
\(\alpha=\gamma_0\alpha'\). The factor \(\rho(\gamma_0)^{-1}\) cancels
\(\rho(\gamma)\), \(f\) is evaluated at \(\alpha'_{\rm fin}\), and
left multiplication permutes the entire shell. The residual sum is
independent of \(\gamma\). This proves (3) for arbitrary \(f,v,w\).
[All-ideal normalization, EisPower (5)](EISENSTEIN_POWER_SHRINKING_WINDOWS.md#uniform-archimedean-weight-and-finite-level).

For \(k\ge1\), the whole-space norm bound is
\(\|T_m\|\le\tau(m)\sqrt m\), with \(\tau\) the divisor function.
The cuspidality argument uses that \(\operatorname{Sym}^{2k}\) is an
irreducible nontrivial archimedean representation, a property of all
\(k\ge1\). Evaluation on a coset has norm at most one. For \(k=0\) the
scalar cuspidal complement obeys the same bound and the norm-character
contribution gives \(12\sigma(m)\) times the uniform mean of the finite
test conditioned on the joint target norm squareclasses. Averaging
\(f(\gamma^{-1}g)\) leaves this conditional mean unchanged, since every
lift of \(\gamma\) has norm one. A product of local tests consequently
retains the product of its local conditional means.
[Positive-weight cuspidality and spectral bound](EISENSTEIN_POWER_SHRINKING_WINDOWS.md#cuspidality-at-positive-archimedean-weight),
[composite-level scalar projection](EISENSTEIN_BOUNDED_PRIME_FACTORS.md#the-composite-level-estimate).

For the norm coordinate \(t=(N(z)-N(y))/m\), choose \(v=w\) proportional
to \(X^kY^k\). Its matrix coefficient in (3) is the Legendre polynomial
\(P_k(t)\), normalized by \(P_k(1)=1\), by the binomial identity in
EisPower (4). For odd \(k\) this vector changes sign under the unit \(j\).
Equation (3) still estimates its moment against a prescribed finite
orientation. [Matrix coefficient formula](EISENSTEIN_POWER_SHRINKING_WINDOWS.md#uniform-archimedean-weight-and-finite-level).

In the pure residue the test at \(3\) is identically one. In the mixed
residue the target norm is \(1\bmod3\). Among the nine residues of
\(\mathbf Z[\omega]/3\) there are six units, two nonzero nonunits and zero.
A nonzero nonunit has norm of \(3\)-adic valuation exactly one for every
lift. A zero residue has norm divisible by \(9\). Thus among the \(36\)
raw unit-quaternion residues in the specified norm class, exactly
\(6\cdot2=12\) have the specified component a unit and the other a
nonzero nonunit. The conditional mean is \(1/3\) for either orientation.
The events are invariant under scalar projectivization. This imposes the
orientation directly even when reflecting the window would move its center.

## A uniform estimate for the angular test

Fix once for all a smooth \(0\le\phi\le1\) with compact support in
\((-1/2,1/2)\) and \(I_\phi=\int_{\mathbf R}\phi(u)\,du>0\). Put
\[
 h=H/n,\qquad b=M/n,\qquad \lambda=x/n,\qquad
 W(t)=\phi\!\left(\frac{(1+t)/2-\lambda}{h}\right).
\]
The constant Legendre coefficient of \(W\) is \(I_\phi h\), by the
change of variables \(s=(1+t)/2\). On its support
\(|s-\lambda|<h/2\), so \(b/2\le\min(s,1-s)\le3b/2\) and
\(1-t^2=4s(1-s)\asymp b\). The same bounds hold on the support of
its derivatives. The supporting interval has length at most \(2h\).

For \(k\ge1\) we use the uniform estimate
\[
 |P_k(t)|\ll k^{-1/2}(1-t^2)^{-1/4},\qquad -1<t<1.
 \tag{4}
\]
Here is a proof including the behavior as \(t\) approaches either pole.
For \(u\in(0,\pi)\) put
\[
 Y(u)=\sqrt{\sin u}\,P_k(\cos u),\qquad
 Q(u)=(k+1/2)^2+\frac1{4\sin^2u}.
\]
The Legendre equation gives \(Y''+QY=0\). On \((0,\pi/2]\),
\[
 \left(Y^2+\frac{(Y')^2}{Q}\right)'
       =-\frac{(Y')^2 Q'}{Q^2}\ge0 .
\]
At \(\pi/2\), the even central values have magnitude
\(\binom{2j}{j}4^{-j}\le(j+1)^{-1/2}\). Indeed the ratio of successive
values is \((2j+1)/(2j+2)\le\sqrt{(j+1)/(j+2)}\), and the initial
value is one. Odd central values vanish, and
\(P'_k(0)=kP_{k-1}(0)\). These identities follow from the Legendre
generating function and derivative recurrence. Since \(Q(\pi/2)\ge k^2\),
the displayed energy there is \(O(1/k)\). By monotonicity it bounds
\(Y^2\) on the entire half interval. Reflection proves (4).

**Lemma (translated bump).** With
\[
 c_k=\frac{2k+1}{2}\int_{-1}^1 W(t)P_k(t)\,dt ,
\]
the translation and width above satisfy
\[
 \sum_{k\ge1}|c_k|\ll_\phi\sqrt{b/h}=\sqrt{M/H}.
 \tag{5}
\]
The constants are uniform over \(0<h\le b\le1/2\). To prove this let
\(\mathcal L=-\frac d{dt}((1-t^2)\frac d{dt})\). In the variable
\(v=(t-(2\lambda-1))/(2h)\), this operator is \(b/h^2\) times
\[
 -\frac{1-t^2}{4b}\frac{d^2}{dv^2}
       +\frac{th}{b}\frac d{dv}.
\]
On the support in question the coefficients, and all their derivatives
in \(v\) of any fixed order, are bounded by absolute constants. For
example the first two derivatives of the leading coefficient cost
\(O(h/b)\) and \(O(h^2/b)\). Thus
\(\|\mathcal L^j W\|_\infty\ll_{\phi,j}(b/h^2)^j\). Integration
by parts \(j\) times against \(P_k\), whose eigenvalue is \(k(k+1)\),
introduces no boundary terms. Using (4) and support length \(2h\)
both before and after these integrations gives
\[
 |c_k|\ll_{\phi,j}
 h\sqrt{k}\,b^{-1/4}
 \min\!\left(1,\left(\frac{b}{h^2k^2}\right)^j\right).
\]
Split at \(K=\sqrt b/h\ge1\). Summing the first bound through \(K\)
and the second with \(j=2\) above \(K\) costs
\(O_\phi(h b^{-1/4}K^{3/2})\), which is (5). Since \(|P_k|\le1\)
on \([-1,1]\), the series converges absolutely uniformly. Its Legendre
coefficients and polynomial density identify its sum with \(W\).
Consequently (3) applies termwise, including all odd degrees.

## Applying the cofactor sieve

Use the squarefree cofactors \(d,u,v\) of
[EisSemigroupScale (2)](EISENSTEIN_POINTWISE_SEMIGROUP_SCALE.md#the-coordinate-count-and-the-cofactors).
Their primes are split and at most \(w_0=n^{\lambda_0}\). They are
pairwise coprime, \(d\mid n\), \((uv,n)=1\), and \(duv\le n^\kappa\).
On the shell \(m=n/d\) use the same \(\lambda=x/n,\ h=H/n\).
Its component norm \(A_1=N(z)\) satisfies \(W(t)=\phi((dA_1-x)/H)\).
Thus after multiplication by \(d\), a positive weight implies membership
in the requested actual interval. The other residual norm is \(A_2=N(y)\).

Keep odd component norms. At \(p\mid u\) prescribe \(v_p(A_1)=1\).
At \(p\mid v\) prescribe \(v_p(A_2)=1\). Include the prescribed test
at \(3\) above and divide the coordinate weights by \(2^{\omega(uv)}\).
Sieve all other primes \(5\le p\le z_0=n^\eta\), including the primes
of \(d\), from \(A_1A_2\). Let \(\ell\) be a squarefree product of
primes outside \(6uv\). Before this last sieve, impose \(p\mid A_1A_2\)
at every \(p\mid\ell\). The resulting weighted count is
\[
 X_m\prod_{p\mid uv}(t_p/2)\,g_m(\ell)
 + O_\zeta\!\left((uv)^3\ell^{3/2}n^{1/2+\zeta}
                                      \sqrt{M/H}\right)\qquad(\zeta>0).
 \tag{6}
\]
Here \(g_m\) is the multiplicative bad-prime density of
EisSemigroupScale (5)--(6), \(t_p=2(p-1)^2/(p^2(p+1))\) is its exact
split valuation-one probability, and
\[
 X_m=I_\phi h R_o(m)\quad\hbox{in the pure case},\qquad
 X_m=I_\phi h R_o(m)/3\quad\hbox{in either prescribed mixed orientation}.
\]
The odd-shell cardinality is \(R_o(m)=12\{\sigma(m)-\sigma(m/4)\}\), with
a term of nonintegral argument taken as zero.

To verify (6), first take the Hecke index coprime to all finite levels.
The raw finite test is the product of the indicated directional factors,
each bounded by one. Equation (3) estimates this exact product. At
\(p\mid uv\) the prime-square reduction group has order
\(p^3|\mathrm{PGL}_2(\mathbf F_p)|<p^6\). The group size at the
squarefree coprime part of \(\ell\) is at most \(\ell^3\) before the
fixed factor at \(3\). The Hilbert norm bound therefore costs at most
a constant times \((uv)^3\ell^{3/2}\), while (5) gives the angular
cost in (6). At degree zero the character projection gives the scalar
factorization in (6) and its cusp error is absorbed there since
\(\sqrt{M/H}\ge1\).
[Prime-square local computation and level normalization](EISENSTEIN_POINTWISE_SEMIGROUP_SCALE.md#the-joint-level-estimate).

For singular \(p\mid(\ell,m)\) the same local elimination applies.
At inert \(p\) a bad pair is empty or obtained from a shell at \(m/p^2\)
by multiplying both coordinates by \(p\). At split \(p\),
EisBoundedFactors (6) resolves the bad count into terms at \(m/p^j\),
\(1\le j\le v_p(m)\), with the remaining modulo-\(p\) test only at
the \(p\)-unit target. It scales the two component norms by the same
factor, preserving this translated \(W\). Divisibility tests at other
primes, including the specified orientation at \(3\), are preserved.
Any remaining test of divisibility by one of the primes above a
different split prime is compatible by unique factorization of the
Eisenstein coordinates. All \(p\mid uv\) stay coprime to the reduced
targets. Extracted split prime powers and inert prime squares are
\(1\bmod3\), so the scalar mixed mean persists. The combined absolute
coefficient loss is at most \(\tau(m)^4\). Odd subtraction removes
the scaled shell at \(m/4\) with identical angular test and preserved
local valuations. These operations prove (6) for arbitrary target
valuations, absorbing the divisor powers into \(n^\zeta\).
[Singular-prime recurrence and bounds](EISENSTEIN_BOUNDED_PRIME_FACTORS.md#singular-primes).

Choose fixed \(\kappa,\xi>0\) satisfying
\[
                 4\kappa+\tfrac52\xi<\varepsilon/2 .
 \tag{7}
\]
Use the lower block polynomial at divisor level \(D=n^\xi\), with
\(\eta=\xi/S\) for its fixed block exponent \(S\). The upper bound for
the absolute coefficients is \(O(\log\log n)\). Its main prediction for
each triple is at least
\(\tfrac12 X_m V_m(z_0)\prod_{p\mid uv}a_p\), where
\(V_m(z_0)=\prod_{5\le p\le z_0}(1-g_m(p))\) and
\(a_p=(p-1)/(p(p-2))\). Summing the errors in (6) over \(\ell\le D\)
costs \(D^{5/2}\). The number of triples with \(duv\le n^\kappa\)
is \(n^{\kappa+o(1)}\), and \((uv)^3\le n^{3\kappa}\).
Consequently the aggregate error is
\[
 \ll_\zeta n^{1/2+4\kappa+(5/2)\xi+\zeta}\sqrt{M/H},
\]
where \(\zeta>0\) may be chosen arbitrarily small. Condition (1) says
\(\sqrt{nM/H}\le H n^{-\varepsilon/2}\). By (7) the error is
\(o(H/\log n)\), uniformly in \(x,H\).
[Block weights and their level](EISENSTEIN_POINTWISE_SEMIGROUP_SCALE.md#sifting-outside-the-prescribed-cofactors).

The scalar calculation uses \(h m=H/d\) for the current moving window.
The unit-target sieve product is comparable to \(1/(\log z_0)^2\).
At split \(p\le w_0\) with \(p\nmid n\) the three cofactor options
supply \(1+2a_p\). At split \(p\le w_0\) dividing \(n\), the option
of a common divisor supplies, together with the target enhancement, at
least
\[
 \frac{p-1}{p-2}+\frac1p=1+2a_p .
\]
The product over these split primes is comparable to \(\log w_0\).
Under the normalized cofactor weights the expected \(\log(duv)\) is
at most \(C'\log w_0\), uniformly in \(n\). Taking
\(0<\lambda_0<\min(\eta,\kappa/(3C'))\) retains at least half the
model mass on (2) of EisSemigroupScale by Markov's inequality. At
inert target primes the combined shell and roughness factors produce
\(J(n)\) up to absolute factors uniformly in their valuations. These
are precisely the scalar estimates (9)--(12) of EisSemigroupScale,
applied with the above \(X_m\). Since \(I_\phi\) is fixed, they give
total retained weight \(\gg_\varepsilon HJ(n)/\log n\).
[Target-divisor calculation and cofactor truncation](EISENSTEIN_POINTWISE_SEMIGROUP_SCALE.md#uniform-contribution-from-the-target-divisors).

## Shifted inert tail and distinct value pairs

For each cofactor \(d\) choose an Eisenstein integer \(\rho_d\) with
\(N(\rho_d)=d\). The map \((z,y)\mapsto(\rho_d z,\rho_d y)\) embeds the
residual coordinates into the original shell. At primes at most \(w_0\),
the original small split primes of valuation one in both norms, only
in the first, or only in the second recover \(d,u,v\), respectively.
All other such valuations are excluded by the sieve. Images of distinct
triples are therefore disjoint. The map preserves inert divisibility
at primes above \(z_0\), the prescribed orientation and the fractional
norm coordinate.

If an inert \(p\) divides the first original component norm, the
component equals \(p\beta\). Its squared radius \(N(\beta)\) lies between
\((x-H/2)/p^2\) and \((x+H/2)/p^2\). For the second component use
\(n-x\) in place of \(x\). In either annulus the squared-radius difference
is \(H/p^2\), its inner radius is nonnegative, and its outer radius is
at most \(\sqrt n/p\). The disjoint translated fundamental cells of
\(\mathbf Z[\omega]\) have fixed positive area and radius at most one.
Enlarging the two annulus boundaries by one bounds the number of possible
\(\beta\) by \(O(H/p^2+\sqrt n/p+1)\), including the case of an inner
radius less than one. The complementary component is positive and has at
most \(6\mathcal D(n)\) choices, where
\(\mathcal D(n)\le\max_{k\le n}\tau(k)=n^{o(1)}\).
Summing over \(z_0<p\le\sqrt n\), the loss in the entire weighted
family, whose weights do not exceed one, is
\[
 \ll \mathcal D(n)\{H/z_0+\sqrt n(1+\log n)\}=o(H/\log n).
\]
Indeed \(z_0=n^\eta\), and (1) together with \(M\ge H\) implies
\(H\ge n^{(1+\varepsilon)/2}\). This gives uniformity even for the
smaller permitted summands.
[Fundamental-cell and complementary norm bounds](EISENSTEIN_POWER_SHRINKING_WINDOWS.md#a-large-prime-tail-inside-the-central-window).

After removing these large inert primes, the \(3\)-condition and oddness
give exactly the allowed semigroup types in their stated orientation.
For residual norms \(A_1=uk,A_2=vj\), every split prime of \(k,j\)
exceeds \(n^\eta\), hence each has at most \(1/\eta\) such prime factors
with multiplicity. The signed coordinate count for a given positive norm
is \(6g\), \(g=1*\chi_3\), with \(g(p^e)=e+1\) at split primes and
\(g(3^e)=1\). After division by \(2^{\omega(uv)}\) the fibre of a
residual value pair has weight at most
\(36g(k)g(j)\le36\cdot2^{2/\eta}\). Each original value pair determines
its residual pair and triple. This proves the lower bound of (2).
[Coordinate fibres for the cofactor sieve](EISENSTEIN_POINTWISE_SEMIGROUP_SCALE.md#removing-the-remaining-inert-primes-and-counting-value-pairs).

For the upper bound sieve the first actual summand in the interval of
length \(H\) centered at \(x\). At each inert \(p\ge5\), residues \(0,n\)
are excluded, a single residue when \(p\mid n\). Counting each CRT residue
class in an interval gives its length divided by the modulus with
\(O(1)\) error independently of its start. The upper block sieve from
EisSemigroupScale therefore applies at this center too. Its fixed-power
cutoff can be chosen uniformly since \(H\ge n^{(1+\varepsilon)/2}\).
The prime product has size \(O(J(n)/\log n)\), with total counting error
\(o(H/\log n)\), proving (2).
[Upper sieve and singular factor](EISENSTEIN_POINTWISE_SEMIGROUP_SCALE.md#upper-bound-for-the-same-singular-factor).

## Consequences for prescribed summand sizes

- For fixed \(0<\lambda<1\), \(x=\lambda n\) and \(H=n^\theta\) with
  \(2/3<\theta<1\), (1) holds for any fixed
  \(0<\varepsilon<3\theta-2\) once \(n\) is large.
- For \(x=n^\gamma\) or \(n-n^\gamma\) with \(1/2<\gamma<1\), every
  \((1+\gamma)/3<\theta\le\gamma\) gives (2) for \(H=n^\theta\), by
  taking \(\varepsilon<3\theta-\gamma-1\).
- For \(1/2<\theta<1\), put \(Y=n^\theta\) and fix
  \(0\le a<b<\infty\). Restrict the first actual summand to \([aY,bY]\).
  Take \(x=(a+b)Y/2\) and a width \(H\asymp_{a,b}Y\) small enough
  that its interval lies in \([aY,bY]\) and \(H\le x\). For large \(n\),
  \(M=x\), and (1) holds for \(\varepsilon<2\theta-1\). This gives the
  lower bound. The direct upper sieve applies to the full interval of
  length \((b-a)Y\). Thus the count is
  \(\asymp_{a,b,\theta}YJ(n)/\log n\), including \(a=0\).
- For any fixed interval \(I\subset(0,1)\) of positive length, the count
  with \(s_1/n\in I\) is \(\asymp_I nJ(n)/\log n\). For the lower bound
  take a fixed interior subinterval of \(I\) of width at most its
  midpoint's distance to the nearest endpoint of \([0,1]\), then apply
  (1) at that scale. The upper sieve again applies to the full interval.

These consequences hold for each family of the theorem. At the center,
a fixed-power width still requires \(\theta>2/3\) in this argument. The
small-summand result has strict threshold \(\theta>1/2\). The theorem
gives an order of magnitude. The exact leading constant and an explicit
numerical onset remain open.

The [all-valuation joint coset theorem](EISENSTEIN_BOUNDED_PRIME_FACTORS.md#all-target-valuations-and-joint-cofactor-levels)
also applies to the directional bump above. Its aggregate remainder allows
\(\kappa+\xi<\varepsilon/2\) in place of the sufficient condition (7).
It leaves the location criterion (1) and central exponent range unchanged.
The proof above uses the earlier bound and does not require this improvement.

## Dependencies and review

The angular input is the class-number-one coefficient model, the
cuspidal all-ideal Hecke bound and its scalar norm-character projection
from EisPower and EisBoundedFactors. Equation (3) permits independent
choices of the finite test and the archimedean vectors. Equations (4)--(5)
supply a localization estimate uniform up to the stated moving endpoint.
Combining these with the cofactor, singular-prime and interval-sieve
arguments of EisSemigroupScale proves (2). The earlier fixed-region
angular theorem concerns normalized coordinate witnesses. Here the count
is of distinct summand values with a width depending on both \(n\) and
the position of the window.
[Coefficient-valued Hecke formula, Greenberg--Voight, Section 2](https://jvoight.github.io/articles/orthohecke-240123.pdf),
[angular spectral inputs and primary references](EISENSTEIN_POWER_SHRINKING_WINDOWS.md#cuspidality-at-positive-archimedean-weight),
[fixed-progression prime estimates, Loughran--Matthiesen, Lemma 2.4](https://arxiv.org/html/1904.12845v2#S2),
[fixed-region angular theorem EisAngle](EISENSTEIN_SIFTED_ANGULAR_DISTRIBUTION.md).

A separate internal reviewer reconstructed the directional projection,
the translated angular estimate, the level inequality and the shifted
annulus bound. External mathematical review and priority assessment
remain outstanding.
