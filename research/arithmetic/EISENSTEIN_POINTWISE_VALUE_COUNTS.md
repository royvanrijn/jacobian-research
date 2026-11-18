# Pointwise counts of distinct Eisenstein decompositions

The coordinate supply in the power-window theorem gives a logarithmic lower
bound for the number of distinct value pairs. The additional input is an
upper bound for a fixed higher moment of the coordinate multiplicity. The
argument applies to every sufficiently large admissible target, including
targets with arbitrarily large prime-power factors. For central power windows,
the separate [cofactor theorem EisSemigroupScale](EISENSTEIN_POINTWISE_SEMIGROUP_SCALE.md)
now gives the stronger matching order \(n^\theta J(n)/\log n\) throughout
\(2/3<\theta<1\). The moment argument here also covers arbitrary fixed
macroscopic intervals as stated below. The separate
[moving-window theorem EisLocationWindows](EISENSTEIN_MOVING_WINDOWS.md)
now supplies the matching order \(nJ(n)/\log n\) in every fixed positive
macroscopic interval, and \(HJ(n)/\log n\) in its stated translated
shrinking windows.

Let
\[
 A=\{a\ge1:p\mid a\Longrightarrow p\equiv1\pmod3\},
\]
including \(1\), and put
\[
\begin{split}
 r_{0,H}(n)&=\#\{(a,b)\in A^2:a+b=n,\ |a-b|\le H\},\\
 r_{1,H}(n)&=\#\{(a,b)\in A^2:a+3b=n,\ |a-3b|\le H\}.
\end{split}
\]
These counts give each ordered value pair weight one. The mixed count refers
to the actual summand orientation \((a,3b)\); reversing that orientation has
the same count.

**Theorem EisValueCounts.** For every fixed \(2/3<\theta<1\) and
\(\nu>0\), there are \(c_{\theta,\nu}>0\) and
\(n_0(\theta,\nu)\) such that, for every integer \(n\ge n_0\),
\[
 r_{0,n^\theta}(n)\ge
 c_{\theta,\nu}\frac{n^\theta}{(\log n)^{\log4+\nu}}
 \quad(n\equiv2\pmod6),
\]
\[
 r_{1,n^\theta}(n)\ge
 c_{\theta,\nu}\frac{n^\theta}{(\log n)^{\log4+\nu}}
 \quad(n\equiv4\pmod6).
\]
Here \(\log4=2\log2\simeq1.38629\) is a natural logarithm. The constants
are uniform in the individual target. No numerical onset is supplied for this
value-count theorem; the coordinate theorem retains its explicit cutoffs.

## A uniform higher-moment upper bound

Let \(r_E(m)\) count ordered integer representations by
\(x^2-xy+y^2\), and put \(g(m)=r_E(m)/6\) for positive \(m\).
Its prime-power coefficients are
\[
 g(p^e)=
 \begin{cases}
 e+1,&p\equiv1\pmod3,\\
 \mathbf1_{2\mid e},&p\equiv2\pmod3,\\
 1,&p=3.
 \end{cases}
\]
For a fixed \(1<s\le2\), write \(A_s=2^s\) and
\[
 \mathcal T_s(n)=
 \prod_{\substack{p\mid n\\p\ge5,\ p\equiv1\ (3)}}
       \left(1+\frac{(A_s-1)^2}{p}\right)
 \prod_{\substack{p\mid n\\p\ge5,\ p\equiv2\ (3)}}
       \left(1+\frac1p\right).
\]

**Lemma EisMoment.** For every fixed \(0<\theta<1\), with \(H=n^\theta\),
\[
 \sum_{\substack{a+b=n\\a,b>0,\ |a-b|\le H}}
       g(a)^s g(b)^s
 \ll_{s,\theta} H(\log n)^{A_s-2}\mathcal T_s(n).
 \tag{M}
\]
This sums over all positive values, with zero weight at a value not represented
by the norm. It has no prime-factor sieve restriction and majorizes the
higher moment on either family in the theorem.

Apply Henriot's **New Theorem 5**, p.377, in the
[2014 erratum](https://doi.org/10.1017/S0305004114000280), using the growth-class
and coefficient-majorant definitions in the
[original paper, Section 2](https://arxiv.org/html/1102.1643v1#S2).
Take \(Q_1(t)=t\), \(Q_2(t)=n-t\) and
\(F(u,v)=g(u)^s g(v)^s\). The primitive irreducible linear factors have
primitive product, total degree two, and coefficient norm \(O(n)\).
For every fixed \(\zeta>0\), multiplicativity and the divisor bound give
\[
 g(p^e)^s\le(e+1)^s\le A_s^e,\qquad
 g(m)^s\le\tau(m)^s\ll_{s,\zeta}m^\zeta.
\]
Thus \(F\) belongs to Henriot's \(\mathcal M_2\) with fixed parameters,
and is its own coefficient majorant for these two factors.

Enlarge the closed summation window by one at each end and take
\[
 x=(n-H)/2-1,\qquad y=H+2.
\]
For sufficiently large \(n\), \(x\asymp n\), \(y\asymp H\), and
\(y\le x\). Choose short-interval exponent \(\alpha=\theta/2\),
coefficient-height exponent \(\Delta=1/2\), and growth exponent
\(\zeta=\theta/1600<\theta/800\). These satisfy the actual erratum's
restriction \(\zeta<\alpha/[50\cdot2(2+\Delta^{-1})]\), as well as
\(x\ge C_0\|Q_1Q_2\|^\Delta\) and \(y>x^\alpha\).
The enlargement avoids any half-open endpoint loss.

For clarity, evaluate the corrected Euler majorant directly. New Theorem 5
bounds the interval sum by a constant times
\[
 y\prod_{2<p\le x}(1-\rho_Q(p)/p)
 \sum_{d_1d_2\le x}g(d_1)^s g(d_2)^s
 \frac{\breve\rho(d_1,d_2)}
 {[d_1\operatorname{rad}(d_1),d_2\operatorname{rad}(d_2)]}.
\]
Drop the product bound, keeping primes at most \(x\), and factor the larger
nonnegative sum. Its local factor is
\[
 E_p=\sum_{a,b\ge0}g(p^a)^s g(p^b)^s P_p(a,b).
\]
The corrected \(\breve\rho\) gives \(P_p(0,0)=1\) with no conditions.
Every nonzero vector prescribes exact positive valuations, and a zero
coordinate opposite a positive coordinate prescribes a unit. The nonzero
valuation events are mutually disjoint.

For \(p\ge5\), put \(c_p=g(p)^s\), equal to \(A_s\) at split primes
and zero at inert primes. If \(p\nmid n\), the two roots are distinct,
\(\rho_Q(p)=2\), and
\[
 E_p=1+2(1-p^{-1})\sum_{e\ge1}g(p^e)^sp^{-e}
     =1+2c_p/p+O_s(p^{-2}).
\]
If \(p\mid n\), the corrected cross-coprimality kills vectors with just
one positive coordinate, and \(\rho_Q(p)=1\). The vector \((1,1)\)
has probability \((p-2)/p^2\) when \(v_p(n)=1\), or \((p-1)/p^2\)
otherwise. All terms with \(R=\max(a,b)\ge2\) contribute at most
\[
 \sum_{R\ge2}(2R+1)(R+1)^{2s}p^{-R}\ll_s p^{-2},
\]
since a unit-slope linear congruence prescribes one residue modulo \(p^R\).
Consequently
\[
 E_p=1+c_p^2/p+O_s(p^{-2})\qquad(p\mid n),
\]
uniformly in every target valuation. At \(p=2,3\), all prime-power
weights \(g(p^e)^s\) are at most one. Disjointness of the nonzero events
therefore gives \(E_p\le2\), independently of \(v_p(n)\). These finite
factors are absorbed into a constant.

Use the baseline \(B_p=(1-2/p)(1+c_p/p)^2\) for \(p\ge5\).
At primes not dividing \(n\), the ratio to the actual Euler factor is
\(1+O_s(p^{-2})\). At primes dividing \(n\),
\[
 \frac{(1-1/p)E_p}{B_p}
 =1+\frac{(c_p-1)^2}{p}+O_s(p^{-2}).
\]
The error products are bounded uniformly in the target. The finite small-prime
factors are also uniformly bounded by the convergent valuation sums above.
Adding the factors for primes of \(n\) above \(x\) only increases
\(\mathcal T_s(n)\).

The fixed-progression Mertens estimate
[Loughran--Matthiesen, Lemma 2.4(2)](https://arxiv.org/html/1904.12845v2#S2)
gives
\[
 \log\prod_{5\le p\le x}B_p
 =2\sum_{5\le p\le x}\frac{c_p-1}{p}+O_s(1)
 =(A_s-2)\log\log x+O_s(1),
\]
because the split and inert prime classes each have density \(1/2\).
This proves (M), including arbitrarily deep target valuations.

## From coordinate mass to value pairs

First note the uniform elementary estimate
\[
 \mathcal T_s(n)\ll_s(\log\log(3n))^{C_s},\qquad
 C_s=\max\{1,(A_s-1)^2\}.
 \tag{T}
\]
Indeed, majorize by \(\prod_{p\mid n}(1+C_s/p)\). Below \(\log n\),
ordinary prime Mertens bounds this by a constant times
\((\log\log n)^{C_s}\). Above \(\log n\), at most
\(\log n/\log\log n\) primes occur, each with reciprocal at most
\(1/\log n\), giving a bounded extra factor.

Fix \(2/3<\theta<1\) and put \(h=n^{\theta-1}\), \(H=n^\theta\).
[Theorem EisPower](EISENSTEIN_POWER_SHRINKING_WINDOWS.md) gives at least
\(hR_o(n)/20\) coordinate witnesses per prescribed good orientation,
where \(R_o(n)\ge18n\). For every actual value pair \((x,y)\), its
coordinate multiplicity is exactly \(36g(x)g(y)\). Thus the good value-pair
family \(\mathcal S\) satisfies
\[
 \sum_{(x,y)\in\mathcal S}g(x)g(y)\ge H/40.
\]
For the mixed orientation \((x,y)=(a,3b)\), the identity \(g(3b)=g(b)\)
holds on \(b\in A\). Hölder and (M) give
\[
 \#\mathcal S\gg_{s,\theta}
 \frac H{(\log n)^{(2^s-2)/(s-1)}
             \mathcal T_s(n)^{1/(s-1)}}.
 \tag{V}
\]
For any fixed \(\nu>0\), choose \(s>1\) close enough to one that
\((2^s-2)/(s-1)<\log4+\nu/2\). Equation (T) absorbs the remaining fixed
power of \(\log\log n\) into \((\log n)^{\nu/2}\). This proves
EisValueCounts for every sufficiently large admissible target.

The exponent \(\log4\) is the infimum available from these moments:
convexity gives \((2^s-2)/(s-1)\ge2\log2\) for every \(s>1\).
The expected \(H/\log n\) pointwise scale therefore needs additional
information beyond this moment-transfer argument. The separate cofactor
amplification provides that information for the central windows in
EisSemigroupScale and the translated windows of EisLocationWindows.

## Macroscopic intervals

**Corollary EisMacroscopicCounts.** For every fixed interval
\(I\subset[0,1]\) of positive length and every \(\nu>0\), the ordinary
or prescribed mixed family has
\[
 \gg_{I,\nu}\frac n{(\log n)^{\log4+\nu}}
\]
distinct ordered value pairs with first actual summand in \(nI\), for all
sufficiently large admissible targets. This moment bound is weaker than
the matching macroscopic order proved separately in EisLocationWindows.

Choose a fixed compact subinterval \(J\subset I\cap(0,1)\), short enough
that its length is less than its left endpoint. The
[angular theorem EisAngle](EISENSTEIN_SIFTED_ANGULAR_DISTRIBUTION.md)
supplies coordinate mass \(\gg_I n\) there. The proof of (M) applies to
this interval with \(x\asymp n\), \(y\asymp_I n\), and
\(\alpha=\Delta=1/2\). The Euler calculation is unchanged. The same
Hölder argument proves the corollary. In particular it applies with both
actual summands between \(n/3\) and \(2n/3\).

## Scope and attribution

The [bounded prime factor theorem](EISENSTEIN_BOUNDED_PRIME_FACTORS.md)
adds roughness and a fixed bound on the number of prime factors in each
permitted part, with a lower count of order \(n^\theta/(\log n)^2\).

The lower bound uses EisPower or EisAngle. Henriot's corrected theorem supplies
only the upper moment; no sieve lower-bound theorem is imported from it.
The independent-form correlation theorem is not applied to the proportional
linear parts of \(t\) and \(n-t\). The moment transfer itself gives the
weaker logarithmic exponent shown here; the cofactor theorem attains the
expected Eisenstein scale in the central and translated windows proved in
EisSemigroupScale and EisLocationWindows. The simultaneous modulus-12
every-target lower bound remains open.

The corrected local factors, target-valuation uniformity, finite factors at
2 and 3, interval parameters, coordinate fibers and Hölder conversion received
separate internal derivation and source audits. External specialist review,
formal verification and priority remain open. No new finite coverage search
is used in this proof.
