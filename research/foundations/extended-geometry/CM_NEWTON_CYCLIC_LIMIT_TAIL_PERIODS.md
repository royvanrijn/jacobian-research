# Escaping residue periods at the cyclic limit

Consider \(Q=(z^2+1)^4/z\) and a conductor-18 phase
\(f=t^{18}Q+t^{17}B+t^{16}D\) with \(\deg B\le5\) and
\(D=-z^{-1}+c+ez\). Use fibre value \(\lambda\), put
\(\varepsilon=\lambda^{-1/18}\), \(u=\varepsilon t\), and
take marker \(H=t^{12}A+t^{11}C\). A normalized residue
period is an integral of \(H\omega/df\) multiplied by
\(\lambda^{5/18}\). The [compact cycle limit](CM_NEWTON_SPHERE_PERIOD_CYCLES.md)
gives its leading differential away from the branch disks.

## The disks at zero and infinity

Near \(z=0\) the scaled fibre equation is \(z=K(u,z)\), where
\[
 K(u,z)=u^{18}P(z)+\varepsilon u^{17}zB(z)
                         +\varepsilon^2 u^{16}zD(z),\qquad P=zQ.
\]
It has a smooth graph \(z=z(u,\varepsilon)\) in small disks,
because \(K_z-1\) is a unit. The normalized residue there is
\[
 \pm\frac{u^{12}A(z)+\varepsilon u^{11}C(z)}{K_z-1}\,du .
\]
It has a holomorphic primitive uniformly small on a shrinking
\(u\)-disk. In particular the distinguished path endpoints
at \(t=\pm1\) have \(u=\pm\varepsilon\) on this graph.
Differences of primitive values there tend to zero.

At infinity write \(z=w^{-18}\), \(u=w^7U\). The scaled
equation is
\[
    q(w)U^{18}+\varepsilon O(w^{29})U^{17}
            +\varepsilon^2 O(w^{94})U^{16}=1,\quad q(0)\ne0.
\]
The 18 roots in \(U\) extend as regular branches on a small
disk, uniformly near \(\varepsilon=0\). On these branches
the coherent differential has order at least zero for
\(\deg A\le5\). The correction term has order at least 47
for \(\deg C\le2\), before its factor of \(\varepsilon\).
These orders follow from the Jacobian factor \(w^6\)
in \((dz/z)\wedge du\). The two disks therefore have
holomorphic primitives with controlled endpoint values
for the three-correction-input problem.

## Divergence despite central exactness

**Proposition (a divergent correction period at a finite root).**
On a phase in this strip there are compact cycles in fibres
escaping near a finite root on which the normalized correction
period grows as \(\varepsilon^{-8/7}\) times a nonzero constant,
even though the correction differential on the central cover
is globally exact.

Take \(\alpha=i\), \(v=z-\alpha\), \(B=v^4\), \(D=v^2/z\),
and \(C=1-4z^2\), with \(A=0\). The central primitive is
\(\frac32 Y^{-12}(1+z^2)\), as differentiation using \(Y^{18}=Q\)
shows. In the scaled fibre equation set
\(u=\varepsilon^{-2/7}U,\ v=\varepsilon^{9/7}V\).
The \(B\)-term vanishes at the limit and the leading level is
\[
        qU^{18}V^4+dU^{16}V^2=1,\qquad q=16/i,\quad d=1/i.
\]
The normalized two-form before taking residue has an overall
factor \(\varepsilon^{-8/7}\), obtained from
\(\varepsilon u^{11}(dz/z)\wedge du\).
With \(W=UV\) the limiting level is \(U^{14}P_*(W)=1\),
where \(P_*=qW^4+dW^2\), and the residue after removing
the factor is, up to sign, \(5U^{11}dW/(14i)\).

The differential \(U^{11}dW\) has orders \(-5,2,14\)
at \(W=0\), the nonzero roots of \(P_*\), and infinity,
respectively. The expansions at zero increase in steps
divisible by seven and give no residue. If this differential
had a meromorphic primitive, its projection to the deck
character 11 would be \(U^{11}R(W)\) for rational \(R\), with
\[
 R'-\frac{11}{14}\frac{P_*'}{P_*}R=1 .
\]
At a finite pole of \(R\) the smallest Laurent exponent
cannot cancel. At a root of \(P_*\) the equation then
forces \(R\) to vanish, so \(R\) is a polynomial with
three distinct zeros. At infinity the leading coefficient
of the left side has multiplier \(\deg R-22/7\), forcing
degree one. This contradiction shows there is a nonzero
compact period, since vanishing of all periods of a
residue-free meromorphic differential gives a meromorphic
primitive by path integration. Choose a representative
away from the branch points. The implicit function
theorem in the chart \(W,U\) lifts it to the actual fibres,
because the limiting \(U\)-derivative of the level equation
is \(14/U\). Its rescaled integral converges to the
nonzero limiting period.

This \(C\) is a formal correction marker. Actual
inverse-Frobenius correction inputs in this example have
an additional root zero. In characteristic \(p\), the
deficit-\(p\) numerator is
\([\xi^p](Q+\xi B+\xi^2D)^m\). At \(z=i\) each summand with
\(v_2\) factors of \(D\) has order \(4m-2v_2\ge4m-p\).
Applying Cartier after multiplication by \(dz/z\) extracts
the inverse-Frobenius correction polynomial times \(dz/z\).
The order is at least \(\lfloor(4m-p)/p\rfloor\ge1\) for
the sufficiently large good primes in the strip class.
Thus the example concerns geometric control for the full
space of correction markers.

## The coherent exponent on an escaping cycle

**Proposition (a finite-root period at the coherent order).**
There is a phase in the same strip for which compact cycles
supported in a finite-root scaling region have a nonzero
normalized coherent period limit with \(C=0\). The limit
survives projection to the same cyclic character as the
coherent period on the central cover.

Let \(\alpha=i\), \(v=z-\alpha\), \(\tau=v/\alpha\), and choose
\[
 c=-4\alpha/3,\qquad e=2/3,\qquad B=v^2 B_2,\qquad
 B_2=-\tfrac2{\sqrt3}(4+5\tau^2-\tau^3).
\]
Then \(zD=-(1+2\tau^2)/3\). Put \(h=1/t=v^2(c_0+k)\), with
\(c_0=-B_2/(2D)\). Completing the square gives
\[
 k^2+R(v)=\lambda v^{32}(c_0(v)+k)^{18}/D(v),
 \qquad R(v)=Q/(v^4D)-B_2^2/(4D^2).
\]
Writing \(R(v)=v^5(r_0+r_1v+\cdots)\) and
\(\beta(v)=c_0(v)^{18}/D(v)\), the local coefficients satisfy
\[
 r_0\ne0,\qquad r_1/r_0=1/(4\alpha),\qquad
 \beta_1/\beta_0=19/\alpha.
\]
For the first two identities, expand
\((2+\tau)^2\sqrt{1+2\tau^2}/(1+\tau)\) through order six.
Its difference from \(4+5\tau^2-\tau^3\) starts with
\(-2\tau^5+(7/2)\tau^6\). The resulting \(R\) is a nonzero
constant times \(-16\tau^5-4\tau^6+O(\tau^7)\).

Put \(\mu=\lambda^{-1/54}\), \(v=\mu^2V\), \(k=\mu^5K\).
The limiting level is \(K^2=\beta_0V^{32}-r_0V^5\).
Along paths away from \(K=0\), the implicit lift at fixed
\(V\) has denominator \(K\) varying analytically in \(\mu\).
To order \(\mu^2\) its square changes by
\(\mu^2(\beta_1V^{33}-r_1V^6)\).

For the protected marker \(A=z(z^2+1)^2\), the original
residue, up to sign, comes from
\(h^4A\,dv\wedge dh/[z(Q+hB+h^2D-\lambda h^{18})]\).
After substituting and dividing by the level derivative,
it is
\[
 \frac{v^8(c_0+k)^4(A/v^2)\,dv}
      {z(2Dk-18\lambda v^{32}(c_0+k)^{17})}.
\]
Write \(a(v)=c_0^4(A/v^2)/(2zD)\), so \(a_1/a_0=6/\alpha\).
The leading coefficient at \(\mu^{13}\) is
\(a_0V^8dV/K\), and the coefficient at \(\mu^{15}\) is
\[
 \left[a_1-\frac{a_0}{2}
       \frac{\beta_1V^{27}-r_1}{\beta_0V^{27}-r_0}\right]
       \frac{V^9\,dV}{K}.
\]
Differentiating \(V^{10}/K\) shows that in periods the factor
\(\beta_0V^{27}/(\beta_0V^{27}-r_0)\) has the same effect
as \(5/9\) on \(V^9 dV/K\). The coefficient at
\(\mu^{15}\) therefore has periods equal to
\(2a_0/(3\alpha)\) times those of \(V^9dV/K\). The latter
is holomorphic and nonzero on the compact normalization.
Its orders at zero, the simple branch points, and infinity
are 14, 0, and 5, respectively.

Let \(\rho=\exp(2\pi i/54)\). Put
\(\delta(V,K)=(\rho^2V,\rho^5K)\) on the limit. A
counterclockwise turn of \(\lambda\) replaces \(\mu\) by
\(\rho^{-1}\mu\). Continuing a fixed path in the scaled
chart then gives \(\delta^{-1}\) in the original chart.
The pullback by \(\delta\) has characters \(\rho^{13}\)
and \(\rho^{15}\) on the two coefficient forms.
On a compact cycle with nonzero \(V^9dV/K\) period take the
weighted sum of its 54 translates with weights
\(\rho^{-15j}\), then continue that chain along the fixed
\(V\) paths. The \(\mu^{13}\) term cancels and the
\(\mu^{15}=\lambda^{-5/18}\) term survives. Its normalized
limit is nonzero and comes entirely from the chart shrinking
to \(z=i\). The transported character \(\rho^{-15}\)
is \(\exp(-2\pi i\,5/18)\), the same as for the central
coherent order.

The examples first construct compact fibre cycles. The
[sphere variation calculation](CM_NEWTON_SPHERE_ESCAPING_VARIATIONS.md)
shows that their relevant character components lie in the
span of actual sphere differences for both phases and
computes their joint ranks with the core. The
[coefficient checker](verify_cm_newton_cyclic_limit_tail_periods.py)
replays the exact expansion and differential factors with
SageMath 10.9.
