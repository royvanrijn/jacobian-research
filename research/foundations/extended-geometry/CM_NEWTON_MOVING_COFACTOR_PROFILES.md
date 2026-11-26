# Moving cofactor profiles at split real-cubic primes

Let \(B/K\) be an abelian threefold whose geometric rational endomorphism
algebra is a totally real cubic field \(F\), with all endomorphisms defined
over \(K\). At a good degree-one prime \(v\) with characteristic \(p\) split
in \(F\), write \(a_v\in\mathcal O_F\) for the relative Frobenius trace and put
\[
 b_v=\frac{N_{F/\mathbb Q}(a_v)}p,
 \qquad d_v=-\frac{N_{F/\mathbb Q}(a_v^2-4p)}p.
\]
Both are integral on the nonordinary locus. One failed factor gives
\(p\mid N(a_v)\) and \(p\mid N(a_v^2-4p)\). The Weil bounds are
\(|b_v|\le8\sqrt p\) and \(0\le d_v\le64p^2\).
See the [trace cofactor sieve](CM_NEWTON_TRACE_COFACTOR_SIEVE.md) for the
residual image and effective Chebotarev inputs.

## Square-root phase in a fixed Chebotarev class

**Theorem 1.** Let \(L/K\) be a fixed finite Galois extension of number
fields, \(C\) a conjugacy class in \(G=\operatorname{Gal}(L/K)\), and
\(c\ne0\) real. Among degree-one unramified primes \(v\) of \(K\) with
Frobenius in \(C\), the fractional parts \(\{c\sqrt{Nv}\}\) are
equidistributed. In particular the count in any fixed arc of length
\(\delta\) and with \(Nv\le X\) is
\((\delta |C|/|G|+o(1))\operatorname{Li}(X)\).

The cyclic reduction and ideal Vaughan identity below are the mechanisms
of Akbal and Güloğlu, [*Piatetski-Shapiro meets Chebotarev*, Lemmas 8–10](https://www.impan.pl/shop/en/publication/transaction/download/product/83752).
Their main theorem concerns Piatetski–Shapiro primes in a range of powers
near one. Here the phase is a fixed square root, and the target is a fixed
arc. The uniform growing-field theorem of Malik and Prabhu uses a zeta
zero-free hypothesis. [Malik–Prabhu, Theorem 1.1 and Remark 2](https://sites.math.rutgers.edu/~am2365/JointDistribution.pdf).

**Cyclic reduction.** Write \(e(t)=\exp(2\pi i t)\). Fix \(\sigma\in C\),
\(H=\langle\sigma\rangle\), and \(J=L^H\). A prime \(v\) of absolute degree
one whose Frobenius belongs to \(C\) has exactly
\(|C_G(\sigma)|/|H|\) primes of relative degree one in \(J\) with Frobenius
\(\sigma\) in \(L/J\). For another class there are none. Indeed these primes
correspond to cosets \(gH\) satisfying \(g^{-1}\operatorname{Frob}_v g=\sigma\).
Use character orthogonality in \(H\), and view its characters as finite
order ray class characters of \(J\). For the weighted Weyl criterion it
suffices to bound
\[
 \sum_{X<N\mathfrak a\le Y}
 \chi(\mathfrak a)\Lambda_J(\mathfrak a)e(c_0\sqrt{N\mathfrak a}),
 \qquad X\le Y\le2X,\quad c_0\ne0,
 \tag{1}
\]
by \(O_{J,\chi,c_0}(X^{1-\epsilon_0})\) for some \(\epsilon_0>0\).
Here \(\Lambda_J(\mathfrak p^m)=\log N\mathfrak p\) and is zero otherwise.
Higher-degree primes and higher prime powers contribute
\(O_J(\sqrt X(\log X)^2)\), and the characters may be extended by zero
on ideals meeting the modulus. Partial summation of (1) removes the
logarithm, and Chebotarev supplies the total count for the Weyl criterion.

**Bilinear sums.** For Dirichlet convolution on ideals let \(1\) be the
constant function, \(\mu\) the Möbius function, and \(\log\) the logarithm
of the norm. Subscripts \(\le U\) and \(>U\) truncate by norm. For
\(U=X^{1/10}\) use
\[
 \Lambda_J=\mu_{\le U}*\log+\Lambda_{J,\le U}
 -\mu_{\le U}*1*\Lambda_{J,\le U}
 +(\mu_{>U}*1)*\Lambda_{J,>U}.
\]
The second term vanishes on (1). The last term yields sums on dyadic
ranges \(m\asymp M\), \(n\asymp N\), \(MN\asymp X\), with \(m,n>U\),
subject to \(X<mn\le Y\), of the form
\(S=\sum_{m,n}\alpha_m\beta_n e(c_0\sqrt{mn})\).
The weights are obtained by aggregating ideals of equal norm and are
\(O_{J,\epsilon}(X^\epsilon)\) for every \(\epsilon>0\). Indeed the number
of ideals of norm \(n\) is bounded by the ordinary \([J:\mathbb Q]\)-fold
divisor function, and an ideal of norm \(n\) has at most
\(\tau(n)^{[J:\mathbb Q]}\) ideal divisors. These bounds control the
convolution weights and include the logarithmic factor in \(\Lambda_J\).

Exchange the two variables if necessary to take \(M\ge N\). Apply
Cauchy–Schwarz in \(m\) and expand the square. For \(n_1\ne n_2\), with
\(h=|n_1-n_2|\), the phase in \(m\) has second derivative of size
\(A_h/M^2\) with \(A_h\asymp \sqrt X h/N\). The second-derivative test
gives \(O(\sqrt{A_h}+M/\sqrt{A_h})\) on the interval where both products
qualify. [Graham–Kolesnik, Theorem 2.2](https://assets.cambridge.org/97805213/39278/excerpt/9780521339278_excerpt.pdf).
Summing the first term over pairs gives \(O(N^2X^{1/4})\), and summing the
second gives \(O(MN^2X^{-1/4})\), since
\(\sum_{1\le h\ll N}h^{-1/2}=O(\sqrt N)\). The diagonal contributes
\(O(MN)\). As \(M\gg\sqrt X\), this gives
\[
 |S|\ll_{J,c_0,\epsilon}
 X^{1+\epsilon}(N^{-1/2}+X^{-1/8}).
\]
There are only logarithmically many dyadic blocks, and \(N>U\) gives a
power saving for sufficiently small \(\epsilon\).

**The short-factor sums.** The remaining terms use an outer ideal of
norm \(m\le U^2\), and an inner sum
\[
 \sum_{X/m<N\mathfrak a\le Y/m}
    \chi(\mathfrak a)e(c_0\sqrt{m N\mathfrak a}).
 \tag{2}
\]
The \(\log N\mathfrak a\) in the first Vaughan term costs only a factor
\(O(\log X)\) by partial summation. We prove a bound
\(O((X/m)X^{-\eta_0})\) for (2), uniformly in \(m\le U^2\) and \(Y\),
with \(\eta_0>0\) depending on the fixed field.

Enlarge the modulus of \(\chi\) to include all real places and work
within one narrow ray class. Write \(D=[J:\mathbb Q]\) and choose an
integral inverse ray class representative \(\mathfrak c\), coprime to the
finite modulus \(\mathfrak f\). Write
\(\mathfrak a\mathfrak c=(\gamma)\) with \(\gamma\gg0\) and
\(\gamma\equiv1\bmod\mathfrak f\). The generators lie in one affine
coset of the lattice \(\mathfrak c\mathfrak f\) under the Minkowski
embedding. For \(\gamma\) restrict the normalized logarithms
\(\log|\gamma_i|-D^{-1}\log|N\gamma|\) to a fixed half-open
fundamental parallelotope for the logarithms of ray units. Each ideal
then occurs a fixed number of times equal to the number of torsion ray
units. The weighted sum of the normalized logarithms, with weights one
at real and two at complex places, is zero. For
\(T=(X/m)^{1/D}\), the generator domain scaled by \(T^{-1}\) lies in a
fixed compact region with all \(|\gamma_i/T|\) bounded above and away
from zero. The norm window and the upper endpoint \(Y\) only vary the
intersection with this region.

Take a fixed even integer \(k\ge D+2\). At each point in this region there
is a direction with nonzero \(k\)-th directional derivative of
\(g(z)=|N(z)|^{1/2}\). At a real place vary just that coordinate radially
to obtain \(g(z)(1+t)^{1/2}\). If there are only complex places, vary
one coordinate in its perpendicular direction to obtain
\(g(z)(1+t^2)^{1/2}\). The derivatives of order \(k\) at zero are
nonzero. Approximate the direction by that of a primitive vector in
the lattice \(\mathfrak c\mathfrak f\). By continuity a neighborhood of
the point has the absolute derivative bounded above and away from
zero along this lattice direction. Compactness and a sufficiently fine
fixed grid partition the scaled region into finitely many boxes, each
with one such direction.

Complete the chosen primitive vector to a lattice basis in each box.
The number of rows parallel to it that meet the box of scale \(T\) is
\(O(T^{D-1})\). Each row intersects the generator domain in a bounded
number of intervals. To see this, after removing the finitely many
coordinate zeros in the row, a logarithmic inequality defining the
unit parallelotope or norm band has derivative a linear combination
of rational functions. Their denominators are real linear functions
or the squared moduli of complex affine functions, which are real
quadratics. The degree bound is fixed by \(D\). The inequality is
constant or changes monotonicity only a bounded number of times.
The box and positivity constraints also give finitely many intervals.

The phase in (2) on one such interval has absolute \(k\)-th derivative
\(\asymp\lambda_k=X^{1/2}/T^k\), since the unscaled generator norm equals
\(N\mathfrak c\cdot N\mathfrak a\). In particular
\(\lambda_k\ll X^{-3/10}\). For a sufficiently small fixed \(\eta>0\),
intervals of length at most \(TX^{-\eta}\) cost
\(O(1+TX^{-\eta})\) each. An interval of length \(R>TX^{-\eta}\) has
\(R\gg X^{4/(5D)-\eta}\) and \(R^k\lambda_k\gg X^{1/2-k\eta}\).
Heath-Brown's \(k\)-th derivative estimate gives
\[
 O_{\epsilon}\!\left(R^{1+\epsilon}
 \left(\lambda_k^{1/(k(k-1))}+R^{-1/(k(k-1))}
 +(R^k\lambda_k)^{-2/(k^2(k-1))}\right)\right).
\]
Apply the estimate with the phase or its negative so the derivative
is positive. With \(\eta<\min(1/(4k),2/(5D))\), all three terms save a
fixed power of \(X\), and \(\epsilon\) can be chosen small enough to
retain a power saving. [Heath-Brown, Theorem 1](https://arxiv.org/pdf/1601.04493).
Summing rows, boxes and the finitely many ray classes gives the claimed
bound for (2). Finally the outer Vaughan weights aggregated by \(m\)
are \(O_{J,\epsilon}(m^\epsilon)\) as above. Summing
\((X/m)X^{-\eta_0}m^\epsilon\) over \(m\le X^{1/5}\), with the
optional logarithmic weight, gives \(O(X^{1-\eta_0/2})\) on reducing
\(\epsilon\). Together with the bilinear bound this proves (1). \(\square\)

## Bounded-width square-root profiles

**Theorem 2.** Fix \(\lambda\in\mathbb R\setminus\{0\}\) and \(H\ge0\).
The degree-one split nonordinary primes satisfying
\[
                   |b_v-\lambda\sqrt p|\le H
\]
have natural density zero, also relative to the degree-one split primes.
For fixed \(B,K,\lambda\), this holds for some nondecreasing
unbounded width \(H=H(X)\) when counting primes up to \(X\).

**Proof.** Fix an auxiliary rational prime \(\ell\) split in \(F\).
On the appropriate locus, \(b_v\bmod\ell\) is the conjugacy-invariant
function \(d^{-1}t_1t_2t_3\) on the residual torsion image, where
\(d\) is the common determinant and \(t_i\) the three traces. Discard
the finitely many ramified primes, including those above \(\ell\).
For a residue \(j\), the claimed archimedean inequality implies that
\(\{\lambda\sqrt p/\ell\}\) lies in the arc of radius \(H/\ell\) about
\(j/\ell\). Theorem 1 applies to the union of conjugacy classes
with residue \(j\). Summing the finitely many arcs and class densities
gives upper density at most \(2H/\ell\) for \(\ell>2H\). Endpoints
can be included by enlarging the arc by arbitrarily small amounts.
Let \(\ell\) tend to infinity through split auxiliary primes.
The split degree-one primes have positive density by Chebotarev.
For the growing width choose increasing thresholds beyond which
the counts with width \(j\), divided by \(\operatorname{Li}(X)\), are
at most \(1/j\). Use width \(j\) between successive thresholds.
This gives a diagonal width with no specified rate. \(\square\)

## Explicit widths for polynomial strips

**Proposition 3.** Let \(P\in\mathbb Q[z,s,b,d]\) be fixed. Assume the
substitution
\[
 Q(z,t_1,t_2,t_3)=
 P\!\left(z,\sum_i t_i,z^{-1}\prod_i t_i,
                  -z^{-1}\prod_i(t_i^2-4z)\right)
\]
is nonconstant as a Laurent polynomial in \(z\) and a polynomial in
the \(t_i\). Among degree-one split nonordinary primes in \((X,2X]\),
the count satisfying
\[
 |P(p,\operatorname{Tr}_{F/\mathbb Q}(a_v),b_v,d_v)|\le H
\]
is
\[
 O_{B,K,P}\!\left((H+1)\frac{X(\log\log X)^{1/6}}{(\log X)^{7/6}}\right).
\]
Examples are strips around \(d_v=\lambda p^2\) and
\(b_v=\lambda\operatorname{Tr}_{F/\mathbb Q}(a_v)\) for fixed rational
\(\lambda\). These strips can have growing width
\(H=o((\log X/\log\log X)^{1/6})\) and still have relative density zero.

**Proof.** Clear the fixed denominators of \(P\). At a sufficiently
large split auxiliary prime \(\ell\), the full common-determinant
residual image gives at most \(O(\ell^{-4})\) probability to each
tuple \((z,t_1,t_2,t_3)\) with \(z\ne0\). Indeed \(z\) is uniform on
\(\mathbb F_\ell^\times\); at fixed determinant each factor is uniform
in its determinant fiber, and given a trace there are at most
\(2\ell^2\) matrices. For the latter count write the matrix rows as
\((a,b),(c,t-a)\) and count solutions to \(bc=a(t-a)-z\).
For every prescribed value \(r\), the reduction of \(Q-r\) after
clearing the Laurent denominator is a nonzero polynomial of bounded
degree, for all sufficiently large \(\ell\). Schwartz–Zippel bounds
its fiber by \(O_P(\ell^3)\). An integer of absolute size \(O_P(H)\)
gives \(O_P(H+1)\) candidate residue values, so the corresponding
union has group mass \(O_P((H+1)/\ell)\).

Apply the classwise unconditional Chebotarev upper bound of the
[trace cofactor sieve](CM_NEWTON_TRACE_COFACTOR_SIEVE.md#counting-small-trace-cofactors).
Its choice \(\ell\asymp c(\log X/\log\log X)^{1/6}\) and the abelian
subgroups used for the classwise estimates depend on the residual
image, independent of the choice of the union of classes. This gives
the stated bound. \(\square\)

## Power widths under GRH

**Proposition 4.** Assume GRH for the Dedekind zeta functions of the
auxiliary torsion fields \(K(B[\ell])\). Fix
\(\lambda\in\mathbb R\setminus\{0\}\) and \(0\le\nu<1/42\). The count
of degree-one split nonordinary primes in \((X,2X]\) with
\(|b_v-\lambda\sqrt p|\le p^\nu\) is \(O(X^{1-\delta})\) for some
\(\delta>0\).

**Proof.** Choose \(\ell\asymp X^\beta\), split in \(F\), with \(\beta\)
strictly between \(\nu\) and \(1/42\). As in Theorem 2, for each
residual class the allowed arc has length \(O(H/\ell)\), now with
\(H=(2X)^\nu\). Cover it by at most two intervals of lengths comparable
to \(H/\ell\), enlarging endpoints as needed. Apply
[Malik–Prabhu, Theorem 1.1 and Remark 2](https://sites.math.rutgers.edu/~am2365/JointDistribution.pdf)
to the class with \(\theta=1/2\), \(\omega=1\), and
\(\alpha=|\lambda|/\ell\). For negative \(\lambda\), use \(-b_v\).
Here the absolute auxiliary degree is \(n_L\ll_{B,K}\ell^{10}\). The
parameter hypothesis follows from
\(\alpha^{1/4}(n_L/(H/\ell))^{1/2}\ll H^{-1/2}\ell^{21/4}\),
where the exponent of \(X\) is below \(1/8\). The main term and
the relative error of order \((H/\ell)X/\log X\) give a saving since
\(\beta>\nu\). Summing the other relative errors uses that the sum of
class masses is one, and yields
\[
 O\!\left(H^{1/2}\ell^{17/4}X^{7/8}\log X
       +H\ell^{19/2}X^{3/4}(\log X)^2\right).
\]
The remaining error is at most
\(O(H^{1/2}\ell^{57/4}X^{3/8}\log X)\), even on summing over at most
\(O(\ell^{10})\) separate conjugacy classes. All exponents of \(X\)
are strictly below one. \(\square\)

These restrictions control specific moving profiles, including profiles
with full Weil-size cofactor. They leave the density of the full
one-failure locus unresolved.
