# Leading constants for primitive quadratic norm sums

Fix an imaginary quadratic order \(\mathcal O\), of discriminant
\(\Delta=f^2D_0<0\). Use \(P,\Lambda,S,\chi,J\) from
[the all-orders theorem](PRIMITIVE_NORM_ALL_ORDERS.md):
\(S=\{p:p\mid2\Delta\}\), \(\chi=\chi_{D_0}\), and \(J\) includes
only inert primes outside \(S\). All counts are ordered numerical
value counts. The order and every coefficient/modulus below are fixed.

This note proves the one-variable constants, averaged additive leading
term, auxiliary-weight removal and conductor-prefix normalization.
Its pointwise implication uses BinaryCoreCorrelation from Section 2.
That correlation is now supplied by
[the pointwise binary theorem](PRIMITIVE_NORM_BINARY_CORRELATION.md),
whose Section 6 discharges its two imported analytic inputs.
Consequently the explicit all-order pointwise formula below holds.
The class group and conductor bookkeeping introduce no further unknown
constant. The earlier positive-restoration argument alone does not
supply this leading term.

## 1. Unconditional one-variable and averaged constants

Let
\[
 L=8\prod_{\substack{p\in S\\p\ {\rm odd}}}p,\quad
 \alpha_L=\varphi(L)/L,\quad H=\ker(\chi\bmod L),\quad
 \mathcal C=\langle q:q\notin S,\ \chi(q)=1\rangle.
\]
The positive number
\[
 c_{\mathcal C}^2=
 \frac{L(1,\chi)}{\pi}
 \prod_{p\in S}(1-p^{-1})(1-\chi(p)p^{-1})
 \prod_{\substack{q\notin S\\\chi(q)=-1}}(1-q^{-2})
 \tag{LC.1}
\]
is fixed by the actual value indicator, with no ideal/coordinate weights.

Indeed its Dirichlet series satisfies the exact Euler identity
\[
 F_{\mathcal C}(s)^2=
 \zeta(s)L(s,\chi)
 \prod_{p\in S}(1-p^{-s})(1-\chi(p)p^{-s})
 \prod_{\substack{q\notin S\\\chi(q)=-1}}(1-q^{-2s}).
 \tag{LC.2}
\]
Thus \(c_{\mathcal C}=\lim_{s\downarrow1}
(s-1)^{1/2}F_{\mathcal C}(s)/\Gamma(1/2)\).
Selberg--Delange gives
\(\#(\mathcal C\cap[1,X])\sim c_{\mathcal C}X/\sqrt{\log X}\).
For any fixed set \(A\subset H\) of unit residues,
\[
 \#\{n\le X:n\in\mathcal C,\ n\bmod L\in A\}
 =\left(c_{\mathcal C}\frac{|A|}{|H|}+o(1)\right)
                          \frac X{\sqrt{\log X}}.
 \tag{LC.3}
\]
For nonempty \(A\) this is an asymptotic equivalent; for empty \(A\)
the count is identically zero.
To check this normalization, expand by all characters modulo \(L\).
The two twists \(1,\chi\) equal one on \(\mathcal C\).
Every other twist has both associated Dirichlet character factors
nonprincipal and hence no half-order pole at one; fixed-character
Selberg--Delange gives a smaller-order contribution. In particular
the denominator is \(|H|\), not \(\varphi(L)\).
The imported one-variable theorem and its Gamma normalization are
de la Bretèche--Tenenbaum, *Remarks on the Selberg--Delange method*,
Theorem 1.2 and equations (1.9)--(1.12), in the inspected author edition.
For each character twist the multiplicative coefficients have absolute
value at most one; the prime-power condition holds for every
\(\sigma>1/2\). The fixed-progression prime number theorem gives its
prime-mean hypothesis with arbitrary fixed logarithmic power.
The parameter is \(\rho=1/2\) for \(1,\chi\) and \(\rho=0\) otherwise;
in the latter case all displayed leading coefficients vanish.

Give \(\mathbb Z_p\) Haar measure one and define the local VALUE densities
\[
 d_p=\mu(\mathscr V_p),\qquad
 d_p(N)=\int_{\mathbb Z_p}
       \mathbf1_{\mathscr V_p}(a)
       \mathbf1_{\mathscr V_p}(N-a)\,da\qquad(p\in S).
 \tag{LC.4}
\]
These count possible numerical values, not coordinate preimages.
OrderLocalCosets makes each a finite rational computation:
if \(m=p^{2v_p(\Delta)+1}\), form the complete primitive value residue
set \(V\subset\mathbb Z/m\); then
\(d_p=|V|/m\) and
\(d_p(N)=m^{-1}\#\{a\bmod m:a,N-a\in V\}\).
Consequently \(d_p(N)\) is periodic and is positive exactly when
\(N\in\mathscr V_p+\mathscr V_p\).

**Theorem OrderMarginalAndAverage.** For every fixed order,
\[
 \#(\Lambda\cap[1,X])\sim\#(P\cap[1,X])
       \sim c_\Lambda\frac X{\sqrt{\log X}},\qquad
 c_\Lambda=\frac{2c_{\mathcal C}}{\alpha_L}\prod_{p\in S}d_p.
 \tag{LC.5}
\]
Moreover
\[
 \sum_{N\le X}R_P(N)\sim\sum_{N\le X}R_\Lambda(N)
                  \sim c_\Lambda^2\,\frac{X^2}{2\log X}.
 \tag{LC.6}
\]
If \(2\operatorname{Pic}(\mathcal O)\ne0\), then for every fixed \(0<\delta<\delta_*(2\operatorname{Pic}(\mathcal O))\), with \(\delta_*\) defined in the all-orders theorem, the individual hidden-class failures obey
\(\#((\Lambda\setminus P)\cap[1,X])
       \ll_{\mathcal O,\delta}X/(\log X)^{1/2+\delta}\);
if \(2\operatorname{Pic}(\mathcal O)=0\) that set is empty. In particular \(\delta=1/5\) is always available, with order-dependent implied constants.

**Proof.** Factor \(n=r c\), where \(r\) is supported on \(S\) and
\((c,L)=1\). The exact shadow decomposition from (O.11)--(O.13) is
\[
 \Lambda=\coprod_{r\in\langle S\rangle}
      r\{c\in\mathcal C:c\bmod L\in A_r\}.
\]
There are no good inert factors in a primitive local norm.
Let \(f_r=|A_r|/|H|\). Equation (LC.3) gives the coefficient
\(c_{\mathcal C}\sum_r f_r/r\).
This sum converges absolutely. To justify the limit, the core upper
bound is \(O(Y/\sqrt{\log(2Y)})\); for \(r\le\sqrt X\) apply dominated
convergence with \(\sqrt{\log X/\log(X/r)}\le\sqrt2\).
For \(r>\sqrt X\), use the trivial bound \(X/r\) and the finite-support
Rankin bound
\[
 \sum_{\substack{r>\sqrt X\\r\in\langle S\rangle}}r^{-1}
 \le X^{-1/8}\prod_{p\in S}(1-p^{-3/4})^{-1}.
\]
This also makes that scaled tail negligible.

Under independent local Haar measure, fixing all bad valuations to
\(r\) has mass \(\alpha_L/r\), and the unit cofactor satisfies the local
norm conditions with probability \(|A_r|/\varphi(L)\).
Cross-prime unit factors in \(r\) merely permute uniform unit residues.
Since \(A_r\subset H\) and \(|H|=\varphi(L)/2\),
\[
 \prod_{p\in S}d_p=\frac{\alpha_L}{2}\sum_r\frac{f_r}{r},
\]
proving the shadow mean (LC.5).

The triangle \(a,b>0,\ a+b\le X\), with these two regularly varying
one-variable counting functions, gives (LC.6) for \(\Lambda\).
Explicitly, subdivide \(a/X\) into a fixed mesh, apply the means away
from the two endpoints, and then shrink the endpoint strips; their
mass is \(O(\varepsilon X^2/\log X)\). The remaining integral is
\(\int_0^1(1-t)\,dt=1/2\).

The target factor has bounded average. Indeed
\[
 J(N)=\sum_{d\mid N}a(d),\qquad
 a(q)=\frac1{q-2}\quad(q\notin S,\ \chi(q)=-1),
 \quad a(q^k)=0\ (k\ge2),
\]
with \(a\) multiplicative and \(\sum a(d)/d<\infty\).
Thus \(\sum_{N\le Y}J(N)\ll Y\).
Summing OrderClassError, splitting \(N\le\sqrt X\) off trivially, gives
\[
 \sum_{N\le2X}(R_\Lambda(N)-R_P(N))
           \ll X^2/(\log X)^{1+\delta}.
 \tag{LC.7}
\]
This transfers (LC.6). If \(D(X)=\#((\Lambda\setminus P)\cap[1,X])\),
the same sum is at least \(D(X)\#(\Lambda\cap[1,X])\): every such
ordered pair has at least one hidden failure and sum at most \(2X\).
The factor \(\#(\Lambda\cap[1,X])\asymp_{\mathcal O}X/\sqrt{\log X}\)
is the shadow mean already proved above, without any principal-value
count. Dividing by it gives \(D(X)\ll_{\mathcal O,\delta}
X/(\log X)^{1/2+\delta}\), with the same \(\delta\) as in
OrderClassError, and then the principal mean follows. Thus the improved
range \(\delta<\delta_*(2\operatorname{Pic}(\mathcal O))\) incurs no
additional loss. This reasoning is independent of a pointwise leading
asymptotic and does not use the principal mean circularly.
\(\square\)

There is no extra \(\mathrm B(1/2,1/2)=\pi\) in (LC.6) or the pointwise
candidate. The halves are powers of a logarithmic singularity.
The numerical counting functions have power \(X^1\); the additive
archimedean integral for a fixed target is \(\int_0^1dt=1\).
The two Gamma factors are already in \(c_{\mathcal C}^2\).

## 2. The precise binary correlation

Set
\[
 C_{\rm in}=
 \prod_{\substack{q\notin S\\\chi(q)=-1}}
             \left(1-\frac1{(q-1)^2}\right).
 \tag{LC.8}
\]
For coprime \(S\)-supported fixed \(A,B\), and fixed subsets
\(U,V\subset H\) that are unions of unit-square states, put
\[
 \Gamma_{A,B;U,V}(M)=
 \frac L{|H|^2}\#\{(u,v)\in U\times V:Au+Bv=M\bmod L\}.
\]

**Hypothesis BinaryCoreCorrelation.** Uniformly as \(M\) tends to infinity
through the positive-density targets, for each such fixed data,
\[
 \#\{Ax+By=M:x,y\in\mathcal C,\ x\bmod L\in U,\ y\bmod L\in V\}
 =
 \left(c_{\mathcal C}^2\Gamma_{A,B;U,V}(M)C_{\rm in}+o(1)\right)
                  \frac{M J(M)}{AB\log M}.
 \tag{LC.9}
\]
A zero finite density has count identically zero.

This is the exact pointwise correlation input to the normalization
argument below. It is supplied by PointwiseBinaryCorrelations in
PRIMITIVE_NORM_BINARY_CORRELATION.md. The older FiniteGenusFilter
theorem supplies comparability on the same targets.

At a good inert prime the two values forbid
\(\nu_q(M)=2-\mathbf1_{q\mid M}\) residues. Relative to the two single
means, its factor is
\[
 \frac{1-\nu_q(M)/q}{(1-1/q)^2}
 =
 \begin{cases}
 q(q-2)/(q-1)^2,&q\nmid M,\\
 q/(q-1),&q\mid M.
 \end{cases}
\]
Their product is exactly \(C_{\rm in}J(M)\); good split primes have
factor one. This computes the candidate singular series but does
not prove (LC.9).

Two unconditional mean checks fix the normalization.
The average of \(\Gamma\) over \(M\bmod L\) is
\((|U|/|H|)(|V|/|H|)\).
The Cesàro mean of \(C_{\rm in}J(M)\) is one, since
\[
 \sum_d a(d)/d
 =\prod_{\substack{q\notin S\\\chi(q)=-1}}
       \left(1+\frac1{q(q-2)}\right)=C_{\rm in}^{-1}.
\]
The same average holds jointly with \(\Gamma\) by CRT.
One-variable means prove the corresponding cumulative pair count
\(c_{\mathcal C}^2(|U||V|/|H|^2)X^2/(2AB\log X)\).
This verifies a necessary average of (LC.9), without differentiating
an average into a pointwise theorem.

## 3. Conditional full leading term, with a proved conductor-tail bound

**Theorem ConditionalOrderLeadingTerm.** If BinaryCoreCorrelation holds
for the core of the fixed order, then on all sufficiently large
locally admissible targets,
\[
 R_P(N)\sim R_\Lambda(N)\sim
 \frac{4c_{\mathcal C}^2}{\alpha_L^2}
 \left(\prod_{p\in S}d_p(N)\right)
 C_{\rm in}J(N)\frac N{\log N}.
 \tag{LC.10}
\]
Equivalently the coefficient of \(J(N)N/\log N\) is
\[
 \frac{4L(1,\chi)}{\pi}
 \prod_{p\in S}\frac{1-\chi(p)/p}{1-1/p}
 \prod_{\substack{q\notin S\\\chi(q)=-1}}
       \left(1-\frac2{q(q-1)}\right)
 \prod_{p\in S}d_p(N).
 \tag{LC.11}
\]
All infinite products in this expression converge absolutely.
This theorem proves the transfer from the stated correlation.
Its premise is supplied by PointwiseBinaryCorrelations in the separate
binary-correlation note.

**Proof.** For exact bad prefixes \(r,s\in\langle S\rangle\), put
\(g=(r,s)\), \(A=r/g\), \(B=s/g\), \(M=N/g\).
Use the original sets \(A_r,A_s\). For fixed \(r,s\) with \(g\mid N\),
(LC.9) predicts a contribution, relative to \(J(N)N/\log N\), of
\[
 c_{\mathcal C}^2 C_{\rm in}
       \Gamma_{r,s}(N)\frac g{rs}.
 \tag{LC.12}
\]
There is no loss in \(J\), since \(g\) is bad-supported.
For \(g\nmid N\) the contribution is zero.

The infinite sum can be justified, rather than assumed.
Parametrize the integer line \(rx+sy=N\). Its parameter interval has
length at most \(N g/(rs)\). At every good inert prime both linear
coefficients are units, so the usual two-root Selberg upper bound
applies uniformly in \(r,s\), even when these prefixes grow:
\[
 R_{r,s}(N)
 \ll_{\mathcal O}
       \frac{N g J(N)}{rs\log N}+O_\epsilon(N^{0.45}).
 \tag{LC.13}
\]
For example use prime-index level \(z=N^{1/10}\) and the exact
root-count/endpoint argument of MixedPrimeSquareSieve with no marked
primes. The interval may be shorter than \(z\); the displayed endpoint
covers that case. All local residue sets can be dropped for this upper
bound. The good-prime target roots merge precisely when \(q\mid N\).
No dependence on the sizes of the bad-supported coefficients occurs
in CRT interval counting.

There are \(O_{\mathcal O}((1+\log N)^{2|S|})\) possible \(r,s\le N\).
Also
\[
 \sum_{r,s\in\langle S\rangle}\frac{(r,s)}{rs}
 =\prod_{p\in S}\sum_{e,f\ge0}p^{-\max(e,f)}
 =\prod_{p\in S}\frac{1+p^{-1}}{(1-p^{-1})^2}<\infty.
 \tag{LC.14}
\]
Thus the total prefix tail beyond any fixed bound is arbitrarily small
relative to \(NJ(N)/\log N\), after \(N\) tends to infinity;
the sum of all endpoints in (LC.13) is negligible.
The fixed-prefix asymptotics can therefore be summed by a finite-first,
tail-second argument.

To identify their sum, partition the local Haar line \(a+b=N\) by
the two exact valuation prefixes. Its stratum has measure
\[
 \frac g{Lrs}\#
   \{(u,v)\in A_r\times A_s:Au+Bv=M\bmod L\}
 =\frac{\alpha_L^2}{4}
       \Gamma_{r,s}(N)\frac g{rs}.
\]
This follows prime by prime by scaling out the common valuation;
one of \(A,B\) is a unit, so a unit-pair root has the usual unique
linear lift. Summing disjoint valuation strata gives
\[
 \sum_{r,s}\Gamma_{r,s}(N)\frac g{rs}
       =\frac4{\alpha_L^2}\prod_{p\in S}d_p(N).
 \tag{LC.15}
\]
The right side is positive with a fixed minimum on the finite set of
admissible bad residue classes. Equations (LC.12)--(LC.15) prove the
shadow assertion with a relative error tending to zero.
OrderClassError transfers it to \(P\). Finally combine (LC.1) and (LC.8);
\((1-q^{-2})(1-(q-1)^{-2})=1-2/(q(q-1))\).
\(\square\)

The finite lower-profile cover used for sharp order bounds cannot by
itself compute the leading coefficient. Formula (LC.15) includes every
valuation stratum; (LC.13) is what permits that passage for split
conductor tails.

### 3.1. Mixed primitive forms at one fixed discriminant

Use \(P_Q,\Lambda_Q\) and the conductor-safe ideal-class convention of
Section 4.1 in the all-orders note. First fix one discriminant
\(\Delta=f^2D_0\), and two primitive positive definite integral forms
\(Q_1,Q_2\) of that discriminant. They may be in different genera.
Keep \(S,L,\alpha_L,H,\mathcal C,\chi,J\) as above and define
\[
 d_{p;1,2}(N)=\int_{\mathbb Z_p}
       \mathbf1_{\mathscr V_{p,Q_1}}(a)
       \mathbf1_{\mathscr V_{p,Q_2}}(N-a)\,da,
 \qquad p\in S.
 \tag{FM.1}
\]
Haar measure gives \(\mathbb Z_p\) mass one. These are densities of
numerical values, not of coordinate preimages. Set
\[
 \mathfrak b_{D_0,S}=\frac{4L(1,\chi)}{\pi}
  \prod_{p\in S}\frac{1-\chi(p)/p}{1-1/p}
  \prod_{\substack{p\notin S\\\chi(p)=-1}}
                         \left(1-\frac2{p(p-1)}\right).
 \tag{FM.2}
\]

**Theorem FormBinaryLeading.** With BinaryCoreCorrelation, as supplied
by the proved fixed-modulus binary theorem under its stated ZF/G inputs,
uniformly through sufficiently large locally admissible targets,
\[
 R_{Q_1,Q_2}(N)\sim R_{\Lambda_{Q_1},\Lambda_{Q_2}}(N)\sim
 \mathfrak b_{D_0,S}
       \left(\prod_{p\in S}d_{p;1,2}(N)\right)
                       J(N)\frac N{\log N}.
 \tag{FM.3}
\]
Every nonadmissible target has zero count. Each bad density is computable,
rational and periodic. The coefficient depends only on the two genera;
there is no class-number divisor. The full asymptotic has a relative
\(o(1)\) remainder, with no new quantitative remainder inferred from
the class-error bound.

**Proof.** The form-filter lemma gives the separate sets
\(A_{r,1},A_{s,2}\subset H\), retaining their individual genus biases.
For fixed prefixes let \(g=(r,s)\), \(A=r/g\), \(B=s/g\) and, when
\(g\mid N\), \(M=N/g\). Put
\[
 \Gamma_{r,s;1,2}(N)=\frac L{|H|^2}
  \#\{(u,v)\in A_{r,1}\times A_{s,2}:Au+Bv=M\bmod L\},
 \tag{FM.4}
\]
and set it to zero otherwise. The original two filters remain after
cancelling \(g\). The core asymptotic contributes
\(c_{\mathcal C}^2C_{\rm in}\Gamma_{r,s;1,2}(N)g/(rs)\)
relative to \(J(N)N/\log N\).

The uniform upper bound (LC.13) still applies after dropping these two
filters: the coefficients are good-prime units, and inert roots merge
exactly when the prime divides \(N\). There are only polynomially
many logarithmic prefix pairs, whereas (LC.14) is summable. Thus the
fixed-prefix limits and summable tail may be taken in that order.

The mixed Haar identity is
\[
 \sum_{r,s}\Gamma_{r,s;1,2}(N)\frac{(r,s)}{rs}
       =\frac4{\alpha_L^2}\prod_{p\in S}d_{p;1,2}(N).
 \tag{FM.5}
\]
At a bad prime write \(i=v_p(r)\), \(j=v_p(s)\), \(k=v_p(L)\).
Each compatible unit-pair root modulo \(p^k\) has Haar measure
\(p^{-\max(i,j)-k}\) on the line. If \(i\le j\), prescribing the
second value modulo \(p^{j+k}\) fixes the required first residue;
the other case is symmetric. Therefore each global root has measure
\(g/(Lrs)\). Summing its roots gives
\((\alpha_L^2/4)\Gamma_{r,s;1,2}(N)g/(rs)\), because
\(|H|=\varphi(L)/2\). Sum the disjoint valuation strata to obtain
(FM.5). Finite value cosets give a positive minimum for the bad product
on its positive residue classes, so the resulting \(o(1)\) is relative
uniformly on admissible targets.

The mixed class error (F.3) transfers the shadow asymptotic to the
actual values. Expanding (LC.1) and \(C_{\rm in}\) gives (FM.2) by
\((1-p^{-2})(1-(p-1)^{-2})=1-2/[p(p-1)]\). No factor \(1/h\)
arises: saturation makes a chosen compatible class attainable for almost
all numerical values in its genus shadow; these value sets do not form
a disjoint class partition. \(\square\)

### 3.2. Different fixed conductors in the same quadratic field

Fix \(K\) of fundamental discriminant \(D_0\), orders
\(\mathcal O_i\) of discriminants \(\Delta_i=f_i^2D_0\), and primitive
positive definite forms \(Q_i\) of discriminants \(\Delta_i\), for
\(i=1,2\). Put
\[
 S=\{p:p\mid2\Delta_1\Delta_2\},\qquad
 J_S(N)=\prod_{\substack{p\mid N\\p\notin S,\ \chi_{D_0}(p)=-1}}
                               \frac{p-1}{p-2}.
 \tag{FM.6}
\]
Use this common \(S\) in (FM.1)--(FM.2) and in the definitions of
\(L,H,\mathcal C\). Let \(G_i=\operatorname{Pic}(\mathcal O_i)\),
\(B_i=2G_i\), and
\[
 \eta_{\rm pair}=\min_i\delta_*(B_i),\qquad\delta_*(0)=+\infty.
 \tag{FM.7}
\]
This convention records an exact zero class error in that slot.

**Theorem CommonFieldGenusInvariance.** Without ZF/G, uniformly on all
sufficiently large locally admissible targets, the mixed shadow and
actual counts have order \(NJ_S(N)/\log N\). For every fixed
\(0<\eta<\eta_{\rm pair}\), the all-target class error is
\[
 0\le R_{\Lambda_{Q_1},\Lambda_{Q_2}}(N)-R_{Q_1,Q_2}(N)
 \ll\frac{NJ_S(N)}{(\log N)^{1+\eta}},\qquad N\ge3.
 \tag{FM.8}
\]
If \(Q'_i\) is in the same genus as \(Q_i\) in its own discriminant,
then on their common sufficiently large admissible targets,
\[
 R_{Q_1,Q_2}(N)=R_{Q'_1,Q'_2}(N)
                            (1+O((\log N)^{-\eta})).
 \tag{FM.9}
\]
Nonadmissible targets have zero count. If both \(B_i=0\), the counts
are equal for every target and the class error is zero. Otherwise
\(\eta=1/5\) is a common permissible exponent. Constants and onsets
depend on the fixed forms/orders and \(\eta\), and no endpoint or
optimality claim is made.

**Proof.** Apply FormIdealShadow and FormGenusCoreFilters to each order,
with the common enlarged \(S\). Their separate local multipliers remain
units with product of Hilbert symbols one; hence both actual filters
lie in the same \(H\). No equality of their local genus biases is used.
Deleting the extra bad primes removes finitely many terms from each
class-prime family and does not change its harmonic density.

For an exceptional summand apply saturation in its own Picard group.
Mark only that summand's outside-maximal-subgroup family. Both shadows
exclude the same good inert primes; the marked family is split in the
common field. The same one-sided exact-one sieve proves the two class
errors at the two respective thresholds, and their sum gives (FM.8).
When one \(B_i=0\), that side contributes zero. There is no simultaneous
marking of two class partitions, and no joint ring-class prime
counting theorem is required.

For the classical lower bound use
\(e_p=1+2\max_i v_p(\Delta_i)\). Both value sets are unions of full
cosets modulo \(p^{e_p}\). The common shift of the two local values by
\(\pm p^{e_p}t\) bounds both valuations by \(e_p+1\), so there are
finitely many covering profiles. The same classical fixed-filter theorem
for the common split core gives a positive lower bound in a positive
profile. The unmarked inert sieve gives the upper bound. Finally,
componentwise genus equivalence makes the two shadows identical.
Subtract (FM.8) for the two pairs and divide by the classical positive
lower bound, exactly as in (F.6). This proves (FM.9), with no use of the
pointwise binary asymptotic. \(\square\)

**Corollary CommonFieldFormPairs.** With the same existing binary inputs
as FormBinaryLeading, the explicit asymptotic (FM.3) holds for the two
possibly different conductors, using (FM.6), the mixed local densities,
and \(\mathfrak b_{D_0,S}\). Classical order bounds and componentwise
genus invariance are supplied independently by the preceding theorem.
There is no identification of genera across discriminants and no
assertion for different quadratic fields or varying orders.

**Proof.** Use the common core and the separate prefix masks just proved.
For every fixed prefix pair the existing BinaryCoreCorrelation applies.
The bound (LC.13), the summable series (LC.14), and the mixed root measure
\(g/(Lrs)\) are unchanged. Hence (FM.5) holds with the enlarged \(S\),
and the Euler expansion gives precisely (FM.2). The one-sided class
errors in (FM.8) transfer the shadow asymptotic to actual values.
The bad product has a fixed positive minimum on the finite positive
residue classes; its vanishing makes the count identically zero.
The leading remainder remains relative \(o(1)\). \(\square\)

## 4. A proved auxiliary-weight removal for the binary problem

Averaging many small split-prime divisor weights can remove them with
a relative error tending to zero. This isolates a different possible
analytic route from formal inversion of the coordinate weights.

Write \(\mathcal L=\log M\). Fix \(K\ge1\).
For \(1\le j\le J\), \(1\le\nu\le K\), choose disjoint prime bands
\[
 \mathcal P_{j,\nu}=
 \{q:\exp(\mathcal L^{\theta_{j,\nu}})<q
        \le\exp(2\mathcal L^{\theta_{j,\nu}}),\
             q\notin S,\ \chi(q)=1\},
\]
where \(1/10<\theta_{j,\nu}<1/5\).
One may use the equally spaced grid
\(\theta_{j,\nu}=1/10+((\nu-1)J+j)/(10(KJ+1))\)
and \(J\le\log\mathcal L/(20K\log2)\).
For large \(M\) these bands are pairwise disjoint.
Let
\[
 V_{j,\nu}=\sum_{q\in\mathcal P_{j,\nu}}\frac1q,\quad
 W_j(n)=\prod_{\nu=1}^K
       \frac{\sum_{q\in\mathcal P_{j,\nu}}\mathbf1_{q\mid n}}{V_{j,\nu}},
 \qquad \overline W=J^{-1}\sum_jW_j .
 \tag{LC.16}
\]
Fixed-character prime number theory gives
\(V_{j,\nu}=(\log2)/2+o(1)\), uniformly over this grid.
In particular these normalizers are bounded above and below.

**Lemma BinaryAuxiliaryRemoval.** For fixed \(A,B,U,V\) as in (LC.9),
let \(R(M)\) denote its actual unweighted count and \(R_{W_j}(M)\)
the same count weighted by \(W_j(x)\). Uniformly in all targets \(M\),
and throughout the displayed range of \(J\),
\[
 \left|R(M)-\frac1J\sum_{j=1}^J R_{W_j}(M)\right|
 \ll_{\mathcal O,A,B,K}
           \frac{M J(M)}{\log M\,\sqrt J}.
 \tag{LC.17}
\]
This is unconditional. Therefore (LC.9) would follow from its version
with the small auxiliary weights \(W_j(x)\), with uniform \(o(1)\)
error over the \(J\) families; choosing \(J\asymp\log\log M\) removes
the average weight at cost \(O((\log\log M)^{-1/2})\).

**Proof.** Let \(\mathbb E\) be the product model in which
\(\mathbf1_{q\mid x}\) at the auxiliary good split primes is Bernoulli
with probability \(1/q\). There is no good split-prime restriction on
either value. The distinct bands are independent, so
\[
 \mathbb EW_j=1,\quad \mathbb EW_j^2\ll_K1,\quad
 \mathbb E(\overline W-1)^2\ll_KJ^{-1}.
 \tag{LC.18}
\]
The identities are exact except for the harmless bounded second moment.
Set \(Z(x)=(\overline W(x)-1)^2\). Expanding the square gives
\[
 Z(x)=\sum_Q c_Q\mathbf1_{Q\mid x},\quad
 Q\le\exp(O_K(\mathcal L^{1/5})),\quad
 \sum_Q|c_Q|\le\exp(O_K(\mathcal L^{1/5})),
 \tag{LC.19}
\]
with squarefree \(Q\) composed only of auxiliary split primes.
The bounds are independent of \(J\): before collecting primes the
absolute scalar coefficients are those of \(1-2\overline W+\overline W^2\),
and their total is four. Also
\(\sum_Q c_Q/Q=\mathbb EZ\).

Apply the direct Selberg upper sieve to the line \(Ax+By=M\), dropping
the fixed residue filters for this upper estimate.
For a squarefree product \(d\) of good inert primes, the forbidden
root set \(\Omega_d\) has \(\rho(d)=\prod_{q\mid d}
(2-\mathbf1_{q\mid M})\) residues modulo \(d\).
The interval for the integer-line parameter has length \(X=O_{A,B}(M)\);
congruence \(Q\mid x\) still has density \(1/Q\) in that parameter,
because \(Q\) is prime to \(B\).
Since \((Q,d)=1\), exact CRT counting gives
\[
 \sum_{\substack{\text{line parameter in its interval}\\
                   \text{parameter}\bmod d\in\Omega_d}}Z(x)
  =X\frac{\rho(d)}d\,\mathbb EZ
       +O\!\left(\rho(d)\exp(O_K(\mathcal L^{1/5}))\right).
 \tag{LC.20}
\]
Here the equality uses the nonnegative polynomial \(Z\), while its
signed coefficients are retained in the exact main term.

Use the same Selberg weights and index level \(z=M^{1/10}\) as in
MixedPrimeSquareSieve. The quadratic main term is at most
\(O(J(M)/\log M)\), uniformly in target valuations. Its total endpoint
is \(O_\epsilon(z^{4+\epsilon}\exp(O_K(\mathcal L^{1/5})))=O(M^{1/2})\)
for sufficiently large \(M\). Equations (LC.18)--(LC.20) imply
\[
 \sum_{\substack{Ax+By=M\\x,y\in\mathcal C\\
                          \text{fixed filters}}}Z(x)
 \ll \frac{M J(M)}{J\log M}+O(M^{1/2}).
\]
Finally the unweighted direct sieve gives
\(R(M)\ll M J(M)/\log M\).
Cauchy--Schwarz applied on this actual pair set proves (LC.17);
the endpoint is absorbed uniformly for \(J=O(\log\log M)\).
\(\square\)

This argument is inspired by the auxiliary-weight averaging in
Jiang--Lau--Lü, Section 8, but the above CRT/sieve proof is given
explicitly for this binary problem. Their prime-plus-ordinary-norm
population has a different local divisibility model \(1/(q-1)\);
our split-prime model is \(1/q\). The lemma alone does not evaluate
\(R_{W_j}\). Expanding \(W_j\) leaves an unweighted core indicator on
both remaining cofactors. The earlier \(\tau\mathbf1_{\mathcal C}\)
weighted theorem does not supply that correlation.

## 5. Gaussian calibration

For \(K=\mathbb Q(i)\), let \(\mathcal C=C_4\) be the positive integers
supported on primes \(1\bmod4\). If \(K_{\rm LR}\) is the
Landau--Ramanujan constant, (LC.1) gives
\[
 c_{\mathcal C}=\frac1{4K_{\rm LR}},\qquad
 P=C_4\sqcup2C_4,\qquad c_P=\frac3{8K_{\rm LR}}.
\]
At two,
\(\mathscr V_2=\{a:a=1\bmod4\}\cup\{a:a=2\bmod8\}\).
Its marginal density is \(3/8\), and
\[
 d_2(N)=
 \begin{cases}
 1/4,&N\bmod8\in\{2,3,6,7\},\\
 1/8,&N\bmod8=4,\\
 0,&N\bmod8\in\{0,1,5\}.
 \end{cases}
\]
Thus the leading coefficient of \(J(N)N/\log N\) is
\[
 A_G=\frac12\prod_{q=3\bmod4}
                 \left(1-\frac2{q(q-1)}\right)
\]
on \(\{2,3,6,7\}\bmod8\), and \(A_G/2\) on \(4\bmod8\).
For the odd core alone it is \(A_G\) on \(N=2\bmod4\).
The full coefficient, including \(J(N)\), has Cesàro mean \(c_P^2\),
in agreement with (LC.6).
These explicit pointwise formulas follow from the now proved
BinaryCoreCorrelation input, as stated in Section 6 of the binary note.

## 6. Comparison with the older analytic inputs

The exact deweighting identity for the coordinate weight
\(V=\tau\mathbf1_{\mathcal C}\) is
\(\mathbf1_{\mathcal C}=V*(\mu\mathbf1_{\mathcal C})\).
Its correction has no absolutely summable tail at exponent one.
Applying the weighted theorem only to small divisors and discarding
the long signed tail is therefore invalid.
The positive-restoration proof retains fixed positive fractions;
it establishes comparability, without a coefficient approaching one.

In Fourier language, the minor-arc contribution for two factors is
a signed integral of \(S_X(t)^2e(-Nt)\).
Its absolute \(L^2\) bound has size \(X/\sqrt{\log X}\), larger than
the desired \(X/\log X\). The extra factor available in a ternary
supremum-times-\(L^2\) argument is absent.
This describes a limitation of these inputs, not an impossibility
theorem for binary asymptotics.

The targeted primary review distinguishes the following statements.

- Ramaré--Viswanadham, arXiv:2509.23260v1:
  the set defined on p.2 is exactly the odd primitive Gaussian core.
  Theorem 1.5 is a Fourier bound and Theorem 1.6 is ternary.
  Neither is (LC.9). The displayed single-mean constant in their
  p.5 equation (6) is incompatible with (LC.2); it is not used here.
- Iwaniec, *The half dimensional sieve*, Acta Arith.29 (1976),
  pp.69--73, Theorem 3 and Corollary 2:
  the indicator correlation formulas retain explicit progression-error
  remainders. Dropping those remainders would assume the missing step.
- Freiberg--Kurlberg--Rosenzweig, arXiv:1701.01157v1,
  Conjecture 1.1, and Smilansky, arXiv:1208.2419v2:
  ordinary shifted sum-of-two-squares indicator correlations are
  conjectural inputs, with proved averaged singular-series results.
  They concern related shifted populations, not a proved equivalence
  with our primitive positive-sum problem.
- Dietmann--Elsholtz, *Sums of two numbers having only prime factors
  congruent to one modulo four* (2019): eventual existence, without
  an unweighted leading term.
- Jiang--Lau--Lü, arXiv:2611.12467v2, Theorem 1.1, p.2:
  an explicit unweighted prime-plus-ordinary-norm asymptotic at fixed
  shift one. It is a recent adaptation lead, not (LC.9).
  Its auxiliary-weight removal motivates Section 4.
  The variable-target determinant issue in its earlier bilinear
  argument must be resolved before that method is imported.

The [source ledger](primitive_norm_leading_constant_source_ledger.json)
records exact inspected versions and scope. This comparison does not
claim historical exhaustiveness, priority or a numerical effective
cutoff. The actual pointwise theorem is proved in the separate binary
note using its two precisely identified recent analytic inputs.


## 7. The variable-target issue in the prime-plus-norm method

This section explains a limitation of directly adapting that particular
proof. The separate dilation-graph argument resolves the binary theorem
through a different target-grid estimate; it does not assert that the
source theorem already covers the enlarged determinant.

The full parameter chain in Jiang--Lau--Lü Sections 3--8 was inspected
for this adaptation, not just its introductory theorem.
Their Theorem 4.1 controls a mean square with
\(\alpha_m=\mu(m)m^{iv}\) on a rough interval, bounded rough
\(\beta_n\), and extracted auxiliary products \(a,c\asymp Y=x^{o(1)}\).
It uses the paired equations
\(mn-1=ah,\ rs-1=ch\), giving
\[
                         cmn-ars=c-a.
\]
For a binary target \(T\asymp x\), they become
\[
 T-mn=ah,\quad T-rs=ch
       \quad\Longrightarrow\quad cmn-ars=T(c-a).
 \tag{LC.21}
\]
A direct substitution expands their small determinant/path-coordinate
bounds by a factor \(T\), outside the stated moment/error budget.
Endpoint primitivity, previously forced by roughness and a nonzero
determinant below the roughness threshold, also needs a new argument.

An index-\(T\) projective lattice
\(\{v:T\mid\det(P,v)\}\) is a possible way to normalize determinants.
This repair is not proved here: it must retain ambient primitivity,
the correct continuous frame phase, and all target-prime gcd cases.
Even granting its geometric part, its major-arc input changes.
Let \(D=cmn-ars\); normalized determinant Fourier analysis contains
\[
 \mathbf1_{T\mid D}e(\theta D/T)
       =\frac1T\sum_{h\bmod T}e((h+\theta)D/T).
\]
At \(\theta\) near \(u/q\), \(q\le(\log x)^B\), these frequencies
have denominators dividing \(Tq\), generally of target size.

More precisely, for a smooth determinant cutoff \(\psi\), the
unprojected bilinear model produces the major-arc expression
\[
 \frac1{T\mathcal V^2}\sum_{h\bmod T}
 \int_{\theta\in\mathfrak M}\int_{\mathbb R}\widehat\psi(\xi)
          |S_h(\theta,\xi)|^2\,d\xi\,d\theta,\qquad
 S_h=\sum_{c,m,n}\eta(c/Y)\alpha_m\beta_n
  e\!\left(\frac{h+\theta+\xi/Y}{T}cmn-\theta c\right),
 \tag{LC.22}
\]
where \(\mathcal V\) is the auxiliary-band normalizer.
This is the exact finite Fourier expansion of the divisibility
condition and the smooth cutoff, before the additional endpoint
projections. The required estimate concerns this actual \(h\)-average;
replacing it by a supremum over characters would be stronger.

Source Section 7.2, equations (7.3)--(7.8), uses Möbius cancellation
only for character conductors at most a fixed power of \(\log x\),
as specified in Lemma 7.2. It supplies no estimate for (LC.22) with
\(T\asymp x\). Retaining the original determinant Fourier variable
moves the resonances \(h/T\) to other arcs instead of removing them.
This is a concrete additional major-arc/bilinear obligation, without
an impossibility claim. The source's recent zero-free inputs are not
imported here; regardless of their status they do not extend its
stated conductor range to \(Tq\).

There is a second independent interface. Direct deweighting of our
core introduces \(\mu_{\mathcal C}=\mu\mathbf1_{\mathcal C}\), while
Theorem 4.1 specifically uses \(\mu\). Principal twists of
\(\mu_{\mathcal C}\) have a fractional-power model and do not satisfy
the source's arbitrarily strong logarithmic cancellation without
extracting that model. A complete Type I/II identity with its actual
main terms is therefore required as well. The auxiliary-removal
lemma (LC.17) solves neither of these two evaluation problems.
