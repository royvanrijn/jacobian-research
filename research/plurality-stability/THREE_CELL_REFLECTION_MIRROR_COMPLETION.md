# An exact mirror score identity forces label transposition

Let \(A\) be an actual balanced planar Gaussian three-cell
noise-stability maximizer at \(0<r<1\). Suppose \(A_i\) is invariant
under reflection \(R\) in a line \(L\) through the origin, up to null
sets. For the other two labels \(j,k\), write \(U_l=T_r1_{A_l}\).
If \(U_j=U_k\) on \(L\), reflection exchanges \(A_j,A_k\) up to
null sets. In particular, if the invariant cell is bounded, the
reflection mismatches of the other two labels have positive Gaussian
mass arbitrarily far in both directions parallel to the mirror.

## The mass-preserving completion

Choose coordinates \(x=(u,v)\) with \(L=\{v=0\}\), let \(H=\{v>0\}\),
and use the restriction of \(d\gamma_2\) on \(H\) without renormalizing
it. Set \(g=1_{A_j}-1_{A_k}\) and define its even and odd parts on \(H\)
by
\[
 e(x)=\tfrac12[g(x)+g(Rx)],\qquad
 s(x)=\tfrac12[g(x)-g(Rx)].
\]
An orbit in \(A_i\) contributes \(e=s=0\). On the complement, an
orbit with two copies of the same label contributes \(e=\pm1,s=0\).
An orbit using both labels contributes \(e=0,s=\pm1\). Thus
\(|e|+|s|=1_{A_i^c}\) on \(H\). The mismatch set here is
\(E=\{x\in H:|e(x)|=1\}\), the orbits in the complement that fail
to exchange labels. Balance gives \(\int_H e\,d\gamma_2=0\), so the
two signs of \(e\) occupy equal Gaussian mass on \(H\).

With \(M_r\) the Mehler density relative to the Gaussian marginal,
write \(K^\pm(x,y)=M_r(x,y)\pm M_r(x,Ry)\). Use
\(\langle p,Kq\rangle_H\) for the corresponding quadratic or
bilinear integral over \(H^2\). Reflection parity gives
\[
 \tfrac12\langle g,T_rg\rangle
     =\langle e,K^+e\rangle_H+\langle s,K^-s\rangle_H .          \tag{MC1}
\]
The left side is the only term in total stability that changes when
the labels \(j,k\) are reassigned with \(A_i\) fixed. Indeed their
combined indicator is \(1_{A_i^c}\), and their combined stability is
half the energy of this indicator plus half the energy of \(g\).

The kernel \(K^-\) is pointwise positive on \(H^2\). On the orbits
already using both labels choose \(j\) on \(H\) and \(k\) on \(RH\).
This preserves each label's mass and replaces \(s\) by \(|s|\) in
(MC1), weakly increasing stability. Now use the same assignment on
every orbit in \(E\). This also preserves each label's total mass,
since the previous \(j,j\) and \(k,k\) orbits have equal mass. The
completed partition puts half of the Gaussian mass of \(A_i^c\) in
each of \(j,k\). Its folded contrast has even part zero and odd part
\(|s|+|e|\). Its gain over the sorted partition is exactly
\[
 \langle |e|,K^-|e|\rangle_H-\langle e,K^+e\rangle_H
       +2\langle |s|,K^-|e|\rangle_H .                        \tag{MC2}
\]
The unadjusted mirror score \(U_j-U_k\) is identical before and
after sorting the orbits using both labels. Its odd part vanishes
on the mirror, and its even part is unchanged by sorting.

## The mirror identity removes the constant term of the even kernel

Put \(\sigma=\sqrt{1-r^2}\), \(a=r/(1-r^2)\), \(b=r^2/(1-r^2)\),
and write \(M_r^{[1]}\) for the scalar relative Mehler density.
Factorization across the two coordinates gives, for
\(x=(u,v),y=(z,w)\in H\),
\[
 \begin{aligned}
 W(x,y)&=\frac{2}{\sigma}M_r^{[1]}(u,z)
                                  e^{-b(v^2+w^2)/2},\\
 K^+(x,y)&=W(x,y)\cosh(avw),&
 K^-(x,y)&=W(x,y)\sinh(avw).
 \end{aligned}                                                \tag{MC3}
\]
On the mirror, the same factorization gives
\[
 (U_j-U_k)(u,0)=\frac{2}{\sigma}T_r^{[1]}f(u),\qquad
 f(z)=\int_0^\infty e(z,w)e^{-bw^2/2}\,d\gamma_1(w).
                                                                    \tag{MC4}
\]
The transverse conditional variance on the mirror is \(\sigma^2\),
so the weight \(e^{-bw^2/2}d\gamma_1(w)\) is \(\sigma\) times the
centered Gaussian density of that variance.

The hypothesis \(U_j-U_k=0\) on \(L\) makes \(f=0\) almost
everywhere because \(T_r^{[1]}\) has strictly positive Hermite
eigenvalues \(r^m\) on \(L^2(\gamma_1)\). Consequently the folded
integral with kernel \(W\) alone vanishes by Fubini. Thus
\[
 \langle e,K^+e\rangle_H
   =\iint_{H^2} W(x,y)[\cosh(avw)-1]e(x)e(y)
                                          \,d\gamma_2(x)d\gamma_2(y).
                                                               \tag{MC5}
\]
For \(t>0\), \(0<\cosh(t)-1<\sinh(t)\) since
\(\sinh(t)-\cosh(t)+1=1-e^{-t}\). The integrals are finite by
boundedness of the label contrasts and integrability of the Mehler
density. Replacing \(e(x)e(y)\) in (MC5) by \(|e(x)e(y)|\) and then
applying the strict kernel inequality gives
\[
       \langle e,K^+e\rangle_H<\langle |e|,K^-|e|\rangle_H
\]
whenever \(E\) has positive measure. Equation (MC2) would then give
a strict improvement, which contradicts maximality. Hence \(E\) is
null. Sorting is also strict in (MC1) if both signs of \(s\) occur
on sets of positive measure. After a possible interchange of
labels, \(j\) therefore occupies \(H\setminus A_i\) and \(k\)
occupies \(RH\setminus A_i\). This proves transposition symmetry of
the original partition.

## Constant mirror pressure and a bounded invariant cell

Let the actual adjusted scores be \(V_l=U_l-\beta_l\) with
constant volume prices. On any open mirror segment of the \(j,k\)
interface they give \(U_j-U_k=\beta_j-\beta_k\). Analyticity extends
that equality to all of \(L\). Write \(\lambda=\beta_j-\beta_k\).
By (MC4) the weighted fibre mean is \(f=\sigma\lambda/2\) almost
everywhere. In this case the constant term removed in (MC5) equals
\(\langle e,We\rangle_H=(2/\sigma)\langle f,T_r^{[1]}f\rangle
=\sigma\lambda^2/2\). The same pointwise comparison and (MC2) give
the necessary inequality
\[
 \frac{\sigma\lambda^2}{2}\ \ge\
 \iint_{H^2} W(x,y)(1-e^{-avw})|e(x)e(y)|\,d\gamma_2(x)d\gamma_2(y)
   +2\langle |s|,K^-|e|\rangle_H.                              \tag{MC6}
\]
[Fixed-mass score assignment](MAXIMIZER_PERTURBATIVE_RIGIDITY.md#assignment-and-the-exact-information-supplied-by-global-maximality)

If \(E\) is supported in \(\{u\le B\}\) for some finite \(B\),
then \(f\) has the same support bound, up to null sets.
Since \(f\) is constant, it vanishes. The same holds for a bound
\(\{u\ge B\}\). In either case the exact unadjusted mirror identity
follows once the adjusted scores agree on a mirror segment.

Suppose \(A_i\) is bounded and \(E\) has such a one-sided bound.
Choose the parallel direction so both are absent almost everywhere
from a half-plane \(\{u>B_0\}\). Sorting at the maximizer puts \(j\)
on one side of \(L\) and \(k\) on the other throughout this
half-plane, up to null sets. At every point of the mirror ray
\(\{(u,0):u>B_0\}\), approach from both sides by points of the
respective almost-everywhere winning labels. Continuity of the
adjusted scores makes both labels maximal there. The mirror
identity then forces exact global transposition. The hypotheses
of (RF1) hold because the partition is an actual balanced maximizer
in its active plane and is invariant under \(R\) combined with
the transposition of \(j,k\). The fixed cell \(i\) is therefore an
unbounded monotone graph region, contradicting its boundedness.
[Reflection reduction (RF1)](THREE_CELL_REFLECTION_REDUCTION.md#orientation-and-gaussian-rearrangement)

The completion changes actual Gaussian partitions and preserves all
three label masses. The mirror kernel comparison needs only the
reflection of the fixed cell, the two remaining masses to be
equal, and the exact unadjusted score identity. The bounded-cell
consequence uses the balanced reflection reduction (RF1), which
uses the Gaussian two-set comparison and its parallel-halfspace
equality cases. [Mossel and Neeman, Theorems 1.1 and 1.3](https://ems.press/content/serial-article-files/32017)

## A bisected cell and a one-sided binary exterior

There is a version with no reflection symmetry of \(A_i\).
Suppose the centered line \(L\) bisects \(A_i\) in Gaussian
mass. If for some finite \(B\) the half-plane \(\{u>B\}\)
has label \(j\) above \(L\) and \(k\) below \(L\), up to null
sets, then \(j,k\) are separated by \(L\) globally and every
strict-winner component is unbounded. In particular a bounded
cell bisected by \(L\) cannot have such an exterior.

Indeed the exterior mirror ray makes \(U_j-U_k\) constant on
\(L\) by adjusted-score equality and analyticity. Formula (MC4)
holds for every contrast \(g\), regardless of how \(A_i\) meets
reflection orbits. Its \(f\) vanishes for \(u>B\) because
\(g(u,v)=\operatorname{sgn}(v)\) there. Hence the constant is
zero and \(f=0\) everywhere up to null sets. For the split of the
complement of \(A_i\), apply the orientation comparison from
(CD5). In detail, the unadjusted line equality lets one subtract
the transverse constant term in the full-plane Mehler energy of
\(g\). After absorbing \(\operatorname{sgn}(v)\) into \(g\), the
remaining kernel has positive factor
\((e^{avw}-1)\operatorname{sgn}(vw)\). Replacing
\(\operatorname{sgn}(v)g\) by its absolute value increases this
energy strictly if both signs occur on positive mass. The
original constant term is zero, and the replacement
\(g_*=\operatorname{sgn}(v)|g|\) has nonnegative constant-term
energy. The exterior supplies positive mass of the correct
sign. Since \(L\) bisects \(A_i\), it also bisects \(A_i^c\),
so \(g_*\) gives balanced cells. Maximality makes the
incorrect-sign set null. The affine-separation theorem then gives
the claimed unboundedness.
[Constant-kernel subtraction and orientation (CD5)](THREE_CELL_COMPACT_SIMPLEX_DISCREPANCY.md#a-spectator-bisector-forces-the-simplex)
[Affine separation and unbounded strict winners](THREE_CELL_AFFINE_SEPARATION_RIGIDITY.md)
