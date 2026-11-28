# Finite-swap bounds at a Gaussian centroid mirror

Let \(A\) be an actual planar Gaussian three-cell noise-stability
maximizer at \(0<r<1\), with positive prescribed masses. Write
\(m_i=\int_{A_i}x\,d\gamma_2(x)\). Fix a reflection \(R\) through
the Gaussian origin, let \(H\) be an open half-plane of \(R\),
and let \(n\) be its inward unit normal. The measure on \(H\)
is the restriction of Gaussian probability, without renormalization.

If \(m_k\cdot n=0\), the other two labels \(i,j\) occur on opposite
members of reflected pairs on a set of positive Gaussian mass.
Moreover, either \(A_k\) is invariant under \(R\) up to null sets,
or all three unordered pairs of distinct labels occur on positive
Gaussian mass. The first occurrence bound is uniform over
maximizers on compact positive-correlation bands and compact sets
of positive mass triples. The other two bounds are uniform when
the reflection defect of \(A_k\) is bounded below.

## A finite-swap inequality for any mirror

Let \(F=(1_{A_i})_{i=1}^3\), \(\Delta(x)=F(x)-F(Rx)\) on \(H\),
and write \(W(x,y)=M_r(x,y)-M_r(x,Ry)\) with \(M_r\) the
Gaussian-relative Mehler density. The reflection-dependent part
of stability is
\[
     \frac12\iint_{H^2} W(x,y)\Delta(x)\cdot\Delta(y)
                                      \,d\gamma_2(x)d\gamma_2(y).
                                                               \tag{CM1}
\]
Indeed, summing the four kernel interactions of the two mirror
orbits splits the original integrand into this odd part and an
even part depending only on \(F(x)+F(Rx)\). Swapping the labels
at \(x,Rx\) leaves that sum and each label's total mass unchanged.
[The folded even and odd identity (RM2)](THREE_CELL_REFERENCE_MIRROR_MISMATCHES.md#a-centered-reflection-with-at-most-two-mismatch-types)

Distinguish any unordered pair of distinct labels and let \(E\)
be its mismatch set in \(H\), of Gaussian mass \(\delta\). Orient
the other two potential mismatch types out of their common label.
On their union write \(\Delta(x)=\varepsilon(x)d(x)\), where
\(d(x)\) is the oriented label difference and
\(\varepsilon(x)\in\{-1,1\}\). Let \(D_+,D_-\) be the two sign
sets. The inner product \(d(x)\cdot d(y)\) is at least one
whenever both points use these retained types. The absolute inner
product of either retained type with the distinguished type is
one, in either orientation.

Swapping only on \(D_-\) increases the retained cross interaction
in (CM1) by at least \(2\int_{D_+\times D_-} W\), and loses
at most \(2\int_{E\times D_-}W\) in interactions with \(E\).
Swapping only on \(D_+\) gives the same retained gain and loses
at most \(2\int_{E\times D_+}W\). Products of sets in these
integrals use \(\gamma_2\otimes\gamma_2\). The two ordered
cross terms and the sign reversal give the factor two after
the factor \(1/2\) in (CM1). The kernel \(W\) is positive for
points in \(H\), and \(\int_H W(x,y)\,d\gamma_2(y)\le1\) since
\(\int_{\mathbb R^2} M_r(x,y)\,d\gamma_2(y)=1\). Hence
\[
 \int_{D_+\times D_-}W
    \le \min\!\left(\int_{E\times D_+}W,\int_{E\times D_-}W\right)
    \le \frac{\delta}{2}.                                     \tag{CM2}
\]
If the first inequality fails for an ordinary partition, one
of the two specified swaps strictly improves stability at
exactly its original masses.

For \(R_0>0\), put \(v(x)=x\cdot n\) on \(H\), and
\(\omega_\pm(R_0)=\int_{D_\pm\cap B_{R_0}(0)}v\,d\gamma_2\).
Set \(\sigma^2=1-r^2\), \(a=r/\sigma^2\), \(b=r^2/\sigma^2\),
and \(c_{r,R_0}=2a e^{-(a+b)R_0^2}/\sigma^2\). With \(u\)
the signed coordinate along the mirror,
\[
 W(x,y)=\frac{2}{\sigma^2}
     e^{-b(|x|^2+|y|^2)/2+a u(x)u(y)}
                  \sinh(a v(x)v(y)).
\]
On \(B_{R_0}(0)^2\cap H^2\) this is bounded below by
\(c_{r,R_0}v(x)v(y)\). Thus
\[
             \delta\ge 2c_{r,R_0}\omega_+(R_0)\omega_-(R_0).
                                                               \tag{CM3}
\]

## Consequences at the centroid mirror of a label

Suppose \(m_k\cdot n=0\). Let \(D_k\) be all mismatch orbits
involving \(k\), put \(q=\gamma_2(D_k)=\gamma_2(A_k\triangle
R A_k)/2\), and write \(w_k=\int_{D_k}v\,d\gamma_2\).
The centroid equation is
\[
        0=\int_H v(x)[1_{A_k}(x)-1_{A_k}(Rx)]\,d\gamma_2(x).
                                                               \tag{CM4}
\]
The two signs of the difference thus each carry normal
moment \(w_k/2\). For any set \(B\subset H\) of measure \(p\),
\(\int_B v\,d\gamma_2\ge p^2/(2\phi(0))\), with \(\phi\)
the scalar standard Gaussian density. Indeed
\(\gamma_2(B\cap\{v>t\})\ge(p-t\phi(0))_+\), and integration
in \(t\) proves the claim. Write
\(\tau_{R_0}=\int_{H\setminus B_{R_0}(0)}v\,d\gamma_2\).
An arbitrary set of mass \(\delta\) in \(H\) has normal moment
at most \(\sqrt{\delta/2}\) by Cauchy–Schwarz and
\(\int_H v^2\,d\gamma_2=1/2\).

First distinguish the type \(\{i,j\}\), so the retained
types share \(k\). Put \(t=|m_i\cdot n|=|m_j\cdot n|\).
Full planar centroid rank gives \(t>0\). The identity
\(m_i\cdot n=\int_H v\Delta_i\,d\gamma_2\) gives
\(w_k\ge t-\sqrt{\delta/2}\), since \(\Delta_i\) vanishes
off \(D_k\cup E\) and has magnitude at most one. The
orientation in (CM2) makes the two normal masses each
\(w_k/2\). Therefore
\[
     \delta\ge 2c_{r,R_0}
           \left[\tfrac12(t-\sqrt{\delta/2})-\tau_{R_0}
           \right]_+^2.                                       \tag{CM5}
\]
Choose \(R_0\) so that \(\tau_{R_0}\le t/8\). If
\(\delta<t^2/2\), the bracket is at least \(t/8\).
In all cases
\(\delta\ge t^2\min(1/2,c_{r,R_0}/32)>0\).
On the stated compact parameter sets, full centroid rank
has a uniform lower singular value. Since \(m_k\cdot n=0\)
and \(\sum_i m_i=0\), this bounds \(t\) below. The Gaussian
tail \(\tau_{R_0}\) is independent of mirror orientation and
\(c_{r,R_0}\) has a positive lower bound on a compact
correlation band for fixed \(R_0\). This makes the lower
bound uniform.
[Uniform full moment rank](GAUSSIAN_NEAR_MAXIMIZER_COMPACTNESS.md#uniform-full-moment-rank-for-three-cells)

Next distinguish one of the two types involving \(k\).
The other type involving \(k\) remains among the retained
types, and on it the retained orientation gives the
negative of the sign in (CM4), since the common label is
the other label of this type. Removing \(E\) costs at most
\(\sqrt{\delta/2}\) of either of the two normal masses in
(CM4), so both retained signs have normal moment at least
\(w_k/2-\sqrt{\delta/2}\). Thus
\[
     \delta\ge 2c_{r,R_0}
          \left[\frac{q^2}{4\phi(0)}-\sqrt{\delta/2}
                         -\tau_{R_0}\right]_+^2.              \tag{CM6}
\]
This forces positive \(\delta\) when \(q>0\). If \(q\ge\eta>0\),
take \(R_0\) with \(\tau_{R_0}\le\eta^2/(16\phi(0))\).
Equation (CM6) gives
\[
 \delta\ge \frac{\eta^4}{\phi(0)^2}
                    \min\!\left(\frac1{128},\frac{c_{r,R_0}}{32}\right).
\]
The bound applies to either distinguished type involving \(k\).

In particular an actual maximizer cannot have the two labels
other than \(k\) separated in the coordinate along a centroid
mirror of \(k\), even if the separating sets in that coordinate
are disconnected. Such separation would give no \(\{i,j\}\)
mismatch because reflection preserves the coordinate along
the mirror. More generally the exact finite-swap test is
(CM2) for the reflected orbits themselves. The statements
apply to balanced maximizers at a right contact and use
neither the simplex nor a score equality on the mirror.
