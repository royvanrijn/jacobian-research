# A uniquely exposed point of a compact simplex discrepancy

Let \(A\) be an actual balanced planar Gaussian noise-stability
maximizer at \(0<r<1\), agreeing with a centered regular simplex \(T\)
outside a compact set and differing from \(T\) on positive measure.
Write \(w_1,w_2,w_3\) for the unit sector centers of \(T\). Let \(D\)
be the closure of the open set of strict winners whose labels differ
from the reference label, with the reference rays removed from that
open set. Define
\[
     M=\max_{x\in D}\max_i w_i\cdot x,\qquad
     F_k=\{x\in D:w_k\cdot x=M\}.
\]
If \(F_k\) is a singleton \(\{p\}\), then \(p\) lies strictly inside
\(T_k\) and is a genuine four-arm junction. In normal and tangent
coordinates at the supporting line, the two boundary graphs of its
middle winner have the same even contact order \(m\ge2\) with that
line. Inward from the supporting line, their positive leading
coefficients have ratio \(2^m\).

Use the weighted slice identities for the compact discrepancy.
In particular on almost every line perpendicular to the outward
direction of a reference \(i,j\) ray, the Gaussian-weighted integral
of the change in \(1_{A_i}-1_{A_j}\) from the reference contrast is
zero. The Gaussian slice density is positive and smooth.
[Compact discrepancy slice identities, (CD1)](THREE_CELL_COMPACT_SIMPLEX_DISCREPANCY.md#every-label-must-receive-a-bounded-discrepancy)

## Location and linear geometry

The compact set \(D\) lies in \(\{x:w_i\cdot x\le M\ \forall i\}\).
Since \(D\) is nonempty and has interior, \(M>0\). Thus \(p\) is
in the closed reference sector \(T_k\), and uniqueness of its
exposed point confines every change on slices \(w_k\cdot x\) near
\(M\) to an arbitrarily small neighborhood of \(p\).

Suppose first that \(p\) is on the reference \(k,l\) ray and denote
the third label by \(j\). Along the ray just beyond \(p\) there is
no discrepancy, so it is an actual \(k,l\) interface there. The
regularity bound (G20) makes the pair gradient nonzero at \(p\).
The adjusted \(k,l\) score difference vanishes on the whole
supporting line by the exterior agreement and analyticity. Thus
near \(p\) the two sides have the reference strict order of these
two scores. The only possible changes there are \(k\to j\) and
\(l\to j\), both giving strictly negative change in the \(l,j\)
contrast. But the reference \(l,j\) ray points along \(-w_k\),
so (CD1) forbids that one sign on any slice close to the exposed
line. There are changes on a positive measure of slices arbitrarily
near the exposed line by the definition of \(D\). This contradicts
uniqueness of the exposed point on the ray.
[Regularity at finite arc accumulations, (G20)](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#regular-interface-arcs-and-isolated-triple-junctions)

Thus \(p\) lies inside \(T_k\). On slices sufficiently close to the
exposed line, all changes come from source label \(k\). If \(l,j\)
are the two other labels, (CD1) for their pair gives equal positive
weighted amounts of \(k\to l\) and \(k\to j\) on any such slice
with positive change. Hence both wrong strict winners accumulate
at \(p\), and the original label \(k\) wins in the outside local
half-plane. The point is a genuine triple junction. Its three
gradient differences must be collinear. Indeed the score of \(k\)
dominates each other score to first order in every direction into
the outside open half-plane. A nonzero linear form nonnegative
there must be a positive multiple of the outward normal form.
The gradient differences are nonzero by (G20) and the genuine
junction classification. In particular \(k\) has the largest
outward derivative, and one of the two wrong labels has the
smallest derivative while the other is strictly between.

## Matched first terms of the slice lengths

Choose local coordinates \((t,x)\) at \(p\), tangent and outward
normal to the supporting line, so the discrepancy lies in \(x\le0\).
Denote the middle-gradient wrong label by \(m_0\) and the low-gradient
wrong label by \(l_0\). There are three analytic pair roots
\[
  x=u(t)\quad(k,m_0),\qquad x=g(t)\quad(k,l_0),\qquad
                     x=v(t)\quad(m_0,l_0),
\]
all zero at \(t=0\) with horizontal tangents. The derivatives of the
three pair scores in increasing gradient order remain strictly
positive in the outward variable. Thus label \(k\) wins above
\(\max(u,g)\). If \(u>g\), the middle label wins between \(v\)
and \(u\), with \(v<g<u\), and the low label wins below \(v\).
If \(u<g\), the middle label is absent and the low label wins
below \(g\). These assertions follow by testing the ordered
contrasts at \(g\) and using their strict outward monotonicity.

Both \(u,g\) are nonpositive near zero, since otherwise the
discrepancy would reach the outside half-plane. Neither vanishes
identically, since a zero upper wall along an interval would
put that interval into \(D\) on the supporting line. Thus they
have negative leading coefficients of even orders. Also
\(u-g\) is not identically zero because the middle label wins
arbitrarily near \(p\). Its first term gives either a one-sided
middle strip or a strip on both sides. In the one-sided case its
order is odd. Both \(u,g\) then have the same earlier even leading
term \(-A t^m\), with \(A>0\), since different orders or leading
coefficients would give an even leading term for their difference.
The same leading term holds for \(v\). Indeed at \(x=g(t)\) the
two outer scores agree, and the strictly monotone transverse
derivatives give \(|v-g|\le C|u-g|\).

In this one-sided case the wrong labels together occupy length
asymptotic to \(2(\delta/A)^{1/m}\) on a slice \(x=-\delta\) as
\(\delta\downarrow0\). The middle strip has length
\(o(\delta^{1/m})\), since all its boundary roots share that
same leading term. To see these estimates, on each side of
\(t=0\) a nonpositive root with leading term \(-A t^m\) is
monotone away from zero in a sufficiently small neighborhood.
Its intersection with \(x=-\delta\) consequently occurs at
\(|t|=(\delta/A)^{1/m}(1+o(1))\). The slice Gaussian density
is its value at \(p\) times \(1+o(1)\) on all these intervals.
Equality of the two weighted wrong-label amounts is impossible.

In the two-sided case, \(v<g<u<0\) on both sides, so the two
middle-strip roots have leading terms \(u=-A t^m+o(t^m)\)
and \(v=-B t^n+o(t^n)\), with \(A,B>0\) and \(m,n\) even.
The wrong labels together have slice length asymptotic to
\(2(\delta/A)^{1/m}\), and the low label has length asymptotic
to \(2(\delta/B)^{1/n}\). The same smooth density factor applies
to both. Since the low label occupies exactly half the weighted
wrong-label length on almost every such slice, one must have
\(n=m\) and \(B=2^m A\). The middle winner persists on both
sides, giving four actual incident arcs.
