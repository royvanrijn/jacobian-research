# The orthogonal nullcone in the balanced two-pair Image problem

The geometric verification question in the [two-pair counterexample, Section 5](TWO_PAIR_IMAGE_MATHIEU_COUNTEREXAMPLE.md#5-propagation-to-every-balanced-degree-dgeq4)
asks for a precise group action and its implication for fixed markers.
Use the full orthogonal action on the Hopf quadric, as follows.

Let \(k\) have characteristic zero, let \(R=\xi_1z_1+\xi_2z_2\), and let
\(V_d=k[\xi_1,\xi_2,z_1,z_2]_{(d,d)}\) for \(d\ge0\). Put
\(q=t^2+2xy\) and \(A=k[x,y,t]/(q-1)\). Dehomogenizing at \(R=1\) sends
\((\xi_1z_2,2\xi_2z_1,\xi_1z_1-\xi_2z_2)\) to \((x,y,t)\) and identifies
\(V_d\) with the degree-at-most-\(d\) part \(A_{\le d}\). The inverse
homogenizes in \(R\). This is the filtration isomorphism of the
[rank-one trace proof, Section 5](RANK_ONE_TRACE_AUTOMORPHISM_INVARIANCE_ALL_IMAGE_PAIRS.md#5-primitive-original-lift-and-full-ring-power-descent-in-every-number-of-pairs).
Transport the pullback action of \(\mathrm O(q)\) on \(A_{\le d}\) to
\(V_d\). The geometric nullcone \(\mathcal N_d\) consists of forms whose
orbit closure over \(\overline{k}\) contains zero. Write \(\mathcal P_d\)
for the forms \(H\in V_d\) with \(\mathcal E_2(H^m)=0\) at every \(m\ge1\).

**Proposition (nullcone implies all-marker safety).** For every \(d\ge0\),
\(\mathcal N_d\subseteq\mathcal P_d\), and for every \(H\in\mathcal N_d\)
and fixed original polynomial marker \(Q\), \(\mathcal E_2(QH^m)=0\) for
all sufficiently large \(m\).

Over \(\overline{k}\), the \(\mathrm O(q)\)-orbit closure is a finite
union of translates of the \(\mathrm{SO}(q)\)-orbit closure. The group
elements fix zero, so zero already lies in the connected orbit closure.
For nonzero \(H\), affine Hilbert–Mumford for the connected reductive
group \(\mathrm{SO}(q)_{\overline{k}}\) gives a cocharacter taking \(H\)
to zero. Apply the weight argument of Lemma `lem:unstable` in
[*The Mathieu property for quadratic complements*](../../../results/mathieu-quadratic-complements/main.tex)
to the Hopf algebra. If the image \(h\) of \(H\) has least cocharacter
weight \(\delta>0\) and a fixed nonzero regular marker \(a\) has least
weight \(b\), every weight in \(ah^m\) is positive once
\(m\delta+b>0\). The [sphere trace](TWO_PAIR_IMAGE_CUBIC_STRATIFICATION.md#the-full-coefficient-space-and-its-spherical-functional)
is orthogonally invariant, so it vanishes on every nonzero weight space.
Taking \(a=1\) gives purity at every positive order. The normalization
and marker balancing in Sections 3 and 5 of the rank-one trace proof
give eventual vanishing for every original polynomial marker. These
equalities descend to \(k\). The zero cases hold directly.

**Corollary (failure of the moment–nullcone equality).** For every
\(d\ge4\), \(\mathcal N_d\subsetneq\mathcal P_d\).

Equation (5.3) of the two-pair counterexample gives a pure balanced form
\(F_d=R^{d-4}F\) in every \(d\ge4\), with its contraction against the
fixed marker \(\xi_1z_2\) nonzero at every positive order. The proposition
therefore excludes \(F_d\) from \(\mathcal N_d\).
