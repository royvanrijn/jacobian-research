# Jacobian and generalized Dixmier transfer through plane graph modules

Write \(A=A_1(\mathbb C)\) for the first Weyl algebra. In
[Kanel-Belov–Fedorov v2](sources/kanel-belov-fedorov-2610.33919v2.pdf),
Conjecture 1.2, the assertion \(DC'_1\) says that every unital
\(\mathbb C\)-algebra map \(A\to D\) to a complex domain Morita equivalent
to \(A\) is an isomorphism. Theorem 2.3 identifies this assertion with
invertibility of all endomorphisms of \(\operatorname{Mat}_2(A)\).
In particular \(DC'_1\) implies ordinary \(DC_1\) by taking \(D=A\).

## From the plane Jacobian conjecture

**Theorem PlaneJacobianToGeneralizedDixmierTransfer.** The complex plane
Jacobian conjecture \(JC_2\) implies \(DC'_1\).

This is Proposition 4.3, printed pages 23–24. For an endomorphism of
\(\operatorname{Mat}_2(A)\), the authors choose a smooth arithmetic
model after one localization. In positive characteristic the reductions
act on the Weyl center and preserve its commutator bracket, so the
induced plane maps have Jacobian one. If the original endomorphism has
generator-image degrees bounded by \(d\) in the Bernstein filtration,
the center maps have ordinary degree at most \(d\). Under \(JC_2\), the
bounded-degree inverse test in §4.1 makes all such center maps invertible
in sufficiently large characteristic. Their inverse degrees are at most
\(d\). Invertibility of the center map implies invertibility of the
Azumaya algebra map. The matching inverse degrees give finite module
generation for the twisted arithmetic action by Lemma 3.5. The
endomorphism's cokernel has zero generic fiber by the finite-module
specialization argument in §3.3. This proves surjectivity in
characteristic zero, and simplicity gives injectivity.

## From bounded faithful families

**Theorem FaithfulBoundedGraphsGeneralizedDixmierTransfer.** Assume
\(DC'_1\) and suppose each complex plane Keller map of Jacobian one
admits a faithful arithmetic family in which the Berest–Wilson invariants
of its actual splitting graph modules are bounded on an ultrafilter-large
set. Then \(JC_2\) holds.

The hypothesis is the one used in the proof of Corollary 1.10, §5.4,
printed page 31. For a given map \(f\), choose a finitely generated
integral \(\mathbb Z\)-algebra \(R_0\subset\mathbb C\) of definition,
reductions \(s_i:R_0\to k_i\) to finite fields with
\(p_i=\operatorname{char}k_i\to\infty\), and a nonprincipal ultrafilter
\(\mathcal U\) such that \(R_0\to\prod_{\mathcal U} k_i\) is injective.
The reductions have \(p_i>\deg f\) after discarding finitely many
indices. The prescribed modules are \(M_{f_{k_i},k_i}\), including
their right actions and coefficientwise Frobenius twists on the centers.

Definition 1.7 and Theorem 1.8 give an arithmetic realization of those
same modules along \(\mathcal U\). Its integral finite-type coefficient
ring \(R\) lies over \(R_0\) inside \(\prod_{\mathcal U} k_i\). The
specialization maps \(t_i:R\to k_i\) extend \(s_i\) and represent each
element of \(R\) in the ultraproduct. The characteristic-zero bimodule
over \(\operatorname{Frac}R\) can be chosen projective of rank one on the
left. After embedding \(\operatorname{Frac}R\) into \(\mathbb C\) over
\(R_0\), this is a left generator by simplicity of \(A\), and its right
action gives a map to the opposite endomorphism domain, which is Morita
equivalent to \(A\). By \(DC'_1\) that map is an isomorphism, making the
bimodule invertible. Faithfully flat descent and Lemma 4.1 spread
invertibility after localization. Every nonzero inverted coefficient
remains nonzero at \(\mathcal U\)-almost all indices. On them the reductions of \(f\) are
automorphisms, because their Frobenius-twisted maps are the center
isomorphisms of the preserved splitting bimodules. The finite inverse
identities (22) give invertibility of \(f\) in characteristic zero.

**Corollary UniformPlaneGraphBoundSuppliesFaithfulFamilies.** The
uniform all-prime bound of v2 Conjecture 1.9 yields the faithful
bounded-family hypothesis above.

Indeed, for each \(f\), enumerate the nonzero elements of a finitely
generated model \(R_0\). At the \(i\)-th stage one can pick a finite
closed point of characteristic greater than \(i\) and \(\deg f\)
avoiding the zero sets of the first \(i\) entries. The surviving open
subset meets the generic fiber, and its image in \(\operatorname{Spec}
\mathbb Z\) contains a nonempty open set by Chevalley's theorem.
Any nonprincipal ultrafilter gives faithfulness on these reductions.
The uniform bound supplies boundedness on a set in the ultrafilter.
The [two-degree fullness criterion](PLANE_KELLER_GRAPH_BW_FILTRATION.md#exact-two-degree-criterion)
gives the equivalent route from a common two-degree bound on the
actual graph modules. The faithful-family hypothesis only asks for
one such family per complex map, with a bound allowed to depend on
that map and family. The multiplicity-one lift with constant arithmetic
support on a dense open belongs to the separate implication in
Theorem 1.6(ii), which assumes Conjecture 1.5.
