# Fibres and monodromy of cotangent assignments

## Polynomial product structure

Let \(k\) have characteristic zero and \(G:\mathbb A^n_k\to\mathbb A^n_k\)
be a polynomial Keller map. For \(M(x)\in\operatorname{GL}_n(k[x])\)
and \(b(x)\in k[x]^n\), put
\[
 T(x,v)=(G(x),M(x)v+b(x)).
\]
The source change \(\Psi(x,v)=(x,M(x)v+b(x))\) is a polynomial
automorphism, with inverse \((x,t)\mapsto(x,M(x)^{-1}(t-b(x)))\).
Consequently \(T=(G\times\mathrm{id}_{\mathbb A^n})\circ\Psi\).

For any field extension \(E/k\), the scheme fibre of \(T_E\) over
\((y,t)\in\mathbb A^{2n}(E)\) is isomorphic to the scheme fibre of
\(G_E\) over \(y\). The inverse fibre map gives each base point \(x\)
the unique coordinate \(v=M(x)^{-1}(t-b(x))\), an identity in the
fibre algebra. Thus empty fibres, full-degree finite étale fibres,
and their residue field degrees transfer in both directions.

The same formula identifies the entire maps after a source change.
The maximal target open where the map is finite is the base open
times \(\mathbb A^n\). For one inclusion, use stability of finiteness
under base change. For the reverse inclusion on geometric points,
if the product map is finite over an open containing \((y,t)\), pull
back along the section with second coordinate \(t\) to obtain a
neighborhood of \(y\) over which \(G\) is finite. The finite normalization
also commutes with this product. If \(B\) is the integral closure of
\(k[y]\) in \(k(x)\) via \(y=G(x)\), then \(B[t]\) is normal and finite
over \(k[y,t]\), with the required fraction field \(k(x)(t)\).
It is therefore the integral closure for the product. Over
\(\mathbb C\) the nonproperness locus of the quasi-finite map is thus
the planar locus times the added affine factor.

Writing \(L=k(x)\) and \(K=k(G(x))\), the generic extension is
\(L(t)/K(t)\) with \(t\) independent. If \(N/K\) is the Galois closure of
\(L/K\), then \(N(t)/K(t)\) has the same Galois group and the same
point stabilizer. Indeed a \(K\)-basis remains a \(K(t)\)-basis, and the
\(K(t)\)-conjugates of \(L(t)\) generate \(N(t)\). Applying the argument
over \(\overline k\) also preserves geometric monodromy with its
permutation action. The induced generic degree is \([L:K]\).

## Hessian gradients and Weyl position assignments

Write \(\Phi(x,v)=R(x)+v^{\mathsf T}G(x)\). After swapping the two
output blocks, its gradient has the form above with
\(M=DG^{\mathsf T}\) and \(b=
abla R\). Conversely, taking \(R=0\)
gives such a gradient from each \(G\). The block Hessian determinant
is \((-1)^n(\det DG)^2\), so it is a nonzero constant.
For a complex four-variable potential, Tang's structural alternative
reduces every noninvertible gradient to this form with \(n=2\).
For the linear source change \(x=Lu\), the identity
\(
abla(\Phi\circ L)(u)=L^{\mathsf T}(
abla\Phi)(Lu)\) gives the
corresponding invertible target change. Under Tang, the possible
generic degrees and geometric permutation actions of noninvertible
four-variable gradients therefore coincide exactly with those of
noninvertible plane Keller maps. The direction from a planar map to
its doubled gradient uses the displayed determinant and product
structure alone.
[Tang, Theorem 1.1 and Lemma 1.2](https://hexagonmath.org/pdf/2611.00156v2).

For a finite Weyl endomorphism sending each \(x_i\) to \(G_i(x)\),
the [positional classification](../../dixmier/dixmier-formal-lifting/IDENTITY_LIFTS_AND_UNBOUNDED_INVERSES.md#converse-for-assignments-preserving-the-multiplication-algebra)
gives momentum images \(\sum_j ((DG)^{-1})_{ji}p_j+b_i(x)\), where
\(b_i\) has momentum
order zero. The polynomial map obtained by making these ordered
expressions commutative has \(M=(DG)^{-\mathsf T}\). Its fibres and
generic monodromy are therefore those of the positional map \(G\)
in the sense just proved. The Weyl endomorphism is surjective exactly
when \(G\) is polynomially invertible. Conversely every Keller map
gives an exact finite Weyl assignment by taking \(b=0\).
In rank two this yields the same bidirectional correspondence of
degrees and permutation actions with plane Keller maps, without a
Hessian structural assumption.

## Low generic degrees and the \(A_4\) exclusion

A whole-affine Keller map of generic degree one is an automorphism.
For each source coordinate write its rational inverse as a reduced
fraction \(a(y)/b(y)\). If \(b\) has a prime divisor downstairs,
[the affine divisor pullback lemma](THEORY.md#3-every-affine-divisor-is-met-without-ramification)
supplies a source divisor dominating it. The fraction has a pole there
since its numerator is a unit at the generic point downstairs.
Its pullback is a regular source coordinate, a contradiction. Each
inverse coordinate is thus polynomial, and the generic inverse
identities hold everywhere.

Orevkov's Theorem 1.1, §1, excludes generic degrees two and three
for complex polynomial plane Keller maps. Its setup defines
multiplicity as the cardinality of a generic fibre. Sections 3–5
apply the compactification Euler characteristic and the topology
of the resulting branch curve in degree three. The complete
four-sheet exclusion is the separately retained
[Domrina theorem](PLANE_A4_MONODROMY_EXCLUSION.md#the-unconditional-four-sheet-theorem).
The same note's finite-coefficient-field descent transfers both
published exclusions to every characteristic-zero field. These inputs
therefore give the lower bound of five on the generic degree of any
noninvertible plane Keller map. Through the product descriptions
they give the same bound for the polynomial map associated to a
nonsurjective positional rank-two Weyl assignment, and under Tang
for a noninvertible complex four-variable gradient. A proof of
inversion through any given generic degree on either recipient
also transfers back to the plane by the exact cotangent maps.
[Orevkov, Theorem 1.1](https://www.math.univ-toulouse.fr/~orevkov/jc86.pdf).

In both of the reverse Hessian transfers, first apply the doubling over
the complex numbers. The finite-coefficient-field descent then extends the
plane exclusion to characteristic zero.

Likewise the [plane \(A_4\) exclusion](PLANE_A4_MONODROMY_EXCLUSION.md)
rules out every faithful transitive generic geometric \(A_4\) action
in the Weyl positional class in rank two, and under Tang for complex
four-variable gradients. The correspondence concerns generic
permutation actions. Special-fibre groups and polynomial component
degrees have their own invariants.
