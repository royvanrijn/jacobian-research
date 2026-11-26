# A finite free moment scheme from a proved support deletion

Let \(m_0,\ldots,m_r\) be distinct monomials in any finite
number of variables, with \(r\ge1\), and let \(\Lambda\) be
the usual factorial functional. Assume the family spanned
by \(m_1,\ldots,m_r\) has no nonzero phase with \(r\)
consecutive zero moments. For example this follows from
the Strong Factorial window for every support in that
family. The assumption is only about the family with
\(m_0\) deleted.

Fix a start \(n\ge1\), and write
\[
 P_j(W,Z)=\Lambda\left((Wm_0+\textstyle\sum_{i=1}^r Z_i m_i)^j\right),
 \qquad F_j(Z)=P_j(1,Z),\qquad H_j(Z)=P_j(0,Z).
\]
All have integer coefficients. Put
\[
 R_n=\operatorname{Res}_Z(H_n,\ldots,H_{n+r-1})\in\mathbb Z,
 \qquad d_n=\prod_{j=n}^{n+r-1}j .
\]

**Theorem.** The integer \(R_n\) is nonzero. Over
\(A=\mathbb Z[1/R_n]\) the affine algebra
\[
                C_n=A[Z_1,\ldots,Z_r]/(F_n,\ldots,F_{n+r-1})
                                                               \tag{1}
\]
is free of rank \(d_n\). In particular all solutions of
these \(r\) zero moment equations in the normalized chart
are algebraic and integral at every prime outside \(R_n\).
Using the unchanged moment polynomial \(F_{n+r}\), the
specific length-\((r+1)\) window is equivalent to
\[
   \det\nolimits_A\bigl(\text{multiplication by }F_{n+r}
                                    \text{ on }C_n\bigr)\ne0.   \tag{2}
\]
For each residue prime outside \(R_n\), invertibility
of that same multiplication map modulo the prime is
an exact sufficient certificate for (2).

If \(m_0=1\), every prime dividing any of the \(r+1\)
indices of an actual zero window divides \(\gcd(R_n,R_{n+1})\).
Use both adjacent length-\(r\) subwindows, since the normalized
coefficient tuple would belong to both algebras.
This restricts the denominator-prime budget using only
the deleted family, before elimination of the inhomogeneous
system. No bound on the size or height of the remaining
coefficients is assumed.

## Finite freeness and the determinant

The leading homogeneous forms of the \(F_j\) are the \(H_j\).
They have no projective common zero by the deletion
hypothesis, so their homogeneous resultant is nonzero.
After inverting it the same projective nonvanishing holds
in every residue field. Fix \(e_j=j\) for
\(n\le j\le n+r-1\), and put \(d=1+\sum_j(e_j-1)\).
Over every residue field and over \(\mathbb Q\), the
forms \(H_j\) are a homogeneous system of parameters in
a polynomial ring in \(r\) variables, hence a regular
sequence by Cohen–Macaulayness.
[Stacks Project, Lemma 10.104.2](https://stacks.math.columbia.edu/tag/02JN) Successively taking the
quotient by each homogeneous nonzerodivisor gives the
Hilbert series \(\prod_j(1-t^{e_j})/(1-t)^r\).
In particular the degree-\(d\) homogeneous quotient is zero.
It follows that the degree-\(d\) map over \(A\)
\[
          \bigoplus_j A[Z]_{d-e_j}\longrightarrow A[Z]_d,
       \qquad (V_j)_j\longmapsto\textstyle\sum_j V_j H_j
\]
is surjective. Its cokernel is finite and has zero
fiber at every maximal ideal, by the degree calculation.

Replacing \(H_j\) by \(F_j\) in a representation of a
degree-\(d\) monomial makes its residue in (1) a sum
of terms of smaller degree. Multiplication by further
monomials and induction reduce every degree to below
\(d\). Thus (1) is finite as an \(A\)-module.
The projective homogenizations \(P_j\) have no zero at
infinity in any fiber. Each projective intersection is
zero-dimensional, since a positive-dimensional projective
component would meet the infinity hyperplane. Its length
is \(d_n\) by the regular sequence Hilbert series in
\(r+1\) variables, whose eventual value is \(\prod e_j\).
This is the same affine fiber length. Since \(A\) is a
localization of \(\mathbb Z\), the structure theorem for
finite modules over a PID implies freeness: torsion at
a prime would enlarge the fiber dimension there beyond the
generic rank \(d_n\).

The fixed chart loses no possible common zero of the full
length-\((r+1)\) system, since \(W=0\) would contradict
the deletion hypothesis. In the finite algebra \(C_n\)
over an algebraically closed field of characteristic zero,
an element is a unit exactly when it vanishes in no
residue field of a maximal ideal. It is a unit exactly
when its multiplication map is invertible. This gives (2)
and its residue-prime criterion; no separability or
multiplicity-one property is required. The equality with the
homogeneous resultant is
\[
 \operatorname{Res}_{Z,W}(P_n,\ldots,P_{n+r})
     =R_n^{n+r}\,
        \det\nolimits_A(\text{multiplication by }F_{n+r}
                                      \text{ on }C_n).
\]
Use the variable order \(Z_1,\ldots,Z_r,W\), normalized by
the resultant of the corresponding variable powers.
The homogeneous Poisson formula compares the last polynomial
with \(W^{n+r}\); restriction to \(W=0\) then gives the factor
\(R_n^{n+r}\). The determinant is the product of the values
at \(W=1\) with their local algebra lengths.
[Busé–D’Andrea–Sombra–Weimann, Propositions 2.1 and 2.2](https://www.maia.ub.edu/~sombra/publications/Flexes/Flexes.pdf)

Integrality of the
affine coordinates also follows from finiteness over \(A\).
For constant \(m_0=1\), the
[Frobenius argument for a multiple of a good prime](FACTORIAL_CONSTANT_TERM_DENOMINATOR_SMOOTH_WINDOWS.md#a-prime-dividing-the-moment-index)
then excludes zero at any index divisible by \(p\nmid R_n\).

## Application to even binary quadratics

Take \(m_0=1\), with \(m_1,m_2,m_3=x^2,xy,y^2\).
The [homogeneous quadratic three-window](BINARY_HOMOGENEOUS_QUADRATIC_STRONG_FACTORIAL.md)
gives the deletion hypothesis. Thus the first three
consecutive moments, at any prescribed start \(n\),
give a finite algebra of rank \(n(n+1)(n+2)\). This
includes the sparse coefficient loci. Exclusion of a
fourth moment is one explicit determinant on this
algebra.

In this case \(R_n\) has no prime factor greater than
\(2n+3\). To verify this, let \(p>2n+3\), so \(p>2n+4\).
Over the residue field normalize the three homogeneous
moments \(H_j\) by \(H_j/(2j)!\), \(j=n,n+1,n+2\),
and denote these forms by \(J_j\). The weighted
[homogeneous recurrence](BINARY_HOMOGENEOUS_QUADRATIC_STRONG_FACTORIAL.md#1-the-interval-identity)
holds after reduction, with the inward terms
\((c-b/2)c^j+(a-b/2)a^j\) and derivative term
\(2j(ac-b^2/4)J_{j-1}/(2j-1)\).
In particular the two inward terms at \(n+1,n+2\)
vanish at any common zero. When \(a,c\) are distinct
and nonzero, their geometric determinant is nonzero,
forcing \(b=2a=2c\), impossible. If \(a=c\ne0\),
the inward terms force \(b=2a\); the phase is then
\(a(x+y)^2\) and \(H_n=(2n+1)!a^n\) is a unit up to
\(a^n\). If \(c=0\ne a\), one instead obtains \(b=2a\).
On that plane \(H_n=4^n(n!)^2 a^n\) by the beta
identity for \(a(x^2+2xy)\). The exchanged case is
identical. If \(a=c=0\), the phase is \(bxy\), whose
moment is \(b^n(n!)^2\). All these formulas are valid
in residue characteristic \(p>2n+3\). This rules out
a common projective zero and proves the prime bound.
At smaller primes the resultant need not be a unit,
so the prime bound alone does not resolve the unrestricted
four-window.
