# The one-operator Image theorem

## Theorem (the image of a constant-shift derivation)

Let \(A\) be a commutative \(\mathbb Q\)-algebra, let \(\delta\) be a
locally nilpotent derivation on \(A\), and let \(a\in\ker\delta\) be a
non-zero-divisor in \(A\). The image \(I\) of \(D=\delta-a\) is a
Mathieu subspace. For integers \(d\geq0\) and \(N\geq1\) there is a
computable \(q=q(d,N)\geq N\) with the following property.
If \(\delta^{d+1} f=0\) and
\[
 f^j\in I\qquad (N\leq j\leq q),
\]
then for each \(g\) with \(\delta^{e+1}g=0\), where \(e\geq0\),
\[
 gf^m\in I\qquad \bigl(m\geq q(d(q-1)+e+1)\bigr).
\]
For \(D=c\partial_z-a\) on \(B[z]\), let \(a,c\) be arbitrary elements
of a commutative \(\mathbb Q\)-algebra \(B\). With \(d,e\) as degree
bounds in \(z\), the same finite window implies \(gf^m\in I\) for all
sufficiently large \(m\). The displayed cutoff holds if \(a\) is a
non-zero-divisor or if \(c=1\). In particular the image is Mathieu.

This includes the one-dimensional regular-sequence Image assertion.
Theorem 2.8 of van den Essen, Wright and Zhao proves it for a
non-zero-divisor generating a radical ideal. Their Proposition 2.10
gives the polynomial case of the image criterion, and their Theorem 4.9
supplies the univariate factorial vanishing theorem. Van den Essen and
Zhao's Theorem 7.1 gives the image result for \(\partial_z-a\) over a
unique factorization domain. The polynomial result here concerns the
constant-coefficient case of their Conjecture 6.1.
[On the Image Conjecture, §§2 and 4](https://arxiv.org/html/1008.3962v2#S2) ·
[Mathieu Subspaces of Univariate Polynomial Algebras, §7](https://arxiv.org/pdf/1012.2017#page=17)

## Recognizing the image after inverting \(a\)

Set \(C=A[a^{-1}]\), with \(A\subseteq C\). The derivation extends to
\(C\) with \(\delta(h/a^r)=\delta(h)/a^r\) and is still locally nilpotent.
The finite exponential
\[
 \Psi(h)=\sum_{i\geq0}\frac{\delta^i h}{i!a^i}t^i\in C[t]
\]
is an algebra homomorphism because the Leibniz rule for \(\delta^i\)
gives the coefficient of \(t^i\) in a product. Let \(L\) be the
\(C\)-linear functional with \(L(t^i)=i!\). On \(C\) the inverse of
\(D=\delta-a\) is \(-a^{-1}\sum_{i\geq0}a^{-i}\delta^i\), a
geometric sum that terminates on each element. Since
\(L(\Psi(h))=-aD^{-1}h\), for every \(h\in A\) we obtain
\[
 h\in I\quad\Longleftrightarrow\quad L(\Psi(h))\in aA.                \tag{1}
\]
In the reverse direction the inverse belongs to \(A\), and its image
equals \(h\) in \(A\) by the injection into \(C\). For
\(\delta=c\partial_z\) on a polynomial ring with constant coefficients
\(a,c\), \(\Psi(h)=h(z+(c/a)t)\) by the polynomial Taylor formula.

## Finite homogeneous moment identities

Let \(S=\mathbb Q[X_0,\ldots,X_d]\), \(\mathfrak m=(X_0,\ldots,X_d)\)
and define the homogeneous degree-\(j\) polynomial
\[
 M_j(X)=L\left(\left(\sum_{i=0}^{d} X_i t^i\right)^j\right),
 \qquad L(t^i)=i!.
\]
The polynomials \(M_j\) for \(j\geq N\) have common zero set
\(\{0\}\) over \(\mathbb C\). At such a zero the corresponding
polynomial \(P(t)\) satisfies \(L((P^N)^s)=0\) for all \(s\geq1\).
The univariate factorial theorem gives \(P^N=0\), hence \(P=0\).
[On the Image Conjecture, Conjecture 4.2 and Theorem 4.9](https://arxiv.org/html/1008.3962v2#S4)

By the Nullstellensatz the radical of the ideal of these moments in
\(\mathbb C[X]\) is \(\mathfrak m\mathbb C[X]\). Choose exponents
\(b_i\geq1\) with \(X_i^{b_i}\) in that ideal. For \(q\geq N\) with
\(q>\sum_i(b_i-1)\), every monomial of total degree \(q\) is divisible
by one of these powers and belongs to the ideal. Its degree-\(q\)
part is the span of the polynomials \(X^\beta M_j\) with
\(N\leq j\leq q\) and \(|\beta|=q-j\). These have rational
coefficients, so the corresponding rational matrix has full rank as
well. Consequently we may find \(q\) by testing that rank for
\(q=N,N+1,\ldots\). The argument proves termination and gives
\[
 \mathfrak m^q\subseteq (M_N,\ldots,M_q)S.                              \tag{2}
\]

Moreover the monomials of degree less than \(q\) generate \(S\) as a
module over \(\mathbb Q[M_N,\ldots,M_q]\). To prove this by induction,
take any monomial of degree \(s\geq q\), factor a degree-\(q\) monomial
out of it, and use the degree-\(q\) identities above. The resulting
terms have the form \(M_j\) times monomials of degree \(s-j<s\).
The induction starts with the monomials of degree less than \(q\).

## Descending marked powers to the original ring

Write \(p(t)=\Psi(f)=\sum_{i=0}^{d}p_i t^i\) in \(C[t]\), and put
\(B=A[p_0,\ldots,p_d]\subseteq C\). By (1) and the hypotheses
\(M_j(p_0,\ldots,p_d)\in aA\) for \(N\leq j\leq q\). Evaluating the
graded module argument at the \(p_i\) shows that \(B\) is spanned over
\(A\) by monomials in the \(p_i\) of degree less than \(q\). Since
\(a^i p_i\in A\), these generators give
\[
 a^{d(q-1)}B\subseteq A.                                                 \tag{3}
\]
For the ideal \(J=(p_0,\ldots,p_d)B\), (2) gives \(J^q\subseteq aB\).
Hence \(J^{qr}\subseteq a^r B\) for every \(r\geq0\).

For any \(R\geq1\) take \(r=d(q-1)+e+R\) and \(m\geq qr\).
Every coefficient of \(p^m\) belongs to \(J^m\subseteq a^rB\).
Also \(a^e\Psi(g)\in A[t]\). By (3) the coefficients of
\(\Psi(gf^m)=\Psi(g)p^m\) therefore belong to
\(a^{r-e}B\subseteq a^RA\). For \(R=1\) their factorial moment lies in
\(aA\), so (1) proves \(gf^m\in I\). Taking \(N=1\) gives the
Mathieu property.

## Removing the coefficients killed by powers of \(a\) when \(c=1\)

For \(D=\partial_z-a\) on \(B[z]\) let
\(T=\{b\in B:a^k b=0\text{ for some }k\geq1\}\), an ideal.
Write \(H\) for integration with constant term zero. Given \(h\in T[z]\),
take \(k\) such that \(a^k h=0\). Then
\[
 D\left(\sum_{j=0}^{k-1} a^j H^{j+1}h\right)
   =h-a^k H^k h=h.
\]
Thus \(T[z]\subseteq I\). Under the quotient to \((B/T)[z]\), the
inverse image of the image of \(\partial_z-\bar a\) is exactly \(I\).
A polynomial preimage over \(B/T\) lifts to a polynomial over \(B\),
with an error in \(T[z]\). If \(\bar a\bar b=0\), some power of \(a\)
kills \(ab\), so \(b\in T\). The assertion is immediate if \(B/T=0\).
Otherwise \(\bar a\) is a non-zero-divisor, and the quotient preserves
the degree upper bounds and the finite list of image hypotheses.
The result with the same cutoff follows from the regular case.

## Two arbitrary constant coefficients

Fix \(f,g\) and the polynomial preimages witnessing the finite list
\(f^j\in I\) for \(N\leq j\leq q\). Replace \(B\) by the
\(\mathbb Q\)-subalgebra generated by \(a,c\) and the coefficients of
these polynomials and preimages. This preserves every hypothesis used
and the resulting preimages lift by inclusion. The subalgebra is
Noetherian, so the ascending chain of annihilators of \(a^j\)
stabilizes. For the ideal \(T\) of coefficients killed by some power of
\(a\), choose \(k\geq1\) such that \(a^kT=0\). Then
\(T\cap a^kB=0\), because if \(a^k b\in T\) then \(b\in T\)
and \(a^k b=0\).

Use \(M=2kN\). Since every coefficient of \(D(v)\) belongs to
\((a,c)B\), the hypothesis on \(f^N\) gives
\(f^M\in(a,c)^{2k}B[z]\). In particular \(f^M T[z]\subseteq c^kT[z]\)
by expanding the ideal power and using \(a^kT=0\).
Integration \(H\) with constant term zero satisfies \(DH=c-aH\), so
\[
 D\left(H\sum_{j=0}^{k-1}c^{k-1-j}(aH)^j h\right)=c^k h
 \qquad (h\in T[z]).
\]
It follows that \(f^M T[z]\subseteq I\).

Over \(B/T\), either the ring is zero or \(\bar a\) is a
non-zero-divisor. For
\[
 m\geq
 \max\{q(d(q-1)+e+k+1),\ M+q(d(q-1)+e+k)\},
\]
apply the coefficient bound for \(\Psi\) in the regular case at \(R=k+1\)
to \(h=gf^m\). By the inverse formula, its reduction has a preimage
divisible by \(\bar a^k\), since the factorial moment lies in
\(\bar a^{k+1}(B/T)[z]\). Lift such a preimage to
\(w\in a^k B[z]\). Then \(h-Dw\in T[z]\).
Apply the coefficient bound at \(R=k\) to \(gf^{m-M}\). Its reduction
lies in \(\bar a^k(B/T)[z]\), as seen from the constant term in \(t\)
of \(\Psi(gf^{m-M})\). Write \(gf^{m-M}=a^k v+h_0\) with
\(h_0\in T[z]\). If \(B/T=0\), both lifting steps hold with \(w=v=0\).
The polynomial \(h-Dw-f^M h_0\) belongs to \(T[z]\) and to
\(a^kB[z]\), and is therefore zero. Since \(f^M h_0\in I\), this gives
\(h\in I\) and completes the proof.

## An example over a nonnormal ring

The intermediate ring in the regular case can be strictly larger than
the original ring. For example, take the polynomial ring
\(B[z]\) with \(B=\mathbb Q[u^2,u^3]\), and let \(a=u^2\) and
\(f=u^3z\). Then \(\Psi(f)=u^3z+ut\) has a coefficient outside
\(B[z]\). For \(\partial_z-a\) the localized inverse of \(D\) applied
to \(h\in B[z]\) has constant term \(-a^{-1}L(h(t/a))\).
Coefficient recursion \((i+1)u_{i+1}=h_i+a u_i\) shows that
\(h\in I\) exactly when this constant term lies in \(B\).
For \(f^m\) the moment \(L(f^m(t/a))=m!u^m\) belongs to \(aB\)
exactly for \(m=2\) and \(m\geq4\). For the marker \(z^k\) the moment
is \((m+k)!u^{m-2k}\) and belongs to \(aB\) exactly when \(m-2k=2\)
or \(m-2k\geq4\). The sharp tail cutoff for this marker and phase is
\(2k+4\).

## Corollary (exact operator count for the regular-sequence assertion)

For each positive integer \(n\), consider the assertion that for every
commutative \(\mathbb Q\)-algebra \(A\) and every regular sequence
\(a_1,\ldots,a_n\) in \(A\), the sum of the images of
\(\partial_{z_i}-a_i\) is a Mathieu subspace of
\(A[z_1,\ldots,z_n]\). This assertion holds exactly when \(n=1\).
The same exact count applies to the Special Image assertion over every
characteristic-zero field.

The theorem gives the positive case. For each \(n\geq2\), take
\(A=K[\xi_1,\ldots,\xi_n]\) over a characteristic-zero field \(K\),
and use the regular sequence \(a_i=\xi_i\). The retained two-pair quartic
\(F\) and its marker \(Z\), viewed as polynomials in the first two
pairs, have the same contractions in any larger ring because all
their exponents in the extra pairs are zero. The image-kernel identity
then gives \(F^m\) in the summed image and \(ZF^m\) outside it for every
\(m\geq1\).
[Two-pair witness, Theorem 1.1 and the extension to additional pairs](TWO_PAIR_IMAGE_MATHIEU_COUNTEREXAMPLE.md#1-statement)
