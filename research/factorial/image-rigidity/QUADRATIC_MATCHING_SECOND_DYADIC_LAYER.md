# The second dyadic matching layer for arbitrary quadratic supports

Consider a zero-constant quadratic
\[
 f=\sum_{i=1}^m(\alpha_i X_i^2/2+d_i X_i)
                        +\sum_{i<j}\beta_{ij}X_iX_j,\qquad
 E(T)=\sum_{n\ge0}\Lambda_m(f^n)T^n/n!,
\]
with natural coefficients integral at a place above two.
Write \(\mathcal O\) for its valuation ring, with \(v(2)=1\).
The following formula gives \(E\) coefficientwise modulo
\(4\mathcal O\), for arbitrary graph, dimension and coefficient
support. It refines the
[first matching layer](QUADRATIC_STRONG_FACTORIAL_WINDOWS_FROM_MATCHINGS.md#2-the-universal-matching-generating-function)
even in ramified coefficient rings. Its new graph terms
are a matching with one local defect. The pole orders of
the correction are controlled by the least common multiple
of the vertex denominators.

## The single-vertex correction

At vertex \(i\), set
\[
 U_i=1-\alpha_iT-d_i^2T^2,\qquad F_i=(1+d_iT)/U_i,\qquad
 J_{r,i}(T)=\frac1{r!}\Lambda_1
                \big(X_i^r\exp(T(\alpha_i X_i^2/2+d_iX_i))\big).
\]
The constant terms of the denominators are one. For a
single vertex, suppress the subscript \(i\). The
[univariate factorial recurrence](QUADRATIC_STRONG_FACTORIAL_WINDOWS_FROM_MATCHINGS.md#2-the-universal-matching-generating-function)
gives the exact formal equation
\[
                      UJ_0=1+dT+2\alpha T^2 J_0'.
\]
Consequently \(J_0\equiv F+(2\alpha T^2/U)F'\pmod4\).
The exact identity \(J_1=(1+T\partial_T+\alpha\partial_\alpha)J_0\)
follows by weighting the \(\alpha^a d^b T^{a+b}\) coefficient
by \(1+2a+b\). This is a universal polynomial
differentiation before substituting the coefficients.
In particular it preserves the coefficientwise congruence.
The two vertex formulas simplify to
\[
 \begin{split}
 J_0&\equiv \frac{1+dT}{U}
              +2\frac{\alpha(\alpha+d+d^3T^2)T^2}{U^3},\\
 J_1&\equiv -\frac1U+2\frac{1+dT}{U^2}
                             +2\frac{\alpha^2T^2}{U^3}
                                     \pmod{4\mathcal O[[T]]}.
 \end{split}                                                 \tag{1}
\]
Indeed \(F'\equiv(\alpha+d+d^3T^2)/U^2\pmod2\).
Also \((1+T\partial_T+\alpha\partial_\alpha)F
=-1/U+2(1+dT)/U^2\) exactly.
On \(\alpha(\alpha+d+d^3T^2)T^2/U^3\) modulo two, the
denominator derivatives of \(T\partial_T+\alpha\partial_\alpha\)
cancel, while the \(\alpha d^3T^4\) contribution cancels
with its own numerator derivative. The combined numerator
is \(\alpha^2T^2\).
This proves (1).

## The graph terms of valuation one

For each array of cross counts \(j_e\ge0\), put
\(r_i=\sum_{e\ni i}j_e\), \(J=\sum_e j_e\), and
\[
                W(\mathbf j)=\frac{\prod_i r_i!}{\prod_e j_e!}.
\]
The exact contribution to \(E\) is
\[
 W(\mathbf j)\,T^J\prod_e\beta_e^{j_e}
                                   \prod_i J_{r_i,i}(T).    \tag{2}
\]
The multiplier \(W(\mathbf j)\) is an integer by the
vertex multinomial coefficients. If \(s_2\) is binary
digit sum, Legendre's formula gives
\[
 v_2(W)=\sum_e v_2(j_e!)
      +\sum_i\left(\sum_{e\ni i}s_2(j_e)-s_2(r_i)\right).    \tag{3}
\]
Each vertex term in parentheses is the total carry
count in adding its incident numbers in base two,
including carry multiplicities.

The patterns of value zero in (3) are exactly the
ordinary matchings. The value-one patterns are:
- Every count is zero or one, and there is precisely one
  vertex of degree two or three. The other nonzero
  edges, away from this star, form a matching.
- Exactly one edge has count two or three, and all
  other nonzero counts are one. The count-one edges
  form a matching. At a count-three edge neither
  endpoint may meet a count-one edge. At a count-two
  edge each endpoint may meet at most one count-one
  edge.

For the first type this follows from the factorial
valuation of the vertex degree. For the second type
the sole external factorial factor of \(2\) is already
accounted for. Therefore all vertex additions must
be carry-free. The count-two edge uses bit one and the
count-three edge uses bits zero and one, yielding the
stated rule. This classification uses integer weights,
so nonunit cross coefficients simply raise valuations
further. It is valid in ramified fields.

All other patterns vanish modulo \(4\mathcal O\).
In a value-one pattern every \(r_i\le3\).
By the
[vertex digit identity](BINARY_ZERO_CONSTANT_QUADRATIC_WEIGHTED_MATCHING_LOCAL.md#exact-cross-weights-below-the-unit-valuation),
one needs only
\[
        J_{r_i,i}\equiv
          \frac{(1+d_iT)^{1-(r_i\bmod2)}}{U_i^{\,1+\lfloor r_i/2\rfloor}}
                                                \pmod{2\mathcal O}.
                                                               \tag{4}
\]
The digit proof is a polynomial identity in characteristic
two, so (4) is valid in \(\mathcal O/2\mathcal O\),
not just after reduction modulo the maximal ideal.

## The rational correction and its denominator

Put \(D=\prod_i U_i\), and let
\[
 N=\sum_{M\ {\rm matching}} T^{|M|}
          \prod_{e\in M}\beta_e\prod_{i\ {\rm unmatched}}(1+d_iT).
\]
Then
\[
                 E\equiv N/D+2(R_{\rm mat}+R_{\rm graph})
                                               \pmod{4\mathcal O[[T]]}.
                                                               \tag{5}
\]
Here the correction terms on the right need only be
specified in \((\mathcal O/2\mathcal O)[[T]]\).
For \(R_{\rm mat}\), form the sum over matchings of
the product of the base vertex factors \(F_i,-1/U_i\)
for unmatched and matched vertices, with exactly one
vertex factor replaced by its coefficient of \(2\) in
(1). Retain \(T^{|M|}\prod_{e\in M}\beta_e\).
For \(R_{\rm graph}\), sum over the two value-one types
in (3), using (2) without \(W\) and using (4) for their
vertices. The omitted multiplier \(W/2\) is odd.
Equations (1)–(4) prove (5), coefficient by coefficient,
with no infinite summation exchange: every coefficient
uses only finitely many cross arrays.

Now make a second reduction, into the residue field
\(k=\mathcal O/\mathfrak m\). Put
\[
       L=\operatorname{lcm}_i\bar U_i
\]
with constant coefficient one. The residue denominator
of \(R=R_{\rm mat}+R_{\rm graph}\) divides
\[
                              \bar D L^2.                   \tag{6}
\]
In a corrected matching, only one vertex has up to
two extra denominator powers. A simple star defect has
only one vertex with one extra power. A repeated-edge
defect has at most two such vertices. Since
\(\bar U_i\bar U_j\mid L^2\) for every \(i,j\),
the displayed denominator works for every term.

In particular, suppose \(N\equiv QD\pmod{2\mathcal O[T]}\)
for an integral polynomial \(Q(T)\), and put
\(Z=(N-QD)/2\). At indices beyond \(\deg Q\), the residue
of \(P_n/2\) is the coefficient of the rational function
\(\bar Z/\bar D+\bar R\), with denominator dividing (6).
If its reduced denominator is nonconstant of degree \(r\)
and its polynomial part is \(I\), then
\[
                       r+\max(0,\deg Q,\deg I)
\]
is a sufficient window length at every positive start.
After discarding indices at which \(Q\) or \(I\) can contribute,
the proper fraction cannot have \(r\) consecutive zero
coefficients, by the
[backward recurrence](QUADRATIC_MATCHING_RESIDUE_DENOMINATOR_WINDOWS.md#proof-of-the-rational-series-criterion).
In particular this excludes an actual-support block when
the given length is at most the support size.

## Fewer extra pole powers without mixed vertex units

Assume \(\bar N/\bar D\) is a polynomial, and
\(\bar\alpha_i\bar d_i=0\) at every vertex. Then the
denominator of \(\bar R\) in (5) divides
\[
 \bar D\,\operatorname{lcm}
      \bigl(\{\bar U_i\}_i,\
                  \{\bar U_i\bar U_j:\bar\beta_{ij}\ne0\}\bigr).
                                                               \tag{7}
\]
Indeed, a corrected matching has a term of denominator
\(\bar D\bar U_i^2\) only when \(\bar\alpha_i\ne0\).
Then \(\bar d_i=0\), so the original matching factor
at that vertex has numerator one whether it is matched
or unmatched. The two high-denominator corrections in
(1) have the same numerator \(\bar\alpha_i^2 T^2\).
Summing all matchings at that vertex therefore makes
their high term
\[
                    \frac{\bar\alpha_i^2 T^2\bar N}
                                       {\bar D\bar U_i^2}.
\]
Since \(\bar N/\bar D\) is polynomial,
\(\bar D\bar U_i\) is already a denominator for this term.
The remaining matching and star terms add at most
one extra \(\bar U_i\). Repeated-edge terms add
\(\bar U_i\bar U_j\) at the endpoints of the edge,
and contribute only when that edge has unit residue.
This proves (7). In particular the extra denominator
has only one vertex power per pole when the nonconstant
\(\bar U_i\) are pairwise coprime.

## Restoring a pole for a quadratic star

There is an explicit second-layer test even when the
entire first residue is the polynomial \(1\).
Consider a star with unit diagonal at its center and
scale the phase so that
\[
                       f=X^2/2+dX+X\sum_i\beta_iY_i.
\]
Let \(h_j=h_j(\beta_1,\ldots,\beta_m)\) denote the
complete homogeneous coefficients,
\(\sum_j h_j T^j=\prod_i(1-\beta_iT)^{-1}\).
Assume \(d^2\in2\mathcal O\) and \(1+d-h_1\in2\mathcal O\).
Set \(\varepsilon=(1+d-h_1)/2,\ \tau=d^2/2\).
The residue of the exact coefficient series \((E-1)/2\) is
\[
  \frac{(\bar\varepsilon T+\bar\tau T^2)(1+T)
              +T+(\bar h_2+1)T^2+\bar h_3T^3}{(1+T)^2}.
                                                               \tag{7}
\]
To obtain it, group the edges by their total count \(j\).
The leaf factorials cancel the count factorials in (2), so
the exact expansion is \(\sum_j j!h_j T^jJ_j\) at the
center. Only \(j\le3\) remain modulo four. Use (1) for
\(j=0,1\), and (4) for \(j=2,3\).
After the base term \((1+(d-h_1)T)/U\) is subtracted,
the residue of the correction equals the numerator
\(T+(\bar h_2+1)T^2+\bar h_3T^3\) over \((1+T)^2\);
here \(\bar d=0\) and \(\bar h_1=1\).
Subtracting \(1\) from the base term contributes
\((\bar\varepsilon T+\bar\tau T^2)/(1+T)\).
This gives (7).

Its numerator and its derivative at \(T=1\) have values
\[
                    \bar h_2+\bar h_3,\qquad
                    \bar\varepsilon+\bar\tau+1+\bar h_3 .
\]
If the first is nonzero, the reduced denominator
has degree two and the numerator degree is at most
three, proving a three-window. If the first is zero
and the derivative is nonzero, one factor cancels and
the resulting numerator has degree at most two with
denominator of degree one, proving a two-window.
This applies also when some leaf coefficients vanish.

For example take \(d=0\) and two leaf coefficients
\(\beta_1=\zeta,\beta_2=\zeta^2\), where
\(\zeta^2+\zeta+1=0\).
Then \(h_1=-1,\ h_2=0,\ h_3=1\), so the first
layer is exactly \(1\) in the residue field and the
second layer has a degree-two pole. The three displayed
monomials satisfy their actual-support window.
The star test follows from the same coefficient identities
in every order.
