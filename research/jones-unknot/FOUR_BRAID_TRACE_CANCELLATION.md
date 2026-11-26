# Jones unknot detection in a four-braid twist family

The Jones Unknot Conjecture asks whether \(V_K(t)=1\) implies that \(K\) is
the unknot. There is no counterexample in the family

\[
\widehat{\sigma_1^a\sigma_2^b\sigma_3^c
          \sigma_2^d\sigma_1^e\sigma_3^f\sigma_2^{2k}},
\qquad
a,b,c\in\{-7,-5,-3,-1,1,3,5,7\},\quad
d,e,f\in\{-6,-4,-2,0,2,4,6\},\quad k\in\mathbb Z.
\]

The proof combines an exact trace reduction in the unrestricted exponent \(k\),
an integer Laurent-polynomial computation, and the published 24-crossing
theorem. A separate exact computation rules out all nonzero powers of
Bigelow's explicit specialized Burau-kernel commutator when followed by
\(\sigma_1\sigma_2\sigma_3\).

## Primary literature

- Bigelow's Proposition 1.5 gives the equivalence of the faithfulness of the
  Burau, Temperley–Lieb, and Jones representations of \(B_4\). Section 3
  explains how an element of a hypothetical kernel can yield a closure with
  unit Jones polynomial. Nontriviality of the resulting knot is a separate
  issue, expressed in Conjectures 3.1 and 3.2. Corollary 6.2 supplies a
  nontrivial braid in the Burau kernel **at the specialization \(q=2\)**.
  [Bigelow, *Does the Jones polynomial detect the unknot?*](https://arxiv.org/pdf/math/0012086)
- The Main Theorem and Corollary 1.1 in the September 2026 preprint of
  Bharathram, Birman, and Brendle report faithfulness of the generic Burau and
  Jones representations on \(B_4\). Faithfulness concerns equality of
  representations of braids. The closed-knot invariant retains a weighted
  trace, so equality of that invariant can still occur.
  [Bharathram, Birman, and Brendle, arXiv:2607.05283v2](https://arxiv.org/html/2607.05283v2)
- Korzun, Lanina, and Sleptsov's Theorem 3.1 excludes a four-strand knot
  with unit Jones and Alexander polynomials and nonunit HOMFLY-PT
  polynomial. Their Theorem 3.2 and Section 4.1 give necessary HOMFLY-PT
  forms for the nonunit HOMFLY-PT branch, with writhe \(2m+5\) up to
  mirroring, \(m\ge 0\). These forms constrain a hypothetical braid
  closure. Realizing a form by a knot and treating the unit HOMFLY-PT
  branch are separate questions.
  [Korzun, Lanina, and Sleptsov, *Closed 4-braids and the Jones unknot conjecture*](https://ems.press/content/serial-article-files/53933?nt=1)
- Kauffman's Theorem 2.8 gives the normalized bracket formula used below.
  Champanerkar and Kofman's Proposition 2.1 gives the general twist-bracket
  expansion that underlies the single-site calculation.
  [Kauffman, *State models and the Jones polynomial*](https://webhomes.maths.ed.ac.uk/~v1ranick/papers/kauffmanjones.pdf),
  [Champanerkar and Kofman, *On links with cyclotomic Jones polynomials*](https://msp.org/agt/2006/6-4/agt-v6-n4-p03-s.pdf)
- Tuzun and Sikora's Theorem 1 proves that any knot with a diagram having at
  most 24 crossings and unit Jones polynomial is the unknot.
  [Tuzun and Sikora, *Verification of the Jones Unknot Conjecture up to 24 Crossings*](https://arxiv.org/html/2003.06724)

## Exact reduction for one even twist

Let \(B\) be a braid whose permutation is a single cycle. Write \(E_i\) for
the Temperley–Lieb diagram connecting the \(i\)-th and \((i+1)\)-st strands
by a cap and a cup with all other strands vertical. With
\(\delta=-A^2-A^{-2}\), one has \(E_i^2=\delta E_i\). The writhe-normalized
braid image \(X(B)\) is multiplicative and sends \(\sigma_i\) to
\(R_i=-A^{-2}I-A^{-4}E_i\). The trace \(\tau\) closes a diagram and evaluates
\(\ell\) circles as \(\delta^{\ell-1}\). Thus
\(\tau(X(B))=V_{\widehat B}(A^{-4})\).

Set \(F_k(A)=\tau(X(B)R_i^{2k})\) and \(C(A)=\tau(X(B)E_i)\).
Over \(\mathbb Q(A)\), the two projectors \(I-E_i/\delta\) and \(E_i/\delta\)
give eigenvalues \(-A^{-2}\) and \(A^{-6}\), respectively, for \(R_i\).
For every integer \(k\) this gives

\[
F_k(A)=A^{-4k}P(A)+A^{-12k}Q(A),\qquad
Q=C/\delta,\quad P=F_0-Q.
\]

All closures in this sequence are knots because \(\sigma_i^{2k}\) is pure.
For a knot the Jones skein relation gives \(V(1)=1\) and \(V'(1)=0\).
Indeed, at \(t=1\) crossing changes preserve the value, so a link of \(r\)
components evaluates as \((-2)^{r-1}\) via the \(r\)-component unlink.
At a knot crossing both signed diagrams have value \(1\) and the oriented
smoothing has value \(-2\). Differentiating the skein relation at \(t=1\)
therefore gives equality of the two knot derivatives. The unknot derivative
is zero. Applying these facts to \(F_k\) shows

\[
P(1)=3/2,\quad Q(1)=-1/2,\quad
P'(1)+Q'(1)=0 .
\]

The first two equalities follow by taking the value at \(k=0\) and
subtracting the derivative identities for \(k=0\) and \(k=1\). Since
\(\delta(1)=-2\) and \(\delta'(1)=0\), the coefficient of \((A-1)^2\) is

\[
[(A-1)^2]F_k=[(A-1)^2]F_0+4C'(1)k-24k^2 .
\]

Every unit Jones value must be one of the at most two integer roots of
this quadratic. If two distinct exponents give unit Jones, they are
consecutive. To see this, shift the smaller exponent to \(0\), let the
larger be \(h>0\), and put \(z=A^{-4h}\). Then \(P+Q=1\) and
\(zP+z^3Q=1\) imply \(Q=-1/(z(1+z))\). The Laurent polynomial \(C=\delta Q\)
can exist only if \(1+z\) divides \(\delta\) in
\(\mathbb Q[A,A^{-1}]\). Their respective spans in \(A\) are \(4h\) and \(4\).
Spans add in this integral domain, so \(h=1\). The bound is attained
by \(B=\sigma_1\sigma_2\sigma_3\), twisting on \(\sigma_1\) at \(k=-1,0\).

## Seven-block computation

The first three odd blocks give the permutation of
\(\sigma_1\sigma_2\sigma_3\), a four-cycle. Each of the four even blocks
leaves that permutation unchanged. The calculation covers every tuple of
the first six exponents. For each tuple it computes the quadratic and all
of its integer roots by the discriminant and divisibility by 48. It then
computes the full normalized bracket in \(\mathbb Z[A,A^{-1}]\) at each root.

The diagram computation starts from the four vertical strands. Appending
\(E_i\) pairs the two chosen bottom points with each other and joins their
previous partners. If the bottom points were already paired, it instead
adds one circle. Iterating these rules reaches 14 pairings. Closing a
pairing joins corresponding top and bottom points, and connected
components give the circle count. Appending individual signed crossings
uses \(R_i^{\pm1}=-A^{\mp2}I-A^{\mp4}E_i\) in the Laurent calculation. The
quadratic calculation uses the same pairings and works modulo
\((A-1)^3\), with a twist of \(m\) crossings given by

\[
R_i^m=(-A^{-2})^m I
       + \frac{A^{-6m}-(-A^{-2})^m}{\delta}\,E_i.
\]

The quotient has integer Laurent coefficients because it is the
coefficient of \(E_i\) in the power of the local signed-crossing formula.
The code expands it by truncated division near \(A=1\) and checks
integrality. The full Laurent computation uses the two signed generators
directly.

| Count or bound | Value |
| --- | ---: |
| Six-block prefixes | 175,616 |
| Candidate tuples after the quadratic screen | 43,256 |
| Largest \(\lvert k\rvert\) among candidate tuples | 17 |
| Tuples with unit Jones polynomial | 2,610 |
| Largest original word length among unit tuples | 29 |
| Unit tuples with original word length above 24 | 272 |
| Largest shortened length among those 272 tuples | 16 |

The counts refer to parameter tuples. For the 272 longer words, cyclic
rotation, commutation of distant generators, inverse cancellation, and
outer-strand Markov destabilization supply shorter diagrams of the same
knots. For a destabilization the outer generator occurs once, so a cyclic
rotation makes its single signed crossing the last letter and all other
letters act on the remaining strands. The left version follows by
reflecting the strand indexing before and after the usual right-end
destabilization. The 24-crossing theorem applies to all 2,610 unit tuples,
so each is an unknot. The other word families in \(B_4\) lie outside the
enumeration.

## Powers of Bigelow's specialized witness

Write \(a=\sigma_1,\ b=\sigma_2,\ c=\sigma_3\). Use ordered products and
commutator convention \([x,y]=xyx^{-1}y^{-1}\). In the notation of
Bigelow's Corollary 6.2 set

\[
\begin{aligned}
\psi &= a^{-3}b^{-2}c^{-1}bc^4b^{-1}cbabc^2ba^{-1}b^{-1}c^{-2},\\
D &= abc,\qquad \beta=[(ba)^3,\psi^{-1}b\psi].
\end{aligned}
\]

The braid \((ba)^3\) is pure, so \(\beta\) is pure and every
\(\widehat{\beta^jD}\), \(j\in\mathbb Z\), is a knot. At \(A=1\) the
normalized braid image factors through the permutation group because
\(R_i(1)^2=I\). Hence \(X(\beta)-I\) is divisible by \(A-1\). In the algebra
modulo \((A-1)^3\),

\[
X(\beta)^j=I+j(X(\beta)-I)
 +\binom{j}{2}(X(\beta)-I)^2
\]

for all integers \(j\). In particular
\([(A-1)^2]\tau(X(\beta)^j X(D))\) is a polynomial of degree at most two
in \(j\). Exact Laurent trace evaluations at \(j=-1,0,1\) give
\(384,0,-384\), respectively. The coefficient is therefore \(-384j\).
Every nonzero exponent has nonunit Jones polynomial. At \(j=0\) the
closure of \(D\) is the unknot by successive destabilizations.

## Replay and scope

From the repository root, run the
[integer trace computation](search_four_braid_twists.py) with ordinary
Python 3. The recorded local environment was Python 3.12.3 under WSL.

```sh
python3 research/jones-unknot/search_four_braid_twists.py
```

The script generates its exact inputs deterministically and prints the
three Bigelow interpolation coefficients, the search counts, and the
list of unit tuples still above 24 crossings after shortening. That list
is empty at the stated bounds. Assertions also check the unknot and
trefoil normalizations and the second coefficient of every Laurent
polynomial surviving the quadratic screen. This computation and the
published crossing theorem establish the family exclusion. The unrestricted
conjecture and the unrestricted four-strand scalar trace problem remain
open.
