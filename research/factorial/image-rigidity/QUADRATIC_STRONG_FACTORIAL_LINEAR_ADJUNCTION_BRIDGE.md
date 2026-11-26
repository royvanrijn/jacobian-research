# Reversible linear adjunction for dense quadratic Strong Factorial windows

Let \(K\) have characteristic zero, and let \(\Lambda_K\) send a monomial to the product of its exponent factorials. A polynomial with \(s\) nonzero monomials has its Strong Factorial property if no \(s\) consecutive positive-index moments vanish.

The three classes considered here are
\[
 g(x)=ax^2+bx+c,\qquad
 F(x,y)=g(x)+dy,\qquad
 H(x,y,z)=g(x)+dy+ez,                         \tag{1}
\]
with every displayed coefficient nonzero. Their support sizes are respectively \(3,4,5\).

**Theorem 1 (equivalence of the three universal assertions).** The following assertions are equivalent:

1. Every dense \(g\) in (1) has its three-moment Strong Factorial property.
2. Every four-term \(F\) in (1) has its four-moment property.
3. Every five-term \(H\) in (1) has its five-moment property.

The equivalence holds over any one fixed algebraically closed characteristic-zero field, in particular over \(\mathbb C\). It also holds when each assertion is universally quantified over all characteristic-zero fields. The reverse construction over an arbitrary fixed field can require an algebraic extension; no equivalence over each fixed non-algebraically-closed field is asserted.

These are restricted classes of genuine binary and ternary polynomials. The theorem does not identify either assertion with full \(SFC(2)\), arbitrary six-term binary quadratics, or full \(SFC(3)\).

## 1. The exact adjunction identity

For a polynomial \(P\) independent of a new variable \(y\), write
\[
 G_j=\Lambda(P^j),\quad
 A_j(d)=\Lambda((P+dy)^j),\quad G_0=A_0=1.
\]
The binomial theorem and \(\Lambda(y^i)=i!\) give
\[
 A_j(d)=\sum_{i=0}^j\frac{j!}{(j-i)!}d^iG_{j-i},
 \qquad
 A_j(d)=jd\,A_{j-1}(d)+G_j\quad(j\ge1).       \tag{2}
\]
In particular \(A_j(d)\), regarded as a polynomial in \(d\), has degree exactly \(j\), leading coefficient \(j!\), and constant coefficient \(G_j\).

If \(A_m,\ldots,A_{m+s}\) vanish, then (2) gives
\[
 G_{m+1}=\cdots=G_{m+s}=0.                    \tag{3}
\]
Thus adjoining a nonzero independent linear monomial carries an \(s\)-moment property to an \((s+1)\)-moment property. This forward result is classical: Rocks' author preprint, Theorem 3.12, gives linear adjunction, and Corollary 3.13 gives its iteration. Equation (2) supplies the complete argument needed here.

The primary statement and supplied proof were read in the [Rocks author preprint](https://people.math.sc.edu/rocks/BradyRocksSFC2015.pdf). The reverse implication below is proved locally; no novelty assertion is made.

## 2. The two initial windows cannot vanish

For the dense quadratic \(g=ax^2+bx+c\), let \(G_j=\Lambda(g^j)\). Its native factorial recurrence is
\[
 G_j=2aj(2j-1)G_{j-1}
       +j(j-1)(b^2-4ac)G_{j-2}
       +c^j+bj c^{j-1}\quad(j\ge2).           \tag{4}
\]
This identity follows by integration by parts over \(\mathbb C\), or by comparing factorial coefficients, and is a polynomial identity over \(\mathbb Q\), hence valid over every \(K\).

If \(G_n=G_{n+1}=G_{n+2}=0\), equation (4) at \(j=n+2\) gives
\[
 c^{n+1}(c+(n+2)b)=0.
\]
Since \(bc\ne0\), this forces \(c=-(n+2)b\). Divide the phase by \(b\) and put \(A=a/b\).

For an initial three-zero window \(n=1\), the normalized phase is \(Ax^2+x-3\). Its first moment \(2A-2\) forces \(A=1\), but its second moment there is \(29\ne0\). Therefore
\[
 \text{the dense quadratic class has no three-zero window starting at }1. \tag{5}
\]

Suppose \(F=g+dy\) has its initial four moments zero. Equations (2)--(3) force \(G_2=G_3=G_4=0\), and (4) forces the normalization \(Ax^2+x-4\). Its second and third native moments are
\[
 P(A)=24A^2-4A+10,\qquad
 Q(A)=720A^3+72A^2+24A-34.
\]
The explicit polynomial identity
\[
 (-10980A^2+4032A+5577)P(A)+(366A-232)Q(A)=63658             \tag{6}
\]
excludes a common root in every characteristic-zero field. Thus
\[
 \text{the four-term class has no four-zero window starting at }1.       \tag{7}
\]
These initial exclusions do not assume any unresolved Strong Factorial assertion.

## 3. The reverse counterexample lift

The following observation makes the field scope explicit. Suppose a polynomial \(P\) has a block of \(s\) consecutive zero moments. Move its start left to the beginning \(r\) of that zero run. Then
\[
 G_{r-1}\ne0,\qquad G_r=\cdots=G_{r+s-1}=0,                \tag{8}
\]
where \(G_0=1\) if \(r=1\). If the class has no zero window starting at \(1\), then \(r\ge2\).

By (2), \(A_{r-1}(d)\) is nonconstant of degree \(r-1\), with nonzero leading coefficient \((r-1)!\) and nonzero constant \(G_{r-1}\). It therefore has a nonzero root \(d\) in an algebraic closure. For this choice,
\[
 A_{r-1}(d)=0.
\]
Equation (2) and the successive zero values in (8) then give
\[
 A_{r-1}(d)=A_r(d)=\cdots=A_{r+s-1}(d)=0.                  \tag{9}
\]
This is an \((s+1)\)-zero block with a positive start. All coefficients of \(P\) are retained and the new coefficient is nonzero.

Apply this first to dense \(g\), using (5). Any three-window counterexample lifts to a four-window counterexample \(g+dy\). Apply it again to that four-term polynomial, using the unconditional exclusion (7). It lifts to a five-window counterexample \(g+dy+ez\), with \(e\ne0\). Conversely, the forward implication (3), once or twice, sends any four- or five-window counterexample back to a three-window counterexample for \(g\). This proves Theorem 1.

Over an algebraically closed \(K\), both choices of roots stay in \(K\). Under universal quantification over all characteristic-zero fields, a lift in a finite algebraic extension is enough to disprove the corresponding universal assertion. This proves that formulation without assuming that an arbitrary large coefficient field embeds in \(\mathbb C\).

## 4. Complete binary and ternary consequences

**Corollary 2.** Over every characteristic-zero field, every \(F\) in (1)
has a nonzero moment in every positive block of length four, and every
\(H\) in (1) has one in every positive block of length five. The
[dense univariate quadratic theorem](DENSE_UNIVARIATE_QUADRATIC_STRONG_FACTORIAL.md#5-the-dense-univariate-quadratic-conclusion)
gives the three-window property for \(g\) over that field. Applying
(2) once or twice gives the claimed windows for \(F\) and \(H\).

These results allow arbitrary nonzero coefficients in (1). General binary
quadratics, including six-term supports, and unrestricted two- and
three-variable Strong Factorial remain separate problems.

### Earlier resonance criterion

The separately proved dense univariate resonance criterion gives an explicit set \(\mathcal R\) of possible ratios for three-zero blocks: if such a block starts at \(n\), then
\[
 -c/b=n+2\in\mathcal R.                                      \tag{10}
\]
For an exact definition and its complete local proofs, see [the dense quadratic resonance criterion](UNIVARIATE_QUADRATIC_STRONG_FACTORIAL_RESONANCE.md), Sections 1--5.

Equations (2)--(3) now imply, without a converse for a specified coefficient choice,
\[
 \begin{array}{ll}
 \Lambda(F^m)=\cdots=\Lambda(F^{m+3})=0
   &\Longrightarrow -c/b=m+3\in\mathcal R,\\
 \Lambda(H^m)=\cdots=\Lambda(H^{m+4})=0
   &\Longrightarrow -c/b=m+4\in\mathcal R.
 \end{array}                                                 \tag{11}
\]
Consequently every phase in either restricted class whose ratio \(-c/b\) is outside that explicit set satisfies its complete support-length property. The earlier criterion first left \(q=91\), corresponding to starts \(88\)
for \(F\) and \(87\) for \(H\). Corollary 2 covers those starts and every
other positive start.

The exact equivalence in Theorem 1 concerns universal assertions for the whole three classes, not an assertion that every particular \(d\), or every particular pair \(d,e\), preserves a counterexample. Its reverse constructs suitable nonzero coefficients in an algebraic closure. No actual quadratic counterexample, full binary or ternary Strong Factorial theorem, or implication from ordinary Factorial is asserted.

## 5. Supplementary native identity checks

The [portable standard-library checker](verify_quadratic_linear_adjunction_bridge.py) reconstructs the finite native adjunction identities and the two small initial-window calculations. The [retained receipt](quadratic_linear_adjunction_bridge_verification.json) records these exact checks. From the directory containing this note, run a fresh verification:

    python verify_quadratic_linear_adjunction_bridge.py --output fresh_adjunction_verification.json

An existing output file is never overwritten, and optimized Python is refused. These finite controls verify the displayed coefficient transcriptions and literal identity (6); the universal root-existence and counterexample-lift arguments are the written proofs above. The checker does not construct a quadratic counterexample or infer a universal assertion from tested starts.
