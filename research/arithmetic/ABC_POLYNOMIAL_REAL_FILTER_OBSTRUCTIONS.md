# Polynomial obstructions to real-filtered target tests

A small formal radical limits how much a polynomial discrepancy can
vanish at a real angle. For a binary model whose forms are formally
highly powerful, the near-equality window in (SI3) therefore has
only boundedly many primitive parameters per height shell. The
argument concerns the ordinary absolute value of the discrepancy.
It places no condition on the actual prime divisibility of the
specialized values.

## Degree and multiplicity of the discrepancy

Let \(A,B,C\) be pairwise coprime nonzero polynomials over a field
of characteristic zero, with at most one constant. Suppose
\(D=A+B-C\ne0\), and write
\[
                      M=\max\{\deg A,\deg B,\deg C\},
          \qquad S=\deg\operatorname{rad}(ABC).
\]
Then
\[
                       \deg D\ge M-2S+3.                       \tag{PR1}
\]
The bound has content when the formal radical is small relative to
the largest degree.

If \(A,B,C\) are independent over the constants, let \(W\) be
their Wronskian with derivative orders zero, one, and two. It is
nonzero. At an irreducible divisor \(P\) of one letter, the
Wronskian is divisible by at least \(P^{(v_P(ABC)-2)_+}\),
since the column belonging to that letter loses at most two
orders of divisibility. The prime factors of the three letters
are disjoint. Consequently
\[
                  \deg W\ge\deg A+\deg B+\deg C-2S.
\]
Replace a letter of largest degree in the Wronskian by \(D\).
This changes the Wronskian by a sign because the coefficient
of that letter in \(D\) is a unit. The degree of the new
Wronskian is at most
\[
                  \deg A+\deg B+\deg C-M+\deg D-3.
\]
These two inequalities prove (PR1).

In constants rank two, any two of \(A,B,C\) are independent by
coprimality and the hypothesis on constants. Their two-column
Wronskians differ only by nonzero constant factors. For a divisor
of any one letter use a basis containing that letter, giving
its exponent minus at most one in the Wronskian. The degree is
therefore at least \(\deg A+\deg B+\deg C-S\). Using the basis
that omits a letter of maximal degree bounds the degree above
by \(\deg A+\deg B+\deg C-M-1\). It follows that \(M\le S-1\),
which already implies (PR1) because \(M\ge1\) and \(D\ne0\).
Rank one is impossible under these hypotheses.

This is the same formal Wronskian as in the support-depth
argument, with the discrepancy inserted for the upper degree
bound instead of the support-wise lower divisibility bound. Generalized polynomial abc bounds also use Wronskians of linear combinations, with more general support and subsum hypotheses.
[De Bondt, Another generalization of Mason’s ABC-theorem](https://arxiv.org/abs/0707.0434v3)
· [Polynomial support-depth bound](ABC_POLYNOMIAL_SUPPORT_DEPTH.md)

Now let \(A,B,C\) be pairwise coprime homogeneous binary
forms of the same degree \(M\ge1\), and let \(S\) denote their
combined formal radical degree. Every projective root of the
nonzero binary form \(D=A+B-C\) has multiplicity at most
\[
                                  2S-4.                         \tag{PR2}
\]
Indeed, work over an algebraic closure and send that root
to infinity by a linear change of coordinates. If its
multiplicity is \(m\), the resulting univariate discrepancy
has degree \(M-m\). At most one of the letters vanishes
at the selected projective point, so the maximum of the
univariate degrees is still \(M\) and at most one letter
can be constant. The radical degree can only decrease
when the common affine chart discards infinity. If a letter
vanishes at the selected point, the affine radical degree is at
most \(S-1\), so (PR1) gives \(m\le2S-5\). Otherwise all
three affine degrees are \(M\). In rank three replace one
column by \(D\) as in (PR1), and cancel the leading coefficient
between the two remaining degree-\(M\) columns. This gives
one further degree saving in the Wronskian upper bound and
hence \(m\le2S-4\). In rank two, \(m\le M\le S-1\le2S-4\).

## Primitive specializations in a near-equality window

Fix \(M\ge1\), \(\theta\in(0,1)\), and an upper bound
\(S_0\ge3\) for the formal radical degree. In each shell
let the three forms be integral, with coefficient heights
\(H^{o(1)}\) uniformly in the shell. They may vary with
\(H\), subject to the degree, coprimality, and nonzero
discrepancy hypotheses above. Use primitive integer
parameters \((s,t)\) with \(H\le\max(|s|,|t|)\le2H\).
The number satisfying
\[
                       |D(s,t)|\le H^{\theta M+o(1)}
\]
is at most
\[
 O_M\!\left(1+H^{\,2-\zeta+o(1)}\right),
               \qquad \zeta=\frac{(1-\theta)M}{2S_0-4}.
                                                               \tag{PR3}
\]
The error in the window and the subpower coefficient
bounds are taken uniformly. The window may include exact
zeros, which give at most a bounded number of additional
primitive pairs.

Choose one of the two charts with denominator of magnitude
at least \(H\). Its parameter \(x\) is a reduced rational
of height at most \(2H\) and \(|x|\le1\). In that chart the
univariate discrepancy \(d(x)\) satisfies
\[
                      |d(x)|\le H^{-(1-\theta)M+o(1)}.
\]
If \(d\) is a nonzero constant, the inequality is impossible
for large \(H\). Otherwise its leading coefficient is a
nonzero integer, and the product over its at most \(M\)
roots first gives a root within \(H^{-(1-\theta)+o(1)}\).
Distinct roots have separation \(H^{-o(1)}\). For example,
the integral squarefree part has subpower coefficient
height by the fixed-degree polynomial factor bound. Its
nonzero discriminant is an integer, and its other root
differences and leading coefficient are bounded above
by subpower factors. The nearest root must therefore
be unique and real. All other factors in \(d(x)\) are
bounded below by subpower factors, and its nearest
root has multiplicity at most \(2S_0-4\) by (PR2).
The displayed inequality now gives
\(|x-\alpha|\le H^{-\zeta+o(1)}\) for that real root \(\alpha\).

Distinct reduced rationals of denominators at most \(2H\)
have separation at least \(1/(4H^2)\). Each of the
boundedly many resulting intervals contains
\(O(1+H^{2-\zeta+o(1)})\) chart parameters. Primitivity gives
at most two pairs per parameter in a chart. This proves (PR3).

In particular, if the three forms are formally \(k\)-full,
then \(S_0\) can be taken at most \(3M/k\). For
\(k>12/(1-\theta)\), (PR3) gives \(O_M(1)\) parameters per
shell. Here formally \(k\)-full means that every irreducible
factor has exponent at least \(k\). If a real-filtered
countermodel has value height \(X\le H^{M+o(1)}\), its
window \(|A+B-C|\le X^\theta\) is covered by (PR3).
This uses the formal support of the polynomials and
allows arbitrary actual factorizations at the corresponding
integer parameters. It supplies no such bound for candidate
families lacking a polynomial representation.
[Short-interval real filter (SI3)](ABC_SHORT_INTERVAL_RADICAL_FILTER.md#consequences-for-real-filtered-popularity-tests)
