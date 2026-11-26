# Finite height reduction for opposite unit-gap rows

Fix an ordinary core with positive factorial coefficients
\(c_{j,k}=[z^k]Q_j^C(z)\). Append \(r\) copies of \((A,A+1)\)
and \(s\) copies of \((B+1,B)\), with \(r,s\) fixed positive
integers and independent nonnegative integer heights \(A,B\).
For two coprime orders \(n,N\ge2\) put
\[
 F_j(A,B,z)=\sum_{k=0}^j c_{j,k}(jA+1)_k^r
                                      (jB+1)_{j-k}^s z^k.
\]
This strips only the factorials of the smaller endpoints.
Either order may be labeled \(N\).

For \(1\le h\le N\), put \(A_0=-1/N,\ B_h=-h/N\) and
\[
 \begin{split}
 P_h(z)&=z^{-\boldsymbol 1_{h=N}} F_n(A_0,B_h,z),\\
 S_h(z)&=\sum_{k=N-h+1}^N c_{N,k}
        [N(k-1)!]^r (1-h)_{N-k}^s z^k .
 \end{split}                                               \tag{1}
\]
The polynomial \(P_h\) has nonzero endpoints. Denote by \(P_h^*\)
its monic squarefree part. The two exact tests are
\[
 \begin{split}
 &\operatorname{Res}_z(P_h^*(z),T-S_h(z))
                                      \text{ is squarefree in }T,\\
 &\gcd(P_h,P_h',S_h')=1 .
 \end{split}                                               \tag{2}
\]
Equivalently \(S_h\) separates distinct roots of \(P_h\), and
its derivative is nonzero at every repeated root of \(P_h\).

**Theorem.** If (2) holds at every \(h\) in (1), the complex
height locus where the two polynomials have a common factor of
degree at least two is zero-dimensional on the chart where
both degrees are full and both constant coefficients nonzero.
If one order is even, there are therefore only finitely many
ordinary height pairs with common roots. These pairs admit an
exact finite height cutoff by algebraic elimination. For
consecutive orders, passing the tests at every pair below the
[absolute terminal cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md)
gives a finite decision for this two-height family in every
degree. The choice of which order degenerates, or which row
increases after exchanging columns, can be made at each pair.

The limiting test here only involves polynomials of degrees
at most \(N\) from the fixed core. Unlike the
[same-orientation reduction](MULTIVARIATE_TWO_HEIGHT_FACTOR_SPECIALIZATION.md),
the opposite row can put a whole block of nonconstant terms
at the same height as the constant at the degeneration.
The evaluations of \(S_h\) record this block.

## Degeneration of a possible curve

Use the endpoint chart in the statement. Its closed rank
condition in that chart comes from the Sylvester minors of
size \(n+N-1\). A curve there cannot have \(A\) constant:
specializing \(B=-1/N\) on that line makes the upper
polynomial a nonzero multiple of \(z^N\) and the lower
constant nonzero. The upper leading coefficient is independent
of \(B\), and the lower constant is nonzero by coprimality
of the orders. The same specialization at generic \(A\)
also rules out a two-dimensional component.

On an irreducible curve with nonconstant \(A\), take a sequence
of points in the chart with \(A\to A_0\). Such a sequence
exists since the projection of the plane curve onto \(A\)
contains all but finitely many values, and leaving the chart
removes at most finitely many points of the curve. After
a subsequence, \(B\) has a limit in the projective line.
If \(B\to\infty\), evaluate \(B^{-js}F_j(A,B,B^s x)\).
Its order-\(n\) limit has nonzero leading coefficient,
and its order-\(N\) limit is a nonzero constant.
Homogeneous resultant continuity excludes this case.

At finite \(B\), the upper limit is constant since
\(N A_0+1=0\). The lower leading coefficient stays
nonzero. Thus the upper constant must vanish, so \(B\to B_h\)
for one of (1). For \(h<N\), \(nB_h\) is not an integer
and the lower limit has nonzero constant term.
For \(h=N\), the lower limit has a simple root at zero,
and all other roots have nonzero finite limits. In that
case the small lower root \(z_0\) is of order
\((B-B_N)^s\) with a nonzero leading ratio, by its
constant and linear coefficients. The upper constant
has the same order with nonzero leading ratio. At \(z_0\),
all nonconstant terms have the additional vanishing
factor \((A-A_0)^r\) and are of smaller magnitude than
the upper constant. Thus \(z_0\) cannot be common.
All potential common roots in either case stay bounded
and away from zero.

## The simultaneous limiting factor

Put \(V=(B-B_h)^s/(A-A_0)^r\), which is defined along
the sequence in the chart. Divide the upper coefficients
by \((A-A_0)^r\). The positive indices less than
\(N-h+1\) now tend to zero due to their residual \(B\)
factor, and indices at least \(N-h+1\) tend to
the coefficients of \(S_h\). The constant is \(V\)
times a factor tending to a nonzero number \(C_h\),
since the corresponding rising product at \(B_h\) has
a simple zero before taking its \(s\)-th power.
If \(V\to\infty\) the constant dominates on the bounded
root region. Take a subsequence with \(V\to V_0\) finite
instead. The upper limit is \(S_h(z)+V_0C_h\).

If there are two common roots counted with multiplicity,
the monic quadratic they determine divides both parents.
Its coefficients stay bounded. After taking a subsequence
its limit is monic and divides both limiting polynomials
by monic division. Its roots stay nonzero, so it divides
\(P_h\) also when \(h=N\). But (2) implies that
\(P_h\) and \(S_h(z)+V_0C_h\) have at most one common
root, and if that root is repeated in \(P_h\) it is simple
in \(S_h(z)+V_0C_h\). This contradiction rules out the
curve component.

The degree-at-least-two locus is therefore finite after
saturating by the endpoint products. At ordinary heights an
even moment is strictly positive on the real line, so a common
root gives a common conjugate pair. The univariate elimination
ideals in the two heights then give effective finite upper
bounds by the Cauchy root bound.

## Ordinary coefficients at unbounded heights

At orders \(14,15\), take fixed core rows \((2,5),(7,3)\),
and take \(r=s=2\). All fifteen tests (2) pass. In fact
reduction modulo \(101\) keeps the degree of each \(P_h\),
makes it squarefree, and gives a squarefree norm polynomial.
To compute the latter, use the matrix of multiplication by
\(S_h\) on \(\mathbb F_{101}[z]/(P_h)\). Its characteristic
polynomial is the monic normalization of the displayed
resultant, and its degree is \(\deg P_h\). A squarefree
reduction of this degree certifies the same property in
characteristic zero. Thus for this family, with the two
heights free and the directions opposite, any common root
at this pair must occur at one of finitely many height pairs.
The [exact replay](verify_multivariate_opposite_height_factor_specialization.py)
checks the original core factorials and every specialization.
