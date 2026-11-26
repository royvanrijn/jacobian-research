# Weighted dyadic matching controls for zero-constant binary quadratics

Let \(f=A X^2+BXY+C Y^2+dX+eY\) be nonzero, with
\(P_n=\Lambda_2(f^n)/n!\) and \(\Lambda_2(X^iY^j)=i!j!\).
Use the natural coefficients \(\alpha=2A,\ \beta=B,\ \gamma=2C\).
Work at any place above \(2\), with \(v(2)=1\).
The following are necessary local conditions for a block
of five zero moments, and give four-window exclusions
in the stated sectors. The residue field is arbitrary
of characteristic two.

After division of the phase by a coefficient of least
valuation, any block of four or more zero moments
requires one of two configurations:
- the cross coefficient is the only unit among the five
  natural coefficients;
- both linear coefficients are nonunits, exactly one axis
  coefficient is a unit, and the cross coefficient has the
  same residue as that axis coefficient.

We give exact refinements for both cases. In the first,
the relative gap between the cross and the other valuations
determines a sharp leading residue pattern with a possible
tie at a dyadic threshold. In the second, a unit pole gives
a congruence modulo \(2\) in the valuation ring, including
when that ring is ramified.

## The two exceptional residues

The [universal matching identity](QUADRATIC_STRONG_FACTORIAL_WINDOWS_FROM_MATCHINGS.md#2-the-universal-matching-generating-function)
gives the coefficientwise congruence, for integral natural
coefficients,
\[
 \sum_{n\ge0}P_nT^n\equiv N(T)/D(T)\pmod{2\mathcal O[[T]]},
 \quad
 N=1+(d+e+\beta)T+deT^2,\quad
 D=(1+\alpha T+d^2T^2)(1+\gamma T+e^2T^2).                 \tag{1}
\]
Indeed both coefficient arrays consist of integer polynomials
in the natural parameters. Their polynomial congruence in
characteristic two can therefore be reduced modulo the
ideal \(2\mathcal O\).

Reduce (1) into the residue field. If both linear
coefficients are units, \(\deg\bar D=4\) and
\(\deg\bar N\le2\). If exactly one is a unit,
\(2\le\deg\bar D\le3\) and \(\deg\bar N\le1\).
If both are nonunits and both axes are units, the
respective degrees are \(2\) and at most \(1\).
In each case the proper rational series has no four
consecutive zeros at positive indices. For a formal
series \(P/D\), with \(D(0)=1,\ \deg D=r\) and
\(\deg P<r\), a block of \(r\) zeros propagates
backward by the coefficient recurrence, contradicting
a nonzero constant coefficient. This argument works
over any residue field and does not need coprime
numerator and denominator.

If both linear coefficients are nonunits and exactly
one axis is a unit, say \(\bar\alpha\ne0\), then
\(\bar N/\bar D=(1+\bar\beta T)/(1+\bar\alpha T)\).
Its coefficient at \(n\ge1\) is
\((\bar\beta+\bar\alpha)\bar\alpha^{n-1}\),
so even one positive zero requires
\(\bar\beta=\bar\alpha\). Finally, if both linears
and both axes are nonunits, normalization forces
the cross to be the only unit. This proves the classification.

## A unit-factor congruence at the full modulus

Let \(D=UV\) with \(U(0)=V(0)=1\), where \(U\)
contains the factors \(1-\mu T\) with \(\mu\) a unit
and \(V\) those with \(v(\mu)>0\), counted with
multiplicity. All \(\mu\)'s are integral, since
they are roots of the reversed monic denominator.
Work in a finite splitting extension when necessary.
Suppose \(U\) has positive degree. Any block of four
consecutive positive zero moments implies
\[
                         U\mid N\pmod{2\mathcal O[T]}.       \tag{2}
\]
Polynomial division by \(U\) is integral, since its
leading coefficient is a unit.
Indeed, by (1) a zero block for \(P_n\) is a zero block
modulo \(2\mathcal O\) for \(N/D\).
Multiplication by \(V\) of degree \(r_V\) transports it
to a block of at least \(4-r_V\ge\deg U\) zeros for \(N/U\),
starting at index \(n+r_V\). For binary \(D,N\) as in
(1), whenever \(U\) has positive degree one has
\(\deg N\le\deg D\). The quotient polynomial in the
division of \(N\) by \(U\) therefore has degree at most
\(r_V\), and has no effect in the transported block.
The remaining proper series has denominator \(U\)
and its recurrence propagates backward modulo \(2\mathcal O\),
because the leading coefficient of \(U\) is a unit.
All its initial coefficients are zero modulo \(2\mathcal O\).
Multiplying by \(U\) gives (2).

In the equal-axis-residue configuration with \(\alpha\) the
unit axis, let \(r\) be the unique unit root of
\(1+\alpha T+d^2T^2\). Its reduction is \(1/\bar\alpha\)
and it is simple there, so it exists by Hensel's lemma.
The other vertex polynomial has no unit root.
Evaluating (2) at \(r\) gives \(v(N(r))\ge1\).
Using the equation at \(r\) and canceling one factor \(r\)
gives the concrete necessary congruence
\[
               v\big(\beta-\alpha+d+e+d(e-d)r\big)\ge1 .   \tag{3}
\]
The symmetric condition holds with the other axis.
In particular, a failure of (3) proves the four-window
locally in the exceptional equal-residue sector.
The congruence works at every dyadic place of the
coefficient field.

## Exact cross weights below the unit valuation

Assume now that the cross is the only unit, and that
some noncross natural coefficient is nonzero.
Let \(c_0\) be a noncross coefficient of smallest
valuation \(s\), and assume \(0<s<1\).
Scale the entire phase by \(c_0\). This keeps the
\(P_n\) zero sets and makes its noncross coefficients
integral with at least one unit. Its cross coefficient
\(\beta'=\beta/c_0\) has valuation \(-s\).
Use the scaled noncross coefficients to form the two
vertex polynomials \(D_X,D_Y\) over the residue field,
and put \(D_0=D_XD_Y\), \(r_0=\deg D_0\).
At least one vertex coefficient is a unit, so
\(1\le r_0\le4\).

For a single vertex with scaled natural coefficients
\(\alpha',d'\), write
\[
 J_j(T)=\frac1{j!}\Lambda_1
        \big(X^j\exp(T(\alpha'X^2/2+d'X))\big).
\]
These series are integral. In the literal coefficient
with \(p\) quadratics and \(l\) linears the integer
coefficient is
\((2p-1)!!\binom{j+2p+l}{j,2p,l}\).
Its reduction is one exactly when the binary bits of
\(j,2p,l\) are disjoint, by the multinomial theorem
in characteristic two. At bit zero the local factor
is \(1+\bar d'T\) if the bit of \(j\) is unused;
at bit \(b\ge1\) it is \(D_X(T)^{2^{b-1}}\) if
the bit is unused. The product for all bits \(b\ge1\)
equals \(1/D_X(T)\) as a formal series, by multiplying
the identity \((1+U)\prod_{i=0}^{m}(1+U^{2^i})
=1+U^{2^{m+1}}\) with \(U=D_X-1\).
Removing the factors for the bits of \(j\) therefore gives
\[
 \overline{J_j}(T)=
   \frac{(1+\bar d' T)^{1-(j\bmod2)}}{D_X(T)^{\,1+\lfloor j/2\rfloor}} .
                                                               \tag{4}
\]
Use the corresponding formula for the second vertex.

The cross expansion of the scaled series is
\(\sum_{j\ge0} j!(\beta')^j T^j J_{j,X}J_{j,Y}\).
Its coefficient valuations in a fixed \(j\)-summand
are bounded below by
\[
     v(j!)-js=(1-s)j-s_2(j)
           =\sum_{b\ge0}j_b\bigl((1-s)2^b-1\bigr).           \tag{5}
\]
This lower bound is attained at the summand's first term
because each \(J_j(0)=1\). Normalize the entire series
by a scalar of valuation the minimum in (5). The residue
is the sum of precisely the minimizing cross groups;
there are at most two. The minimum exists because \(s<1\).
If \(s\) is between two successive thresholds
\(1-2^{-(b-1)}<s<1-2^{-b}\), with \(b\ge1\),
the unique minimizing count is \(j_0=2^b-1\).
Writing \(Q=2^{b-1}\), (4) gives the normalized residue
\[
                         \text{unit}\cdot T^{j_0}/D_0(T)^Q.  \tag{6}
\]
It has support on \(j_0+Q\mathbb Z_{\ge0}\).
Its sampled coefficients, starting at \(j_0\), obey the
proper recurrence of \(1/D_0^{(Q)}\) of degree \(r_0\),
where \({}^{(Q)}\) is coefficient Frobenius.
Hence a window of length \(\max(j_0,Qr_0)\) at any
positive start contains a nonzero coefficient and a
nonzero exact moment. In particular \(0<s<1/2\) gives
a four-window. The range \(1/2<s<3/4\) gives a
four-window whenever \(r_0\le2\).

## The tied residue and its exact limitations

At \(s=1-2^{-b}\) with \(b\ge1\), the minimizing
counts are \(j_0=2^b-1\) and \(j_1=2^{b+1}-1\).
Their difference is \(2Q\), with \(Q=2^{b-1}\).
After removing a unit, the normalized series has residue
\[
 T^{j_0}\frac{D_0(T)^Q+\xi T^{2Q}}{D_0(T)^{2Q}},
                              \qquad \xi\ne0.               \tag{7}
\]
Here \(\xi\) is the residue of the ratio of the two
cross multipliers \(j!(\beta')^j\).
By Frobenius the sampled series is
\((D(S)+\xi S^2)/D(S)^2\), where \(S=T^Q\) and
\(D=D_0^{(Q)}\).

For clarity, consider \(T(D(T)+\xi T^2)/D(T)^2\),
where \(D=1+D_1T+D_2T^2+D_3T^3+D_4T^4\).
Its even and odd coefficient streams, indexed by
powers \(S=T^2\), are respectively
\[
 \frac{S(D_1+D_3 S)}{D^{(2)}(S)},\qquad
 \frac{1+(D_2+\xi)S+D_4S^2}{D^{(2)}(S)}.                    \tag{8}
\]
When \(\deg D=1\), the even stream after removing \(S\)
has a nonzero numerator of degree zero and denominator
of degree one, so every positive even coefficient is nonzero.
When \(\deg D=2\), if \(D_1\ne0\) use the even stream
with denominator of degree two. If \(D_1=0\), the
odd stream has numerator of degree at most one,
constant one and denominator of degree two.
Thus \(r_0\le2\) gives a four-window at \(s=1/2\),
improved to a two-window when \(r_0=1\).

At \(\deg D=3\), a block of five starting at an odd
index has three consecutive odd coefficients, impossible
by the second proper recurrence in (8).
A block starting at an even index has three consecutive
even coefficients. Remove the initial \(S\) in the
first stream, and if \(D_1=0\) remove another \(S\);
the remaining proper series has nonzero constant and
degree-three denominator. If the block starts before
the shifted first nonzero even coefficient, it contains
that coefficient. Thus a five-window holds in this case.

There is also a four-window when \(\deg D=3\) and
\[
                         D_3\ne D_1(D_2+\xi).              \tag{9}
\]
Let \(a_m\) denote the coefficients of \(1/D^{(2)}\).
The two stream values at index \(m\) after removing
the initial \(S\) from the first are
\(D_1a_m+D_3a_{m-1}\) and \(a_m+(D_2+\xi)a_{m-1}\).
The determinant in (9) makes their simultaneous zero
force \(a_m=a_{m-1}=0\).
A four-zero block starting odd has both stream values
at two consecutive \(m\)'s, and a block starting even
has both at one \(m\) and the second stream at \(m+1\).
In either case it forces three consecutive zeros in
\(1/D^{(2)}\), a contradiction. The initial indices
are covered by its constant coefficient one.
In particular (9) holds if one axis residue is zero
and \(\deg D=3\). The degree-two factor then has
the form \(1+d^2 T^2\) and the other factor is
\(1+\gamma T\) with \(\gamma\ne0\).
The determinant simplifies to \(\gamma\xi\ne0\).

When \(\deg D=4\), the numerator of
\((D+\xi S^2)/D^2\) has degree at most four,
so the eight-window is immediate. Seven consecutive
leading-residue zeros can occur over \(\mathbb F_4\).
With \(\omega^2+\omega+1=0\) take
\[
       \bar\alpha'=\bar\gamma'=0,\qquad
       \bar d'=\omega,\quad\bar e'=\omega^2,\quad \xi=1.
\]
Then \(D=1+T^2+T^4\) and the shifted series in (8) is
\[
                    \frac{T(1+T^8)}{1+T^{12}},
\]
whose coefficients at \(2,\ldots,8\) vanish.
These data arise from ordinary dyadic inputs: take
\(c_0=\sqrt2,\ B=1\), zero axes, and \(d,e\) equal to
\(\sqrt2\) times unramified lifts of \(\omega,\omega^2\).
The minimizing cross counts are \(1,3\), with ratio
\(6(B/c_0)^2=3\), of residue one. This shows that
the leading tied residue by itself leaves a gap for the
five-window.

For arbitrary \(b\) in (7), the same parity arguments
apply to its sampled series. They give a window length
\(Qm\) for the original positive-index series with
\(m=2,4,5,8\) for \(r_0=1,2,3,4\), respectively.
Condition (9) allows \(m=4\) for \(r_0=3\).
In these bounds \(Qm\ge j_0\), so a starting index
before the first coefficient is covered as well.
Cases beyond the stated local exclusions still require
a higher-precision comparison or an independent place.

## Later exclusions in the two residue types

In one-axis resonance, [a missing linear term](ZERO_CONSTANT_ONE_AXIS_MISSING_LINEAR.md)
gives the four-window at every positive start.
The [combined one-axis exclusions](ZERO_CONSTANT_ONE_AXIS_THREE_HALVES_BOUNDARY.md)
force \(e\ne0\) and \(v(e)\ge7/4\) for a five-zero block when
\(\alpha=1\) is the unit axis. Both statements allow ramification.
In the pure-cross case with \(s>1\), the
[integer-anchor argument](ZERO_CONSTANT_PURE_CROSS_INTEGER_ANCHORS.md)
gives unit normalized moments when the second-order coefficient
vanishes, and exclusions in deep integer congruence classes
when it is nonzero.
