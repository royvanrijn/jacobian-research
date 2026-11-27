# The order-thirty-four resultant in the surviving order-twelve row

In the remaining reduced orders \((12,42)\), the
[order-twelve reduction](POLYNOMIAL_WEYL_ORDER_TWELVE_RESIDUAL.md#the-remaining-differential-order-pair)
gives a common leading factor \(h\) with root multiplicities
three or fifteen and at least one fifteen. The resultant
of the two operators gives an additional constraint at every
multiplicity-three root. Its lowest Euler component commutes
with the two indicial components at shifts six and twenty-one,
and its own shift is seventeen. These three shifts give a finite
list for the local indicial polynomials.

## The resultant at a multiplicity-fifteen root

**Proposition PolynomialWeylOrderTwelveFortyTwoResultant.**
Let \(P,Q\) be a proper polynomial Weyl pair in reduced
differential orders \((12,42)\) over a characteristic-zero
field with \([Q,P]=\kappa\ne0\).
After extending constants, scaling the leading coefficients
to \(h^2,h^7\), and subtracting a polynomial in \(P\) from \(Q\),
there is a monic cubic \(A\) such that
\[
             Q=(U A(P))_+,\qquad U=P^{1/2}.
\]
For some constant \(c\), put
\[
      T=U A(P)-Q-cU^{-1},\qquad
      W=Q^2-P A(P)^2+2c A(P).
\]
Then \(T\) has differential order \(-8\) and leading
coefficient proportional to \(h^{-4/3}\).
The polynomial differential operator \(W\) has order
thirty-four, and \([W,P]=2\kappa Q\).
At a translated root of multiplicity fifteen, its weight
at \((-1,3)\) is seventeen.

At the high root, translate it to zero. The
[local \((6,15)\) crossing](POLYNOMIAL_WEYL_ORDER_TWELVE_RESIDUAL.md#the-crossings-at-gcd-six)
and its
[two-root edge calculation](POLYNOMIAL_WEYL_BERNSTEIN_FORTY_TWO_REAL_BOUND.md#the-two-roots-and-real-coefficients)
give highest symbols proportional to \(R^2,R^7\), where
\[
    R=x^{-3}r(z),\qquad z=x^3y,\qquad r(z)=z L(z),
    \qquad L(0)\ne0 .
\]
Here \(L\) has two roots of multiplicities four and one,
so \(r\) has degree six and multiplicities \(1,4,1\).
In particular \(r\) has a simple root, and a homogeneous
symbol commuting with \(R^2\) has weight a multiple of
three and differential order equal to twice its weight.
Indeed its expression \(x^{-w}U_0(z)\) satisfies
\(w U_0r'-3U_0'r=0\), so \(w/3\) is a nonnegative
integer and \(U_0\) is proportional to \(r^{w/3}\).

Apply the
[fractional-power expansion by differential order](POLYNOMIAL_WEYL_LEADING_FACTOR_RESIDUE.md#fractional-powers-by-differential-order).
Write \(e=\gcd(6,\{\operatorname{mult}_r h\})=3\).
After removing the integer powers of \(P\), the terms
with index divisible by \(e\) in the expansion give the
nonnegative part of \(F=U A(P)\), where \(A\) is monic
of degree three. Suppose there is another term, and
take its largest index \(j\), using powers \(P^{j/(2e)}\).
Then \(Q^2-F^2\) is polynomial and has differential order
\(42+6j/e\), with leading coefficient proportional to
\(h^{7+j/e}\).
Its weight is at least \(21+3j/e\) by the translated root.
Since \([Q^2-F^2,P]=2\kappa Q\), the top must commute
with \(R^2\): a weight above seventeen would otherwise
give a commutator of weight above twenty-one.
Such a commuting top has weight at least the next multiple
of three above \(21+3j/e\) and hence order strictly above
\(42+6j/e\), because \(e\nmid j\).
This is impossible, proving \(Q=F_+\).

The order of \(F-Q\) is negative. At leading order
\(\ell>-11\), its coefficient is proportional to
\(h^{\ell/6}\) by \([F-Q,P]=-\kappa\).
If \(-7\le\ell<0\) and \(\ell\ne-6\), then \(Q^2-F^2\)
has order \(42+\ell\) and weight at least \(21+\ell/2>17\).
The commuting-top comparison again gives a strictly larger
order. Remove order \(-6\) if present by subtracting
\(cU^{-1}\), giving \(T\) of order \(\ell\le-8\).
The polynomial \(W\) then has order at most
\(\max(42+\ell,-12)\le34\), with \(\ell=-\infty\) if \(T=0\).
Its weight is exactly seventeen: anything lower cannot
give its commutator with \(P\), and anything higher
would give a commuting top of order at least thirty-six.

Write the top of \(W\) as \(x^{-17}w(z)\). The top bracket,
after dividing by \(R\), gives
\[
                        3w'r-17wr'=c_0 r^6,\qquad c_0\ne0 .
\]
At a root of \(r\) of multiplicity \(m\in\{1,4\}\),
the leading coefficient has nonzero multiplier
\(3v-17m\) for the integer order \(v\) of \(w\).
Hence \(v=5m+1\). The three distinct roots give
\(\deg w\ge33\). Degree thirty-three would make the
left side have nonzero leading term of degree thirty-eight,
higher than the right side of degree thirty-six.
Thus \(\deg w=34\), forcing \(\ell=-8\).
The leading-coefficient comparison at this order
gives the asserted \(h^{-4/3}\).

## A third indicial shift at a small root

**Theorem PolynomialWeylOrderTwelveFortyTwoSmallRootIndicialList.**
At any multiplicity-three root in the setting of the
preceding proposition, translate that root to \(x=0\)
and put \(E=xD_x\). The lowest Euler components of \(P,Q,W\)
are
\[
                  x^{-6}f(E),\qquad x^{-21}g(E),\qquad
                  x^{-17}w_0(E).
\]
They have degrees \(12,42,34\) in \(E\), respectively,
and the last two commute with the first.
Their root divisors have a common signed divisor \(D\)
of mass two supported on the integers:
\[
 \operatorname{div}_0(f)=A_6D,\qquad
 \operatorname{div}_0(g)=A_{21}D,\qquad
 \operatorname{div}_0(w_0)=A_{17}D .
\]
Here \(A_\ell=1+\mathcal S+\cdots+\mathcal S^{\ell-1}\),
where \(\mathcal S\) translates a root by one.
The admissible \(D\)'s under these divisor and mean
constraints are given by the following finite list.

Let \(\mathcal S_0\) be the additive semigroup generated by
\(6,17,21\) in the nonnegative integers.
Take any subset \(I\) of its gaps such that if \(j\in I\)
and \(j-a\ge0\) for \(a\in\{6,17,21\}\), then \(j-a\in I\).
Define \(C_k=0\) for \(k<0\) and for \(k\in I\), and \(C_k=1\)
otherwise. The list consists of
\[
                   D=[0]+\sum_{k\ge0}(C_k-C_{k-1})[k-|I|].
                                                                  \tag{1}
\]
There are \(1694\) such candidates. This counts the
indicial divisor conditions only.

The
[small-root regular-singular theorem](POLYNOMIAL_WEYL_SMALL_ROOT_SPECTRAL_CONSTRAINTS.md#exact-subprincipal-coefficients-at-small-roots)
gives the first two lowest components, and gives root mean
\(5/2\) for \(f\). To see the third lowest component,
use rational pseudodifferential operators in \(E\), whose
derivation \(x\,d/dx\) preserves \(k[[x]]\).
The root \(U\) has the form \(x^{-3}\) times a series in
descending powers of \(E\) with coefficients in \(k[[x]]\).
Indeed solve \(U^2=P\) in that form: the leading
coefficient after factoring \(x^{-3}\) is \(h/x^3\), a
unit, and each next coefficient has invertible linear
factor \(2h/x^3\). Differentiating any new coefficient, or
shifting \(E\) past \(x^{-3}\), contributes at lower
differential orders. The inverse \(U^{-1}\) similarly has
lowest Euler exponent at least three. Thus \(T\) is
bounded below in Euler exponent, and the coefficient of
each power of \(x\) is a formal Laurent series in \(E^{-1}\)
of degree at most \(-8\). Its coefficient of \(D_x^{-8}\)
is proportional to \(h^{-4/3}\), of \(x\)-order \(-4\)
at this root. Since the leading term of \(D_x^{-8}\) in
Euler notation is \(x^8 E^{-8}\), the lowest Euler exponent
\(l\) of \(T\) is at most four.

If \(l<4\), write the corresponding term as
\(x^l t(E)\) with \(t\) of degree \(j\le-8\).
The component of Euler exponent \(l-6<0\) in
\([T,P]=-\kappa\) vanishes, yielding
\[
                         t(E-6)f(E)=f(E+l)t(E).
\]
Comparison one degree below the top gives \(-6j=12l\),
impossible for \(j\le-8\) and \(l<4\).
Thus \(l=4\), and this lowest term has degree \(-8\) in \(E\).
The expansion \(W=-F T-TF+(cU^{-1}+T)^2\) now gives
lowest Euler exponent \(-21+4=-17\) and degree \(42-8=34\)
there. The two terms at this top degree in \(E\) add,
since their leading coefficients commute.
Its leading bracket with \(P\) at Euler exponent
\(-17-6=-23\) must vanish, because the lowest Euler
exponent of \(2\kappa Q\) is \(-21\).
The commutation with \(Q\)'s lowest component follows
at exponent \(-6-21\) from the scalar bracket of \(P,Q\).

Write \(F_0,G_0,W_0\) for the root divisors of \(f,g,w_0\).
The commutation for the steps six and seventeen gives
\(A_{17}F_0=A_6 W_0\).
On each integer-translation orbit use the Laurent
polynomial ring in the shift \(\mathcal S\).
Since \(A_6,A_{17}\) are coprime monic polynomials,
their divisibility gives integral \(D\) with
\(F_0=A_6D,\ W_0=A_{17}D\), by Gauss's lemma.
Commutation for the steps six and twenty-one then gives
\(G_0=A_{21}D\) by cancellation of \(A_6\).
Thus \(D\) has mass two and has nonnegative mass on each
translation orbit. Its mass on the integer orbit is at
least one because \(F_0\) contains the falling-factorial
divisor of length six.

To see that no other orbit occurs, use the
[small-root sector congruence](POLYNOMIAL_WEYL_SMALL_ROOT_SPECTRAL_CONSTRAINTS.md#indicial-constraints-for-each-monodromy-eigenvalue).
The primitive indicial multiset modulo \(k=6-3=3\)
is \(A_3D\) modulo three, since \(F_0=A_6D\) reduces to
twice \(A_3D\). On a noninteger orbit of mass one it
has one occurrence each of \(\eta,\eta+1,\eta+2\)
modulo three. The congruence would require
\(3\eta+3\equiv0\pmod3\), which forces \(\eta\) integral.
No orbit of mass two outside the integers is possible,
since the integer orbit already has positive mass.
A zero-mass orbit cannot contain a nonzero contribution:
its \(F_0\) part would be effective of mass zero, forcing
its \(D\) part to vanish by cancellation of \(A_6\).

Polynomiality of the three Euler components makes
\(A_\ell(D-[0])\) effective for \(\ell=6,17,21\).
In fact an Euler term \(x^{-\ell}p(E)\) in the original
polynomial Weyl algebra requires \((E)_\ell\mid p(E)\).
Put \(C_k\) equal to the cumulative sum of \(D-[0]\)
through the integer \(k\). The effectivity conditions
are \(C_k\ge C_{k-\ell}\) for these three steps.
Along each residue class modulo six, the sequence is
nondecreasing with end values zero and one, so \(C_k\)
is zero or one. The set where it is one is closed under
addition of the three steps.
After translating its smallest index to zero, its
complement is precisely a gap subset closed under the
stated subtractions. The mass of \(D\) is two, and
the exact mean of \(f\) gives
\(\operatorname{moment}(F_0)=12\cdot5/2=30
=6\operatorname{moment}(D)+2(0+\cdots+5)\).
Thus \(D-[0]\) has first moment zero.
In the translated sequence with smallest index zero,
the first moment of the successive differences is
the number of omitted nonnegative indices, by telescoping.
The required translation is therefore \(-|I|\), giving (1).

The [exact gap checker](verify_weyl_forty_two_indicial.py)
finds the twenty-five gaps of \(\mathcal S_0\) and enumerates
the downward-closed subsets under the three predecessors.
This yields \(1694\) distinct normalized choices. The same normalized closed sets are
the numerical sets for simultaneous \((6,17,21)\)-core
partitions.
[Constantin, Houston-Edwards and Kaplan, Proposition 2.4](https://arxiv.org/pdf/1509.06077)
The sector and mean constraints were used over \(\mathbb C\).
For characteristic zero they apply by embedding the finitely
generated coefficient field, with the selected roots, into
\(\mathbb C\).


## The number of small roots

**Corollary PolynomialWeylOrderTwelveFortyTwoRootCount.**
In the reduced \((12,42)\) row let \(n_H,n_L\) be the
numbers of roots of multiplicities fifteen and three.
Then \(n_L\le n_H+2\).
If \(n_H=1\), the degree is twenty-four and the common
leading factor has the form
\[
             h=\gamma (x-x_0)^{15}
                    \bigl((x-x_0)^3-\lambda\bigr)^3,
                    \qquad \gamma\lambda\ne0 .
\]

Use the
[rational-map form of the leading-factor residue](POLYNOMIAL_WEYL_LEADING_FACTOR_RESIDUE.md#a-rational-map-from-the-residue-equation).
Here \(M=(5n_H+n_L)/2\), and its divisor exponent gcd
is three. If \(n_L=0\), the extra roots of the residue
polynomial still contribute six to that gcd, since
\(M>n_H\). The residue map has degree \(N=3n_H\) and
ramification index \(M-1\) at infinity; hence \(N\ge M-1\)
gives the inequality.

For \(n_H=1\), divisibility of the degree by six and its
lower bound twenty-four force \(n_L=3\).
The residue map then has a single pole of order three
at \(x_0\), and ramification index three at infinity.
As a polynomial in \(u=(x-x_0)^{-1}\), it must have
the form \(a u^3+b\), where \(ab\ne0\).
Its three simple zeros give the three low roots in
the displayed form of \(h\).
