# A nontorsion branch point independent of the rational elliptic group

The degree-two zero of \(F(x)=128x-25865\) on the existing elliptic
slice is nontorsion. Over its quadratic field it contributes a new
Mordell--Weil direction, independent of the whole rational group.
This gives a fixed branch-point input for a future divisibility
argument; it does not establish an all-primes measure or an infinite
Hasse subfamily.

## 1. The fixed point and an explicit quadratic twist

Retain the models
\[
 E_0:\ y^2=x^3-10656x^2+27091730x,\qquad X=x-3552,
\]
\[
 E:\ y^2=X^3-10758382X+6600759744.
\]
Identify these two models by the displayed translation.
Put \(d=51730=2\cdot5\cdot7\cdot739\) and
\(L=\mathbb Q(\sqrt d)\). The branch identity from
[EHD1](ELLIPTIC_EXCEPTION_AMBIENT_HAAR_DIAGNOSTIC.md) gives
\[
 R_\pm=\left(\frac{25865}{128},
               \ \pm\frac{639735\sqrt d}{2048}\right)\in E_0(L),
 \qquad \operatorname{div}(128x-25865)=R_++R_--2O.   \tag{1.1}
\]
The nontrivial \(\sigma\in\operatorname{Gal}(L/\mathbb Q)\)
sends \(R_+\) to \(R_-=-R_+\).

The short quadratic twist is
\[
 E^{(d)}:\quad
 Y^2=Z^3-10758382\,d^2 Z+6600759744\,d^3.
\]
Its isomorphism over \(L\) to \(E\) is
\[
 (Z,Y)\longmapsto
       \left(\frac Zd,\frac{Y}{d\sqrt d}\right).
\]
The rational point mapping to \(R_+\) in short coordinates is
\[
 Q=\left(
 d\left(\frac{25865}{128}-3552\right),
 \frac{639735d^2}{2048}\right)\in E^{(d)}(\mathbb Q). \tag{1.2}
\]
Its curve identity follows identically by multiplying the branch
identity of (1.1) by \(d^3\); no point search is involved.

## 2. A torsion-free formal neighborhood

**Theorem EBT1 (nontorsion branch and rank increase).**
The point \(Q\) is nontorsion, hence so are \(R_+,R_-\).
Moreover
\[
 \boxed{\operatorname{rank}E^{(d)}(\mathbb Q)\ge1,\qquad
        \operatorname{rank}E(L)\ge3.}               \tag{2.1}
\]
More precisely, if \(n\in\mathbb Z\), \(P\in E(\mathbb Q)\),
and \(nR_++P\) is torsion in \(E(L)\), then \(n=0\) and \(P\)
is torsion.

**Proof of nontorsion.** The short twist has integral coefficients.
For the coordinates of (1.2),
\[
 v_2(Z(Q))=-6,\qquad v_2(Y(Q))=-9,\qquad
              v_2(-Z(Q)/Y(Q))=3.                  \tag{2.2}
\]
Thus \(Q\) belongs to the formal subgroup \(E^{(d)}_3(\mathbb Q_2)\)
and has exact depth three.

For an integral short Weierstrass equation, Milne's parameter estimate
gives
\[
 t(U+V)=t(U)+t(V)+O(2^{5n})\qquad(U,V\in E_n(\mathbb Q_2)).
\]
Doubling a point of exact depth \(n\ge1\) therefore gives exact
depth \(n+1\), since \(n+1<5n\). Odd integer multiplication preserves
its depth. Consequently every nonzero integer multiple of \(Q\)
has finite depth, and cannot equal \(O\). Hence \(Q\) is nontorsion.
The isomorphism over \(L\) proves the same for \(R_\pm\).

The exact formal-group inputs are Milne, *Elliptic Curves* (2006),
Chapter II, Theorem 4.1 and Corollary 4.2, printed pp.62--63, and
Lemma 5.9, printed p.67:
[primary book](https://www.jmilne.org/math/Books/ectext6.pdf),
[retained copy](sources/milne-elliptic-curves-2006.pdf).
This proof uses the integral short twist over \(\mathbb Q_2\);
it does not require extending the estimate to fractional valuations
in \(L\).

**Proof of independence.** Suppose \(nR_++P\) is torsion. Its conjugate
\(-nR_++P\) is torsion as well, so their difference \(2nR_+\)
is torsion. Since \(R_+\) is nontorsion, \(n=0\), and then \(P\)
is torsion. The [known rational rank is two](ELLIPTIC_EXCEPTION_ARITHMETIC.md),
so adjoining this independent direction proves the second inequality
in (2.1).
The rational nontorsion point \(Q\) proves the first. \(\square\)

## 3. Scope of this branch-point input

The result concerns the same fixed elliptic curve, its fixed branch
field and one fixed quadratic twist. It does not produce another
moving-factor family or identify a basis of either Mordell--Weil
group. The ranks in (2.1) are lower bounds, not exact ranks.

The independent point supplies the nondegeneracy input to
[the Kummer branch-exclusion theorem](ELLIPTIC_EXCEPTION_KUMMER_BRANCH_EXCLUSION.md).
It gives a positive-density set of inert primes at which neither
finite branch residue belongs to the rational reduction image, even
though all four quadratic classes pass the earlier eligibility tests.
The local sieve there has only its pole branch. Neither result gives
a lower bound for the actual conditional prime losses of
[EFA4](ELLIPTIC_EXCEPTION_FINITE_LOCAL_ASYMPTOTICS.md#6-conditional-prime-losses-and-the-limiting-measure)
nor a value of the intrinsic all-primes Haar measure.

The coordinate and valuation identities are exact rational arithmetic.
The nontorsion and independence conclusions follow from the written
formal-group and conjugation arguments. There is no new finite search,
formal verification or external review.
