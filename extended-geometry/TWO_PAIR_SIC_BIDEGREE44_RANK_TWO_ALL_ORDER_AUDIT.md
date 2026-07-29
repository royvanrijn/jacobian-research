# All-order audit of the bidegree-\((4,4)\) rank-two moment fiber

## 1. Outcome

The proposed all-order calculation cannot yet be specialized to a known
rank-two survivor.  The current exact input is weaker:

1. the first thirteen moments have a semistable common zero on the
   rank-at-most-two determinantal variety;
2. one squarefree rank-one Rabinowitsch chart remains open, so exact
   coefficient rank two has not been forced; and
3. the Hilbert-series proof is existential and records neither coordinates
   nor a residue field for its semistable point.

Consequently there is no exact coefficient point at which to derive and
solve a scalar recurrence.  Neither an all-order rank-two witness nor an
exact tail obstruction is proved here.  The first moment not imposed by the
known prefix is \(\mu _{14}\), but its value on the existential fiber is
unknown.

This corrects a possible misreading of equation (6.1) in the
[rank frontier](TWO_PAIR_SIC_BIDEGREE44_RANK_FRONTIER.md).  That displayed
rank-two matrix is a transversality point used to prove Jacobian rank
thirteen.  It is not a moment zero:
\[
\begin{aligned}
 \mu _1&=7414,\\
 \mu _2&=3675739680,\\
 \mu _3&=12167497410877440,\\
 \mu _4&=148010006143680629760000.
\end{aligned} \tag{1.1}
\]

## 2. Exact rank-two period

The determinantal parameterization can nevertheless be put into the form
needed by creative telescoping.  Write
\[
 C=UW,\qquad
 U=(u_{iq})\in\operatorname {Mat}_{5\times2},\qquad
 W=(w_{qj})\in\operatorname {Mat}_{2\times5}. \tag{2.1}
\]
With
\[
 \Phi_C(x,y)=\sum_{i,j=0}^4c_{ij}x^iy^j,
\]
define the Laurent polynomial
\[
\begin{aligned}
 A_q(u)&=\sum_{i=0}^4u_{iq}u^{4-i},\\
 B_q(u,t)&=\sum_{j=0}^4
   w_{qj}u^{j-4}t^j(1-t)^{4-j},\\
 P_{U,W}(u,t)&=A_1(u)B_1(u,t)+A_2(u)B_2(u,t).
\end{aligned} \tag{2.2}
\]
Direct multiplication gives
\[
 P_{U,W}(u,t)
 =\Phi_C\left(1,u,t,\frac{1-t}{u}\right). \tag{2.3}
\]
The formal beta identity from the full-rank witness therefore gives
\[
 \boxed{\frac{\mu_m(C)}{(4m+1)!}
 =\operatorname {CT}_u\int_0^1P_{U,W}(u,t)^m\,dt.} \tag{2.4}
\]
Equivalently, the factorial-normalized ordinary generating function is
\[
 \boxed{
 {\cal G}_{U,W}(s)
 =\sum_{m\geq0}\frac{\mu_m(C)}{(4m+1)!}s^m
 =\operatorname {CT}_u\int_0^1
   \frac{dt}{1-sP_{U,W}(u,t)}.} \tag{2.5}
\]

Formula (2.5) is the exact holonomic starting point.  Closure of holonomic
functions under coefficient extraction and definite integration implies
that every exact algebraic specialization of \(U,W\) has a P-finite
moment sequence.  This existence statement is not a recurrence
certificate: an all-order result still requires an explicit scalar
operator, telescoping certificates including the endpoint terms, its
singular-step audit, and enough exact initial values.

## 3. Rank two does not supply a small generic cutoff

For a generic rank-two factor point, the exponent support of (2.2) has
Newton polygon
\[
 \operatorname {conv}\{(-4,0),(0,0),(4,4),(-4,4)\}. \tag{3.1}
\]
Its Euclidean area is \(24\), hence its normalized two-dimensional volume
is \(48\).  Thus the determinantal rank condition does not by itself
collapse the Laurent support to the six-point, normalized-volume-eight
support of the known full-rank witness.

The number \(48\) is not asserted to be the scalar recurrence order.
The period realization is relative to the endpoints \(t=0,1\), and the
coefficient scaling curve passes through the zero polynomial at \(s=0\).
The ordinary-point and scalar cyclic-vector gates from
[the holonomic algorithm note](HOLONOMIC_HYPERGEOMETRIC_ALGORITHMS.md)
therefore remain necessary.  In particular, the thirteen zero moments do
not propagate merely from holonomicity or from (3.1).

## 4. Exact gates before recurrence solving

The required order is:

1. certify the remaining squarefree rank-one target-only membership
   \[
   \lambda^4(\lambda-1)^4
   \bigl(p(8c-3d^2)\bigr)^5
   \in(f_3,f_4,f_5,f_6);
   \]
   the exponent \(5\) is the common least exponent modulo
   \(101,103,107\), but the rational lift remains open;
2. compute an explicit closed point, or an explicit positive-dimensional
   component and its function field, in the semistable
   \(\mu_1=\cdots=\mu_{13}=0\) fiber;
3. verify exact coefficient rank two at that point;
4. specialize (2.5) and produce a checkable creative-telescoping
   recurrence; and
5. evaluate \(\mu_{14}\) and the later bridge values required by the
   recurrence.

If the recurrence propagates zero, a fixed mixed multiplier must still be
tested before the point becomes an SIC witness.  If a tail moment is
nonzero, that exact value is the desired obstruction for that closed
point; excluding rank two globally would require treating every
semistable component.

## 5. Reproduction

Run

```bash
python3 scripts/verify_two_pair_sic_bidegree44_rank_two_all_order_audit.py
```

The dependency-free checker verifies the factor identity (2.3), the
beta/constant-term identity (2.4) through order four at the displayed
exact rank-two chart point, the four nonzero values in (1.1), and the
Newton polygon and normalized volume in (3.1).  The finite replay is an
audit of the formulas and of the proposed starting point; it is not an
all-order recurrence certificate.
