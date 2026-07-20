# Exact finite-field value distribution

Let

\[
F(x,y,z)=\left(
(1+xy)^3z+y^2(1+xy)(4+3xy),
y+3x(1+xy)^2z+3xy^2(4+3xy),
2x-3x^2y-x^3z
\right).
\]

For a finite field \(k=\mathbb F_q\), write \(N_j\) for the number of
targets in \(k^3\) whose fiber in \(k^3\) has cardinality \(j\).

## Theorem (characteristic different from 2 and 3)

If \(\operatorname{char} k>3\), every fiber has cardinality \(0,1\), or
\(3\), and

\[
N_3=\frac{(q-1)(q^2+2)}6,\qquad
N_0=\frac{(q-1)(q^2+2)}3,
\]

\[
N_1=\frac{q^3+q^2-2q+2}{2}.
\]

Consequently

\[
|F(k^3)|=N_1+N_3
=\frac{2q^3+q^2-2q+2}{3},
\]

and the image density is

\[
\frac{|F(k^3)|}{q^3}
=\frac23+\frac1{3q}-\frac2{3q^2}+\frac2{3q^3}
\longrightarrow \frac23.
\]

### Proof

Write a target as \((a,b,c)\). On the chart \(x\ne0\), put
\(t=y+1/x\). The inverse equation and reconstruction identity are

\[
P_{a,b,c}(T)=cT^3-2T^2+bT-2a=0,
\qquad P'_{a,b,c}(t)=\frac2x.
\]

Thus every root which comes from a finite source point is simple. If
\(c\ne0\), every simple root in \(k\) reconstructs a unique source point,
whereas a repeated root has \(P'=0\) and gives no finite point. A cubic has
\(0,1\), or \(3\) simple roots in its ground field (a repeated-root cubic has
at most one simple root), so these are the only possible fiber sizes over
\(c\ne0\).

A three-point fiber over \(c\ne0\) is determined by an unordered triple
\(\{t_1,t_2,t_3\}\) of distinct elements. Comparing the coefficient of
\(T^2\) gives

\[
c(t_1+t_2+t_3)=2.
\]

Hence such targets are in bijection with the three-element subsets of \(k\)
having nonzero sum. Conversely, their nonzero sum determines \(c\), and the
remaining elementary symmetric functions determine \(b\) and \(a\).

If \(c=0\), the plane \(x=0\) contributes exactly one point to every fiber,
while the other chart is governed by the quadratic

\[
-2T^2+bT-2a.
\]

The fiber has three points precisely when this quadratic has two distinct
roots. These targets are therefore in bijection with the unordered pairs of
distinct elements of \(k\), of which there are \(\binom q2\).

It remains to count zero-sum triples. In characteristic different from \(3\),
the number of ordered triples of distinct elements satisfying
\(r+s+t=0\) is

\[
q^2-(3q-2)=(q-1)(q-2).
\]

Division by \(6\) gives \((q-1)(q-2)/6\) unordered triples. Therefore

\[
N_3=\binom q2+\binom q3-\frac{(q-1)(q-2)}6
=\frac{(q-1)(q^2+2)}6.
\]

Since \(\det DF=-2\ne0\), the preceding description accounts for every
rational point and every fiber has size \(0,1\), or \(3\). Counting targets
and then source points gives

\[
q^3=N_0+N_1+N_3,
\qquad q^3=N_1+3N_3.
\]

Thus \(N_0=2N_3\), and the remaining formulas follow. This proof uses only
finite-field arithmetic and therefore applies without change to every
\(\mathbb F_{p^m}\) with \(p>3\).

## Characteristic 3

The map is still etale because \(\det DF=-2=1\). The fiber sizes are again
\(0,1,3\), but the histogram is

\[
N_3=\frac{q^2(q-1)}6,\qquad
N_0=\frac{q^2(q-1)}3,\qquad
N_1=\frac{q^2(q+1)}2,
\]

and hence

\[
|F(k^3)|=\frac{q^2(2q+1)}3.
\]

Indeed, in characteristic \(3\), if two entries of a zero-sum triple agree,
then all three agree. There are therefore \(q^2-q=q(q-1)\) ordered distinct
zero-sum triples, or \(q(q-1)/6\) unordered ones. Substitution in the preceding
argument gives the formulas. In particular, over \(\mathbb F_3\) the
histogram is

\[
(N_0,N_1,N_3)=(6,18,3).
\]

The image density is \(2/3+1/(3q)\), so it has the same limit \(2/3\).

## Characteristic 2

Here \(\det DF=0\), and the map is not etale. For \(q=2^m\), its only fiber
sizes are \(0,1,q-1,2q-1\), subject to coincident labels when \(q=2\), and

\[
N_0=2(q-1)^2,
\]

\[
N_1=q^3-2q^2+2q-1,\qquad
N_{q-1}=2(q-1),\qquad N_{2q-1}=1.
\]

For \(q=2\), the \(N_{q-1}\) targets are also singleton fibers, so the actual
histogram is \(N_0=2,N_1=5,N_3=1\). For \(q=4\), \(q-1=3\), but there is no
additional three-point class: the displayed classes already give the whole
histogram.

To prove the claim, retain \(t=y+1/x\) on \(x\ne0\). Direct simplification in
characteristic \(2\) gives

\[
b=ct^2,
\qquad a+tb+t^2=\frac{t}{x}.
\]

Because squaring is bijective on \(k\), these two equations classify all
points on this chart. On \(x=0\), one has simply

\[
F(0,y,z)=(z,y,0),
\]

so every target with \(c=0\) has one additional point.

For \(c\ne0\), the targets \((0,0,c)\) have \(q-1\) preimages; for every
nonzero \(b\), exactly one value of \(a\) is omitted and all other values have
one preimage; and \((a,0,c)\) is omitted when \(a\ne0\). For \(c=0\), the
origin has \(2q-1\) preimages, each \((a,0,0)\) with \(a\ne0\) has \(q-1\),
and all \((a,b,0)\) with \(b\ne0\) have one. Summing these cases gives the
stated histogram and

\[
|F(k^3)|=q^3-2(q-1)^2.
\]

Thus the characteristic-two image density tends to \(1\), rather than
\(2/3\).

## Four refinements in characteristic greater than 3

The elementary proof contains more information than the aggregate histogram.
The following refinements make the boundary contribution and the \(S_3\)
factorization law explicit.

### 1. The strata \(c=0\) and \(c\ne0\)

On the target plane \(c=0\), every target has its unique source point on
\(x=0\). The quadratic on the other chart either contributes two points or no
points. Consequently

\[
N^{c=0}_0=0,\qquad
N^{c=0}_1=\frac{q(q+1)}2,\qquad
N^{c=0}_3=\frac{q(q-1)}2.
\]

All omitted values therefore lie in \(c\ne0\). On that stratum,

\[
N^{c\ne0}_3=\frac{(q-1)^2(q-2)}6,
\]

\[
N^{c\ne0}_1=\frac{(q-1)^2(q+2)}2,
\qquad
N^{c\ne0}_0=\frac{(q-1)(q^2+2)}3.
\]

For the point count used here, fix \(c\ne0\). The source equation for the
third coordinate has exactly one \(z\) for each \(x\ne0\) and \(y\), hence
\(q(q-1)\) source points above the target plane with that fixed \(c\). Combining
this with the three-fiber count determines the singleton and empty counts.

### 2. Complete cubic factorization law

For fixed \(c\ne0\), normalization makes \(P_{a,b,c}\) a monic cubic whose
root sum is the fixed nonzero element \(2/c\). As \((a,b)\) ranges over
\(k^2\), its exact factorization distribution is

| Factorization over \(k\) | Number for each fixed \(c\ne0\) | Finite fiber |
|---|---:|---:|
| three distinct linear factors \(1+1+1\) | \((q-1)(q-2)/6\) | 3 |
| one linear and one irreducible quadratic \(1+2\) | \(q(q-1)/2\) | 1 |
| one double and one simple linear factor \(2+1\) | \(q-1\) | 1 |
| one triple linear factor \(3\) | \(1\) | 0 |
| irreducible cubic | \((q^2-1)/3\) | 0 |

For the \(2+1\) row, choose the double root \(r\); the simple root is forced
by the prescribed sum, and precisely one value of \(r\) would make all three
roots equal. For the \(1+2\) row, choose the linear root and then a monic
irreducible quadratic with prescribed trace. The remaining two rows follow by
subtraction, or from the standard counts with prescribed nonzero trace.

A repeated root contributes no finite source point because \(P'(t)=0\) would
force \(2/x=0\). This explains why a \(2+1\) factorization gives only the one
simple-root source, while the triple-root polynomial gives an empty fiber.

### 3. Discriminant square classes

For every fixed \(c\), including \(c=0\), the \(q^2\) targets \((a,b,c)\)
split as

\[
\#\{\operatorname{Disc}P=0\}=q,
\]

\[
\#\{\operatorname{Disc}P\text{ a nonzero square}\}
=\#\{\operatorname{Disc}P\text{ a nonsquare}\}
=\frac{q(q-1)}2.
\]

For \(c\ne0\), a separable cubic has square discriminant exactly when its
Frobenius permutation is even. Thus the square class consists of the split
\(1+1+1\) cubics and the irreducible cubics (identity and 3-cycle), while the
nonsquare class consists of the \(1+2\) cubics (transpositions). The zero class
is the union of the \(2+1\) and \(3\) rows. The table above gives the displayed
equalities immediately.

For \(c=0\),

\[
\operatorname{Disc}P=4(b^2-16a),
\]

which runs uniformly through \(k\) as \(a\) varies. A nonzero square gives the
two quadratic roots and hence a three-point fiber; a nonsquare or zero gives
only the boundary point.

### 4. Double and triple fiber products

Let

\[
X^{[2]}=\{(u,v)\in(k^3)^2:F(u)=F(v),\ u\ne v\}
\]

and let \(X^{[3]}\) be the analogous space of ordered, pairwise distinct
triples in one fiber. Then

\[
\#X^{[2]}(k)=\sum_y(|F^{-1}(y)|)_2=6N_3
=(q-1)(q^2+2),
\]

\[
\#X^{[3]}(k)=\sum_y(|F^{-1}(y)|)_3=6N_3
=(q-1)(q^2+2).
\]

The equality is structural: an ordered pair of distinct points can occur only
in a three-point fiber, where it determines the remaining third point. The
full double and triple fiber products, with diagonals retained, therefore have

\[
\#(k^3\times_F k^3)=q^3+(q-1)(q^2+2),
\]

\[
\#(k^3\times_F k^3\times_F k^3)
=N_1+27N_3=q^3+4(q-1)(q^2+2).
\]

Because the off-diagonal formula holds over every \(\mathbb F_{p^m}\), its
local zeta function is explicitly

\[
Z(X^{[2]}/\mathbb F_p,T)=Z(X^{[3]}/\mathbb F_p,T)
=\frac{(1-p^2T)(1-T)^2}{(1-p^3T)(1-pT)^2}.
\]

The limiting normalized fiber law is consequently

\[
(\Pr(0),\Pr(1),\Pr(3))\longrightarrow
\left(\frac13,\frac12,\frac16\right),
\]

exactly the fixed-point distribution of the natural action of \(S_3\).

## Context

The result is an exact value distribution for a multivariate polynomial map,
not merely an asymptotic obtained from monodromy. Its elementary explanation is
the cubic inverse model plus the affine boundary chart. It therefore provides a
concrete comparison point for general work on value sets of polynomial maps and
for current work deriving finite-field distributions from monodromy.

Two useful points of comparison are:

- Gary L. Mullen, Daqing Wan, and Qiang Wang, [*Value sets of polynomial maps
  over finite fields*](https://arxiv.org/abs/1210.8119), which treats
  multivariate polynomial maps and general value-set bounds.
- David Kumallagov, [*Affine monodromy and exact value distributions over
  finite fields*](https://arxiv.org/abs/2607.09799), which computes complete
  fiber laws for an explicit family of one-variable covers using monodromy.
