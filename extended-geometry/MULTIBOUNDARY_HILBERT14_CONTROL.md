# A multiboundary Hilbert--14 control

## Result

The one-divisor saturation ladders on the normalized \((2,3)\) and
\((2,4)\) factorization slices both terminate.  This note supplies the
missing two-divisor comparison: an explicit finitely generated domain with
two commuting locally nilpotent derivations whose common invariant ring is
not finitely generated.

The example is a multiboundary control, not a Keller map.  It identifies the
precise conductor and bidegree pattern that a Keller-attached construction
would need to reproduce.

## 1. Why the obvious factorization boundary cannot work

For unequal-degree binary factors, write the top coefficients as

\[
 A=aT^r+bT^{r-1}S+\cdots,\qquad
 B=pT^q+dT^{q-1}S+\cdots.
\]

The tangent normalization used in the factorization slices is

\[
 [AB]_{T^{r+q-1}S}=ad+bp=1.                         \tag{1}
\]

Hence

\[
 1=ad+bp\in(a,p).
\]

The two leading divisors \(a=0\) and \(p=0\) are therefore disjoint on every
such normalized slice.  They do not define a normal-crossings corner and
cannot support a genuine bifiltration
\(C_{i,j}a^{-i}p^{-j}\).  Merely increasing the factor degrees does not
create the multiboundary geometry sought here.

This is compatible with the independent obstruction in
[Controlled boundary suspensions](../cancellation/CONTROLLED_BOUNDARY_SUSPENSIONS.md#4-a-first-independent-two-boundary-ansatz):
the smallest one-reconstruction-variable two-boundary ansatz acquires an
unwanted third divisor.

## 2. Two cusp boundaries and two additive actions

Put

\[
 A_0=k[s^2,s^3,t^2,t^3],\qquad
 \mathcal R=A_0[X,Y,U,V].
\]

Its normalization is

\[
 \widetilde{\mathcal R}=k[s,t,X,Y,U,V].
\]

Define two commuting locally nilpotent derivations

\[
 D_s=s^3{\partial\over\partial X}
     -s^2{\partial\over\partial Y},\qquad
 D_t=t^3{\partial\over\partial U}
     -t^2{\partial\over\partial V}.                  \tag{2}
\]

They preserve \(\mathcal R\), and in the normalization they fix

\[
 P=X+sY,\qquad Q=U+tV.                              \tag{3}
\]

Since multiplication by the nonzero factors \(s^2,t^2\) does not change a
derivation kernel in a domain,

\[
 \widetilde{\mathcal R}^{D_s,D_t}=k[s,t,P,Q].
\]

Thus the common invariant ring under study is the explicit intersection

\[
 \mathcal K
 =\mathcal R^{D_s,D_t}
 =k[s,t,P,Q]\cap\mathcal R.                          \tag{4}
\]

## 3. The bivariate invariant ladder

For every \(m,n\geq0\), set

\[
 F_{m,n}=s^2t^2P^mQ^n.                              \tag{5}
\]

These functions are killed by both derivations.  They lie in the original
ring because

\[
\begin{aligned}
F_{m,n}
={}&\sum_{i=0}^m\sum_{j=0}^n
 {m\choose i}{n\choose j}
 s^{2+i}t^{2+j}
 X^{m-i}Y^iU^{n-j}V^j,
\end{aligned}
\]

and every exponent \(2+i\) and \(2+j\) belongs to the numerical semigroup
\(\langle2,3\rangle\).

The two axes have conductor ideals

\[
 I_s=(s^2,s^3)\mathcal R,\qquad
 I_t=(t^2,t^3)\mathcal R.                            \tag{6}
\]

Give the normalization the bidegree

\[
 \deg P=(1,0),\qquad \deg Q=(0,1).
\]

Every bihomogeneous common invariant of bidegree \((m,n)\) has the form
\[
 h(s,t)P^mQ^n.
\]
If \(m>0\), its expansion contains the pair of coefficients
\[
 h(s,t)X^m,\qquad m\,s\,h(s,t)X^{m-1}Y.
\]
The cusp ring contains no term of \(s\)-order one.  Therefore the least
\(s\)-order of \(h\) cannot be zero, and the invariant belongs to \(I_s\).
Similarly,

\[
 m>0\Longrightarrow F\in I_s,\qquad
 n>0\Longrightarrow F\in I_t.                       \tag{7}
\]

This is the two-axis version of Maubach's conductor lemma.

## 4. Uniform non-finite-generation proof

Reduce \(\mathcal K\) modulo

\[
 J=(s^4,t^4)\mathcal R.                             \tag{8}
\]

By (6)--(7), a product of two invariants of positive \(P\)-degree lies in
\(I_s^2\subset(s^4)\) and vanishes modulo \(J\).  Likewise a product of two
positive-\(Q\)-degree invariants vanishes modulo \(J\).

Suppose \(\mathcal K\) had finitely many bihomogeneous generators.  Let
\(d_P\) and \(d_Q\) be the largest \(P\)- and \(Q\)-degrees among them.  A
nonzero monomial in those generators modulo \(J\) can use at most one factor
of positive \(P\)-degree and at most one factor of positive \(Q\)-degree.
Every surviving bidegree therefore lies in the finite rectangle

\[
 [0,d_P]\times[0,d_Q].
\]

But

\[
 F_{d_P+1,d_Q+1}
\equiv s^2t^2X^{d_P+1}U^{d_Q+1}+\cdots\not\equiv0
\pmod J,
\]

and has bidegree \((d_P+1,d_Q+1)\).  This is a contradiction.  Hence

\[
 \boxed{\mathcal R^{D_s,D_t}\text{ is not finitely generated}.} \tag{9}
\]

Unlike a bounded search, (9) is an arbitrary-bidegree proof.  The finite
checker only replays a configurable rectangle of the uniform identities.

## 5. Consequence for the Keller search

The comparison now separates the observed behaviors cleanly.

| model | boundary geometry | invariant outcome |
|---|---|---|
| normalized \((2,3)\) slice | one effective leading divisor | finite |
| normalized \((2,4)\) slice | one effective leading divisor | finite |
| tangent leading pair \(a=0,p=0\) | disjoint by (1) | no bifiltration |
| two-cusp control | intersecting conductor axes | non-finite |

The next Keller-attached candidate must therefore contain two intersecting
normalization/conductor primes.  The current direct two-boundary suspension
with one reconstruction variable is unsuitable because it introduces a
third divisor.  The concrete construction choices left by the existing
boundary audit are:

1. use two reconstruction variables;
2. accept the third divisor and run a three-index saturation ledger; or
3. extract an LND action on a completed two-boundary normalization chart
   already present in the plane boundary compiler.

In each case the decisive test is now explicit: reduce modulo the squares of
the two conductor ideals and look for invariant bidegrees escaping every
finite rectangle.

Run

```bash
.venv/bin/python scripts/verify_multiboundary_hilbert14_antichain.py
```

to check the commuting LNDs, the cusp-semigroup expansions, the
modulo-\((s^4,t^4)\) product vanishings, the finite rectangle replay, and the
factorization leading-divisor obstruction.
