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

## 5. The exact finite-generation ideal and SAGBI ladder

For a \(k\)-subalgebra \(E\) of a finitely generated domain, its
finite-generation ideal is

\[
 \mathfrak f_E
 =\{0\}\cup\{g\in E\setminus\{0\}:E_g
                   \text{ is a finitely generated }k\text{-algebra}\}.
                                                               \tag{10}
\]

This is a radical ideal.  It was introduced as an algorithmic invariant by
Derksen--Kemper, and recent Hilbert--14 calculations use conductors and
infinite SAGBI bases to determine it explicitly; see
[Derksen--Kemper](https://sites.lsa.umich.edu/hderksen/wp-content/uploads/sites/614/2018/05/A.I.a.28.pdf)
and
[Hart](https://arxiv.org/abs/2203.15569).

The present control admits a particularly transparent exact computation.
First consider one cusp.  Inside \(C_s=k[s,P]\), put

\[
 K_s=k+s^2C_s,\qquad \mathfrak m_s=s^2C_s.          \tag{11}
\]

The homogeneous coefficient test in Section 3 gives

\[
 k[s,P]\cap k[s^2,s^3,X,Y]=K_s.
\]

Moreover, \(\mathfrak m_s=[K_s:C_s]\) is the conductor.  If
\(0\ne g\in\mathfrak m_s\), then
\((K_s)_g=(C_s)_g\), so \(g\in\mathfrak f_{K_s}\).  Conversely, write
\(g=\lambda+h\) with \(\lambda\in k^\times\) and
\(h\in\mathfrak m_s\).  Modulo \(\mathfrak m_s^2=s^4C_s\), \(g\) is a unit
with inverse

\[
 g^{-1}\equiv\lambda^{-1}-\lambda^{-2}h
       \pmod{\mathfrak m_s^2}.                      \tag{12}
\]

Products of two positive-\(P\)-degree classes vanish in this quotient.
Thus any finite list of generators of \((K_s)_g\) has a bounded surviving
\(P\)-degree, while \(s^2P^d\) survives for every \(d\).  Hence

\[
 \boxed{\mathfrak f_{K_s}=\mathfrak m_s=s^2k[s,P].} \tag{13}
\]

The two-cusp ring factors as

\[
 \mathcal K=K_s\otimes_kK_t,\qquad
 K_t=k+t^2k[t,Q].                                  \tag{14}
\]

Its conductor to the invariant algebra
\(C=k[s,t,P,Q]\) on the normalized ambient ring is

\[
 \mathfrak c=[\mathcal K:C]
 =\mathfrak m_s\otimes_k\mathfrak m_t
 =s^2t^2C.                                         \tag{15}
\]

Every nonzero \(g\in\mathfrak c\) makes
\(\mathcal K_g=C_g\), so \(\mathfrak c\subset\mathfrak f_{\mathcal K}\).
For the converse, suppose \(\mathcal K_g\) is finitely generated.  Specialize
the \(t,Q\) factor at a \(k\)-point for which the image of \(g\) modulo
\(\mathfrak m_s\) is nonzero.  Such a point exists for every nonzero
polynomial because \(k\) is infinite.  The resulting quotient is a
localization of \(K_s\) at an element outside \(\mathfrak m_s\), contradicting
(13).  Therefore
\(g\in\mathfrak m_s\otimes K_t\).  Interchanging the two factors gives
\(g\in K_s\otimes\mathfrak m_t\), and the intersection of these two
subspaces is (15).  Consequently

\[
 \boxed{\mathfrak f_{\mathcal K}
        =[\mathcal K:k[s,t,P,Q]]
        =s^2t^2k[s,t,P,Q].}                        \tag{16}
\]

Thus the finite-generation locus is exactly the locus \(D(\mathfrak c)\)
where the comparison \(\mathcal K\subset C\) becomes an isomorphism.  This
does not assert that the generally nonfinite extension
\(\mathcal K\subset C\) is an integral normalization map.  The bad locus
\(V(\mathfrak c)\) is the union of the two pinched axes, not merely their
set-theoretic corner.

There is also a literal SAGBI interpretation.  In the normalization
coordinates the algebra is monomial, with infinite SAGBI basis

\[
\begin{aligned}
 \mathcal G_s={}&\{s^2,s^3\}
 \cup\{s^2P^m,s^3P^m:m\ge1\},\\
 \mathcal G_t={}&\{t^2,t^3\}
 \cup\{t^2Q^n,t^3Q^n:n\ge1\}.                     \tag{17}
\end{aligned}
\]

Their union is a SAGBI basis of \(\mathcal K\).  The finite-generation ideal
is generated as an ideal by the four return ladders

\[
 s^\epsilon t^\delta P^mQ^n,\qquad
 \epsilon,\delta\in\{2,3\},\quad m,n\ge0.           \tag{18}
\]

In particular the displayed \(F_{m,n}\) are not just witnesses to
non-finite generation: they are indispensable low-conductor-order
generators of \(\mathfrak f_{\mathcal K}\).  No finite subfamily generates
that ideal, because multiplication by a positive \(P\)- or \(Q\)-degree
element of \(\mathcal K\) raises the corresponding conductor order.

## 6. Arbitrarily many independent cusp boundaries

The two-axis calculation is not exceptional.  For \(r\ge1\), put

\[
 C_i=k[t_i,P_i],\qquad
 K_i=k+t_i^2C_i,\qquad
 \mathfrak m_i=t_i^2C_i,
\]

and define

\[
 K^{(r)}=\bigotimes_{i=1}^rK_i
 \subset
 C^{(r)}=k[t_1,P_1,\ldots,t_r,P_r].                 \tag{19}
\]

Then

\[
 \boxed{
 \mathfrak f_{K^{(r)}}
 =[K^{(r)}:C^{(r)}]
 =\bigotimes_{i=1}^r\mathfrak m_i
 =\left(\prod_{i=1}^rt_i^2\right)C^{(r)}.
 }                                                   \tag{20}
\]

Indeed every nonzero element of the last ideal makes
\(K^{(r)}\subset C^{(r)}\) an equality after localization.  Conversely,
fix \(i\).  If the image of a putative
\(g\in\mathfrak f_{K^{(r)}}\) in

\[
 K^{(r)}/\mathfrak m_iK^{(r)}
 \simeq\bigotimes_{j\ne i}K_j
\]

were nonzero, specialize every other factor at a \(k\)-point where that
image is nonzero.  The resulting quotient would be a localization of
\(K_i\) at an element outside \(\mathfrak m_i\), contradicting (13).
Thus \(g\) lies in the extension of every \(\mathfrak m_i\).  The vector
space splittings \(K_i=k\oplus\mathfrak m_i\) identify the intersection of
those \(r\) extended ideals with their tensor product, proving (20).

The union of the \(r\) one-axis systems

\[
 \{t_i^2,t_i^3\}
 \cup\{t_i^2P_i^n,t_i^3P_i^n:n\ge1\}
                                                               \tag{21}
\]

is an infinite monomial SAGBI basis of \(K^{(r)}\).  Its
finite-generation ideal has \(2^r\) return ladders

\[
 \prod_{i=1}^r t_i^{\epsilon_i}P_i^{n_i},
 \qquad
 \epsilon_i\in\{2,3\},\quad n_i\ge0.                \tag{22}
\]

Modulo \((t_1^4,\ldots,t_r^4)\), products of two positive-\(P_i\)-degree
classes vanish in the \(i\)-th direction.  Any finite generating set
therefore occupies a bounded \(r\)-dimensional degree box, while (22)
escapes it.  The localization/specialization proof of (20), rather than
this bounded replay, computes the exact finite-generation ideal.

This theorem is the control for the three-index option below.  It applies
when the boundary pinches and their generic invariant coordinates separate
as the tensor product (19).  It does **not** apply automatically to the
symmetric Cox coupling \(xyz=p\): there the three axes share one conductor
parameter.  That coupled construction is instead closed by the
[three-boundary Cox-fill obstruction](CONDUCTOR_THREE_BOUNDARY_COX_FILL_OBSTRUCTION.md).

## 7. Consequence for the Keller search

The comparison now separates the observed behaviors cleanly.

| model | boundary geometry | invariant outcome | finite-generation ideal |
|---|---|---|---|
| normalized \((2,3)\) slice | one effective leading divisor | finite | unit ideal |
| normalized \((2,4)\) slice | one effective leading divisor | finite | unit ideal |
| tangent leading pair \(a=0,p=0\) | disjoint by (1) | no bifiltration | not applicable |
| two-cusp control | intersecting conductor axes | non-finite | exact conductor (16) |

The next Keller-attached candidate must therefore contain two intersecting
normalization/conductor primes.  The current direct two-boundary suspension
with one reconstruction variable is unsuitable because it introduces a
third divisor.  The concrete construction choices left by the existing
boundary audit are:

1. use two reconstruction variables;
2. accept the third divisor and run a three-index saturation ledger,
   using (20) as the exact independent-boundary control; or
3. extract an LND action on a completed two-boundary normalization chart
   already present in the plane boundary compiler.

In each case the decisive test is now explicit: reduce modulo the squares of
the two conductor ideals and look for invariant bidegrees escaping every
finite rectangle.  In parallel, compute the conductor to the normalized
generic quotient and test both containments in the finite-generation ideal:
localization gives the easy inclusion, while one-axis specializations give
the exclusion outside the conductor.  If all one-axis specializations
separate, (20) computes the answer.  If they do not, retain the mixed Rees
relations: their failure to tensor-separate is precisely the extra data
that distinguishes a coupled Keller candidate from the control.

Run

```bash
.venv/bin/python scripts/verify_multiboundary_hilbert14_antichain.py
```

to check the commuting LNDs, the cusp-semigroup expansions, the
modulo-\((s^4,t^4)\) product vanishings, the finite rectangle replay, the
monomial conductor characterization, and the factorization leading-divisor
obstruction.
