# Geometric-degree and generic-fiber sieve for F2 `(75,125)`

> **Status.** Exact consequence of two published geometric-degree theorems
> and the certified F2 carrier vertex, followed by conditional generic-slice
> identities.  A hypothetical F2 counterexample has geometric degree
> `d_geo<=28`.  Nguyen Van Chau's semigroup theorem and the known exclusion
> through degree five sharpen the possibilities to
> `d_geo in {6,8,9,...,28}`.  The squarefree terminal row therefore has 22
> possible degrees, while the doubled terminal row has the 17 possibilities
> `12,...,28`.  This is a finite reduction, not an exclusion.
>
> For affine nonproperness components with normalization degrees `(3k_j,5k_j)`,
> the generic coordinate-fiber formulas below are exact.  They deliberately
> retain the unknown pole counts `p,q`; the terminal Belyi passport does not
> imply `p=q=2`.

The arithmetic is replayed by
[`verify_f2_geometric_degree_sieve.py`](../scripts/verify_f2_geometric_degree_sieve.py).

## 1. Strict geometric-degree ceiling

The F2 standard pair has leading carrier vertices

\[
 v_P=(15,60),\qquad v_Q=(25,100),
\]

and coprime leading ratio `(a_0,b_0)=(3,5)`.  These are the unramified
coordinates corresponding to the transformed vertices `(75,60)` and
`(125,100)` after adjoining `X=x^(1/5)`.

Makar-Limanov proves that, for a normalized Jacobian pair with leading
vertex `(m,n)` and coprime vertical-degree ratio `a_0/b_0`, the field degree
`N=[C(x,y):C(P,Q)]` satisfies the strict inequality

\[
 N<(n-m)\frac{n b_0}{n(a_0+b_0)-a_0}.
\]

See page 10 of
[Makar-Limanov, *On the Newton polytope of a Jacobian pair*](https://arxiv.org/abs/2106.06869).
For F2 this is

\[
 N<45\frac{60\cdot5}{60\cdot8-3}
   =\frac{1500}{53}<29.
\]

In characteristic zero the field degree is the geometric degree, hence

\[
 \boxed{d_{\rm geo}\le28.}
\]

This application uses the certified F2 standard-pair normalization.  It is
not an assertion that an arbitrary degree pair `(75,125)` with no chosen
Newton corner automatically has the same vertex.

## 2. Exact candidate degrees

Order the coordinate degrees as `125=25*5` and `75=25*3`.  Theorem B of
[Nguyen Van Chau, *Non-zero constant Jacobian polynomial maps of C^2*](https://matwbn.icm.edu.pl/ksiazki/apm/apm71/apm7135.pdf)
gives

\[
 d_{\rm geo}=5r+3s,
 \qquad r,s\ge0,\quad r+s\ge1.
\]

Geometric degrees at most five are excluded.  Intersecting this fact, the
semigroup, and the strict ceiling gives

\[
 \boxed{d_{\rm geo}\in\{6,8,9,10,\ldots,28\}.}
\]

In particular `7` is impossible.  The terminal cover gives the independent
lower bounds

\[
\begin{array}{c|c}
\text{terminal stratum}&\text{remaining degrees}\\ \hline
\text{squarefree}&6,8,9,\ldots,28,\\
\text{double}&12,13,\ldots,28.
\end{array}
\]

The former Bezout interval `6..9375` has therefore collapsed to 22 integer
degrees; the doubled interval has collapsed to 17.

## 3. Safe generic-slice identities

Let the affine nonproperness components be `C_j`, with polynomial
normalizations of degrees `(3k_j,5k_j)`.  Over the generic point of `C_j`,
write the boundary rows as `(e_ji,f_ji)` and put

\[
 A_j=\sum_i e_{ji}f_{ji},\qquad
 F_j=\sum_i f_{ji}.
\]

Thus `A_j` is the moved-sheet degree and `F_j` is the number of boundary
punctures counted with residue degree.  A generic vertical coordinate line
meets `C_j` in `3k_j` points and a generic horizontal line meets it in
`5k_j` points.  Euler integration of the finite normalization cover gives

\[
 \boxed{\chi(P^{-1}(c))=d-3\sum_jk_jA_j,}
 \qquad
 \boxed{\chi(Q^{-1}(c))=d-5\sum_jk_jA_j.}       \tag{3.1}
\]

Let `p` and `q` be the still-unresolved pole-branch counts over target
infinity, away from the affine nonproperness punctures.  Then

\[
 \boxed{n_P=p+3\sum_jk_jF_j,}
 \qquad
 \boxed{n_Q=q+5\sum_jk_jF_j.}                  \tag{3.2}
\]

Set

\[
 T=\sum_jk_j(A_j-F_j)
   =\sum_jk_j\sum_i(e_{ji}-1)f_{ji}.           \tag{3.3}
\]

Using `chi=2-2g-n`, equations (3.1)--(3.2) become

\[
 \boxed{2g_P=2-d-p+3T,\qquad
        2g_Q=2-d-q+5T.}                        \tag{3.4}
\]

They recover Chau's Theorem C identically:

\[
 75(\chi_Q-d)=125(\chi_P-d).
\]

Consequently Theorem C adds no new equation after the complete affine
branch ledger has been entered.  The useful new constraints are instead
integrality, nonnegative genus, and puncture bounds.  For example,

\[
 T\ge\left\lceil\frac{d+p-2}{3}\right\rceil
 \ge\left\lceil\frac{d-1}{3}\right\rceil.      \tag{3.5}
\]

Chau's puncture theorems also imply, for a nonautomorphism,

\[
 n_P\ge4,\qquad n_Q\ge6,                       \tag{3.6}
\]

where `Q` is the degree-125 coordinate.  Equations (3.2), (3.4), and (3.6)
are the correct finite genus/puncture sieve.  Replacing `p,q` by `2,2`
without a separate boundary proof would be invalid.

## 4. Immediate computational target

The degree ceiling changes the global search qualitatively.  It is now
finite to:

1. enumerate transitive complement actions only through index 28;
2. require a fixed affine sheet for every nonproperness meridian;
3. require the terminal `A_6` section in the global image (and use
   nonsolvability plus divisibility of the image order by `360` as weaker
   necessary filters when an exact section test is unavailable);
4. apply (3.2)--(3.6) to the surviving action rows; and
5. run the localized-`ch_2` packet budget only on that finite list.

An `A_6` section need not be a composition factor of the ambient group;
`A_7\supset A_6` is the smallest immediate counterexample.  Therefore a
composition-factor filter cannot be used as a necessary-condition
exclusion here.

The accompanying low-index research driver is
[`research_f2_affine_k1_low_index_actions.py`](../scripts/research_f2_affine_k1_low_index_actions.py).
It is exploratory until its complete index-28 output is pinned and checked.
