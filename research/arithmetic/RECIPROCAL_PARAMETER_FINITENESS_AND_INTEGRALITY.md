# Finiteness and denominator control for reciprocal paired pencils

For every degree \(n=4h+2\) and every aggregate reciprocal sign split,
the normalized reciprocal equations define a finite scheme. Its exact
length, including degenerate solutions and multiplicities, is

\[
             2h\binom{2h}{p}\binom{2h}{q},
             \qquad p+q=h.
\]

Every normalized good complex model is therefore algebraic. Over any
number field, its parameter coefficients and numerator coefficients are
integral outside primes at most \(n-1\); the reciprocal parameter and
constant numerator coefficient are units there.

The key is a uniform description of the equations at weighted infinity:
every putative nonzero limit produces a polynomial with only one finite
critical value. An elementary root count forces that polynomial to be a
monomial and makes the limit zero. This applies to all sign splits,
including the balanced and nonsquare-descent cases.

## The raw normalized equations

Fix \(h\ge1\), \(p,q\ge0\), \(p+q=h\), and \(L=4h+1=n-1\).
Use the normalization of the
[mixed coupled-ODE note](MIXED_RECIPROCAL_COUPLED_ODE.md).
Write

\[
 P=w^p+\sum_{i=0}^{p-1}a_iw^i,\qquad
 Q=w^q+\sum_{j=0}^{q-1}b_jw^j,\qquad s=P Q.
\]

There are \(h+1\) coefficient variables: all \(a_i,b_j\), and \(r\).
Set

\[
 \Sigma(v)=v^h s(v+v^{-1}),\qquad
 F(v)=L\sum_{\ell=0}^{2h+1}
       \frac{\Sigma_{\ell-1}-r\Sigma_\ell}{2\ell-1}v^\ell,
                                                        \tag{1}
\]

where the out-of-range coefficients of \(\Sigma\) are zero.
In particular \(F\) is monic, \(F(0)=Lr\), and
\(2vF'-F=L(v-r)\Sigma\).

Let

\[
 X=\frac{F(v)+vF(1/v)}{v+1},\qquad
 Y=\frac{F(v)-vF(1/v)}{v-1},\qquad
 G=\frac{F(r)}{2r}.                                      \tag{2}
\]

Here \(X,Y\) are polynomials in \(w=v+v^{-1}\), monic of degree
\(2h\). The apparent division by \(r\) in \(G\) cancels identically,
so \(G\) is a polynomial in the coefficient variables. Moreover

\[
                          G|_{r=0}=L.                   \tag{3}
\]

Define the **raw scheme** by

\[
       \operatorname{rem}_P X=0,\qquad
       \operatorname{rem}_Q Y=0,\qquad G=0.               \tag{4}
\]

There are exactly \(p+q+1=h+1\) scalar equations. A degree-zero
factor contributes no remainder equations. These equations include
degenerate solutions; no automatic \(w=\pm2\) factors or
discriminant factors have been divided out. They are necessary for
every normalized good reciprocal model. On the good open set of
the mixed note they also suffice.

All coefficients of (4) lie in

\[
       A=\mathbb Z[1/\ell:\ell\text{ prime},\ \ell\le L].
\]

The constructions use only integer operations, monic division,
and division by 2 or by odd integers at most \(L\).

## Theorem 1. Weighted infinity consists only of the origin

Assign weights

\[
       \operatorname{wt}(a_i)=p-i,\qquad
       \operatorname{wt}(b_j)=q-j,\qquad
       \operatorname{wt}(r)=1.
\]

Over every algebraically closed field of characteristic zero or
characteristic greater than \(L\), the common zero of the
weighted-leading forms of (4) is only the origin.

To prove this, write \(s(w)=\sum_{k=0}^h s_kw^k\), \(s_h=1\),
and give \(w\) weight one. The common top part of \(X,Y\) is

\[
 Z(w)=L\sum_{k=0}^h s_k
       \left(\frac{w^{h+k}}{2h+2k+1}
           -\frac{r\,w^{h+k-1}}{2h+2k-1}\right).          \tag{5}
\]

It is monic of degree \(2h\), and its coefficients have weights
\(2h-j\) in degree \(w^j\). This follows either directly from (1),
or from the coupled operators: at top weight \(C_h,D_h\) both
become \(w^h\), and the two operators both become
\(T\mapsto T+2wT'\). Thus

\[
          Z+2wZ'=L w^{h-1}(w-r)s(w).                    \tag{6}
\]

Monic remainder division respects these weights. The top equations
in (4) are consequently

\[
                    P\mid Z,\qquad Q\mid Z,\qquad G_{\rm top}=0,
\]

where

\[
 G_{\rm top}
  =-L\sum_{k=0}^h
       \frac{s_k r^{h+k}}{(2h+2k-1)(2h+2k+1)}
  =\frac{Z(r)}2.                                        \tag{7}
\]

Now form the monic degree-\(L\) polynomial \(g(t)=tZ(t^2)\).
Equations (6)--(7) give

\[
                 g'(t)=L t^{2h-2}(t^2-r)s(t^2).         \tag{8}
\]

Every zero of this derivative is a zero of \(g\):
\(g(0)=0\); the roots above \(r\) have \(Z(r)=0\); and every root
of \(s=P Q\) is a root of \(Z\). No coprimality or squarefreeness
of \(P,Q\) is needed for this statement.

A polynomial of degree \(L\) whose derivative has all its roots
among the polynomial's roots has only one distinct root when
the characteristic is zero or greater than \(L\). Indeed, if
its distinct root multiplicities are \(m_1,\ldots,m_d\), then
the derivative has multiplicities \(m_i-1\) there. The equality
of its full degree \(L-1\) with \(\sum_i(m_i-1)=L-d\) forces
\(d=1\). The characteristic restriction ensures this derivative
multiplicity formula and its leading degree.

Since \(g\) is monic and \(g(0)=0\), it follows that \(g=t^L\),
hence \(Z=w^{2h}\). Equation (6) then reduces to

\[
                            (w-r)s(w)=w^{h+1}.
\]

Thus \(r=0\), \(P=w^p\), \(Q=w^q\). Every coefficient variable
is zero, proving the theorem.

## Theorem 2. Exact finite raw length

Over a field of characteristic zero or characteristic greater than
\(L\), the raw scheme (4) is a zero-dimensional complete intersection
of length

\[
                  D_{h,p,q}=2h\binom{2h}{p}\binom{2h}{q}.       \tag{9}
\]

These are scheme lengths, not counts of distinct good covers.

Here is the filtered-algebra justification. The polynomial ring has
\(h+1\) variables with weights
\(1,\ldots,p;\ 1,\ldots,q;\ 1\). Theorem 1 shows that its
\(h+1\) homogeneous leading equations cut out only the origin.
They therefore have height \(h+1\) and form a regular sequence
in the polynomial ring. This is the standard dimension criterion
in a Cohen--Macaulay polynomial ring; see
[Stacks, Lemma 10.104.2](https://stacks.math.columbia.edu/tag/02JN)
and [Section 10.104](https://stacks.math.columbia.edu/tag/00N7).

For clarity, the passage to the inhomogeneous equations does not
require an unstated flatness assumption. For a regular sequence
of leading forms, their relations are generated by the elementary
Koszul relations \(H_j e_i-H_i e_j\). This follows inductively:
in a relation \(\sum c_iH_i=0\), reduce modulo the first \(m-1\)
forms; regularity puts \(c_m\) in their ideal, and induction treats
the remaining relation. If the highest-weight terms cancel in
a representation \(\sum d_i f_i\) of a raw ideal element, use
these relations, with \(f_i\) in place of \(H_i\), to lower the
maximum weight of the representation without changing its sum.
Iterating proves that the initial ideal of the raw equations is
exactly their homogeneous leading ideal.

The associated graded quotient therefore has the Hilbert series
of that homogeneous complete intersection. Each homogeneous
nonzerodivisor of weight \(d\) multiplies the preceding Hilbert
series by \(1-T^d\), by its multiplication exact sequence. The
equation weights are

\[
 (2h,2h-1,\ldots,2h-p+1),\quad
 (2h,2h-1,\ldots,2h-q+1),\quad 2h.
\]

Taking the limit at \(T=1\) divides their product by \(p!q!\),
which is (9). This also proves finiteness of the ungraded quotient.

For degree 26 the raw lengths are:

| Pair counts \((p,q)\) | Raw sign degrees | Raw scheme length |
|---|---|---:|
| \((1,5)\) | \((2,10)\) | 114048 |
| \((2,4)\) | \((4,8)\) | 392040 |
| \((3,3)\) | \((6,6)\) | 580800 |

Removing automatic factors, restricting to the good locus,
or eliminating variables can reduce the lengths. For example,
\(h=2,p=q=1\) has raw length 64, while dividing the automatic
\(c+2,d-2\) factors gives the degree-ten mixed system of length
36. The four-variable degree-26 \((2,10)\) system is similarly
a smaller necessary locus; its
[global rational obstruction](DEGREE_TWENTY_SIX_MIXED_2_10_GLOBAL_OBSTRUCTION.md)
does not assert that the raw complex scheme is empty.

Every complex point of (4) is algebraic over \(\mathbb Q\):
the raw coordinate algebra over \(\mathbb Q\) is finite-dimensional.
Thus there are only finitely many normalized good reciprocal
models in each sign split. This includes models whose initial
coefficient field was not assumed algebraic.

## Theorem 3. Uniform integrality and unit constraints

Let a solution of (4) lie in any finite extension \(K/\mathbb Q_\ell\),
with \(\ell>L\). Then all \(a_i,b_j,r\) are integral, and \(r\)
is a unit. This includes residue extensions and ramification
of arbitrary degree.

Theorem 1 holds over the algebraic closure of the residue field.
Apply the weighted valuation argument from the
[global obstruction note](DEGREE_TWENTY_SIX_MIXED_2_10_GLOBAL_OBSTRUCTION.md):
a nonintegral coordinate tuple can be rescaled in a finite totally
ramified extension to a nonzero common residue zero of the top
forms. Their only such zero is the origin, so the tuple is integral.
If \(r\) were not a unit, reducing \(G\) and using (3) would give
\(0=L\) in residue characteristic \(\ell>L\), a contradiction.

Consequently, for a normalized good model over a number field,
write \(S\) for the prime ideals above rational primes at most \(L\).
Then

\[
 \begin{gathered}
 a_i,b_j\ \text{are }S\text{-integers},\qquad
 r,\rho=1/r,\lambda=\rho^2\ \text{are }S\text{-units},\\
 N(u)=\rho^{2h+1}F(u/\rho)
 \ \text{has }S\text{-integral coefficients},\qquad
 N(0)=L\rho^{2h}=L\lambda^h\ \text{is an }S\text{-unit}.
 \end{gathered}
\]

If the sign factors require a quadratic coefficient extension,
apply this statement there and descend the conclusions for
coefficients already in the original field.

This does not prove good reduction: critical points or values
can still collide modulo \(\ell\). A finite-field exclusion must
include all residue solutions of an appropriate necessary system,
or separately exclude every lift of its degenerate residue points.

## Replay, sources, and scope

~~~text
.venv/bin/python research/arithmetic/verify_reciprocal_parameter_finiteness.py --write-receipt research/arithmetic/reciprocal_parameter_finiteness_verification.json
~~~

The checker verifies actual coupled-operator coefficients, the
weighted \(Z\) limits, \(Z(r)=2G_{\rm top}\), equation (8), denominator
support, and equation (3) in retained symbolic cases. It independently
computes full raw Groebner quotient lengths in degrees six and ten,
including the mixed length 64. The written arguments supply the
uniform quantifiers and exact-length formula.

The primary commutative-algebra inputs are retained as
[Stacks tag 02JN](sources/stacks_02JN_regular_sequence_2026_11_18.html)
and [Stacks Section 00N7](sources/stacks_00N7_cohen_macaulay_2026_11_18.html).
The filtered-ideal and Hilbert-series steps are explained above.
No additional cover classification is imported.

This is a structural theorem about the normalized reciprocal ansatz.
It does not construct rational degree-26 models, classify the remaining
rational points, or give finiteness for nonreciprocal normalized
paired-critical models. The raw lengths count many degenerate points
and must not be reported as counts of Hurwitz orbits or good covers.
No novelty comparison or external human review is asserted.
