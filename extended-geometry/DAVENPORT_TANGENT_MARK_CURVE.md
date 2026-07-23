# The Davenport tangent-mark curve

The global relative weighted Sunada construction adjoins a tangent mark
\(r\) satisfying

\[
g_T(r)-g_T(0)-r g_T'(0)=0,\qquad r\ne0.               \tag{1}
\]

This note determines that marking curve exactly.  It is rational, but its
nonzero-mark open is punctured at two rational points.  Consequently it
cannot be nontrivially parametrized by affine one-space.  This rules out the
most direct attempt to turn the relative weighted construction into an
absolute polynomial family by merely replacing \(T,r\) with polynomial
functions of one affine parameter.

Work over \(K=\mathbb Q(a)\), \(a^2+a+2=0\).

## 1. The marking equation

After removing the double factor \(r^2\), put

\[
J(T,r)=
\frac{g_T(r)-g_T(0)-r g_T'(0)}{r^2}.                  \tag{2}
\]

Direct reduction gives

\[
\begin{aligned}
7J(T,r)={}&
\bigl(14ar+28a-21r-14\bigr)T^2\\
&+7(a+1)r^2(r+1)T+r^5.                               \tag{3}
\end{aligned}
\]

Thus \(J\) has degree five in \(r\), as used by the finite tangent-mark base
change, but only degree two in \(T\).  Its \(T\)-discriminant is

\[
\operatorname {disc}_T(J)=r^4Q(r),                   \tag{4}
\]

where

\[
Q(r)=
\frac{(5-a)r^2+(-2a-6)r+7a-7}{7}.                    \tag{5}
\]

The quadratic \(Q\) is separable:

\[
\operatorname {disc}_r(Q)=-\frac{16(11a-7)}{49}\ne0. \tag{6}
\]

Away from \(r=0\), completing the quadratic in \(T\) identifies the
normalization of (3) birationally with the smooth conic

\[
\boxed{v^2=Q(r).}                                    \tag{7}
\]

More precisely, if

\[
A_2=(2a-3)r+4a-2,
\]

then one branch of the inverse birational map is

\[
T=\frac{r^2(-(a+1)(r+1)+v)}{2A_2}.                   \tag{8}
\]

Changing \(v\) to \(-v\) gives the other quadratic branch.

## 2. Rational parametrization

At \(r=0\),

\[
Q(0)=a-1=(1+a)^2.                                    \tag{9}
\]

Hence the conic has the two \(K\)-points

\[
P_\pm=(0,\pm(1+a)).                                   \tag{10}
\]

Lines \(v=1+a+mr\) through \(P_+\) give the explicit parametrization

\[
r(m)=
-\frac{2(7am+a+7m+3)}{a+7m^2-5},\qquad
v(m)=1+a+mr(m).                                      \tag{11}
\]

Substitution in (5) verifies \(v(m)^2=Q(r(m))\).
Combining (8) and (11) gives rational functions \(T(m),r(m)\) which
globalize the weighted family over a punctured rational base.  This is a
useful simplification: the former degree-five coefficient extension is a
rational curve, not a positive-genus obstruction.

## 3. Why rationality is not enough

The weighted normalization requires \(r\ne0\).  On the smooth normalization
(7), this removes both points \(P_+\) and \(P_-\).  Since the conic has a
\(K\)-point, its projective completion is \(\mathbb P^1_K\).  The
nonzero-mark locus is therefore contained in

\[
\mathbb P^1_K\setminus\{P_+,P_-\}\simeq\mathbb G_{m,K}. \tag{12}
\]

Any morphism

\[
\mathbb A^1_K\longrightarrow\mathbb G_{m,K}
\]

is constant because the image of the Laurent coordinate must be a unit in
\(K[t]\), hence an element of \(K^*\).  It follows that:

> There is no nonconstant affine-line reparametrization of the present
> nonzero Davenport tangent-mark family.

Equivalently, the denominators in (11) cannot all be removed by changing the
single affine parameter.  The relative weighted maps can be written over a
rational boundary complement, but not over an affine-line base by this
marking construction.

This is independent of the Cox-boundary unit obstruction.  The
[Cox audit](DAVENPORT_COX_BOUNDARY_OBSTRUCTION.md) rules out stable
straightening of the direct common-target construction; the result here
rules out polynomial affine-line descent of the current relative weighted
construction.

## 4. Moving both marks: the affine ansatz is empty

One possible escape is to abandon the fixed tangent point \(Y=0\).  Put

\[
\mathcal T(T,q,r)=
\frac{g_T(r)-g_T(q)-g_T'(q)(r-q)}{(r-q)^2}.           \tag{13}
\]

This is polynomial because the numerator has a double zero along \(r=q\).
Test the full affine-linear moving-mark ansatz

\[
q(T)=q_0+q_1T,\qquad r(T)=r_0+r_1T.                   \tag{14}
\]

Substitution in (13) gives a polynomial of degree five in \(T\).  Take its
six coefficient equations in

\[
K[q_0,q_1,r_0,r_1].
\]

Exact Gröbner reduction gives the unit ideal:

\[
\boxed{
\left([T^0]\mathcal T,\ldots,[T^5]\mathcal T\right)=(1).
}                                                     \tag{15}
\]

Therefore no tangent pair has both marks affine-linear in \(T\).  This
includes constant marks and proportional linear marks.  A moving-mark
escape, if one exists, must use higher polynomial degree, rational functions
with poles, more base variables, or a different normalization mechanism.

Such an escape does exist after writing \(T\) as a quadratic function of a
new parameter: the
[proportional-section audit](DAVENPORT_PROPORTIONAL_TANGENT_SECTIONS.md)
classifies four sections with marks linear in \(s\) and \(T=\tau s^2\).
The statement above remains the exact obstruction for affine-linear
functions of \(T\); it is not an obstruction to this ramified base change.

## 5. Remaining scope

The theorem concerns the normalization based at the explicit tangent point
\(Y=0\).  It does not rule out:

1. a higher-degree or multivariable family of tangent pairs \((q,r)\);
2. an affine modification using more than one base variable;
3. a non-weighted suspension of the \((1,1,2)\) Cox boundary ledger; or
4. an absolute construction whose inverse core is a different
   \(GL_3(\mathbb F_2)\) realization.

It does remove the simplest hoped-for descent: rationality of the
degree-five tangent-mark extension does not by itself produce an absolute
polynomial parameter family.

## 6. Verification

Run

```bash
.venv/bin/python scripts/verify_davenport_tangent_mark_curve.py
```

The checker verifies (3)--(11), smoothness of the conic, the two distinct
rational punctures, and the unit Gröbner basis (15).
