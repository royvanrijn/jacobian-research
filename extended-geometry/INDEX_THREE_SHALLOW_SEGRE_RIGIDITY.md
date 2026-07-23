# Index-three rigidity for the one-stage shallow Segre circuit

## Result

Let \(k\) be a characteristic-zero field, let

\[
Q:k^m\longrightarrow k^m,\qquad
C:k^m\longrightarrow k^m
\]

be homogeneous of degrees two and three, and put

\[
\mathcal F(x,y,t)=(x,y,t)+\mathcal H(x,y,t),
\qquad
\mathcal H=(tQ(x)-t^2y,\ C(x),\ 0).
\tag{1}
\]

This is the one-stage full-output Segre circuit used to turn
\(K=x+Q+C\) into a cubic-homogeneous map.  It has \(2m+1\) variables and
only one bilinear homogenization layer.

> **Theorem (shallow Segre rigidity).**  
> If \((J\mathcal H)^3=0\), then both \(K=x+Q+C\) and \(\mathcal F\) are
> polynomial automorphisms.  In particular, no noninjective
> cubic-homogeneous index-three map occurs in the class (1), in any
> dimension.

This is an exhaustive classification of a well-defined shallow-circuit
class, not a failed search.  It covers the odd dimensions \(7,9,11\) in the
proposed first test range by taking \(m=3,4,5\).

The exact checker is
[`verify_index_three_shallow_segre_rigidity.py`](../scripts/verify_index_three_shallow_segre_rigidity.py).
It verifies the polynomial block identities, the integrability/commutator
identity, and a seven-variable rank-four exact-index-three calibration.

## 1. Collision gauge

There is a small but important correction to the proposed quotienting.
Translations do not preserve the form \(X+H\) with \(H\) cubic homogeneous.
They should not be used to replace a collision by \((0,e_1)\).

No translation is needed.  Suppose

\[
F=X+H,\qquad H(\lambda x)=\lambda^3H(x),
\qquad (JH)^3=0.
\]

Two distinct collision points cannot be linearly dependent.  Indeed, if
\(b=\lambda a\), then \(F(a)=F(b)\) gives

\[
(1-\lambda)\bigl(a+(1+\lambda+\lambda^2)H(a)\bigr)=0.
\]

The cases \(a=0\) and \(1+\lambda+\lambda^2=0\) are immediate
contradictions.  Otherwise Euler's identity gives

\[
JH(a)a=3H(a)
=-\frac{3}{1+\lambda+\lambda^2}a,
\]

contradicting nilpotence of \(JH(a)\).  Thus \(a,b\) are independent, and a
linear conjugacy sends them to \(e_1,e_2\).  This is the correct homogeneous
collision gauge.

On their span, write

\[
JH(se_1+te_2)=s^2A+stB+t^2C.
\]

The necessary pencil problem is

\[
(s^2A+stB+t^2C)^3=0,
\tag{2}
\]

together with

\[
Be_1=2Ae_2,\qquad Be_2=2Ce_1,\qquad
Ae_1-Ce_2=3(e_2-e_1).
\tag{3}
\]

Equations (2)--(3) are the recommended finite-field discovery frontend for
the unrestricted tensor search.  Any survivor must still be extended to a
full cubic tensor and checked away from the collision plane.

## 2. Block calculation

Put

\[
A=JQ(x),\qquad B=JC(x).
\]

In the variable order \((x,y,t)\),

\[
N:=J\mathcal H=
\begin{pmatrix}
tA&-t^2I&Q-2ty\\
B&0&0\\
0&0&0
\end{pmatrix}.
\tag{4}
\]

The \((x,y)\)-block of \(N^3\) is

\[
\begin{pmatrix}
t^3(A^3-AB-BA)&-t^4(A^2-B)\\
t^2(BA^2-B^2)&-t^3BA
\end{pmatrix},
\tag{5}
\]

and its last column is

\[
\binom{t^2(A^2-B)(Q-2ty)}{tBA(Q-2ty)}.
\tag{6}
\]

Consequently

\[
N^3=0
\quad\Longleftrightarrow\quad
B=A^2\ \text{ and }\ A^3=0.
\tag{7}
\]

The reverse implication in (7) follows immediately from (5)--(6).

## 3. Integrability forces a constant triangular flag

Write

\[
A(x)=\sum_{i=1}^m x_iA_i.
\]

Since \(A=JQ\), equality of mixed second derivatives gives

\[
A_i e_j=A_j e_i.
\tag{8}
\]

Since \(A^2=JC\), equality of mixed second derivatives of \(C\) gives

\[
\partial_i(A^2)e_j=\partial_j(A^2)e_i.
\tag{9}
\]

Using (8), the difference between the two sides of (9) is exactly

\[
(A_iA_j-A_jA_i)x.
\tag{10}
\]

Hence all coefficient matrices \(A_i\) commute.  Each \(A_i\) is nilpotent
because \(A_i^3=0\), obtained from \(A(x)^3=0\) by setting \(x=e_i\).
A commuting family of nilpotent matrices over an algebraic closure is
simultaneously strictly triangular.  In that constant linear gauge, both
\(JQ=A\) and \(JC=A^2\) are strictly triangular.  Therefore

\[
K=x+Q+C
\]

is a triangular polynomial automorphism.

Finally, on the fiber \(t=\tau\), eliminate \(y\) from (1):

\[
\mathcal F(x,y,\tau)=(u,v,\tau)
\quad\Longleftrightarrow\quad
x+\tau Q(x)+\tau^2C(x)=u+\tau^2v,\qquad
y=v-C(x).
\]

The first map is the scaled triangular automorphism
\(\tau^{-1}K(\tau x)\) (with the evident value at \(\tau=0\)).
Thus \(\mathcal F\) is a polynomial automorphism as well.

## 4. Exact in-range calibration

The theorem is not the assertion that index three collapses to index two.
Let

\[
Q=(0,\tfrac12x_1^2,x_1x_2),\qquad
C=(0,0,\tfrac13x_1^3).
\tag{11}
\]

Then

\[
A=
\begin{pmatrix}
0&0&0\\
x_1&0&0\\
x_2&x_1&0
\end{pmatrix},
\qquad
JC=A^2,
\qquad
A^2\ne0,\quad A^3=0.
\]

The resulting map (1) has seven variables,

\[
\operatorname{rank}J\mathcal H=4,\qquad
(J\mathcal H)^2\ne0,\qquad
(J\mathcal H)^3=0.
\]

It lies exactly in the proposed range \(6\le n\le12\),
\(3\le\operatorname{rank}JH\le6\), but it is triangularizable and
invertible.  It is therefore a useful positive control for collision-aware
search code.

## 5. Rank-compressed successor

For the repository's rectangular rank-compressed circuit, write

\[
\mathcal H=(tQ+t^2UY,\ -c,\ 0),
\]

where \(U\) is an \(n\)-by-\(r\) constant matrix,
\(A=JQ\), \(D=Jc\), \(R=UD\), and \(E=A^2-R\).
Direct multiplication gives the exact residual system

\[
\begin{aligned}
EU&=0,\\
DE&=0,\\
DAU&=0,\\
A^3-RA-AR&=0,\\
EQ&=0,\\
DAQ&=0.
\end{aligned}
\tag{12}
\]

When \(r=n\) and \(U\) is invertible, (12) reduces to the theorem above.
For \(r<n\), the only possible escape is the rectangular defect

\[
E=A^2-UD,
\]

which must simultaneously vanish on \(\operatorname{im}U\), land in
\(\ker D\), and annihilate \(Q\).  This is a sharper finite-field search
target than the full cubic coefficient space.  A survivor of (12) with an
exact collision would be the next candidate to promote; absence of
survivors in an unclassified beam is not a lower bound.
