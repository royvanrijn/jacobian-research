# Seven support targets suffice for the global Rodrigues gate

## Statement

In the normalized complete cyclic-frame model with \(D=1/8\), let
\[
\mathcal Z=\{p,s_0,s_1,U_{000},U_{002},U_{022},U_{222}+128s_2\}.
\]
The first six moments, together with the seven equations \(z=0\) for \(z\in\mathcal Z\), have exactly the three normalized Rodrigues points as their complex support. Adding moments through nine makes this intersection scheme-theoretically equal to those three reduced points. In particular, if \(P\) is their full defining ideal in the cyclic coordinate ring, then
\[
\boxed{I_9+(\mathcal Z)=P.}
\]

Consequently, given the retained exact local lengths seventeen at the three Rodrigues points,
\[
\dim_{\mathbb C}A_{14}=51
\quad\Longleftrightarrow\quad
z^7\in I_{14}\text{ for every }z\in\mathcal Z.
\]
Any separately verified smaller local nilpotence exponent for a target may replace seven. No one of these global power memberships is asserted here.

## 1. Reduction to a full four-parameter family

Set the seven targets to zero. The Gram relation gives
\[
s_2^3=1/64,
\]
so \(s_2\ne0\). The three trace equations reduce to
\[
U_{011}+2U_{002}=0,\qquad U_{111}+2U_{012}=0,
\qquad U_{112}+2U_{022}=0.
\]
Thus \(U_{011}=U_{112}=0\) and \(U_{111}=-2U_{012}\).

Scaling the phase by \(\lambda\), with \(\lambda^6=1\), preserves \(D=1/8\) and the seven target equations. It acts transitively on the three choices of \(s_2\). Indeed \(s_2\) scales by \(\lambda^4\), while \(U_{ijk}\) scales by \(\lambda^{-2-i-j-k}\); in particular \(U_{222}\) and \(s_2\) have the same weight modulo six. It suffices to take \(s_2=1/4\).

The cyclic frame \((v,Bv,B^2v)\) then has Gram matrix
\[
\begin{pmatrix}0&0&1/4\\0&1/4&0\\1/4&0&0\end{pmatrix}
\]
and oriented determinant \(1/8\). There is a determinant-one isometry sending it to
\[
(e_z,\tfrac12e_t,\tfrac14e_w).
\]
Cayley-Hamilton, with \(p=0\), determines the quadratic operator in these coordinates:
\[
B=\begin{pmatrix}0&0&4r\\1/2&0&0\\0&1/2&0\end{pmatrix}.
\]
This argument includes all values of \(r\), including the nilpotent case. There is no eigenvalue or nonisotropic-vector assumption.

Write
\[
R=4r,\qquad A=\tfrac32U_{001},\qquad
B_*=-\tfrac14U_{012},\qquad C=\tfrac3{32}U_{122}.
\]
The coefficient change has determinant \(-9/64\), a nonzero constant. The full remaining phase is
\[
F=w+zt-\tfrac12z^3+Rw^2+Atw^2+B_*(t^3-3ztw)+Cz^2t.
\]
For the rest of this note, write \(B\) for the scalar \(B_*\); it is no longer the quadratic operator.

The replay verifies this adapter using every symmetric cubic tensor slot, the three trace identities, the Gram matrix and the oriented determinant.

## 2. Exact first moments

Apply the sphere rule
\[
h(z^i t^j w^k)=0\quad\text{unless }i=k\text{ and }j\text{ is even},
\quad
h(z^i t^{2d}w^i)=\frac{i!(2d-1)!!}{(2i+2d+1)!!}.
\]
After multiplication by nonzero rational constants, the first two relevant equations are
\[
M_2=4AC+15B^2=0,
\qquad M_3=-4BCR-9B-6R=0.
\]
The exact primitive integer polynomials \(M_4,M_5,M_6\) are exported in the JSON certificate. Two evaluations of \(M_5\) are
\[
\left.M_5\right|_{B=C=R=0}=-78A,\qquad
\left.M_5\right|_{A=B=R=0}=286C,
\]
with the substitutions indicated explicitly.

If \(B=0\), then \(R=0\) and \(AC=0\); the displayed fifth-moment evaluations force \(A=C=0\). If \(C=0\), the second moment first forces \(B=0\). Therefore every possible nonzero solution has \(BC\ne0\).

## 3. A univariate contradiction on the complementary chart

Put \(x=BC\) and \(T=B^4\). Then \(x\ne0\), and the third-moment equation shows \(2x+3\ne0\). The first two equations give
\[
A=-\frac{15B^3}{4x},\quad
R=-\frac{9B}{2(2x+3)},\quad C=\frac{x}{B}.
\]
After substitution and multiplication by a unit, the fourth moment becomes
\[
3h(x)T+26x^2g(x)=0,
\]
where
\[
h=128x^3+1684x^2+6528x+6435,\qquad
g=176x^2+843x+990.
\]
These polynomials are coprime. Thus \(h(x)=0\) would force \(g(x)=0\), an impossibility. Consequently
\[
T=-\frac{26x^2g(x)}{3h(x)}.
\]

Substituting this into the fifth and sixth moments, and removing powers of the units \(x,B,2x+3,h(x)\), gives
\[
P_8(x)=Q_{10}(x)=0,
\]
with
\[
\begin{aligned}
P_8={}&139722752x^8+2124898304x^7+16357588736x^6\\
&+78221171808x^5+234755107902x^4+422954870481x^3\\
&+416469786462x^2+178673103180x+11925856800,
\end{aligned}
\]
and
\[
\begin{aligned}
Q_{10}={}&1339467759616x^{10}+31049215901696x^9+274113475043328x^8\\
&+1153787427753216x^7+2148687162230400x^6+401397850161216x^5\\
&-1874858578107984x^4+11224150637836302x^3+40751130406017735x^2\\
&+47989543125273750x+20140211802818250.
\end{aligned}
\]
The replay exports and verifies rational Bezout identities for both \((h,g)\) and \((P_8,Q_{10})\). Their gcds are one.

A smaller independently checkable certificate works modulo 101 without any degree drop. For \((h,g)\), the Bezout multipliers are
\[
-48x-46,\qquad -11x^2+10x+1.
\]
For \((P_8,Q_{10})\), they are
\[
\begin{aligned}
u={}&-34x^9-24x^8+23x^7+17x^6+23x^5-43x^4+12x^3+6x^2-2x+6,\\
v={}&-39x^7-22x^6-13x^5-3x^4-30x^3+26x^2-40x+10.
\end{aligned}
\]
In each case the corresponding linear combination is one modulo 101. Since the reductions preserve both polynomial degrees, this also certifies coprimality over \(\mathbb Q\). Therefore there is no nonzero solution on the complementary chart.

The origin \(R=A=B=C=0\) is the Rodrigues phase and satisfies all moments. Undoing the three possible scalar normalizations proves the stated support theorem.

## 4. The first-nine intersection is reduced

In the parameter order \((R,A,B,C)\), the linear parts of the primitive odd-moment equations of indices three, five, seven and nine form the matrix
\[
\begin{pmatrix}
-6&0&-9&0\\
0&-78&0&286\\
646&0&6460&0\\
0&5244&0&0
\end{pmatrix}.
\]
Its determinant is \(49411883664\ne0\). Thus the first-nine intersection is reduced at the origin. The first-six support result shows there are no other points, so its entire affine coordinate ring is the residue field.

This argument transports to all three values of \(s_2\); the roots of \(s_2^3-1/64\) are simple. Therefore \(I_9+(\mathcal Z)=P\) scheme-theoretically in the full normalized cyclic model. The replay separately verifies the four linear derivatives directly from \(n h(F_0^{n-1}Q)\), avoiding any reliance on a truncated high-order multiplication table.

## 5. Consequence for the global power-membership gate

If the seven powers \(z^7\) belong to \(I_{14}\), every complex closed point of \(A_{14}\) has all seven target coordinates zero. It therefore lies at one of the three Rodrigues points by the support theorem, which used only the first six moments. A finite-type complex algebra with this finite support is Artinian. Its local lengths are the retained values seventeen, so its total length is fifty-one.

Conversely, global length fifty-one leaves exactly the three known Rodrigues local algebras. Their maximal ideals have seventh power zero; hence all seven target seventh-powers vanish globally.

This reduction neither proves any global power membership nor assumes that the full global quotient is finite. It turns seven such memberships into a sufficient and necessary certificate.

## Replay

Run `python verify_odd_slice.py`. The script uses exact rational arithmetic, verifies the complete cyclic adapter, computes moments by multinomial enumeration, independently recomputes them by direct spatial expansion, and verifies both rational and mod-101 Bezout identities. The exported data are `odd_slice_support_data.json`; the run log is `verification.log`.
