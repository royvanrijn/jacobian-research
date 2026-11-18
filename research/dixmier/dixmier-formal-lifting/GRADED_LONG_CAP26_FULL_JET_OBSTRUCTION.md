# The complete arbitrary-parity graded cap-26 Long obstruction

There is no scalar graded canonical Long jet of common spatial degree
at most **26** satisfying all six Weyl commutators through \(h^5\),
over any characteristic-zero field. Every coefficient of every image,
including both odd and even first-image deformations, is allowed.
This is a field classification; no equality of nonreduced parameter
schemes is asserted.

Together with the actual rational even degree-27 positive point in
[the minimal-jet proof](MINIMAL_EVEN_LONG_JET_AND_FIXED_LOWER_OBSTRUCTION.md),
the minimum degree for an arbitrary-parity graded full jet through
\(h^5\) is exactly **27**. Neither assertion constructs an all-orders
finite lift or treats jets without the stated grading.

The [checker](check_graded_long_cap26_full_jet.py) uses only the Python
standard library and exact rational arithmetic. The
[witness](graded_long_cap26_full_jet_witness.json) supplies nonzero
modular minors, independent kernel minors, every coefficient of the
complete lower families, and rational obstruction functionals. The
[receipt](graded_long_cap26_full_jet_controls.json) pins all inputs.
Discovery used exact rational linear algebra; replay requires no
SymPy or FLINT.

## Coefficients and the complete first kernel

Use source order \((X,Q,P,Z)\), brackets
\(\{P,X\}=\{Z,Q\}=1\), weights \((-1,1,4,2)\), and \(h\)-weight 3.
Let \(f=(R,T,D,S)\) be the full polynomial Long symbols from the
[canonical transcription](../dixmier-audit/verify_long_pbw.py).
Their weights are \(w=(-1,2,4,1)\). Products below are products of
commuting Weyl-symbol variables, followed by Weyl quantization.

Write \(F_i=f_i+\sum_{m\ge1}h^mF_{m,i}\). The ENTIRE coefficient
basis at index \(m\) and cap 26 is
\[
 X^\alpha Q^\beta P^\gamma Z^\delta,\qquad
 \alpha=\beta+4\gamma+2\delta+3m-w_i\ge0,\qquad
 \alpha+\beta+\gamma+\delta\le26.                       \tag{1}
\]
Define
\[
 (\mathcal A_m U)_{ij}=\{U_i,f_j\}+\{f_i,U_j\},\qquad i<j.
\]
The complete dimensions and rational ranks are
\[
\begin{array}{c|rrrr|r|r|r}
m&R&T&D&S&\text{columns}&\operatorname{rank}\mathcal A_m&
 \dim\ker\mathcal A_m\\ \hline
1&114&156&188&141&599&588&11\\
2&80&114&141&102&437&431&6\\
3&54&80&102&71&307&306&1\\
4&34&54&71&47&206&-&-
\end{array}                                                   \tag{2}
\]
The square minors of sizes 588, 431, and 306 have determinants
\(476290,181171,525733\) modulo the prime \(1000003\).
The matching independent exact kernels prove the upper rank bounds.

Here are compact generators. Put
\[
\begin{split}
t&=XQ,\quad B=3X^2P+2(1-3t)Z,\quad
 A_0=Q+3XQ^2-XB/3,\quad Y=RA_0-3t,\\
H_4&=Y^4+4tY^2-\tfrac{16}{3}tY+10t^2,\\
N_5&=9Y^5+45tY^3-60tY^2+45t^2Y+110tY-720t^2,\\
H_6&=R^4D+N_5/25,\\
H_5&=9Y^5+45tY^3-60tY^2+135t^2Y-10tY-555t^2,\\
E&=X^2Z,\qquad J=(RS)^2,\qquad L_0=X^4P,\qquad N_0=X^3QZ .
\end{split}                                                   \tag{3}
\]
The first kernel is precisely the span of the eleven tuples
\[
 K_G=\{G,f\},\qquad
 G\in(t,Y,Y^2,R^2T,H_4,H_6,E,H_5,J,L_0,N_0).              \tag{4}
\]
Their independence minor is \(10888888320\).

Only the uniform source \(t\) and target \(RS\) gauges are removed.
The source gauge sends \(P\mapsto P-\lambda hQ\),
\(Z\mapsto Z-\lambda hX\); the target gauge sends
\(T\mapsto T+\mu hR,\ D\mapsto D-\mu hS\).
They preserve every common cap. The relation
\[
                       RS=-6t-3Y^2-8Y
\]
allows removal of the \(t,Y^2\) directions. No other source or target
flow is used for normalization.

## The \(h^3\) equations, including every complex branch

Before any field elimination, the normalized first Hamiltonian has
NINE independent coefficients:
\[
 G=aY+bR^2T+cH_4+dH_6+e_{\rm old}E+fH_5
             +\eta J+jL_0+kN_0.                         \tag{5}
\]
All 437 second-coefficient columns are retained. Seven exact
rational functionals on the six \(h^3\) equations give
\[
\begin{gathered}
d^2-250df+35625f^2=0,\qquad e_{\rm old}j-k^2/6=0,\\
e_{\rm old}k+\tfrac23jk+\tfrac23k^2=0,\qquad fj=fk=0,\\
\eta j+\tfrac12\eta k+\tfrac1{36}jk+\tfrac1{72}k^2=0,
 \qquad (j+k/2)^2=0.                                   \tag{6}
\end{gathered}
\]
Their supports are \(413,329,329,329,329,81,81\) rows.

Over a FIELD the last square forces \(j=-k/2\). Set
\[
 k=-6z,\quad j=3z,\quad e_{\rm old}=e+2z,\quad
 K_*=3L_0+2E-6N_0=X^2B .
\]
Then (5) becomes
\[
 G=aY+bR^2T+cH_4+dH_6+eE+fH_5+\eta J+zK_* .            \tag{7}
\]
The necessary and sufficient FIELD conditions through \(h^3\) are
\[
 d^2-250df+35625f^2=0,\qquad zd=ze=zf=0.                 \tag{8}
\]
For example \(zf=0\) follows directly from (6); multiplying the conic
by \(z^2\) then proves \(zd=0\) over fields. Conversely substituting
\(j=3z,k=-6z,e_{\rm old}=e+2z\) in (6) proves all seven equations.

The witness retains a complete rational polynomial second family
modulo (8). Replay reduces \(d^2\) to \(250df-35625f^2\), and kills
monomials divisible by \(zd,ze,zf\). It checks every coefficient of
\[
 \mathcal A_2F_2
   =-\mathcal P_1(F_1,F_1)-\mathcal P_3(f,f)/24.          \tag{9}
\]
The six homogeneous freedoms, with parameters
\((u,v,\gamma,\delta,\chi,\psi)\), are respectively \(K_G\) for
\[
 R^3,\ R^4S,\ X^3,\ R^3(Y+2),\ X^4Q,\ R^3J_4,\qquad
 J_4=9Y^4+66Y^3+36tY^2+150tY+171t^2-380t.              \tag{10}
\]
Their independence determinant is \(2597584896\). Thus no
second-coefficient freedom has been omitted.

The conic in (8) includes both
\(d=(125+100q)f\) and \(d=(125-100q)f\), \(q^2=-2\),
after a splitting field extension. It is not discarded for lack of
rational points.

## Complete field cover through \(h^4\)

The next equation uses every one of the 307 third-coefficient columns:
\[
\begin{split}
\mathcal A_3F_3={}&-\mathcal P_1(F_1,F_2)-\mathcal P_1(F_2,F_1)\\
 &-[\mathcal P_3(f,F_1)+\mathcal P_3(F_1,f)]/24.          \tag{11}
\end{split}
\]
Exact necessary functionals give \(f^2(d+375f)=0\).
If \(f\ne0\), then \(d=-375f\), and the conic becomes
\(270000f^2=0\), a contradiction. Thus \(f=d=0\), for both complex
signs and their intersection. Further exact functionals give
\[
 e^2\eta=0,\quad \psi z=0,\quad
 \psi(22c+45\eta)=0,\quad e\psi=0,\quad
 e(18ce+\chi)=0.                                      \tag{12}
\]
The latter three functionals are applied AFTER \(d=f=0\).
Since we work over fields, (8) and (12) imply
\(e\eta=ez=e\psi=0\).

The remaining cover is transparent in three cases.

If \(e\ne0\), then \(\eta=z=\psi=0,\ \chi=-18ce\).
Put \(K_E=\delta+(976c^2+216ce)/27\).
Necessary functionals give
\[
                 (a-544c/27)K_E=0,\qquad cK_E=0.
\]
Hence either \(K_E=0\) (Efull), or \(a=c=0\) (Ezero).

If \(e=\psi=0\), put
\[
\begin{split}
L&=a-544c/27,\quad t_*=z-6\eta,\\
K&=\delta+976c^2/27+192c\eta+488cz/9+144\eta z+33z^2,\\
\widehat\chi&=\chi+100cz/3-\tfrac92z(49z-288\eta),\\
M_*&=2c+72\eta-9z .
\end{split}
\]
Three necessary equations are
\[
 KL+576zt_*^2=0,\qquad M_*K=0,\qquad
                   8t_*\widehat\chi+3Kz=0.             \tag{13}
\]
They are constant invertible row combinations of the remaining
three large coefficient constraints; no division by a parameter
enters (13). For example, with \(T_*=KL+576zt_*^2\), the original
three constraints at positions 1, 2, and 6 are
\[
 \frac{297T_*+832KM_*+7488Kz+19968t_*\widehat\chi}{19968},
 \quad-\frac{297T_*+832KM_*}{7488},\quad\frac{9T_*}{4160}.
\]
If \(K=0\), the first equation gives \(z=0\) or \(t_*=0\).
When \(z=0\), (13) gives \(\eta\chi=0\): \(\chi=0\)
belongs to A0, and \(\eta=0\) belongs to HZres.
When \(t_*=0\), one obtains HZres. If \(K\ne0,t_*=0\),
(13) forces \(z=\eta=c=L=0\), giving Bzero.
Finally \(K t_*\ne0\) gives Bmain by solving (13).

If \(\psi\ne0\), then \(e=z=0,\ \eta=-22c/45\). Put
\[
 Q_*=\delta-7792c^2/135+908\psi/3,\qquad L=a-544c/27.
\]
Necessary equations give
\[
 LQ_*=0,\qquad
       138112cQ_*=187200L\psi .                         \tag{14}
\]
Because \(\psi\ne0\), these force \(L=0\).
The remaining equations give \(c\chi=0\).
If \(c=0\) the point lies in Bzero; otherwise
\(Q_*=\chi=0\), giving Cspec.

Here is the resulting exhaustive cover. Every row has \(d=f=0\).
Unmentioned variables among \(a,b,c,e,\eta,z,\delta,\chi,\psi\)
are free, subject only to the displayed formulas. In every row
\(u,v,\gamma\) and the third-kernel parameter \(w\) are free.
\[
\begin{array}{c|l}
\text{family}&\text{conditions}\\ \hline
\mathrm{Efull}&\eta=z=\psi=0,\quad
 \delta=-(976c^2+216ce)/27,\quad\chi=-18ce\\
\mathrm{Ezero}&a=c=\eta=z=\chi=\psi=0\\
\mathrm{A0}&e=z=\chi=\psi=0,\quad
 \delta=-976c^2/27-192c\eta\\
\mathrm{HZres}&e=\psi=0,\quad z=6\eta,\quad
 \delta=-976c^2/27-1552c\eta/3-2052\eta^2\\
\mathrm{Bzero}&a=c=e=\eta=z=0\\
\mathrm{Cspec}&e=z=\chi=0,\quad a=544c/27,\quad
 \eta=-22c/45,\quad\delta=7792c^2/135-908\psi/3\\
\mathrm{Bmain}& e=\psi=0,\quad t_*K\ne0,\quad
 \eta=(z-t_*)/6,\quad c=6t_*-3z/2,\\
 & a=544c/27-576zt_*^2/K,\\
 & \delta=K-976c^2/27-192c\eta-488cz/9-144\eta z-33z^2,\\
 & \chi=-3Kz/(8t_*)-100cz/3+\tfrac92z(z+48t_*)
\end{array}                                                   \tag{15}
\]
Rows overlap intentionally. Efull and Ezero are allowed also at
\(e=0\); Cspec and Bzero also at \(\psi=0\).
Bmain uses precisely the declared nonzero denominators \(t_*,K\);
all their zero loci are covered by the other rows.

The witness supplies an actual particular \(F_3\) on EVERY row of
(15), and replay verifies all six equations (11), coefficient by
coefficient, including every Laurent parameter coefficient in
Bmain. Since
\[
 \ker\mathcal A_3
    =\Bbbk K_{R^6},\qquad K_{R^6}=(0,0,-6R^5,0),
\]
the retained \(F_3+wK_{R^6}\) is the ENTIRE third family for the
specified \(F_1,F_2\). The rank-306 minor and this nonzero kernel
prove completeness. Thus (15), rather than only a necessary
superset or a collection of sampled points, classifies all
field-valued lower jets through \(h^4\).

## Every row fails at \(h^5\)

All 206 fourth-coefficient columns are allowed in
\[
\begin{split}
\mathcal A_4F_4={}&-\mathcal P_1(F_1,F_3)
 -\mathcal P_1(F_2,F_2)-\mathcal P_1(F_3,F_1)\\
 &-[\mathcal P_3(f,F_2)+\mathcal P_3(F_1,F_1)
                  +\mathcal P_3(F_2,f)]/24\\
 &-\mathcal P_5(f,f)/1920 .                             \tag{16}
\end{split}
\]
In particular \(\mathcal P_3(F_1,F_1)\) is INCLUDED.
For each row of (15), exact rational row functionals annihilate
ALL 206 columns and have these values on its complete right side:
\[
\begin{array}{c|l}
\text{family}&\text{necessary values, all required to vanish}\\ \hline
\mathrm{Efull,Ezero,A0,Cspec}&1\\
\mathrm{Bzero}&\chi^2,\quad1-768\chi\gamma/49\\
\mathrm{Bmain}&K^2\\
\mathrm{HZres}&\eta^4,\quad(27a-544c)^4
\end{array}                                                   \tag{17}
\]
Bzero is impossible: \(\chi^2=0\) forces \(\chi=0\), leaving 1.
Equivalently, with \(q=768\chi\gamma/49\),
\(1=(1+q)(1-q)+(768\gamma/49)^2\chi^2\).
Bmain is impossible since \(K\ne0\).

For HZres, (17) first forces
\(\eta=0,\ a=544c/27\). Recomputing the ENTIRE right side after
only these justified field substitutions gives further functionals
\[
                  c^4,\qquad\chi^2,\qquad
                              1-768\chi\gamma/49 .      \tag{18}
\]
These again contradict field solvability. No parameter other than
a proved field consequence is removed at this stage.

The supports of the first obstruction functionals are respectively
164, 162, 162, \((160,161)\), \((167,168)\), 163, and 161
for Efull, Ezero, A0, HZres, Bzero, Cspec, Bmain.
The three specialized HZres supports are \((158,159,160)\).
Replay checks every coefficient, not only values at parameter
samples. In particular all \(b,u,v,\gamma,\delta,\chi,\psi,w\)
freedoms that occur on their rows remain present.

As an independent positive control, the checker converts the actual
HZres lower point
\[
 (a,b,c,\eta,u,v,\gamma,\chi,w)=(0,0,0,1,0,0,0,0,0)
\]
to normal PBW order and verifies all six relations through \(h^4\)
using separate PBW multiplication. This includes \(z=6\) and the
new \(K_*\) direction; it is not the proof of (17).

This proves the universal cap-26 FIELD obstruction through \(h^5\).
The positive cap-27 point cited at the start then proves the sharp
minimum for this finite graded stage. No evenization theorem,
cap-growing source normalization, unrestricted scalar theorem, or
matrix obstruction is inferred.
