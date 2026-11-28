# The posterior lift can lower noise stability for balanced partitions

For each \((\sqrt7-1)/3<\rho<1\) there is an ordinary balanced
planar three-cell partition \(A\) such that
\(S_\rho(\mathcal L_s A)<S_\rho(A)\) for all \(s<1\) sufficiently
close to one. The construction modifies a translated tripod in a
distant region and retains exact balance. It concerns comparison on
all ordinary partitions and makes no maximizing claim about \(A\).
[Definition of the lift](GAUSSIAN_WEIGHTED_CENTROID_INEQUALITY.md#iteration-improves-first-moments)

## Finite-line partitions and the derivative at the identity

Let the interfaces of a planar partition consist of finitely many
segments or rays. Orient each \(e=ij\) from \(i\) to \(j\), let
\(N_e\) be the unit normal, and set \(d_e=e_i-e_j\) and
\(q_e=d_e/\sqrt2\). Use Gaussian arclength \(d\nu_e\).
Allow distinct segments of the same interface type to have different
normals. Define
\[
 \begin{split}
  \mathcal E_\rho(A)
     &=\sum_{e,f}(d_e\cdot d_f)(q_e\cdot q_f)
                     \iint M_\rho(x,y)\,d\nu_e(x)d\nu_f(y),\\
  \mathcal N_\rho(A)
     &=\sum_{e,f}(d_e\cdot d_f)
                     \iint M_\rho(x,y)N_e(x)\cdot N_f(y)\,
                                      d\nu_e(x)d\nu_f(y).
 \end{split}
\]
Differentiating the joint Gaussian density in \(\rho\) equals the
sum of its mixed spatial derivatives. Two integrations by parts
give \(S'_\rho(A)=\mathcal N_\rho(A)\).
[Correlation derivative and interface normalization](THREE_CELL_COUPLED_TRANSLATION_TEST.md#a-lower-bound-using-unoriented-interface-measures)

Write \(\epsilon=\sqrt{1-s^2}\) and view \(A\) as a cylinder in
\((x,z)\), with \(z\) the model label-plane coordinate. Put
\(e_s=1_{\mathcal L_sA}-1_A\). Averaging over \(z\) gives \(T_s1_A-1_A\).
Expanding the objective at the cylinder gives
\[
 S_\rho(\mathcal L_sA)-S_\rho(A)
    =2(S_{s\rho}(A)-S_\rho(A))
          +\sum_i\langle e_{s,i},T_\rho e_{s,i}\rangle .
\]
The last operator acts on both coordinate blocks at correlation
\(\rho\). The normal flux limit below gives
\[
 \lim_{s\uparrow1}
     \frac{S_\rho(\mathcal L_sA)-S_\rho(A)}{1-s^2}
             =\rho\bigl(\mathcal E_\rho(A)-S'_\rho(A)\bigr).
                                                               \tag{PC1}
\]

On a compact interior segment of an \(i,j\) interface put
\(x=a+\epsilon uN_e\), with \(a\) on the wall. The conditional
masses tend to \(\Phi(-u)\), \(\Phi(u)\), and zero for the spectator.
The model \(i\) indicator differs from the binary winner
\(1_{\{z_i-z_j\ge\beta_i-\beta_j\}}\) only where the spectator wins.
The conditional mass of this set tends to zero, so the binary
Gaussian CDF gives \(\beta_i-\beta_j\to\sqrt2 u\).
The integrated normal change therefore tends to
\[
  d_e\int_{\mathbb R}
       [1_{\{u<q_e\cdot z\}}-1_{\{u<0\}}]\,du
                      =d_e(q_e\cdot z).                         \tag{PC2}
\]
Consequently the limiting signed measure for
\(e_s(x,z)d\gamma(x)d\gamma(z)/\epsilon\) is
\(\sum_e d_e(q_e\cdot z)d\nu_e(x)d\gamma(z)\).

To justify this limit and its Mehler pairing, put
\(m_s(x)=1-(T_s1_A(x))_{A(x)}\). The flip from \(x\) to
\(sx+\epsilon G\) must cross one of the finitely many supporting
lines. For a line \(\ell(x)=n\cdot x-c=0\) use
\(\ell(sx)=s\ell(x)-(1-s)c\) and the Gaussian tail bound. This gives
\(m_s(x)\le C\sum_\ell \exp(-c_0\ell(x)^2/\epsilon^2)\)
for \(s\) near one. In the model lift
\(\int |e_s(x,z)|_1d\gamma(z)=2m_s(x)\).
Inserting \(e^{\lambda'|z|^2}\) with \(\lambda'<1/4\) bounds
that integral by \(C_{\lambda'}\sqrt{m_s(x)}\), by
Cauchy–Schwarz. The line envelopes integrate against
\(e^{\lambda'|x|^2}d\gamma(x)\) to \(O(\epsilon)\), by the
normal Gaussian integral. On radius-\(\delta\) disks around the
finitely many endpoints and corners the cost is
\(O(\epsilon\delta)\), since their tangential projections have
length at most \(2\delta\). On compact sets away from the actual
walls the flip probability is exponentially small. This gives
the required domination on regular wall charts and removes the
endpoints.

The Mehler density is bounded by
\(C_\rho e^{\lambda_\rho(|x|^2+|y|^2)}\) with
\(\lambda_\rho=\rho/[2(1+\rho)]<1/4\), in either coordinate block.
Taking \(\lambda_\rho<\lambda'<1/4\) gives uniform tails in the
pairing after division by \(\epsilon^2\). On compacts the kernel
is continuous and bounded, so the weak flux limit applies to the
product of the signed measures. Its \(z\)-pairing gives
\(\rho(q_e\cdot q_f)\). The linear term divided by \(\epsilon^2\)
tends to \(-\rho S'_\rho(A)\), which proves (PC1).

## Exact balance and the distant interfaces

Use coordinates \((x,y)\). Let label \(2\) of the translated
tripod occupy \(x>t,\ |y|<\sqrt3(x-t)\). For \(x<R\) use that
tripod. For \(x\ge R\) use label \(1\) above \(y=0\) and label
\(3\) below it. For \(R>1\) the remaining label \(2\) mass is
a continuous strictly decreasing function of \(t<R\).
At \(t=0\) it is below \(1/3\) and at \(t\to-\infty\) tends
to \(\Phi(R)>1/3\). Thus a unique \(t_R<0\) gives mass \(1/3\).
Reflection across \(y=0\) ensures the other two masses are
\(1/3\). Also \(t_R\to0\) since the truncation removes at
most \(\overline\Phi(R)\) from the strictly decreasing
untruncated sector mass.

Of the original interfaces two slanted segments between the
vertex and \(x=R\), and the leftward \(1,3\) ray, are retained.
The new interfaces are
\[
       V_+=\{(R,y):0<y<H_R\},\qquad
       V_-=\{(R,y):-H_R<y<0\},\qquad
       H=\{(x,0):x>R\},\quad H_R=\sqrt3(R-t_R).
\]
In the cyclic orientations \(12,23,31\), the normals on
\(V_+\) and \(V_-\) are opposite. Their two ordered terms in
\(\mathcal E_\rho-\mathcal N_\rho\) contribute
\(-m(V_+,V_-)\), where \(m(B,C)=\iint M_\rho\,d\nu_B\,d\nu_C\).
Pairs of old tripod pieces have zero difference, since their
normal Gram matrix equals the label-normal Gram matrix. The
difference is zero for two points on the same new piece.

With \(p_\rho\) the standard bivariate Gaussian density at
correlation \(\rho\),
\[
       m(V_+,V_-)=p_\rho(R,R)\,
                   \mathbb P(0<Y<H_R,\ -H_R<Y'<0),
\]
where \((Y,Y')\) has correlation \(\rho\).
The probability tends to a positive constant, and the density
is a positive constant times \(e^{-R^2/(1+\rho)}\).
The other terms are smaller on this scale:

- A pairing of \(H\) with a vertical piece is at most
  \(C_\rho\int_R^\infty p_\rho(R,x)\,dx
                 =O_\rho(p_\rho(R,R)/R)\).
  The first bound integrates the remaining vertical variable
  over the whole line. The second is the tail bound for \(x\)
  normal with mean \(\rho R\) and variance \(1-\rho^2\).
- A vertical piece and a retained slanted piece are bounded
  by the joint density of \(X_x=R,\ Y\cdot n=c_R\), extending
  both pieces to their supporting lines. Here \(X,Y\) are the
  correlated planar observations,
  \(|n_x|=\sqrt3/2\), and \(c_R=O(t_R)\).
  This gives exponential bound
  \(\exp[-R^2/(2(1-3\rho^2/4))+o(R^2)]\).
  For \(H\) and a slanted piece the corresponding conditions
  are \(X_y=0,\ X_x\ge R,\ Y\cdot n=c_R\).
  Given \(X_y=0,Y\cdot n=c_R\), the conditional \(X_x\)
  variance is at most \(1-3\rho^2/4\) and its mean is
  \(O(t_R)\). This gives the same exponential bound.
- A pairing of a new piece with the old leftward ray is
  bounded by \(\exp[-R^2/(2(1-\rho^2))+o(R^2)]\).
  For a vertical piece the horizontal variables give the
  density at \(X_x=R\) with \(Y_x\le t_R\). For \(H\)
  they give \(\mathbb P(X_x\ge R,Y_x\le t_R)\), bounded
  by the tail \(X_x-\rho Y_x\ge R-\rho t_R\).

For \(\rho>(\sqrt7-1)/3\) both exponential rates exceed
\(1/(1+\rho)\). The two conditions are
\(3\rho^2/2+\rho>1\) and \(\rho>1/2\), respectively.
Therefore \(\mathcal E_\rho(A)-S'_\rho(A)
=-m(V_+,V_-)+o(p_\rho(R,R))<0\) for all sufficiently large
\(R\). Fix such \(R\) and apply (PC1) as \(s\uparrow1\).

For a fixed \(s<1\) and the same nonsimplicial partition,
the separate perimeter comparison gives strict improvement
when the tested correlation tends to one. The two statements
take different limits.
[Perimeter improvement of the lift](GAUSSIAN_POSTERIOR_LIFT_PERIMETER_IMPROVEMENT.md)
