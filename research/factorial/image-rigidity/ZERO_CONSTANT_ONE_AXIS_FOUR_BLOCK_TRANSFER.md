# Four-block residue transfer for a zero-constant binary quadratic

Use the one-axis normalization
\[
 f=x^2/2+\beta xy+\gamma y^2/2+dx+ey,\qquad
 \bar\beta=1,\qquad \bar\gamma=\bar d=\bar e=0
\]
over a valued extension of \(\mathbb Q_2\), with \(v(2)=1\). Put \(\sigma=v(e)\). Consider a candidate block of five zero factorial moments starting at \(n\ge1\).

**Theorem.** The five-zero block is impossible for \(3/2<\sigma<7/4\). At \(\sigma=3/2\) it is impossible unless
\[
                   v(\gamma)=v(e),\qquad
                  n\equiv0\pmod4,\qquad \overline{\gamma^2/8}=1.
\]
There is no restriction on the residue-field degree. The equal-axis residue torsion constraint from [unit-pole removal](ZERO_CONSTANT_ONE_AXIS_UNIT_POLE_REMOVAL.md) continues to apply in the surviving case.

## One more valuation in the symmetric approximation

Use the [symmetric approximation](ZERO_CONSTANT_ONE_AXIS_SYMMETRIC_APPROXIMATION.md), writing \(\theta=v(\gamma)-\sigma\), \(u=\beta-1\), \(t=d-e\), and \(w=(n+3)\sigma+\theta-1\). The candidate gives \(v(d)=\sigma\), \(0\le\theta\le1\), and \(v(u)=v(t)=w\). Use \(r_i,z_i\) for the axis ratio and normalized weighted forcing from the [unit-pole proof](ZERO_CONSTANT_ONE_AXIS_UNIT_POLE_REMOVAL.md#exact-interior-and-weighted-axis-equations). In its last force equation, at \(i=n+4\), one has
\[
                u+t=(\gamma-\beta)z_i-u(r_i+d-1).
\]
Here \(v(z_i)\ge i\sigma+\theta>w+1\), from the axis bounds in the symmetric approximation. Since \(\sigma\ge1\), the recurrence for \(r_i\) gives \(v(r_i-1)\ge1\). It follows that \(v(u+t)\ge w+1\).

Write \(P_i\) and \(P_i^0\) for the normalized moments of \(f\) and the equal-linear phase \(f_0=x^2/2+xy+\gamma y^2/2+e(x+y)\). Their coefficient polynomials in the natural parameters are integral. At \(f_0\) their two first derivatives in \(\beta,d\) are both one modulo the integer two for every \(i\ge1\). Indeed all \(\gamma,d,e\) there are in \(2\mathcal O\), and the two derivatives of the universal matching series at \((\beta,\gamma,d,e)=(1,0,0,0)\) both give \(T/(1+T)\) modulo two. Taylor expansion in \(u,t\) consequently gives
\[
                       v(P_i-P_i^0)\ge w+1,                  \tag{1}
\]
since the linear term agrees with \(u+t\) to this precision and \(2w\ge w+1\). This uses divisibility by the integer two and also holds in a ramified field.

## The paired transfer

The exact symmetric recurrence is
\[
 P_i^0=2iP_{i-1}^0+e^2P_{i-2}^0+B_i+eB_{i-1},
 \qquad P_0^0=1,\quad P_1^0=2+\gamma+2e ,
\]
where \(B_i\) is the \(y\)-axis moment divided by \(i!\). Put
\[
 E_m=P_{2m}^0/e^{2m},\qquad O_m=P_{2m+1}^0/(2e^{2m}),\qquad
 b_m=8m/e^2,\qquad q_m=2m+1 .
\]
The recurrence on the vector \((E_m,O_m)^{\mathsf T}\) has matrix and forcing
\[
 A_m=\begin{pmatrix}1&b_m\\q_m&1+q_m b_m\end{pmatrix},\qquad
 F_m=\begin{pmatrix}C_m\\q_m C_m+D_m\end{pmatrix},\quad
 C_m=\frac{B_{2m}+eB_{2m-1}}{e^{2m}},\quad
 D_m=\frac{B_{2m+1}+eB_{2m}}{2e^{2m}} .
\]
The \(C_m\) are integral and the \(D_m\) have positive valuation. With \(g=\overline{\gamma/e}\), let \(c_i\) denote the reductions of \((B_i+e B_{i-1})/e^i\) at positive indices. Thus \(\bar C_m=c_{2m}\). These \(c_i\) are all zero when \(g=0\). When \(g\ne0\) they are \(\xi^i+\xi^{-i}\), with \(\xi+\xi^{-1}=g\), by the axis-series identity in the symmetric approximation.

Take \(3/2<\sigma<7/4\), put \(z=e^2/8\) and \(\delta=v(z)=2\sigma-3\in(0,1/2)\). At even \(m\) use the state \((z E_m,O_m)\). Its two-step transition is \(A_mA_{m-1}\) conjugated by \(\operatorname{diag}(z,1)\). Here \(v(b_m)\ge1-\delta\), \(v(b_{m-1})=-\delta\) and \(v(q_m+q_{m-1})=v(4m)\ge3\). The top-left entry reduces to one, and the top-right reduces to one from \(z b_{m-1}=m-1\). The bottom-left is
\[
        (q_m+q_{m-1}+q_m b_m q_{m-1})/z ,
\]
of valuation at least \(1-2\delta>0\). The bottom-right is
\[
      1+(q_m+q_{m-1})b_{m-1}+q_m b_m(1+q_{m-1}b_{m-1}),
\]
which reduces to one, since each correction has positive valuation. The scaled forcing \(F_m+A_m F_{m-1}\) reduces to \((0,c_{2m})^{\mathsf T}\): the first coordinate gains \(z\), and in the second all terms from the preceding even forcing combine with \(q_m+q_{m-1}\) or \(b_m\). All odd forcing gains disappear in reduction.

Consequently the reductions \((U_r,V_r)=(\overline{z E_{2r}},\overline{O_{2r}})\) obey
\[
            U_r=U_{r-1}+V_{r-1},\qquad V_r=V_{r-1}+c_{4r},
                      \qquad U_0=0,\quad V_0=1.             \tag{2}
\]
In the intermediate odd step use \((z E_{2r+1},z O_{2r+1})\). Its reduction is \((U_r+V_r,U_r+V_r)\), because \(z b_{2r+1}=2r+1\) and the factors of \(z\) on all remaining terms have positive valuation. Thus the sequence of scaled residues in the \(r\)-th block is
\[
                           U_r,\ V_r,\ U_r+V_r,\ U_r+V_r.     \tag{3}
\]
The initial values and the integral transitions justify taking all these residues. At each index \(k\) the normalization in (3) divides \(P_k^0\) by a scalar of valuation strictly below \(k\sigma\). By (1), at \(n,n+1,n+2,n+3\) the base moments have valuation at least \(w+1=(n+3)\sigma+\theta\). All four scaled residues there must vanish.

## Exclusion of the four residue positions

If \(g=0\), all \(c_{4r}=0\) and \(V_r=1\); each four-index interval contains a \(V\) position, a contradiction.

If \(g\ne0\), the residue torsion condition is \(\xi^{n+2}=1\). The exact ratios of \(\xi,\xi^{-1}\) have the same odd order as \(\xi\). Equation (2) gives
\[
                  V_r=\sum_{s=-r}^r \xi^{4s},\qquad
                  U_r=\sum_{i=0}^{r-1} V_i .                 \tag{4}
\]
In particular \(V_r=0\) implies \(\xi^{2r+1}=1\). Under this last equality, \(U_r=1/(\xi^2+\xi^{-2})\ne0\). To see this, put \(\eta=\xi^2\) and use
\[
 V_i=\frac{\eta^{-2i}+\eta^{2i+2}}{1+\eta^2}.
\]
Over \(0\le i<r\), the two geometric sums in its numerator are \(\eta^2/(1+\eta)\) and \(\eta/(1+\eta)\), whose sum is \(\eta\). This handles \(n=4r\). At \(n=4r+1\), vanishing forces \(V_r=0\), incompatible with the torsion \(\xi^{4r+3}=1\). At \(n=4r+2\), the last of the four residues is \(V_{r+1}\), whose vanishing is incompatible with \(\xi^{4r+4}=1\). At \(n=4r+3\), \(V_{r+1}=0\) is incompatible with \(\xi^{4r+5}=1\). In each incompatible pair the two exponents are coprime, while \(\xi\ne1\). This proves the open band.

## The boundary \(\sigma=3/2\)

Now \(z\) is a unit and the same two-step matrix reduces to \(\begin{pmatrix}1&1\\0&1\end{pmatrix}\). Its forcing gains a first coordinate \(\bar z(c_{4r-2}+c_{4r})\), and the initial first coordinate is \(\bar z\). In an intermediate odd step the two residues are \(F_r\) and \(\bar z V_r+F_r\), where \(F_r=U_r+V_r+\bar z c_{4r+2}\). Thus the first three residues of the proposed block, which vanish by (1), already imply \(g\ne0\): otherwise \(V_r=1\), and the two intermediate residues cannot both vanish.

Formula (4) for \(V_r\) is unchanged. At \(n=4r+1,4r+3\), the same incompatible-order argument uses a \(V\) among the first three residues. At \(n=4r+2\), both intermediate residues vanish, giving \(V_r=0\), incompatible with \(\xi^{4r+4}=1\). Thus \(n=4r\). Now \(V_r=0\) and (4) still gives \(\sum_{i<r}V_i=1/(\xi^2+\xi^{-2})=1/g^2\), whereas
\[
                 U_r=\sum_{i<r} V_i
                        +\bar z\left(1+\sum_{i=1}^{2r} c_{2i}\right).
\]
The parenthesis is \(\sum_{s=-2r}^{2r}\xi^{2s}=1\) when \(\xi^{2r+1}=1\), by another finite geometric sum. The condition \(U_r=0\) therefore gives \(\bar z g^2=1\), which is exactly \(\overline{\gamma^2/8}=1\).

## Exact checks

The [finite replay](verify_zero_constant_one_axis_four_block_transfer.py) checks 2,352 paired transfer entries, 7,376 residue sums and 22 matching first-order coefficients.
