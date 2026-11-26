# Unit-pole removal in one-axis zero-constant quadratic windows

Let \(f=x^2/2+\beta xy+\gamma y^2/2+dx+ey\) over a valued extension of \(\mathbb Q_2\), with \(v(2)=1\). Assume \(\bar\beta=1\) and \(\bar\gamma=\bar d=\bar e=0\), where bars denote residues. Include zero coefficients by assigning them infinite valuation.

**Theorem.** Five consecutive factorial moments of \(f\) at positive indices cannot all vanish in any of the following cases.

- \(v(\gamma)<v(e)\).
- \(e\ne0\) and \(v(\gamma)>v(e)+1\). This includes \(\gamma=0\).
- \(e=\gamma=0\).

If \(v(\gamma)=v(e)<\infty\), let \(\zeta\) be the ratio of the two roots of \(Z^2-\overline{\gamma/e}Z+1\) in the algebraic residue field. A five-zero block starting at \(n\ge1\) requires
\[
                     \zeta^{n+2}=1\quad\text{or}\quad\zeta^{n+3}=1.
\]
The conclusions allow ramified coefficient fields. When the residue field is \(\mathbb F_2\), the equal-valuation condition has \(\operatorname{ord}(\zeta)=3\); in larger residue fields the ratio need not have order three.

## Exact interior and weighted axis equations

Define the normalized moments
\[
 P_j=\Lambda_2(f^j)/j!,\quad
 A_j=\Lambda_1((x^2/2+dx)^j)/j!,\quad
 B_j=\Lambda_1((\gamma y^2/2+ey)^j)/j!.
\]
Put \(\Delta=\gamma-\beta^2\), \(D=1+\gamma-2\beta\), and \(V=\gamma d^2-2\beta de+e^2\). Both \(D,\Delta\) are units. The exact recurrence for positive \(j\ge2\) is
\[
 \begin{split}
 D P_j={}&2\Delta j P_{j-1}+ V P_{j-2}\\
        &+(\gamma-\beta)B_j+(\gamma d-\beta e)B_{j-1}
              +(1-\beta)A_j+(e-\beta d)A_{j-1}.                \tag{1}
 \end{split}
\]
To derive it, write \(P,A,B\) for the three formal series of normalized moments, and \(N_x,N_y\) for the series of factorial moments with \(x,y\) inserted. The two factorial Stein identities are
\[
  P-B=T(N_x+\beta N_y+dP),\qquad
  P-A=T(\beta N_x+\gamma N_y+eP).
\]
Euler's identity gives \(N_x+N_y=2P+2TP'-T(dN_x+eN_y)\). Apply the adjugate of \(\begin{pmatrix}1&\beta\\ \beta&\gamma\end{pmatrix}\) in the Stein system, then substitute in the Euler identity. The cross terms linear in \(T\) cancel in the two products involving \(P\), giving (1).

The axis recurrences from the [matching calculation](QUADRATIC_STRONG_FACTORIAL_WINDOWS_FROM_MATCHINGS.md#2-the-universal-matching-generating-function) give \(A_0=B_0=1\), \(A_1=1+d\), \(B_1=\gamma+e\), and
\[
 A_i=(2i-1)A_{i-1}+d^2A_{i-2}, \qquad
 B_i=\gamma(2i-1)B_{i-1}+e^2 B_{i-2}.
\]
All \(A_i\) are units. Let \(r_i=A_i/A_{i-1}\) and \(q_i=r_{i+1}-r_i\). Then
\[
             q_1=2-d/(1+d),\qquad
             q_{i+1}=2-\frac{d^2q_i}{r_i r_{i+1}}.            \tag{2}
\]
In particular, if \(v(q_i)<1\), then \(d\ne0\) and \(v(q_i)=(2i-1)v(d)\). This is obtained by reading (2) backwards, since no cancellation with the term of valuation one is possible at a smaller valuation.

Set \(\rho=(\gamma d-\beta e)/(\gamma-\beta)\) and \(C_i=B_i+\rho B_{i-1}\). Five consecutive zero \(P_i\), starting at \(n\), imply for \(j=n+2\) the three equations
\[
            (1-\beta)r_i+(e-\beta d)
                +(\gamma-\beta)z_i=0,\quad i=j,j+1,j+2,
                \qquad z_i=C_i/A_{i-1}.                     \tag{3}
\]
If the two consecutive differences of the \(z_i\) have valuations \(L\) and \(L+\sigma\), with \(v(z_{j+1})\ge L+\sigma\) and \(\sigma>0\), then \(\beta\ne1\) and
\[
 v(q_{j+1})=v(q_j)+\sigma,\qquad
 v(1-\beta)= v(e-\beta d)=L-v(q_j).                           \tag{4}
\]
Indeed the differences of (3) give the first equality and \(v(1-\beta)=L-v(q_j)<L+\sigma\). The middle equation gives the remaining equality since \(r_{j+1}\) is a unit. By (2), \(v(q_j)<1\), since at valuation at least one the next difference has valuation exactly one. Thus \(\lambda=v(d)>0\) is finite and
\[
                       v(q_j)=(2j-1)\lambda.                \tag{5}
\]

## Two axis valuations with \(\sigma=v(\gamma)\le v(e)\)

If \(v(\gamma)<v(e)\), then \(v(B_i)=i\sigma\) by its recurrence and \(v(\rho)>\sigma\). Thus \(v(C_i)=i\sigma\) for all \(i\ge1\).

If \(v(\gamma)=v(e)\), normalize \(B_i\) by \(e^i\), and put \(g=\overline{\gamma/e}\). The residual generating function of \(C_i/e^i\) is
\[
        \frac{(1+T)^2}{1+gT+T^2}
                 =1+\frac{gT}{1+gT+T^2},
\]
taking the constant coefficient to be one. In fact
\(\rho-e=\gamma(d-e)/(\gamma-\beta)\) has valuation greater than \(\sigma\). The polynomial \(Z^2-gZ+1\) has distinct unit roots in the algebraic residue field. The coefficient at \(i>0\) is zero precisely when their ratio has \(i\)-th power one. If the claimed torsion condition fails, \(v(C_j)=j\sigma\) and \(v(C_{j+1})=(j+1)\sigma\), while \(v(C_{j+2})\ge(j+2)\sigma\).

In either situation (4) applies with \(L=j\sigma\). We claim \(\sigma=2\lambda\) with \((2j+1)\lambda\le1\). If \((2j+1)\lambda<1\), (2) proves this immediately. If \((2j+1)\lambda>1\), write \(\nu=(2j-1)\lambda\). Then \(\sigma=1-\nu\), and
\[
 j\sigma-\nu<\lambda,\qquad j\sigma-\nu<\sigma .
\]
The first follows from \(j\sigma-\nu-\lambda=j[1-(2j+1)\lambda]<0\), and the second from \(\nu>(2j-1)/(2j+1)>(j-1)/j\). This contradicts the valuation of \(e-\beta d\) in (4), since \(v(e)\ge\sigma\). In the equality case \((2j+1)\lambda=1\), (2) gives \(\sigma\ge2\lambda\); then \(v(e-\beta d)=\lambda\) and (4) gives \(\sigma=2\lambda\).

Use the exact binary [matching series](QUADRATIC_STRONG_FACTORIAL_WINDOWS_FROM_MATCHINGS.md#2-the-universal-matching-generating-function)
\[
 P(T)\ \equiv\
 \frac{1+(\beta+d+e)T+deT^2}
      {(1+T+d^2T^2)(1+\gamma T+e^2T^2)}
                                                    \pmod{2}. \tag{6}
\]
Here the congruence means valuation at least one, even in a ramified field. The matching proof is a universal polynomial reduction modulo the integer two, so it applies over the quotient of the valuation ring by \(2\), which can have nilpotents. It does not merely give agreement of residues modulo the maximal ideal.

Factor \(1+T+d^2T^2=(1+aT)(1+uT)\), with \(a\) a unit and \(v(u)=2\lambda=\sigma\), using the simple residue root for \(a\). Factor \(1+\gamma T+e^2T^2=(1+bT)(1+cT)\) over a valued extension. If \(v(e)>\sigma\) the two valuations are \(\sigma\) and greater than \(\sigma\), including infinity if \(e=0\). When \(v(e)=\sigma\), both are \(\sigma\). This follows by comparing the valuations of their sum and product.

Multiply (6) by \(1+aT\). Of the remaining poles, either two or three of \(u,b,c\) have valuation \(\sigma\), with any other having higher valuation. At coefficient \(n+i\), \(i=1,2,3\), the numerator term of degree one has unit coefficient and yields valuation at least \((n+i-1)\sigma\). After division by \(u^{n+i-1}\) its residue is, up to a unit, the complete homogeneous polynomial of degree \(n+i-1\) in the minimum-valuation pole parameters divided by \(u\). The other numerator terms are higher, since \(v(de)>\sigma\). With two or three unit parameters these complete homogeneous polynomials have no two or three consecutive zeros, respectively. Their geometric generating function has denominator of that degree with nonzero leading coefficient, so consecutive zeros would propagate back to the constant coefficient one.

But \((n+2)\sigma=2j\lambda<1\). Thus among \(n+1,n+2,n+3\) one coefficient on the right has valuation below one. On the left those coefficients vanish exactly by the moment hypothesis. This proves the strict gamma-dominant exclusion and the equal-valuation torsion condition.

## The strongly linear \(y\)-axis

Suppose \(e\ne0\) and \(v(\gamma)>v(e)+1\), and put \(\sigma=v(e)\). In \(B_i/e^i\), the axis recurrence differs from that of the constant unit sequence by terms of valuation greater than one. Since \(\rho/e=1+\) terms of valuation greater than one, it follows that \(v(C_i)=1+i\sigma\) at every \(i>0\). Thus (4) applies with \(L=1+j\sigma\).

Again (5) holds. If \((2j+1)\lambda\le1\), (2) gives \(\sigma\ge2\lambda\), forcing \(v(e-\beta d)=\lambda\), whereas \(L-v(q_j)\ge1+\lambda\). If \((2j+1)\lambda>1\), then \(\sigma=1-(2j-1)\lambda\) and
\[
                  v(e-\beta d)=1+j\sigma-(2j-1)\lambda
                               =(j+1)\sigma.
\]
Unless \(\sigma=\lambda\), the left side is the smaller of \(\sigma,\lambda\), a contradiction. Hence \(\sigma=\lambda=1/(2j)\).

In (6), remove the unit pole as above. The small root from the first denominator has valuation \(2\lambda\), and the two roots from the second denominator have valuation \(\lambda\). Both normalized residues there are one since \(\gamma/e\) reduces to zero. Their complete homogeneous coefficients have no adjacent zeros. The numerator's quadratic coefficient has valuation \(2\lambda\), strictly higher than the pole scale \(\lambda\). Among the coefficients at \(n+1,n+2\) one thus has valuation \(n\lambda\) or \((n+1)\lambda<1\). This contradicts the exact zeros after removing the unit pole.

Finally if \(e=\gamma=0\), all positive \(B_i\) vanish. Three force rows (3) would give two adjacent \(q_i=0\) if \(\beta\ne1\), contrary to (2). For \(\beta=1,d\ne0\) the forcing is already nonzero in a single row. For \(\beta=1,d=0\) the phase is \((x^2+2xy)/2\) and all moments are \(2^j(j!)^2\).

The unresolved one-axis valuations include \(v(e)\le v(\gamma)\le v(e)+1\) outside the exclusions above. These statements concern the five-window. A support-specific shorter window requires its own argument.

## Exact checks

The [finite replay](verify_zero_constant_one_axis_unit_pole_removal.py) checks the interior and axis identities, the consecutive slope formulas, and matching reduction modulo two in ramified finite coefficient rings.
