# Symmetric approximation in the remaining one-axis residue

Use the notation of [unit-pole removal](ZERO_CONSTANT_ONE_AXIS_UNIT_POLE_REMOVAL.md) for
\[
 f=x^2/2+\beta xy+\gamma y^2/2+dx+ey,\qquad
 \bar\beta=1,\quad\bar\gamma=\bar d=\bar e=0 .
\]
Let \(n\ge1\) be a proposed start of five zero factorial moments, \(j=n+2\), \(\sigma=v(e)\), and \(\theta=v(\gamma)-v(e)\).

**Theorem.** A five-zero block in this class must satisfy \(e\ne0\), \(0\le\theta\le1\), and
\[
               v(d)=\sigma,\qquad
               v(\beta-1)=v(d-e)=(n+3)\sigma+\theta-1.          \tag{1}
\]
For \(\theta<1\) it also requires \(n\sigma\ge1\), while for \(\theta=1\) it requires \((n+1)\sigma\ge1\).

There is no such block for \(\sigma<3/2\), except possibly when \(\sigma=1\), \(\theta=0\), \(n\) is odd, and \(\overline{\gamma/2}=1\). At \(\sigma=3/2\) and \(\theta>0\) there is no block either. In the equal-valuation case \(\theta=0\), the residue-root order constraint from [unit-pole removal](ZERO_CONSTANT_ONE_AXIS_UNIT_POLE_REMOVAL.md) is still necessary.
The [valuation-one lift](ZERO_CONSTANT_ONE_AXIS_UNIT_LINEAR_SCALE.md)
excludes the exceptional case at \\(\\sigma=1\\).

## Valuations forced by the three axis equations

By the exclusions in the cited note, the remaining possibilities are \(\theta=0\) with its residue torsion condition, \(0<\theta<1\) with \(n\) even, or \(\theta=1\) with \(n\) odd and \(\overline{\gamma/(2e)}=1\). In each such case the weighted \(y\)-side sequence \(C_i\) satisfies
\[
       v(C_j)>j\sigma+\theta,\qquad
       v(C_{j+1})=(j+1)\sigma+\theta,\qquad
       v(C_{j+2})>(j+1)\sigma+\theta .
\]
These follow from the equal-valuation generating function or the even and odd first-order axis formula in the cited proof. In the equal-valuation case, if \(\zeta^j=1\) then \(\zeta^{j+1}\ne1\) because \(\zeta\ne1\).

Let \(z_i=C_i/A_{i-1}\) and \(q_i=r_{i+1}-r_i\) as in the cited proof. Write \(A=v(z_{j+1}-z_j)>j\sigma+\theta\) and \(B=v(z_{j+2}-z_{j+1})=(j+1)\sigma+\theta\). Since \(B\) is finite, the weighted axis equations imply \(\beta\ne1\). If \(A<B\) then \(v(q_j)<1\), and the slope recurrence gives \(v(q_j)=(2j-1)\lambda\) for \(\lambda=v(d)\). The weighted equation at \(j+1\) then yields
\[
                v(e-\beta d)=A-(2j-1)\lambda.
\]
If \(\sigma\ge2\lambda\), the right side is strictly greater than \(\lambda\), while the left side is \(\lambda\). Thus \(\sigma<2\lambda\). Remove the unit pole from the matching generating function as in the cited proof. The two poles of the \(y\)-axis alone have minimum valuation \(\sigma\), and \((n+1)\sigma<(2j-1)\lambda<1\). Their complete homogeneous coefficients do not vanish at consecutive indices, contradicting the exact zero coefficients at \(n+1,n+2\) after pole removal.

Consequently \(A\ge B\). The two axis differences force \(v(q_{j+1})\le v(q_j)\), so the slope recurrence gives \(v(q_j)\ge1\) and \(v(q_{j+1})=1\). This also holds if \(A=+\infty\), since then \(q_j=0\) and \(q_{j+1}=2\). Taking the difference at \(j+1,j+2\) and substituting in the equation at \(j+1\) now gives
\[
                  v(1-\beta)=v(e-\beta d)=B-1.                \tag{2}
\]
If \(\lambda\ne\sigma\), then \(B-1=\min(\lambda,\sigma)\). When \(\lambda>\sigma\), this gives \(j\sigma=1-\theta\), and the two \(y\)-axis poles remain visible before precision one at \(n+1,n+2\). When \(\lambda<\sigma\), it gives \(j\sigma<1\), and in the matching argument every minimum pole has valuation at most \(\sigma\). There are at most three such poles after removing the unit one. At least one of the coefficients at \(n+1,n+2,n+3\) is again visible before precision one. These contradictions give \(\lambda=\sigma\), and (2) then gives (1), since \((\beta-1)d\) has strictly higher valuation.

With \(\lambda=\sigma\), the \(y\)-axis poles alone have minimum valuation \(\sigma\). For \(\theta=0\), their complete homogeneous coefficient at degree \(n\) is nonzero in the residue field by \(\zeta^{n+2}=1\) and \(\zeta\ne1\). For \(0<\theta<1\) the normalized poles both have residue one and this coefficient is nonzero since \(n\) is even. In either case the coefficient at \(n+1\) after unit-pole removal gives \(n\sigma\ge1\). For \(\theta=1\) and odd \(n\), use the coefficient at degree \(n+1\) in the normalized two-pole sequence, obtaining \((n+1)\sigma\ge1\).

## Comparison with equal linear coefficients

Let \(P_i\) denote the normalized original moments and let \(P_i^0\) denote the normalized moments of
\[
         f_0=x^2/2+xy+\gamma y^2/2+e(x+y).
\]
The exact recurrence (1) of the preceding note specializes to
\[
 P_i^0=2i P_{i-1}^0+e^2 P_{i-2}^0+B_i+e B_{i-1},
 \quad P_0^0=1,\quad P_1^0=2+\gamma+2e .                     \tag{3}
\]
By integrality of the natural coefficient polynomials, \(P_i-P_i^0\) has valuation at least \(w=(n+3)\sigma+\theta-1\). We need a stronger comparison at \(\sigma\le1\).

Write \(R,R_0\) for the right side of the binary matching identity modulo two, evaluated at \(f,f_0\) respectively, and regard them as actual rational series with integral coefficients. Their numerators \(U,U_0\) satisfy
\[
  U-U_0=[(\beta-1)+(d-e)]T+e(d-e)T^2 .
\]
The errors \(P-R\) and \(P^0-R_0\) are universal coefficient polynomials divisible by the integer two. Their difference has coefficient valuations at least \(1+w\), because the two parameter sets differ only in \(\beta,d\), each by valuation \(w\).

Factor the \(x\)-axis denominators as \((1+aT)(1+uT)\) and \((1+a_0T)(1+u_0T)\) respectively, with \(a,a_0\) units and \(v(u)=v(u_0)=2\sigma\). Choosing the simple unit branches, \(v(a-a_0)=v(u-u_0)\ge w+\sigma\); subtract the equations for the two unit roots, whose divided difference is a unit, and use \(v(d^2-e^2)\ge w+\sigma\). The common denominator \(1+\gamma T+e^2T^2\) on the \(y\)-side factors into parameters of valuation \(\sigma\).

Multiply \(R-R_0\) by both unit factors. The resulting denominator has only parameters of valuation at least \(\sigma\). The numerator is
\[
      (1+a_0T)U(1+u_0T)-(1+aT)U_0(1+uT).
\]
Its constant coefficient is zero, and its coefficients of degrees \(1,2,3,4\) have valuations at least \(w,w,w+\sigma,w+3\sigma\), respectively. These follow by using the bounds on the numerator difference, the root differences, the values \(v(u)=v(u_0)=2\sigma\), and the valuation \(2\sigma\) of the quadratic coefficient of \(U_0\). Thus for \(k\ge2\), the coefficient of the filtered difference of rational series has valuation at least \(w+(k-2)\sigma\). Adding the error above gives valuation strictly greater than \((k-2)\sigma\) for \(k=n+2,n+3,n+4\), since \(1+w=(n+3)\sigma+\theta\). The corresponding coefficients of the filtered \(P\) itself vanish, by the original five-zero block.

When \(\sigma\le1\), equation (3) gives \(P_k^0/e^k\) integral at every \(k\). In the filtered \(P^0\), after division at \(k\) by \(e^{k-2}\), only the term \(aa_0P^0_{k-2}/e^{k-2}\) can survive reduction. The comparison therefore forces three consecutive zero residues of \(P^0_i/e^i\) at \(i=n,n+1,n+2\).

## The base residue recurrences

Put \(g=\overline{\gamma/e}\). The normalized axis forcing \((B_i+e B_{i-1})/e^i\) has residue \(c_i\), with generating function \((1+T)^2/(1+gT+T^2)\) if its constant coefficient is set to one. For \(0<\sigma<1\), the residue series of \(P_i^0/e^i\) is therefore \(1/(1+gT+T^2)\), since \(2i/e\) reduces to zero in (3). Its coefficients have no adjacent zeros. For \(\sigma=1,\theta>0\), the forcing residues at positive indices are zero, and the homogeneous residue recurrence has coefficient one on the second predecessor. Again adjacent zeros are impossible.

Suppose \(\sigma=1,\theta=0\), and let \(\xi,\xi^{-1}\) be the two roots with sum \(g\ne0\). Then \(c_i=\xi^i+\xi^{-i}\) for \(i>0\). Write \(p_i=\overline{P_i^0/e^i}\) and \(u=\overline{2/e}\). Equation (3) gives
\[
  p_{2m}=p_{2m-2}+c_{2m},\qquad
  p_{2m+1}=p_{2m-1}+u p_{2m}+c_{2m+1}.
\]
Here \(p_0=1,\ p_1=u+g\). Thus \(p_{2m}=\sum_{s=-m}^{m}\xi^{2s}\), which vanishes exactly when \(\xi^{2m+1}=1\). The ratio torsion from the preceding note gives \(\xi^{n+2}=1\), because a root of unity in characteristic two has odd order. Even \(n\) is impossible since then \(p_n=0\) would also give \(\xi^{n+1}=1\). For odd \(n=2m+1\), the same torsion gives
\[
     \sum_{s=0}^{m}c_{2s+1}=1,\qquad
     \sum_{s=0}^{m}p_{2s}=1/g
\]
by summing the two geometric series. Consequently \(p_n=1+u/g\), so a zero requires \(u=g\), or equivalently \(\overline{\gamma/2}=1\).

Finally consider \(1<\sigma\le3/2\). In this case \(w>(n+2)\sigma\), so the first three unfiltered base moments already vanish below their coefficient precisions. Normalize the even base moments by \(e^{2m}\) and the odd ones by \(2e^{2m}\); the resulting values are integral by (3). The even-normalized recurrence has coefficient \(8m/e^2\) on the preceding odd value, and the odd-normalized recurrence has a unit on the preceding even one. The odd axis forcing gains \(\sigma-1>0\). For \(\sigma<3/2\), the reduced recurrences for these even and odd values \(E_m,O_m\) are
\[
                 E_m=E_{m-1}+c_{2m},\qquad
                 O_m=O_{m-1}+E_m,\qquad E_0=O_0=1 .
\]
If \(\theta>0\), all \(E_m=1\). If \(\theta=0\), \(E_m\) is the same geometric sum as \(p_{2m}\) above. At even \(n\) it cannot vanish under \(\xi^{n+2}=1\), and at odd \(n=2m+1\) the value \(O_m\) is the same sum \(1/g\) computed above. This excludes the interval \(1<\sigma<3/2\). At \(\sigma=3/2,\theta>0\), the reduced even recurrence instead has the additional term \(\overline{8m/e^2}\,O_{m-1}\), and still has \(c_{2m}=0\). Any pair \(E_m=O_m=0\) then propagates backwards to \(E_0=O_0=0\), a contradiction. The three consecutive base residue zeros contain such a pair.

## Exact checks

The [finite replay](verify_zero_constant_one_axis_symmetric_approximation.py) checks 240 base and normalized recurrences, 192 filtered coefficient bounds and 3,728 finite-field geometric sums.
