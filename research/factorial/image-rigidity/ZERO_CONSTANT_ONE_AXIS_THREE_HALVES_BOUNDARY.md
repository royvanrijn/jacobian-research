# The three-halves boundary for a zero-constant one-axis quadratic

In the one-axis normalization
\[
 f=x^2/2+\beta xy+\gamma y^2/2+dx+ey,\qquad
            \bar\beta=1,\quad\bar\gamma=\bar d=\bar e=0 ,
\]
over a finite extension of \(\mathbb Q_2\), put \(v(2)=1\).
**Theorem.** When \(v(e)=3/2\), five consecutive positive-index factorial
moments of \(f\) cannot all vanish. The proof works also after
passing to finite extensions to split the residue polynomials.

**Corollary.** Any five-zero block in one-axis resonance requires
\(e\ne0\) and \(v(e)\ge7/4\). The
[symmetric approximation](ZERO_CONSTANT_ONE_AXIS_SYMMETRIC_APPROXIMATION.md)
gives the lower range with one exception at valuation one, which is
excluded by the [valuation-one lift](ZERO_CONSTANT_ONE_AXIS_UNIT_LINEAR_SCALE.md).
The theorem above supplies \(3/2\), and the
[paired transfer](ZERO_CONSTANT_ONE_AXIS_FOUR_BLOCK_TRANSFER.md)
covers \(3/2<v(e)<7/4\).

The [paired residue comparison](ZERO_CONSTANT_ONE_AXIS_FOUR_BLOCK_TRANSFER.md#the-boundary-sigma32)
reduces a candidate start to \(n=4r\), \(r\ge1\), with
\(v(\gamma)=v(d)=v(e)=\sigma=3/2\) and
\(\overline{\gamma^2/8}=1\). We use its improvement of the
symmetric approximation by one integer valuation digit.
The next nonzero term comes from the paired odd moment.

## The forced value of one axis term

Put \(j=n+2\), and retain the notation of
[unit-pole removal](ZERO_CONSTANT_ONE_AXIS_UNIT_POLE_REMOVAL.md#exact-interior-and-weighted-axis-equations):
\[
 \rho=\frac{\gamma d-\beta e}{\gamma-\beta},\qquad
 z_i=\frac{B_i+\rho B_{i-1}}{A_{i-1}},\quad
 r_i=\frac{A_i}{A_{i-1}},\quad q_i=r_{i+1}-r_i .
\]
Here \(A_i,B_i\) are the normalized one-variable axis
moments, and \(A_i\) are units. The three weighted axis
equations hold at \(j,j+1,j+2\). Write
\(w=(j+1)\sigma-1=v(\beta-1)=v(d-e)\), by the
[symmetric approximation](ZERO_CONSTANT_ONE_AXIS_SYMMETRIC_APPROXIMATION.md).
One has \(v(q_i)\ge1\) and \(v(q_i-2)\ge2\sigma+1\)
for \(i\ge2\), by the exact slope recurrence in that proof.
Taking the second difference of the three equations gives
\[
                      v(z_{j+2}-2z_{j+1}+z_j)
                                \ \ge\ w+2\sigma+1.
\]
The equal-valuation axis estimates give
\(v(z_{j+2})\ge(j+2)\sigma\) and \(v(z_{j+1})=(j+1)\sigma\).
Consequently \(2z_{j+1}\) determines
\(v(z_j)=(j+1)\sigma+1\). The first estimate in this
comparison is strict relative to this value because
\(\sigma>1\).

Put \(g=\gamma/e\), and at each positive index set
\[
                 c_i(g)=\frac{B_i+eB_{i-1}}{e^i}.
\]
Replacing \(\rho\) by \(e\) in \(A_{j-1}z_j\) adds
valuation at least \(j\sigma+w\).
Indeed \(\rho-e=\gamma(d-e)/(\gamma-\beta)\) and
\(v(B_{j-1})\ge(j-1)\sigma\). Since \(w>\sigma+1\),
it follows that
\[
                          v(c_j(g))=\sigma+1=5/2.        \tag{1}
\]

Write \(S_i(g)=\lambda^i+\lambda^{-i}\), where
\(\lambda+\lambda^{-1}=g\). These are integer polynomials
in \(g\). For reference, the reduction of the exact axis
recurrence gives
\[
                 c_{2m}(g)\equiv(-1)^m S_{2m}(g)\pmod4,
 \qquad          c_i(g)\equiv S_i(g)\pmod2\quad(i\ge1).
                                                               \tag{2}
\]
To see the first congruence, write \(b_i=B_i/e^i\), with
\(b_0=1,\ b_1=g+1\) and
\(b_i=(2i-1)g b_{i-1}+b_{i-2}\).
In the even and odd steps the coefficients reduce to
\(-g,g\) modulo four, respectively. The two-step matrix
has determinant one and trace \(2-g^2\), so the sequence
of even-index sums \(b_{2m}+b_{2m-1}\) has that trace in its
recurrence. Its starting values are \(2,2-g^2\), using
the auxiliary value \(b_{-1}=1\). These are the starting
values and trace of \((-1)^m S_{2m}\). The second
congruence follows from the one-step recurrence modulo
two and the initial values \(2,g+2\) for the same sums.

Choose \(\lambda\) with residue \(\xi\), where
\(\xi+\xi^{-1}=\bar g\), and the two residue roots are
distinct since \(\bar g\ne0\). Such a root is in a finite
extension. The residue torsion condition gives
\(\xi^{j}=1\), with \(\xi\ne1\). Set \(d_*=2r+1\),
so \(j=2d_*\) and \(\xi^{d_*}=1\) because the residue
order is odd. By (1) and (2),
\(S_{2d_*}(g)\equiv0\pmod4\). Hence
\[
                      S_{d_*}(g)^2\equiv 2\pmod4 .
\]
For \(a=\lambda^{d_*}-1\), the identity
\(S_{d_*}(g)=2+a^2/\lambda^{d_*}\) gives \(v(a)=1/4\)
and \(\overline{a^4/2}=1\).
Thus, putting \(W=\sum_{i=-r}^r\lambda^{4i}\), the
geometric sum gives
\[
       v(W)=1,\qquad
       \overline{W/2}=\frac{\xi^{-4r}}{\xi^4-1}
             =\frac{\xi^2}{\xi^4-1}
             =\frac1{\bar g^{\,2}} .                       \tag{3}
\]
In fact \(\lambda^{4d_*}-1=(1+a)^4-1\), whose term \(a^4\)
has valuation one and whose other terms are higher.
The denominator \(\lambda^4-1\) is a unit.

## Paired odd recurrence at precision above one

For the base phase
\(f_0=x^2/2+xy+\gamma y^2/2+e(x+y)\), write \(P_i^0\)
for its moments divided by \(i!\), and put
\[
 E_m=P_{2m}^0/e^{2m},\quad O_m=P_{2m+1}^0/(2e^{2m}),\quad
 z=e^2/8,\quad o=e/2,\quad b_m=m/z,\quad q_m=2m+1 .
\]
In this section \(o_v(1)\) means terms of valuation strictly
above one. The recurrence from the
[paired comparison](ZERO_CONSTANT_ONE_AXIS_FOUR_BLOCK_TRANSFER.md#the-paired-transfer)
is
\[
 E_m=E_{m-1}+b_m O_{m-1}+C_m,\qquad
 O_m=O_{m-1}+q_m E_m+D_m,
 \qquad C_m=c_{2m}(g),\quad D_m=o\,c_{2m+1}(g).
\]
Both state components are integral at this boundary, and
\(v(D_m)\ge1/2\). At \(m=2\ell\), \(b_m\) has valuation
at least one, \(q_m\equiv1\pmod4\), and
\(q_{m-1}\equiv-1\pmod4\). Eliminating the intermediate
odd state gives
\[
\begin{split}
 O_{2\ell}={}&O_{2\ell-2}+C_{2\ell}
                   +D_{2\ell-1}+D_{2\ell}\\
 &+b_{2\ell}\bigl[
       (1-b_{2\ell-1})O_{2\ell-2}-E_{2\ell-2}-C_{2\ell-1}
                          \bigr]+o_v(1).                    \tag{4}
\end{split}
\]
Indeed \(O_m=(1+q_m b_m)O_{m-1}+q_m E_{m-1}+q_m C_m+D_m\).
Substituting the odd equation for \(O_{m-1}\) leaves
\(-b_m\) times \(E_{m-1}\) at this precision. Substituting
the even-component recurrence for \(E_{m-1}\) then
gives (4). The extra term \(b_mD_{m-1}\) has valuation
at least \(3/2\).

The square bracket in (4) has zero residue. To verify this,
put \(t=\ell-1,\ \eta=\xi^2,\ u=\eta+\eta^{-1}=\bar g^{\,2}
=\bar z^{-1}\), and use the even-state residues from the
paired comparison:
\[
 \begin{split}
 V_t&=\overline{O_{2t}}=\sum_{i=-t}^{t}\eta^{2i},\\
 U_t&=\overline{zE_{2t}}
       =\sum_{i<t}V_i+\bar z W_t,\qquad
 W_t=\sum_{i=-2t}^{2t}\eta^i.
 \end{split}
\]
The residue in question is
\[
             (1+u)V_t+u\sum_{i<t}V_i+W_t+
                    \eta^{2t+1}+\eta^{-2t-1}.
\]
Since \(uV_t=\eta^{2t+1}+\eta^{-2t-1}\) in characteristic
two, this becomes \(V_t+W_t+u\sum_{i<t}V_i\).
The first two terms leave the odd powers of \(\eta\)
between \(-2t\) and \(2t\). The sum of the endpoint
identities \(uV_i=\eta^{2i+1}+\eta^{-2i-1}\) over \(i<t\)
is exactly the same set of powers. This proves the
cancellation for each \(\ell\).

In (4), equation (2) gives
\(C_{2\ell}=S_{4\ell}(g)+o_v(1)\). For the odd forces,
the mod-two coefficientwise congruence suffices since
\(v(o)=1/2\). Therefore
\[
 D_{2\ell-1}+D_{2\ell}
    =o\bigl(S_{4\ell-1}(g)+S_{4\ell+1}(g)\bigr)+o_v(1)
    =o\,g S_{4\ell}(g)+o_v(1).
\]
Also \(O_0=1+og+o_v(1)\).
Summing (4) yields
\[
                 O_{2r}=(1+og)W+o_v(1).
\]
By (3), \(v(O_{2r})=1\).

Finally \(P_{n+1}^0=2e^n O_{2r}\). The first precision
improvement in the
[paired comparison](ZERO_CONSTANT_ONE_AXIS_FOUR_BLOCK_TRANSFER.md#one-more-valuation-in-the-symmetric-approximation)
and the assumed zero at \(n+1\) instead give
\(v(P_{n+1}^0)\ge(n+3)\sigma\). This requires
\(v(O_{2r})\ge3\sigma-1=7/2\), the contradiction.
