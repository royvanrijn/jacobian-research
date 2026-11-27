# PBW transport and the final coefficient comparison

Two different orders must be kept separate. An edge-centroid mismatch
changes a coefficient at head minus \(N\). Under the monomial-face
hypothesis it does not create a coefficient at head minus one, for
\(N>1\). Separately, exact PBW multiplication contributes
\(N\)-dependent constants at head minus \(N+1\). Retaining those
constants changes the last rational comparison in the primary proof.

The exact first-component and elementary Schur calculations below
require no DC-pair. The final comparison is conditional on the
source's additional shared-head, body-tail and scalarization
constructions; those are not established here.

## Source conventions and raw transported components

The pinned source is Zheglov, arXiv:2410.06959v5,
[the cached PDF](sources/zheglov-2410.06959v5.pdf).
Definition (2.1), printed page 15, gives
\(\Phi_{N,1}(x)=x+\partial^N,\ \Phi_{N,1}(\partial)=\partial\).
The definitions on printed pages 10 and 17 give
\(\Gamma_1=x\partial\) and
\(\partial^p\Gamma_1=(\Gamma_1+p)\partial^p\).
Normal PBW ordering has \([\partial,x]=1\).
The homogeneous order in this calculation is \(j-i\) on
\(x^i\partial^j\).

Suppose fixed polynomials \(P_0,Q_0\) have rectangular support with
unique top-\(x\) faces \(x^a\partial^b\) and
\(\epsilon x^c\partial^e\), respectively. Here \(a,b,c,e\) are positive
integers and \(\epsilon\ne0\). Define fixed scalars
\[
 \rho_P=\frac{[x^{a-1}\partial^b]P_0}{a},\quad
 \rho_Q=\frac{[x^{c-1}\partial^e]Q_0}{\epsilon c},\quad
 A=[x^{a-1}\partial^{b-1}]P_0,\quad
 B=\frac{[x^{c-1}\partial^{e-1}]Q_0}{\epsilon}.
\]
No equality between the centroids \(\rho_P,\rho_Q\) is assumed.
For \(N>1\), put
\[
                         p_N=aN+b,\qquad q_N=cN+e.
\]
The two raw heads after \(\Phi_{N,1}\) are \(1,\epsilon\).
Every homogeneous component strictly between head and head minus
\(N\) is zero. The next two components are exactly
\[
\begin{split}
 (\Phi P_0)_{p_N-N}&=a\rho_P\,\partial^{p_N-N},\\
 (\Phi P_0)_{p_N-N-1}
   &=\left(a\Gamma_1+A+N\binom a2\right)\partial^{p_N-N-1},\\
 (\Phi Q_0)_{q_N-N}&=\epsilon c\rho_Q\,\partial^{q_N-N},\\
 (\Phi Q_0)_{q_N-N-1}
   &=\epsilon\left(c\Gamma_1+B+N\binom c2\right)
                                      \partial^{q_N-N-1}.     \tag{1}
\end{split}
\]

For the proof, the top monomial becomes
\((x+\partial^N)^a\partial^b\).
Its first \(x\)-term is \(a x\partial^{p_N-N}\).
Each ordered pair of factors contributes one contraction
\([\partial^N,x]=N\partial^{N-1}\), yielding
\(N\binom a2\partial^{p_N-N-1}\).
Terms with at least two \(x\)'s or two contractions have lower
homogeneous order. The original adjacent monomial contributes \(A\);
the original head-minus-\(x\) monomial contributes the centroid term.
All columns two steps left are at most head minus \(2N\), below the
last order in (1). This proves (1), and the same proof applies to
\(Q_0\).

Thus the raw scalar accompanying \(a\Gamma_1\) is
\(A_N^{\rm raw}=A+N\binom a2\), and its terminal counterpart is
\(B_N^{\rm raw}=B+N\binom c2\).
These replace the scalars described as \(N\)-independent \(a,b\)
in the first displayed input of Step 4 on printed page 72.
They are not absorbed into fixed coefficients of \(P_0,Q_0\).

## Centroids and the actual head-minus-one issue

An elementary Schur factor of order \(-N\) clears the first centroid
coefficient. Its coefficient is linear in \(\Gamma_1\); the change
of the terminal coefficient at head minus \(N\) is
\[
             \epsilon\left(c\rho_Q
                         -\frac{a q_N}{p_N}\rho_P\right).     \tag{2}
\]
It does not change the next component in (1): repeated products of
an order-\(-N\) term first contribute at head minus \(2N\), which
is below head minus \(N+1\) when \(N>1\).

If the terminal corner satisfies
\[
                              ae-bc=a-b,                      \tag{3}
\]
and \(\rho_Q=\rho_P+\eta\), expression (2) is
\[
                   \epsilon\left(c\eta
                         -\frac{a-b}{p_N}\rho_P\right).
\]
The one-root control's nonzero centroid difference therefore lives
at head minus \(N\). Meanwhile (1)'s gap already gives zero at
head minus one.

Conjugations \(1+O(\partial^{-N})\) preserve head minus one for
\(N>1\). The filtered Schur, scalarization and tail conditions
spelled out in the
[conditional Laurent note](CONDITIONAL_LAURENT_COEFFICIENT_INVARIANCE.md)
preserve that zero if their asserted constructions exist. Under these
hypotheses the final scalar Laurent data therefore have
\(\mathscr A_N=\epsilon,\ \mathscr B_N=0\), even with
\(\rho_P\ne\rho_Q\).
The [corrected exceptional invariant](EXCEPTIONAL_LAURENT_CORRECTION.md)
then equals the original distinguished coefficient: its
\(\mathscr A_N\mathscr B_N\) correction vanishes. No actual-pair
centroid equality follows from this observation.

## Elementary next Schur factor

At order \(-N-1\), use the polynomial part
\(\gamma_2\Gamma_1^2+\gamma_1\Gamma_1\) of the Schur coefficient.
The shift identity for \(\Gamma_1\) gives exactly
\[
 \gamma_2=-\frac{a}{2p_N},\qquad
 \gamma_1=\frac a2-\frac{A_N^{\rm raw}}{p_N}.                   \tag{4}
\]
Indeed these cancel \(a\Gamma_1+A_N^{\rm raw}\) in the first
component. The terminal component becomes
\[
 \epsilon\left[
 \left(c-\frac{a q_N}{p_N}\right)\Gamma_1+
 B_N^{\rm raw}+\gamma_2q_N^2+\gamma_1q_N
 +(\hbox{shift of the vector part})\right].                   \tag{5}
\]
Under (3), its \(\Gamma_1\)-coefficient is \((b-a)/p_N\).
Equations (4)--(5) are elementary Schur algebra. An asserted
vector-part normalization and its averaging remain separate inputs.

## Conditional corrected final comparison

Assume now that the source's subsequent body-tail decomposition and
auxiliary averaging have the form stated in Lemma 6.4 and Step 7,
printed pages 74--75. In its notation put \(a=dn,\ b=ln\) and
\(c=\widetilde M_{I-1}+1\). The terminal edge gives (3).
Let \(c_0\) denote its fixed shared-head constant. With the raw
components (1) retained, its normalized distinguished coefficient
must then be
\[
 \frac{\mathcal C_N}{\epsilon}
 =(b-a)\left(\frac12-\frac1{2p_N}-c_0\right)
    -\frac{a q_N^2}{2p_N}
    +\left(\frac a2-\frac{A+N\binom a2}{p_N}\right)q_N
    +B+N\binom c2.                                         \tag{6}
\]
This is a conditional substitution into those later formulas, not
a proof of the body-tail construction.

Direct algebra using (3) reduces (6) to
\[
 \boxed{\quad
 \frac{\mathcal C_N}{\epsilon}
     =C_\infty+
        \frac{(a-b)(ab-2A)}{2a(aN+b)}.\quad}                 \tag{7}
\]
One explicit \(N\)-independent constant is
\[
 C_\infty=B-\frac ca A+(b-a)(\tfrac12-c_0)
             +\frac{bc}{2}-c+\frac12
             -\frac{b(c-1)^2}{2a}.
\]
The \(N^2\) numerator terms from the two binomial contractions cancel
the apparent linearly growing dependence in the printed formula.
They cannot be dropped when \(N\) varies.

For \(a>b>0\), the remaining expression is independent of \(N\)
on an infinite set exactly when
\[
                                  A=ab/2.                    \tag{8}
\]
In original parameters the possible residual is
\[
             -\frac{(d-l)(2A-dln^2)}
                     {2dn(dN+l)} .
\]
Consequently, if all later application hypotheses and Laurent
invariance are supplied, this comparison would force (8).
It does not by itself give a contradiction or force the earlier
homogeneous-branch parameter to vanish. Condition (8) is precisely
the adjacent PBW coefficient of the Weyl-ordered leading corner.
For actual centered rectangular Weyl pairs, the
[exact indicial-mean theorem](POLYNOMIAL_WEYL_FIXED_AND_MOVING_SPECTRAL_SPLITTINGS.md#exact-indicial-means-and-the-mixed-weyl-coefficient)
now supplies the zero mixed Weyl corner directly. The Fourier change
exchanging the two rectangle directions preserves this zero condition.
Thus the numerical corner value is settled in that centered class.
The uncentered source data and the subsequent transport, body-tail and
scalarization hypotheses of the comparison remain separate obligations;
zero corner alone gives no contradiction.

## Exact controls and scope

[The Fraction checker](check_pbw_transport_final_coefficient.py)
performs independent normal PBW multiplication for tame shears,
including unequal centroid and adjacent coefficient controls. It
also verifies (7) as a seven-variable polynomial identity after
clearing denominators, and checks the constant and nonconstant cases
of (8). The [receipt](pbw_transport_final_coefficient_controls.json)
pins the original edition of this note, checker, primary PDF and
corrected Laurent source. Its historical hashes are retained.

The raw first-component assertion and the elementary Schur
calculations are exact. The passage to (6) is conditional on the
source's body-tail, common-head, commuting-curve, trace and auxiliary
conjugation bridges. This note constructs no actual DC-pair and
asserts no whole-paper refutation or DC1/DC2 resolution.
