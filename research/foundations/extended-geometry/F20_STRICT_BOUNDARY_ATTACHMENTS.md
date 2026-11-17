# Strict-boundary attachments for the corrected \(F_{20}\) Cox atlas

## Status and outcome

The exceptional Cox atlas now attaches to three strict-boundary families.
Exact controlled transforms certify

1. the final cusp divisor \(E_4\) meeting the strict \(r\)-boundary;
2. the triple-orbit divisor \(E_2\) meeting the strict \(d\)-boundary; and
3. the \(A\)-packet over the \(q\)-\(r\) tangency meeting the strict
   \(q\)-boundary.

The same calculation proves a general weighted Taylor--Cox attachment
theorem that treats unramified, simply ramified, and index-four boundary
scales uniformly.  It also isolates the \(q\)-node failure exactly: the
scaled \(h_q\) incidence fixes the exceptional slope but loses its root
coordinate there.  That edge must therefore be constructed in the
conductor-normalized chart, not in the present exceptional chart.

The exact checker is
[`verify_f20_strict_boundary_attachments.py`](../scripts/verify_f20_strict_boundary_attachments.py),
with generated certificate
[`f20_strict_boundary_attachments.json`](../artifacts/generated-results/f20_strict_boundary_attachments.json).

This result does not define the global Čech class.  The strict \(r\)
attachments at the triple and \(q\)-\(r\) packets, together with the
conductor-to-node/triple/transverse transitions, remain open.

## 1. Weighted Taylor--Cox attachment theorem

Let \(A\) be a regular local domain of characteristic zero, let \(\lambda\)
be a boundary scale on a finite normalized root-cover chart, and let
\(F\in A[Y]\).  Suppose a Cartier-compatible root centre \(c\) is defined
in a finite local incidence algebra over \(A\).  Fix \(m\geq 1\) and assume

\[
 \partial_Y^kF(c)\in(\lambda^{m-k}),
 \qquad 0\leq k<m.                              \tag{1.1}
\]

Taylor expansion after \(Y=c+\lambda W\) gives

\[
 F(c+\lambda W)=\lambda^m\Phi,\qquad
 F_Y(c+\lambda W)=\lambda^{m-1}\Phi_W.          \tag{1.2}
\]

Indeed, the term of order \(k<m\) in the expansion contains
\(\lambda^k\partial_Y^kF(c)\), hence is divisible by \(\lambda^m\); every
term of order at least \(m\) already contains that power.  Differentiating
the first identity with respect to \(W\) gives the second.

Now suppose a controlled transform satisfies

\[
 P=e^NF,
 \qquad X=\xi+e^rY.                              \tag{1.3}
\]

Then

\[
 P_X=e^{N-r}\lambda^{m-1}\Phi_W.                \tag{1.4}
\]

> **Weighted Taylor--Cox attachment theorem.** If the compact Cox monomial
> has local exceptional/boundary order \((N-r,m-1)\), it divides \(P_X\)
> literally on the attachment chart, with quotient \(\Phi_W\).

The theorem is expressed in the normalized scale \(\lambda\), not merely
in the base boundary equation.  This is essential here:

| boundary | pullback equation | root scale |
|---|---|---|
| \(q\) | \(q=\lambda\) on each colored branch | unramified |
| \(r\) | \(r=\lambda^2\) up to exceptional factors | simple ramification |
| \(d\) | \(d=\lambda^4\) up to exceptional factors | index four |

Thus one theorem replaces three separate derivative-order ledgers.

## 2. Generic strict \(q\) chart

On \(s\ne1\), the linear incidence equation \(h_q=0\) has centre

\[
 A_q=\frac{2s^2t+2s^2+3s-4}{2(s-1)}.            \tag{2.1}
\]

Exact substitution gives

\[
 P(A_q)\in(q^2),\qquad P_X(A_q)\in(q)            \tag{2.2}
\]

in \(\mathbf Q[s,t,(s-1)^{-1}]\).  Hence
\(X=A_q+qY\) gives

\[
 P=q^2\Phi,\qquad P_X=q\Phi_Y.                  \tag{2.3}
\]

This proves the generic unramified colored-branch attachment used below.
It does not extend through \(s=1\), where the leading coefficient of
\(h_q\) vanishes.

## 3. The cusp \(E_4\)-to-\(r\) attachment

Use

\[
\begin{aligned}
 s&=11+50e^2z+180e^4z^2+e^5z^2,\\
 t&=-\frac12+e^2z,\\
 X&=-1+eA.
\end{aligned}                                    \tag{3.1}
\]

The exact controlled equations are

\[
 r=e^{10}z^4R,\qquad h_r=e^8H,\qquad
 P=e^5F,\qquad P_X=e^4F_A.                       \tag{3.2}
\]

The strict boundary meets \(E_4\) at

\[
 e=0,\qquad z=\frac1{2560},qquad R_z=2560.      \tag{3.3}
\]

At that corner the incidence residual is

\[
 H(0,1/2560,A)
 =-\frac{256A^2-16A-1}{2^{30}},                 \tag{3.4}
\]

which is separable and therefore gives two geometric attachment colours.
Exact Singular saturation by the locally invertible coefficient
\(\operatorname{lc}_A(H)\) proves

\[
 F,F_A\in (R,H):\operatorname{lc}_A(H)^\infty.  \tag{3.5}
\]

On either incidence branch, put \(R=\lambda^2\) and
\(A=A_0+\lambda Y\).  The theorem gives

\[
 P=e^5\lambda^2\Phi,\qquad
 P_X=e^4\lambda\Phi_Y.                           \tag{3.6}
\]

The local compact section is \(Z_r=e^4\lambda T_r\), exactly the derivative
order in (3.6).

## 4. The triple \(E_2\)-to-\(d\) attachment

Over \(\mathbf Q(i)\), put

\[
\begin{aligned}
 s&=2i+e^2z,\\
 t&=-\frac34+\frac i2+e,\\
 X&=1+i+eZ.
\end{aligned}                                    \tag{4.1}
\]

Then

\[
 d=e^2D,\qquad P=e^4F,\qquad P_X=e^3F_Z,         \tag{4.2}
\]

and the exact \(d\)-root centre is \(Z_c=ez/2\).  Polynomial division proves

\[
 \partial_Z^kF(Z_c)\in(D),\quad 0\leq k\leq3,   \tag{4.3}
\]

whereas \(\partial_Z^4F(Z_c)\notin(D)\).  The strict boundary is transverse
at \(e=z=0\), since \(D_z=4i\).  With \(D=\lambda^4\) and
\(Z=Z_c+\lambda Y\), (1.2) gives

\[
 P=e^4\lambda^4\Phi,\qquad
 P_X=e^3\lambda^3\Phi_Y.                         \tag{4.4}
\]

Here \(Z_d=e\lambda T_d\), so \(Z_d^3\) has exactly the derivative order
\(e^3\lambda^3\).

## 5. The \(q\)-\(r\) \(A\)-packet-to-\(q\) attachment

Work over

\[
 8\alpha^3+16\alpha^2+2\alpha-7=0               \tag{5.1}
\]

with

\[
\begin{aligned}
 s_0&=4\alpha^2-5,\\
 m&=-8\alpha^2-4\alpha+10,\\
 X_A&=-2-2\alpha,\\
 s&=s_0+me+ze^2,\quad t=\alpha+e,\quad X=X_A+eY.
\end{aligned}                                    \tag{5.2}
\]

There are exact polynomials \(Q,H,F\) with

\[
 q=e^2Q,\qquad h_q=eH,\qquad P=e^3F,\qquad
 P_X=e^2F_Y.                                      \tag{5.3}
\]

The strict \(q\)-corner has

\[
 z=-\frac8{25}(166\alpha^2+139\alpha-93).       \tag{5.4}
\]

The coefficient \(H_Y\) is nonzero there.  Exact pseudo-remainder and
Gröbner calculations over the cubic residue field prove

\[
 \operatorname{prem}_Y(F,H)\in(Q^2),\qquad
 \operatorname{prem}_Y(F_Y,H)\in(Q).             \tag{5.5}
\]

Since \(H_Y\) is a unit locally, this is precisely the weighted Taylor
criterion with \(\lambda=Q\) and \(m=2\).  Therefore

\[
 P=e^3\lambda^2\Phi,\qquad
 P_X=e^2\lambda\Phi_W.                           \tag{5.6}
\]

The compact sections are \(Z_q=e\lambda T_q\) and \(Z_r=eT_r\), whose
product has order \(e^2\lambda\), as required.

## 6. Remaining strict-\(r\) Cartier centres

The full family saturation is still open at three further strict-\(r\)
packet types, but their corner fibres are no longer unknown.

At triple \(E_1\), over \(b^2+i=0\), the strict \(r\)-corner is

\[
 z=-\frac1{20}-\frac i{10}.                      \tag{6.1}
\]

After the exact scalings \(r=e^4R\), \(h_r=e^5H\), \(P=e^{10}F\), and
\(P_X=e^7F_Y\), the residual incidence has centre \(Y=0\) and nonzero pivot

\[
 H_Y=b\left(-\frac{527}{32}-\frac{21}{2}i\right). \tag{6.2}
\]

The residual of \(F\) has an exact double root there.  Thus the two
conjugate triple colours pass the Cartier-centre and double-root fibre gates.

At the \(q\)-\(r\) tangency, the strict \(r\)-corner is

\[
 z=-\frac6{25}(52\alpha^2+58\alpha-21).          \tag{6.3}
\]

The \(A\)- and \(B\)-packet root centres are, respectively,

\[
 c_A=\frac{8\alpha^2+2\alpha-9}{5},\qquad
 c_B=\frac{36\alpha^2+24\alpha-23}{5}.           \tag{6.4}
\]

For the \(B\)-packet, the first transform of \(h_r\) shares the residual
strict-boundary equation.  The required saturated incidence transform is
obtained by subtracting

\[
 -\frac32(434\alpha^2+353\alpha-317)R            \tag{6.5}
\]

and dividing once more by the exceptional parameter.  Both resulting
incidences have nonzero root-coordinate pivots, and \(F,F_Y\) have exact
double-root fibres at (6.4).  This certifies six further candidate colours,
three conjugates for each packet.

These eight residual certificates do not yet prove that \(F,F_Y\) belong
to the strict-\(r\) incidence ideal throughout a neighbourhood.  That
saturated family membership is the next algebraic gate.

## 7. Exact conductor-degenerate frontier

At the \(q\)-node use

\[
 s=1+e,\qquad t=-\frac12+ze,\qquad X=a+eY,
\qquad a^2-3a+1=0.                              \tag{7.1}
\]

Writing \(q=e^2Q\) and \(h_q=eH\), exact reduction gives

\[
 Q=4e^2z^2+8ez^2+4z^2+8z-1,
 \qquad H|_{e=0}=2a-2z-5.                        \tag{7.2}
\]

In particular, \(H|_{e=0}\) has no \(Y\)-term.  It selects a base slope,
not a root centre, so the Cartier-centre hypothesis of the theorem fails.
This is a precise conductor-normalization gate, not evidence against the
existence of the attachment.

The next exact targets are:

1. saturate the triple-\(E_1\)-to-strict-\(r\) incidence family from
   (6.1)--(6.2);
2. saturate the \(A\)- and \(B\)-packet strict-\(r\) families at the
   \(q\)-\(r\) tangency;
3. normalize (7.2) together with the conductor equalizer and compute its
   Cox-frame transition;
4. compute the conductor transitions at the triple and transverse packets.

Only after those edges pass is there a global compact
degree-\((3,1,1)\) Čech class to test for inverse-adjugate polynomiality and
affine-space recognition.

## 8. Reproduction

Run

```bash
.venv/bin/python scripts/verify_f20_strict_boundary_attachments.py \
  --output artifacts/generated-results/f20_strict_boundary_attachments.json
```

Singular is required for the saturated cusp-incidence membership.  The
other attachment gates use exact SymPy arithmetic over \(\mathbf Q\),
\(\mathbf Q(i)\), and the cubic field (5.1).
