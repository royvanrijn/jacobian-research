# Exceptional Cox corner charts for the corrected \(F_{20}\) cover

## Status and outcome

The one-parameter exceptional atlas has now been upgraded at every genuine
positive exceptional--exceptional intersection.  Seven exact two-parameter
corner types cover the ramphoid-cusp and \(q\)-\(r\)-tangency resolution
graphs.  On each corner the compact \((3,1,1)\) Cox monomial divides \(P_X\)
literally, with quotient the derivative of the bivariate strict transform.

The complementary blowup charts at the \(q\)-node also have an exact
transition law preserving \(P_X\).  At the triple orbit, exact residual-root
data proves that the positive \(E_1\) and \(E_2\) colored primes are
disjoint, even though their base exceptional rays are adjacent.  Thus no
transition map between those two positive colors should be imposed.

This closes the exceptional--exceptional corner gate represented by the
current colored resolution graph.  It does not yet construct the attachment
corners where strict \(d\)-, \(q\)-, or \(r\)-boundary components meet the
exceptional locus, nor the conductor-to-corner transitions.  Consequently a
global Čech class and global Cox algebra are still not defined.

The exact checker is
[`verify_f20_exceptional_cox_corners.py`](../scripts/verify_f20_exceptional_cox_corners.py),
with generated certificate
[`f20_exceptional_cox_corners.json`](../artifacts/generated-results/f20_exceptional_cox_corners.json).

## 1. General two-parameter Cox theorem

Let \(A\) be a regular corner chart with exceptional primes

\[
 E_u=(u=0),\qquad E_v=(v=0).
\]

Suppose

\[
 P=u^Nv^MF(u,v,z,Y),\qquad
 X=\xi+u^rv^sY.                                    \tag{1.1}
\]

Then

\[
 P_X=u^{N-r}v^{M-s}F_Y.                            \tag{1.2}
\]

If

\[
 D_j=a_jE_u+b_jE_v,
\]

its local multi-Rees generator is

\[
 Z_j=u^{a_j}v^{b_j}T_j.                            \tag{1.3}
\]

> **Exceptional Cox corner theorem.** A multidegree \(n=(n_j)\) has exact
> derivative order at the corner precisely when
> \[
>  \sum_jn_j(a_j,b_j)=(N-r,M-s).                   \tag{1.4}
> \]
> In that case its local Cox monomial divides \(P_X\) literally and the
> quotient is \(F_Y\).  Coordinatewise smaller orders leave an additional
> monomial factor; an order exceeding either coordinate cannot divide
> \(P_X\).

This is the two-dimensional version of the controlled-transform theorem.
It replaces two separate valuation ledgers by one exact corner identity.

## 2. The cusp corner graph

Use tangent coordinates

\[
 V=t+\frac12,\qquad W=(s-11)-50V.                  \tag{2.1}
\]

The embedded-resolution graph relevant to the colored root cover is

\[
 E_1\;--\;E_2\;--\;E_4\;--\;E_3.                  \tag{2.2}
\]

The divisor \(E_4\) is the blowup of the \(E_2\)-\(E_3\) corner; its base
valuation vector \((2,5)\) is the sum of \((1,2)\) and \((1,3)\).

The normalized source-corner substitutions are:

| corner | \(V\) | \(W\) | \(X+1\) |
|---|---|---|---|
| \(E_1E_2\) | \(u^5v^5\) | \(u^5v^{10}\) | \(uv^2Y\) |
| \(E_2E_4\) | \(u^5v^2\) | \(180V^2+u^{10}v^5\) | \(u^2vY\) |
| \(E_4E_{3,\mathrm{unr}}\) | \(u^2v\) | \(180V^2+u^5v^3\) | \(uvY\) |
| \(E_4E_{3,\mathrm{ram}}\) | \(u^2v^2\) | \(180V^2+u^5v^6\) | \(uvY\) |

Exact substitution produces:

| corner | \((N,M)\) | \((r,s)\) | \(D_r\) order | \(P_X\) order |
|---|---:|---:|---:|---:|
| \(E_1E_2\) | \((5,10)\) | \((1,2)\) | \((4,8)\) | \((4,8)\) |
| \(E_2E_4\) | \((10,5)\) | \((2,1)\) | \((8,4)\) | \((8,4)\) |
| \(E_4E_{3,\mathrm{unr}}\) | \((5,3)\) | \((1,1)\) | \((4,2)\) | \((4,2)\) |
| \(E_4E_{3,\mathrm{ram}}\) | \((5,5)\) | \((1,1)\) | \((4,4)\) | \((4,4)\) |

The corner residuals are, respectively,

\[
 Y^5+\frac{25}{2},\quad
 Y^5+\frac{25}{2},\quad
 \frac{25}{2}(40Y+1),\quad
 Y(Y^4-50Y^2+500),                                 \tag{2.3}
\]

and are separable.  Since only \(D_r\) occurs at these cusp colors, (1.4)
is immediate and exact on all four corners.

## 3. The \(q\)-\(r\) tangency corners

Work over

\[
 8\alpha^3+16\alpha^2+2\alpha-7=0
\]

and put

\[
 V=t-\alpha,\qquad
 U=(s-s_0)-mV,qquad
 m=-8\alpha^2-4\alpha+10.                          \tag{3.1}
\]

For the ramified \(A\)- and \(B\)-corners use

\[
 V=u^2v,\qquad U=u^2v^2,qquad X-X_{A/B}=uvY.       \tag{3.2}
\]

For the unramified \(A\)-corner use

\[
 V=uv,qquad U=uv^2,qquad X-X_A=uvY.              \tag{3.3}
\]

The exact order table is:

| corner | \((N,M)\) | \(D_q\) | \(D_r\) | derivative order |
|---|---:|---:|---:|---:|
| \(E_{1,A,\mathrm{ram}}E_{2,A}\) | \((3,3)\) | \((1,1)\) | \((1,1)\) | \((2,2)\) |
| \(E_{1,A,\mathrm{unr}}E_{2,A}\) | \((2,3)\) | \((1,1)\) | \((0,1)\) | \((1,2)\) |
| \(E_{1,B,\mathrm{ram}}E_{2,B}\) | \((2,2)\) | \((0,0)\) | \((1,1)\) | \((1,1)\) |

Thus \(Z_qZ_r\), \(Z_qZ_r\), and \(Z_r\), respectively, have the exact
two-variable order of \(P_X\).  The cubic, linear, and quadratic corner
residuals are separable over the cubic residue field, so all three strict
transforms are generically regular at the corner.

## 4. The \(q\)-node complementary transition

Let \(a^2-3a+1=0\) be either repeated root at the node.  The two ordinary
blowup charts are

\[
\begin{array}{lll}
 s-1=u,&t+\frac12=zu,&X-a=uY_L,\\
 t+\frac12=v,&s-1=\zeta v,&X-a=vY_R.
\end{array}                                        \tag{4.1}

On their overlap,

\[
 v=zu,qquad \zeta=z^{-1},\qquad Y_R=z^{-1}Y_L.    \tag{4.2}

If \(F_L,F_R\) are the strict transforms and \(G_L,G_R\) their derivative
quotients, exact calculation gives

\[
 F_R=z^{-2}F_L,qquad G_R=z^{-1}G_L.                \tag{4.3}

The \(D_q\) frames satisfy

\[
 Z_{q,R}=zZ_{q,L}.                                  \tag{4.4}

Equations (4.3)--(4.4) make \(P_X=Z_qG\) invariant on the overlap.  This is
the first literal Čech transition in the exceptional Cox construction.

## 5. Colored root-center separation

Base-fan adjacency does not imply adjacency of colored primes on the
normalized root cover.

> **Root-center separation lemma.** Let two candidate colored primes over
> adjacent base divisors have normalized residual root centers \(c_1,c_2\)
> in a common residue extension.  If \(c_1-c_2\) is a unit, their strict
> transforms are disjoint and no Cox transition edge joins them.

At triple \(E_1\),

\[
 \frac{X-X_0}{\tau^2}=b+O(\tau),\qquad b^2+i=0,    \tag{5.1}
\]

whereas the positive \(E_2\) cluster is centered at \(X-X_0=0\).  Since
\(b\) is a unit modulo \(b^2+i\), the two positive colored primes do not
meet.  The apparent transition requested by the base fan is therefore a
false edge and must be deleted from the colored overlap graph.

## 6. Remaining global attachments

The exceptional--exceptional part of the overlap graph is now exact.  The
remaining edges are of three kinds:

1. strict \(d\)-, \(q\)-, and \(r\)-boundary components meeting their final
   exceptional divisors;
2. generic ramification charts attaching to those strict transforms;
3. conductor frames attaching to the node, triple, and \(q\)-\(r\) corner
   charts.

Until those edges are computed, the compact multidegree-\((3,1,1)\) Čech
class is not defined globally.  Full inverse-adjugate polynomiality and
affine-space recognition therefore remain gated.

## 7. Reproduction

Run

```bash
.venv/bin/python scripts/verify_f20_exceptional_cox_corners.py \
  --output artifacts/generated-results/f20_exceptional_cox_corners.json
```

The checker verifies the seven bivariate strict-transform identities, their
separable corner residuals, the compact Cox orders, the \(q\)-node transition
law, and the triple root-center separation certificate.
