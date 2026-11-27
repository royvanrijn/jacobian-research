# A generic conic-argument pencil for the rank-two double conic

**Theorem QuinticKellerRankTwoDoubleConicGenericPencil.** Let \(q\) be
a geometrically irreducible homogeneous ternary quadratic and
let \(z\) be a source linear coordinate. A degree-five
three-variable Keller map cannot have leading parts
\[
 P_5=z(q^2+\delta z^4),\qquad Q_5=z^3q+\nu z^5,\qquad R_4=q^2
\]
if
\[
                              \delta(\delta+5\nu^2/9)\ne0 .
\]
The theorem works over \(\mathbb C\) and hence over characteristic
zero by embedding the finitely generated coefficient field.

For the more general \(P_5=z(q^2+\lambda z^2q+\mu z^4)\) and the same
\(Q_5,R_4\), subtract \(\lambda Q\) from \(P\) to obtain the
displayed form with \(\delta=\mu-\lambda\nu\). The alternatives
\(\delta=0\) and \(\delta=-5\nu^2/9\) are outside this theorem,
as are pencil spans without a \(q^2\) term in either quintic part.
If the line in the quintic span with zero \(q^2\) coefficient has
no \(z^3q\) term, it is the line spanned by \(z^5\).
That branch is outside this theorem. The special span
\(\langle z^5,zq^2\rangle\) is covered by the
[quartic-argument double-conic exclusion](QUINTIC_KELLER_RANK_TWO_DOUBLE_CONIC_QUARTIC_PENCIL.md).

## The mixed homogeneous equation

Write \(r=q(x,y,1)\), an irreducible binary polynomial of degree two.
The polynomial constants of \(\{r,\cdot\}_{x,y}\) are \(\mathbb C[r]\)
by the closed-polynomial criterion in Wang and Yan,
[Proposition 2.2 and Lemma 2.3](https://arxiv.org/pdf/2610.36495v1).
Denote the dehomogenizations of \(P_j,Q_j,R_j\) at \(z=1\)
by \(p_j,q_j,t_j\).

The degree-ten equation is
\[
  2r\big((3r+5\nu)p_4-(r^2+5\delta)q_4\big)
          +5(-r^2-2\nu r+\delta)t_3\in\mathbb C[r].
 \tag{1}
\]
To derive it, use rational source coordinates \((z,q,x)\) on the
chart \(q_y\ne0\). The determinant of this change of coordinates
is \(-q_y\), and
\[
          -q_y\,\partial_x|_{z,q}=\{q,\cdot\}_{x,y}.
\]
The top three rows in these coordinates depend only on \(z,q\).
The degree-ten contributions use the normal derivative of
\(P_4,Q_4\), or \(R_3\), with the minor of the other two leading
rows in \(z,q\). Set \(\lambda_0=q/z^2\) on \(z\ne0\).
The three minors, with cofactor signs, are
\[
 \begin{aligned}
   &2z^6\lambda_0(3\lambda_0+5\nu),\\
   &-2z^6\lambda_0(\lambda_0^2+5\delta),\\
   &5z^7(-\lambda_0^2-2\nu\lambda_0+\delta).
 \end{aligned}
\]
They depend only on \(q,z\), so summing the three contributions
and setting \(z=1\) gives the bracket of the left side of (1)
with \(r\), up to sign. Its vanishing gives (1).

## The kernel and its removal

Choose binary source coordinates so that \(q\) has nonzero
\(x^2\) coefficient and normalize that coefficient to one
by rescaling \(x\). Monic division in \(x\) expresses every
polynomial uniquely in the free \(\mathbb C[q,z]\)-basis
\(x^i y^j\), \(i\in\{0,1\}, j\ge0\). In this assertion the
monic equation is \(x^2+\text{lower terms in }x-q=0\).
The expansion preserves homogeneous degrees when \(q,z\)
have degrees two and one. The same basis is free over
\(\mathbb C[r]\) after dehomogenization.

For a basis monomial of positive degree \(d=i+j\), the respective
coefficient degrees in \(r\) in (1) are at most
\(\lfloor(4-d)/2\rfloor,\lfloor(4-d)/2\rfloor,\lfloor(3-d)/2\rfloor\),
where a negative bound means zero. Their contribution to (1)
must vanish. For \(d\ge2\) it suffices to use the largest
bounds, two affine polynomials and one constant. Write these
coefficients as \(a_1r+a_0,\ b_1r+b_0,\ c_0\).
The coefficients of \(r^4,r^3,r^0\) successively give
\[
                         b_1=0,\qquad b_0=3a_1,\qquad c_0=0 .
\]
Here the last equation uses \(\delta\ne0\).
The \(r^2\) coefficient then gives \(a_0=-5\nu a_1/3\),
and the \(r\) coefficient gives
\[
                           -30(\delta+5\nu^2/9)a_1=0 .
\]
Thus the kernel is zero on these bounds.

For \(d=1\) the third coefficient may also be affine.
The triple \((2r,1,2r)\) satisfies (1) identically.
Subtract a multiple to remove the \(r\) coefficient of the
third entry. The remaining triple is on the preceding bounds
and is zero. Thus the only nonconstant-basis terms in
\((P_4,Q_4,R_3)\) arise from a binary linear form \(L\) times
\((2zq,z^3,2q)\), the three derivatives of the leading forms
with respect to \(q\).

The homogeneous quadratic \(q\) is nondegenerate because its
projective conic is irreducible over \(\mathbb C\).
Choose a source translation for which the linear correction
to \(q\) is \(-L\). In the three leading forms its first-order
effect removes the terms just identified. The additional
effect through the translation of \(z\) already depends only
on \(q,z\). Hence \(P_4,Q_4,R_3\) are now polynomials in \(q,z\).

Proceed with \(j=3,2,1\), in that order. Assume all higher pieces
\(P_i,Q_i\) for \(i>j\) and \(R_i\) for \(i>j-1\) depend only
on \(q,z\). At degree \(j+6\) of the determinant, the terms
involving a normal derivative are precisely the same cofactor
terms as in (1), with \(p_j,q_j,t_{j-1}\) in place of
\(p_4,q_4,t_3\). Contributions using only higher pieces
vanish, since they factor through \(q,z\); contributions with
a lower piece or smaller cofactor degrees lie below \(j+6\).
Now all positive basis degrees have their coefficient bounds
within the zero-kernel bounds above. Thus the next three
pieces also depend only on \(q,z\). Induction makes all
components polynomials in \(q,z\), contradicting the
nonzero determinant.
