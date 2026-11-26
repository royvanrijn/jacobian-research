# Unconditional factorial nonvanishing for quadratic phases

Let \(\Lambda(X^\alpha)=\prod_i\alpha_i!\) over a field of
characteristic zero.

**Theorem.** Every nonzero polynomial \(F\) of total degree at most
two in \(m\) variables has infinitely many nonzero factorial moments.
More precisely, put \(R=2^m\). Every window at positive orders of
length
\[
                     W_m=R+\frac{(R+1)R(R+3)}2
\]
contains a nonzero factorial moment, uniformly over the coefficients.
In every residue class of moment orders, the moments are either
identically zero on that class or nonzero infinitely often.
In particular every positive power of a nonzero quadratic phase
has infinitely many nonzero factorial moments.

The calculation uses the exact exponential moment series on all
coordinate faces. It eliminates some singular faces algebraically
and gives differential equations for the others, including ones whose
quadratic parts are singular. At infinity the only zero residue comes
from constant phases. Its contribution to any nonconstant phase is
of strictly positive order in the reciprocal parameter. This excludes
every nonzero polynomial solution of a scalar equation for the series
and gives the window bound.

## The initial coefficient-moment ideal

**Corollary.** Fix \(m\ge1\) and a characteristic-zero field \(K\).
Put \(A=K[c_\alpha:|\alpha|\le2]\) and
\(F=\sum_{|\alpha|\le2}c_\alpha X^\alpha\). The ideal
\(I=(\Lambda(F),\ldots,\Lambda(F^{W_m}))\) has radical
\((c_\alpha:|\alpha|\le2)\). In particular \(A/I\) is
finite-dimensional over \(K\).

Indeed, if \(\mathfrak p\) is a prime of \(A\) containing \(I\), then
the coefficient images in \(\operatorname{Frac}(A/\mathfrak p)\)
give a quadratic phase with its first \(W_m\) moments zero. The
theorem in that characteristic-zero field forces every coefficient
image to vanish. Thus \(\mathfrak p\) equals the coefficient maximal
ideal. Conversely the generators of \(I\) vanish at the origin.
Intersecting the primes over \(I\) gives the radical equality. Each
coefficient has a power in \(I\), so finitely many monomials span
the quotient.

## Exact coordinate-face equations

Work over \(\mathbb C\) and write
\[
                F=c+\ell^{\mathsf T}X+X^{\mathsf T}QX,\qquad
                Q=Q^{\mathsf T}.
\]
For \(I\subseteq\{1,\ldots,m\}\) let \(F_I\) be the result of
setting the variables outside \(I\) to zero. Write \(Q_I\) for
the principal submatrix and \(\ell_I\) for the restricted column.
Let \(\mathbf 1_I\) be the all-ones column on \(I\), and set
\[
 \Phi_I(s)=\Lambda_I(e^{sF_I})
           =\sum_{n\ge0}\frac{\Lambda_I(F_I^n)}{n!}s^n,\qquad
 V_I(s)=\Lambda_I(X_Ie^{sF_I}),\qquad
 U_I(s)=(\Phi_{I\setminus\{i\}}(s))_{i\in I}.
\]
These are formal power series. The empty face has
\(\Phi_{\varnothing}=e^{cs}\), as does every face on which
\(F_I=c\).

For polynomials \(h\) in the variables of \(I\), integration by
parts is the factorial identity
\[
      \Lambda_I((\partial_i-1)h)
                  =-\Lambda_{I\setminus\{i\}}(h|_{X_i=0}).
 \tag{1}
\]
For a monomial \(h\) this follows from
\(k(k-1)!=k!\) when its \(X_i\)-exponent \(k\) is positive.
If \(k=0\) the left side is minus the factorial value on the
remaining coordinates. Thus (1) acts coefficientwise on all
the exponential series. Applied to \(e^{sF_I}\) it gives
\[
                2sQ_I V_I+(s\ell_I-\mathbf 1_I)\Phi_I=-U_I.
 \tag{2}
\]
Applied to \(X_i e^{sF_I}\), the right side of (1) is zero.
Summing over \(i\in I\) gives
\[
 2s\Lambda_I(X_I^{\mathsf T}Q_I X_I e^{sF_I})
       +(s\ell_I-\mathbf 1_I)^{\mathsf T}V_I
       =-|I|\Phi_I.
\]
Substitute this into \(\Phi_I'=\Lambda_I(F_I e^{sF_I})\) to get
\[
 \Phi_I'
     =\left(c-\frac{|I|}{2s}\right)\Phi_I
                   +\tfrac12(\ell_I+\mathbf 1_I/s)^{\mathsf T}V_I.
 \tag{3}
\]
Only formal power or Laurent series occur in these identities.

The face equation (2) is the factorial version of the integration
by parts in Koyama and Takemura's Lemma 8. Their Theorem 9 obtains a
Pfaffian system by inverting the principal submatrices in the
convergent negative-definite setting. Here the equations must also
handle singular principal submatrices and arbitrary complex
coefficients. The formal calculation and the eliminations below
supply these cases.
[Koyama–Takemura, Lemma 8 and Theorem 9](https://arxiv.org/pdf/1211.6822)

## Eliminations on singular faces

Keep one coordinate \(\Phi_{\varnothing}=e^{cs}\) for all constant
faces. For a nonconstant face \(I\), use the following three cases.

1. If some \(k\in\ker Q_I\) has \(k^{\mathsf T}\ell_I\ne0\), use (2)
   to write
   \[
           \Phi_I=-\frac{k^{\mathsf T}U_I}
                    {s\,k^{\mathsf T}\ell_I-k^{\mathsf T}\mathbf 1_I}.
   \tag{4}
   \]
   Each coefficient on a proper face is \(O(1/s)\) at infinity.

2. If \(\ell_I\) annihilates the kernel and some \(k\in\ker Q_I\)
   has \(k^{\mathsf T}\mathbf 1_I\ne0\), use
   \[
                     \Phi_I=\frac{k^{\mathsf T}U_I}
                                         {k^{\mathsf T}\mathbf 1_I}.
   \tag{5}
   \]
   A child face \(J=I\setminus\{i\}\) with \(F_J=c\) has
   \(Q_{J,J}=0\) and \(\ell_J=0\). If any \(Q_{ji}\ne0\) for
   \(j\in J\), the \(j\)-th kernel equation gives \(k_i=0\).
   If all those entries vanish and \(Q_{ii}\ne0\), the \(i\)-th
   equation gives \(k_i=0\). In the remaining case \(Q_I=0\)
   and the nonconstant face has \(\ell_i\ne0\), so
   \(k^{\mathsf T}\ell_I=0\) gives \(k_i=0\). Thus the constant
   coefficients in (5) are supported on nonconstant children.

3. In the remaining case both \(\ell_I\) and \(\mathbf 1_I\)
   annihilate the kernel. Since \(Q_I\) is symmetric for the
   transpose pairing, both belong to its image. Choose columns
   \(u_I,v_I\) solving
   \(Q_Iu_I=\ell_I,\ Q_Iv_I=\mathbf 1_I\).
   Multiplying (2) by \((u_I+v_I/s)^{\mathsf T}\) and substituting
   into (3) gives
   \[
   \begin{split}
    \Phi_I'&=
      \left(a_I-\frac{|I|}{2s}
             +\frac{\mathbf 1_I^{\mathsf T}v_I}{4s^2}\right)\Phi_I
          -\frac1{4s}(u_I+v_I/s)^{\mathsf T}U_I,\\
    a_I&=c-\frac{\ell_I^{\mathsf T}u_I}{4}.
   \end{split} \tag{6}
   \]
   The mixed terms in the diagonal coefficient cancel because
   \(u_I^{\mathsf T}\mathbf 1_I=v_I^{\mathsf T}\ell_I\).
   The formula uses any such pair of solutions.

There is one additional cancellation in (6). Suppose \(a_I=c\)
and \(F_{I\setminus\{i\}}=c\) for some \(i\). Then
\[
                                (u_I)_i=0. \tag{7}
\]
Indeed, outside the \(i\)-th row and column all entries of \(Q_I\)
are zero, and only the \(i\)-th entry of \(\ell_I\) can be nonzero.
A nonzero off-diagonal entry \(Q_{ji}\) gives
\(Q_{ji}(u_I)_i=(\ell_I)_j=0\).
If those off-diagonal entries all vanish and \(Q_{ii}\ne0\), the
equality \(a_I=c\) gives \(\ell_i^2/Q_{ii}=0\), so
\((u_I)_i=\ell_i/Q_{ii}=0\).
If \(Q_I=0\), \(\ell_I\in\operatorname{im}Q_I=0\) would make
the face constant. This proves (7) for all faces of the third type.

Starting with the smallest faces, express every eliminated face
by substituting (4) and (5). Keep \(Y_0=e^{cs}\) and a coordinate
\(Y_I=\Phi_I\) for each nonconstant face of type 3.
In particular \(\Phi_{[m]}=b(s)Y\) for a rational row \(b(s)\).
All face rows in this basis are bounded at infinity. For each
nonconstant face their coefficient on \(Y_0\) is \(O(1/s)\).
This holds on type 3 rows by the choice of coordinates, on
type 1 rows by (4), and on type 2 rows by induction because
(5) has zero coefficients on constant child faces.

The equation for \(Y_0\) and the equations (6) with the substituted
children give a triangular rational system. In increasing order
of face size its expansion at infinity is
\[
                Y'=A(s)Y,\qquad A(s)=A_0+s^{-1}A_1+O(s^{-2}).
 \tag{8}
\]
Here \(A_0=\operatorname{diag}(c,(a_I)_I)\), and \(A_1\) is
lower triangular with diagonal entries \(0\) on \(Y_0\) and
\(-|I|/2\) on \(Y_I\). In any row with \(a_I=c\) the
\(Y_0\)-column entry of \(A_1\) is zero. To see this, an
eliminated nonconstant child has \(Y_0\)-coefficient \(O(1/s)\),
which multiplies \(O(1/s)\) in (6). A constant child
has its \(1/s\) coefficient in (6) equal to zero by (7).

## Formal outputs at infinity

Put \(z=1/s\). A formal change \(Y=G(z)Z\) with
\(G=I+O(z)\) makes (8) block diagonal by the equal values of
the diagonal entries \(\alpha\) of \(A_0\). To see this, seek
the block diagonal form as \(G^{-1}AG-G^{-1}G'\).
At order \(z^k\), the new off-block coefficient of \(G\) enters
only through its commutator with \(A_0\), since
\(\partial_s=-z^2\partial_z\) raises \(z\)-order.
For entries joining distinct blocks this commutator is
multiplication by the nonzero difference of their \(\alpha\)'s.
Thus each order can be set to zero recursively. Within each block
the residue \(R_\alpha\), the coefficient of \(z\), is the
restriction of \(A_1\).

A block with constant part \(\alpha\), after removing that
constant part from its \(\partial_s\)-equation, has
\(z\partial_z\)-matrix \(C(z)\) with \(C(0)=-R_\alpha\).
It has a full formal solution basis with entries in
\[
                     E=\mathbb C((z^{1/2}))(\log z).
\]
Here \(\log z\) is a formal variable with \(z\partial_z\log z=1\).
For completeness, for each generalized eigenvector \(v\) of
\(C(0)=C_0\) with eigenvalue \(\lambda\), take
\[
 z^\lambda\sum_{k\ge0}z^k p_k(L),\qquad
 L=\log z,\qquad p_0(L)=e^{(C_0-\lambda I)L}v.
\]
The last expression is polynomial in \(L\). Writing
\(C(z)=\sum_{j\ge0}C_j z^j\), successive equations are
\[
  \bigl((\lambda+k)I-C_0+\partial_L\bigr)p_k
                       =\sum_{j=1}^k C_j p_{k-j}\quad(k\ge1).
 \tag{9}
\]
On a generalized eigenspace of \(C_0\) for eigenvalue
\(\mu\ne\lambda+k\), the operator is a nonzero scalar plus a
nilpotent operator on any finite space of polynomials of bounded
degree. Its finite inverse preserves the degree bound.
If \(\mu=\lambda+k\), it is \(\partial_L-N\) for a nilpotent
matrix \(N\). Multiplication by \(e^{-NL}\), polynomial
antidifferentiation, and multiplication by \(e^{NL}\) give a
polynomial solution. Only finitely many \(k\) are resonant,
so the degrees of the \(p_k\)'s are bounded. The
\(\lambda\)'s in (8) are among \(0\) and \(|I|/2\). The
resulting columns therefore lie in \(E\). Their determinant is
nonzero after factoring \(z^\lambda\) from each column, because
the leading columns in \(L\) are the original generalized
eigenbasis acted on by invertible block matrices \(e^{NL}\).

Combine these columns across blocks and multiply by \(G(z)\).
This gives an invertible matrix of columns \(w_j\) over \(E\)
with \(w_j'=A(s)w_j-\alpha_j w_j\), each associated to its own
constant part \(\alpha_j\). If \(\alpha_j=0\) and its initial
exponent \(\lambda\) is positive, \(b(s)w_j\) has positive
\(z\)-order because \(b\) is bounded.
For the block with \(\alpha_j=0\) to have exponent \(\lambda=0\),
it must contain \(Y_0\), so \(c=0\).
Its residue is lower triangular, has only one zero diagonal entry,
and its \(Y_0\)-column is zero by (7) and (8).
Hence the zero generalized eigenspace is exactly the line spanned
by the coordinate vector \(e_0\) for \(Y_0\).
The corresponding column starts with \(e_0\) and
\(b(s)e_0=O(z)\) for a nonconstant full phase.
Since \(G=I+O(z)\), \(b(s)w_j\) again has positive \(z\)-order.
This proves the required positive order for every zero-constant
block column.

## Excluding polynomial solutions and bounding the window

For a nonconstant \(F\) let \(b_0=b\) and \(b_{j+1}=b_j'+b_jA\).
Since \(bY=\Phi_{[m]}\) has constant coefficient one,
\(b\ne0\). Let \(r\) be the first dependent row index among
\(b_0,b_1,\ldots\) over \(\mathbb C(s)\).
Here \(1\le r\le\dim Y\le2^m\).
The relation expresses \(b_r\) in the first \(r\) rows and
gives a scalar equation for \(\Phi_{[m]}\), with an equivalent
first-order companion system on its \(r\) initial derivatives.
Denote the companion matrix by \(T(s)\).

For each of the formal columns \(w_j\), form
\(h_j=(b_iw_j)_{0\le i<r}\). Using the derivative recursion
and the relation for \(b_r\) gives \(h_j'=T h_j-\alpha_j h_j\).
The matrix of all the \(h_j\)'s has rank \(r\) over \(E\),
because the first \(r\) rows \(b_i\) are independent and the
matrix of \(w_j\)'s is invertible. Keep \(r\) independent columns
as a matrix \(H\), with diagonal \(\mathcal D\) of their
\(\alpha_j\)'s. Then \(H'=TH-H\mathcal D\).

Suppose a polynomial \(P(s)\) solves the scalar equation.
The column \(J\) of its first \(r\) derivatives solves
\(J'=TJ\). Thus \(\gamma=H^{-1}J\in E^r\) satisfies
\(\gamma'=\mathcal D\gamma\).
Every element of \(E\) expands as a Laurent series in
\(z^{1/2}\) with rational functions of \(L\) as coefficients.
This follows by expanding the denominator at its first nonzero
coefficient in \(z^{1/2}\), over \(\mathbb C(L)\).
For a nonzero \(\gamma_j\), \(\partial_s=-z^2\partial_z\)
raises its order by at least one, so \(\gamma'_j=\alpha_j\gamma_j\)
forces \(\alpha_j=0\).
Also the constants of \(E\) are exactly \(\mathbb C\).
Indeed a Laurent coefficient \(q(L)\) of \(z^\nu\) in a
constant satisfies \(q'(L)+\nu q(L)=0\).
A nonzero rational function has logarithmic derivative tending
to zero at \(L=\infty\), so this forces \(\nu=0\) and \(q\)
constant. It follows that \(P\), the first entry of \(H\gamma\),
is a constant linear combination of outputs \(b(s)w_j\) from
zero-constant blocks. All have positive \(z\)-order.
A nonzero polynomial in \(s=1/z\) has order at most zero.
This contradiction excludes every nonzero polynomial solution of
the scalar equation.

To bound its width, set \(B(s)\) equal to \(s^2\) times the
product of the linear denominators in (4), one factor for each
face of type 1. In the recursive substitutions each chain has
strictly descending face sizes, so each linear face denominator
occurs at most once on that chain. The product clears all the
denominators of every face row. Thus \(B A\) and \(B b\) are
polynomial, of degree at most \(q=\deg B\) because \(A,b\)
are bounded at infinity. We have \(q\le2^m+1\).
By induction \(b_j=L_j/B^{j+1}\) with a polynomial row \(L_j\)
of degree at most \((j+1)q\). In the induction the numerator
is \(B L_j'-(j+1)B' L_j+L_j B A\).

The signed minors of the \(r+1\) rows \(L_0,\ldots,L_r\)
in any \(r\) columns of full rank on the first \(r\) rows
give a relation \(\sum_{j=0}^r c_j L_j=0\) with \(c_r\ne0\).
Indeed the relation holds in the selected columns, and there
it determines uniquely the linear combination of the first
\(r\) rows equal to the last row. The same combination holds
on the full rows by their rank \(r\). Thus an annihilator is
\(\sum_{j=0}^r c_j B^j\partial_s^j\).
The minor for \(c_j\) omits the row of degree at most
\((j+1)q\), so
\[
               \deg(c_j B^j)\le
                        K=q\left(\frac{(r+1)(r+2)}2-1\right).
\]
Multiplication by \(s^r\) and
\(\partial_s^j=s^{-j}\theta(\theta-1)\cdots(\theta-j+1)\),
\(\theta=s\partial_s\), put the equation in the form
\(\mathscr L=\sum_{k=0}^{K+r}s^k P_k(\theta)\).
The transformations preserve its polynomial solution space.

If a block of \(K+r\) zero coefficients started at order \(n\ge1\),
truncate \(\Phi_{[m]}\) strictly below \(n\) to a polynomial \(P\).
Then \(\mathscr L P\) has degree at most \(n+K+r-1\).
Since \(\Phi_{[m]}-P\) starts at degree \(n+K+r\)
and the Euler form of \(\mathscr L\) preserves this lower order
bound, \(\mathscr L P=-\mathscr L(\Phi_{[m]}-P)\) also
starts at degree at least \(n+K+r\).
Hence \(\mathscr L P=0\), contradicting its nonzero constant
coefficient \(P(0)=1\). The bounds \(r\le R\), \(q\le R+1\)
give \(K+r\le W_m\). Dividing the factorial moments by \(n!\)
preserves their zeros. The nonzero constant phases have every
moment nonzero directly.

## Powers and coefficient fields

The polynomial exclusion argument persists for a direct sum of
finitely many systems and the sum of their scalar outputs, as long
as each scalar output on the block with constant part zero has
positive reciprocal order. Apply it to the systems for
\(\Phi_P(\zeta s)\) on the full face, for the \(d\)-th roots of
unity \(\zeta\). The \(s\)-rescaling changes each \(\alpha_j\)
to \(\zeta\alpha_j\), and preserves the reciprocal order.
For \(0\le a<d\), their weighted average is
\[
 \frac1d\sum_{\zeta^d=1}\zeta^{-a}\Phi_P(\zeta s)
        =\sum_{n\ge0}\frac{\Lambda(P^{dn+a})}{(dn+a)!}s^{dn+a}.
\]
The scalar exclusion shows that this series is either identically
zero or has infinitely many nonzero coefficients. At \(a=0\)
its constant coefficient is one. This proves the assertion on
\(P^d\) and the stated alternative in every residue class.
For a nonzero constant \(P\), both assertions follow directly.

Over any characteristic-zero field, the subfield generated over
\(\mathbb Q\) by the finitely many phase coefficients embeds in
\(\mathbb C\). Ordinary factorial identities and nonvanishing
therefore descend to the given field. The quadratic window holds
at every start, with singular principal forms and complex
coefficients included. The window length here depends on \(m\)
and can exceed the actual monomial support size.
