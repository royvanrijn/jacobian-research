# Higher transfer blocks

For an allocation change

\[
(i,j)\longleftrightarrow(i-3k,j+2k)
\]

at one collision root, the stronger factorization correspondence is the
stronger factorization scheme

\[
\mathfrak Z_k={(U,V):U^2=V^3},
\]

where `U,V` are monic of degrees `3k,2k`, completed along

\[
U=S^3,\qquad V=S^2,qquad \deg S=k.
\]

Its reduced locus is the affine space of monic degree-`k` polynomials `S`.
At `S=Z^k`, the tangent equation

\[
2Z^{3k}\delta U-3Z^{4k}\delta V=0
\]

gives `2k` tangent directions, hence `k` excess infinitesimal directions over
the `k`-dimensional reduced locus.

## The two-transfer theorem

Let

\[
S=Z^2+pZ+q.
\]

Use the two coefficients transverse to the square locus to write

\[
V=S^2+XZ+Y.
\]

Then the coefficient equations of `U^2=V^3` eliminate all six nonleading
coefficients of `U` triangularly.  The remaining ideal is exactly

\[
\boxed{
X^3=0,\qquad 2XY=pX^2,\qquad Y^2=qX^2.}          \tag{1}
\]

Consequently

\[
\boxed{
\mathfrak Z_2=
\operatorname{Spec}
k[p,q,X,Y]/(X^3,2XY-pX^2,Y^2-qX^2)}             \tag{2}
\]

formally along the reduced `S`-plane.

The relations have monic leading terms `X^3,XY,Y^2`.  Polynomial division
therefore gives the free `k[p,q]`-basis

\[
1,\quad X,\quad Y,\quad X^2.
\]

Thus `Z_2` is finite flat of rank four over its reduced locus.  This proves
that the generic length `2^2`, obtained by separating the two roots of `S`,
does not jump when those roots collide.

## The coincident-root fiber

At `S=Z^2`, or `p=q=0`, the transverse algebra is

\[
\boxed{
k[X,Y]/(X^3,XY,Y^2).}                            \tag{3}
\]

It has basis `1,X,Y,X^2`, Hilbert function `(1,2,1)`, and nilpotency index
three.  Its socle is spanned by `X^2,Y`, so it has socle dimension two and is
not Gorenstein.

This is an important correction to the naive model

\[
k[\epsilon_1,\epsilon_2]/(\epsilon_1^2,\epsilon_2^2).
\]

The two algebras have the same length and Hilbert function, and they agree
after separating the two roots, but their multiplications differ at the
coincident-root point.  Higher transfer blocks are therefore nontrivial flat
degenerations of the squarefree subset algebra, not simply tensor products of
dual numbers at every point.

## Consequence for branch intersections

Normalized seed coefficients discard the constant and linear coefficients of
the collision polynomial.  The correct two-transfer equation is therefore

\[
U^2-V^3=\lambda Z+\mu.                            \tag{4}
\]

After the same triangular elimination, equation (4) retains only the
coefficient equations in degrees `5,4,3,2`; the degree-one and constant
equations are omitted.  Their exact Groebner basis is nevertheless still

\[
(X^3,2XY-pX^2,Y^2-qX^2).
\]

Moreover, the discarded degree-one and constant equations reduce to zero
modulo this ideal.  Thus affine difference does not enlarge the two-transfer
block:

\[
\boxed{\mathfrak Z_2^{\mathrm{aff}}=\mathfrak Z_2.} \tag{5}
\]

In particular, the actual isolated two-transfer correspondence has length
four:

\[
\operatorname{length}\mathfrak Z_2^{\mathrm{aff}}=4.
\]

For a global pair of allocations, transfers at different collision roots are
coupled by the fixed total `(a,b)` and by the single global affine difference.
The local result (5) supplies the exact two-transfer factor, but assembling
several compensating blocks still requires a Hensel-product/equalizer theorem.

## The all-`k` structure

The conjecture is now a theorem.  Over ordered roots `S=prod(Z-r_i)`, attach
commuting square-zero variables `epsilon_i` and set

\[
V=\prod_i((Z-r_i)^2+\epsilon_i),\qquad
U=\prod_i((Z-r_i)^3+	frac32\epsilon_i(Z-r_i)).
\]

The transfer block is the symmetric-group quotient of this Boolean
thickening.  Consequently it is finite flat of rank `2^k`, affine difference
adds nothing for every `k`, and the coincident-root Hilbert series is
`(1+t)^k`.  The complete invariant-theoretic proof is in
[ALL_K_TRANSFER_BLOCK_THEOREM.md](ALL_K_TRANSFER_BLOCK_THEOREM.md).  The
explicit results below are retained as coordinate models and independent
low-degree audits of the structural theorem.

## The three-transfer theorem

Let

\[
S=Z^3+pZ^2+qZ+r,\qquad V=S^2+XZ^2+YZ+T.
\]

The nine highest coefficient equations eliminate the nine nonleading
coefficients of `U`.  Over the coefficient ring `k[p,q,r]`, the remaining
equations have the monic relative Groebner basis

\[
\begin{aligned}
T^3={}&6T^2r^2+6X^2pr^3-12XYr^3,\\
T^2X={}&2T^2pr+2X^2p^2r^2-4XYpr^2,\\
TX^2={}&2T^2q+2X^2pqr-4XYqr,\\
X^2Y={}&4T^2p+4X^2p^2r-8XYpr,\\
X^3={}&6T^2+6X^2pr-12XYr,\\
2TY={}&-X^2pq+X^2r+2XYq,\\
Y^2={}&-2TX-X^2p^2+X^2q+2XYp.                 \tag{6}
\end{aligned}
\]

Its standard monomials are

\[
1,\quad T,\quad T^2,\quad X,\quad TX,\quad X^2,
\quad Y,\quad XY.
\]

Therefore

\[
\boxed{\operatorname{rank}_{k[p,q,r]}\mathfrak Z_3=8.} \tag{7}
\]

As for `Z_2`, retaining only the coefficient equations in degrees at least
two gives exactly the same ideal: the discarded linear and constant equations
reduce to zero modulo (6).  Thus

\[
\boxed{\mathfrak Z_3^{\mathrm{aff}}=\mathfrak Z_3.}    \tag{8}
\]

At `S=Z^3`, the special fiber is

\[
\begin{aligned}
T^3=T^2X=TX^2=X^2Y=TY=0,\\
X^3=6T^2,qquad Y^2=-2TX.                       \tag{9}
\end{aligned}
\]

It has length eight and Hilbert function `(1,3,3,1)`.  Both `T^2` and `XY`
lie in its socle, so the socle has dimension two and the algebra is not
Gorenstein.  Once again, collision preserves the generic length `2^3` but
changes the multiplication of the squarefree subset algebra.

## The four-transfer theorem

Let

\[
S=Z^4+pZ^3+qZ^2+rZ+t,
\qquad
V=S^2+AZ^3+BZ^2+CZ+D.
\]

The twelve highest coefficient equations eliminate the twelve nonleading
coefficients of the monic degree-twelve polynomial `U`.  The coefficient
equations in degrees eleven through two have a thirteen-element monic
relative Groebner basis over `k[p,q,r,t]`.  Its leading monomials are

\[
\begin{gathered}
A^2CD,D^3,CD^2,BD^2,B^2D,B^3,AD^2,ABD,AB^2,\\
A^2B,A^3,C^2,BC.
\end{gathered}                                                    \tag{10}
\]

Consequently its standard monomials are

\[
\begin{gathered}
1,A,B,C,D,A^2,AB,B^2,AC,AD,BD,CD,D^2,\\
A^2C,A^2D,ACD,
\end{gathered}                                                    \tag{11}
\]

and the quotient is free over the full monic-quartic parameter space:

\[
\boxed{\operatorname{rank}_{k[p,q,r,t]}\mathfrak Z_4=16.}        \tag{12}
\]

The discarded linear and constant coefficient equations reduce to zero
modulo the same basis.  Hence affine difference again introduces neither
extra components nor an embedded correction:

\[
\boxed{\mathfrak Z_4^{\mathrm{aff}}=\mathfrak Z_4.}              \tag{13}
\]

At the coincident quartic `S=Z^4`, a convenient specialized Groebner basis is

\[
\begin{gathered}
A^2CD, D^3, CD^2, BD^2, D(AC-2B^2), B^3-3A^2D,\\
AD^2, ABD, A(AC+B^2), A^2B-2D^2,\\
A^3-12CD, C^2+2BD, AD+BC.
\end{gathered}                                                    \tag{14}
\]

The powers of its maximal ideal have dimensions `15,11,5,1,0`, so the
Hilbert function is `(1,4,6,4,1)`.  Its socle has dimension four.  Thus the
length remains the generic value `2^4`, while the increasingly large socle
shows that the coincident fibers are not approaching a complete-intersection
or Gorenstein pattern.

Run:

```bash
python scripts/classify_transfer_block_k2.py
python scripts/classify_transfer_block_k3.py
python scripts/classify_transfer_block_k4.py
```
