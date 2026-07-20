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

## General program

The result supports the finite-flat conjecture

\[
\operatorname{rank}_{k[S]}\mathfrak Z_k=2^k.
\]

For `k>=3`, the correct question is not whether the special fiber is a tensor
product of dual numbers—it need not be—but whether one can find a monic
Groebner basis with `2^k` standard monomials and describe its multiplication
as a flat degeneration of the subset algebra over the squarefree-root locus.

Run:

```bash
python scripts/classify_transfer_block_k2.py
```
