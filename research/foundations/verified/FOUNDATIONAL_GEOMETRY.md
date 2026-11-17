# Foundational Keller map and rational collision

Put `u=1+xy` and define

\[
F(x,y,z)=\left(
 u^3z+y^2u(4+3xy),
 y+3xu^2z+3xy^2(4+3xy),
 2x-3x^2y-x^3z
\right).
\]

## Theorem

Over every characteristic-zero field,

\[
\det DF=-2
\]

and

\[
F(0,0,-1/4)=F(1,-3/2,13/2)=F(-1,3/2,13/2)=(-1/4,0,0).
\]

Thus `F` is an everywhere-etale, noninjective polynomial map.  Appending
identity coordinates gives the same phenomenon in every dimension at least
three.

## Exact proof

On the dense open `x!=0`, set

\[
t=y+1/x,\qquad r=2/x,\qquad c=F_3(x,y,z),
\]

and write the first two target coordinates as `(a,b)`.  Direct substitution
gives

\[
a=t^2+\frac{rt}{2}-ct^3,
\qquad
b=r+4t-3ct^2.
\]

The two coordinate Jacobians are

\[
\det\frac{\partial(a,b,c)}{\partial(t,r,c)}=\frac r2,
\qquad
\det\frac{\partial(t,r,c)}{\partial(x,y,z)}=-2x.
\]

Since `r=2/x`, their product is `-2`.  Both sides are polynomials, so equality
on the dense open proves the determinant identity everywhere.  The collision
is exact substitution, and the three displayed source points are distinct.

## Role and assurance

This is the established Keller counterexample used as an explicit input to
the retained constructions. The determinant identity and rational collision
are proved directly above. The
[Archive of Formal Proofs entry](https://isa-afp.org/entries/Jacobian_Counterexample.html)
provides the separately authored formal reference.

The current applications include
[bounded forward lifts with unbounded inverses](../../dixmier/dixmier-formal-lifting/IDENTITY_LIFTS_AND_UNBOUNDED_INVERSES.md)
and the [integral cubic constructions](../../jacobian/cubic-reduction/PROOF12.md).
