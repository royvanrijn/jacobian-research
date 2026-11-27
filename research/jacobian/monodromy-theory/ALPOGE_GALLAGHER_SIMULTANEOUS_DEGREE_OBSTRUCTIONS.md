# Simultaneous coordinate changes in two degree-seven families

**Theorem AlpogeGallagherFlagTriangularSourceDegreeObstruction.**
Let \(F=(A,B,C)\) be the retained Alpöge–Gallagher Keller map
over a characteristic-zero field. Let \(\Theta\) be any target
polynomial automorphism. For every source automorphism of the form
\[
 \Psi(x,y,z)=(a x+b,\ c y+p(x),\ d z+q(x,y)),
 \qquad acd\ne0 ,
\]
one has \(\deg(\Theta\circ F\circ\Psi)\ge7\).

**Theorem AlpogeGallagherIndependentSourceLeadingFormsObstruction.**
Let \(\Psi=(x_1,y_1,z_1)\) be any dominant polynomial source map.
Set \(g_1=2-3x_1y_1-x_1^2z_1\).
If the ordinary highest forms of \(x_1,y_1,g_1\) are
algebraically independent, then for any polynomial target
automorphism \(\Theta\) one has
\(\deg(\Theta\circ F\circ\Psi)\ge7\).

Both assertions use the same target basis as the
[separate-coordinate obstruction](ALPOGE_GALLAGHER_SEPARATE_DEGREE_OBSTRUCTIONS.md).
For the triangular source, the leading \(z\)-coefficient of
each pullback gives an ordinary degree bound even after a
polynomial shear of \(y\) or \(z\).

## A basis in the localized target ring

Use \(u=1+xy,\ g=2-3xy-x^2z,\ w=gu\).
Put \(X=BC,\ Y=AC^2\), and let
\[
            \Delta=4X^2-4X^3-64Y-108Y^2+72XY .
\]
The source identities in the cited separate-coordinate proof
give \(C=xg\) and
\[
            X=-3w^2+4w+2g,\qquad Y=-w^3+w^2+gw,\qquad
            \Delta=4g^2(9w^2-12w+4-8g).
\]
These follow also from
[Gao, Theorem 3.3 and the displayed expansions](sources/gao-2608.00222v2.pdf).

After inverting \(C\), use coordinates \(C,X,Y\) on the
target and group a polynomial \(V(A,B,C)\) by its power
of \(C\) as \(\sum_\lambda C^\lambda J_\lambda(X,Y)\).
For each monomial \(X^iY^j\) in \(J_\lambda\) one has
\(i+2j\ge-\lambda\), because it comes from
\(B^i A^j C^{\lambda+i+2j}\).
In particular \(J_{-2}\in(Y,X^2)\).
The coefficient of \(Y\) in \(J_{-2}\) is exactly that of
\(A\) in the target linear part.

The ring \(k[X,Y]\) is free over \(k[X,\Delta]\)
with basis \(1,Y\), by the monic quadratic relation for
\(Y\) after dividing by its leading scalar.
Expand all \(J_\lambda\) in this basis:
\[
                         C^\lambda X^iY^j\Delta^k,
                \qquad i,k\ge0,\qquad j\in\{0,1\}.
 \tag{1}
\]
Every expression below is a finite sum. Degrees in the
localized ring can be computed after multiplication by
a power of \(C\) and then subtracting the degree of that
power.

## Leading \(z\)-coefficients and flag-triangular shears

The leading \(z\)-coefficients of \(C,X,Y,\Delta\) in the
source are nonzero scalar multiples of
\[
                  x^3,\quad x^4u^2,\quad x^6u^3,\quad x^8u^2
\]
and their \(z\)-degrees are \(1,2,3,4\).
Thus a term of (1) has \(z\)-degree \(n\) and leading
\(z\)-coefficient proportional to \(x^a u^m\), where
\[
       n=\lambda+2i+3j+4k,\qquad
       a=2n+\lambda,\qquad m=n-\lambda-2k .
 \tag{2}
\]
For fixed \(n\), the exponent pair \((a,m)\) determines
\(\lambda,k,j,i\): first recover \(\lambda,k\) and then use
the parity of \(n-\lambda-4k\) to determine \(j\).
The functions \(x,u\) are algebraically independent.
Consequently the maximum \(n\) of a term present in (1)
is the actual \(z\)-degree of the pullback \(V\circ F\),
and its leading \(z\)-coefficient is the corresponding
nonzero sum of the displayed \(x^a u^m\).

Let \(m_0\) be the largest \(u\)-exponent in that sum.
The highest \(y\)-coefficient there is a sum of distinct
nonzero monomials \(x^{a+m_0}\), since \(u=1+xy\).
These exponents are nonnegative because the coefficient
is a polynomial in \(x\), although some individual
\(a\) in the localized expansion can be negative.
The substitution \(\Psi\) in the first theorem preserves
\(z\)-degree and replaces its highest coefficient by
a nonzero scalar times the same polynomial evaluated
at \(a x+b,c y+p(x)\).
It preserves the highest \(y\)-degree of that coefficient
and the \(x\)-degree of its highest \(y\)-coefficient.
Hence any term at the maximum \(n,m_0\) gives an ordinary
degree lower bound \(n+a+2m_0\) after the source change.
By (2),
\[
                            a+2m=3n+2i+3j .
\]
If \(n\ge2\), this gives degree at least \(4n\ge8\).

Suppose now that \(V\) has nonzero coefficient of \(A\)
in its target linear part and \(\deg(V\circ F\circ\Psi)\le6\).
Then \(n\le1\). In the \(\lambda=-2\) group, (2)
allows only \(1,X,Y\), since \(2i+3j+4k\le3\).
As \(J_{-2}\in(Y,X^2)\) and its \(Y\)-coefficient is
nonzero, that group is a nonzero scalar multiple of
\(A=C^{-2}Y\). In particular \(n=1\) and \(m_0\ge3\).

If \(m_0=3\), the \(A\) term gives the monomial \(x^3\)
in the top \(y\)-coefficient. This gives total degree at
least \(1+6=7\) after the source change.
If \(m_0\ge5\), (2) at \(n=1\) gives \(a+m_0=3-2k\)
for every highest-\(y\) term. Nonnegativity forces
\(k\le1\), so \(a+m_0\ge1\) and the total degree is at least
\(1+m_0+1\ge7\).

In the remaining case \(m_0=4\), a highest-\(y\) term with
\(k=0\) already forces \(1+4+3=8\).
The only term with \(n=1,\ m=4,\ k=1\) is
\(C^{-5}X\Delta\). This term cannot occur with \(n\le1\):
the \(\lambda=-5\) group must obey \(i+2j\ge5\) in its
ordinary \(X,Y\)-monomials, and its allowed basis terms
have \(2i+3j+4k\le6\). Among those terms only \(X\Delta\)
has an \(X^4\) contribution, with coefficient \(-4\);
the constraint forbids \(X^4\).
Thus this case is also impossible.

Every polynomial automorphism of the target has invertible
linear part at the target origin. At least one of its
components has a nonzero coefficient of \(A\) there, proving
the first theorem.

## Ordinary leading forms after a source change

For the second theorem write
\(d_x=\deg x_1,\ d_y=\deg y_1,\ d_g=\deg g_1,\)
\(d_C=\deg(C\circ\Psi)\), and \(d_w=\deg(w\circ\Psi)\).
Dominance implies \(3y_1+x_1z_1\) is nonconstant, so
\(d_g\ge d_x+1\), and
\(d_C=d_x+d_g,\ d_w=d_g+d_x+d_y\).
The leading forms of \(C\circ\Psi,g_1,w\circ\Psi\) are
algebraically independent under the stated hypothesis:
in their rational function field one recovers the
leading forms of \(x_1,y_1\) by dividing the leading
form of \(C\) by that of \(g\), and that of \(w\) by
that of \(C\), respectively.

The ordinary leading terms of the basis (1) are
proportional to independent Laurent monomials in those
three leading forms. Their powers are
\[
                        (\lambda,\ 2k,\ 2i+3j+2k)
\]
in \(C,g,w\), respectively. They determine the original
indices by parity as before. After clearing a power of
\(C\), algebraic independence thus gives the maximum
ordinary degree over the terms of that basis, with no
cancellation of the leading degree.

If \(V\) has \(A\)-coefficient nonzero and its pullback
after \(\Psi\) has ordinary degree at most six,
\(J_{-2}\) must have weighted degree at most \(6+2d_C\).
In its normalized basis expansion, \(Y,\Delta,X^2\)
have respective weights \(3d_w,\ 2d_g+2d_w,\ 4d_w\),
all strictly larger than \(6+2d_C\). Indeed
\[
 \begin{aligned}
  3d_w-2d_C &= d_g+d_x+3d_y\ge7,\\
  2d_g+2d_w-2d_C &=2d_g+2d_y>6,\\
  4d_w-2d_C &=2d_g+2d_x+4d_y>6.
 \end{aligned}
\]
For the first two inequalities, if \(d_y=1\), then
\(d_g=2d_x+\deg z_1\ge3\) because the \(x_1^2 z_1\)
term strictly dominates \(3x_1y_1\). If \(d_y\ge2\),
use \(d_g\ge d_x+1\). Thus the only basis terms
permitted in \(J_{-2}\) are \(1,X\), contrary to
\(J_{-2}\in(Y,X^2)\) with nonzero \(Y\)-coefficient.
This proves the second theorem.

Together the two theorems impose a necessary condition on a
simultaneous reduction of this map to ordinary degree at most
six. The source automorphism must lie outside the displayed
flag-triangular class, and the highest forms of
\(x_1,y_1,2-3x_1y_1-x_1^2z_1\) must be algebraically dependent.

The two source classes need not coincide. For example
\(\Psi=(x,y-xz/3,z)\) satisfies the second criterion:
the highest forms of \(x_1,y_1,g_1\) are
\(x,-xz/3,-3xy\).
