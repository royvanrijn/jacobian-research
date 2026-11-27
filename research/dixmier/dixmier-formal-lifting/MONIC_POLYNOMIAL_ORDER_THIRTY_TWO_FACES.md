# Necessary faces for a constant-leading Weyl member of order thirty-two

**Theorem MonicPolynomialOrderThirtyTwoFaces.** Let \(k\) have
characteristic zero. Suppose \(P,Q\in A_1(k)\) have nonzero scalar
commutator, generate a proper subalgebra, and \(P\) has differential
order \(32\) and constant nonzero leading coefficient. After extending
constants to an algebraic closure, polynomial changes of Weyl
generators and subtraction of powers of \(P\) from the mate give
a pair with orders \(2g,gb\), with odd \(b\ge3\). Its leading
differential coefficients are proportional to \(x^{2d},x^{bd}\).
The possibilities are
\[
                      (g,d,\rho)=(4,12,4),\ (6,10,2).
\]
In the Fourier coordinates \(X=D,\ Y=-x\) write \(z=XY^\rho\)
in the commutative symbol ring. The leading symbols in direction
\((\rho,-1)\) are proportional to \(R^2,R^b\), where
\[
 R=X L(z),\qquad
 L(z)\ \text{is proportional to}\quad
       (z-u)^{g-1}\quad\text{or}\quad(z-u)^{g-2}(z-2u),
                   \qquad u\ne0 .
\]
The same pair has Bernstein degrees \(32,16b\), and has leading
terms proportional to \(D^2,D^b\) for the weight that assigns
\(D\) weight \(1\) and \(x\) weight \(-1\).

These are necessary faces of a putative pair, with one parameter
up to scaling. The [successor analysis](MONIC_POLYNOMIAL_ORDER_THIRTY_TWO_EXCLUSION.md)
rules out all four faces.

## Saturating the leading-degree bound

Extend constants to an algebraic closure. A proper generated
subalgebra stays proper because its nonzero vector-space
quotient stays nonzero after extension of scalars. Apply the
[order-thirty-one bound and its first residual](POLYNOMIAL_WEYL_WEIGHTED_DESCENT_AND_ORDER_THIRTY_ONE.md#the-monic-bound-and-the-first-residual).
After mate descent and polynomial depression the Bernstein leading
symbols have two independent linear factors. For the order-\(32\)
member these occur with unequal multiplicities, since in the
primitive factor the multiplicities have gcd one and its degree
is \(4,8\), or \(16\). Use a linear Weyl change so the leading
symbol becomes a nonzero multiple of
\[
                      x^{2d_0}\xi^{2g_0},
               \qquad d_0+g_0=16,\qquad 0<g_0<d_0 .
\]
The new differential order \(i\) is at least \(2g_0\) and
less than \(32\), since the leading symbol contains no
\(\xi^{32}\) term. The differential leading coefficient
cannot be constant by the order-thirty-one theorem. The
[individual leading-degree bound](POLYNOMIAL_WEYL_WEIGHTED_DESCENT_AND_ORDER_THIRTY_ONE.md#nonconstant-leading-factors-and-bernstein-degree)
then gives \(i+\deg\operatorname{lc}P\ge32\).
Equality must hold, since the Bernstein degree is \(32\).
The unique Bernstein leading monomial therefore lies in the
differential leading term, and \(i=2g_0\).

Reduce the mate by subtracting powers of this fixed member
whenever its order is divisible by \(i\). If the remaining
order \(j\) divides \(i\), then \(j<i\) and the leading
commutator identity makes the leading coefficient of that
mate nonconstant. The same identity gives
\(\deg\operatorname{lc}P=(i/j)\deg\operatorname{lc}Q\).
The individual bound for that mate would give
\[
            i+\deg\operatorname{lc}P
             = (i/j)(j+\deg\operatorname{lc}Q)\ge64 ,
\]
contrary to the equality above. Thus the remaining orders
are nondividing. If \(g\) is their gcd, their leading
coefficients have common factor \(h\) with exponents \(i/g\)
and \(j/g\). The general nonconstant-leading bound gives
\(g+\deg h\ge16\). Since \(i/g\ge2\), the equality
\[
                     32=(i/g)(g+\deg h)
\]
forces \(i=2g\), \(g=g_0\), and \(\deg h=d_0\).
The leading coefficient proportional to \(x^{2d_0}\) makes
\(h\) proportional to \(x^{d_0}\). Write \(d=d_0\) and
\(j=gb\), where \(b\ge3\) is odd.

## The standard minimal pair and its first edge

Here \(d>g\). The infinity shear and degree calculation
in the nonconstant-leading theorem give a presentation with
\[
               \deg_x p_i\le 2(d-g)+i,\qquad
               \deg_x q_j\le b(d-g)+j
\]
for the differential coefficients. These inequalities are
the coefficient bounds there with \(R=-s\). The leading
differential coefficients remain pure powers of \(x\).
In Fourier coordinates \(X=D,\ Y=-x\), with \([Y,X]=1\),
the supports lie in the rectangles with corners
\((2g,2d),(bg,bd)\) and contain those corners. Reordering
\(Y^jX^i\) only introduces terms with both indices smaller.
The only point with \(X\)-exponent \(2g\) in the first support
is \((2g,2d)\), and similarly for the mate. Hence
\(v_{1,-1}(\operatorname{st}_{1,0})<0\) for both.
Their total-degree gcd is \(16\). By
[Guccione–Guccione–Valqui, Corollary 7.4](sources/guccione-guccione-valqui-1111.6100v3.pdf),
this attains the minimum possible gcd for a proper pair.
After a scalar rescaling to \([Q,P]=1\), this is a standard
minimal pair in their Definition 6.10, printed page 36.

Propositions 7.2 and 7.3 of that source, printed pages 38–42,
give a direction \((\rho,\sigma)\) with \(\rho>-\sigma>0\)
and normalized endpoints \(C_0=(u_0,v_0),C_1=(r_0,s_0)\).
The two normalized endpoints of the first leading face are
\(C_0=\operatorname{en}_{\rho,\sigma}(P)/2\) and
\(C_1=\operatorname{st}_{\rho,\sigma}(P)/2\).
In particular \(u_0<v_0\), \(u_0\ge3\),
\(\gcd(u_0,v_0)>1\), \(r_0>s_0\ge0\), and \(r_0<u_0\).
Both endpoints have the same \((\rho,\sigma)\) weight.
Their Proposition 7.3 also gives an integer vector
\((f_1,f_2)=\tau C_0\), where \(f_1\ge2,\ 0<\tau<1\), and
\((\rho,\sigma)\) is the primitive direction of
\((f_2-1,1-f_1)\).

The numerical check in the proof of Corollary 7.4 rules
out every tuple satisfying these constraints and Proposition
7.3(6) with \(u_0+v_0\le15\). The support rectangle gives
\(u_0\le g,\ v_0\le d,\ u_0+v_0\le16\).
Thus \(u_0+v_0=16\) and \(C_0=(g,d)\). Since \(3\le g<8\)
and \(\gcd(g,16)>1\), the only corner choices are \((4,12)\)
and \((6,10)\).

For \((4,12)\), the vector \((f_1,f_2)\) is \((2,6)\) or
\((3,9)\). The first would give direction \((5,-1)\).
There is no integer \(0\le s_0<r_0<4\) with
\(5r_0-s_0=8\). The second gives \((4,-1)\), and
\(4r_0-s_0=4\) gives only \(C_1=(1,0)\).
For \((6,10)\), the only vector is \((3,5)\), giving
\((2,-1)\); \(2r_0-s_0=2\), with \(0\le s_0<r_0<6\),
again gives only \((1,0)\).
The leading face endpoints in either case are therefore
\(2(g,d),2(1,0)\) for \(P\) and \(b(g,d),b(1,0)\)
for \(Q\). Increasing the penalty on positive \(Y\)
exponents from \(\sigma=-1\) to \(\sigma=-\rho\), keeping
the \(X\) weight fixed, makes the left endpoint the sole
leading term. Every term involving \(Y\) loses positive
weight relative to that endpoint, while a term of
\(Y\)-degree zero already had \(X\)-degree at most the
endpoint's \(X\)-degree. After rescaling the direction
this is the asserted \((1,-1)\) statement.

## The polynomials on the two edges

The common face factor in the proof of Proposition 7.3,
equation (7.13), has endpoints \(C_0=(g,d)\) and
\(C_1=(1,0)\). Since \(d=\rho(g-1)\), it is
\(R=X L(XY^\rho)\) with \(\deg L=g-1,\ L(0)\ne0\).
The homogeneous element \(F\) preceding Proposition 7.3
has weight \(\rho-1\) and \([P,F]_{\rho,-1}=\ell_{\rho,-1}(P)\).
The leading commutator is the corresponding Poisson bracket,
with a choice of overall sign. Dividing the bracket identity
by a nonzero scalar times \(2R\) gives
\(\{R,F\}=cR\) with \(c\ne0\), writing \(F\) also for its
symbol. A polynomial of weight \(\rho-1\) in this direction
has the form \(F=XYG(z)\). Indeed a monomial \(X^iY^j\)
of that weight has \(j=\rho(i-1)+1\) and \(i\ge1\).
Using \(\{X,Y\}=1\) and absorbing a sign into \(c\) yields
\[
     \{X L(z),XYG(z)\}
       = X\bigl(LG+\rho z LG'-(\rho-1)z L'G\bigr).
\]
Consequently
\[
          \rho zG'+
          \left(1-(\rho-1)z\frac{L'}{L}\right)G=c.         \tag{1}
\]
In particular \(G(0)=c\). At any root of \(L\), \(G\) must
vanish to cancel the pole. It must vanish simply for the
left side to have the nonzero value \(c\). A root of \(L\)
of multiplicity \(m_0\) is therefore impossible if
\(\rho-(\rho-1)m_0=0\). Since \(\deg L>0\), \(G\) is
nonconstant, and comparison at infinity gives
\[
               \deg G=\frac{(\rho-1)(g-1)-1}{\rho}=2 .
\]
Thus \(L\) has at most two distinct roots, all nonzero.
When \(\rho=2\) a root of multiplicity two is impossible.

If \(L\) has one root \(u\), write \(L=L_0(z-u)^{g-1}\)
and \(G=\gamma(z-u)(z-v)\). In (1) the linear term
after cancellation is proportional to
\(-(\rho+1)u+\rho v\). Hence \(v=(\rho+1)u/\rho\).
If \(L\) has two distinct roots, their multiplicities
are \(g-2,1\) in both corner cases. For \(g=4\) this
is the only partition of three into two positive parts;
for \(g=6\) the partition \(3,2\) is excluded.
Write \(L=L_0(z-u)^{g-2}(z-v)\) and
\(G=\gamma(z-u)(z-v)\). The linear term is proportional
to
\[
              -2u+\bigl((\rho-1)(g-2)-(\rho+1)\bigr)v
                          = -2u+v .
\]
Therefore \(v=2u\). This gives the four stated face profiles.

The local normalization at \(x=0\) is consistent with
these exact directions. Put \(x=t^{-e}\), where
\(e=g/\gcd(g,d)\), and take \(E\) proportional to
\(x^{d/g}D\). Then \([E,t]\) is proportional to
\(t^{-1}\) in both cases, since \(ed/g-e=2\).
With \(\operatorname{wt}t=\operatorname{wt}E=1\) the
original derivative has weight \(1+ed/g=e\rho\),
and \(x\) has weight \(-e\). This gives precisely the
edge direction above. Writing \(\eta=\xi_E/t\), the
common symbol \(R\) becomes a scalar multiple of
\(t^g\eta L(c_1\eta)\) with \(c_1\ne0\).
Its \(\eta=0\) root is simple, so it is primitive.
Depression in \(E\) only translates \(\eta\), preserving
the multiplicities. The normalized \(r=1\) is attained:
the translated common polynomial has at least two
distinct roots and hence has a nonzero lower term.

## Primary comparisons and remaining equations

The closing discussion in the introduction of
[Guccione–Guccione–Valqui, *The Dixmier conjecture and the
shape of possible counterexamples II*, arXiv:1205.6827v1](https://arxiv.org/pdf/1205.6827)
already identifies \((4,12)\) and \((6,10)\) as the two
corners to analyze for \(B=16\). Its reduction of other
corners does not eliminate either of these.
[Han–Tan, Corollary 3.10](https://arxiv.org/pdf/2210.00257)
requires the centralizer of the lowest homogeneous component
to be generated by that component. Here the lowest component
is proportional to \(D^2\), whose centralizer contains \(D\).
The homogeneous-component bound of
[Guccione–Guccione–Valqui, arXiv:2402.11135](https://arxiv.org/abs/2402.11135)
excludes pairs with at most four components in one member.
The residual above gives no such upper bound.
