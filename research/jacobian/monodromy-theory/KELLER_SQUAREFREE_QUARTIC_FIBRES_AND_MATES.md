# Fibres and source-direction mates for the squarefree quartic scalar

Consider
\[
 a=x(xy-1),\qquad
 R=az+b,\qquad
 b=\lambda_0y+\mu_0+
       xy\bigl((\lambda_1-\lambda_0)y+\mu_1-\mu_0\bigr),
 \qquad \lambda_0\lambda_1\ne0 .
\]
This is the [squarefree-coefficient quartic form](QUARTIC_AFFINE_SOURCE_COMPONENTS_THROUGH_DEGREE_SEVEN.md)
left by the affine-source reduction. Write
\(\Delta=(\lambda_1-\lambda_0)y+\mu_1-\mu_0\), \(Y=1-xy\).
The identity \(R-\mu_1=\lambda_1 y-Y(xz+\Delta)\)
will be useful.

**Theorem KellerSquarefreeQuarticFibreModel.** Over \(\mathbb C\),
for \(t\ne\mu_1\) the fibre \(R=t\) is the Pełka–Raźny surface
\(S_{T+\theta(t),1}\), where
\[
                  \theta(t)=-\frac{\lambda_1(t-\mu_0)}
                                      {\lambda_0(t-\mu_1)}.
\]
These fibres are factorial, have constant units and have logarithmic
Kodaira dimension zero. They are pairwise nonisomorphic when
\(\mu_0\ne\mu_1\), and the family off \(\mu_1\) is algebraically
trivial when \(\mu_0=\mu_1\). The exceptional fibre \(R=\mu_1\)
is a principal open subset of the affine plane. Each ordinary
fibre contains exactly three images of nonconstant polynomial
maps from the affine line, with one additional line when its
parameter \(\theta\) is \(0\) or \(-1\).

**Theorem KellerSquarefreeQuarticNoPlanarMate.** Over a
characteristic-zero field, this scalar has no polynomial
Keller completion \((P,Q,R)\) with either mate independent
of \(z\). This obstruction has no degree bound. In particular,
the three polynomials \(P_z,Q_z,R_z\) are linearly independent
over the coefficient field.

**Theorem KellerSquarefreeQuarticNoSimultaneousAffineMates.**
The same squarefree scalar has no polynomial Keller completion
with both mates affine in \(z\), over any characteristic-zero
field and with no ordinary degree bound.

## The fibre model

Fix \(t\). Put \(y_0=(t-\mu_0)/\lambda_0\).
On the fibre, \(L=(y-y_0)/x\) is regular because the equation
\(R=t\) gives
\[
           \lambda_0(y-y_0)=x(Yz-y\Delta).
\]
For \(t\ne\mu_1\), put \(s=t-\mu_1\),
\(x_1=\lambda_1/s,\ X=x/x_1\) and \(\theta=-x_1y_0\).
There are regular functions on the fibre
\[
  U=\frac{Y-1-\theta X}{X^2}=-x_1^2 L,\qquad
  V=\frac{X-1}{Y}=-1-\frac{x(xz+\Delta)}{\lambda_1}.
\]
They give a map to the surface defined by
\[
                X^2 U=Y-1-\theta X,\qquad
                YV=X-1 .                                   \tag{1}
\]
The opens \(X\ne0\) and \(Y\ne0\) cover (1).
Recover \(x=x_1X,\ y=y_0+xL,\ L=-U/x_1^2\).
On \(X\ne0\) recover \(z\) using the formula for \(V\);
on \(Y\ne0\) use \(z=(\lambda_0 L+y\Delta)/Y\).
These formulas agree on the overlap using (1), and substitution
gives the two inverse morphisms. In particular (1) is the
affine ring obtained by adjoining the two displayed fractions
to \(\mathbb C[X,Y]\).

By [Pełka–Raźny, Theorem 1.1 and Construction 1.3](sources/pelka-razny-2002.10995v2.pdf),
this is exactly \(S_{p_1,p_2}\) for the monic polynomials
\(p_1(T)=T+\theta\), \(p_2(T)=1\). Their theorem gives the
stated algebraic invariants and the classification by the
unordered pair of polynomials. In our case their different
degrees distinguish the two members of the pair. Hence a
nonconstant Möbius parameter \(\theta(t)\) gives pairwise
nonisomorphic fibres. When \(\mu_0=\mu_1\), the formulas for
\(X,Y,U,V\) and the inverse are regular over \(\mathbb C[t,(t-\mu_1)^{-1}]\)
and give the product with the fixed surface.

For \(t=\mu_1\), the function \(Y\) is a unit on the fibre.
Indeed \(Y=0\) would give \(xy=1\) and \(\lambda_1y=0\).
The coordinates \(x,L\) identify this fibre with
\[
                  1-y_0x-x^2 L\ne0
\]
in \(\mathbb A^2\). The inverse uses \(y=y_0+xL\),
\(Y=1-xy\), and \(z=(\lambda_0 L+y\Delta)/Y\).

For the affine-line assertion, let \(X,Y,U,V\) in (1) be
polynomials in a single parameter. If \(X,Y\) are both
nonconstant, then \(Y\mid X-1\) gives \(\deg Y\le\deg X\).
Hence \(X^2\mid Y-1-\theta X\) implies \(U=0\).
Now \(Y=1+\theta X\) divides \(X-1\), forcing \(\theta=-1\)
and \(V=-1\). This gives the extra line \(U=0,\ Y=1-X\)
in that case. If \(X\) is nonconstant and \(Y\) is constant,
the same degree comparison forces \(\theta=0,\ Y=1,\ U=0\),
giving the extra line \(V=X-1\). If \(X\) is constant and
\(Y\) is nonconstant, (1) forces \(X=1,\ V=0\).
When \(X,Y\) are both constant, a nonconstant map is possible
precisely on the two exceptional lines \(X=0,Y=1\) and
\(X=1,Y=0\). Each listed locus is a copy of the affine
line, and a nonconstant polynomial map into such a line
is surjective over \(\mathbb C\).

## Exclusion of a mate in the base variables

Suppose \(P=P(x,y)\) in a Keller triple with determinant
\(\kappa\ne0\). Work over \(\mathbb C\) by embedding the
finitely generated coefficient field. The differential \(dP\)
is nowhere zero on the base plane. By the closed-generator
decomposition of [Wang–Yan, Proposition 2.2 and Lemma 2.3](sources/wang-yan-2610.36495v1.pdf),
\(P=h(p)\) with \(p\) closed. If \(h'\) had a zero, the
nonconstant polynomial \(p\) would attain that value and
\(dP\) would vanish. Thus \(h\) is affine linear and the
rational constants of \(\{P,\cdot\}_{x,y}\) are \(\mathbb C(P)\).

With \(t\) an indeterminate set
\(\widehat Q=Q(x,y,(t-b)/a)\).
It belongs to \(\mathbb C[x,y,a^{-1},t]\). On the open
\(a\ne0\), changing source coordinates to \(x,y,t\) gives
\[
                  \{P,\widehat Q\}_{x,y}=\kappa/a.          \tag{2}
\]
For either divisor \(D\) defined by \(x=0\) or \(xy=1\), the
Hamiltonian derivation in (2) must be tangent to \(D\).
Otherwise at its generic point the derivative of a pole
of order \(m\ge1\) has order \(m+1\), while the derivative
of a regular function stays regular. Neither gives the
simple pole in (2). Tangency means that \(P|_D\) is constant.
For the line this follows from \(P_y=0\) there, and on
the hyperbola it follows from
\(x\frac{d}{dx}P(x,1/x)=xP_x-yP_y=0\).

Writing \(\widehat Q=\sum Q_j(x,y)t^j\) in (2) gives
\(Q_j=f_j(P)\) for rational \(f_j\) when \(j>0\).
The constant coefficient obeys
\[
            Q_0=Q(x,y,0)-\sum_{j>0} f_j(P)b^j .
\]
Let \(p_D=P|_D\). If all \(f_j\) were regular at \(p_D\),
then \(\widehat Q\) would be regular along \(D\).
The derivation in (2) preserves regularity there, so at
least one \(f_j\) has a pole at \(p_D\). Since every
\(Q_j\) is regular off \(a=0\), all components of
\(P-p_D=0\) must be among the components of \(a=0\).
The fibres of \(P\) are reduced by the nowhere-vanishing
differential. If the two values \(p_D\) differed, the
fibre containing \(x=0\) would force \(P=p_D+\lambda x\),
\(\lambda\ne0\), which is not constant on \(xy=1\).
The two values are therefore equal, and
\(P=P_0+\lambda a\) with \(\lambda\ne0\).

At fixed \(a,t\), the coordinates are rational in \(x\)
since \(y=1/x+a/x^2\) on \(a\ne0\).
Equation (2) becomes
\(-\lambda x^2\partial_x\widehat Q=\kappa/a\), so
\[
                    \widehat Q=\frac{\kappa}{\lambda a x}
                                      +\Psi(a,t).
\]
In the generic \(t\)-fibre, \(a=0\) consists of the line
\(x=0,y=y_0\) and the line \(xy=1,x=x_1\),
both reduced. The first term has pole orders two and one
on them respectively, while \(\Psi(a,t)\) has the same order
on both, since \(a\) is a local parameter at both.
This contradicts regularity of \(\widehat Q\) along the two
lines of the generic fibre. The same argument interchanging
the mates gives the other exclusion. A constant target-linear
dependence among the three \(\partial_z\)-entries would
produce a mate independent of \(z\), after an invertible
target-linear change fixing \(R\), and is therefore impossible.

## Exclusion of a pair of affine mates

Suppose all three components are affine in \(z\) and write
\(F=zA(x,y)+B(x,y)\). Work over \(\mathbb C\).
The three entries of \(A\) are linearly independent by
the preceding theorem, and have no common zero by the
Keller determinant. The coefficient of \(z^2\) in that
determinant is \(\det(A_x,A_y,A)\), so the projective
direction \([A]\) has image of dimension one. Indeed, on
the chart \(A_3\ne0\) the two coordinate ratios have zero
Jacobian and are algebraically dependent. Restricting
\([A]\) to an affine line on which it is nonconstant
shows that its image curve is rational by Lüroth's
theorem.

Choose a rational parameter \(f/g\) of that direction curve,
with relatively prime \(f,g\in\mathbb C[x,y]\). Its map to
\(\mathbb P^2\) is given by binary forms \(H_i\) of a common
degree without a common zero on \(\mathbb P^1\).
Then \(A_i\) are common scalar multiples of \(H_i(f,g)\).
To see that the multiplier is constant, at any irreducible
divisor at least one of \(f,g\) is nonzero in its residue
field. Thus the specialized \(H_i(f,g)\) have no common
divisor. The multiplier has no pole, and cannot have any
divisor of zeros since \(A\) has no common zero.

Factor \(H_3(f,g)\) into members of the pencil generated by
\(f,g\). Its value is a constant times the squarefree
polynomial \(x(xy-1)\). If the two irreducible factors
belong to different members, a linear change of pencil
parameter makes \(f=x,\ g=xy-1\). There is then no
constant nonzero linear combination of \(f,g\), so
the binary forms have degree two. Since the \(A_i\) are
linearly independent, an invertible linear target change
fixing \(R\) makes \(A=(x^2,(xy-1)^2,x(xy-1))\).
With \(v=y-1/x,\ Z=x^2z\), the map takes the form
\[
               F=Z A_0(v)+B_0(x,v),\qquad A_0=(1,v^2,v).
                                                               \tag{3}
\]
The Jacobian determinant for the change \((x,y,z)\mapsto
(x,v,Z)\) equals \(x^2\).

If the two factors belong to one member of the pencil,
that member is a multiple of \(a=x(xy-1)\).
The other factors of \(H_3(f,g)\) must give constants,
and at least one such factor exists: forms of degree one
would put \([A]\) in a line. Thus the pencil has basis \(a,1\).
With \(v=a,\ Z=z\) we again have (3), with \(A_0(v)\) a
vector of three linearly independent polynomials in \(v\).
The Jacobian determinant of this change is again \(x^2\).

In both cases the new determinant is \(\kappa/x^2\),
with \(\kappa\ne0\) a scalar after the target change.
Let \(N=A_0\times A_0'\), with prime denoting \(v\)-derivative.
The coefficient of \(Z\) makes \(N\cdot B_{0,x}=0\).
Hence \(N\cdot B_0=h(v)\) for some rational \(h\), and
\(B_{0,x}=\alpha A_0+\beta A_0'\) over the rational
function field. Put \(W=N'\cdot B_0\) and \(C=N'\cdot A_0'\).
Then
\[
       W_x=C\beta,\qquad \beta(h'-W)=\kappa/x^2 .
\]
Here the second equality is the constant determinant,
using \(\det(A_0',B_{0,v},A_0)=N\cdot B_{0,v}=h'-W\).
Also \(C\ne0\) by linear independence of the three
polynomials in \(A_0\): up to sign it is their
Wronskian, whose leading coefficient is nonzero after
a constant basis change to distinct polynomial degrees.
Consequently
\[
                   (h'-W)^2=2\kappa C/x+J(v)
\]
for some \(J(v)\in\mathbb C(v)\). The right side has a
simple pole at \(x=0\) over \(\mathbb C(v)\), which is
impossible for a square. This excludes the two affine
mates after embedding any finitely generated
characteristic-zero coefficient field in \(\mathbb C\).
