# Nonbinary quartic factors on a sextic line pencil

Write \(\nu(x,y)=(x^2,xy,y^2)\). In this note the ordinary
highest tuple is \(g(x,y,z)\nu(x,y)\), where \(g\) is a
nonbinary quartic, and the basepoint of the line pencil is
\([0:0:1]\).

**Theorem SexticKellerHighTransverseQuarticFactorReduction.**
For a characteristic-zero three-variable Keller completion,
the \(z\)-degree of \(g\) can only be one or three. If it
is three, then the whole map has \(z\)-degree three. In
an affine binary source frame write
\[
               F=z^3 A(x,y)+z^2 B(x,y)+zC(x,y)+D(x,y).
\]
Then \(A=\ell(x,y)\nu(x,y)\) with \(\ell\) a nonzero
homogeneous linear form, \(B\) has no terms of degree
less than two, and \(C(0,0)=0\).

**Theorem SexticKellerNonbinaryQuarticLineWeightTwoReduction.**
Suppose \(\deg_z g=1\), write \(g=z\beta(x,y)+g_4(x,y)\),
and suppose the full map has a term with \(z\)-exponent at
least two. Increase the weight of \(z\) from one, keeping
\(x,y\) of weight one, until the first such term ties the
term \(z\beta\,\nu(x,y)\). If the resulting highest tuple,
divided by \(z\), factors through a conic with parameter
forms of weight two, then the weight of \(z\) is two.
In a binary linear and constant target frame the weighted
highest tuple is
\[
                      zy\,\nu(z+xy,y^2).
\]
In particular \(\beta\) is a cube of a linear form.

## The highest \(z\)-coefficient

Work over an algebraic closure. The
[leading-pencil differential](SEXTIC_KELLER_INDEPENDENT_LEADING_PENCIL_RESTRICTIONS.md#the-smooth-line-and-conic-cases)
requires a triple or quadruple zero of \(g\) on a generic
pencil line. If that zero is at the vertex, \(g\) has
multiplicity at least three there, giving \(\deg_z g\le1\).
Otherwise \(g\) contains a line not through the vertex
with multiplicity at least three, giving
\(\deg_z g\ge3\). Thus the only possible positive
\(z\)-degrees are one, three and four.

Suppose \(s=\deg_z g\ge3\). If the map has a term with
\(z\)-exponent greater than \(s\), increase the weight
\(w\) of \(z\) to its first tie, i.e. the minimum of
\((6-s-a-b)/(c-s)\) over the occurring terms
\(x^a y^b z^c\) with \(c>s\). This is greater than one.
The highest tuple has weight \(6-s+sw\) and is
divisible by \(z^s\). Its quotient restricts at \(z=0\)
to the leading binary \(z\)-coefficient of
\(g\,\nu(x,y)\). The vanishing of the highest determinant
forces its projective image to be the same conic, since
the restriction has dense image in that conic.
Factor the quotient as \(G\nu(P,Q)\) with coprime
weighted homogeneous \(P,Q\). Their common weight \(D\)
is a positive integer by restriction to \(z=0\), so
\(6-s=\operatorname{wt}G+2D\) gives \(D=1\) and
\(\operatorname{wt}G\le1\). Thus the quotient is binary,
contradicting the tying \(z\)-dependent term.

Write \(F=z^s A+z^{s-1}B+z^{s-2}C+\cdots\). The coefficient
of \(z^{3s-1}\) in the determinant implies that the
projective image of \(A\) is contained in the same conic:
the equation is \(\det(A_x,A_y,A)=0\), and the highest
binary parts already give dense image in that conic.
Writing \(A=G\nu(P,Q)\) with coprime binary polynomials
gives \(6-s=\deg G+2\max(\deg P,\deg Q)\).
The parameter forms are affine linear with independent
linear parts, since their highest ratio is nonconstant.
After centering and normalizing, \(A=\nu(x,y)\) if \(s=4\),
or \(A=(\ell(x,y)+c)\nu(x,y)\) for \(s=3\).
The use of the leading \(z\)-coefficient is related to the
quartic argument of [Wang–Yan, Lemma 3.6](sources/wang-yan-2610.36495v1.pdf).

Let \(A_d\) be the lowest nonzero homogeneous part of \(A\),
of degree \(d=2\) or \(3\), and let \(B_i\) be the
degree-\(i\) part of \(B\). The coefficient of \(z^{3s-2}\)
first gives \(B_0=0\). For \(i=0\), and once \(B_0=0\)
also for \(i=1\), its relevant binary homogeneous equation
contains a factor
\[
                        (s-1)d-si
\]
multiplying the normal contraction of \(B_i\) with the
two derivative columns of \(A_d\).
Indeed take \(x=u,\ y=ut\) and write \(A_d=u^d a(t)\),
\(B_i=u^i b(t)\). The coefficient comes from
\[
 (s-1)\det((A_d)_u,(A_d)_t,B_i)
 +s\det((A_d)_u,(B_i)_t,A_d)
 +s\det((B_i)_u,(A_d)_t,A_d).
\]
The middle term is zero, and the radial factors of
the other two yield the displayed factor.
It is nonzero in both degrees for the indicated
\((s,d)\). A constant vector whose normal contraction
vanishes is zero, as the conic normals span the three-dimensional
target dual space.
The degree-one vector tangent to the conic has the form
\[
                         B_1=a\nu_x+b\nu_y
\]
with \(a,b\in k\). This follows by solving the conic
tangent equation \(y^2 B_{1,1}-2xy B_{1,2}+x^2B_{1,3}=0\).

## Exclusion with a quadratic lowest coefficient

Suppose \(d=2\) and rescale so \(A_2=\nu\). Substitute
\(x=\epsilon u,\ y=\epsilon v,\ z=\epsilon^{-1}\) and write
\(F=\epsilon^{-\alpha}H(u,v,\epsilon)\), where
\(\alpha=s-2\in\{1,2\}\). Because \(B_0=0\), the
series \(H\) is polynomial in \(u,v,\epsilon\).
Its coefficient at \(\epsilon=0\) is
\(\nu(u,v)+B_1(u,v)+C_0\).
A constant shift of \(u,v\) removes the linear tangent
part, leaving \(\nu+c_0\) for some constant vector \(c_0\).
The Keller identity in these coordinates is
\[
           \det(H_u,H_v,\alpha H-\epsilon H_\epsilon)
                              =\kappa\epsilon^{3\alpha+1}.
 \tag{1}
\]
The coordinate substitution has determinant \(-1\), and
the constant shift of \(u,v\) leaves that determinant
unchanged. At order zero (1) forces \(c_0=0\): the
contraction with the conic normal is zero.

Use a reference of the form
\[
                  H_{\rm ref}
       =\nu(U,V)+\epsilon^{\alpha/2} M(U,V)+\epsilon^\alpha c
\]
where \(U,V\) are series with leading \(u,v\),
\(c\) is a constant vector, and \(M\) is linear when
\(\alpha=2\), and absent when \(\alpha=1\).
Its determinant on the left of (1) is zero because all
three columns are combinations of its derivative columns
in \(U,V\). For the third column this follows by applying
\(\alpha-\epsilon\partial_\epsilon\) to the three displayed
terms, first with \(U,V\) fixed, and then by the chain rule.

If the first discrepancy is \(\epsilon^r E(u,v)\),
its linear contribution to the determinant at that order is
\[
   (\nu_u\times\nu_v)\cdot
       \left((\alpha-r)E-\tfrac{\alpha}{2}(uE_u+vE_v)\right).
 \tag{2}
\]
This uses the Euler identity \(2\nu=u\nu_u+v\nu_v\).
For \(1\le r<3\alpha+1\), the right side of (1) vanishes.
In each binary degree \(i\) the factor in (2) is
\(\alpha-r-\alpha i/2\). When it is nonzero, the
contraction vanishes and the error is a polynomial
combination of \(\nu_u,\nu_v\), by the coprime conic minors.
It can be absorbed by corrections to \(U,V\).
The possible zero factors have \((r,i)=(\alpha,0)\)
or \((1,1)\) for \(\alpha=2\), and are absorbed by
\(c\) and \(M\) respectively. This raises the first
error order to \(3\alpha+1\). At that order the left
side is again (2), which vanishes at \(u=v=0\) because
the conic normal vanishes there. The right side is
\(\kappa\ne0\), a contradiction.

Thus \(s=4\) is impossible, and \(s=3\) requires \(c=0\)
in the centered \(A=(\ell+c)\nu\).
Take \(\ell=x\) by a binary source change and a compensating
target change preserving the conic. For this cubic \(A\)
use the coefficient of \(z^6\).
Its binary degree-three part comes only from \(A,B_1,B_1\)
and is
\[
 2\det(A_x,(B_1)_y,B_1)
 +2\det((B_1)_x,A_y,B_1)
 +3\det((B_1)_x,(B_1)_y,A)
                    =-6x(bx-ay)^2.
\]
The equality follows by substituting
\(A=(x^3,x^2y,xy^2)\) and
\(B_1=(2ax, ay+bx,2by)\).
Hence \(B_1=0\). The degree-four part of the same coefficient
now reduces to \(\det(A_x,A_y,C_0)=0\), so \(C_0=0\).
This proves the first theorem.

## Weight-two parameters when the quartic factor is linear in \(z\)

Let \(s=1\) and use the threshold from the second theorem.
After division by \(z\) the leading tuple has weight five,
and at \(z=0\) it is \(\beta\nu(x,y)\).
The projective image argument at the threshold applies as
above. If the parameter weight is two, the common factor
has weight one and is a binary linear form \(G\).
At least one parameter form must involve \(z\), and it has
weight two with a nonzero constant \(z\)-coefficient.
Thus \(w=2\). The parameters restrict to \(xh,yh\) on
\(z=0\) for some binary linear \(h\), by equality of their
projective ratio there. Rotating \(x,y\), then rescaling
\(z\) and compensating on the target, gives the entire
weighted tuple as
\[
                      z G\,\nu(P,Q),\qquad
                      P=xh+z,\quad Q=yh .
\]
In particular \(\beta=Gh^2\).

Use the [weighted leading-pair slice](KELLER_LEADING_PAIR_RATIONAL_SLICE.md)
with \(S=zG P^2\) and \(S T\), \(T=Q/P\), as the pair.
At generic fixed nonzero \(T\), \(z=h(y/T-x)\), so
\(S\) becomes the binary form \(G h^3 y^2(y-Tx)/T^3\)
of degree seven. The relative differential on this
plane level equals its area differential multiplied
by \(Q\) up to a nonzero scalar on the level:
division by \(dT\) contributes \(P^2/Q=Q/T^2\),
and the other factor from \(d(ST)\) is constant on
the chosen level.
Using coordinates radial \(R\) and projective \(u\)
on the binary plane, its level is \(R^7 f(u)=a\),
where \(f\) is the dehomogenization of
\(G h^3 y^2(y-Tx)\) and \(a\ne0\).
Since the binary area differential is proportional
to \(R^2\,du\) on this level, the differential to
be integrated is proportional to \(R^4 M(u)\,du\),
where \(M\) is the dehomogenization of \(yh\).

Choose the affine chart in \(u\) so all linear factors
are finite. The \(u\)-root from \(y-Tx\) is simple and
distinct from the others for generic \(T\). The extension
by \(R\) has degree seven, since \(a/f\) has a simple
pole there. If the differential has a rational primitive,
take its component under the cyclic action \(R\mapsto
\zeta R\) with character \(\zeta^4\), where \(\zeta\)
is a primitive seventh root of unity. This gives a
primitive \(R^4 L(u)\) with \(L\) rational, so
\[
                      L'-\tfrac47(f'/f)L=M
\]
after multiplying by a nonzero scalar.
At a finite root of \(f\) of multiplicity \(m\), \(1\le m\le6\),
the coefficient on the leading local power of \(L\) is
its order minus \(4m/7\), which is never zero.
It follows that \(L\) has no pole there and its
vanishing order must be one more than the order of \(M\)
there. It has no finite pole away from \(f=0\) either,
since differentiation would give a pole of higher order.
At infinity, a polynomial \(L\) of degree \(d>4\)
gives a term of degree \(d-1\ge4\) with coefficient
\(d-4\), incompatible with the degree-two right hand side.
Thus \(\deg L\le4\).

The zero orders at the distinct roots of \(f\) sum to
their number plus two, since both factors of \(M\) occur
among them. Hence \(f\) has at most two distinct roots.
As \(y-Tx\) is independent of all the fixed factors,
\(G,h,y\) are proportional. Rescaling the nonzero
coefficients produces \(zy\,\nu(z+xy,y^2)\) as stated.
