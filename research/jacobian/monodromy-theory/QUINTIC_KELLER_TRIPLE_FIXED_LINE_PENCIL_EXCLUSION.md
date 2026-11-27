# Excluding a triple fixed line in top rank one

**Theorem QuinticRankOneTripleFixedLineExclusion.** A noninvertible
degree-five three-variable Keller map of top rank one cannot have
primitive cancelling-pencil degree one with a triple fixed line
in the common factor of its quartic leading components.

We use
[QuinticRankOneLinePencilCubicRestriction](QUINTIC_KELLER_RANK_ONE_LINE_PENCIL_CUBIC_RESTRICTION.md)
and the
[primitive pencil representation](QUINTIC_KELLER_REDUCED_BOUNDARY_PENCIL_RESIDUAL.md#rank-one-pencil).
After linear changes choose
\[
                  Q_4=z^4,\qquad R_4=xz^3 .
\tag{1}
\]
The polynomial \(P_5\) is binary in \(x,z\), since all its
irreducible factors belong to members of the line pencil.
Both \(Q_3\) and \(R_3\) are divisible by \(z\).
We work over \(\mathbb C\) and transfer the conclusion to
characteristic zero by coefficient-field embedding and descent
of the polynomial inverse.

## A square restriction from the line blowup

For a component in the cancelling pencil let \(U\) be its cubic
term divided by \(z\) and restricted to \(z=0\). Blowing up the
line \(z=w=0\) in the projective fibre, use the chart \(w=zW\),
where \(w\) homogenizes the affine source. In a chart on the
line the original fibre equation has terms
\[
              z^3 L+wzU+w^2 V+\text{higher terms in }z,w .
\]
Here \(L\) is the other linear leading factor and \(V\) is the
restriction of the homogeneous quadratic term to the line.
The strict-transform equation after division by \(z^2\)
restricts to \(W(U+WV)=0\) at \(z=0\).
If \(U\ne0\), the section \(W=0\) of this exceptional curve
is smooth in the surface over the open \(U\ne0\). The
quartic dualizing form there has a simple pole with residue
a nonzero constant times \(dx/U(x,1)\). Indeed the ambient
form gains a factor \(z\) from the blowup while the fibre
equation gains \(z^2\), and its derivative in \(W\) on the
section is \(U\). This residue is a rational differential
on the projective line. If the binary quadratic \(U\) has
distinct roots, its simple poles give nonzero periods on
small loops. The shrinking normal tubes over these loops
contradict the
[exact fibre form identity](KELLER_QUARTIC_COMPONENT_BOUNDARY_OBSTRUCTION.md#residues-on-a-resolution).
The loops avoid the remaining exceptional components and
strict infinity, so the tubes lie in the original affine
fibre. Thus \(U\) is zero or a square up to scalar.

This conclusion holds for every general cancelling component.
A pencil contained in the cone of binary squares has
one-dimensional span or zero span, since the sum of two
independent squares of linear forms has rank two.
The two cubic restrictions therefore have span at most one.
We can choose the components so that
\[
                    Q_4=z^4,\qquad R_4=xz^3,\qquad z^2\mid R_3,
\tag{2}
\]
although the pure fourth power requires consideration if its
own restriction initially vanishes. If the restriction of the
pure fourth power member is nonzero, subtract a multiple of
it from a second member and adjust \(x\) by a multiple of \(z\).
If its restriction and the second restriction are both zero,
(2) is already available.

The outstanding possibility is that the pure fourth power
member \(Q\) has zero restriction while the other restriction
is nonzero. In this case \(z^2\mid Q_3\) and \(Q_4=z^4\).
A nonzero constant \(y^2\) term in \(Q_2\) gives
[quadratic-source completion](KELLER_QUARTIC_QUADRATIC_SUSPENSIONS_AND_HIGHER_CONTACT.md#the-rank-one-normal-jet-of-a-quartic-component).
Otherwise \(Q\) is affine linear in \(y\). A nonzero reduced
coefficient gives inversion by
[reduced quadratic completion](KELLER_QUARTIC_TRIPLE_BOUNDARY_AND_LINEAR_COMPONENTS.md#invertibility-with-a-reduced-quadratic-coefficient),
and a nonzero repeated coefficient gives the repeated-line
residual in the
[linear-source classification](KELLER_LINEAR_SOURCE_COEFFICIENT_CUBIC_CLASSIFICATION.md).
Its leading \(z^4\) forces
\(Q=\mu+\lambda z+z^2(z^2+\alpha x+\beta y+\gamma z+\delta)\)
after translating \(z\), with \(\lambda\beta\ne0\).
The case that \(Q\) is independent of \(y\) has the same residual
with \(\beta=0,\ \alpha\ne0\): a nonzero \(x^2\) term in \(Q_2\)
gives quadratic-source completion, and otherwise the coefficient
of \(x\) depends only on \(z\) and is at most quadratic.
Noncriticality excludes a simple root, while a nonzero constant
gives a graph coordinate. This proves the asserted quartic
exceptional form whenever \(Q\) has zero restriction.

We first treat (2), and then return to this outstanding case.

## The resulting scalar \(R\)

Under (2), a nonzero \(y^2\) coefficient in \(R_2\) gives
quadratic-source completion. Otherwise \(R\) is affine
linear in \(y\), with coefficient of degree at most two
in \(x,z\). If the coefficient is nonzero, the only
noninvertible possibility from the linear-source
classification has a repeated square coefficient in \(z\)
and a constant restriction to that line. The quadratic
part of the factor \(H\) in the repeated-line residual
\(R-\mu=\lambda z+z^2H\) is \(xz\) by (2).
The \(y\) coefficient in \(H\) is nonzero. An affine change
in the binary variables and a target scaling then make
\[
                              R=z+z^2(y+zx).
\tag{3}
\]
The change preserves (1), up to the nonzero scalar
coefficient in \(R_4\), and \(P_5\) remains binary.

If \(R\) is independent of \(y\), its constant \(x^2\)
coefficient, if nonzero, again gives quadratic-source
completion. Otherwise \(R=a(z)x+b(z)\) with \(a\) monic cubic.
Noncriticality excludes every simple root of \(a\), since
one could solve \(a'(z)x+b'(z)=0\) there. Thus \(a\) is a cube
and \(b'\) is nonzero at its root. Translating \(z\) and \(x\)
and scaling gives the planar scalar
\[
                              R=z+d z^2+xz^3.
\tag{4}
\]

## The exceptional bracket restrictions

We will use the following elementary determinant restrictions
for \(r=z+c z^2+X z^j\) with \(j=2\) or \(3\), and a free
variable \(Y\). Put \(D=-z^j\partial_z+r_z\partial_X\).
The restriction of \(D\) on \(\mathbb C(r,z)\) is
\(-z^j\partial_z\) at fixed \(r\). Its rational kernel is
\(\mathbb C(r)\), and no polynomial \(D\)-primitive of a
nonzero constant exists. Indeed a rational primitive of
one is \(z^{1-j}/(j-1)\), so a polynomial primitive would
differ from a nonzero multiple of this function by
\(A(r)/C(r)\). In a reduced fraction \(C(r)\) would have
to divide \(z^{j-1}\), which is incompatible with the degree
\(j+1\) of \(r\), unless \(C\) is constant. This would leave
an uncanceled pole at \(z=0\).

The bracket \(D(S)T_Y-S_YD(T)\) cannot be a nonzero constant
if one component is independent of \(Y\), if one has nonzero
constant leading \(Y\)-coefficient, or if both are affine in
\(Y\). The first case gives a polynomial \(D\)-primitive of a
constant immediately. For the second, normalize a monic
\(S\) of degree \(m\) in \(Y\). Use the formal series
\(W=S^{-1/m}\) with leading term \(Y^{-1}\) in
\[
 D(TS_YW)-\partial_Y(TD(S)W)=-\kappa W .
\]
Taking the \(Y^{-1}\)-coefficient gives a forbidden
polynomial primitive. Finally, affine leading coefficients
\(p,q\) must have ratio in \(\mathbb C(r)\) by the top
bracket equation. They cannot have a common zero, since
the bracket is a nonzero constant. Thus their ratio in
reduced form forces \(p,q\in\mathbb C[r]\), up to a common
constant. The constant equation then gives the same
forbidden primitive. These are the degree-independent
arguments of
[Wang–Yan, Proposition 3.17](https://arxiv.org/pdf/2610.36495v1)
for the exceptional cubic, which depend here only on the
kernel and primitive calculations just given.

## The planar scalar

For (4) use free variable \(y\) and the preceding bracket
with \(j=3\). Write \(m=\deg_y P,\ n=\deg_y Q\), with
top coefficients \(p_m,q_n\). Then \(1\le m\le4\),
\(1\le n\le2\), and
\[
 \deg p_m\le4-m,\quad \deg q_n\le3-n,\quad
                    p_m^n/q_n^m\in\mathbb C(R).
\tag{5}
\]
The degree bounds use that \(P_5,Q_4\) are independent of
\(y\), and that \(Q_3\) is divisible by \(z\). The ratio
is constant: a nonconstant rational function would require
a full fibre factor \(R-c\) to divide one of the leading
coefficients, by a reduced numerator-denominator identity.
Each fibre is squarefree by submersion, and its degree
four exceeds both coefficient degrees.

For \(n=2\), \(q_2\) has degree one by the nonconstant
leading-coefficient requirement. Equation (5) forces
\(m=2\), since \(m\) must be even and \(p_4\) would be
constant. Equal degrees allow cancellation by a constant
multiple of \(Q\) to get \(m=1\), which violates parity in
(5). For \(n=1\), \(q_1\) is nonconstant and degree
comparison forces \(m\le2\). At \(m=2\), subtract a constant
multiple of \(Q^2\) from \(P\) using (5). This gives two
affine components or a component independent of \(y\),
both forbidden. This closes (4).

## The nonparallel scalar in (2)

For (3) shear \(Y=y+zx\) and keep the free variable \(x\).
Now use \(r=R=z+z^2Y\) with \(j=2\) in the exceptional
bracket restrictions. Measure \(m=\deg_x P,\ n=\deg_x Q\)
after the shear and write \(p_m,q_n\) for their leading
coefficients in \(z,Y\). We have \(n\le2\). For any \(m\ge1\),
\[
                        \deg_z p_m\le4.
\tag{6}
\]
Indeed the original \(P_5\) is binary in \(x,z\), and
a term of smaller ordinary degree can gain at most one
power of \(z\) for each substituted \(y\) replaced by
\(-zx\). Also \(m=5\) would make \(p_m\) constant, and
\(p_4\) is independent of \(Y\).

If \(Q_3/z\) restricted to \(z=0\) is a nonzero square involving
\(y\), write it as \((u x+v y)^2\) with \(v\ne0\), absorbing
its scalar into the linear form. Then \(n=2\) and \(q_2\)
is a polynomial in \(z\) of degree three, since its leading
term is \(v^2 z^3\). The top bracket condition gives
\(p_m^2=q_2^m H(r)\) for a polynomial \(H\). A denominator of
\(H\) would require a full factor \(r-c\) to divide \(q_2^m\),
and every such fibre has a factor depending on \(Y\).
For odd \(m\), the highest \(Y\)-coefficient of the right
side has odd degree in \(z\), since that degree is
\(3m+2\deg H\). It cannot be a square. For \(m=4\),
\(p_4\) independent of \(Y\) forces \(H\) constant, and
\(\deg_z p_4=6\) contradicts (6). For \(m=2\), taking the
ratio \(p_2/q_2\) directly gives \(p_2=q_2 K(r)\), with \(K\)
polynomial by the same denominator argument. The bound
(6) forces \(K\) constant. Cancel the proportional leading
term using \(Q\) to reduce to the already impossible case
of odd or zero \(m\).

It follows that \(Q_3/z\) restricts either to zero or to a
nonzero multiple of \(x^2\). If \(Q_2\) has nonzero \(y^2\)
coefficient, quadratic-source completion applies.
Otherwise \(Q\) is affine linear in \(y\), with coefficient
of degree at most two whose degree-two term is a multiple
of \(z^2\). If that coefficient is nonzero, the
linear-source classification makes \(Q\) an exceptional
repeated-line scalar
\[
 Q=\nu+\eta(z-\xi)+(z-\xi)^2
                    (z^2+\mu x+\rho y+\alpha z+\beta).
\tag{7}
\]
Here we used \(Q_4=z^4\) to fix the quadratic part, and
\(\rho\ne0\). In particular the \(zx^2\) coefficient
vanishes. If \(Q\) is independent of \(y\), it has the
form
\[
             Q=f(z)+A(z)x+(kz+d)x^2 ,
             \qquad \deg A\le2,\quad f\text{ monic quartic}.
\tag{8}
\]
If \(k\ne0\), translate \(z\) to absorb \(d\).
For \(A(0)=0\), \(Q_x\) vanishes all along \(z=0\)
and \(Q_z\) is a quadratic in \(x\) with nonzero leading
coefficient there, giving a critical point.
For \(A(0)\ne0\), solve \(Q_x=0\) on \(z\ne0\). Then the
remaining equation is the derivative of
\(f(z)-A(z)^2/(4kz)\). Its numerator after multiplication
by \(z^2\) has degree five and nonzero constant coefficient.
It has a nonzero root, again giving a critical point.
Thus \(k=0\). A nonzero \(d\) gives quadratic completion.
If \(d=0\), \(A\) must have degree two and a double root
for noninvertibility, by the graph and simple-root tests
used above. We get (7) with \(\rho=0,\ \mu\ne0\).

For (7), after the shear \(y=Y-zx\) the free-variable
degree is \(n=1\) with coefficient
\(q_1=(z-\xi)^2(\mu-\rho z)\) of degree two or three.
The highest determinant coefficient now gives
\(p_m=q_1^m H(r)\) with \(H\) polynomial by the same
denominator argument. The bound (6) forces \(m\le2\).
At \(m=2\) it forces \(H\) constant, and subtracting a
multiple of \(Q^2\) yields the forbidden affine or
independent case. This excludes (3) under (2).

## When the pure fourth power already has zero cubic restriction

In the outstanding case \(z^2\mid Q_3\), and \(R_3/z\) restricts
to a nonzero square on \(z=0\). As shown above, the pure fourth-power member
has the repeated form
\[
                Q=\mu+\lambda z+z^2(S+g(z)),
\]
after translating \(z\), with \(\lambda\ne0\), a nonzero binary
linear form \(S\), and \(g\) a univariate quadratic. Rescaling
\(Q\) and \(S+g\), shear \(X=S+g(z)\) so that \(Q=z+z^2X\),
after dropping the target constant. Use the free variable
\(T=x\) if \(S\) has a \(y\)-coefficient, and \(T=y\) otherwise.
The original \(x\) either is \(T\) or is independent of \(T\)
and affine in \(S\).

Write \(m=\deg_T P,\ n=\deg_T R\) after the shear, with leading
coefficients \(p_m,r_n\) in \(z,X\). The exceptional bracket
restrictions give \(m,n\ge1\). The binary leading quintic
contributes to a positive-degree \(T\)-coefficient only if
\(T=x\), and that contribution is independent of \(X\).
All other contributions come from terms of ordinary degree
at most four in the original frame. Replacing \(S\) by
\(X-g(z)\) in such a term gives
\[
 \deg_z p_m\le\max(5-m,2(4-m)),\qquad
 \deg_X p_m\le4-m .
\tag{9}
\]
Only the binary quintic can give \(m=5\), and then its
coefficient is constant and impossible. For \(m=2\), the
coefficient of \(X\) in \(p_2\) has degree in \(z\) at most
two, since after consuming two powers of \(T\) and one
of \(X\) at most one ordinary degree remains before the
quadratic substitution.

The leading quartic \(xz^3\) of \(R\) contributes at most
one power of \(T\). Its homogeneous cubic is \(z\) times
the square on \(z=0\), plus a term divisible by \(z^2\).
Hence \(n\le2\), and when \(n=2\) its coefficient \(r_2\)
is affine in \(z\) and independent of \(X\). It must have
degree one by the constant-leading-coefficient restriction.
When \(n=1\), the \(T\) coefficient in the linear form being
squared must vanish, since otherwise that square would give
a nonzero \(zT^2\) term. Therefore
\[
                       r_1=\gamma X+h(z),
                       \qquad \deg h\le3 .
\tag{10}
\]
In (10) the \(X\) coefficient can only come from a term of
ordinary degree at most two before shear.

Neither \(r_n\) can contain a full fibre factor of \(Q-c\).
For \(r_2\) this follows from independence of \(X\).
For (10), if \(\gamma\ne0\) the coefficient of \(X\) is
a nonzero constant, whereas \(Q-c\) has \(X\)-coefficient
\(z^2\). Thus the highest coefficient relation
\(p_m^n/r_n^m\in\mathbb C(Q)\) is polynomial in \(Q\).
If \(n=2\), parity of the \(z\)-degree of the leading
\(X\)-coefficient on the right in \(p_m^2=r_2^m H(Q)\)
excludes odd \(m\). For \(m=4\), the \(z\)-bound one in
(9) contradicts the required degree at least two on
the left-hand square root. For \(m=2\), use
\(p_2=r_2 H_1(Q)\), since the ratio itself is invariant.
The bound \(\deg_z p_2\le4\) gives \(\deg H_1\le1\). A
nonconstant \(H_1\) would give the \(X\)-coefficient of
\(p_2\) degree three in \(z\), contrary to its bound two.
So \(p_2\) is proportional to \(r_2\) and cancellation
reduces to the already forbidden odd or zero \(m\).

If \(n=1\), we have \(p_m=r_1^m H(Q)\).
For \(m\ge3\), any positive \(\deg_z r_1\) contradicts
(9), and if \(\deg_z r_1=0\) then \(r_1\) must have degree
one in \(X\), also contradicting (9).
For \(m=2\) we show \(H\) is constant. If \(\deg_z r_1\)
is at least two, this follows from the bound four in
(9). If \(\deg_z r_1=0\), the degree-two bound in \(X\)
gives the result. If \(\deg_z r_1=1\), the \(z\)-bound
gives \(\deg H\le1\). A nonconstant \(H\) would give
degree three in \(X\) when \(\gamma\ne0\); when \(\gamma=0\)
it would give the \(X\)-coefficient degree four in \(z\).
Both violate the bounds above. Thus \(H\) is constant
and a multiple of \(R^2\) cancels the \(T^2\) term of
\(P\). The exceptional bracket now has both components
affine in \(T\), or one independent of \(T\), a contradiction.
This finishes the outstanding case and the theorem.
