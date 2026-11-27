# Quadratic and affine mates of the squarefree quartic scalar

**Theorem SexticKellerSquarefreeQuarticQuadraticAffineExclusion.**
In characteristic zero, the squarefree quartic normal form \(R=az+b\),
where
\[
 v=xy-1,\quad a=xv,\quad
 b=\lambda_0y+\mu_0+xy\bigl((\lambda_1-\lambda_0)y+\mu_1-\mu_0\bigr),
 \qquad \lambda_0\lambda_1\ne0,
\]
has no polynomial Keller completion of ordinary degree at most six
whose mate degrees in \(z\) are \(2,1\).

Work over \(\mathbb C\) by embedding the finitely generated
coefficient field. Write \(P=pz^2+p_1z+p_0,\ Q=qz+d\), with
\(p,q\ne0\) and \(\deg p\le4,\ \deg q\le5\).
In the chart \(z=(r-b)/a\) at fixed \(r=R\), the coefficient
of \(r^3\) in the mate bracket gives
\(\{p/a^2,q/a\}=0\). If \(q\) is a multiple of \(a\) by a
constant, the corresponding subtraction of \(R\) gives a
[forbidden planar mate](KELLER_SQUAREFREE_QUARTIC_FIBRES_AND_MATES.md#exclusion-of-a-mate-in-the-base-variables).
Thus \(q/a\) is nonconstant. A leading coefficient proportional
to \(x\) is also impossible by the
[bounded coordinate argument](SEXTIC_KELLER_SQUAREFREE_REDUCED_CUBIC_CASES.md#the-mate-with-linear-source-coefficient).

## An integration constraint at fixed affine mate

On \(q\ne0\), use \(W=Q\) as third coordinate and substitute
\(z=(W-d)/q\). At fixed \(W\), suppose that the two transformed
mates have the form
\[
 \widehat R=c(t)W+r_0,\qquad
 \widehat P=F(t)W^2+GW+H,\qquad c=a/q,\quad F=p/q^2,
\]
where \(t\) is a rational function of the base variables.
The underlying Jacobian equation is
\(\{\widehat R,\widehat P\}=\kappa/q\), with
\(\kappa\ne0\) a constant of the appropriate sign.
In each application below, \(\mathbb C(x,y)=\mathbb C(t,\xi)\)
with \(\xi=x\) or \(\xi=y\), and \(c\) is nonconstant.
Divide the equation by \(\{t,\xi\}\) to use the bracket in
\(t,\xi\). Put \(M=F'/c'\), with prime denoting \(d/dt\).
The equations at \(W^2\) and \(W\) integrate at fixed \(t\) to give
\[
 G=Mr_0+g(t),\qquad
 H=\frac{M'}{2c'}r_0^2+\frac{g'}{c'}r_0+h(t).
\]
Here the rational functions with zero \(\xi\)-derivative are
exactly \(\mathbb C(t)\). In the constant equation the terms
containing \(\partial_t r_0\) cancel. Its left side is
\[
 -\partial_\xi\!\left(
 \frac{1}{3}\Big(\frac{M'}{2c'}\Big)'r_0^3+
 \frac{1}{2}\Big(\frac{g'}{c'}\Big)'r_0^2+h' r_0\right).
 \tag{1}
\]
Its right side is \(\kappa/(q\{t,\xi\})\).

Two applications will suffice. If the right side has a simple
pole as a rational function of \(\xi\) over \(\mathbb C(t)\),
it cannot be a rational derivative. If it is a nonzero element
of \(\mathbb C(t)\) times \(1/\xi^2\), integration of (1) gives
a polynomial in \(r_0\) equal to a rational function of degree
one in \(\xi\). The polynomial must have degree one, since
degrees of nonconstant rational functions multiply under
composition. In particular \((M'/2c')'=0\). Integrating this
equality twice in \(t\), using \(M=F'/c'\), gives
\[
                            F\in\mathbb C[c]_{\le2}.          \tag{2}
\]

In particular, if \(t=v/x=y-1/x\), then \(\{t,x\}=-1\).
For \(q=x^2 f(t),\ p=x^4\), a quadratic polynomial \(f\) with
nonzero leading coefficient gives \(c=t/f,\ F=1/f^2\).
Equation (1) applies with \(\xi=x\) and right side proportional
to \(1/x^2\). Equation (2) would imply \(1\) is a scalar linear
combination of \(f^2,tf,t^2\), whose respective degrees are
\(4,3,2\). This is impossible.

Two other applications use \(q=v f(x)\) or \(q=x f(v)\).
For \(q=v f(x)\), with \(t=x,\xi=y\), the right side has a
simple pole at \(y=1/x\). It excludes this case whenever \(p/q^2\)
depends on \(x\) alone. For \(q=x f(v)\), take \(t=v,\xi=x\).
Then \(\{t,x\}=-x\), and (1) has right side proportional to
\(1/x^2\). If \(p=x^2(\gamma v+\delta)\), equation (2) requires
\[
          \gamma v+\delta \in
                \operatorname{span}_{\mathbb C}\{f^2,vf,v^2\}.
 \tag{3}
\]
For quadratic \(f\) this is impossible by the distinct degrees
\(4,3,2\). For nonconstant affine \(f\) with \(v\mid f\), a
constant change of the mate gives a planar mate. For affine \(f\)
with \(v\nmid f\), subtract a constant multiple of \(R\) to
obtain a mate with leading coefficient proportional to \(x\),
which is also impossible.

## Leading-coefficient cases

We use the
[closed rational pencil criterion](sources/petravchuk-iena-math-0701588v2.pdf),
Lemma 2 and Theorem 1, in the same form as in the
[high source-degree argument](SEXTIC_KELLER_SQUAREFREE_MATE_DEGREES.md#leading-coefficients-in-the-r-fibre-chart).
A quotient \(u\) is closed when \(\mathbb C(u)\) is relatively
algebraically closed in \(\mathbb C(x,y)\). Then
\(\{u,p/a^2\}=0\) implies \(p/a^2\in\mathbb C(u)\).
If \(u=q/a\) is closed, has a simple pole along a component of
\(a=0\), and every finite fibre meets \(a\ne0\) in a curve, then
\(p/a^2=\sum_{i=0}^2 k_i u^i\). Indeed, regularity off \(a=0\)
rules out finite poles and the simple pole of \(u\) bounds the
degree by two. Subtract \(\sum_{i=0}^2 k_i R^{2-i}Q^i\) from
\(P\). This lowers its \(z\)-degree to at most one and invokes
the [exclusion of two affine mates](KELLER_SQUAREFREE_QUARTIC_FIBRES_AND_MATES.md#exclusion-of-a-pair-of-affine-mates),
which has no ordinary degree bound.

**Neither factor of \(a\) divides \(q\).** If \(q\) and \(a\)
are algebraically dependent, the kernel of \(a\) is
\(\mathbb C[a]\), and the degree bound on \(q\) makes it
constant or affine in \(a\). The cases of constant leading
coefficient and nonconstant affine leading coefficient with
nonzero constant term are already excluded by the
[formal root exclusions](SEXTIC_KELLER_SQUAREFREE_MATE_DEGREES.md#constant-leading-coefficient)
and the
[two-divisor obstruction](SEXTIC_KELLER_SQUAREFREE_REDUCED_CUBIC_CASES.md#a-leading-coefficient-nonvanishing-on-both-components).
In the independent case, \(q-ca\) is nonconstant for every
scalar \(c\) and coprime to \(a\). Thus if \(u=q/a\) is closed,
the preceding subtraction applies.

If \(u\) is nonclosed, write it as a quotient of coprime binary
forms of common degree \(D\ge2\) evaluated at coprime
polynomials. The evaluated numerator and denominator are
coprime. No nonzero member of the polynomial pencil can be
constant: then \(q,a\) would both be polynomials in one
generator, contrary to independence. The denominator is
squarefree with only the factors \(x,v\). Thus \(D=2\), and
the denominator comes from two distinct pencil members
proportional to \(x,v\). After subtracting a multiple of \(R\)
the leading coefficient becomes \(q=\eta x^2+\nu v^2\), with
\(\eta\nu\ne0\). The rational constants of \(u\) are
\(\mathbb C(t)\), \(t=v/x\), since the full base field is
\(\mathbb C(t)(x)\). Thus \(p/a^2\) is Laurent polynomial in
\(t\) with exponents from \(-2\) to \(2\): finite nonzero fibres
meet \(a\ne0\), and \(t\) has simple zero and pole on the two
divisors of \(a\). The corresponding terms of \(p\) have
ordinary degrees \(6+i\) for exponent \(i\). Only \(i=-2\)
is allowed, so scale \(p\) to \(x^4\). This is the excluded
case \(q=x^2 f(t)\) of (2).

**Both factors divide \(q\).** Write \(q=ah_0\) with \(h_0\)
nonconstant of degree at most two. A polynomial of degree at
most two that is a nontrivial composition is a polynomial
in an affine linear polynomial, by the
[closed polynomial generator theorem](sources/wang-yan-2610.36495v1.pdf),
Proposition 2.2 and Lemma 2.3. Thus choose a closed polynomial
generator \(h\) of degree at most two so that \(h_0\in\mathbb
C[h]\), and \(p/a^2\in\mathbb C(h)\) by the bracket. This
rational function must have a finite pole, since otherwise
\(p\) is divisible by \(a^2\). All components of the corresponding
fibre of \(h\) must be supported on \(a=0\). A polynomial of
degree at most two supported on \(xv=0\) is a scalar times
\(x,x^2\) or \(v\). The \(x^2\) case contradicts closedness of
\(h\). Thus after affine change of the parameter \(h\) it is
\(x\) or \(v\). If \(h=x\), then \(q=v f(x)\) and \(p/q^2\)
depends only on \(x\), excluded by the simple pole in (1).
If \(h=v\), then \(q=x f(v)\), with \(\deg f\le2\), and
\(p=x^2(\gamma v+\delta)\). In fact \(p/a^2\) can only have
a pole of order at most two at \(v=0\); the degree bound on
\(p\) then gives the displayed form. This case is excluded by
(3) or the constant and affine \(f\) cases.

**Only \(x\) divides \(q\).** Write \(q=xh\), \(v\nmid h\).
If \(h\) and \(v\) are dependent, then \(h\in\mathbb C[v]\).
Thus \(q=x f(v)\) with \(\deg f\le2\), and the same kernel
argument gives \(p=x^2(\gamma v+\delta)\). Apply (3).
Otherwise \(h/v\) is closed by the irreducible denominator
criterion. We can use the earlier subtraction unless some
finite fibre is supported on \(x=0\) or is empty, in which
case \(h-cv\) is a nonzero scalar times \(x^j\) with \(0\le
j\le4\). Subtract \(cR\) and scale so \(q=x^{j+1}\). The
case \(j=0\) is excluded. For \(j\ge1\) the quotient
\(u=x^j/v\) is closed since \(\mathbb C(x,y)=\mathbb C(u)(x)\).
At \(x=0\) it vanishes to order \(j\), and at \(v=0\) it
has a simple pole. Every finite nonzero fibre meets \(a\ne0\).
Thus \(p\) is a linear combination of
\[
       a^2u^i=x^{2+ji}v^{2-i},\qquad
                     -\lfloor 2/j\rfloor\le i\le2 .
\]
The degrees are \(6+(j-2)i\). If \(j=2\), the tied degree-six
terms have distinct leading monomials. The bound \(\deg p\le4\)
leaves only \(j=1,i=2\), giving \(q=x^2,\ p\propto x^4\).
Subtract a multiple of \(Q^2\) to obtain two affine mates,
which is impossible.

**Only \(v\) divides \(q\).** Write \(q=vh\), \(x\nmid h\).
If \(h\) and \(x\) are dependent, then \(q=v f(x)\) and
\(p/q^2\) depends only on \(x\), excluded by (1). Otherwise
\(h/x\) is closed by the irreducible denominator criterion,
and the earlier subtraction applies unless some \(h-cx\) is a
nonzero scalar times \(v^j\). Here \(\deg h\le3\) gives
\(j=0\) or \(1\). After subtracting \(cR\), the constant case
is already covered, so scale to \(q=v^2\). For \(u=v/x\)
the bracket again gives a Laurent polynomial \(p/a^2\) with
exponents \(-2\) through \(2\). The degree bound leaves
\(p\propto x^4\). In \(t=v/x\) this gives \(q=x^2t^2\),
excluded by (2).

The four divisibility cases exhaust \(q\).
