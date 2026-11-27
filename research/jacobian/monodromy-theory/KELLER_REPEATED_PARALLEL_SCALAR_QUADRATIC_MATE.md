# Quadratic mates of a polynomial scalar linear in a source direction

**Theorem KellerRepeatedParallelScalarQuadraticMate.**
Let \((P,Q,r)\) be a complex polynomial Keller map in \(x,v,y\)
with \(r=g(x)+A(x)v\), where \(A\in\mathbb C[x]\) is nonconstant.
Suppose \(Q=qy^2+By+C\), with \(q\ne0\) and \(q,B,C\in\mathbb C[x,v]\).
Then every prime factor of \(q\) is a component of some fibre of \(r\).
There are \(e\in\mathbb C(r)^\times\) and \(w\in\mathbb C(x)^\times\),
with the zeros and poles of \(w\) supported on the roots of \(A\)
and infinity, such that \(q=e w^2\) in \(\mathbb C(r,x)\).

More precisely, put \(k=\mathbb C(r)\) and use
\[
   Z=w\left(y+\frac{B}{2q}\right),\qquad
   K=C-\frac{B^2}{4q},\qquad Q=eZ^2+K .
\]
There exist polynomials \(\Psi\in k[Q]\) and \(F\in k[T]\), with
\(F\ne0\), such that in the function field
\[
        P=\Psi(Q)+\int_0^Z F(Q-eu^2)\,du .
                                                        \tag{1}
\]
The integral denotes a polynomial antiderivative. If \(J'=F\),
there are rational functions \(V\in\mathbb C(x)\) and \(U\in k\)
with
\[
        J(K)=V(x)+U(r),\qquad V'(x)=\frac{\kappa}{A(x)w(x)} ,
                                                        \tag{2}
\]
where \(\kappa\ne0\) is the signed constant Jacobian determinant.
In particular the generic sheet degree of \((P,Q,r)\) is
\(\deg_x K\,(2\deg F+1)\), where \(\deg_x K\) is the degree
of the rational function \(K\) over \(k\).

## Parity after centering

With the [scalar derivation](KELLER_REPEATED_PARALLEL_SCALAR_AFFINE_MATE.md#localization-by-the-affine-mate)
\(D=-A(x)\partial_x\) in \(k(x)\), the determinant condition is
\[
                         D(P)Q_y-P_yD(Q)=\kappa .
\]
Write the equation over \(k(x)\) and set \(Y=y+B/(2q)\),
so \(Q=qY^2+K\). Translation of \(y\) does not change the
bracket. The bracket of the even part of \(P\) with \(Q\)
is odd in \(Y\), and that of the odd part with \(Q\) is even.
The even part therefore commutes with \(Q\). Expressed as
a polynomial in \(Q\) over \(k(x)\), its coefficients are
\(D\)-invariant, so the even part is \(\Psi(Q)\) for \(\Psi\in k[Q]\).

The odd part has bracket \(\kappa\), so it is nonzero.
Write \(l\) for its odd \(Y\)-degree and \(h\) for
its leading coefficient. The bracket coefficient at
\(Y^{l+1}\) is \(2qDh-lhDq\), hence \(h^2/q^l\in k^\times\).
It follows that \(q=e w^2\) with \(e\in k^\times\) and
\(w\in k(x)^\times\), by taking \(w=h/q^{(l-1)/2}\).
One can use any such representation for the subsequent
calculation.

Set \(Z=wY\), so \(Q=e Z^2+K\). Denote the odd part in
\(Z\) by \(O=\sum_{i\ge0} f_i Z^{2i+1}\).
The new bracket is \(w(2e Z D(O)-DK\,O_Z)\), with \(D\)
in that expression taken at fixed \(Z\).
Comparison at \(Z=0\) gives \(-w f_0 DK=\kappa\), so \(DK\ne0\).
The top coefficient \(f_s\) is invariant and the successive
equations are
\[
                         2e Df_{i-1}=(2i+1)f_i DK
                         \qquad (1\le i\le s).
\]
Integrating downward gives polynomials in \(K\) with
coefficients in \(k\), since the kernel of \(D\) is \(k\).
In particular \(f_0=F(K)\) for a polynomial \(F\in k[T]\)
of degree \(s\), with nonzero leading coefficient.
In coordinates \((Q,Z)\), the bracket of \(O,Q\) equals
\(-w DK\) times \(\partial_Z O\) at fixed \(Q\). This
derivative equals \(F(K)=F(Q-eZ^2)\), and the odd part
vanishes at \(Z=0\). This proves (1).

## Regularity on the generic scalar fibre

The ring \(R=k[x,A^{-1}]\) contains the original coefficients
\(q,B,C\) via \(v=(r-g)/A\), and contains all coefficients of
\(P\) as a polynomial in \(y\). At a prime of \(R\) dividing \(q\),
\(w\) has positive order. If \(K\) is regular at this prime
then \(DK\) and \(F(K)\) are also regular there. This
contradicts \(-w F(K)DK=\kappa\).
Thus \(K\) has a pole at that prime. At \(y=0\), write
\(Z_0=B/(2e w)\). Since \(eZ_0^2=C-K\), \(Z_0\) has a pole.
Evaluation of (1) at \(y=0\) gives
\[
                         P|_{y=0}=\Psi(C)
                                  +\int_0^{Z_0}F(C-eu^2)\,du .
\]
The integral has a pole: as a polynomial in \(Z_0\)
it has degree \(2\deg F+1\), its leading coefficient is
a nonzero element of \(k\), and every other coefficient
is regular at the prime. But \(P|_{y=0}\) and \(\Psi(C)\)
are both regular there. This proves that \(q\) is a unit of \(R\).

Equivalently \(q\) divides a nonzero polynomial in \(r\)
in \(\mathbb C[x,v]\). Indeed \(R\) is the localization
of \(\mathbb C[x,v]\) by all nonzero polynomials in \(r\).
In this localization \(A\) is a unit since each of
its roots \(a\) satisfies \(x-a\mid r-g(a)\), and after
inverting \(A\) the inclusion and its inverse follow
from \(v=(r-g)/A\). Thus each prime of \(q\) divides
some \(r-c\). As in the
[affine-mate theorem](KELLER_REPEATED_PARALLEL_SCALAR_AFFINE_MATE.md#the-components-of-a-scalar-fibre),
factoring these \(r-c\) gives
\(q=q_r(r)q_x(x)\) with
\(q_r\in k^\times\) and \(q_x\in\mathbb C(x)^\times\)
supported at the roots of \(A\) and infinity.
The square relation already proved forces all the
orders of \(q_x\) at these roots even. Thus we may
take \(e=q_r\) and \(w\in\mathbb C(x)^\times\) as stated,
absorbing a nonzero scalar into \(e\). The parity calculation
applies in these coordinates.

The constant coefficient equation now reads
\(\partial_x J(K)=\kappa/(Aw)\). The right side is independent
of \(r\), so we can specialize \(J(K)\) at a regular constant
\(r\) value to obtain a rational primitive \(V(x)\). Subtraction
of that primitive gives the invariant \(U(r)\) in (2).

Finally \(k(x,Z)\) is the source function field over the generic
\(r\)-line. With \(Z\) fixed, \(Q-eZ^2=K(x)\), so the extension
\(k(x,Z)/k(Q,Z)\) has degree \(\deg_x K\).
With \(Q\) fixed, (1) expresses \(P\) as a polynomial in \(Z\)
of degree \(2\deg F+1\). Hence
\([k(Q,Z):k(Q,P)]=2\deg F+1\), giving the product formula.
