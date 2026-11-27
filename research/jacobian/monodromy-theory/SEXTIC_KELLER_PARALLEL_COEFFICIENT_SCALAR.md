# Sextic Keller maps and scalar components with parallel coefficients

**Theorem SexticKellerParallelCoefficientScalarExclusion.**
Let \(P,Q,r\) be polynomials over a characteristic-zero field,
in affine source coordinates \(x,v,y\), of ordinary total degree
at most six. If
\[
                        r=g(x)+A(x)v
\]
with \(A\) nonconstant, then \(\det J(P,Q,r)\) cannot be a nonzero
constant. Here \(g\) and \(A\) are otherwise arbitrary univariate
polynomials.

**Corollary SepticKellerParallelCoefficientScalarExclusion.** The
same exclusion holds for mates of ordinary degree at most seven
when \(5\le\deg r\le7\).

Embed the finitely generated coefficient field in \(\mathbb C\).
The Keller condition would make \(r\) a submersion. For every
root \(a\) of \(A\), this requires \(A'(a)=0\) and
\(g'(a)\ne0\), by the formula \(r_x=g'(x)+A'(x)v\).
The fibre \(r=g(a)\) is reducible: \(x-a\) divides its defining
polynomial and the quotient has positive \(v\)-degree.
If \(\deg r\le3\), the
[sextic cubic-target theorem](SEXTIC_KELLER_CUBIC_COMPONENT.md)
would make \(r\) a coordinate of an automorphism, whose fibres
are irreducible. For the rest of the proof, let \(d_0=6\) for the
theorem and \(d_0=7\) for the corollary. In both cases the
ordinary total degree is at most \(d_0\) and
\(\deg r\ge d_0-2\).

## Reduction to quadratic mates

Put \(B=\mathbb C[x,v]\) and
\[
 D=-A(x)\partial_x+(g'(x)+A'(x)v)\partial_v
   =-A(x)\partial_x|_r .
\]
Then \(D\) fixes \(r\), has rational kernel \(\mathbb C(r)\),
and the Keller equation, up to orientation, is
\[
                         [P,Q]=D(P)Q_y-P_yD(Q)=\kappa\ne0 .
 \tag{1}
\]
No polynomial in \(B\) has \(D\)-derivative a nonzero
constant. Indeed, expanding such a polynomial over
\(\mathbb C(x)[r]\) would give \(F(x)+H(r)\) with \(H\)
polynomial in \(r\). Polynomiality in \(x,v\) then gives
\(F\in\mathbb C[x]\), and \(-AF'\) cannot be a nonzero constant.

A leading coefficient \(p\) of a mate of \(y\)-degree \(m\)
cannot be a nonzero constant. More generally, after adjoining
an \(m\)-th root of \(p\) and inverting \(p\), the binomial
series \(S=P^{-1/m}\) has leading term \(p^{-1/m}y^{-1}\).
The identity
\[
                   D(QP_yS)-\partial_y(QD(P)S)=-\kappa S
 \tag{2}
\]
shows that \(p^{-1/m}\) has a primitive for \(D\), up to a
nonzero scalar. Indeed, take the coefficient of \(y^{-1}\);
the formal \(y\)-derivative contributes zero.
For constant \(p\) this would give a polynomial primitive
of a nonzero constant. This is the formal residue argument
of [Wang and Yan, Proposition 3.17](https://arxiv.org/pdf/2610.36495v1).

By the [affine-mate exclusion](KELLER_REPEATED_PARALLEL_SCALAR_AFFINE_MATE.md),
both degrees in \(y\) are at least two. Write these degrees
as \(m\ge n\), with respective leading coefficients \(p,q\).
Thus \(m\le d_0-1\), \(\deg p\le d_0-m\), and
\(\deg q\le d_0-n\).
The highest coefficient of (1) gives
\[
                              p^n/q^m\in\mathbb C(r).
 \tag{3}
\]
A nonconstant reduced rational function on the right has a
finite zero or pole. This forces a full fibre polynomial \(r-c\)
to divide the corresponding numerator or denominator
coefficient, by substitution of the Bézout identity and
squarefreeness of \(r-c\) from submersion.
If \(m=n\), relative algebraic closedness of \(\mathbb C(r)\)
in \(\mathbb C(r,x)\) gives \(p/q\in\mathbb C(r)\).
For \(m=n\ge3\) neither coefficient can contain a full fibre
factor because its degree is at most \(d_0-3<\deg r\).
For \(m=n=2\), a nonconstant ratio would force such a factor
in one coefficient. Its degree bound of \(d_0-2\) makes it a
scalar multiple of \(r-c\), so the other is also invariant.
Both are then linear polynomials in \(r\) by their degree
bounds. A constant target combination would give a nonzero
constant leading coefficient, impossible by (2).
Thus equal degrees always give constant ratio, and a constant
target shear makes them unequal without increasing ordinary
degree. Reduction to an affine mate is already excluded.

For \(m>n\ge2\), the coefficient \(p\) has degree at most
\(d_0-3\). A full fibre factor in \(q\) would force \(n=2\) and
\(q\) a scalar multiple of \(r-c\), since \(\deg q\le d_0-2\).
Then (3) makes \(p\) invariant, contrary to its being
nonconstant of degree less than \(\deg r\).
Hence (3) is constant. If \(d=\gcd(m,n)\), factorization
gives \(p=\lambda\ell^{m/d}\), \(q=\mu\ell^{n/d}\) with
nonzero scalars \(\lambda,\mu\).
The inequalities
\[
                     (m/d)\deg\ell\le d_0-m,\qquad
                     (n/d)\deg\ell\le d_0-n
\]
leave only \((m,n)=(3,2),(4,2)\), and in both \(\ell\)
is nonconstant affine linear.

## The root of the leading coefficient

The [quadratic-mate theorem](KELLER_REPEATED_PARALLEL_SCALAR_QUADRATIC_MATE.md)
makes every prime factor of \(q\) a component of some fibre of \(r\).
In either possible degree pair, \(q\) is a power of the affine
linear \(\ell\). If \(\ell=\alpha x+\beta v+\gamma\) with
\(\beta\ne0\), the restriction of \(r\) to \(\ell=0\) is
\(g(x)-A(x)(\alpha x+\gamma)/\beta\). Its derivative at every
root of \(A\) is \(g'(a)\ne0\), which contradicts constancy
on a fibre component.
Hence \(\ell\) is univariate, and its root \(a\) must also be
a root of \(A\), since otherwise \(r\) varies with \(v\)
along \(x=a\).

Translate the shared root to zero. In either pair,
\(q=q_0 x^e\) with \(q_0\ne0\), where \(e=2\) for
\((3,2)\) and \(e=1\) for \((4,2)\).

## Centered bracket at the vertical component

In both pairs write \(Y=y+b/(2q)\), \(K=C-b^2/(4q)\).
The odd part of \(P\) is \(fY+hY^3\) and has bracket
\(\kappa\) with \(Q=qY^2+K\), by parity.
The even part commutes with \(Q\). It is a polynomial in
\(Q\) over \(B[q^{-1}]\), and its coefficients have zero
\(D\)-derivative, by expressing its bracket in \(Q\).
Since \(B[q^{-1}]\cap\mathbb C(r)=\mathbb C[r]\), the
even part is polynomial in \(x,v,y\). The intersection
holds because a denominator factor \(r-c\) in a reduced
fraction would divide a power of the univariate \(q\),
which is impossible. Thus the odd part is polynomial
in \(x,v,y\) as well, and \(h\in B\) as its \(y^3\)
coefficient. If \((m,n)=(3,2)\), then \(h=p\)
is a nonzero scalar multiple of \(x^3\).

The constant and quadratic coefficients of its bracket
are
\[
                      -f DK=\kappa,\qquad
                      2q Df-fDq-3h DK=0 .
 \tag{4}
\]
So \(f=f_0 x^k\) is a unit of \(B[q^{-1}]\), with
\(f_0\ne0,\ k\in\mathbb Z\), and \(DK\) is univariate.
Expanding \(K\) over \(\mathbb C(x)[r]\) gives
\(K=T(x)+H(r)\) with \(H\) polynomial in \(r\).
At any other root of \(A\), the rational function \(T\)
is regular because \(K=C-b^2/(4q)\) and \(H(r)\)
are regular there. Then \(DK=-A T'\) vanishes there,
while \(f\) is a unit there, contradicting (4).
Thus \(A=a_0 x^M\) with \(a_0\ne0\) and \(M\ge2\).

If \(T\) is regular at zero, then \(-k\ge M\) from
(4). But the formula
\[
                         h=-\frac{f(2qDf-fDq)}{3\kappa}
 \tag{5}
\]
then gives \(h\) a pole, because its nonzero leading
coefficient has a factor \(2k-e\) and its \(x\)-exponent
is \(2k+e+M-1<0\). Thus \(K=T+H(r)\) has a pole.
From \(K=C-b^2/(4q)\) and \(e\in\{1,2\}\), that
pole has order \(e\) and \(b_0=b(0,v)\ne0\).
The first equation of (4) gives
\[
                  k=1+e-M,\qquad
                  b_0^2=\frac{4q_0\kappa}{e\,a_0 f_0}.
 \tag{6}
\]
In particular \(b_0\) is a nonzero scalar.

For \(e=2\), formula (5) and \(h=p\) of order three
force \(M=4\) and \(k=-1\). The coefficient of \(x^{-3}\)
in \(fY+hY^3\) at \(y=0\) is
\(f_0 b_0/(6q_0)\) by (5) and (6).
This contradicts its polynomiality.
For \(e=1\), formula (5) gives \(h\) of order \(4-M\),
so \(M\le4\). If \(M>2\), the \(y\)-coefficient of
the odd part,
\[
                         f+\frac{3h b^2}{4q^2},
\]
has order \(2-M\) with leading coefficient
\((4-2M)f_0\), a pole. If \(M=2\), the odd part at
\(y=0\) instead has nonzero \(x^{-1}\)-coefficient
\(f_0 b_0/(3q_0)\). Each possibility contradicts
polynomiality, proving the theorem.

The [exact coefficient replay](verify_sextic_parallel_coefficient.py)
checks the finite degree comparison, (4), and the pole
coefficients used in the last step.
