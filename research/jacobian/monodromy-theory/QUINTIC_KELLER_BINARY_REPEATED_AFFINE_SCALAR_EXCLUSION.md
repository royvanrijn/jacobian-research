# The binary multiplicity-two affine scalar

**Theorem QuinticKellerBinaryRepeatedAffineScalarExclusion.** There is
no three-variable Keller triple of characteristic zero with
\(\deg P,\deg Q\le5\), with homogeneous degree-five terms of \(P,Q\)
in \(x,y\), and with third component
\[
                    R=\mu+\lambda x+x^2(z+h(x,y)),
             \qquad \lambda\ne0,\qquad \deg h\le2 ,
\]
where the coefficient of \(y^2\) in \(h\) is nonzero. Here \(P_5\) or \(Q_5\) may be zero. In particular, the repeated
affine-source case of the
[binary multiplicity-two quartic restriction](QUINTIC_KELLER_BINARY_QUARTIC_REPEATED_CUBIC_VANISHING.md#excluding-the-planar-multiplicity-two-case)
cannot occur.

## Degree bounds after the shear

Work over \(\mathbb C\) by embedding the finitely generated
coefficient field. An affine change in \(y\) completes the square
in \(h\).
Absorb the part of \(h\) affine in the source variables into \(z\),
and normalize \(R\) and \(z\) to make \(\lambda=1,\mu=0\).
The resulting form is
\[
 R=x+x^2(z+g(x,y)+b x^2),
 \quad g=a y^2,\quad a\ne0 .
\]
All affine changes preserve the binary
degree-five parts. Set \(V=z+g(x,y)\), leaving \(y\) free.
The polynomial \(r=x+x^2 V+b x^4\) is now the third component.

In \(\mathbb C[x,V,y]\) the determinant is, up to sign,
\(D(P)Q_y-P_y D(Q)=\kappa\ne0\), where
\(D=-x^2\partial_x+(1+2xV+4b x^3)\partial_V\).
The field \(\mathbb C(x,V)=\mathbb C(r,x)\), with
\(D=-x^2\partial_x\) at fixed \(r\), has constants \(\mathbb C(r)\).
Its polynomial constants are \(\mathbb C[r]\) by a
univariate Bézout identity after substitution.
There is no polynomial \(D\)-primitive of a nonzero constant:
such a primitive would differ from a multiple of \(1/x\) by
a rational function of \(r\), and a reduced denominator after
substitution would divide \(x\), although \(\deg r\ge3\).
As in the
[exceptional scalar argument](QUINTIC_KELLER_CUBIC_COMPONENT.md#the-exceptional-scalar-through-degree-five),
these facts forbid a component independent of \(y\), a nonzero
constant leading coefficient in \(y\), or two components
affine in \(y\). They are the degree-independent portion of
[Wang–Yan, Proposition 3.17](https://arxiv.org/pdf/2610.36495v1).

If the resulting component has coefficient
\(p_m\) of \(y^m\), then
\[
                \deg_{x,V}p_m\le\max(5-m,\,4-\lceil m/2\rceil).
 \tag{1}
\]
For terms from the binary quintic part this follows directly.
For a monomial \(x^i y^j z^k\) of degree at most four, a term in
its substitution with \(l\) factors \(-a y^2\) has
\(m=j+2l\) and coefficient degree \(i+k-l\le4-m+l\).
Since \(l\le\lfloor m/2\rfloor\), this gives (1).
The same bounds hold for \(Q\).
In particular the nonconstant leading coefficient restriction
makes both \(y\)-degrees at most six.

## Unequal leading coefficients

Let \(m\ge n\) be the two \(y\)-degrees after the shear, and
\(p_m,q_n\) their leading coefficients. The highest determinant
coefficient gives
\[
                        p_m^n/q_n^m \in\mathbb C(r).
 \tag{2}
\]
Each fibre \(r-c\) is squarefree by submersion and has degree
three or four. If a rational function \(H(r)\) is nonconstant,
every finite zero or pole of \(H\) forces a full fibre factor
in the corresponding numerator or denominator of (2).
This uses coprime univariate numerator and denominator
after substitution. Also each fibre has a factor involving
\(V\) because it has degree one in \(V\).

If \(m=n\), the ratio \(p_m/q_n\) itself is in \(\mathbb C(r)\),
as follows by differentiating (2). For \(m\ge3\) the degree
bounds are at most two, so the ratio is constant.
For \(m=n=2\), a nonconstant ratio requires one coefficient
to be proportional to \(r-c\), since both degrees are at most
three. The other is then also in \(\mathbb C[r]\), of \(r\)-degree
at most one. If these leading coefficients are not proportional,
a constant linear combination has nonzero constant leading
coefficient in \(y\), which is forbidden. Thus for \(m=n\ge2\)
one can cancel a leading term by constant target subtraction.
This preserves the original degree constraints. The case \(m=n=1\)
is already forbidden. Reorder to obtain \(m>n\).

For \(m>n\ge3\), both coefficients have degrees at most two,
so (2) is constant. The degree comparison \(n\deg p_m=m\deg q_n\)
is impossible when \(m=5,6\), because \(\deg p_m\le1\).
For \(m=4,n=3\), it would require \(4\mid \deg p_4\),
again impossible.

For \(m\ge3,n=2\), a fibre factor in \(q_2\) would make it a
scalar multiple of that factor by its degree bound three.
Equation (2) would then make \(p_m\) invariant and hence
constant by its smaller bound. Otherwise (2) is constant.
The degree comparison shows the only candidate is
\(m=4,n=2,\deg p_m=2,\deg q_2=1\).

For \(m\ge3,n=1\), if \(q_1\) contains a fibre factor, write
\(q_1=(r-c)\ell\) with \(\deg\ell\le1\).
The ratio \(p_m/\ell^m\) is invariant, and is constant because
neither part has a full fibre factor. This is impossible
since \(p_m\) is nonconstant of degree at most two.
If \(q_1\) has no fibre factor, (2) is already constant and
forces \(\deg p_m=m\deg q_1\), also impossible.

## The two unequal-degree cancellations

No such bracket can have a \(y\)-degree-two component with
linear leading coefficient \(q_2\). Indeed the analogue of (2)
holds for a mate of any \(y\)-degree \(j\). Every finite zero
or pole of the rational ratio \(p_j^2/q_2^j\) has even order:
test an irreducible factor involving \(V\) of the corresponding
squarefree fibre. Such a factor does not divide \(q_2\).
In fact if \(q_2\) involved \(V\) and divided \(r-c\), the
other factor would depend only on \(x\), and comparison of the
\(V\)-coefficient forces that factor to be proportional to \(x^2\).
The term of degree one in \(x\) of \(r-c\) rules this out.
Thus the rational ratio is a square in \(\mathbb C(r)\), and
unique factorization forces \(j\) even. Then \(p_j/q_2^{j/2}\)
is polynomial in \(r\), since a reduced denominator would
require a full fibre factor in \(q_2\). Subtract that multiple
of \(Q^{j/2}\) from the mate. The coefficient in \(r\) is
\(D\)-constant and independent of \(y\), so the bracket is
preserved. Descending this way would yield an odd-degree or
independent mate, both impossible. This excludes the candidate
\((m,n)=(4,2)\).

Finally let \(m=2,n=1\). A constant ratio in (2) allows
subtraction of a multiple of \(Q^2\), giving an affine
or independent mate. Otherwise one of \(p_2,q_1\) contains
a full fibre factor. If \(p_2\) does, its degree bound makes
it a scalar multiple of that fibre, and (2) makes \(q_1\)
invariant and nonconstant, so \(q_1\) contains a full fibre
factor in either case. Write \(P=p_2 y^2+c_1y+c_0\) and
\(Q=q_1 y+b_0\). Reducing the bracket modulo the fibre
factor of \(q_1\) gives
\(-\overline{P_y}\,\overline{D b_0}=\kappa\).
On every component of this reduced fibre the identity forces
\(\overline{p_2}=0\). Thus \(r\) must have degree three and
\(p_2,q_1\) are scalar multiples of the same fibre polynomial,
using invariance of \(q_1\) and its bound four.
Writing \(p_2=\delta q_1\), the middle coefficient gives
\(c_1-2\delta b_0=L(r)\), and the constant equation is
\[
                 D(q_1 c_0-\delta b_0^2-L(r)b_0)=\kappa .
\]
This is the forbidden polynomial primitive. That finishes
the case that \(h\) has nonzero \(y^2\) coefficient.
