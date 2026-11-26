# A matched tanh and Newton polygon gate for torsion common factors

Consider an ordinary two-monomial profile with exactly one
private unit row on each side and arbitrary shared rows.
Put
\[
                         A=\sum_{\rm shared\ rows}
                                    \frac{(b_i-a_i)^2}{\min(a_i,b_i)}.
\]
This is determined by the exact matched tanh data. For a shared
Gauss factor at adjacent orders \(n,N=n+1\), put
\[
   v_i=\left(\frac{n}{T_i'}-\frac{n-1}{T_i}\right)^{-1}.
\]
The [exact matched grid](MULTIVARIATE_MATCHED_TANH_INTEGER_MULTIPLICITY.md#recovering-the-gauss-grid)
has \(d\) factors with \(v_i=d/(2\ell+d)\) for a shared row
of smaller exponent \(\ell\) and gap \(d\). Therefore
\(A=\sum_i 2v_i/(1-v_i)\) over the shared factors.

Fix an adjacent pair \(n,N=n+1\), \(n\ge2\).
Suppose an exact argument bounds the rational common factor
degree by an integer \(D\ge2\). At each prime let \(S_p\)
be the common set of root valuations of the two Newton polygons.
Assume the following two conditions.

- At each prime, any two distinct values in \(S_p\) have reduced
  denominators whose sum exceeds \(D\).
- The shared curvature estimate satisfies
  \[
                         (D-2)A<\frac{\log 2}{4D}.            \tag{1}
  \]

Then every irreducible common factor would have one of the forms
in the table. The column on the right gives a necessary valuation
at the stated prime for its roots. If all rows of the table are
excluded, the actual moments are coprime.

| Degree | Positive root arguments | Lower order \(n\) | Necessary valuation |
|---:|---|---|---|
| \(2\) | \(\pi/4\) | \(n\equiv5\pmod 8\) | Half-integer outside \(\mathbb Z\) at \(p=2\) |
| \(2\) | \(\pi/6\) | \(n\equiv7,8,9\pmod {12}\) | Half-integer outside \(\mathbb Z\) at \(p=3\) |
| \(4\) | \(\pi/12,5\pi/12\) | \(n\equiv17\pmod {24}\) | Half-integer outside \(\mathbb Z\) at \(p=2\) |

For instance \(D=4\) permits \(A<\log2/32\).
The first hypothesis needs checking only at primes dividing
the endpoint products of both integer moment polynomials.
The [Newton edge width bound and delayed marks](MULTIVARIATE_QUADRATIC_NEWTON_TRACE_SIEVE.md)
give ways to obtain the degree bound without a root calculation.
For \(D=2\), the two hypotheses reduce to the integral shared
valuation criterion of that proof.

The parameters and valuation tests concern the original support.
They can be combined with the
[absolute angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md)
for all-degree coprimality. The first pair is coprime by positivity
of the quadratic over the reals. The gate itself proves physical
coprimality and uses a matched curvature cap together with local
arithmetic.

## The matched root annulus and algebraic units

The two unit private rows cancel the binomial coefficient up
to a constant. Endpoint-normalize the moments in degree \(m\)
so that the constant and leading coefficients are one.
Extend the log of their coefficient weights to a function \(f(x)\)
on \(0\le x\le m\) using Gamma functions. For a shared row, the
logarithmic second derivative is bounded above by
\(d^2/(m\ell)\), by the positive trigamma series and its integral
comparison. Hence \(0\le f''(x)\le A/m\) and \(f(0)=f(m)=0\).
Integration gives \(-A/2\le f'(0)\le f'(m)\le A/2\).
The ratios of neighboring coefficients lie in
\([e^{-A/2},e^{A/2}]\). The Eneström–Kakeya annulus therefore
bounds the difference between the largest and smallest
log root moduli by \(A\), in the normalized or original variable.

Let \(H\) be a monic irreducible rational common factor of
degree \(d\le D\). All its roots are nonreal by real positivity
of the even moment, so \(d\) is even. At each prime a Newton
polygon edge of \(H\) with slope of reduced denominator \(b\)
has width divisible by \(b\), since its ordinates are integral.
With two distinct root valuations the sum of those denominators
would be at most \(d\). The first hypothesis therefore puts
every root of \(H\) at the same valuation at each prime.

Write \(c=\prod_{H(\alpha_i)=0}\alpha_i>0\), and let
\(\beta=\alpha^d/c\) for one root \(\alpha\).
At every prime and every embedding, \(v_p(\beta)=0\),
so \(\beta\) is an algebraic unit. In fact its minimal
polynomial has integral coefficients at every prime
and a unit constant coefficient. Its degree is at most \(d\).
For a root \(\alpha_i\) of maximal modulus,
\[
       \log|\alpha_i^d/c|=
            \sum_j(\log|\alpha_i|-\log|\alpha_j|)
                       \le (d-2)A .
\]
Both \(\alpha_i\) and its complex conjugate have that
maximal modulus and give zero terms. This bounds the house
of \(\beta\), the greatest modulus of its conjugates.

The [Schinzel–Zassenhaus theorem of Dimitrov](https://annals.math.princeton.edu/articles/22975)
implies that a nonzero algebraic integer of degree \(e\)
which is not a root of unity has house at least
\(2^{1/(4e)}\), by applying the stated polynomial theorem
to the reversed monic minimal polynomial. Condition (1)
therefore implies that \(\beta\) is a root of unity.
All \(\alpha_i\) have the same modulus \(R>0\), and every
quotient \(\alpha_i/\alpha_j\) is a root of unity.

## The cyclotomic orders allowed by the angle boxes

All roots of \(H\) have arguments in the open right half-plane:
the [both-private angle boxes](MULTIVARIATE_PRIVATE_UNIT_ASYMMETRIC_ANGLE_BOXES.md)
put each common upper root in the overlap of the two consecutive
boxes with the same label, and the lower roots are their conjugates.
Their arguments are all distinct. Thus \(\operatorname{Tr}(\alpha)\)
is nonzero, since every root has positive real part.

Let \(L\) be the splitting field of \(H\) and \(L_0\) the field
generated by all root quotients. This is a cyclotomic field,
since the generators are roots of unity. An automorphism of
\(L\) over \(L_0\) acts on the roots by a common multiplier,
which must be one because their sum is a nonzero rational.
Hence \(L=L_0\) and the Galois group is abelian. Its transitive
action on the roots is regular, so \(\mathbb Q(\alpha)=L\).

Put \(q=\alpha/\bar\alpha\), of order \(M\ge3\). Complex
conjugation commutes with the Galois group. Its action on
the conjugates of \(\alpha\) therefore shows that the
conjugates of \(q\) are the ratios of those roots to their
complex conjugates. These are \(d\) distinct values: doubling
an angle is injective modulo \(2\pi\) on the open right
half-plane. Consequently \(\varphi(M)=d\) and
\(L=\mathbb Q(q)\). The positive angles of \(H\) are
exactly \(j\pi/M\) for \(0<j<M/2\) coprime to \(M\).

For any such \(j\), the overlap condition
\((2k-1)\pi/n<\theta<2k\pi/(n+2)\) is equivalent to
\[
                   M<n j\bmod 2M<2M-2j,                     \tag{2}
\]
where the residue is in \(\{0,\ldots,2M-1\}\).
For odd \(M\), use \(j=(M-1)/2\); the interval becomes
\((M,M+1)\), which contains no integer.

For \(M=4m\), use \(j=2m-1\). Its required residue is
\(4m+1\) modulo \(8m\), so \(n\equiv2m-1\pmod {8m}\)
when \(m\) is even and \(n\equiv6m-1\pmod {8m}\)
when \(m\) is odd. This follows from
\((2m-1)^2=1+4m(m-1)\). Even \(m\) is excluded by
\(j=1\). For odd \(m\), substitution into (2) allows
exactly the primitive \(j<2m\) with \(j\equiv1\pmod4\).
If \(m\ge5\) and \(m\equiv1\pmod4\), the primitive
\(j=m-2\) violates this. If \(m\ge7\) and \(m\equiv3\pmod4\),
use \(j=m-4\). Hence only \(M=4,12\) remain, with
\(n\equiv5\pmod8\) and \(n\equiv17\pmod {24}\), respectively.

Finally, let \(M=2m\) with \(m\) odd. The case \(m=3\)
gives \(M=6\) and \(n\equiv7,8,9\pmod {12}\).
For \(m\ge5\), the test at \(j=1\) writes the residue
of \(n\) as \(2m+x\), \(1\le x\le2m-3\).
The primitive index \(j=m-2\) requires residue \(2m+r\)
with \(1\le r\le3\). At \(j=m-4\) the residue is
\(2m+r-2x\) modulo \(4m\), and (2) requires
\(r-2x\bmod4m\in\{1,\ldots,7\}\).
As \(7-4m\le r-2x\le1\), this forces either
\((x,r)=(1,3)\) or \((2m-3,1)\).
Substitution into the \(j=m-2\) test forces \(m=5\)
in both cases. Thus \(M=10\) is the only additional order.

## The remaining norm and valuation restrictions

Complex conjugation fixes \(R^2=\alpha\bar\alpha\).
It is rational, since all conjugates have the same
modulus. For \(M=10\) or \(12\),
\[
                b=\frac{\alpha}{1+q}
                        \in L^+,\qquad
         R^2=b^2(2+q+q^{-1}),
\]
where \(L^+\) is the real quadratic subfield fixed by
complex conjugation. At \(M=10\),
\(\mathrm{Nm}_{L/\mathbb Q}(1+q)=\Phi_{10}(-1)=5\).
It follows that
\[
    \mathrm{Nm}_{L/\mathbb Q}(\alpha)
            =5\,\mathrm{Nm}_{L^+/\mathbb Q}(b)^2
                         =(R^2)^2 ,
\]
which is impossible over the rationals. This excludes order ten.

At \(M=12\) the real field is \(\mathbb Q(\sqrt3)\), and
\(2+q+q^{-1}=2\pm\sqrt3=(1\pm\sqrt3)^2/2\).
The element \(y=(1\pm\sqrt3)b\) has \(y^2=2R^2\in\mathbb Q\).
Writing \(y=a+b_0\sqrt3\) shows \(ab_0=0\), so \(y^2\)
is a rational square or three times a rational square.
In either case \(v_2(R^2)\) is odd. As the roots of \(H\)
have uniform valuation at \(2\), their shared valuation
is a half-integer outside \(\mathbb Z\).

At \(M=4,6\) the factor is quadratic. For its rational
trace \(u\), one has
\(R^2=u^2/(2+q+q^{-1})\). The denominator is \(2\) or \(3\),
respectively, giving the other two half-integer tests.
These are all the possibilities in the table.
