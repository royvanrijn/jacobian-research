# The free-degree pair \((9,6)\) in the sextic tangent-conic scalar

**Theorem SexticTangentConicNineSixExclusion.** A
three-variable Keller map of ordinary degree at most six
in characteristic zero with quartic component
\(R=g(x)+x^2(xz+y^2)\) cannot have free-variable degrees
\((9,6)\) after the localized substitution \(z=v-y^2/x\).

Work over \(\mathbb C\) and use the
[tangent-conic leading restriction](SEXTIC_KELLER_TANGENT_CONIC_HIGH_MATES.md).
The leading free-variable coefficients of \(P,Q\) can
be normalized to \(x^{-3},x^{-2}\). Write \(r=g(x)+x^3v\)
and \(k=\mathbb C(r)\). Adjoin \(t\) with \(x=t^3\).
If \(B_5=[y^5]Q\), take \(Z=t^{-1}(y+B_5/(6x^{-2}))\).
In this variable both mates are monic and \(Q\) is centered,
meaning \([Z^5]Q=0\). Since \(\partial_y Z=t^{-1}\),
the bracket at fixed \(Z\) becomes
\[
              D(P)Q_Z-P_ZD(Q)=\kappa t,\qquad
              D=-\tfrac13 t^7\partial_t,\quad \kappa\ne0 .
 \tag{1}
\]
Over \(k(t)\), its constant field is \(k\).
The automorphism \(\sigma:t\mapsto\zeta t\) for a primitive
cube root of unity \(\zeta\) fixes the original mates and
sends \(Z\) to \(\zeta^{-1}Z\).

## The formal square-root truncation

Use \(W=Q^{1/6}=Z+O(Z^{-1})\) as a formal Laurent
variable. This is the formal-root framework of
[Wang and Yan, Proposition 3.17](https://arxiv.org/pdf/2610.36495v1).
Equation (1) gives \(D(P)|_W=\kappa t/Q_Z\), whose
order is \(-5\) at \(W=\infty\). Every nonnegative
\(W\)-coefficient of \(P\) is therefore invariant.
Since \(\sigma(W)=\zeta^{-1}W\), the only possible
nonnegative powers are nine, six, three and zero.
Consequently \(P\) is the polynomial part in \(Z\)
of \(Q^{3/2}+\lambda Q+\mu Q^{1/2}+\nu\), with
\(\lambda,\mu,\nu\in k\).

Subtract \(\lambda Q+\nu\) from \(P\). Let
\(S=Z^3+aZ+b\) be the polynomial part of \(Q^{1/2}\),
and write \(Q=S^2+U\) with \(\deg_Z U\le2\).
For \(U=cZ^2+dZ+e\), the polynomial part of
\(U^2/S\) is \(c^2Z+2cd\); all terms beyond
it in the binomial expansion of \(Q^{3/2}\) have
negative \(Z\)-degree. Adding \(2\mu/3\) to \(Q\)
only shifts the constant term of \(U\), leaving
the bracket unchanged. Relabel the result as
\[
        Q=S^2+U,\qquad
        P=S^3+\tfrac32 S U+\tfrac38(c^2Z+2cd).
 \tag{2}
\]
In an invariant polynomial, the coefficient of \(Z^j\)
has character \(\zeta^j\). Here \(S\) and \(U\) are
invariant as polynomials in the function field.
Thus \(a,d\) have character \(\zeta\), \(c\) character
\(\zeta^2\), and \(b,e\) are invariant.

## The coefficient obstruction

In the bracket of (2), the coefficients at \(Z^4,Z^3\)
are the \(D\)-derivatives of \(-9I/4,-9J/4\), where
\[
            I=-ac^2+2ce+d^2,\qquad J=-2acd-bc^2+2de .
\]
The [exact replay](verify_sextic_tangent_conic_nine_six.py)
checks these and the remaining identities used below.
Since \(I,J\in k\) by (1) but their characters are
\(\zeta^2,\zeta\), both vanish.
If \(c=0\), then \(d^2=0\) follows from \(I=0\).
Now \(Q=S^2+e\) and
\([P,Q]=-3e(De)S_Z/2\) from (2).
As \(S_Z\) has degree two, this cannot be a nonzero
scalar in \(Z\).

Suppose \(c\ne0\). Put \(u=d/c\), \(v_0=e/c\).
The equations \(I=J=0\) give
\[
                      a=2v_0+u^2,\qquad b=-2u(v_0+u^2).
\]
The coefficient at \(Z\) is now
\(9c^2(cDu+3uDc)/8\), so \(c^3u\in k\).
But \(u\) has character \(\zeta^2\), hence \(u=0\).
In particular \(d=b=0,\ a=2e/c\). The coefficient at
\(Z^2\) gives \(e^2-c^3/6\in k\), and the constant
coefficient equation is
\[
                  -\tfrac38(c^2+8e^2/c)De=\kappa t .
 \tag{3}
\]
If the invariant \(e^2-c^3/6\) is nonzero, the pair
\((c,e)\) gives a rational map from the line over
\(\overline k\) to the smooth elliptic curve
defined by this equation. Such a map is constant
by Riemann–Hurwitz. This would give \(De=0\).

If the invariant vanishes, put \(s=e/c\).
Then \(c=6s^2,\ e=6s^3\), and the left side of
(3) is a nonzero scalar multiple of \(D(s^7)\).
Since \(D=-t^7\partial_t/3\), equation (3) makes
\(s^7\) a nonzero scalar multiple of \(t^{-5}\)
plus an invariant. Its pole order five cannot be
that of a seventh power in \(k(t)\).
This excludes the last case.
