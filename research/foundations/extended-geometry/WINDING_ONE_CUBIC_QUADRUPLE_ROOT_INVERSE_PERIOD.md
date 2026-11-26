# An inverse-period condition for a single quadruple leading root

Work over a characteristic-zero field with the
[regular winding-one presentation](MIXED_CUBIC_SEPTIC_LEADING_FACE_REDUCTION.md)
\(f=x^{-1}(C+Bt+Dt^2+Et^3)\) and the normalized sphere functional
\(\Lambda\). Suppose \(E=e(1+\lambda x)^4\) with \(e\lambda\ne0\).
For a Laurent series at \(x=\infty\), write
\(\operatorname{CT}_{x=\infty}\) for its constant coefficient.

**Proposition.** If the positive pure moments of \(f\) vanish eventually,
then the formal branch at \(t=\infty\) of the inverse equation
\[
 f(x,T(x,w))=w^3,qquad
 T(x,w)=(x/E)^{1/3} w+T_0(x)+\sum_{j\ge1}T_j(x)w^{-j}
\]
satisfies \(\operatorname{CT}_{x=\infty} T_j(x)=0\) for every
\(j\ge1\). Here \((x/E)^{1/3}\) is expanded at \(x=\infty\), after
choosing a cube root of \(e\lambda^4\). Equivalently the constant
coefficient in \(x\) of the inverse branch is independent of \(w\).

## Arbitrary fixed deficits in the prime test

Set \(U=1+\lambda x\) and let \(n\ge2\). For \(n\not\equiv1\pmod3\),
choose arbitrarily large rational primes \(p\equiv1-n\pmod3\), and put
\(m=(p-1+n)/3\). The angular-constant polynomial in \(f^m\) has height
at most \(3m=p-1+n\). For \(p>n+1\), multiplying its interval moment
by \(p\) isolates the coefficient of height \(p-1\), since no other
height has a denominator divisible by \(p\) and \(p-1\) is even.

Put
\[
 W_n(x)=[z^n]\bigl(1+(D/E)z+(B/E)z^2+(C/E)z^3\bigr)^{(n-1)/3}.
\]
It is a rational function with possible finite poles only at \(U=0\).
In characteristic \(p\), the exponent \(m\) can be replaced by
\((n-1)/3\) inside this coefficient of degree \(n<p\), because
each contributing binomial coefficient has factorial denominator
prime to \(p\). Eventual purity therefore gives
\[
 0=[x^m]U^{4m}W_n(x)=[x^m]U^{m+n-1}W_n(x)\pmod{mathfrak p}.
\]
Here scalar factors \(e^m\) are removed at good primes, and the second
identity uses \(4m=p+m+n-1\) and
\(U^p=1+\lambda^p x^p\) in the residue field. For sufficiently large
\(p\), \(U^{m+n-1}W_n\) is a polynomial and \(m<p\), so the
\(x^p\) term cannot contribute.

Read the last coefficient at infinity and remove the unit
\(\lambda^{m+n-1}\). The finite binomial expansion needed at that
coefficient only depends on
\(m+n-1\equiv4(n-1)/3\pmod p\). Thus
\[
 [x^{-(n-1)}]\bigl(1+(\lambda x)^{-1}\bigr)^{4(n-1)/3}
                 W_n(x)=0\pmod{mathfrak p}.
\]
Over a number field the expression is fixed, so vanishing at all
sufficiently large good primes in the progression forces equality
in the coefficient field. Absolute residue degree one is unnecessary.
Over any characteristic-zero field, if a fixed such expression were
nonzero, a number-field specialization preserving it and \(e\lambda\)
would preserve every putative vanishing moment. This is impossible by
the prime test. These arguments use the
[same prime progression and specialization rules](SPHERICAL_LEADING_BINARY_FORM_MOMENT_OBSTRUCTIONS.md#3-a-leading-coefficient-detected-at-good-primes)
as the leading-coefficient test, with Dirichlet's theorem for the
nonzero residue class modulo three.

Choose \(\kappa^3=e\lambda^4\). The last coefficient, times
\(\kappa^{n-1}\), is
\[
 \operatorname{CT}_{x=\infty}[t^{-1}]f(x,t)^{(n-1)/3}.
\]
The power of \(f\) here is expanded with leading cube root
\(\kappa x t\) and power-series corrections at the two infinities.
For \(n\equiv1\pmod3\), the same coefficient vanishes automatically,
because \(f^{(n-1)/3}\) is polynomial in \(t\).

Now apply invariance of the formal residue in the variable \(t\)
under \(t=T(x,w)\). This gives
\[
 [t^{-1}] f(x,t)^{(n-1)/3}
   =[w^{-1}]w^{n-1}\partial_w T(x,w)
   =-(n-1)T_{n-1}(x).
\]
Both residue signs are the same at the two infinities. Taking the
constant coefficient in \(x\) proves the proposition. The leading
coefficient \((x/E)^{1/3}\) has order one at \(x=\infty\), so its
constant coefficient already vanishes.

## The remaining implication in the logarithmic boundary case

Under the [boundary divisibilities](WINDING_ONE_CUBIC_LOGARITHMIC_BOUNDARY.md#a-quadruple-root)
write \(D=U^3 d\), \(B=U^2 b\), and \(C=U c\) with polynomial
\(d,b,c\). Setting \(v=Ut\) gives
\[
 f=\frac{U}{x}\bigl(e v^3+d(x)v^2+b(x)v+c(x)\bigr).
\]
In the residue representation one may take
\(s=(\lambda x/U)^{1/3}\), with \(s=1\) at \(x=\infty\).
If \(a=(x/E)^{1/3}\), then, for compatible cube roots,
\(a\,dx/x=3(e\lambda)^{-1/3} ds\). Thus the condition on the inverse
is the vanishing at \(s=1\) of the local residues of the negative
\(w\)-coefficients of \(T/a\). The two-variable formal inverse depends
on the value parameter \(w\); a fixed rational moment problem in
\(s\) applies only under additional composition hypotheses.

For example, if the two depressed coefficients of the polynomial
in \(v\) above were constant in \(x\), that polynomial would be
\(e(v-S(x))^3+p_0(v-S(x))+q_0\) with \(S\) polynomial. Its sphere
degree would then be \(\max(6,3\deg S)\), by the
[winding-one degree formula](MIXED_CUBIC_ZERO_CENTER_AND_CONDITIONAL_OFFCENTER.md).
This explains the relevance of the inverse-period condition to the
septic branch. The implication from the residue condition to constant
depressed coefficients is unresolved.
