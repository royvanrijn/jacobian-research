# Fixed strip moments with supersingular cyclic covers of unbounded conductor

Let \(k\) be a field of characteristic zero. Write the affine sphere
as \(2xy+t^2=1\) and substitute \(y=(1-t^2)/(2x)\). For a sphere
polynomial \(f\) this gives a literal Laurent polynomial in \(x,t\).
The sphere functional \(\mathcal A\) takes the constant term in \(x\)
and applies the normalized interval integral in \(t\). Thus an angular
constant monomial \(t^j\) has moment zero for odd \(j\) and
\(1/(j+1)\) for even \(j\).
The support bounds are those of the
[fixed strip construction](CM_NEWTON_DIRECT_STRIP_MOMENTS.md), with
the simpler binomial top face replacing the four-branch cover.

**Theorem.** Let \(s\ge3\) be odd and \(e\ge2\), with
\(\gcd(s,n)=1\) and \(se<n\). Choose \(A,c\in k^\times\) and
\(P(z)=A(z-c)^e\). Consider any fixed sphere polynomial \(f\) whose
literal support obeys
\[
 j+2i\le n,\qquad (n-1)j+2ni\ge n(n-3)
\]
and whose angular polynomial on \(j+2i=n\) is \(z^{-1}P(z)\).
Then its least sphere-polynomial degree is \(n+1\), and it is an
absolute nonpower. For all sufficiently large rational primes
\(p\equiv s^{-1}\pmod n\), put \(m=(sp-1)/n\). Then
\[
                  \mathcal A(f^m)\ne0 .
\]
The progression has natural density \(1/\varphi(n)\).

If, in addition, \(s^h\equiv-1\pmod n\) for some \(h\ge1\), the same
certifying primes have supersingular geometric reduction of the
entire Jacobian of the connected auxiliary cover
\[
                         W^n=z^{n-1}P(z).
\]
At these primes a single coherent Frobenius step from the primitive
character \(s\) to character one is nonzero. In particular,
\(n=s^h+1\), \(e=2\), and \(h\ge2\) give unbounded actual deck
conductors, with all compatible lower-weight terms allowed.

## Isolation of the moment

At an angular-constant height \(J\) of \(f^m\), the support
inequalities give
\[
             \frac{n(n-3)}{n-1}m\le J\le nm=sp-1.
\]
Since \(s\le n-1\), the lower bound is strictly greater than
\((s-2)p-1\). Indeed \(n(n-3)/(n-1)\ge n(s-2)/s\), and multiplying
the latter by \(m\) gives \((s-2)p-(s-2)/s\).
The only possible interval denominators divisible by \(p\) are
\((s-1)p\) and \(sp\). The former has odd height and gives zero
under the interval integral. Equality at the latter height requires
every factor to come from the top face. Thus, for \(p>s\) and integral
coefficients at a place above \(p\),
\[
  sp\mathcal A(f^m)\equiv [z^m]P(z)^m
       = A^m\binom{em}{m}(-c)^{(e-1)m}\pmod p .
\]
Here congruence at a coefficient place means reduction in its residue
field. The formula removes every term whose denominator is a unit.
For sufficiently large \(p\), both \(A\) and \(c\) are units and
\(em< (es/n)p<p\), so the right side is nonzero. The prime number
theorem in arithmetic progressions gives the stated density.

For general \(k\), take a finitely generated \(\mathbb Q\)-algebra
containing the coefficients of \(f\), \(A\), \(c\), and their two
inverses \(A^{-1},c^{-1}\). A maximal-ideal specialization to a
number field preserves the top face and both support bounds. It
preserves the nonvanishing of \(A,c\) and sends each moment to the
corresponding specialized moment. Only finitely many rational primes
are excluded by integrality for this one specialization.

## Degree and the power test

The largest negative angle is \(-1\) on the top face. Its monomial
\(x^{-1}t^{n+2}\) has the regular lift
\(-2yt^n=x^{-1}(t^{n+2}-t^n)\); both literal terms satisfy the bounds.
The other top angles range from zero to \(e-1\), with nonnegative
heights since \(se<n\) and \(s\ge3\). Thus such phases exist.

Use the regular sphere normal form with no monomial involving both
\(x\) and \(y\). Passing to it by \(2xy=1-t^2\) cannot increase
degree. At angle \(-v<0\), the literal coefficient is divisible by
\((1-t^2)^v\). If it is nonzero, its highest and lowest \(t\)-powers
differ by at least \(2v\). The support inequalities allow a difference
of at most \(2(n-v)/(n-1)\), so \(v\le1\). At nonnegative angles the
normal monomials have degree at most \(n\), and at angle \(-1\) the
normal term corresponding to the top coefficient is a nonzero multiple
of \(yt^n\). It gives degree exactly \(n+1\).

Put \(z=x/t^2\) in the Laurent ring after inverting \(t\). The leading
part in \(t\)-weight of a power is the power of its leading part.
The leading angular polynomial \(z^{-1}A(z-c)^e\) has order \(-1\)
at \(z=0\), which is indivisible by every power exponent at least two.
The argument also applies to a rational root, by taking the leading
coefficient of its Laurent expansion at \(t=\infty\) over
\(\overline{k}(z)\).

## Coherent Frobenius on a supersingular quotient

The branch exponents at \(0,c,\infty\) are \(n-1,e,n+1-e\).
The first is coprime to \(n\), so the Kummer cover is connected
and has degree \(n\). For the character one, the finite eigensheaf
generator is \(W\), and its generator at infinity is
\(W/z^2\), since \(1<\deg(z^{n-1}P)/n<2\). The corresponding
coherent line is generated by the overlap class of \(W/z\).

For character \(s\), the finite generator is
\(E_s=W^s/z^{s-1}\). The floors of the two finite branch contributions
are \(s-1\) and zero, using \(se<n\). The remaining pole at infinity
has order \(1+s(e-1)/n\) in the base coordinate, strictly between
one and two. Thus the generator at infinity is \(E_s/z^2\) and
\([W^s/z^s]\) generates the source coherent line. Frobenius sends
this overlap function to
\[
       W^{sp}/z^{sp}=W z^{-m-1}P(z)^m .
\]
The coefficient at \(W/z\) is exactly the nonzero binomial scalar
in the moment calculation. This is the selected map of lines.

For supersingularity, fix such a prime and let \(q=p^h\). The
congruences give \(n\mid q+1\). Over the algebraic closure of the
residue field, the Fermat curve \(u^n+v^n=1\) dominates the auxiliary
curve. Take \(z=cu^n\) and \(W=C u^{n-1}v^e\) with
\(C^n=A c^{n-1}(-c)^e\). The Fermat curve in turn is dominated by
the Hermitian curve \(U^{q+1}+V^{q+1}=1\) by raising each coordinate
to the power \((q+1)/n\). This is the classical supersingular
Fermat case. [Aoki, Introduction and Corollary 2.5](https://rikkyo.repo.nii.ac.jp/records/16621).

The Hermitian curve over \(\mathbb F_{q^2}\) has \(q^3+1\) points.
In the affine chart the norm map to \(\mathbb F_q\) counts
\(2(q+1)+(q-2)(q+1)^2=q^3-q\) points. There are \(q+1\) points
at infinity. Its genus is \(q(q-1)/2\), so equality in the upper
Weil bound forces every eigenvalue of Frobenius over
\(\mathbb F_{q^2}\) to be \(-q\). Pullback on first cohomology for a
nonconstant map of smooth projective curves is injective by the
pushforward-degree identity. After taking a finite field over which
both maps are defined, the resulting Frobenius slopes on the auxiliary
Jacobian are therefore all \(1/2\).

For example, take \(n=10,s=3,e=2,A=c=1\) and
\(f=-2yt^{10}-2t^{10}+xt^8\). Its second power has angular constant
term \(6t^{20}-2t^{18}\), hence \(\mathcal A(f^2)=24/133\).
At \(p=7\) the normalized moment is \(21\mathcal A(f^2)\equiv6\pmod7\),
agreeing with \(\binom42\). The auxiliary cover is a quotient of the
Hermitian curve with \(q=49\).

## A density bounded away from zero as the conductor grows

**Corollary.** Let \(n=3^a\) with \(a\ge2\), and take
\(P(z)=A(z-c)^2\) with \(A,c\ne0\). Prescribe the same leading
angular polynomial on \(j+2i=n\), and use the strip
\[
                         j+i\ge n-1 .
\]
Every fixed regular phase in this strip has least degree \(n+1\)
and is an absolute nonpower. It has nonzero moments certified by
nonzero coherent Frobenius maps on a set of rational primes of
natural density
\[
                        \frac{1-3^{1-a}}4 .
\]
The whole auxiliary Jacobian is supersingular at those primes.

**Proof.** The regular lifts of the leading face in angles
\(-1,0,1\) satisfy the strip. In angle \(-v\) the two support bounds
permit a spread of at most \(1+v\). Regularity needs a spread of
at least \(2v\) because of the factor \((1-t^2)^v\), so \(v\le1\).
Nonnegative angles have normal degree at most \(n\), and the leading
negative term still forces degree \(n+1\). The order \(-1\) at
\(z=0\) gives the same power test.

Use all residues \(u\) with \(1\le u<n/2\) and \(u\equiv-1\pmod3\).
At large \(p\equiv u^{-1}\pmod n\), set
\[
 m_0=(up-1)/n,\qquad
 \epsilon=\begin{cases}0&u\text{ odd},\\1&u\text{ even},\end{cases}
 \qquad s=u+\epsilon n,\qquad m=m_0+\epsilon p .
\]
Now \(s\) is odd and \(s<2n\). Consequently the strip bound
\(B=n-1\) is at least \(n(s-2)/s\), and the interval calculation
still isolates \(sp\). Because \(2m_0<p\), Frobenius on
\((1+z)^p\) gives
\[
      \binom{2m}{m}\equiv
      \binom{2m_0}{m_0}\binom{2\epsilon}{\epsilon}\not\equiv0\pmod p .
\]
The scalar factors \(A^m(-c)^m\) are units outside finitely many
primes after specialization. For the coherent interpretation, the
character-\(u\) line has basis \([E_u/z]\) with \(E_u=W^u/z^{u-1}\),
since \(2u<n\). When \(\epsilon=1\), the source overlap function
is \(E_u P(z)/z^2\) and its class is \(P'(0)[E_u/z]\). This is a
unit multiple at the same good primes. The target of its Frobenius
is still the character-one line with entry \([z^m]P(z)^m\).

For every such residue \(u\), the unit \(-u\equiv1\pmod3\) has order
dividing \(3^{a-1}\) modulo \(n\). Raising \(u\) to that odd power
therefore gives \(-1\), and the same holds for \(p\). The Hermitian
quotient proves supersingularity. Write \(n=6v+3\). There are
exactly \(v\) residue representatives \(u=2,5,\ldots,3v-1\) in
the specified half interval, and \(\varphi(n)=4v+2\). Their prime
progressions are disjoint, giving the stated density after removing
the finite exceptions. \(\square\)

For \(n=9,u=2,p=23\), one has \(s=11\) and \(m=28\).
With \(A=c=1\), the isolated residue is
\(\binom{56}{28}\equiv2\binom{10}{5}\equiv21\pmod{23}\).
