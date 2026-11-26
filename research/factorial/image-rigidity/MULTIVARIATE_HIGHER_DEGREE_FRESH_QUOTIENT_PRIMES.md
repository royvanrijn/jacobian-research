# Fresh quotient primes and binomial factors in arbitrary degree

The [quadratic cyclotomic product test](MULTIVARIATE_QUADRATIC_CYCLOTOMIC_PRODUCT_PRIMES.md)
forces a fresh prime of prescribed quotient order for a
noncyclotomic quadratic direction. A related test works in
arbitrary common-factor degree. One auxiliary prime exceeding
an effective degree threshold either leaves a large common
prime divisor in the two pair products or shows that every
surviving common irreducible factor is a factor of a rational
binomial. A separate power transform detects the latter
case. The necessary conditions apply to the unchanged
ordinary factorial arrays, with no bound on the row list.

## A degree-uniform primitive prime for an algebraic trace

Let \(\gamma\ne0\) be a nontorsion algebraic number, meaning
that it is not a root of unity, and set
\(\tau=\gamma+\gamma^{-1}\). Let \(P_\tau\) be the primitive
integral minimal polynomial of \(\tau\), with positive leading
coefficient. Fix \(E\ge1\) with \(e=\deg P_\tau\le E\).
For an odd prime \(q\), define \(C_q\in\mathbb Z[Y]\) by
\[
      C_q(X+X^{-1})=X^{-(q-1)/2}\Phi_q(X).
\]
There is an effectively computable \(Q(E)\) such that
for every prime \(q>Q(E)\) the integer
\[
                     R_q(\tau)=\operatorname{Res}(P_\tau,C_q)
\]
has a rational prime divisor \(p\ne q\). Every such prime
satisfies
\[
                  q\mid p^f-1\ \text{or}\ q\mid p^f+1
                       \quad\text{for some }1\le f\le e.      \tag{1}
\]
In particular \(p\ge(q-1)^{1/e}\). The threshold
does not depend on the heights of the trace and quotient.

Here is a direct proof for prime \(q\), including quotients
whose fractional ideals have nonprincipal numerator or
denominator. Write \(P_\gamma\) for the primitive minimal
polynomial of \(\gamma\), \(l=\deg P_\gamma\le2E\), and
let \(h(\gamma)\) be the absolute logarithmic Weil height.
For any conjugate \(\gamma_i\) take \(x_i=\gamma_i\) or
\(\gamma_i^{-1}\) so that \(|x_i|\le1\). A complex
linear form in two logarithms gives
\[
 \log|1-x_i^q|\ge -K_E\max(1,h(\gamma))\log(3q).
\]
In fact, when \(x_i^q\) is within \(1/2\) of \(1\), its
nonzero principal logarithm is
\(q\log x_i-2k\log(-1)\) with \(|k|\le(q+1)/2\).
The principal logarithms have magnitudes at most
\(2E h(\gamma)+\pi\) and \(\pi\), and their bases lie
in a number field of degree at most \(2E\).
Apply [Matveev's complex Corollary 2.3](https://doi.org/10.1070/IM2000v064n06ABEH000314)
and use \(|\log(x_i^q)|\le2|x_i^q-1|\) in this disk.
The bound outside the disk is immediate on increasing \(K_E\).

The geometric sum for \(\Phi_q\) now gives
\[
 |\Phi_q(\gamma_i)|\ge
       \tfrac12\max(1,|\gamma_i|)^{q-1}|1-x_i^q|.
\]
Multiplication with the leading coefficient contribution
to the resultant yields
\[
  \log|\operatorname{Res}(P_\gamma,\Phi_q)|
       \ge l\Bigl((q-1) h(\gamma)-\log 2
              -K_E\max(1,h(\gamma))\log(3q)\Bigr).             \tag{2}
\]
There is an effective positive lower bound on \(h(\gamma)\)
for nontorsion numbers of degree at most \(2E\).
It suffices to enumerate the primitive integer
polynomials of degree at most \(2E\) with Mahler measure
at most \(\exp(2E)\), using the coefficient bounds from
elementary symmetric functions. This finite list covers
all candidates with \(h(\gamma)\le1\). Kronecker's
criterion gives height zero only for the roots of unity
and zero. A positive lower bound for the remaining
members can be certified by algebraic isolation of
their root moduli. Hence we can take \(Q(E)\) so that
(2) is greater than \(2E\log q\) for all \(q>Q(E)\).

The polynomial \(X^eP_\tau(X+X^{-1})\) is integral and
primitive. For primitivity, if \(k\) is the largest
index with coefficient of \(P_\tau\) nonzero modulo
\(p\), the coefficient at \(X^{e+k}\) in the transformed
polynomial is that unit modulo \(p\). The primitive
minimal polynomial \(P_\gamma\) divides the transform.
Evaluation at primitive \(q\)-th roots, pairing each
root with its inverse, gives
\[
 \left|\operatorname{Res}\bigl(X^eP_\tau(X+X^{-1}),\Phi_q\bigr)\right|
                                    =|R_q(\tau)|^2 .
\]
No resultant on these factors vanishes, since a
cyclotomic value of the reciprocal trace would force
\(\gamma\) to be a root of unity. Gauss's lemma and
(2) now imply \(|R_q(\tau)|>q^E\).

To bound the contribution of the index prime, put
\(P_\tau(2+Z)=\sum_i b_i Z^i\) and let \(k\) be the first
index whose coefficient is a \(q\)-adic unit.
For a primitive \(q\)-th root \(\zeta\) the number
\(Z=\zeta+\zeta^{-1}-2=(\zeta-1)^2/\zeta\) has
\(q\)-adic valuation \(2/(q-1)\), normalized by
\(v_q(q)=1\). This follows from the Eisenstein polynomial
\(\Phi_q(1+X)\). Take \(Q(E)>2E+1\). All terms with
coefficients divisible by \(q\) have valuation at
least one, and the others have distinct valuations
\(2i/(q-1)<1\). Each evaluation therefore has valuation
\(2k/(q-1)\). Summing over the \((q-1)/2\) cyclotomic
traces gives \(v_q(R_q(\tau))=k\le e\).
The size bound forces a prime \(p\ne q\).
At such a prime the two reduced polynomials have a
common root \(t\) of degree \(f\le e\) over \(\mathbb F_p\),
since \(C_q\) is monic even when the degree of \(P_\tau\)
drops. Every root of \(Y^2-tY+1\) in the residue algebraic
closure has order \(q\) by the cyclotomic identity.
The \(p^f\)-power map fixes or inverts it. This proves (1).

## Transfer to the two factorial parents

For an integral polynomial \(F\) of degree \(m\ge2\) with
leading coefficient \(a\), nonzero constant coefficient,
and roots \(\alpha_i\) counted with multiplicity, define
\[
 {\cal T}_F(Y)=a^{2(m-1)}
       \prod_{i<j}(\alpha_i\alpha_jY-\alpha_i^2-\alpha_j^2).
\]
The symmetric product has per-variable degree at most
\(2(m-1)\), so \({\cal T}_F\in\mathbb Z[Y]\) by the
symmetric degree argument in the
[cyclotomic product construction](MULTIVARIATE_QUADRATIC_CYCLOTOMIC_PRODUCT_PRIMES.md#integral-pair-products).
For odd prime \(q\),
\[
                      \operatorname{Res}({\cal T}_F,C_q)
                                          = I_q(F).           \tag{3}
\]
This follows term by term from
\(\Psi_q(x,y)=(xy)^{(q-1)/2}C_q(x/y+y/x)\).
No root separation or squarefreeness is needed for (3).

Suppose the ordinary adjacent parents \(F,G\) have a
common irreducible factor of degree \(d\ge2\), and
some pair of its distinct roots has nontorsion
quotient \(\gamma\). Every conjugate of the pair
trace \(\tau\) is one of the at most \(\binom d2\)
unordered traces from that factor. In particular
\(e\le\binom d2\), and \(P_\tau\) divides both
\({\cal T}_F\) and \({\cal T}_G\) over \(\mathbb Z\).
It follows from (3) that
\[
                 R_q(\tau)\mid \gcd(I_q(F),I_q(G)).
\]
Consequently every prime \(q>Q(\binom d2)\) forces a
prime \(p\) of the gcd with (1) for some
\(1\le f\le\binom d2\). If \(M\) bounds the original
factorial arguments and the two orders, choose \(q\)
also greater than \(M^{\binom d2}+1\). Then \(p>M\),
so every original parent coefficient is a unit
at \(p\), including after primitive normalization.
All roots of the parents are units at such a prime.
The reduction of \(P_\tau\) then comes from actual
reductions of the pair traces of the common factor:
in a splitting field each conjugate of \(\tau\) is
integral, and the leading coefficient of its
primitive minimal polynomial is a unit at \(p\).
Thus the new prime sees two distinct common residue
roots with exact quotient order \(q\), outside the
entire coefficient support. The forced order can be
tested without finding the candidate factor.

## Reduction to rational binomial factors

Fix any a priori degree bound \(d_0\ge2\) for an
irreducible common factor, and let \(E_0=\binom{d_0}2\).
If for one prime \(q>Q(E_0)\) the common divisor of
\(I_q(F)\) and \(I_q(G)\) has no rational prime
\(p\ne q\) satisfying \(q\mid p^f\pm1\) for
\(1\le f\le E_0\), then every surviving factor of
degree at most \(d_0\) has only torsion quotients.
Take
\[
         L=\operatorname{lcm}
            \{h\ge1:\ \varphi(h)\le d_0(d_0-1)\}.
\]
The set is finite since \(h\le2\varphi(h)^2\).
For two different roots \(\alpha,\beta\) of a
rational irreducible factor of degree \(d\), the
field \(\mathbb Q(\alpha,\beta)\) has degree at most
\(d(d-1)\), by adjoining \(\beta\) after removing
\(\alpha\)'s linear factor over \(\mathbb Q(\alpha)\).
Thus the torsion orders of their quotients divide
\(L\). All \(L\)-th powers of the roots of the factor
are one value \(c\). Galois invariance makes
\(c\in\mathbb Q^\times\). The factor therefore
divides \(z^L-c\). In particular \(c\) is a shared
rational root of the parent power transforms
\(F^{[L]},G^{[L]}\), of multiplicity at least \(d\)
in each. If the two transforms have no common
rational root of multiplicity at least two in both,
the nontorsion-prime test then excludes every common
factor of degree at most \(d_0\). Degree one is
excluded by positivity of the even moment. One may
use \(d_0\) equal to the smaller parent degree.

There are also two smaller ways to use this branch.
First, for a factor of degree \(d\) dividing \(z^L-c\),
all its roots have the same valuation at each \(p\),
namely \(v_p(c)/L\). Thus a prime where no one shared
slope has remaining width \(d\) already excludes the
binomial branch. Second, its possible conjugate-pair
arguments are a finite set depending only on \(d_0\).
A fixed nonreal diffuse contact with conjugate quotient
outside that set cannot be approached by this branch,
regardless of the size of the rational scale \(c\).

This is a finite reduction for any specified pair.
The fresh-prime test does not require that the field
be quadratic, the ratio have modulus one at the
real place, or its prime support be fixed. It does
not prove that the required prime is always missing
for unrestricted ordinary factorial profiles.

## Exact examples in degree four

Let \(H(z)=z^4+z+1\). Its reductions modulo \(2\) and
\(3\) have irreducible factor degrees \(4\) and \(1,3\),
respectively. They are squarefree, so the Galois group
contains a four-cycle and a three-cycle. It is the
symmetric group on four letters, since its order is
divisible by twelve and a subgroup of index two cannot
contain a four-cycle. The polynomial \({\cal T}_H\) is
\[
                Y^6+4Y^5+4Y^4-Y^3-3Y^2-3Y-1 .
\]
It is squarefree modulo \(2\). Because the Galois group
is transitive on unordered root pairs, the six distinct
traces are conjugate. Thus \({\cal T}_H\) is an example
of \(P_\tau\) of degree six. At \(q=11\) its resultant
with \(C_{11}\) has magnitude \(89\), and at \(q=31\)
the magnitude is \(5^6\cdot1303\).
At \(5\) the order condition uses \(5^3\equiv1\bmod31\);
at \(1303\) it uses the degree-one condition.
For \(F=H(1+2z)\) and \(G=H(2+3z+z^2)\), each
integer resultant divides both exact self products,
as in (3).

For the other branch, \(H_0=z^4+2\) is irreducible and
all four of its roots have fourth power \(-2\).
The power transforms at exponent four of
\(H_0(1+z)\) and \(H_0(1+z+z^2)\) both contain
\(z+2\) with multiplicity exactly four. These examples
test the arithmetic of the two branches with planted
factors, without asserting that either pair consists of
factorial moments.
