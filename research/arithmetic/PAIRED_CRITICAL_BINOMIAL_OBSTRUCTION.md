# A uniform obstruction to binomial residual critical polynomials

In the normalized even-numerator problem, a binomial residual critical
polynomial never gives paired simple critical values. It forces a specific
four-term numerator with full symmetric monodromy. This excludes one uniform
degree-twenty-six construction strategy; the general rational existence
question remains open.

## Theorem 1. The binomial residual ansatz has symmetric monodromy

Let b>=3, let N be a monic complex polynomial of degree b with a double
zero at u=1, and suppose

\[
               2uN'(u)-N(u)=(u-1)S(u),\qquad
               S(u)=c(u^{b-1}-k),\quad c\ne0.                  \tag{1}
\]

Then, necessarily,

\[
 c=2b-1,\qquad k=-\frac1{(2b-1)(2b-3)},
\]
\[
 N(u)=u^b-\frac{2b-1}{2b-3}u^{b-1}+\frac{u+1}{2b-3}.           \tag{2}
\]

For f(z)=N(z^2)/z, all 2b finite critical points are simple and nonzero.
Exactly two lie over zero. All the remaining 2b-2 critical points have
different nonzero critical values. Its finite branch passport is therefore

\[
                (2,2,1^{2b-4}),\quad
                (2,1^{2b-2})^{\,2b-2},
\]

and its geometric monodromy is S_(2b). Its discriminant is not a nonzero
constant times a square over C[s].

For the rational square-discriminant target b is odd, but this obstruction
holds for every b>=3 under the normalized double-zero hypothesis (1).

**Proof.** Put d=b-1. The operator 2u d/du-1 multiplies the coefficient
of u^j by 2j-1, so is invertible on complex polynomials. Equation (1)
therefore determines

\[
 \frac{N(u)}c=
 \frac{u^b}{2b-1}-\frac{u^{b-1}}{2b-3}-ku-k.
\]

Leading coefficients give c=2b-1. Substituting u=1, where N vanishes,
gives the displayed k and hence (2). The same equation then gives
N'(1)=0.

The binomial S is squarefree, S(0)!=0 and S(1)!=0. The critical identity

\[
                  z^2 f'(z)=(z^2-1)S(z^2)                    \tag{3}
\]

therefore has 2b distinct nonzero roots. In particular N has exactly
a double zero at one.

For a root u of S, reduce N modulo u^d-k to write

\[
 N(u)=\alpha u+\beta,\qquad
 \alpha=-\frac{2ck(b-1)}{2b-1},\quad
 \beta=-\frac{2ck(b-1)}{2b-3}.
\]

Both constants are nonzero, and their ratio is
rho=beta/alpha=(2b-1)/(2b-3)>1. The squared critical value at this
u is

\[
 M(u)=\frac{N(u)^2}{u}
      =\alpha^2u+2\alpha\beta+\frac{\beta^2}{u}.
\]

For distinct roots u,v, equality M(u)=M(v) is equivalent to

\[
                              uv=\rho^2.
\]

But u^d=v^d=k, so this equality would force k^2=rho^(2d).
Here 0<k^2<1<rho^(2d), which is impossible. All the M(u) are distinct.
They are nonzero too: N(u)=0 would give u=-rho, whereas
|(-rho)^d|>1>|k|.

Each u-root lifts to the two critical points z=+-sqrt(u), whose values
are opposite, nonzero, and have square M(u). Thus all 2b-2 of their
values are distinct. The remaining pair z=+-1 maps to zero.
This proves the exact finite passport.

The map is transitive: N(z^2)-sz is irreducible over C(s), because as
a polynomial linear in s its coefficients -z and N(z^2) are coprime
and N(0)=1/(2b-3)!=0. Its poles have orders (2b-1,1). The long
cycle at infinity forces primitivity. A nonzero finite branch has a
transposition as inertia, so the primitive group is the full symmetric
group. For an elementary justification, form the graph of all transpositions
conjugate to this one inside the group. Its components are blocks;
primitivity makes the graph connected. The edge transpositions of a
connected graph generate the symmetric group.

Every nonzero finite branch value contributes a simple zero to the
discriminant. Hence the discriminant cannot be a constant times a square. ∎

## Corollary 2. Consequence for the degree-twenty-six search

In the normalized coefficient criterion
N(u)=(u-1)^2G(u), deg(G)=11, the choice S(u)=c(u^12-k)
forces

\[
 H(z)=z^{26}-\frac{25}{23}z^{24}+\frac1{23}z^2,
 \qquad a=\frac1{23}.
\]

It has one double-transposition branch and twenty-four separate
transposition branches, and full S26. Thus the binomial S strategy
cannot produce the required thirteen double-transposition branches.

The scope is the residual binomial S after the zero critical pair has
been normalized to u=1. General S, general four-term numerators,
and general rational degree-twenty-six pencils are outside this obstruction.

## Exact evidence and dependencies

The proof is elementary except for the usual interpretation of simple
critical points as inertia transpositions. Its normalized coefficient
setting and even-degree discriminant identity are the ones in
[Proposition 4 of the uniformity note](ALTERNATING_HURWITZ_UNIFORMITY_AND_REAL_NODES.md).

The [exact checker](verify_paired_critical_binomial.py) verifies in
degrees ten and twenty-six:

* the complete derivative and normalization identities;
* the reduced linear numerator and the collision identity;
* squarefreeness and nonvanishing of the squared-critical-value polynomial;
* the full polynomial discriminant identity and all nonzero simple factors.

The [receipt](paired_critical_binomial_verification.json) pins the checker,
Python 3.12.3 and SymPy 1.14.0. The written proof supplies the quantifier
for every b>=3. The two exact controls do not prove that quantifier.

Run from the repository root:

~~~sh
.venv/bin/python research/arithmetic/verify_paired_critical_binomial.py
~~~

This is an internal exact obstruction for one ansatz. No rational
degree-twenty-six impossibility theorem, formal verification, external
review or independent external replay is asserted.

