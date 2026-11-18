# Reciprocal critical pairing: a classification and a degree-26 obstruction

The reciprocal pairing that produces the rational degree-ten pencil has a
uniform one-parameter description. Its extension with one reciprocity sign
has no rational point in degree twenty-six: a degree-twelve obstruction
polynomial has no zero in the projective line over F11.

This excludes the uniform-sign reciprocal ansatz. The mixed-sign reciprocal
cases remain possible, with an explicit factor-degree classification below.
The general rational degree-twenty-six existence problem remains open.

## Setting

Let b=2h+1, n=2b=4h+2, and let N in Q[u] be monic of degree b, with
N(1)=N'(1)=0. Put

\[
 2uN'(u)-N(u)=(u-1)S(u),\qquad
                      f(z)=\frac{N(z^2)}z.                    \tag{1}
\]

Assume the normalized paired-simple-critical conditions: S is squarefree,
S(0)S(1)!=0, N is coprime to S, and the squared critical values

\[
                              M(u)=N(u)^2/u
\]

at the roots of S occur precisely twice each. These are the conditions of
[the normalized coefficient criterion](ALTERNATING_HURWITZ_UNIFORMITY_AND_REAL_NODES.md),
with the additional reciprocal hypothesis specified next.

Suppose the roots of S are paired by the fixed-point-free involution
u->lambda/u, with lambda in Q*. Choose rho with rho^2=lambda and set
K=Q(rho). Pairing gives M(u)=M(lambda/u) at every S-root.

The **uniform-sign case** means that, after choosing the sign of rho,

\[
                   N(u)+\frac u\rho N(\lambda/u)=0
                           \quad\text{modulo }S.               \tag{2}
\]

Laurent expressions modulo S are legitimate because S(0)!=0.
An opposite sign in (2) is absorbed by replacing rho with -rho.

## Theorem 1. Uniform-sign reciprocity has an explicit ODE classification

Under these hypotheses rho is rational. Put r=1/rho, v=u/rho, and

\[
 F(v)=\rho^{-b}N(\rho v),\qquad
 \mathsf S(v)=\frac{S(\rho v)}{(4h+1)\rho^{2h}}.
\]

Then F is monic of degree 2h+1, while mathsf S is a monic reciprocal
polynomial of degree 2h with constant coefficient one. They satisfy

\[
 2vF'(v)-F(v)=(4h+1)(v-r)\mathsf S(v),                         \tag{3}
\]
\[
 v^{2h}F(v)+v^{2h+1}F(1/v)=(v+1)\mathsf S(v)^2.                \tag{4}
\]

Equivalently mathsf S is the unique polynomial with coefficients
s_0=s_(2h)=1,

\[
 s_1=-\frac{(4h+1)r+4h-1}{4h-3},\qquad
 s_j=-\frac{4j-4h-3}{4j-4h-1}s_{j-1}
                   \quad(2\le j\le2h-1).                     \tag{5}
\]

Given mathsf S, equation (3) determines F coefficient by coefficient:
divide the coefficient of v^j in its right side by 2j-1.
The zero branch requires F(r)=0.

Conversely, for every nonzero r, construction (5) and (3) gives monic
polynomials satisfying the full identity (4). A value F(r)=0 supplies
the normalized double zero, provided mathsf S(r)!=0. To obtain the
desired passport one must still impose squarefreeness, absence of
fixed roots +-1, nonvanishing of the other critical values, and
distinctness of the values attached to different reciprocal pairs.

**Proof.** If rho were quadratic irrational, apply its conjugation to
(2), using N,S,lambda in Q. This gives the same congruence with the
opposite sign. Adding forces N=0 modulo S, contradicting coprimality.
Thus rho is rational.

The fixed-point-free reciprocal pairing makes the constant-to-leading
coefficient ratio of S equal to lambda^h. After the stated normalization,
mathsf S is therefore reciprocal with constant coefficient one.
Equation (3) follows by changing variables in (1).

Set

\[
                    A(v)=v^{2h}F(v)+v^{2h+1}F(1/v).
\]

It is a monic reciprocal polynomial of degree 4h+1. Equation (2) says
mathsf S divides A. Using (3) and reciprocity gives

\[
 2vA'-(4h+1)A
  =(4h+1)(v^{2h+1}-rv^{2h}+rv-1)\mathsf S.                     \tag{6}
\]

At each simple nonzero root of mathsf S, the right side and A vanish,
so A' vanishes as well. Hence mathsf S squared divides A.
Its quotient is monic of degree one and reciprocal, hence v+1.
This proves (4).

Substitution into (6) gives the first-order equation

\[
 4v(v+1)\mathsf S'
 -((4h-1)v+4h+1)\mathsf S
 =(4h+1)(v^{2h+1}-rv^{2h}+rv-1).                               \tag{7}
\]

Its coefficient equations are exactly (5), with the two boundary
coefficients one. None of the displayed denominators vanishes for
integer h>=1. The middle recurrence is symmetric under reversal,
so s_(2h-j)=s_j.

For the converse, those coefficient equations prove (7). Define F
from (3), which is possible uniquely since 2j-1 never vanishes.
Equation (6) follows. Subtracting the expression (v+1)mathsf S^2
shows that their difference D satisfies 2vD'-(4h+1)D=0.
This operator has no nonzero polynomial kernel, since 4h+1 is odd.
Thus (4) holds identically.

Finally N(1)=0 is F(r)=0, and (3) then gives F'(r)=0.
When mathsf S(r)!=0 the zero has exactly order two. The remaining
conditions in the statement exclude every collision or branch
multiplicity other than the prescribed reciprocal pairs. ∎

This classification concerns a rational one-parameter auxiliary family.
Solving F(r)=0 and the nondegeneracy conditions is a separate arithmetic
step; the family does not automatically have alternating monodromy.

## Theorem 2. No rational uniform-sign reciprocal model in degree 26

For h=2, equation F(r)=0 becomes

\[
 F(r)=\frac{2r(2r-15)(4r^3-6r^2-9r-15)}{25}.
\]

Taking r=15/2 and rho=2/15 recovers exactly the known monic
[rational degree-ten numerator](RATIONAL_DEGREE_TEN_ALTERNATING_PENCIL.md):

\[
 N(u)=u^5-\frac{96}{25}u^4+\frac{3129}{625}u^3
                  -\frac{14006}{5625}u^2+\frac{1804}{5625}u
                  +\frac{16}{5625}.
\]

For h=6, the corresponding identity is

\[
                       F(r)=\frac{2rP_6(r)}{193849929},
\]
\[
 \begin{split}
 P_6(r)={}&3516552r^{12}-5171400r^{11}+8166600r^{10}
 -14469000r^9+32478600r^8\\
 &-207017160r^7-117119800r^6-26033800r^5
 +79781000r^4-234304200r^3\\
 &+1527631560r^2+1769265225r+4846248225.
 \end{split}
\]

It has no rational root. Consequently a rational degree-twenty-six pencil
satisfying the setting above cannot have uniform-sign reciprocal pairing.

**Proof.** The identities follow from the exact coefficient recurrence
(5) and integration rule (3), rather than interpolation.
Modulo 11 the descending coefficient list of P_6 is

\[
                     (6,8,2,4,0,2,5,10,2,2,4,2,5).
\]

Its values at 0,1,...,10 are

\[
                        (5,8,1,6,8,3,3,4,3,5,7).
\]

Every value is nonzero, as is the leading coefficient 6. Thus its
degree-twelve homogenization has no zero in P1(F11), including infinity.
A rational root would reduce to such a projective zero after choosing
coprime numerator and denominator. This is impossible.
Theorem 1 makes r nonzero and rational, so F(r)=0 cannot hold. ∎

The obstruction is over Q. This argument does not exclude algebraic
parameters r or real-algebraic models.

## Corollary 3. Exact remaining reciprocal factor patterns

In degree twenty-six, any rational model satisfying the reciprocal
pairing setting must have S reducible over K=Q(sqrt(lambda)).
More precisely:

* If lambda is a rational square, S splits over Q into nonconstant
  opposite-sign reciprocal factors of degrees (2,10), (4,8), or (6,6).
* If lambda is nonsquare, S splits over K into two conjugate
  opposite-sign reciprocal factors, both of degree six.
  If S is irreducible over Q, K is a quadratic subfield of the
  degree-twelve field defined by any S-root.

These are necessary patterns, without an existence assertion for them.

**Proof.** In the etale algebra K[u]/S, paired equality factors as

\[
 \left(N(u)-\frac u\rho N(\lambda/u)\right)
 \left(N(u)+\frac u\rho N(\lambda/u)\right)=0.
\]

The factors have no common root of S, since N is nonzero there.
They therefore partition its roots into the two signs and give coprime
factors S_-,S_+ over K. A reciprocal pair retains its sign, so each
nonempty factor has even degree. Theorem 2 excludes either factor
being constant.

When rho is rational these factors are rational, giving the three
possible even degree splits of twelve. When rho is quadratic
irrational its conjugation swaps the signs, so the factors are
conjugate and have equal degree six. If S is irreducible over Q,
failure of irreducibility over the quadratic field K forces K to lie
in the S-root field, by the degree formula for adjoining K. ∎

In particular a reciprocal S irreducible over Q(sqrt(lambda)) cannot
produce the rational degree-twenty-six target. The search must retain
genuine mixed signs. The corollary does not cover nonreciprocal
pairings of critical points.

## Exact evidence and scope

The [checker](verify_paired_critical_reciprocal.py) and
[receipt](paired_critical_reciprocal_verification.json) verify the complete
symbolic recurrence, ODE and squared reciprocity identities for h=2,6,
the exact degree-ten recovery, and every point of P1(F11) in the
degree-twenty-six obstruction. It also rejects corrupted finite-root
and infinity-root controls.

The uniform classification and quadratic-field factor arguments are
written proofs. The finite controls do not supply their all-degree
quantifiers. All arithmetic is exact, with the pinned Python 3.12.3
and SymPy 1.14.0 environment.

~~~sh
.venv/bin/python research/arithmetic/verify_paired_critical_reciprocal.py
~~~

This is a complete obstruction and structural reduction for the stated
reciprocal ansatz. Rational degree-twenty-six existence, mixed-sign
solutions, and uniform rational constructions remain open.
No formal verification, external human review or independent external
replay is claimed.

