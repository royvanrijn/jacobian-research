# Complete affine singularity strata of the polynomial gap-21 family

14 November 2026. PF2R7COLL1 classifies all affine singularities in the
polynomial coefficient component L of PF2R7FAM1. It gives the exact
34 exceptional parameters and all their singularity types. The Laurent
component K remains a separate collision calculation.

PF2R7COLL2 applies the triple-point strata to the residual degree-six
Keller ledger and corrects the interpretation of the old sharp braid
sample. No global cover, additional image, whole-r7 exclusion or JC2
theorem is claimed. This is an authored written local-curve proof with
exact algebraic certificates, without formalization or independent review.

## 1. The answer

Use the normalized family L from [PF2R7FAM1](../jc2-r7-image-family/PROOF.md),
with b=1,a=2d,c=d^2,B=-7d/6 and the displayed polynomial D,E. In particular

    u=t^3+d*t, p=u^2+u.

Put z=d^3 and define the following polynomials (the checker uses d):

    C=3z+2,
    U=2187z^3+3402z^2+1512z+196,
    M=243z+64,
    T=3761479876608z^6+15898976317440z^5+25767272444544z^4
      +20621286822727z^3+8479009232730z^2+1615258153728z
      +91424858112.

The four polynomials in d are squarefree, pairwise coprime and nonzero
at d=0. Their d-degrees are3,9,3,18. The complete affine classification is:

| Parameter stratum | Number of d values | Affine singularities |
|---|---:|---|
| d*C*U*M*T nonzero | all other values | 12 ordinary nodes |
| d=0 | 1 | 3 ordinary triple points + 3 nodes |
| U=0 | 9 | 1 ordinary triple point + 9 nodes |
| M=0 | 3 | 1 ordinary triple point + 9 nodes |
| T=0 | 18 | 1 ordinary tacnode A3 + 10 nodes |
| C=0 | 3 | A2 cusp + A5 two-branch tangency + 8 nodes |

Thus exactly34 complex parameter values fail to have twelve distinct
ordinary nodes. An ordinary triple point has three smooth branches with
pairwise different tangent lines (type D4). The final row is the retained
nonimmersed class I; Section6 explains the normalization comparison.

## 2. Unordered collision coordinates preserve the exceptional cases

For two normalization parameters s,t put v=s+t,w=s*t. Define the
symmetric divided differences

    A(v,w)=(p(s)-p(t))/(s-t),
    B(v,w)=(q(s)-q(t))/(s-t).

The identity is polynomial, including on the diagonal. Direct expansion gives

    A=A1*A2,
    A1=v^2-w+d,
    A2=v^3-3v*w+d*v+1.

A1 means that the two u values agree; A2 means that they sum to -1.
On A1=0 we have w=v^2+d and

    B=-H9(v,d)/162,
    H9=162v^9+621dv^7-270v^6+756d^2v^5-702dv^4
       +(240d^3+126)v^3-432d^2v^2+(-64d^4+168d)v
       +32d^3-14.

On A2=0 one necessarily has v nonzero, and

    w=(v^3+dv+1)/(3v), B=-v^6*H3(v,d)/162,
    H3=2v^3+9dv+6.

The leading coefficients of H9,H3 are the constant units162,2.
In particular v=0 on the first branch is retained; it can be a genuine
collision at d^3=7/16. An unfiltered resultant in w has an extraneous
factor v^10 from leading-coefficient degeneration. It is not interpreted
as ten finite collision pairs.

All identities and resultants below are reconstructed from the original
family p,q, rather than from numerical samples. Up to nonzero rational factors,

    Disc_v(H9)=T, Disc_v(H3)=C,
    Res_v(H9,2v^3+2dv-1)=C,
    Res_v(H9,3v^2+4d)=1,
    Res_v(H3,v^3+4dv+4)=C.

The third equation detects overlap of A1 and A2; the last two detect
diagonal pairs s=t in the respective charts. Therefore for C*T nonzero
there are exactly12 distinct off-diagonal unordered collision pairs,
with the charts disjoint and their individual solutions simple.

For an immersed curve, a simple off-diagonal collision pair is
transverse. Indeed the map (s,t) to (v,w) is etale off the diagonal.
The determinant of the derivatives of A,B differs from the determinant
of the two image tangent vectors by a unit there. Along either smooth
A_i chart, the multiplicity of the B zero is the intersection
multiplicity of the two image branches. These facts also identify a
unique double zero with a smooth two-branch tangency of order two.

## 3. Exhaust every possible triple or higher image value

Write the exact remainder of q modulo t^3+d*t-u as

    q=r2(u,d)*t^2+r1(u,d)*t+r0(u,d).

These three explicit polynomials are reconstructed in
[collisions.py](collisions.py). Three distinct parameters with the same
p value have at most two distinct u values, since p=u^2+u. Hence any
triple or higher value contains either three parameters with the same
u, or two with one u and another with -1-u. This is an exhaustive
division, not a generic assumption.

In the same-u case, q is constant on the three roots of the monic
cubic t^3+d*t-u precisely when r2=r1=0. The exact resultant is

    Res_u(r2,r1)=(7/5668704)*d^3*U.

For d nonzero the degrees of r2,r1 are fixed two and three, so this is
an exact existence criterion. At U=0 a frozen modular gcd certificate gives

    gcd(r2,r1) proportional to 126u+405d^6+504d^3+196.

Thus exactly one u value occurs. Its cubic has three distinct roots:
otherwise its repeated parameter would have u'=0 and, since r2=r1=0,
also q'=0. This contradicts the complete immersion criterion C nonzero
from PF2R7FAM2. Consequently every U root gives exactly one genuine
three-branch value. T is nonzero there, so its branches are pairwise
transverse. No higher same-u value is possible in a cubic fibre.

For mixed triples choose the unique same-u pair, with A1=0 and H9=0.
Its u value is u0=-v*(v^2+d), and its common q value is

    q0=r0(u0,d)-r2(u0,d)*(v^2+d).

Define J(v,d) as the exact resultant in a third parameter r of

    r^3+d*r-(-1-u0),
    r2(-1-u0,d)*r^2+r1(-1-u0,d)*r+r0(-1-u0,d)-q0.

No parameter denominator is introduced; the cubic is monic. The next
resultant is, up to a nonzero rational constant,

    Res_v(H9,J)=C^2*M.

Outside the nonimmersed C locus, mixed triples therefore occur exactly
at M=0. A frozen modular gcd gives the unique same-u pair sum

    gcd(H9,J) proportional to 16v-45d^2.

The same-u pair is off the diagonal by Section2. Its two u values
cannot coincide with -1-u0, since that is precisely the A1/A2 overlap
locus C. The resultant supplies at least one distinct third parameter.
There cannot be a fourth: a partition2+2 would give two distinct
same-u pairs and hence two roots of this linear gcd; a partition3+1
would force the disjoint U locus. Thus there is exactly one ordinary
triple point at each M root, with all other image values at most double.

This proves that every triple or higher image value is exhausted by
d=0, U=0, M=0, or the already separated nonimmersed C locus.

## 4. The special d=0 sample and the tacnode stratum

At d=0 the remainder is

    q=t*A0(u), A0=u^3+5u^2/3+7u/9+7/81, u=t^3.

A0 has three distinct nonzero roots. Each supplies three distinct
parameters with q=0 and one common p value. The mixed resultant is
nonzero at d=0, so these three triple values are distinct. The tangent
discriminants are also nonzero there; all are ordinary triple points.
The remaining three unordered collision pairs give three nodes.

At T=0, the frozen modular gcd of H9 and its derivative is monic linear.
Thus H9 has precisely one double root and all other roots simple,
at every root of T. C is nonzero, so the parameters are immersed,
all pairs are off-diagonal, and the two collision charts do not overlap.
The triple loci d,U,M are disjoint from T. Hence this is one two-branch
tangency of intersection multiplicity two, an ordinary A3 tacnode;
the ten remaining pairs give ten ordinary nodes.

At U or M roots all12 pairs are simple and transverse. Three pairs
combine into the unique ordinary triple point, leaving nine nodes.
Outside all exceptional strata the12 pairs give twelve distinct nodes.
No additional finite singularity is possible: the parametrization is
finite and birational, and any singular point either has more than one
normalization parameter or is nonimmersive at one parameter. The latter
parameters are exactly C=0 by PF2R7FAM2.

## 5. Modular certificates use no algebraic-root embedding

For each of U,M,T the certificate bundle stores a monic linear g and
rational polynomial multipliers satisfying

    L0*f+L1*h=g modulo the exceptional polynomial in d.

The checker reconstructs f,h and that exceptional polynomial. It
replays the identity by sparse Fraction multiplication and polynomial
reduction in d, and separately verifies that g divides both generators
by evaluation at its polynomial root in the quotient ring. The leading
coefficient of g is the constant unit1. These two directions prove the
gcd is exactly linear after every exceptional specialization, even if
one does not assume the exceptional polynomial irreducible. Mutating
a multiplier fails the replay.

A direct algebraic-root embedding attempt was capped at45 seconds and
gave no completed result. Early Singular exports with fractional
exponents/coefficients produced parser errors despite exit code zero;
they were rejected. Clearing numeric denominators before export repaired
the input. Those failed attempts are retained in
[failed_attempts.json](failed_attempts.json), and no printed basis is a
proof premise. The final three modular certificates occupy about16KB.

## 6. The nonimmersed row and the degree-six global consequence

At C=0 the line chart is the retained PF2R7CUSP1 class I, with complete
affine type A2+A5+8A1. For its old coefficients b0=-2/3,d0=4/9,
scaling by lambda^3=b0 gives d^3=d0^3/b0^5=-2/3. The weighted relations
a0*b0=2d0 and6B0*b0=-7d0 identify the line chart; its graph functions
then identify c,D,E. The three cube-root choices give all three C roots.
This row uses that earlier exact nonimmersion classification, rather
than inferring its local type from a vanishing derivative resultant.

For a hypothetical degree-six Keller map containing an immersed e2,f1
prime with an L image, PF2OMIT6 gives ambient length two at every point
of that prime. At an ordinary triple value there are three distinct
normalization points of C, giving three disjoint omitted factors of
length two. They consume the full degree-six finite fibre. Therefore:

1. The original source A2 has no point above any such target value.
2. No other omitted image, ramified or unramified, can pass through it.

This is a necessary matching restriction on the unknown B and optional
H. It is not a contradiction: no surjectivity theorem for a hypothetical
Keller map is imported. It supplies three forbidden intersection values
at d=0 and one at each U or M parameter, with explicit equations for
their normalization parameters in Section3.

The old `jc2-normal-r7/braids-generic.json` sample is exactly d=0.
Its stored case label `generic` is historical; it is not generic in
the complete L family. Its saved six braid words and the one-sample
nontransitivity calculation retain their scope. They cannot be promoted
to a result on the twelve-node stratum. This application and corrected
sample interpretation constitute PF2R7COLL2.

The next global inputs remain the Laurent-family collision strata,
actual B,H equations and their intersections, the complement relations
for the actual curve configuration, and a finite normal cover with the
required source A2. No whole normal-r7 or JC2 closure is claimed.
