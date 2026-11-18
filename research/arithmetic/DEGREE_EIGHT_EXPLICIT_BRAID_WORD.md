# A certified braid word for the original degree-eight cover

The literal quadratic-field polynomial

\[
 H(z)=z^8+z^6+ez^2,\qquad
 e=\frac{-1009+94i\sqrt{329}}{10584},\qquad a=-\frac{625}{21952},
\]

has explicitly certified branch cycles that reach the original knot tuple by a
ten-letter braid word. The selected complex embedding is
sqrt(-329)=+i sqrt(329). This strengthens the
[complete Nielsen-orbit result](DEGREE_EIGHT_HURWITZ_ORBIT.md): the calculation
below starts from actual root continuation of this polynomial.

## Theorem 1. Literal model-to-knot braid certificate

Label the eight roots of H(z)-z+a at the base value s=1 by increasing real part,
starting with label zero. The retained certificate gives disjoint rational root
disks proving that these labels are unambiguous.

Let the four branch squares have halfside h=1/500 and respective centers

\[
 q_0=\frac{-1323+752i}{10000},\quad
 q_1=\frac{-811+405i}{10000},\quad
 q_2=\frac{811-405i}{10000},\quad
 q_3=\frac{1323-752i}{10000}.
\]

For each j, take the straight spoke from 1 to q_j+h, traverse the boundary of
its square counterclockwise, and return along the same spoke. Explicitly the
vertices are

\[
 1,\ q_j+h,\ q_j+h+ih,\ q_j-h+ih,\ q_j-h-ih,\
 q_j+h-ih,\ q_j+h,\ 1.
\]

Each square contains exactly one finite branch value. The squares are disjoint,
and the spokes meet only at their common basepoint and avoid the other squares.
The order displayed is counterclockwise around the basepoint. These are four
positive geometric meridians.

Their root permutations, computed by exact continuation certification, are

\[
 g_0=(2\,6)(4\,7),\quad
 g_1=(4\,6)(5\,7),\quad
 g_2=(0\,4)(1\,5),\quad
 g_3=(1\,4)(3\,5).
\]

Use permutation composition pq=p composed with q. Relabel the sheets by the
permutation r given in one-line notation by

\[
 r=[0,4,5,7,3,6,1,2].
\]

Thus label i is sent to r(i), and the relabeled tuple is r g_j r^-1. It is

\[
 T=\bigl((1\,5)(2\,3),\ (1\,3)(2\,6),\
          (0\,3)(4\,6),\ (3\,4)(6\,7)\bigr).
\]

For adjacent involutions define the right Hurwitz action

\[
 (\ldots,x,y,\ldots)\cdot\sigma_i
   =(\ldots,y,yxy,\ldots),\qquad
 (\ldots,x,y,\ldots)\cdot\sigma_i^{-1}
   =(\ldots,xyx,x,\ldots).
\]

Apply the following letters successively in written order:

\[
 w=\sigma_3\,\sigma_2^{-1}\,\sigma_1\,\sigma_3\,\sigma_2\,
       \sigma_1^{-1}\,\sigma_2\,\sigma_2\,\sigma_3\,\sigma_3.
\]

Then the exact equality is

\[
 T\cdot w
 =\bigl((0\,1)(2\,3),\ (1\,4)(2\,3),\
          (0\,5)(1\,6),\ (1\,2)(6\,7)\bigr),
\]

the literal b=4 knot tuple. Both tuples have product
(0 5 1 2 6 7 4)(3). The certificate supplies the sheet relabeling, all ten
intermediate Hurwitz operations, and the final equality.

The result is braid equivalence of branch cycles, with infinity distinguished.
It does not assert equality of branch values or isomorphism of the individual
quadratic-field and rational-coefficient covers. Since the knot tuple generates
A8, the same certificate also gives a direct geometric A8 proof for this model,
using the retained knot-generation theorem.

## The exact continuation argument

The [checker](verify_degree_eight_braid_word.py) uses rational complex arithmetic
throughout verification. Newton iterations are used only by the optional
certificate generator to propose root centers. Default replay neither imports
a numerical root solver nor trusts the proposed permutations.

The coefficient e is enclosed using

\[
 m^2<329\,2^{160}<(m+1)^2,\qquad
 m=\lfloor 2^{80}\sqrt{329}\rfloor.
\]

The midpoint E of the resulting imaginary coefficient interval is rational
complex, and its error bound epsilon is rational. Every root center is an
integer pair divided by 2^70.

Here is the finite inequality checked for each root and each parameter segment.
Write the midpoint-coefficient Taylor expansion at center zeta and the start
value s_0 as

\[
 F_{s_0}(\zeta+w)=\sum_{k=0}^8 b_k w^k,
 \qquad F_s(z)=z^8+z^6+ez^2-sz+a.
\]

For complex rational v put
||v||_1=|Re(v)|+|Im(v)| and
||v||_infinity=max(|Re(v)|,|Im(v)|).
Let delta=||s_1-s_0||_1. The coefficient uncertainties over the entire closed
segment s_0 to s_1 are bounded by

\[
 \epsilon_0=\epsilon\|\zeta^2\|_1+\delta\|\zeta\|_1,\qquad
 \epsilon_1=2\epsilon\|\zeta\|_1+\delta,\qquad
 \epsilon_2=\epsilon,\qquad
 \epsilon_k=0\quad(k\ge3).
\]

For a retained dyadic radius rho, the checker verifies the strict rational
inequality

\[
 \|b_0\|_1+\epsilon_0+
 \sum_{k=2}^8(\|b_k\|_1+\epsilon_k)\rho^k
 <(\|b_1\|_\infty-\epsilon_1)\rho.                  \tag{1}
\]

On |w|=rho the right side bounds the linear term below, and the left side
bounds all remaining terms above, uniformly for every parameter value on the
segment. Rouche's theorem therefore gives exactly one root in the disk at
every point of the segment. All eight tracking disks are pairwise disjoint,
so these are all roots and none collide on the path.

At each endpoint, a disk of radius 2^-40 is independently certified to contain
one root. It is strictly contained in the corresponding tracking disk.
The initial endpoint disk is also contained there. Uniqueness identifies the
two endpoint anchors with the same continued root. At loop closure, each
endpoint anchor is contained in exactly one of the original base disks of
radius 2^-39. This proves the recorded permutation.

There are 6,055 certified parameter segments in total, respectively
1,755, 1,575, 1,444 and 1,281 for the four loops. This is a covering proof of the
entire paths, not a sample of roots at finitely many parameter values.

The four finite branch values are the roots of the exact quartic Q from the
[model identity](DEGREE_EIGHT_ALTERNATING_PENCIL.md).
The same Taylor/Rouche test encloses one Q-root in the radius-1/1000 disk about
each q_j. These disks lie inside the squares and exhaust the quartic's roots.
Exact interval clipping checks that each spoke avoids the other branch squares.
These checks certify the geometric meaning of the four loops in addition to
the root labels along them.

The imported complex-analysis theorem is the standard
[NIST DLMF §1.10(iv), Rouche's theorem](https://dlmf.nist.gov/1.10.iv).
In the polynomial setting used here, it follows by deforming the linear term
to the polynomial while maintaining nonvanishing on the circle; its winding
number stays one, and the argument principle counts the roots. The
[source ledger](degree_eight_braid_word_source_ledger.json) records the exact
inspected source and the import boundary. No formal proof of complex analysis
is claimed.

## Reproduction and verification scope

From the repository root, run

```sh
python3 research/arithmetic/verify_degree_eight_braid_word.py
```

The verifier reads the
[retained exact certificate](degree_eight_braid_word_certificate.json.gz).
The [receipt](degree_eight_braid_word_verification.json) records the root
permutations, relabeling, braid word and intermediate states, segment count,
source hashes, and software version. It was generated with the additional
option
`--write-receipt research/arithmetic/degree_eight_braid_word_verification.json`.

The optional --generate mode uses mpmath to propose the root centers and then
runs the same exact verifier. Replaying the retained certificate requires only
Python's standard library. Both an inadequately small root disk and a truncated
braid word are rejected as negative controls.

This is an internally executed exact analytic-continuation and finite group
certificate. It is not formal verification, independent external replay,
external human review, or a priority claim.
