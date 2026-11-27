# The reflected cubic field of a coherent norm triangle

A coherent cubic prime triangle determines a binary cubic whose
reflected Kummer class is exactly the difference of two edge classes.
The cubic is a line section of a diagonal plane cubic determined by
the coefficient classes. In the common cube-colour case this becomes
a line on the Fermat cubic and gives a cubic field of fundamental
discriminant when the three normalized coefficients are \(1\).

## The exact reflected Kummer class

Fix an imaginary quadratic field \(K=\mathbb Q(\sqrt{-d})\),
\(d>0\) squarefree and \(d\ne3\), and a coherent triangle with distinct
odd prime vertices as in the
[three-coefficient reconstruction](ABC_COHERENT_NORM_PENCILS.md#exact-reconstruction-from-three-signed-coefficients).
Write \((v_{12},v_{23},v_{13})=\delta(x,y,z)\) for its signed
normalized coefficients. Put \(E=d\delta^2\),
\[
   A=p_1^3y^2,\qquad B=p_2^3z^2,\qquad C=p_3^3x^2,\qquad
   H=2(AB+AC+BC)-(A^2+B^2+C^2)=E x^2y^2z^2 .
\]
The positive \(A,B,C\) are pairwise unequal: at \(p_i\) the entry
corresponding to \(p_i\) has odd valuation while the other two
entries have even valuation. Let \(\alpha_{ij}=(u_{ij}+v_{ij}\sqrt{-d})/2\).

**Proposition RF1.** The binary cubic
\[
       G(U,V)=BC U^3+AC V^3-AB(U+V)^3                          \tag{1}
\]
has discriminant \(27(ABC)^2H\). Let \(N\) be its splitting field
and \(L=\mathbb Q(\sqrt{\operatorname{sf}(3d)})\). Then \(N\) contains
\(L\), and
\[
       N(\zeta_3)=K(\zeta_3)
                  \bigl(\sqrt[3]{\,\overline{\alpha}_{12}\alpha_{23}}\bigr).
                                                                  \tag{2}
\]
If the three prime ideals have relative labels \(c_i\) in the common
three-torsion coset, \(G\) is irreducible over \(\mathbb Q\) exactly
when \(c_1+c_2+c_3\ne0\). In that case \(N/L\) is cyclic of degree
three and unramified outside \(3\). For fixed \(K\), two irreducible
triangles give the same cubic field up to isomorphism exactly when
their sums span the same line in \(\mathrm{Cl}(K)[3]\).

For the algebra of (2), write the coefficients of \(G(T,1)\) as
\(g_0,g_1,g_2,g_3\) in descending order. Then \(g_0=B(C-A)\ne0\),
\(g_1=g_2=-3AB,\ g_3=A(C-B)\). In the cubic formula put
\[
 \begin{split}
  D_0&=g_1^2-3g_0g_2=9AB^2C,\\
  D_1&=2g_1^3-9g_0g_1g_2+27g_0^2g_3
       =27AB^2C\bigl((C-A)^2-B(A+C)\bigr).
 \end{split}
\]
Expansion gives
\(D_1^2-4D_0^3=-(27AB^2C)^2(C-A)^2H\), and the identity
\(D_1^2-4D_0^3=-27g_0^2\operatorname{disc}(G)\) gives the
discriminant assertion. Set \(s=xyz\sqrt{-E}\), so \(s^2=-H\).
The cube root in Cardano's formula has cube
\[
       \frac{D_1+\sqrt{D_1^2-4D_0^3}}{2}
       =\frac{27AB^2C}{2}\bigl((C-A)^2-B(A+C)+(C-A)s\bigr).
\]
Over \(K(\zeta_3)\) this cube root generates the splitting field:
the three roots are linear combinations of it and its inverse,
and, conversely, the root combinations with coefficients
\(1,\zeta_3,\zeta_3^2\) recover it.

Define \(\beta_{12}=(A+B-C+s)/2=yz\alpha_{12}\) and
\(\beta_{23}=(B+C-A+s)/2=xz\alpha_{23}\).
Multiplication shows that the displayed Cardano cube is
\(-27AB^2C\overline{\beta}_{12}\beta_{23}\).
The factor \(27AB^2C\,xyz^2=(3p_1p_2^2p_3\,xyz^2)^3\)
is itself a cube. This proves (2) including the cancellation
of all the normalized coefficients.

In \(K\) the cube class of the virtual unit
\(\overline\alpha_{12}\alpha_{23}\) corresponds to the class of
\((\overline{\mathfrak p}_2\mathfrak p_3)/
(\overline{\mathfrak p}_1\mathfrak p_2)\), namely
\(c_3-2c_2+c_1=c_1+c_2+c_3\). The virtual-unit to class-group
map is an isomorphism in this imaginary quadratic field:
its kernel is the units modulo cubes, trivial for \(d\ne3\).
Changing the base label alters the sum by a multiple of three.
The restriction of cube classes to \(K(\zeta_3)\) is injective by
the norm, since the degree is two. Thus reducibility, equivalently
that the splitting field has order at most two, is exactly the
zero class condition. In the irreducible case \(N\) is \(S_3\)
because its positive discriminant has the nonsquare class \(3d\);
in particular the cubic field is totally real. At all primes
away from \(3\), (2) is unramified over \(K(\zeta_3)\) because
its radicand has valuations divisible by three. Ramification
in \(N/L\) would have index three, which cannot disappear after
the quadratic base extension \(L(\zeta_3)/L\).
Finally Kummer theory determines \(N(\zeta_3)\) from the line
in \(\mathrm{Cl}(K)[3]\), and conversely. Its group over
\(\mathbb Q\) is \(S_3\times C_2\), with unique central subgroup
of order two. The fixed field of that subgroup is \(N\).

This is the virtual-unit form of cubic reflection used in
[Cohen–Thorne, Section 3 and Proposition 4.1](https://arxiv.org/pdf/1301.3563).
Here the cubic field is designated by the coherent triangle. The
formula identifies its class direction without choosing any
further class-group generators.

The [cofactor cubic](ABC_COHERENT_NORM_COFACTOR_CUBIC.md#exact-cubic-and-quadratic-forms)
already gives an integral form with the same reflected line.
In its coordinates \(n\in p\times\Lambda\), set
\((U,V)=(p_1n_1,p_2n_2)\), so \(U+V=-p_3n_3\).
Then \(G(U,V)/((p_1p_2p_3 xyz)^3)\) is minus that integral
cubic, by its displayed cofactor formula. This identifies
the plane-cubic line section with the cofactor form and its
marked Hessian without adding a second integral model.

## Line sections and a global restriction on the labels

**Corollary RF2.** Fix the normalized coefficient cube classes
of \(x,y,z\), and choose positive integral representatives
\(x_0,y_0,z_0\) that are cube-free. There are nonzero signed integers
\(b_x,b_y,b_z\) with \(x=x_0 b_x^3,\ y=y_0 b_y^3,\ z=z_0 b_z^3\).
The triangle gives a rational line
\[
    p_1b_y^2 U_1+p_2b_z^2 U_2+p_3b_x^2 U_3=0
\]
on the fixed diagonal plane cubic
\[
                 \frac{U_1^3}{y_0^2}
                  +\frac{U_2^3}{z_0^2}
                  +\frac{U_3^3}{x_0^2}=0.                       \tag{3}
\]
Its intersection has a rational point exactly when the sum of
the three labels vanishes. Otherwise it is a degree-three
closed point defined over the reflected cubic field of RF1.
For fixed \(d\) and the three representatives the map to
projective lines has at most eight choices of signed witnesses
per line.

Indeed divide (1) by \(ABC\) and use \(W=-U-V\) to obtain the line
\(U+V+W=0\) on \(U^3/A+V^3/B+W^3/C=0\). The three changes of
scale yield (3). The line coefficients are positive and
primitive. If a prime divided all three, then, since the
\(b\)'s are pairwise coprime, it would have to divide two
different vertex primes. Each coefficient recovers its vertex
prime as its positive squarefree part, then the corresponding
\(|b|\). The signs give at most eight options. For fixed signed
coefficients the reconstruction gives the traces and common
scale.

When all three unnormalized coefficients have the same cube
class, normalization makes \(x,y,z\) signed cubes. Equation (3)
is the Fermat cubic \(U_1^3+U_2^3+U_3^3=0\). Its rational points
have one coordinate zero and two opposite coordinates, by
Fermat for exponent three. The line cannot meet these points
because the positive line coefficients have distinct prime
square classes. Thus in any such triangle
\(c_1+c_2+c_3\ne0\). In particular the three distinct labels
cannot be the three points of an affine line over \(\mathbb F_3\).
This is a restriction on joint label and coefficient data
beyond the marginal counts for individual edges. More generally
any coefficient profile for which (3) has no rational point
excludes the zero label sum. For \((x_0,y_0,z_0)=(5,3,4)\),
equation (3) is isomorphic to the classical Selmer cubic
\(3u^3+4v^3+5w^3=0\), which is everywhere locally soluble and
has no rational point. This failure is visible through the
label sum. [Van Luijk, introduction](https://arxiv.org/pdf/0708.2752)

For example, the normalized profile \((5,3,4)\) is genuinely
realized on primes \((167,173,179)\). The reconstruction gives
\(E=3762067939459=19\cdot29\cdot6827709509\), which is squarefree
and prime to the three vertices. The reflected cubic field in
this case has discriminant \(27E\), computed from (1) by
PARI/GP `nfdisc`. In contrast, over \(d=43\) choose the
principal oriented primes generated by
\((1+\sqrt{-43})/2,(9-\sqrt{-43})/2,2+\sqrt{-43}\) over
\(11,31,47\), respectively. The products of the conjugate
initial generator and final generator, cubed, give the signed
coefficient triple \((260,-9152,-1508)\).
The normalized triple is \((5,-176,-29)\).
Here (1) has affine root \(-3872/2697\), in agreement with the
zero label sum from principal primes.

## The fundamental discriminant on the Fermat slice

**Proposition RF3.** Assume a triangle of RF1 has
\(x=y=z=1\) and \(p_i>3\). Its reflected cubic field \(M\) has
fundamental discriminant \(\operatorname{disc}(L)\).
If the three \(p_i\) agree modulo \(3\), then \(3\mid d\),
\(\operatorname{disc}(M)=d/3\), and \(3\) splits completely in
\(M\).

The line is \(p_1 X+p_2 Y+p_3 Z=0\) on the Fermat cubic, so
\(M\) is the field of a root of
\[
                  g(T)=p_3^3(T^3+1)-(p_1 T+p_2)^3.             \tag{4}
\]
This has discriminant \(27 p_3^6 E\), either directly or
by changing variables in (1). The nontriviality follows
from RF2 and ramification away from \(3\) from RF1.
If the three primes agree modulo \(3\), their cubes are
equal to the same sign modulo \(9\). Therefore \(E\equiv3\bmod9\).
The binary form of (4) divided by \(3\) has integral
coefficients, and its reduction modulo \(3\) is a nonzero
scalar times \(UV(U+V)\). The three simple projective
roots lift over \(\mathbb Q_3\), giving complete splitting.
Here \(d\equiv3\bmod4\) and \(d/3\equiv1\bmod4\), so
\(\operatorname{disc}(L)=d/3\).

If two congruence signs differ, reorder so \(p_1,p_2\)
have the same residue and \(p_3\) the other. Then
\(g(-1)=-(p_2-p_1)^3\) has valuation at least three at
\(3\), whereas
\(g'(-1)=3p_3^3-3p_1(p_2-p_1)^2\) has valuation one.
Hensel's lemma gives a rational local root, so the local
splitting extension has degree at most two. In either
case \(N/L\) is unramified at \(3\). Wherever \(L/\mathbb Q\)
is ramified, the decomposition group in \(S_3\) is of
order two because its inertia subgroup is a transposition
and is normal in the decomposition group. The local cubic
algebra there is a rational factor and the ramified quadratic
field. At the other primes the cubic algebra is unramified.
This gives \(\operatorname{disc}(M)=\operatorname{disc}(L)\).

For \((17,23,29)\) of the
[progression example](ABC_COHERENT_PRIME_PROGRESSION_NORMS.md#a-prime-progression-and-its-ideal-classes),
exact use of (4) gives \(\operatorname{disc}(M)=61894801=d/3\).
For the mixed residues \((23,29,31)\), \(d=1141198627\) is
squarefree, and the corresponding discriminant is \(3423595881=3d\).
Both are internally reproducible with `nfdisc`.

## A family of cubic points and the fixed-field count

**Corollary RF4.** At least \(\gg X^{1/3}/(\log X)^{B_0}\) distinct
totally real cubic fields of positive fundamental discriminant at
most \(X\), with \(3\) totally split, carry a nontorsion point on
the fixed Fermat elliptic curve \(U_1^3+U_2^3+U_3^3=0\).
They carry such a point from a rational line whose ordered three
coefficients are primes in arithmetic progression, and the
imaginary reflection field has two independent cubic ideal classes
carried by the same line's norm triangle.

Use [PP3](ABC_COHERENT_PRIME_PROGRESSION_NORMS.md#quantitative-count-of-genuine-fields)
and \(p_i\equiv p_j\bmod3\) from that construction. Proposition RF3
gives distinct discriminants \(d/3\). The line section is an
irreducible cubic point, hence determines its field. The fixed
elliptic curve has only finitely many torsion points of degree
at most three, by Northcott and the bounded difference between
the canonical and coordinate heights. Removing their finitely
many fields proves the assertion. Fermat for exponent three
gives rank zero over \(\mathbb Q\), since all rational points
are its three flexes. The count is thus a rank-growth result
with an associated coherent norm description.
[Lemke Oliver–Thorne, Theorem 1.1](https://arxiv.org/pdf/1810.04018)
give a general \(X^{1/3-\varepsilon}\) rank-growth bound in
cubic extensions for a rational elliptic curve, with prescribed
root number. The family here also records the positive
fundamental discriminants, complete splitting at \(3\), and
the cubic prime-line construction.

For \(x=y=z=1\) and fixed \(d\), the squarefree-part equation
\(E=d\delta^2\) over the projective prime directions is a
quadratic twist of the discriminant double cover for the Fermat
cubic. The latter adjoins a square root of \(27E\) on the
chart of (4), after removing the square factor \(p_3^6\). Over an algebraic closure this is the quotient of the
ordered collinear triples by the cyclic permutation, equivalently
\(C\times C\) modulo rotation of the triple.
[Van Luijk, Section 2 and Example 2.13](https://arxiv.org/pdf/0708.2752)
give this double-plane description for a smooth plane cubic.
The reflected class line identifies the field of a cubic point
above the direction. For general coefficients RF2 takes the
corresponding line section on a diagonal cubic. A bound for
points of the fixed genus-two progression slice suffices for PP3.
The full rank-one count would require a uniform count with the
diagonal cubic profiles, fields and marked line coefficients all
varying.
