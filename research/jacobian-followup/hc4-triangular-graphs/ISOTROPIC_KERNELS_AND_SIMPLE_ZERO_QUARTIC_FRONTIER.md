# Isotropic kernels leave only the simple-zero quartic frontier

15 November 2026. Claims **HC4MYJO0/1**. In the affine
quadratic-composite family h=q_A+H(q_B+2v^T A x+c),
every pencil whose kernel contains an isotropic direction
now reduces its universal inverse question to the same
four fixed source-orbit tests. Orthogonality of that
direction to v is not required.

In particular the entire cubic coordinate-inner family
from JN has exactly those four universal tests. Beyond
them, an additional affine mechanism must have a simple
zero eigenvalue with nondegenerate eigenline and nonzero
offset there. Its only remaining generating form is

    H(t)+q_A3(z)+ell(z)+(t-q_B3(z)-beta)^2/(4b^2),
    b!=0,  q_A3,q_B3 nondegenerate ternary quadratics.    (1)

This quartic-plus-univariate class is not resolved here.
Neither h0,h1, the two actual old endpoints nor
unrestricted HC4/JC2 is solved.

## 1. Statements and retained inputs

Work over C in four variables, with A=A^T invertible,
B=B^T, T=A^(-1)B and deg H=d>=2.

**JO0.** For a!=0 and arbitrary b,c, the generic plane curve

    F(t,z)=H(t)+(zt-z^3-2bz^2-cz)/a=kappa               (2)

has finite degree-three t projection, exactly 2d simple
finite ramifications, r=gcd(3,d) actual infinity points,
and infinity ramification 3-r. If it has c_C components
and total compact genus G_C, then

    G_C=d-(3+r)/2+c_C >= d-(1+r)/2 > 0.                (3)

No connectedness theorem is assumed. It is the actual
plane potential q_0+H(g_0) for a marked zero Jordan
block of length two, by the global source coordinate
t=g_0.

**JO1.** If ker T contains any nonzero A-isotropic vector,
the universal full-source inverse question reduces
exactly to the four old/fixed tests. In particular,
every polynomial of the full cubic family

    H(t)+(x/a)t+q(z1,z2)-(x/a)p(x,z1,z2),
    a!=0, q nondegenerate, deg p<=2,                   (4)

is covered by that four-test reduction, with all outer
and source degrees unrestricted. The four tests are
not asserted to be solved.

Outside those tests, the prior closures and JO1 leave
only (1), with T of rank three, zero a simple eigenvalue,
and b nonzero. Both ternary forms in (1) are nondegenerate,
their relative operator may have any Jordan form,
ell is arbitrary linear and beta arbitrary constant.
The normal form is a literal global source equivalence,
not native Hessian invariance.

Retain [HG0](GLOBAL_VALUE_ORBITS.md): under unit Hessian
and gradient nonautomorphy the generic whole profile is
L^3-NL, and native non-Morse or positive-dimensional
critical schemes exclude every unit-Hessian source member.
[HC0](HODGE_TRACE_CLOSURE.md) supplies the additive/product
invariant; [IY0](CONTRACTIBLE_RADICAL_CURVES_AND_PENCIL_DEFECT.md)
retains characteristic-zero Riemann--Hurwitz.

Use the proved constant primary frames in
[JD0](NILPOTENT_QUADRATIC_PENCILS_AND_ENDPOINT_ORBITS.md),
the cyclic moment calculations and nilpotent endpoint
reduction in [JF1/2](AFFINE_QUADRATIC_PENCIL_RESIDUES_AND_ENDPOINTS.md),
the complete finite critical-mass test
[JH0](SIMPLE_STATIONARY_POLES_AND_NONCYCLIC_KRYLOV_CLOSURE.md),
and the literal Witt/eigenvector fibre strata from
[JI](AFFINE_POLAR_CURVES_AND_COMPLETE_INVERTIBLE_SEMISIMPLE_CLOSURE.md)
and [JL](JORDAN_POLAR_CURVES_AND_ALL_INVERTIBLE_PENCIL_CLOSURE.md).
[JM](ZERO_EIGENVECTOR_SHEARS_AND_RANK_ONE_PENCIL_CLASSIFICATION.md)
and [JN](COMPATIBLE_INNER_CENTERS_AND_QUADRATIC_COORDINATE_FRONTIER.md)
already cover rank-at-most-one, compatible centers and
admissible null directions by the same four tests.
[HP0](QUADRATIC_ISOMETRY_SATURATION.md) and
[HM0](AMBIENT_QUADRIC_GRAPHS.md) retain actual ambient transfers.

Write L=uv, U=u+v, epsilon=d mod 2,

    a_d=(d+epsilon-2)/2, rho=2-epsilon,
    M_d=L-a_d U+1-rho,
    Q_d=L^2+(d-1)L,
    V_d=L^3-L^2+L M_d.

These are the complete generic double-curve, split
rank-two-quadric and rank-three-quadric profiles from JL.
Only explicit constant quadratic frames and product
strata are used, not arbitrary finite-cover E multiplicativity.

## 2. The marked zero-Jordan plane has positive genus (JO0)

Use the zero-primary frame

    q_0=2zw,       g_0=z^2+2a w+2b z+c, a!=0.

The polynomial coordinate t=g_0 has inverse
w=(t-z^2-2bz-c)/(2a). The potential q_0+H(g_0)
is exactly (2), with the same global source inverse.
Its level equation is monic cubic in z:

    z^3+2bz^2+(c-t)z+a(kappa-H(t))=0.                  (5)

Thus t is finite free of degree three. Its polar is

    F_z=(t-3z^2-4bz-c)/a=0,
    t=3z^2+4bz+c.

On the polar,

    F=H(3z^2+4bz+c)+(2z^3+2bz^2)/a.

This polynomial has exact degree 2d: 2d>=4 dominates
the degree-three remainder. Generic kappa therefore
meets the smooth parabola transversely in exactly 2d
points. They give all the simple affine ramifications
on the smooth generic level curve.

At infinity set t=tau^(-3), z=tau^(-d)Z.
The leading equation is Z^3-a h_d=0, with three
distinct nonzero simple roots. Every other term in
(5) has strictly smaller pole order for d>=2.
The deck transformation tau->omega tau permutes
these roots by the d-th power of a primitive cube
root omega. It has r=gcd(3,d) orbits, each with
ramification index 3/r. The total degree three
accounts for every branch; no extra infinity point
is suppressed.

Applying Riemann--Hurwitz component by component gives
2G_C-2c_C=-6+2d+(3-r), which is (3).
If r=1 the lower bound is d-1>=1; if r=3 it is
d-2>=1. All affine coefficients are unrestricted.

## 3. A binary-cone trace with a quadratic line is inverse-closed

This lemma will handle a quadratic infinity resonance.
Let b(z1,z2) be a nondegenerate binary quadratic and
a(z)=q(z)+ell(z), with q any quadratic and ell linear.
Consider

    f=s^2+a(z)+tau b(z).                               (6)

The cone b=0 is two distinct lines through zero.
On line i write a=a_i r^2+b_i r. Let j count lines
with a_i!=0, and k count lines with a_i=b_i=0.
For generic kappa!=0, the curve
s^2+a_i r^2+b_i r=kappa has profile L-1, L or 2L
according as a_i!=0, a_i=0,b_i!=0, or both vanish.
The two curves meet at the two points z=0,s^2=kappa.
Hence their union C has

    E(C)=(2+k)L-j-2.

Off b=0, solve tau with s free. On b=0 keep C and
the free tau line. These disjoint strata give

    E(f=kappa)=L^3+k L^2-(j+1)L.                      (7)

If k>0 the L^2 coefficient contradicts HG0 under
hypothetical nonautomorphy. If k=0, a full Morse
critical set has N<=j. Indeed s=0 and b=0.
If ell=0, z=0 leaves a whole critical tau line,
already forbidding unit Hessian. Otherwise z=0 is
not critical. At a nonzero point of cone line i,
the tangential equation is 2a_i r+b_i=0.
It has at most one solution when a_i!=0 and no
solution when a_i=0,b_i!=0. The normal derivative
then fixes tau uniquely because db is nonzero.
Thus N<=j, whereas (7) would require N=j+1.

Every unit-Hessian full-source member of (6) therefore
has automorphic gradient. This is inverse closure,
not a general cancellation or source-stability assertion.

## 4. Exhausting a forced isotropic kernel

After the prior closures, an additional affine mechanism
has rank T two or three, is nonnilpotent, has
v not in im T and no admissible isotropic zero vector.
Suppose nevertheless ker T contains an isotropic vector.

Split orthogonally into the zero primary space and the
nonzero primary space using the JD frames. The zero
primary dimension is two or three; dimension four is
nilpotent and dimension one has nondegenerate eigenline.
If a zero Jordan chain of length at least two is unmarked,
its last isotropic eigenvector is admissible, contrary
to the reduction. Thus that head is marked.

For zero primary dimension two, its form is either
one marked length-two zero block or a semisimple
zero two-plane. In the latter case the zero-space
offset must be nonisotropic: a zero or isotropic
offset would itself admit an orthogonal isotropic
kernel direction. For zero primary dimension three,
the only remaining possibilities are a marked
length-three zero block or a marked length-two zero
block plus a zero line. Three semisimple zero lines
already give rank one and are covered by JM.

It therefore suffices to handle marked zero chains
and a nonisotropically marked semisimple zero two-plane.
This fixed primary split has no dependence on the
outer or source degree.

## 5. Marked zero chains

On a marked length-two zero block, the polynomial part
of the stationary inner function has leading term
3a^2 s^2, a!=0. On a marked length-three block its
leading term is -4a^2 s^3. The top-moment calculation
is exactly the retained JF calculation; other zero
lines only add terms of degree at most one.

Let l>=2 be that degree and m=d-1>=1. At infinity
F=s-H'(R) has degree ml>=2 and the nonzero-primary
determinant P has degree at least one. Thus W=1/(PF)
has decay of order at least three and F is nonzero.

If the nonzero primary dimension is one, every finite
root is simple: it is blocked with a double R pole
or a compatible simple root counted in full by JH0.
The complete critical-mass certificate excludes unit Hessian.
This covers the marked length-three case and a marked
length-two case with an extra zero line.

If the nonzero primary dimension is two, the zero block
is precisely the marked length-two block. The finite
certificate also works for two distinct simple eigenvalues,
a scalar primary with nonisotropic nonzero offset, and
a Jordan primary with marked head: the active R orders
are two or three and determinant multiplicity at most two.
A scalar primary with zero offset gives a full critical
quadric and already forbids unit Hessian.

There remain only a nonzero isotropic scalar offset
or an unmarked nonzero Jordan head. Let f0 be the
zero-primary plane q_0+H(g_0) from Section 2, whose
total genus is positive for every d>=2.

For the isotropic scalar offset, its full-eigenspace
Witt pair gives the exact JI whole-fibre profile

    E=L^3-Q_d+L E(f0=kappa).                           (8)

The subtracted binary quadric is nondegenerate:
mu A_0-B_0 is invertible since B_0 is zero-nilpotent.
For the unmarked nonzero Jordan offset, translate
only that invertible primary to center. Its affine
force is along the first chain coordinate. Cutting
its eigenvector as in JL gives rank three and rank
two quadrics, hence

    E=V_d-Q_d+L E(f0=kappa).                           (9)

Both identities retain the entire closed cut stratum.
Their curve genus enters with a negative sign, and
V_d can only add another nonnegative genus. Therefore
neither profile can be HG0's L^3-NL.

Every marked zero-chain case is thus inverse-closed,
in all affine coefficients and outer/source degrees.
No unaccounted multiple determinant residue is identified
with finite critical mass.

## 6. A nonisotropic offset on a semisimple zero two-plane

Normalize that plane to s^2+w^2 with offset b partial_w,
b!=0. Translate the invertible nonzero two-primary
quadratic to center. Up to a constant the potential is

    h=s^2+w^2+q_A2(z)+ell(z)
                         +H(q_B2(z)+2bw+beta),         (10)

where both binary forms are nondegenerate and ell linear.
Put eta=b^2, m=d-1 and gamma=H'(beta).
The regular stationary inner function at infinity is

    R(s0)=-2eta s0+beta+O(s0^(-2)),    deg P=2.

The parameter here is s0, independent of the source s.
For d>=3 the infinity decay is proper. If d=2,
write H'(r)=u r+v. The leading coefficient is
1+2u eta; if it vanishes the constant is -gamma.
Unless both vanish, W decays at least as s0^(-2).

All the finite good/bad alternatives of Section 5
apply to the nonzero two-primary space. The good
ones therefore satisfy JH0 except for the double
quadratic infinity resonance, treated below.

For the two finite bad types use the actual profiles
(8),(9), now with

    f0=s^2+w^2+H(2bw+beta).

For d>=3 its generic curve has positive genus a_d.
Thus (8),(9) again rule out nonautomorphy.
For a quadratic outer with nonzero leading coefficient
1+2u eta, the curve has profile L-1. For a leading
resonance u=-1/(2eta) but gamma!=0 it has profile L.
The resulting exact profiles and complete Morse-count
bounds are:

| Nonzero primary finite type | Outer case | Whole E | N bound |
|---|---|---|---|
| Nonzero isotropic scalar offset | nonresonant quadratic | L^3-2L | <=1 |
| Unmarked nonzero Jordan head | nonresonant quadratic | L^3-3L | <=2 |
| Nonzero isotropic scalar offset | leading resonance, gamma!=0 | L^3-L | 0 |
| Unmarked nonzero Jordan head | leading resonance, gamma!=0 | L^3-2L | <=1 |

For completeness the bounds retain all singular groups.
In both bad types the nonzero-primary contribution to R
is constant, so globally R=beta-2eta s0.
The nonresonant regular equation is linear, giving
at most one regular point. At leading resonance with
gamma!=0, F=-gamma has no regular root.
The isotropic scalar singular equation is inconsistent.
For the Jordan type, singular stationarity is a line;
its inner value is affine-linear along it or constant.
A quadratic H supplies at most one point or a forbidden
whole critical line. There are no further parameters.
Each table row contradicts its HG0 coefficient.

Finally suppose both infinity coefficients vanish:

    H(r)=H(beta)-(r-beta)^2/(4b^2).

The global source shear

    w=-b tau-q_B2(z)/(4b)

and its polynomial inverse give

    h=s^2+q_A2(z)+ell(z)+tau q_B2(z)+constant.           (11)

The binary-cone trace lemma of Section 3 applies to
every nonzero two-primary Jordan form and affine offset.
It closes this double resonance without omitting
its exceptional cone or stationary points.

This completes the semisimple zero two-plane case.

## 7. The entire isotropic-kernel and cubic families

If the kernel vector was admissible, JN already gives
the four-test reduction. If a compatible center or
rank-at-most-one/nilpotent case occurred, the retained
JN/JM/JF conclusions apply. The remaining forced
isotropic cases are exactly Sections 5 and 6 and
are inverse-closed. This proves the sufficiency of
the four tests for every isotropic-kernel pencil.
All four are already represented in the JN/JM examples
with isotropic kernel, so the reduction is exact.

The full cubic family (4) is the literal global source
normal form of

    2xy+q(z)+H(2a y+p(x,z)).

The vector partial_y is a nonzero isotropic kernel
vector of its inner quadratic. The coordinate t=2a y+p
has polynomial inverse y=(t-p)/(2a). Hence the whole
cubic family is covered, with no degree bound on H
or any subsequent source automorphism.

The four tests actually occur even in this cubic family.
Normalize a=1 and q=z^2+w^2. For p=0,H=t^2 the
potential is a nondegenerate quadratic.
For p=2w,H=-t^2/4+t/2 it is a coordinate source orbit.
For the latter fact and for h0,h1, the constant frame

    X=y,  Y=x-y/2-w,  W=y+w,  Z=z

satisfies 2XY+Z^2+W^2=2xy+z^2+w^2.
The old coordinate example with inner 2W becomes
2y+2w. The actual JM rank-one resonant graph with
inner Z^2+2c Z+2W+c^2 becomes
2y+p(z,w), p=z^2+2cz+2w+c^2.
At c=0 or c=-1/2 it is respectively h0 or h1's
full source orbit. Every map here has an explicit
polynomial inverse. Thus the cubic universal question
has precisely the four retained tests, not an
additional degree or coefficient frontier.

## 8. Only the nondegenerate simple-zero quartic family remains

An additional affine mechanism beyond the four tests
cannot have invertible T, an inner critical point,
rank at most one, an admissible null direction or
even an isotropic vector anywhere in ker T.

Over C a nonzero quadratic space with no nonzero
isotropic vector is a nondegenerate line. Therefore
ker T is such a line. Its orthogonal complement
is T-stable and invertible; zero is a simple eigenvalue
and rank T=3. The zero-space offset is nonzero,
since otherwise v belongs to im T and JN applies.

Normalize that line as w with q_A=w^2+q_A3(z)
and zero-space offset b partial_w, b!=0. Translate
the invertible ternary inner quadratic to center.
The potential becomes

    w^2+q_A3(z)+ell(z)+H(q_B3(z)+2bw+beta)+constant,

with both ternary forms nondegenerate.
The actual global coordinate t=q_B3+2bw+beta
has polynomial inverse w=(t-q_B3-beta)/(2b).
This gives exactly (1), with constant nonzero Jacobian.

It retains every source orbit and every outer degree.
It is not native Hessian invariance or a closure
of the remaining quartic class. The old endpoints,
h0,h1 and unrestricted HC4/JC2 remain open.

We rechecked the retained primary
[Riemann--Hurwitz foundation](https://stacks.math.columbia.edu/tag/0C1B)
and [Soprunov residue context](https://arxiv.org/abs/math/0511684v2).
No new unproved external theorem, stable cancellation
or arbitrary conic-cover multiplicativity is imported.
