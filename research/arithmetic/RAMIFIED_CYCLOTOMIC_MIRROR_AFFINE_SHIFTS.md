# A residual affine shift at a ramified cyclotomic place

Complete splitting of an odd prime-power
extension above the mirror prime can leave an
affine translation in the absolute local action.
For the ramified quadratic completion containing
cube roots of unity, this happens on lower
three-power directions. The largest affected
order gives the exact local discriminant.
A real norm family realizes the shift with a
global unit and an actual moving prime class of
full three-primary order.

Fix \(m=3^n\), and let \(K\) be a quadratic field
different from \(\mathbb Q(\zeta_3)\), with
\(K_3=\mathbb Q_3(\zeta_3)\). Form the diagonal
mirror \(E\) of \(K\mathbb Q(\zeta_m)\) as in
[AFM3](TWO_PRIME_MIRROR_INTERMEDIATE_FIELDS.md).
Let \(S_m(K)\) be the \(m\)-virtual classes
in \(K^\times/K^{\times m}\) which are local
\(m\)-th powers at the prime above three.

## The exact local affine action

**Theorem RAM1.** Scalar abelian extensions of
\(E\) of exponent dividing \(m\), everywhere
unramified and split above three, correspond
to the subgroups \(S\) of \(S_m(K)\) by the
same Kummer bijection as AFM3. Take a type
\[
             G=\operatorname{Hom}(S,\mu_m)
                         \simeq\prod_i C_{3^{a_i}} .
\]
Choose compatible cyclic generators of \(S\).
For a generator of order \(3^{a_i}\), choose
its unique lower class
\(\beta_i\bmod K^{\times3^{a_i}}\), whose
\(3^{n-a_i}\)-th power is the original
class. In the local quotient write
\[
              \beta_i=\zeta_3^{j_i} \pmod{K_3^{\times3^{a_i}}},
                  \qquad j_i\in\mathbb F_3.
 \tag{1}
\]
Here \(j_i=0\) if \(a_i=n\).
On the complement field, the local
Galois generator at three acts on \(G\)
by \(x\mapsto4x+t\) with
\(t_i\equiv j_i\bmod3\).
Put \(c=\max\{a_i:j_i\ne0\}\), with
\(c=0\) when all phases vanish.
For \(0\le s\le n-1\) the number of points
with orbit length at most \(3^s\) is
\[
 B_s=\begin{cases}
          0 &s<c,\\
          \displaystyle\prod_i3^{\min(a_i,s+1)} &s\ge c .
       \end{cases}
 \tag{2}
\]
With \(B_{-1}=0\), the local algebra
has \((B_s-B_{s-1})/3^s\) factors
\(\mathbb Q_3(\zeta_{3^{s+1}})^+\)
for each \(s\).
In particular \(c\) strictly increases
the discriminant exponent at three
for a fixed translation type whenever
it increases among the allowed values.

The global descent in AFM3 still
applies. The additional local point
is that restriction from \(K_3\) to
\(C_n=\mathbb Q_3(\zeta_{3^n})\)
is injective modulo \(m\)-th powers.
If \(x^m\in K_3\), a generator
\(g:\zeta_m\mapsto\zeta_m^4\) of
\(C_n/K_3\) has \(g(x)/x=\zeta_m^z\).
Multiplying the conjugate ratios
shows that \(z\) is killed modulo
\(3^n\) by
\(1+4+\cdots+4^{3^{n-1}-1}\),
of valuation \(n-1\).
Thus \(z\) is divisible by three.
Multiplying \(x\) by
\(\zeta_m^{-z/3}\) makes it invariant,
so the power already exists in
\(K_3\). This includes \(n=1\)
with trivial local extension.
The unramified valuation conditions
at all other finite places and
the central diagonal involution are
as in AFM3. When the quadratic
field ramifies at two, the
quadratic conductor contribution
outside three follows as in AFM5.

The lower class \(\beta_i\) exists
and is unique since \(K\) contains
no nontrivial three-power roots of
unity. Its \(3^{n-a_i}\)-th power
is a local \(m\)-th power, so the
lower class is represented by a
root of unity killed by
\(3^{n-a_i}\). This proves (1):
\(K_3\) contains \(\mu_3\) but
does not contain \(\mu_9\).
On choosing an \(a_i\)-level
radical locally, it has the form
\(u_i\zeta_{3^{a_i+1}}^{j_i}\)
with \(u_i\in K_3\), when \(j_i\)
is nonzero. These roots lie in
\(C_n\) since then \(a_i<n\).
The action of \(g\) changes this
root by the multiplier
\(\zeta_{3^{a_i}}^{j_i}\).
The local decomposition group
of \(E\) has order \(3^{n-1}\)
and cyclotomic generator \(4\):
the mirror is locally
\(\mathbb Q_3(\zeta_{3^n})^+\).
Because the relative extension
splits there, its local action
over the rational prime is
exactly this lifted generator.
These are the stated translations.
A change of local root choices adds
a multiple of three to \(t_i\).

For an integer \(h>0\), the
\(h\)-th power of the affine map
changes \(x\) by
\[
            (1+4+\cdots+4^{h-1})(3x+t).
\]
The first factor has three-adic
valuation \(v_3(h)\). Thus the
orbit length of \(x\) is the
order of \(3x+t\). On a shifted
coordinate, \(3x_i+t_i\) is a
unit and always has order
\(3^{a_i}\). On an unshifted
coordinate one can translate
to make \(t_i=0\).
This gives (2) and the local
field degrees. An orbit of size
\(3^s\) has local discriminant
exponent
\[
                   \frac{(2s+1)3^s-1}{2}.
 \tag{3}
\]
Indeed it is the fixed field
of the order-two tame subgroup
in the full local cyclotomic
field at level \(s+1\).
Dividing (3) by \(3^s\) gives a
strictly increasing function of
\(s\). Formula (2) is obtained
from the scalar orbit exponents
by replacing each \(s\) by
\(\max(s,c)\). Since some
scalar orbits have size one,
the total exponent increases
strictly with \(c\).
The one-coordinate cyclic
action is the high-divisibility
case of
[LD3](LOCAL_DICKSON_NORM_ONE_UNIT_TOWERS_AT_THREE.md).
Here the relative extension is
split at the full ambient level,
and several lower factors can
contribute to the absolute shift.

## The mass at a fixed wild discriminant

**Proposition RAM2.** For \(1\le j\le n\)
define \(\lambda_j:S_m(K)[3^j]\to\mathbb F_3\)
by taking the lower class and its phase
as in (1). Put \(s_j=|S_m(K)[3^j]|\),
\(s_0=1\), and \(z_j=|\ker\lambda_j|\).
For a bound \(c\le d\), set
\[
       h_j(d)=
       \begin{cases}s_j&j\le d,\\z_j&j>d.\end{cases}
\]
The number of subgroups of type
\(C_{3^a}\times C_{3^b}\), \(a\ge b\ge1\),
with \(c\le d\) is
\[
        \frac{(h_a(d)-s_{a-1})(h_b(d)-3s_{b-1})}
                       {|\operatorname{Aut}(C_{3^a}\times C_{3^b})|}.
 \tag{4}
\]
For a cyclic subgroup use
\((h_a(d)-s_{a-1})/\varphi(3^a)\).
Differences in \(d\) give the exact
wild-discriminant mass at fixed
quadratic resolvent and translation type.

An element of order less than \(3^j\)
has zero \(j\)-phase because its lower
class at level \(j\) is a cube of a
lower-level class. For \(i\ge j\),
the \(j\)-phase of \(3^{i-j}x\) equals
the \(i\)-phase of \(x\); the same
lower representative can be used.
Thus \(c\) in RAM1 is intrinsic:
it is the largest \(j\) for which
\(\lambda_j\) is nonzero on
\(S[3^j]\). A pair of generators is
allowed exactly when both phases
at the levels exceeding \(d\)
vanish. Choose the first
generator \(x\) with exact order
\(3^a\) in its allowed subgroup.
Its order-three line lies in
the image of the allowed
second-generator subgroup under
multiplication by \(3^{b-1}\).
For \(b>d\) this is because
\(3^{a-b}x\) has zero \(b\)-phase;
for \(b\le d\) it is automatic.
Each of the three points of
the line has \(s_{b-1}\)
preimages, since all
lower-order elements have
zero phase. Excluding that line
enforces independence and
gives (4) by the same injection
count as AFM3. Equation (3)
and the tame factors identify
the claimed discriminant
mass for the complement fields.

## A genuine moving prime and a shifted unit line

**Theorem RAM3.** For every fixed \(m=3^n\),
\(n\ge2\), there are fixed coefficients and
\(\gg_m X^{1/(2m)}/(\log X)^2\) real
quadratic fields \(|D|\le X\) with
\(K_3=\mathbb Q_3(\zeta_3)\) from
integral prime parameters \(u=P\).
On the same fields, for every
\(1\le b<n\), a genuine moving
prime class has three-primary
order \(m\), and \(S_m(K)\) contains
a subgroup of type \((n,b)\) with
shift \(c=b\).
The first direction comes from
the moving split prime and the
second from a global unit.
Both unramified dyadic types can
be imposed. The affine fields
of degree \(3^{n+b}\) have the
local factorization (2) and
supply \(\gg_m Y^{1/(m(3^{n+b}-1))}
/(\log Y)^2\) fields of their type
with absolute discriminant at
most \(Y\), for fixed coefficients.

Let \(a=3^n+1\) be even,
\(A=-a^m\), and put \(M=2n+5\).
Choose an odd squarefree integer
\(p_0>1\) coprime to \(a\) with
\[
             p_0\equiv 4a^m\pmod{3^{M-1}},
             \qquad q=3p_0+4a^m
                      \ \text{squarefree}.
\]
The squarefree pair exists by
the elementary sieve of ENJ4
applied to \(p_0,3p_0+4a^m\).
The local restrictions make the
two values prime to every
prime dividing the progression
modulus except possibly two,
where both are odd. At all
remaining primes at most two
classes modulo the prime square
are excluded. Both coefficients
are coprime to \(a\), and \(p_0,q\)
are coprime. Use \(p=3p_0\)
and \(B=-(p+q)/2\). Then
\(p=2A-B,\ q=-2A-B\) and
\[
             B/A\equiv14\pmod{3^M}.
\]

Use the integral prime \(P\) with
\(P\equiv a^{-1}\bmod3^M\).
In particular \(P\equiv1\bmod m\).
Put \(t=Au^m\) and
\(R^2=t^2+(B/A)t+1\) at \(u=P\).
Then \(t\equiv-1\bmod3^M\) and
\(R^2\equiv-3(t-1)^2\bmod3^M\).
Hence \(K_3=\mathbb Q_3(\pi)\)
for \(\pi^2=-3\), and choose the
branch \(R/(t-1)\equiv\pi
\bmod\pi^{2M-1}\).
Write \(\gamma_\pm=(t\pm1+R)/
(t\pm1-R)\). The plus ratio is
\(-1\) times a principal unit
of depth at least \(2M-1\)
in the \(\pi\)-valuation.
The minus ratio is
\((1+\pi)/(1-\pi)=\zeta_3\)
times a unit of the same depth.
Principal units of depth at
least \(2n+2\) are \(m\)-th powers
here, by the logarithm on
units of depth at least two.
Thus \(\gamma_+\) is a local
\(m\)-th power and
\(\gamma_-\in\zeta_3 K_3^{\times m}\).

At two choose \(P\) odd with
\(R^2\equiv4\) or \(20\bmod32\)
as in TP5–TP6. At each prime
dividing \(p_0\) choose
\(P\equiv a^{-1}\) so \(t+1\)
vanishes modulo the prime,
and at primes dividing \(q\)
choose \(P\equiv-a^{-1}\)
so \(t-1\) vanishes.
The norm identities give
valuation one there. At other
primes dividing \(a\) require
\(P\) to be a unit.
All these conditions are
compatible with \(P\equiv1
\bmod m\).

At \(P\) the two Kummer residue
frames for \(\gamma_+\) and
\(U=\gamma_+\gamma_-\) on the
branch \(R\) near one are
\(4/p\) and \(q/p\).
Their valuation vectors at
a fixed prime dividing \(p_0\)
and another dividing \(q\)
are unimodular. Both fixed
primes are away from three.
The resulting rank-two
Kummer field has maximal
abelian subfield
\(\mathbb Q(\zeta_m)\), since
inversion acts on both
translations. Chebotarev
can therefore give a basis
vector with \(\gamma_+\)
of full residue order and
\(U\) trivial, compatible
with the stated progressions.
In particular \(\gamma_+\)
is a global noncube.

At \(P\) the valuations of
the two ratios cancel in
\(U\). In fact the product
has the exact formula
\[
 U=\frac{(J+2AR)^2}{pq},\qquad
 J=2A^2P^m+B,\qquad J^2-4A^2R^2=pq.
\]
The auxiliary congruences
give \(pq\mid J\) since \(pq\)
is squarefree and
\(J=2At+B\) reduces to
\(-p\) over \(p\) and
\(-q\) over \(q\). Thus
\(\operatorname{Tr} U=4J^2/(pq)-2\)
is integral. This proves that
\(U\) is a global unit of
norm one. Its local phase
is nontrivial at three,
so it is a global noncube
and generates the image
of the real unit group
modulo \(m\)-th powers.
It is independent of
\(\gamma_+\) modulo cubes
since the latter is a
local cube at three and
a noncube at \(P\).
The two classes generate
\(C_m^2\) modulo \(m\)-th
powers by successively
dividing a putative relation
by three, since \(K\)
contains no nontrivial
three-power roots of unity. The divided ideal
of \(\gamma_+\) is the
ratio of the two split
primes above \(P\), whose
class is twice one of
their classes up to sign.
Independence from the unit
line gives exact order
\(m\) for its three-primary
part. Now \(U^{3^{n-b}}\)
is a local \(m\)-th power
at three and lies in
\(S_m(K)\). Its lower
class is \(U\), of phase
one, so together with
\(\gamma_+\) it yields
type \((n,b)\) and
shift \(c=b\).

The admissible primes \(P\asymp T\)
number \(\gg_m T/\log T\)
for fixed coefficients.
The discriminants are
\(O_m(T^{2m})\).
For any fixed \(K\), choose
the real branch with \(R/t\to1\)
as \(P\to\infty\). Then
\[
    U=\frac{(f_+f_-)^2}{pq P^{2m}}
                    \asymp_m P^{2m}
\]
with \(f_\pm=t\pm1+R\).
The positive norm-one units
of \(K\) in this range number
\(O_m(\log T)\), uniformly
in \(K\), since a positive
norm-one unit greater than one
has integral trace at least three.
The rational function \(U\)
on the fixed curve has finite
degree and is nonconstant,
so each such unit has boundedly
many parameters. This proves
the stated count of fields.
The complement field recovers
\(K\) from its mirror since
\(\ell=3\); the discriminant
has the fixed factor at three
given by RAM1 and the usual
outside-three power of \(D/3\).
This yields the stated count
by absolute field discriminant.

The [exact checker](verify_ramified_cyclotomic_mirror_affine_shifts.py)
replays the higher-rank orbit and
mass formulas and a genuine prime
coefficient and local unit fixture.
