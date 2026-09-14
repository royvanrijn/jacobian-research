# The same-branch orthogonal genus-zero stratum on Q80 is empty

On the fixed direct11952 alternate-Q80 MW17 parent, the full polynomial
genus-zero system

```
x_i^3+A*x_i+B=q*s_i^2,
deg x_i<=4, deg s_i<=5, q squarefree binary quadratic,
x_1 != x_2, gcd(q,Delta)=1
```

has **no rational solution with `q | (x_1-x_2)` and
`<P_1,P_2>=0`**, where `P_i=(x_i,s_i*sqrt(q))` and heights are measured
on the quadratic pullback. This is the entire `(k,j)=(2,1)` stratum below,
including higher contacts, degree drops and infinity. It does not require
the genus-zero cover to have a rational point.

The new finite calculation proves a stronger intermediate statement:
**none of the 49 full minimum-norm-twelve pencils has a rational singular
member**. Together with the retained norm-eight theorem, this excludes
every geometrically integral, rationally defined bisection of arithmetic
genus one and geometric genus zero on this Q80 surface.

The [subsequent branch-field theorem](Q80_RATIONAL_BISECTION_BRANCH_TORSION_2026-09-14.md)
also excludes the same-branch cross heights `+/-2`. Thus all three `k=2`
genus-zero cases are now closed. The later
[reciprocity theorem](Q80_SINGLE_BRANCH_RECIPROCITY_2026-09-14.md)
also closes all four `k=1` cases. The subsequent
[coefficient collision proof](Q80_GENUS_ZERO_COLLISION_CLOSURE_2026-09-14.md)
closes the five `k=0` cases, completing this polynomial genus-zero chart.

The remaining genus-one strata, other parents, higher-height sections,
the positive MW17 construction and a whole-family genus lower bound remain
**OPEN**. In particular, this is not a theorem that every rational bisection
on Q80 has positive genus: arithmetic-genus-zero bisections do exist.

## 1. The full genus-zero coefficient chart

Use binary forms of degrees `x_i:4`, `s_i:5`, `q:2`, `A:8`, `B:12`.
The previous repeated-common-quartic boundary has `D=e^2*q` and
`s_i=e*r_i`. Allowing all quintics s_i removes its extra requirement that
the node factors have a common linear factor. The present chart is broader.

Choose a coordinate where `h=x_1-x_2` has degree four and q has degree two.
With monic `g,d,u,v` and nonzero rational `c,lambda`, put

```
g=gcd(q,h),       k=deg g in {0,1,2},
q=c*lambda*g*d,   h=c*g*u*v,
u=gcd(h/c/g,s_1-s_2),      j=deg u in {0,...,4-k},
s_1=(u*a+v*b)/2,          s_2=(v*b-u*a)/2,
x_2=z,                   x_1=z+h.
```

The twelve degree allocations are

```
deg d=2-k, deg v=4-k-j,
deg a<=5-j, deg b<=1+k+j, deg z<=4.
```

The fixed-parent coefficient equations are

```
A=lambda*d*a*b-3*z^2-3*z*h-h^2,
B=q*(v*b-u*a)^2/4+2*z^3+3*z^2*h+z*h^2-lambda*z*d*a*b.
```

Require `g*d` squarefree, `gcd(d,u*v)=1`, `gcd(v,a)=1`, both s_i
nonzero, and `gcd(q,Delta)=1`. The literal scalar `c*lambda` is retained.
Subtraction of the two section identities gives the factorization exactly
as in [the quartic branch-stratum proof](COMMON_QUARTIC_BRANCH_STRATA_2026-09-14.md#1-a-complete-factor-presentation).
Conversely direct expansion gives both identities. There are 21 coefficient
parameters before one rescaling, and 22 equations for the fixed A,B; this
count alone proves no exclusion.

Seven distinct good rational sites suffice to choose infinity: the binary
form `q*h` has degree six. The previously verified sites `1,...,9` contain
such a set. At a site s not on `q*h`, use `t=s+1/v` with weighted transforms
`x_new=v^4*x`, `s_new=v^5*s`, `q_new=v^2*q`, `A_new=v^8*A`,
`B_new=v^12*B`. Thus the chart covers originally lower degrees and every
infinity contact. Here the letter s in the coordinate change is a scalar
site, distinct from the node polynomials s_i.

The pullback has `chi=4`, irreducible fibres and sections disjoint from zero.
The same local intersection proof as for squarefree quartic branch gives

```
height(P_i)=8,                  P_1.P_2=k+2*j,
<P_1,P_2>=4-k-2*j,
height(P_1+P_2)=24-2*k-4*j,
height(P_1-P_2)=8+2*k+4*j.
```

No ordinary-node or simple-contact assumption is made. These formulas use
the standard intersection height pairing; see
[Schütt--Shioda, Elliptic Surfaces, section11.8](https://arxiv.org/pdf/0907.0298).
The Gram determinant is at least48, so any admissible solution has two
independent anti-invariant directions modulo inherited MW. For `k=2`, the
three cross heights are `2,0,-2`. Only the middle case is closed here.

## 2. All 49 moving norm-twelve singular loci have no rational parameter

Reuse the complete full pencils and their integral models from
[the moving-contact theorem](Q80_ALL_SMOOTH_GENUS_ONE_BISECTIONS_2026-09-13.md).
For each trace `tau=(Nx/h^2,Ny/h^3)` its branch form is

```
g=v*t-u,       M=g*M0-v*c*h^2,       c=coefficient(M0,t^7),
Q(t;u,v)=(M^4-6*M^2*Nx*g^2-8*M*Ny*g^3
           -3*Nx^2*g^4-4*A*h^4*g^4)/h^6.
```

Here Q is binary of degree four in each variable pair. Define
`F(u,v)=disc_t Q(t;u,v)`, a homogeneous form of degree 24. If a member
has geometric genus zero, its branch quartic cannot be squarefree, so
`F(u,v)=0`. This necessary condition suffices for exclusion regardless
of whether the singularity is nodal, cuspidal or more degenerate.

For `Q=a*t^4+b*t^3+c*t^2+d*t+e`, the constructor uses

```
I=12*a*e-3*b*d+c^2,
J=72*a*c*e+9*b*c*d-27*a*d^2-27*b^2*e-2*c^3,
F=(4*I^3-J^2)/27.
```

The coefficients a,...,e are homogeneous quartics in u,v. The **full**
degree 24 discriminant is retained. No square factor or even-multiplicity
singularity is discarded.

At 521 and 523 the retained trace frames have monic degree-four h, coprime
h,Nx, the prescribed norm-twelve trace and the same generic height Gram.
The moving-pencil integrality proof consequently applies to the entire
parameter line at these primes. Each rational `[u:v]` has primitive
`Z_p` coordinates. If F vanishes there, its reduction vanishes on
`P1(F_p)`. Therefore one root-free projective reduction excludes **all**
rational singular parameters for that pencil, with no denominator or
rational-height bound.

| Prime | Pencils tested | Newly excluded | Remaining |
|---:|---:|---:|---:|
|521|49|25|24|
|523|24|24|0|

Every exclusion includes `F(1,0)`, the coefficient of `u^24`. Thus a
parameter at infinity and degree loss of the affine discriminant cannot
escape the check. No comparison between distinct reduced parameters is
used, so coalescence of two rational parameters creates no gap. No
simultaneous-good-reduction assumption across different pencils is needed:
each pencil has its own single-prime nonexistence witness.

The norm-twelve divisor classes are primitive nef isotropic, and their
pencils have no reducible fibres, by the retained minimum-norm proof.
Consequently the preceding genus-zero exclusion says that their rational
members are smooth genus-one curves. All geometric singular fibres lie
over nonrational points of the pencil parameter line.

## 3. Completing the arithmetic-genus-one rational-normalization row

An integral arithmetic-genus-one bisection B has `B^2=0`. In
`NS=U+M(-1)` it has class `(n,2,w)` with `w^2=4n`.
Translate by an inherited rational section to minimize w in its parity
class. The complete parity census gives nonzero isotropic minima
`4:1313, 8:63917, 12:49`; the zero class has minimum zero.
The norm-zero and norm-four cases have `B.O=n-2<0`, impossible for
an irreducible bisection distinct from O. Thus B lies in one of the
63,917 norm-eight or 49 norm-twelve pencils after translation.

This argument uses irreducibility and square zero, **not smoothness**.
Translation is an automorphism of the minimal elliptic K3 and preserves
arithmetic and geometric genus. A rationally defined member corresponds
to a rational point of its complete pencil parameter line.

The retained [norm-eight singular theorem](R17_NORM12_RATIONAL_NORMALIZATION_BOUNDARY_2026-09-04.md#exact-singular-genus-one-certificate)
excludes every nonsplit rational singular member in its 63,917 classes.
It explicitly accounts for the known split members and the entire
even-multiplicity discriminant locus. Section2 supplies the previously
missing full 49 moving norm-twelve layer. Hence **there is no geometrically
integral Q-defined bisection with `(p_a,g_norm)=(1,0)` on Q80**.

The norm-eight theorem is reused with its existing assurance level.
Its merger and retained shards passed the recorded integrity review, but
their polynomial coefficients are not retained for independent replay.
The new49-pencil replay does not repair or upgrade that inherited gap.

## 4. Descent excludes `(k,j)=(2,1)`

Suppose the system in section1 has such a solution. Put
`C:w^2=q`, `L=Q(C)`, `U_+=P_1+P_2`, `U_-=P_1-P_2`.
Both U have height16. At both branch points the P_i choose the same
nonzero 2-torsion point because q divides h and the parent fibre is smooth.

The product Kummer class

```
beta=(x_1-theta)*(x_2-theta),     theta^3+A*theta+B=0
```

has square norm and even valuation at every spectral place. At the branch
points the identical local parity triples cancel; away from them use the
same good/I1 local argument as in the
[Q80 unramified descent proof](COMMON_QUARTIC_BRANCH_STRATA_2026-09-14.md#3-the-unramified-arithmetic-input-on-q80).
That argument depends on the actual branch divisor, not its having four
points, and includes infinity and all singular parent fibres. Its retained
determinant948 arithmetic gives `beta=delta(T)` for an inherited rational
section T.

Kummer exactness over L now gives points

```
R_+=(U_++T)/2,         R_-=(U_-+T)/2.
```

This is actual divisibility in E(L). The torsion-free height formula makes
their traces exactly T and their anti-traces U_+,U_-. In particular neither
point is inherited. Let B be either bisection image. If sigma is the deck
involution and U is its anti-trace, intersection theory on the pullback gives

```
2*B^2=-8+2*(R.sigma(R)),
height(U)=8+2*(R.sigma(R)),
B^2=(height(U)-16)/2=0.
```

Here `R+sigma(R)` in the pullback divisor calculation means addition of
divisors. Since U is nonzero, B is a genuine geometrically integral
bisection with function field L. Its arithmetic genus is one, while its
normalization C has genus zero. This contradicts section3. A single one
of these images is already impossible; no singular-cover collision test
or rational point on C is needed.

## 5. Replay and the next unresolved boundary

The [frozen packet](../artifacts/generated-results/elkies-k3-q80-genus-zero-singular-gate-v1/input.json)
hash-binds the previous frames, generic model, parity and singular theorem.
The [constructor](scripts/certify_q80_genus_zero_singular_gate.sage) records
73 degree 24 discriminant polynomials and their complete projective root
lists, assigns a witness prime to all 49 pencils, and checks the universal
twelve-stratum identities and degrees. Each prime has a separate checkpoint.
No trace, lattice, cover-pair or specialization search is repeated.

The [independent Python checker](scripts/verify_q80_genus_zero_singular_gate.py)
rechecks all 73 used trace frames against the generic basis by the retained
height-bounded finite-fibre method. It verifies each discriminant using the
expanded sixteen-term formula at 25 distinct sites, which is exact by the
degree 24 bound, checks infinity separately and evaluates all 38,154
projective parameters. It uses neither Sage nor the producer's I,J formula.
The [receipt](../artifacts/generated-results/elkies-k3-q80-genus-zero-singular-gate-v1/independent-replay.json)
records a successful 0.208-second replay; the producer used 0.013 CPU seconds
after startup. The limits were 30 CPU seconds and 4 GiB. The
[execution record](../artifacts/generated-results/elkies-k3-q80-genus-zero-singular-gate-v1/execution.json)
retains the preliminary coercion failure and the pre-computation correction
to a helper-source binding. Written geometry and descent are not formal
proof-assistant verification.

```
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 research/elkies-k3/scripts/verify_q80_genus_zero_singular_gate.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m unittest discover -s research/tests -p 'test_q80_genus_zero_singular_gate.py' -q
```

At this singular-pencil gate the genus-zero chart had eleven remaining
`(k,j)` cases. For `k=2`, the
sign-equivalent cases `j=0,2` have cross height `+/-2`; the same descent
gives anti-trace heights12 and20, hence image arithmetic genera zero and
two. They are outside the arithmetic-genus-one exclusion. The
[subsequent branch-field theorem](Q80_RATIONAL_BISECTION_BRANCH_TORSION_2026-09-14.md)
closes them: every smooth rational bisection has a branch field without
nonzero2-torsion on its parent fibre, contradicting the original polynomial
sections. For `k<2`, the
branch factor outside g must split the parent cubic, giving a rational
effective divisor of degree `2-k` on its genus31 splitting curve. No
absence theorem for those degree-one or degree-two divisors is asserted.

These are new arithmetic gates, not instructions to reopen the completed
49 or 63,917 pencil banks. The repeated-common-quartic boundary is narrowed
by this result only where its residual q has `k=2,j=1`; its other cases
and a fixed-parent whole-family genus bound remain unknown.
