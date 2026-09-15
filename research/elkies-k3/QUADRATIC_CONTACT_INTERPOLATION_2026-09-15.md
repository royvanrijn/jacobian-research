# Quadratic contact interpolation for a second singular genus-one layer

**Written theorem and exact interpolation identities; no matching MW17 covers.**
On a fixed24-I1 elliptic K3 over a number field, geometrically integral bisections
of arithmetic genus3 and normalization genus1 are finite modulo inherited
translations. A candidate is determined by a degree-two contact divisor on its
genus-nine halving curve. In suitable quotient charts its construction is a
four-by-four linear system with the norm or trace conditions below.

This extends the [arithmetic-genus-two tangent construction](SINGULAR_BISECTION_TRACE_GEOMETRY_2026-09-15.md)
to the next image genus. Neither finite set has been enumerated or proved empty.
The [two-gain objective](CORRELATED_QUADRATIC_GAINS_2026-09-12.md) remains open.

## 1. Contact rigidity, including higher tangencies

Let B have arithmetic genus a and normalization genus1. Fix its trace T and
quotient pi:S -> Y by eta(P)=T-P. Suppose T is geometrically primitive modulo2.
Then Y is F_n with n<=2 and smooth irreducible branch curve R=D_T of genus9,
as in the cited quotient argument. The image C of B is a section of Y over
the original t-line and

```
C^2=a-1,    R.C=2a+2.
```

At a contact P of multiplicity m_P between C and R, the local double cover has
delta=floor(m_P/2). Define the effective contact divisor on R

```
Z=sum_P floor(m_P/2)*P,     deg Z=a-1.
```

It is defined over the number field whenever B is. This includes ordinary
nodes, cusps and higher contacts, rather than assuming distinct double roots.

**Rigidity:** fixed T, a and Z determine at most one integral section C.
All sections of the ruling with the same self-intersection have the same
numerical class C0+kF. If two distinct ones C,C' had restriction to R vanishing
along2Z, their local graphs relative to smooth R would agree to order at least
2*mult_P(Z) at each P. Hence

```
C.C' >= 2*deg Z = 2(a-1) > a-1=C^2
```

for a>1, a contradiction. Distinct integral sections have no common component.
This proves uniqueness, not existence: the jet constraints may have no solution
or only reducible divisors.

## 2. Why the degree-two contact set is finite

For a=3 the trace must be geometrically primitive. Otherwise translate by a
geometric half of T. The translated integral bisection has trace zero, whereas
the parent intersection identity would give

```
0=4*(B.O)+10-2*3=4*(B.O)+4,
```

contradicting nonnegative intersection of distinct effective curves B and O.
This argument needs no rationality assumption on the proposed half.

The halving curve R has a primitive degree-four ruling and genus9. It has no
degree-two map to a genus0 or genus1 curve. If such a map and the ruling were
independent, Castelnuovo–Severi would give genus at most

```
4*0+2*1+(4-1)*(2-1)=5.
```

If they were not independent, primitivity of the degree-four ruling would
force the degree-two map to factor through it, which is impossible.
By [Kadets–Vogt, Theorem1.2(1)](https://arxiv.org/pdf/2208.01067), together with
finiteness of rational points in genus9, R has only finitely many closed points
of degree at most2 over the fixed number field. Thus there are finitely many
rational effective divisors Z of degree2, including doubles of rational points.
Contact rigidity gives finitely many B for each trace. Finite MW/2 and inherited
translation give the stated global finiteness and finitely many literal covers.
No effective height bound follows.

Here C^2=2 forces n even. The possibilities are C=C0+F on F0 or C=C0+2F on F2.
Each complete linear system has projective dimension3. The branch restriction
to C has degree8; removing its even contact divisor2Z leaves degree4. It must
be squarefree to produce the required genus-one normalization.

## 3. Explicit two-point jet equations

Use an affine chart where the distinct contact points have coordinates
(t_i,u_i), with t_1!=t_2. A section of the ruling cannot contain distinct points
above the same original t. At any eligible contact dt is nonzero on R, since
a section cannot be tangent to a vertical tangent. Set v_i=du/dt on R; for an
implicit equation R(t,u)=0 this is -R_t/R_u.

On F0 a section has equation

```
A+B*t+C*u+D*t*u=0.
```

Vanishing and tangency at each point give rows

```
[1,t_i,u_i,t_i*u_i],    [0,1,v_i,u_i+t_i*v_i].
```

The determinant of the four-row matrix is

```
(u_2-u_1)^2-v_1*v_2*(t_2-t_1)^2.                 (1)
```

Its vanishing is necessary. For rank3 the kernel supplies the candidate;
AD-BC!=0 is required for an integral (1,1) section. A zero determinant with
rank2 can give only reducible divisors and must not be accepted.

On F2 the corresponding equation is

```
A+B*t+C*t^2+D*u=0,
```

with rows [1,t_i,t_i^2,u_i] and [0,1,2*t_i,v_i]. Its determinant is

```
(t_2-t_1)*((t_2-t_1)*(v_1+v_2)-2*(u_2-u_1)).     (2)
```

Here D!=0 is required for an integral section. Coordinate charts and infinity
must be transported with the ruled surface; these affine formulas are not
an exclusion at omitted points. For Z=2P use four consecutive vanishing jets
at P instead of two pairs of jets. The general rigidity and finiteness include
this diagonal case.

## 4. Quadratic fields turn the determinants into norm and trace equations

Suppose the contact divisor is one quadratic point with

```
t^2-s*t+p=0,     u=a+b*t,     v=c+d*t
```

in its residue field. Equality of that field with the original t-field is
necessary for an eligible section. Conjugating t to s-t, equation(1) becomes

```
Norm(v)=c^2+c*d*s+d^2*p=b^2.                     (3)
```

Equation(2) becomes

```
Trace(v)=2*b,    equivalently c=b-s*d/2.         (4)
```

All four linear conditions can be expanded in the basis1,t over the ground
field, so a rank3 kernel gives the candidate over that same ground field.
These equations are evaluated at an actual point of the actual halving curve.
Arbitrary quadratic numbers satisfying(3) or(4) are not incidence evidence.

The checker contains two such linear-algebra controls over Q(sqrt2), recovering
u=(t+1)/(t+2) on F0 and u=t^2+t+3 on F2. They are not points on a K3 halving
curve. A separate same-ruling control has determinant zero and rank2, with
every solution reducible. These controls check the algebra, not the positive
research target.

## 5. Completing a candidate and obtaining two gains

For an actual eligible contact divisor, reconstruct C, restrict the double-cover
branch equation and remove the exact square contact factor. Retain its literal
scalar. Require a squarefree residual binary quartic, an integral inverse image,
and branching at smooth parent fibres for the following height test. The
resulting bisection has anti-trace height24 by the previous genus formula.

If that literal quartic is the same quadratic extension as a verified smooth
genus-one bisection, the two heights24 and16 have nonsquare ratio3/2, proving
two independent new directions. Matching an arithmetic-genus-two bisection
likewise gives ratio24/20=6/5. One must still prove that the covering genus-one
curve has a rational point and a rational nontorsion point. None of these matches
or arithmetic base conditions is supplied by equations(3) and(4) alone.

The [checker](scripts/verify_quadratic_contact_interpolation.py) and
[certificate](../artifacts/generated-results/elkies-k3-quadratic-contact-interpolation-v1/result.json)
verify the determinants, quadratic-field identities and interpolation controls:

```
.venv/bin/python research/elkies-k3/scripts/verify_quadratic_contact_interpolation.py
```

The geometric rigidity, monodromy and finiteness arguments remain written or
inherited mathematics, without formal verification or external review. No
arithmetic point census, new parent construction or curve search was run.
