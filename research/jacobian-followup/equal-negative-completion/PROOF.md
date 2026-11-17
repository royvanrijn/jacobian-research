# Equal-negative weighted strict-profile exclusion through degree seven

Subsequent [KCUB40](../cubic-reduction/EQUAL_NEGATIVE_FULL_DEGREE7.md) closes
the full weight111 degree-seven class by a changed structural proof. This
strict-stage proof and its outside-range unresolved descriptions below are
historical;all original inputs,logs,forty-second timeouts and unreplayed
primary evidence remain unchanged and unused by the new closure.

KCUB23 excludes every three-variable weight(1,-1,-1) Keller counterexample
of ordinary degree<=7 whose affine-invariant degree flag is strictly below
(7,6,4). It also excludes every nonlinear positive first component of degree
at most5 in this class, without a restriction on the flag. It does not claim
that every degree7 map of these weights is an automorphism: a degree7
counterexample with two independent degree7 outputs remains outside the
strict profile range and is not completely resolved here.

Five independently multiplied rational unit certificates support the finite
normal charts:5868 multiplier monomials in total, maximum degree8. The
normalization, leading-form, fibre and residue arguments below are written,
not formalized. The global GLDC5/7/8 and Moh trust boundaries remain explicit.
There is no smaller positive witness or change to the20 remaining entire
unrestricted profiles. Work and source review are dated15 November2026.

## Whole ansatz and exact invariant dictionary

Normalize the tangent linear part of the weighted map by constant output
changes and recenter at the origin. It has the complete form

```
Y=xy, Z=xz,
F=(x*f(Y,Z), x^-1*A(Y,Z), x^-1*B(Y,Z)),
f(0)=1, A(0)=B(0)=0, A1=Y, B1=Z,
deg f<=3, deg A,deg B<=4.
```

All displayed expressions in x,y,z are polynomial: every monomial ofA,B
has invariant degree>=1. Their ordinary degrees are2*deg f+1 and
2*deg A-1,2*deg B-1. No bounded coefficient values are assumed.
Put {u,v}=u_Y*v_Z-u_Z*v_Y. The exact Jacobian condition is

```
Jac(F)=f*{A,B}+A*{f,B}+B*{A,f}=kappa,
```

or equivalently Jac(fA,fB)=kappa*f. The independent Fraction checker builds
all33 free coefficients of this complete ansatz and differentiates the full
three-variable map; it verifies the dictionary coefficient by coefficient.
The negative2-by2 linear block can be normalized to the identity because the
Jacobian at the origin is invertible. Later target changes may replace that
normalization by another invertible constant block, as stated explicitly.

If f is nonconstant, choose a zero(Y0,Z0) of f and take x=1. Its image has
first coordinate0. The plane x=0 maps linearly and bijectively to the whole
output plane F1=0, so it supplies a distinct point with the same image.
Thus every such weighted completion is automatically noninjective. If f=1,
GLDC5/KPAIR1 excludes a collision in this degree range by plane fibres.
For a Keller completion, f=0 must be smooth: an affine singular zero would
make every partial derivative of x*f(xy,xz) vanish at x=1.

## Two elementary leading-form facts

Homogeneous binary forms U_d,V_e with {U_d,V_e}=0 are powers of one common
homogeneous binary form: U=alpha*h^(d/g), V=beta*h^(e/g), g=gcd(d,e),
with h of degreeg. This follows by dehomogenizing and integrating the
zero logarithmic derivative, then unique factorization. In particular,
equal-degree forms are proportional; coprime degrees force h linear.

If A,B have top degrees d,e and f has top degreem, the highest condition is
{f_m*A_d,f_m*B_e}=0. When d=e, this forces A_d,B_d proportional.
If f is linear with topY, the equivalent leading relation is
(d+1)*A_d*(B_e)_Z-(e+1)*B_e*(A_d)_Z=0.
For d=4,e<=3, coprimeness of5 ande+1 forces A4 andB_e pure powers ofY.
Indeed atY=1 their common power would give A4 a fifth power of a
nonconstant polynomial of degree<=4, impossible.

If the actual negative top invariant degree is at most3, the source map
has degree<=5 unless f has degree3. When deg f<=2 and both negative
outputs reach degree5, their top binary cubics are proportional, producing
a target combination of degree<=3. GLDC8 excludes that degree5 map.
If either negative top degree is lower, the same small-component exclusion
applies directly. Hence with deg f<=2 a candidate must have one degree4
negative invariant polynomial. A constant target change makes A4 nonzero
and B4=0. It is essential not to demand tangent identity at the same time
without an additional justified normalization.

## Fixed linear f: complete two-chart reduction

For nonconstant linear f a linear change of negative source coordinates
puts f=1+Y. With A4 nonzero and B4=0, the preceding prime-degree argument
gives A4=a7*Y^4 and B3=b3*Y^3, where a7!=0 and b3 may vanish.
Lower terms of A through degree3 and B through degree2 are all retained.

For any invertible negative linear block, two normalizations cover it.
If B1 has a nonzero Z coefficient, a source Z-shear followed by constant
target shears/scales gives A1=Y,B1=Z. These source changes preserve f=1+Y.
If B1 has zero Z coefficient, it is a nonzero multiple ofY; source and
target changes give A1=Z,B1=Y. The constant determinants become1 and-1,
respectively. The two13-variable models include every remaining lower
coefficient and the saturation rho*a7-1. Independent multiplication gives
unit identities for all18 generators in each chart.

- [ChartI](normal-I/input.json),742 multiplier monomials,max degree7.
- [ChartII](normal-II/input.json),1403 multiplier monomials,max degree6.

The earlier full24-variable34-equation Singular run timed out at40 seconds.
It remains a timeout, not UNSAT. The finite normal-chart certificates and
the complete written reduction, rather than that timeout, prove the result.

## Quadratic f: every first component of degree five

Assume f2 nonzero. After removing B4, consider e=deg B<=3.
If f2 is squarefree and e=3, the top forms f2*A4 and f2*B3 have coprime
degrees6,5 and would be powers of a single line, impossible for squarefree
f2. For e=2, their degrees6,4 give powers h^3,h^2 of a quadratic h.
Squarefreeness forces h proportional to f2, hence B2 proportional to f2.
The cubic source component x^-1*B has leading form x*f2(y,z), three
nonconcurrent distinct lines. The degree-free nodal residue obstruction
reviewed in GLDC7 and explicitly used in GLDC8 excludes any Keller
completion. Its proof uses only a cubic component, exactness of its fibre
form and the boundary periods, not the degree of the other two outputs.
For e=1 the negative component is linear, excluded by KPAIR1.

It remains f2 a square. Set f2=Y^2 by an invertible negative source change.
For e=3, the coprime6,5 argument gives A4 proportional to Y4 andB3 toY3.
For e=2 it gives A4=alpha*Y*L^3,B2=beta*L^2 for a linear L.

On any smooth zero curve C of f, the Jacobian identity restricts to
A*delta_f(B)-B*delta_f(A)=kappa, with delta_f={f,-}. On a graph or a vertical
line this is a constant nonzero polynomial Wronskian in its affine parameter.
Two polynomials with such a Wronskian both have degree at most1: subtract
proportional highest terms when their degrees agree, then compare the
noncancelling highest term for distinct degrees.

If the transverse linear term of f is nonzero, a linear source change gives
f=1+Z+Y^2. On C, Z=-1-Y^2. In the e=2 case, L involving Z makes B2|C
have degree4, while B1|C has degree at most2; this contradicts the Wronskian
bound. Thus L is proportional toY. In every remaining case A4 is pureY4,
B3 is pureY3, and all lower coefficients are arbitrary.
Only target transformations are now used, because the source normalization
must preserve the fixed parabola. If the Y coefficient of B1 is zero,
scale B1 toZ, subtract its multiple fromA and rescale to A1=Y. Otherwise
scale B1 toY+tZ, subtract its multiple fromA to remove A's Y coefficient,
and rescale to A1=Z. These two charts cover every invertible linear block.
Their determinants are1 and-1, and no assumption is made about t.
The two complete models have13/14 variables and23 generators, including
rho*a7-1. Their exact unit identities have179 and3512 multiplier monomials,
maximum degrees3 and5. No parameter t is assumed zero or invertible.

If the transverse linear term vanishes, f=1+aY+Y^2. Smoothness forces two
distinct nonzero roots. For e=2 the same Wronskian on either vertical root
line forces L proportional toY; e=3 is already pureY. The restrictions
of A,B to both root lines are affine inZ. Since A4 andB3 have noZ terms,
the Z^3 coefficient ofA is constant and its Z^2 coefficient has degree<=1
inY; both vanish at the two roots and are zero. The Z^2 coefficient ofB
is constant and likewise zero. Thus A,B are affine inZ. In source coordinates
F is affine in z, with first component independent ofz. The affine reduction
lemma reviewed in GLDC7 yields a plane Keller pair of degree<=13, within
Moh100, and determinesz from the nonvanishing derivative pair. Hence F is
injective, contradicting the automatic collision. This covers the parallel
line case and the entire quadratic-f sector.

## Cubic f: the whole strict-profile branch

If deg f=3 and either negative invariant output has degree4, the positive
and negative source leading forms have different powers x4 andx3. They
are independent, giving two degree7 entries in the invariant flag.
Such a flag is lexicographically above764 and outside this theorem.

In the strict range both negative invariant degrees are at most3. Their
highest equal-degree forms are proportional, so a constant target change
gives negative source degrees at most5 and3. If both were at most3, the
same leading argument would give a linear component, excluded by KPAIR1.
Any surviving cubic negative component has the form
Q=ell(y,z)+x*q2(y,z), with ell nonzero. A squarefree q2 is excluded by the
nodal cubic residue argument. If q2 is a square of a line independent of
ell, Q is a graph coordinate with completion inverse degree3; plane fibres
have degree<=21 and Moh applies. The linear case is already excluded.
Consequently only Q=z+c*x*z^2 remains, after negative linear changes.
Rescale to c=1 and remove its linear part from the other negative output.

Permute the source to(X,Y,Z)=(z,x,y), permute outputs to(Q,F1,P), and reverse
all weights. The weights are now(1,-1,1) and the fixed first component is
r=X+X^2Y. With U=XY,V=YZ, the complete remaining pair is

```
P=Y*f(U,V), Q=Y^-1*B(U,V),
f0=1, B1=V, deg f<=3, deg B<=3.
g=U+U^2,
Jac(r,P,Q)=f*{g,B}+g*{f,B}+B*{g,f}=1.
```

The full opposite-orientation three-variable Jacobian dictionary is checked
independently. The actual top degrees are3,3; if B had lower degree the
original negative component would already have been excluded. The leading
forms g2*f3 and f3*B3 have coprime degrees5,6. They are powers of one line;
the factor U^2 forces that line to beU. Therefore f3=a5*U^3 andB3=b3*U^3.
All lower coefficients are retained. The resulting11-variable16-generator
ideal, with rho*a5-1, has an independently multiplied unit identity with
32 multiplier monomials,max degree8. This closes the cubic-f strict branch.

## Retained outside-range experiments and positive control

The fixed cubic quotient f=1+YZ+Y^3 with one degree7 negative output has
profile77* and is outside the strict range. Both of its small charts have
independently checked unit identities, with6 and5 multiplier monomials.
The wider family f=1+alphaY+betaZ+YZ+Y^3 has a primary unit result in its
first chart and a40-second timeout in the second. The large first-chart
certificate has not been replayed and is not used in this theorem. No
whole cubic-quotient family, all weight111 degree7 maps, or degree6/7
unrestricted profile is declared excluded from those experiments.

A separate exact positive control applies b->b-a*c/3 to the known13-term
three-variable map. It gives degree profile(9,7,3),19 nonlinear occurrences
and an exact collision, with determinant1 independently expanded. It retains weights(1,-1,-2), so it is not an equal-negative upper witness.
This is not a global improvement. Together with GLDC8 it gives only a sufficient
interval[6,9] for the least degree in the extra class with a cubic output
combination; no endpoint is proved attainable optimally.

## Replay and unchanged goal

```
.venv/bin/python research/jacobian-followup/equal-negative-completion/verify.py
```

[Validation](validation.json) records all five required certificates,
the33-parameter dictionary, the opposite dictionary, the two unused finite
controls and both preserved timeouts. The proof is written, with independent
exact certificate arithmetic and imported scalar/plane facts, not formal
or independently externally refereed. Raw13, normalized[10,48], cubic12/36,
cubic-linear89 and global3D degree[5,7] are unchanged. Twenty whole dense
profiles remain unresolved. The full simplicity goal remains open.

## Reused mathematical premises

- [GLDC5/6 small-component/classical source notes](../small-component-gate/CLASSICAL_CUBIC_IMPORT.md).
- [GLDC7 written residue/affine reduction audit](../quartic-three-audit/INVERTIBILITY_AUDIT.md).
- [GLDC8 degree5 component exclusion and explicit degree-free scalar observations](../quartic-three-audit/EXCEPTIONAL_DEGREE5.md).
- [KPAIR1's plane-degree import and exact scope](../sparse-eight/COORDINATE_GATES.md).

The nodal cubic residue statement is used in the degree-free form justified
by its exact fibre2-form argument. Graph completions have inverse degree<=3;
the affine-in-z plane pair has degree<=13. Both lie within the already
imported Moh100 bound for this source degree<=7 problem.
