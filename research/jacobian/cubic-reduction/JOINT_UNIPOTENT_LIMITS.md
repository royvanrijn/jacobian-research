# Complete five-parameter joint-frame limits

KCUB31,15 November2026. Fix the integral twelve-variable cubic K of
[KCUB3](PROOF12.md), and consider S with diagonal1 and the only possible
off-diagonal entries

```
S_(2,7)=u, S_(2,10)=v, S_(3,9)=w,
S_(3,8)=z+w*t, S_(9,8)=t.
```

This parametrizes every matrix with the five off-diagonal positions
proposed after [KCUB30](OUTPUT_FRAME_TRADEOFF.md): the original free
(3,8) entry is z+w*t, and z is unrestricted. There is no omitted parameter
stratum. All parameters lie in an arbitrary characteristic-zero field.

The complete results are:

| Class | Minimum nonlinear occurrences | Additional sharp statement |
| --- | ---: | --- |
|S*K*S^-1, identity linear part|36|Only S=I attains36; every other S needs at least39, and39 occurs|
|Distinct nonlinear source monomials in K*S^-1|30|Only S=I attains30; every other S needs at least33, and33 occurs|
|A*K*S^-1 for arbitrary invertible output A|31|Imported KCUB30 attains31 at S=I; every nonidentity S needs at least34, and34 occurs|

Thus this entire source family, even with arbitrary invertible output
mixing, cannot improve the31-term free-linear-part benchmark or the36-term
identity-linear-part benchmark. No counts37 or38 occur in its prescribed
identity-linear-part conjugacy class. It gives no global affine-orbit,
sparsity, degree or dimension optimum and no smaller headline witness.
It does not settle total-term optimization over all A and S. The global
simplicity objective remains open.

## Matrices, maps and inherited collision

The source inverse has the particularly simple form

```
x2 -> x2-u*x7-v*x10
x3 -> x3-w*x9-z*x8
x9 -> x9-t*x8,
```

with all other source coordinates fixed. The independent checker verifies
both polynomial matrix inverse identities, (S-I)^3=0 and an explicit
factorization into five determinant-one elementary matrices. Hence
G=S*K*S^-1 is cubic, has derivative I at0, and determinant1. Its points
S*p and S*q collide, where p,q are the saved distinct integer collision
points. Their x1 coordinates remain0 and-1 for all parameters, so they
never merge. The complete parameter-polynomial collision identities are
checked directly. The universal determinant inheritance is the written
chain-rule consequence of the verified matrices and pinned KCUB3 premise.

There are90 generic nonlinear coefficient slots in G and19 polynomial
factors. All90 scalar-times-factor-product identities are independently
checked. A count refers to expanded monomial occurrences after collecting,
not the number of terms in a factored expression. Let [f] be1 if f is
nonzero and0 otherwise.

## The normalized minimum36 and next value39

Discarding all but mandatory disjoint coefficient groups gives

```
N(G) >= 24+10[u]+10[v]+5[w]+5[z].
```

If u and v are both nonzero, this is at least44. If at least one is nonzero
and either w or z is nonzero, it is at least39. The remaining cases are
covered completely by the following exact indicator identities.

When u=v=0, the count is

```
N = 31+6[w]+5[z]+[z-4]+[(t-3)*w+z+8]+3[t]
       +[w-2]+[z+2*t]+[w+4]+[t-3].
```

If w and z are both nonzero, the first terms already give42. If w is
nonzero and z=0, use [z-4]=1 and [w-2]+[w+4]>=1 to get at least39.
If w=0,z!=0, the two w-offset indicators contribute2, and
3[t]+[t-3]>=1 gives at least39. Finally, w=z=0 reduces exactly to

```
N=35+4[t]+[t-3].
```

This equals36 at t=0,39 at t=3, and40 at every other nonzero t.

On v=w=z=0 with u!=0, the exact formula is

```
N = 42+4[u+1/3]+[u+3]+[u+1/6]
       +3[t]+[t*(u+1/3)]+3[t-3],
```

so N>=42. On u=w=z=0 with v!=0 it is

```
N = 44+3[v-1]+[v+12]+[v-1/2]+5[t]+3[t-3],
```

so N>=44. This completes every parameter case. The sole36 point is
u=v=w=z=t=0, exactly S=I, and every other point has N>=39.

The independent checker reconstructs each specialization identity up to
nonzero rational factors and verifies the indicator polynomials on all512
Boolean assignments in the first branch and32 in each other branch.
These are exhaustive truth-table identities, not sampled field points.
The first branch also uses necessary compatibility facts such as
(w-2)-(w+4)=-6 and t-(t-3)=3, all checked as exact polynomial identities.
The written field case split and nonzero-constant arguments establish
coverage over every characteristic-zero extension, including complex
parameters. No real-only or generic-parameter assertion is used.

## Counting source monomials before any output change

Put R=K*S^-1. There are68 generic nonlinear source monomials. For each
source monomial e, write its coefficient vector as

```
c_e(P)=g_e(P)*q_e(P),
```

with a supplied rational constant combination of the entries of q_e equal
to1. The checker verifies every multiplication and all68 unit identities.
Consequently the vector vanishes at a parameter exactly when g_e does.
No primary gcd computation, factor irreducibility claim or denominator
stratum is trusted. The scalar factors of g_e are also checked exactly.

The number C of nonzero nonlinear source coefficient vectors satisfies

```
C >= 24+9[u]+9[v]+4[w]+4[z].
```

Thus u!=0 or v!=0 already forces C>=33. When u=v=0, the complete identity is

```
C = 29+4[w]+4[z]+[w-2]+2[t]+[z+2*t].
```

This is30 at w=z=t=0. If w!=0 it is at least33. If w=0,z!=0 it is at least34.
If w=z=0,t!=0 it equals33. Thus every nonidentity S has at least33 nonzero
nonlinear source columns. Invertible output changes cannot erase any such
column. The t=3 controls have exactly33, proving sharpness of that bound.

## Every nontrivial source change needs34 terms even with free outputs

The column count can be strengthened by one where needed. If v!=0 and
u!=0, C>=42. If v!=0 and u=0, the mandatory24+9 is supplemented by three
source columns whose factors u+1/3,u+1/6 and u*w-2*u-2/3 are nonzero
constants on u=0. Hence C>=36 whenever v!=0.

When v=0, the three columns for x1*x2, x1^2*x3 and x1*x10 are always

```
3*e4+e12, e4, e12,
```

independently of all other parameters. They form a pairwise-nonparallel
triple in a plane. After any invertible output change they need at least
four occurrences among their three source monomials. This is the same
elementary projective support argument used by KCUB30, and the three
actual parameter-polynomial vectors at v=0 are reconstructed directly.
It forces an occurrence beyond C. Thus for nonidentity S with v=0,
N(A*R)>=C+1>=34; for v!=0 the stronger36 lower bound already applies.

At S=I, imported KCUB30 gives the sharp31 minimum over all invertible
outputs. Therefore31 remains sharp for the whole five-parameter source
family with arbitrary outputs, and a genuinely changed source needs at
least34. Adding output constants does not change nonlinear counts.

The KCUB26 rigidity theorem gives an additional corollary: every polynomial
target map of degree<=5 which makes K*S^-1 a cubic Keller map is affine.
Indeed composing with the invertible linear source S preserves ordinary
degree, so each target component falls under KCUB26. Its affine linear
part must be invertible because the result is Keller. All count bounds
above therefore also cover that degree<=5 target class. The previously
checked6175-column theorem is an explicit imported premise, not rerun.
Higher target degrees involving all outputs remain outside this corollary.

## Sharp matching controls and normalization

At (u,v,w,z,t)=(0,0,0,0,3), the normalized G has39 nonlinear and12 linear
occurrences. Its source support is33. For the same source, choose the
integer output matrix A with diagonal1 and off-diagonal entries

```
A_(2,7)=-3, A_(2,10)=-12,
A_(3,9)=-4, A_(3,8)=-20,
A_(9,8)=3, A_(11,10)=3.
```

Then A*K*S^-1 has34 nonlinear and17 linear occurrences. Both controls
total51. Every coefficient, both integral collisions and both full
symbolic Jacobian determinants are independently checked. Their Schur
computations keep all12 variables symbolic and verify the9-by9 carrier
inverse in both orders. These controls certify the new class lower bounds;
they are dominated by the saved36/12 and31/17 headline examples and are
not advertised as improved counterexamples.

The zero parameter point is coefficientwise K. The old output-clearing
point (-3,-12,-4,-8,3) is coefficientwise KCUB30's already checked74-term
identity-linear-part conjugate. Both comparisons are model controls;
neither is a new search result or an excuse to repeat an unchanged point.

## Replay, assurance and remaining routes

Run

```
.venv/bin/python research/jacobian/cubic-reduction/joint_unipotent5.py
.venv/bin/python research/jacobian/cubic-reduction/verify_joint_unipotent5.py
```

The primary [producer](joint_unipotent5.py) uses SymPy only for discovery
of coefficients, factors and source-vector quotient units. The independent
[checker](verify_joint_unipotent5.py) reconstructs source substitution and
output multiplication with standard-library rational polynomials. It
checks all matrices, complete maps and generic collisions;90 coefficient
factor products;68 source-column factorizations and residual-vector units;
all branch indicator identities; the forced projective triple; and both
matching maps' full determinants, points and counts.

- [Full parameter model and standalone control coefficients](joint-unipotent5/model.json)
- [Independent receipt](joint-unipotent5/validation.json)

These exact finite checks support the written field case split, support
invariance and chain-rule arguments. This is not formal verification or
external mathematical review. No bounded point search, cap increase,
heuristic score, generic factor assumption or timeout is used as a proof.

For fewer than36 identity-linear-part or31 free-linear-part cubic terms,
leave this source family or change the construction. Larger linear source
blocks, affine source centres, nonlinear graph/source changes, coordinated
target degrees>=6, and other counterexample constructions remain open.
The theorem does not settle total-term optimization over all A,S or the
full affine source/output orbit. The ordinary-cubic dimension interval
[4,12], raw3D13, normalized[10,48], cubic-linear89, global3D degree[5,7]
and strict62/12 profile partition remain unchanged.
