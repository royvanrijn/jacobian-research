# Collision-free classification of the degree-seven negative-weight class

**Theorem (KCUB14).** Let F be a three-variable polynomial map over a
characteristic-zero field, with identity linear part, determinant one,
ordinary degree at most seven and source/output weights (1,-1,-2).
Assume the first-output coefficients on a^2b and a^3c are both nonzero.
Then F is diagonally conjugate to the established 13-term representative.
In particular it has exactly 13 nonlinear monomial occurrences and
has a collision. **No collision, sparsity bound or numerical coefficient
range is assumed.** All 23 compatible nonlinear atoms are allowed.

This proves a matching 13-term minimum throughout this weighted sector
with those two coefficients nonzero. Zero values of either coefficient,
other weights, arbitrary polynomial coordinates, and the global degree,
sparsity and special-form dimension minima remain outside the theorem.

## Normalization and the complete model

If those two coefficients are mu,nu, put
D=diag(1,-3/mu,-1/nu). The map D^{-1}F D still has identity linear
part, determinant one, the same degree and exactly the same nonlinear
support. Its two coefficients become -3,-1. It therefore suffices to
classify this normalized case.

Use u=ab,w=a^2c and F=(af,a^{-1}g,a^{-2}h). As in
[KCUB13](WEIGHTED_FULL_DEGREE7_CLASSIFICATION.md),

\[
J=g(f_u h_w-f_w h_u)+2h(g_u f_w-g_w f_u)
       +f(g_u h_w-g_w h_u)=\det JF. \tag{1}
\]

Here, crucially, the four first-output coefficients are independent:

\[
f=1-3u-w+R u^3+P u^2+Q u w+T w^2. \tag{2}
\]

There is no relation R=P-4Q+16T and no collision equation.
The model has four first-output parameters and all eight g-tail and
nine h-tail parameters, for 21 coefficient parameters in 23 nonlinear
source atoms. Saved parameter names R,P,Q,T are f1,f2,f3,f4.

The smooth-cubic boundary lemma in KCUB13 applies to (2) without using
its former collision relation: a singular point of f=0 contradicts
(1), while RT!=0 and smoothness give a projective cubic with one point
at infinity. The identity d(h/g^2)=-omega/g^3 then has incompatible
zero and pole orders. Thus **RT=0**. The lemma is written mathematical
proof; its projective-curve argument is not formally machine verified.

## The parameter branches are exhaustive

The following three lifted identities have explicit rational polynomial
multipliers. Rational nonzero scalar factors in their targets do not
change the stated vanishing consequences.

| Branch | Cleared target in the necessary ideal | Consequence |
|---|---|---|
| T=0 | Q^16 | Q=0 |
| R=0 | T^2(PT-Q^2/4)^7 | T=0 or 4PT=Q^2 |
| T=Q=0, polynomial boundary constraints | P R^5 | P=0 or R=0 |

The polynomial boundary constraint used in the last row is elementary.
Whenever u=t,w=1-3t+Pt^2+Rt^3 parameterizes f=0, its extension by
v=f has Jacobian partial(u,w)/partial(t,v)=-1. Restricted tails
G(t),H(t) satisfy G H'-2G'H=1. The polynomial Wronskian lemma from
[KCUB11](WEIGHTED_TAIL_CLASSIFICATION.md) forces deg G<=1,deg H<=2.
These are necessary linear coefficient equations.

Four further reduced full-ring branches are needed:

| Branch | Treatment |
|---|---|
| T=Q=P=0, R!=0 | Exact unit identity after all boundary linear equations |
| T=Q=R=0, P!=0 | Exact unit identity after all boundary linear equations |
| R=0,T!=0,4PT=Q^2, kappa!=0 | Exact unit identity below |
| R=0,T!=0,4PT=Q^2, kappa=0 | Elementary parallel-line proof below |

In the last two rows set kappa=Q/(2T)-3, so
P=(kappa+3)^2 T and Q=2(kappa+3)T. For kappa!=0 the polynomial-line
boundary chart is

\[
u=(-1+t-Tt^2)/\kappa,\qquad w=t-(\kappa+3)u.
\]

It extends by replacing -1 with v-1, where v=f, and has Jacobian
-1/kappa. Thus its boundary Wronskian is 1/kappa; the same degree
bounds apply. Clearing coordinate denominators gives 13 necessary
linear equations. Their complete row reduction leaves five tail
variables. Its only additional linear-solution denominator is T,
already nonzero. The exact full-ring certificate includes
1-branch_inv*kappa*T^2. This excludes every kappa*T!=0 solution,
including all roots of the extraneous factor left by the earlier
generic field certificate. No exceptional hypersurface is omitted.

The pure cubic and quadratic charts have linear-solution denominator
one. Their unit identities include 1-branch_inv*R and
1-branch_inv*P respectively. Together with the three cleared identities,
these close T=0 and the nonparallel part of R=0.

## The parallel-line branch

For kappa=0 write L=3u+w and phi(L)=1-L+T L^2, with T!=0.
Then f=phi(L). If T=1/4, f=f_u=f_w=0 along L=2, contradicting (1).
Otherwise phi has two distinct roots, and f=0 is two disjoint affine
lines. On either line its nonzero constant normal derivative makes
the boundary Wronskian a nonzero constant. The degree lemma therefore
applies on each line.

Expand g,h in u,L. For every u-power above the allowed boundary
bounds, its coefficient has degree less than two in L, as direct
source-atom enumeration shows. It vanishes at both roots of phi, so
it vanishes identically. The complete seven-equation linear reduction
therefore gives, globally,

\[
g=A(L)u+B(L),\qquad h=C(L)u^2+D(L)u+E(L),
\]

where deg A<=2, deg C<=1. Identity linear part gives

\[
B(0)=E(0)=0,\qquad D(0)=-3. \tag{3}
\]

In u,L coordinates (1) becomes

\[
J=-\phi'(g h_u-2h g_u)+\phi(g_u h_L-g_L h_u).
\]

If A=0 identically, then J=-(phi B)' h_u=1. This forces
(phi B)' to be a nonzero constant, impossible for nonzero polynomial
B and degree-two phi.

Otherwise the u^2 coefficient gives A C'-2A'C=0, hence C=c A^2
for a constant c. Since deg C<=1, either c=0 or A is constant.
Replace h by h-c g^2. This polynomial target shear preserves (1)
and (3), and makes C=0. It need not preserve the original degree
bound; no further degree bound on the modified D,E will be used.
The bound deg A<=2 is unchanged.

Now the u coefficient is

\[
A^2(\phi D/A)'=0.
\]

Thus phi D=k A for a constant k. Since D(0)=-3, k is nonzero.
The degree-two phi and deg A<=2 force D to be the constant -3
and A=alpha phi with alpha!=0. The constant coefficient in J is then

\[
J=\bigl(\alpha\phi^2 E-D\phi B\bigr)'=1.
\]

Integrating and using (3) at L=0 gives
phi*(alpha phi E-D B)=L. A degree-two polynomial cannot divide L.
This excludes the whole parallel-line branch, including every value
of T!=0, without trusting a timed-out unit lift.

We have now proved R=P=Q=T=0.

## The collision is a conclusion

It remains to classify the tails with f=1-3u-w, without assuming a
collision. On its polynomial boundary, write

\[
g(u,1-3u)=A u+B.
\]

Two exact unit identities exclude A=1 and B=0, each allowing all
17 tail atoms. Therefore beta=B/(A-1) is nonzero and

\[
g(-\beta,1+3\beta)=-\beta.
\]

Take q=(-1,beta,1+3beta). Its first output is zero and second output
is beta. At a=0, the map has first output zero, second output b and
third output c+I b^2 for some coefficient I. Choose c0 so that
p=(0,beta,c0) has the same third output as q. This supplies a collision
in precisely the nonzero-beta family classified by KCUB11. Hence
beta=1, c0=-8, and the tails equal the established 13-term map.
No collision assumption has been reintroduced.

Returning from normalization gives F=D F0 D^{-1}. Its explicit
collision is D(0,1,-8),D(-1,1,4), with common value D(0,1,0).
Diagonal scalings preserve all 13 nonzero nonlinear occurrences.

## Exact replay and retained limits

Run

```sh
.venv/bin/python research/jacobian/cubic-reduction/verify_weighted_unmarked7.py
```

The [independent validation](weighted-unmarked7/validation.json) uses
standard-library Fraction arithmetic. It reconstructs the universal
three-source-variable determinant, verifies that no collision equations
entered the model, and checks all boundary charts and necessary linear
equations. Complete row reductions use only explicitly permitted
nonzero parameter monomials as pivots, exposing every possible
exceptional divisor. It independently multiplies all eight used
certificates: three cleared parameter identities and five full-ring
units, with 30,337 multiplier monomials and maximum multiplier degree 15.
It also checks the polynomial identities and source degree bounds in
the parallel-line proof. The curve, degree and divisibility arguments
are written proof, not proof-assistant formalization.

The unmarked fixed-first tail solve, both before and after linear
reduction, retains its 20-second timeout. The two slope/intercept
strata and the automatic fixed-point collision replace that direct
solve. Generic quadraticgraph and parabola lifts produced unnecessary
exceptional factors; the reduced full-ring units close their entire
nonzero branches instead. The generic parallel lift finished Singular
export but its redundant symbolic postprocessing was stopped and
recorded. Parallel full-ring lifts, simpler saturation, and a combined
basis/transformation attempt all retain their 30-second timeouts;
the elementary proof above is the certificate-independent resolution.
The [Singular manual](https://www.singular.uni-kl.de/index.php/singular.pdf)
documents the combined basis/transformation operation used in that
retained attempt. No timeout or generic result is labeled a global
UNSAT certificate.

The subsequent [KCUB15 classification](WEIGHTED_COMPLETE_DEGREE7_CLASSIFICATION.md)
closes both coefficient-zero boundaries as well: every degree-at-most-seven
weighted counterexample has13 terms. Its linear-first-component branch is
an automorphism sector. Other weights and coordinate forms remain available. The separate
ordinary cubic12/36, cubic-linear89 and normalized three-variable
[10,48] bounds are unchanged; none is a global optimum certificate.
