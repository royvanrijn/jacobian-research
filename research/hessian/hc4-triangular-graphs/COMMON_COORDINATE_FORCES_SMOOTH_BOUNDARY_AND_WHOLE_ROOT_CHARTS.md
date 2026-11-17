# A common regular source coordinate forces smooth boundary and whole root charts

This is a conditional repair of the whole-component root-chart step.
For the actual canonical Keller completion, **given** one regular
polynomial source coordinate on a component open, the entire boundary
curve is algebraically A1 and admits the whole normalized plane chart
proved below. Existence of that coordinate, or one common coordinate for
all components, is not proved.

The argument uses the complete affine component opens and topological
boundary theorem, then removes the singularities that topology alone
could not exclude. It does not use a classification of affine surfaces
or assume a polynomial root chart in constructing the additive action.

## 1. Precise hypotheses and conclusion

Let M be a smooth connected affine complex surface with a dense open
```
X=Spec C[x,y]=A2,   M minus X=D,
```
where D is an irreducible reduced curve homeomorphic to C, with
normalization A1. Suppose M has a nowhere-zero regular two-form Omega
restricting to dx wedge dy on the WHOLE X. Suppose the given polynomial
source coordinate y is regular on M.

For an actual plane Keller map, the
[canonical completion theorem](CANONICAL_ETALE_COMPLETION_BOUNDARY_TOPOLOGY_AND_FREE_DECK_ACTION.md)
provides M=M_i=X union D_i, smooth affine and open in the FULL original
normalization, together with the stated boundary topology. The pulled-back
target form is a suitable Omega. If y is regular on the canonical U,
it is regular on every M_i. Different individual coordinates y_i regular
on different M_i also suffice for the individual conclusion here; this
does not align them globally.

The conclusion is that there are lambda in C and an integer r>=1 with
```
div_M(y-lambda)=L_lambda+rD,   L_lambda={y=lambda} in X,   (1)
```
where L_lambda is a closed source A1 and is disjoint from D.
Moreover D is a SMOOTH algebraic affine line.

Normalize the finite cyclic cover of M given by
```
z^r=y-lambda.
```
Its source inverse image is the ENTIRE A2_(x,z). Over D there are r
disjoint copies D_1,...,D_r, each mapping isomorphically to D. Remove the
whole source line Ltilde={z=0} and all D_k except one D_j. The resulting
open W_j is a smooth affine plane with coordinates (t,z), including the
WHOLE D_j={z=0}, and on its overlap with the source cover
```
y=lambda+z^r,
x=z^(1-r)*t+h(z),   h in C[z,z^-1].                     (2)
```
The Laurent polynomial is finite. The Hamiltonian extension of the
source translation lifts on W_j as
```
delta=z^(r-1)*partial_t.                                (3)
```
No unknown finite point of D is omitted from these charts.

If Omega has a regular primitive on M, which it does for the actual
canonical Keller component open, then
```
[h(z)]_(z^-r)=0.                                        (4)
```
After a polynomial t-shear preserving z and partial_t, the chart may
therefore be written
```
h(z)=sum_(q=r+1)^p a_q*z^-q,   a_p!=0,   p>=r+1.        (5)
```
This is a finite nonzero principal part. Without the primitive assumption,
the analogous shear still leaves a nonzero principal part, with p>=r;
the coefficient at z^-r is then not forced to vanish.

The integer p is a chart pole order, not a proved global bound on source
or target degrees. No claim about common-coordinate existence, target
Ga descent, unique LND kernels, generic deck cyclicity, or Keller inversion
is a conclusion of this conditional theorem.

## 2. Hamiltonian extension and verticality

The [Hamiltonian bridge](HAMILTONIAN_COMMON_COORDINATE_AND_ADDITIVE_COMPLETION.md)
applies to M directly: contraction with its unit form defines
```
contraction_delta(Omega)=dy,   delta|X=partial_x.
```
Its derivation on O(M) is locally nilpotent because O(M) embeds in
C[x,y] and every element has finite x degree. Its complete Ga action
preserves the whole source X and hence D. Also delta(y)=0.

The curve D must be vertical for y. Otherwise y|D is nonconstant.
The induced derivation on C(D) kills this nonconstant element, so is zero:
in characteristic zero C(D) is a finite separable extension of C(y|D),
and derivations extend uniquely through a separable algebraic extension.
Thus the ambient vector field delta vanishes generically along D.
The contraction identity then says dy vanishes generically as a cotangent
vector of M along D. Its pullback to the smooth generic point of D is
zero, forcing y|D to be constant in characteristic zero, a contradiction.

Write this constant as lambda. On X the zero divisor of y-lambda is
the single reduced source line L_lambda. No other component can occur
outside X, because M minus X=D. Thus(1) holds for some r>=1.
The source line is closed in affine M: its restriction algebra in C[x]
contains a nonconstant polynomial, over which C[x] is finite; hence its
image is closed. It is disjoint from D because it is contained in X.

## 3. Normalize the cyclic cover before using smoothness

Let Mtilde be the normalization of
```
Spec O(M)[z]/(z^r-(y-lambda)).
```
The polynomial is irreducible over C(x,y), by the valuation of y-lambda
on the source line. The normalized cover is therefore integral and finite
over M. It is affine. On the whole source X it is already the smooth
plane A2_(x,z), with y=lambda+z^r.

Near a point of D, choose a regular local defining equation f for the
Cartier divisor D. Because L_lambda and D are disjoint, (1) gives
```
y-lambda=u*f^r,   u a unit.
```
In the generic cover set w=z/f. Then w^r=u, and the local normalization
is exactly
```
O_M[w]/(w^r-u),   z=f*w.                                (6)
```
Indeed this algebra is finite etale over the regular local base
(r*w^(r-1) is a unit), hence normal; it has the correct generic field and
is integral over the raw cover. It is therefore its integral closure.
This proof does not assume f cuts out a smooth curve. The RAW equation
z^r=u*f^r need not be smooth when r>1, especially when D has a cusp.
Smoothness comes from its normalization(6).

Over the source line the raw cover is A2_(x,z), hence smooth. Away from
y=lambda it is etale. These loci together with (6) prove that Mtilde is
smooth everywhere.

The restriction above D is a finite etale cover of degree r. Locally
analytically its unit u admits an rth root, so it is a topological cover.
D is homeomorphic to C and hence simply connected. The cover splits as
r disjoint degree-one finite etale covers, algebraically isomorphic to D.
Call them D_1,...,D_r. Its zero divisor is reduced:
```
div_(Mtilde)(z)=Ltilde+D_1+...+D_r.                      (7)
```
All these supports are pairwise disjoint. The source line Ltilde is closed,
and each D_j is closed as a component of the finite inverse image of D.

For a fixed j remove Ltilde and the other r-1 curves from Mtilde.
The result W_j is smooth and affine. Each removed reduced curve is an
effective Cartier divisor on the smooth surface Mtilde. The complement
of their finite sum is affine over Mtilde: on an affine cover trivializing
the Cartier ideal its inverse image is a principal affine open.
Since Mtilde is affine, so is W_j.

There are no further points or fibres:
```
W_j minus D_j = Spec C[x,z,z^-1],
z^-1(0)=D_j reduced,
z^-1(c)=A1_x for every c!=0.                             (8)
```
The restricted map W_j->M is etale EVERYWHERE. The only ramification over
the source line was removed; (6) is etale at D_j, and the nonzero-z
locus was already etale. Thus Omega'=pullback(Omega) is a unit two-form
on W_j. No derivation-extension statement through normalization is
needed in the next step.

## 4. A divided Hamiltonian field removes the possible cusp

Define a regular vector field Delta on W_j directly by
```
contraction_Delta(Omega')=r*dz.                          (9)
```
The unit form makes it regular at every point. On the source overlap
Omega'=r*z^(r-1)*dx wedge dz, so
```
Delta=z^(1-r)*partial_x.                                (10)
```
This is exactly the divided field; (9) avoids assuming that dividing an
arbitrary extended derivation preserves a normalization.

Because
```
O(W_j) subset C[x,z,z^-1],
```
every element is a finite-degree polynomial in x with Laurent-polynomial
coefficients in z. The field (10) kills z; its iterates therefore kill
each element after finitely many steps. Hence Delta is a nonzero LND.
It preserves z and the reduced special fibre D_j.

At a generic point of D_j, z has order one by(7). In smooth local
coordinates its differential is nonzero, so (9) makes the induced
derivation on the curve D_j nonzero. A nonzero LND on ANY irreducible
affine complex curve forces its ring to be C[t], without a smoothness
classification: its kernel is C, because a nonconstant invariant would
make the curve function field finite separable over that invariant field
and force the derivation to vanish. For a moved element, take its last
nonzero derivative. This is a nonzero constant; the preceding derivative,
rescaled, is a slice t with Delta(t)=1. For any element, subtracting its
top derivative times t^n/n! repeatedly proves that it is a polynomial
in t. Thus D_j=A1 and its action is translation.

Since D_j->D is an isomorphism, D is algebraically A1 and is smooth.
In particular topology's previously possible cusps are ruled out by
the GIVEN regular source coordinate and unit volume.

## 5. Why the free action is a torsor and gives a whole plane chart

On z!=0, (10) is a nonzero translation on every A1_x fibre. On z=0,
Section4 gives a nontrivial translation on the whole A1=D_j. Thus
Delta has no zero anywhere on W_j. Equation(9) then makes dz nowhere
zero, so
```
f=z:W_j->A1
```
is a smooth surjective morphism, and every geometric fibre is exactly
one free Ga orbit.

These statements are stronger than freeness alone. They make the
action map
```
Ga times W_j -> W_j times_(A1) W_j,
(s,p) |-> (s*p,p)                                       (11)
```
an isomorphism. Over W_j both sides are smooth of relative dimension one;
the relative derivative of (11) is the nonzero vector Delta, hence (11)
is etale. Its geometric fibres are single points, because each fibre of f
is one free orbit. Thus it is an etale universally bijective morphism,
and therefore an isomorphism. Together with smooth surjectivity of f,
this is precisely the fppf Ga-torsor criterion.

A Ga torsor on the affine line is trivial:
```
H^1_fppf(A1,Ga)=H^1(A1,O)=0.                             (12)
```
The quasi-coherent comparison is
[Stacks03DW](https://stacks.math.columbia.edu/tag/03DW), and affine
vanishing is [Stacks01XB](https://stacks.math.columbia.edu/tag/01XB).
Their exact statements and displayed proofs were inspected and cached.
Equivalently the additive transition functions of local slices are a
quasi-coherent one-cocycle and hence a coboundary.

Choose the resulting global torsor coordinate t with Delta(t)=1.
This gives W_j=Spec C[t,z], the ENTIRE plane, including D_j={z=0}.
On its source overlap, (10) gives
```
Delta(x)=z^(1-r),   x-z^(1-r)*t in C[z,z^-1].
```
This is (2), with a finite Laurent polynomial h. The original source
translation is (3). The construction takes place in the original
normalized cover and includes every point of D_j.

## 6. The primitive kills exactly the logarithmic Laurent coefficient

Now assume Omega=dtheta for a regular one-form theta on M. This holds
in the actual canonical Keller application: for c=Jac(P,Q),
theta=c^-1*P*dQ is regular on M and dtheta=dx wedge dy.

On the whole source, theta-x*dy is a closed polynomial one-form, so
```
theta|X=x*dy+dF,   F in C[x,y].                          (13)
```
This polynomial Poincare statement is elementary: integrate its
dx coefficient in x; the remaining dy coefficient is a polynomial
only in y and can be integrated as well.

Pull (13) to the overlap in (2). It becomes
```
theta'=r*t*dz+r*h(z)*z^(r-1)*dz
       +dF(z^(1-r)*t+h(z),lambda+z^r).                  (14)
```
The left side is a regular polynomial one-form on the entire plane
W_j. Restrict to the whole curve t=0. Its dz coefficient is a polynomial
in z and has zero logarithmic residue at z=0. The last term is the
derivative of a finite Laurent polynomial after this restriction, so
its z^-1 coefficient also vanishes. The middle term therefore has zero
residue, which is exactly
```
r*[h(z)]_(z^-r)=0.
```
This proves(4). Exactness is used on the SAME component open, not merely
on the two separate plane charts. The fact that each plane has exact
two-forms by itself would not give(4).

## 7. Polynomial shears and the unavoidable nonzero principal part

Split the finite Laurent polynomial as
```
h=h_high+h_low,
h_high supported on exponents j>=1-r,
h_low supported on exponents j<=-r.
```
The polynomial coordinate change
```
t_new=t+z^(r-1)*h_high(z)
```
removes h_high from(2). It preserves the global plane chart, z, and
Delta=partial_t. With a primitive, (4) removes the remaining exponent
-r, so h_low has only exponents at most -r-1.

Crucially h_low is NONZERO. Otherwise, after that shear,
```
x=z^(1-r)*t_new.
```
There would be two morphisms A1_z->Mtilde: the source-cover curve x=0
and the boundary-chart curve t_new=0. They agree for every z!=0 by
this formula, but at z=0 their limits lie respectively on the disjoint
curves Ltilde and D_j. This contradicts separatedness of the affine
normalization Mtilde. Thus h_low!=0. Its largest pole order p gives(5).
This argument works for r=1 as well as r>1 and does not assume a
normal form for h before the polynomial shear.

Consequently an exact unit source volume rules out an empty principal
part or a lone order-r logarithmic term; it does not impose an upper
bound on p. For different D_i the individual source coordinates,
root choices, slices and principal parts need not align. No global
intersection-matrix or Xie-bound implication is inferred here.

## 8. Sharp hypothesis control and remaining obligation

The standard smooth affine family
```
M_p: y^p*a=b*(b-1),   p>=1,
b=x*y^p,   a=x*(x*y^p-1)
```
contains the whole source A2_(x,y), with complement
D={y=0,b=1}=A1. It has a unit source volume. Here r=1 and the whole
boundary chart, after removing the source line {y=0,b=0}, is
```
b=1+y^p*t,   a=t*(1+y^p*t),   x=t+y^-p.
```
For p=1 the coefficient h_-1 is nonzero; there is no global regular
primitive of the source volume on M_1. For every p>=2 a regular primitive
is
```
theta=(y*dx+p*x*dy)/(p-1).
```
The relation shows regularity at the whole D, and the Laurent residue
in Section6 explains precisely the distinction. These are controls
of the conditional geometry, not actual Keller target pairs or
counterexamples to JC2.

For the actual completion, the remaining input is still existence of
a suitable regular **source coordinate** on each component open, or
a single common coordinate on all of U. This theorem proves smooth
algebraic boundary and whole root charts once that input is given.
It neither produces the coordinate from valuation negativity nor
establishes compatibility of independently obtained coordinates.

The checker supplies exact local-normalization, LND/slice, torsor,
Laurent/shear/residue and whole affine-control calculations.
The unbounded theorem is the written proof, not its finite models.
The retained construction/topology and Hamiltonian bridge, plus the two
cached cohomological primary statements, are pinned in the receipt.
No older proof or status metadata is changed by this bundle. Internal
checking is not formal verification or external human review.
