# The first deficit: exact monodromy and a determinant component criterion

Use the notation of the
[Hermite threshold](POLYNOMIAL_STABLE_HERMITE_CONNECTIVITY.md).
There are \(k\) non-transposition finite factors, each of index
at least two, \(s\ge1\) literal transposition factors, and \(p\)
nontrivial cycles among the first \(k\) factors. Assume
\[
          \prod g_i=P_N,\qquad \sum\operatorname{ind}(g_i)=N-1,
          \qquad s=p-k-1\ge1.                                \tag{1}
\]
This is one simple branch below the proved threshold.
All statements concern full braids with every factor-type
placement allowed; absolute equivalence quotients by powers
of the fixed product.

The actual group is \(S_N\), with exactly one possible
imprimitive passport family. That family's entire imprimitive
part is one explicitly counted component. In all remaining
cases, the unresolved component count reduces to the
geometric factorization of one explicit moment determinant,
followed by a finite label action. This is a reduction, not
a claim that every remaining determinant is irreducible.

## Primary comparison

Jones--Zvonkin,
[Section 5.4, retained pp155--159](sources/jones-zvonkin-braid-cacti-2002.pdf),
construct imprimitive cacti from copies of a quotient and
attachment data, including the degree-fifteen example on
pp157--158. Their actual construction and branch-placement
counts were inspected. It does not give the unrestricted
first-deficit classification below.

Kochetkov,
[arXiv:1212.6317v1, Sections 3 and 6](sources/kochetkov-zolotarev-trees-2012.pdf),
uses coefficient varieties for three-critical-value polynomials,
including explicit low-degree reducible examples. The actual
eight-page version was inspected. That is useful comparison
for coefficient-space methods, not an imported assertion of
geometric irreducibility or of a general component count here.

The algebraic identities and conormal calculation in K1,
K5 and K6 of the Hermite note are used before its numerical
threshold is imposed; the pertinent argument is repeated
below. Its stable-range connectedness conclusion is not
assumed at the present deficit. Other proof inputs are H2--H3 of the
[polynomial tree note](POLYNOMIAL_PASSPORT_COMPONENTS.md),
and L5--L6 of the
[return-transport note](POLYNOMIAL_PASSPORT_LIFT_TRANSPORT.md).
The positive class count is the corrected primary
Goulden--Jackson formula recorded and audited in K1's source
discussion. No new external component theorem is assumed.

## Theorem U1. The only imprimitive possibility at the first deficit

Under (1), \(k\ge2\). Every imprimitive tuple has \(s=1,k=2,p=4\)
and, for some \(a,b\ge2\), has degree \(N=2m\), \(m=a+b-1\),
and passport
\[
 (a^2\,1^{2b-2}),\qquad(b^2\,1^{2a-2}),\qquad T.              \tag{2}
\]
Its normal transposition blocks have size two. Every tuple
outside this imprimitive part generates the actual \(S_N\),
has only block-count divisors \(1,N\), and has trivial deck group.

**Proof.** If \(k=1\), its only branch has \(I+p\) prescribed
preimages counted with multiplicity, but \(N=I+s+1=I+p-1\),
which is impossible. Use the normal transposition blocks of
L6, of size \(B\ge2\), with quotient degree \(d>1\).
Let \(c\) count the nontrivial quotient cycles, and let
\(I_o,p_o,z\) be the total index, nontrivial-cycle count and
number of factors acting trivially on the quotient.
If \(A\le c\) factors move quotient blocks and \(R\) is
the return index over their nontrivial quotient cycles,
then \(p-k\ge cB-R-A+p_o-z\), with \(p_o\ge z\) and
\(R\le B-1-s-I_o\). Thus the return-budget inequality is
\[
             p-k\ge(c-1)(B-1)+s+I_o.                        \tag{3}
\]
If \(c=1\), all but one of the \(k\ge2\) factors fix quotient
blocks, and its sharper equality gives \(p-k\ge s+2\).
This contradicts \(p-k=s+1\). Thus \(c\ge2\).
Equation (3) forces \(c=2,B=2,I_o=0\).
Every non-transposition factor is nonidentity, so none can
fix the quotient: \(k\le c=2\), hence \(k=2\).
The two-point return budget is one and already contains
the \(s\ge1\) transpositions, so \(s=1\) and all other
returns are identities. Each quotient factor has exactly
one nontrivial cycle, of lengths \(a,b\ge2\).
The quotient index budget gives \(d=a+b-1=m\).
Lifting with identity returns duplicates these cycles
and their fixed points, giving (2).

L6 says that a tuple with one normal transposition block
generates \(S_N\); its block and deck assertions follow.
\(\square\)

## Theorem U2. The complete imprimitive component

For every \(a,b\ge2\), passport (2) has exactly one imprimitive
full component. Put
\[
 \epsilon=2-\boldsymbol1_{a=b},\qquad W=3\epsilon,\qquad
 K=\begin{cases}
       A_m,&a,b\ \text{both odd},\\
       S_m,&\text{otherwise}.
    \end{cases}                                             \tag{4}
\]
The component has actual group \(S_2\wr K\), all block-count
divisors \(\{1,m,2m\}\), deck order two, literal size \(Wm^2\),
and absolute size \(Wm\).

Every other component of (2) has group \(S_{2m}\) and deck
one. Their combined literal size is
\[
               Wm^2\bigl((2a-1)(2b-1)-1\bigr).              \tag{5}
\]
Their separate number and sizes are not asserted uniformly.

**Proof.** The quotient consists of an \(a\)-cycle and a
\(b\)-cycle on \(m=a+b-1\) sheets. Its two cycle vertices
in the incidence tree share exactly one sheet. Their
cyclic orders give a unique absolute two-factor cactus;
reversing the cycle-type order adds the second placement
exactly when \(a\ne b\). Its literal size is \(m\epsilon\).

For completeness, its group is the \(K\) in (4).
Write its cycles as \(A,B\), meeting at \(x\).
Their commutator is a three-cycle with support
\(\{x,A(x),B(x)\}\). Conjugating by powers of \(A\)
gives all triangles consisting of \(B(x)\) and adjacent
vertices of the \(A\)-cycle. These triangles generate the
alternating group on that cycle together with \(B(x)\):
triangles along a connected edge graph with one fixed
apex generate the alternating group. To verify this, start
with one triangle and add vertices along a spanning tree.
An alternating group on a set \(X\), together with a
three-cycle containing two vertices of \(X\) and one new
vertex, generates the alternating group on the enlarged
set: conjugate that three-cycle inside \(A_X\), also using
its inverse when \(\lvert X\rvert=3\).
If \(a\ge3\), this contains the alternating group on
the \(A\)-cycle. Conjugating it by powers of \(B\) gives
alternating groups on sets containing the \(a-1\ge2\)
other \(A\)-vertices and each \(B\)-vertex in turn.
Their overlaps generate \(A_m\).
When \(a=2<b\), exchange \(A,B\) in this argument.
When \(a=b=2\), the two transpositions generate \(S_3\)
directly. The parity of the two original cycles now gives
(4). In particular the quotient is primitive.

Apply L5 with the additional free cyclic branch of degree
two. The transposition uses the entire fibre budget, so
all other returns are identities. L5 gives exactly one
lifted component, group \(C_2\wr K=S_2\wr K\), literal
size \(3m(m\epsilon)=Wm^2\) and deck two.
There is no other imprimitive component by U1.
Every nonsingleton block system coarsens the normal
two-point blocks; the primitive quotient \(K\) gives
exactly the stated block divisors.

The ordered class count of (2) is
\[
             (2m)^2\frac{2b-1}{2}\frac{2a-1}{2}
                         =m^2(2a-1)(2b-1).
 \]
There are \(W\) distinct factor-type placements. Subtract
the unique imprimitive component to obtain (5).
The remainder is positive. U1 identifies every tuple
in it as \(S_{2m}\). \(\square\)

## Theorem U3. One determinant indexes all first-deficit components

Label all nontrivial critical points as in K1, put
\[
 Q(z)=\prod_{i,j}(z-u_{ij})^{q_{ij}},\qquad
 F_{\ell i}(u)=\sum_j\operatorname{res}_{z=u_{ij}}
                                      \frac{z^\ell\,dz}{Q(z)}
                    \quad(0\le\ell\le k-1),                 \tag{6}
\]
and let \(T\) consist of the first \(k-1\) rows of \(F\).
Remove its last column and set
\[
                              \Delta(u)=\det T_{\rm red}(u).
                                                               \tag{7}
\]
These are rational functions on the distinct-point
configuration space. Clearing their denominators there
does not change any locus below.

Retain the geometric irreducible factors of \(\Delta\) which
meet the following open conditions:

- \(\operatorname{rank}F=k-1\) and
  \(\operatorname{rank}T=k-2\);
- \(Fv=(0,\ldots,0,1)^t\) has its unique solution with \(v_k=0\);
- the resulting monic polynomial \(\Phi_v\) satisfies every
  exact ramification and critical-value exclusion (9) of K1.

Let \(H\) permute equal local-degree critical points within
each branch, and also permute branch colors having the same
complete partition. Then full Hurwitz components are exactly
the \(H\)-orbits of these retained geometric factors.
In particular \(\Delta\) is not identically zero, and no
additional component is supported on a lower-rank locus.

This is an exact algebraic component criterion. A uniform
factorization or closed formula for the number of retained
factor orbits is still required to turn it into a numerical
classification of all passports (1).

**Proof.** Here \(h=p-s-2=k-1\). Every row of \(F\) kills
the constant value vector, so its rank is at most \(k-1\).
At a feasible rank \(r\le k-2\), write
\(\delta=k-r\ge2\). The conormal calculation in K5--K6
does not require \(k\ge h+2\): its feasible incidence
dimension bound is still
\[
                              p-\delta+1<p.                 \tag{8}
\]
K5's admissible Jacobian has full row rank \(h+1=k\),
so the good coefficient locus has pure dimension
\(p+k-k=p\). Consequently no component can live at
rank \(r\le k-2\). Every component meets rank \(F=k-1\).

The identities behind that dimension bound can also be
checked directly in this range. For a vector \(w\) with
\(Fw\in\mathbb C e\), CRT gives \(g=\Phi_w\) of degree
at most \(N\), with \(g'=ZB_g\), \(\deg B_g\le s\).
Differentiating its jets at a moving critical point gives
\[
 \beta(dF)w=-\sum_j
                    \frac{\beta(u_j)B_g(u_j)}{P'(u_j)}\,du_j,
 \qquad P(z)=\prod_j(z-u_j),\quad \deg\beta\le k-1.
                                                               \tag{D}
\]
The differentiated denominator has no \(z^{-1}\) term,
since \(\deg g\le M-k\). For an admissible monic \(f\),
\(B_f=NA\) is nonzero at every prescribed point, so (D)
is a surjective Vandermonde Jacobian. At feasible rank
\(r\), left-kernel polynomials have dimension
\(\delta=k-r\) and degree at most \(k-2\), and the
space of \(B_g\)'s has dimension \(k-r=\delta\).
Their products have degree at most \(p-3\) and dimension
at least \(2\delta-1\), by the leading-degree product
argument. Their \(p\)-point evaluations are injective
and conormal to the rank stratum. Adding the solution
fibre of dimension \(\delta\) gives (8). This repeats
the algebraic proof without a stable-range assumption.

Feasibility at rank \(F=k-1\) forces rank \(T=k-2\): if \(T\)
had rank \(k-1\), its only kernel would be the constant
vectors, all killed by the final row as well. Removing
one column from a matrix whose rows kill constants
preserves row rank, so \(\Delta=0\).
Conversely rank \(T=k-2\), rank \(F=k-1\) and \(\Delta=0\)
give a unique nonconstant kernel direction modulo constants;
the final row normalizes it to \(Fv=e\).
The choice \(v_k=0\) supplies a regular section on this
rank open. All other solutions add the same constant
to every value. All good conditions are unchanged by
this translation.

The exact good full-rank coefficient locus is therefore
the product of an affine line with the indicated open
subset of \(\Delta=0\). Each irreducible component of the
full good locus meets this open, and the omitted ranks
have smaller dimension. Smoothness precludes component
intersections in the good locus. Thus the retained
geometric factors give precisely its labeled components.
Nonemptiness follows from the positive cactus count.
The same dimension statement also precludes
\(\Delta\equiv0\): a nonempty good rank open would then
have coefficient dimension \(p+1\), contrary to the
Jacobian dimension \(p\).

Forgetting the finite labels gives the \(H\)-orbit
description. Source and target affine groups are connected,
so they preserve each component; they do not introduce
further mergers. The Hurwitz correspondence and full-twist
argument of K2 identify these components with the stated
absolute and literal full-braid orbits. \(\square\)

## Proposition U4. The degree-six determinant separates the two groups

For \(a=b=2\), normalize the two critical points in the first
branch to \(0,1\) and write the others as \(x,y\). Put
\[
 U=x+y-1,\qquad V=xy,\qquad
 P=-3V^2+U(2U+3)V-(U+1)U^2.
 \]
On the distinct-point configuration,
\[
 \Delta=\operatorname{res}_0\frac{dz}{Q}
             +\operatorname{res}_1\frac{dz}{Q}
         =\frac{2UP}{V^3(V-U)^3},\qquad
 Q=z^2(z-1)^2(z-x)^2(z-y)^2.                                \tag{9}
\]
The two admissible geometric factors are \(U\) and \(P\).
The first is exactly the imprimitive component; the second
is the \(S_6\) component. Their literal sizes are \(27,216\)
and their absolute sizes \(9,36\).

**Proof.** Direct differentiation at the two double poles
gives
\[
 \Delta=2\left(\frac{V+x+y}{V^3}
                -\frac{V-2(x+y)+3}{(V-U)^3}\right),
 \]
whose numerator factors as in (9).
On \(U=0\), both pairs have midpoint \(1/2\). The CRT
interpolant is invariant under \(z\mapsto1-z\), hence
is a polynomial in \((z-1/2)^2\), giving the quadratic
deck symmetry. Conversely U1--U2 show that every
imprimitive polynomial has this symmetry and hence these
paired midpoints.

The other factor is geometrically irreducible. On \(P=0\)
in the configuration open, \(U\ne0\); set \(w=V/U\).
Then
\[
 U=\frac{3w^2-3w+1}{2w-1},\qquad V=wU,\qquad
 (x+y)^2-4xy
       =-\frac{w(w-1)(15w^2-15w+4)}{(2w-1)^2}.              \tag{10}
\]
The quartic numerator in the last expression is squarefree.
It is not a square in \(\mathbb C(w)\), so ordering the two
roots \(x,y\) gives one irreducible quadratic cover of
the rational \(w\)-curve. This proves geometric irreducibility
of \(P(x,y)\). The two factors meet only at \(U=V=0\),
which is outside the configuration open.
The positive total count and U2 show that both have
nonempty good loci. U3 gives two components, and U2's
count gives the asserted sizes. \(\square\)

## Proposition U5. Complete bounded first-deficit controls

For every \(2\le a\le b\le5\), passport (2) has exactly two
components: U2's imprimitive component and one \(S_{2m}\)
component, whose size is (5). This is a bounded exact
statement; the single primitive component is not extrapolated
to all \(a,b\).

Two further first-deficit passports have one component:
\[
 \begin{array}{c|c|c|c}
 N&\text{finite passport}&\text{actual group}&
                         \text{full literal size}\\ \hline
 8&\mathrm{DT},\mathrm{DT},3,T&S_8&38\,400\\
 9&\mathrm{DT},\mathrm{DT},\mathrm{DT},T,T&S_9&1\,771\,470 .
 \end{array}                                                \tag{11}
\]
Here \(3\) denotes one three-cycle with all other sheets fixed.

**Proof and exact evidence.** The
[checker](verify_polynomial_first_deficit.py) generates complete
fixed-pattern absolute braid closures using pure braid
generators, exchanges of equal-type factors, and fixed-product
rotations. These generate the full colored braid group.
For the family controls the two closures are disjoint, their
normal transposition graphs and exact quotient closures give
the asserted actual groups, and literal multiplicities
\(N/\lvert\mathrm{Deck}\rvert\) exhaust an independent ordinary
hook-character total. The two additional closures exhaust
their complete class totals in the same way. The retained
[receipt](polynomial_first_deficit_verification.json) specifies
every seed, domain, orbit hash, group, block, deck and count
convention. The [source ledger](polynomial_first_deficit_source_ledger.json)
separates the uniform written reduction from the bounded
component assertions. \(\square\)

Run from the canonical checkout:

    python3 research/arithmetic/verify_polynomial_first_deficit.py

## Remaining obligation

The subcases \(\min(a,b)=2\) and \(\min(a,b)=3\) of (2)
are now classified, with exact separate sizes, in the
[DT/duplicated-cycle theorem](POLYNOMIAL_DT_DUPLICATED_CYCLE_COMPONENTS.md)
and the [quadratic residue theorem](POLYNOMIAL_DUPLICATED_THREE_CYCLE_COMPONENTS.md).
The latter also gives a residue polynomial of degree
\(\min(a,b)-1\) for arbitrary lengths.
The entire duplicated-cycle component problem, including
\(a,b\ge4\), is closed by the independent
[four-core path normalization](POLYNOMIAL_DUPLICATED_CYCLE_PATH_NORMALIZATION.md):
there is exactly one primitive component for all lengths,
in addition to U2's unique imprimitive component. Its
separate size is therefore the aggregate (5).
The [general two-cycle-pair classification](POLYNOMIAL_TWO_CYCLE_PAIR_COMPONENTS.md)
also resolves unequal within-color lengths and, with the
single-cycle-tail theorem, the complete \(k=2,s=1,p=4\)
sector. Its conclusions hold under pure braids on each
fixed ordered passport as well as under full braids.

For arbitrary passports (1), factor the determinant (7)
geometrically, remove factors with empty good locus, compute
the finite label action, and determine each component's
critical-value-map degree. U1--U2 already determine every
actual monodromy group, block system and deck group at this
level. Uniform primitive connectedness for the other
profiles has not been proved.
The remaining first-deficit profiles include \(k=2,p\ge5\)
with \(s=p-3\ge2\), and \(k\ge3\), outside already
resolved families.
In particular rational factorization alone does not certify
geometric factorization. The general component goal, all
deeper deficits and zero-simple passports outside already
resolved families remain open. Independent review and
replay of the present arguments remain pending.
