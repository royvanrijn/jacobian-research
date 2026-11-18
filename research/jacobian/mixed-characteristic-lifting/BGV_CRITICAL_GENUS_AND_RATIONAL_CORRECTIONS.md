# Positive-genus critical curves exclude rational corrections of the BGV seed

For every odd prime p, the characteristic-zero integer BGV seed has one
geometrically irreducible critical curve. Its smooth projective genus is

    p^2+p-(3+gcd(2p-1,p-2))/2.

This curve is not contracted and has ramification index two. It obstructs
every dominant rational source correction, even when accompanied by an
arbitrary dominant rational change on the target. Thus a different
rational volume normalization of the fixed integer seed cannot produce a
plane Keller map. This is stronger than the radial obstruction, but does
not exclude deformations of the seed or algebraic corrections whose
nonrational terms cancel. It does not resolve JC2, and no novelty or major
impact claim is made.

## 1. A general divisorial obstruction

Work over an algebraically closed field k of characteristic zero. Let
F:A2 -->> A2 be dominant and generically finite. Suppose a prime curve E
in its source has positive geometric genus, F is defined generically on E,
and F maps E dominantly onto a curve Gamma. Let e>=2 be the ramification
index of the corresponding divisorial valuations.

**Theorem.** For any dominant rational maps alpha:A2 -->> A2 and
beta:A2 -->> A2, the rational map beta composed F composed alpha cannot
have a nonzero constant Jacobian.

Neither alpha nor beta is assumed birational. The theorem concerns a
rational Jacobian identity, so it also excludes any polynomial Keller
map obtained by these operations.

Resolve alpha on the projective source plane by finitely many point
blowups pi:S->P2, obtaining a morphism a:S->P2. Such a resolution is the
classical elimination-of-indeterminacy theorem; a primary proof appears
in [Brosnan's surface notes, Theorem 4.24](https://math.umd.edu/~pbrosnan/notes/Surfaces/sect0035.html).
The morphism a is proper, dominant and generically finite. Some irreducible
component D of a^(-1)(Ebar) dominates the projective closure Ebar: otherwise
properness would fail over the generic point of Ebar. Normalize both
curves and pass to their smooth projective models. The induced curve
map is finite and separable. Riemann--Hurwitz implies that D has positive
genus, since its target has positive genus.

All exceptional curves of pi are rational, as is the strict transform
of the original line at infinity. Therefore D is neither exceptional nor
a curve at infinity. Its generic point maps isomorphically to the generic
point of a curve in the original affine source plane. In particular,

    ord_D(pi*(dx wedge dy))=0.                           (1)

This accounts for indeterminacy points inside A2 as well as those at
infinity. One cannot assume that a resolution preserves A2; instead its
additional exceptional curves are excluded by their genus.

For completeness, the differential order formula used next can be checked
locally. If one prime divisor lies above another with ramification index r,
choose a uniformizer t downstairs and u upstairs. Then t=a*u^r with a a
unit. Choose a rational separating parameter along the downstairs residue
curve. Its differential has nonzero pullback along the upstairs residue
curve because the residue extension is separable. Hence for any rational
two-form omega of order m downstairs,

    ord_upstairs(omega pulled back)=r*m+r-1.             (2)

The leading term is r*a*u^(r-1) du wedged with a nonzero residue
differential. Characteristic zero makes its coefficient nonzero.
Thus (2) includes zeros and poles and does not require the form or the
rational map to be regular away from the generic divisors.

Put n=ord_Gamma(beta*(dX wedge dY)), an integer. Applying (2) first
to F and then to alpha, and writing r for the ramification index along D
over E, gives

    ord_D((beta F alpha)*(dX wedge dY))
       =r*(e*n+e-1)+r-1
       =r*e*(n+1)-1.                                    (3)

The right side cannot be zero, because r*e>=2. A nonzero constant
Jacobian would instead make this order zero by (1). This proves the
theorem.

The genus hypothesis is essential. For G=(x^2,y), the ramified curve
x=0 is rational and is not contracted. The dominant rational correction
alpha=(x,y/(2x)) gives G composed alpha=(x^2,y/(2x)), whose rational
Jacobian is one. Thus ramification alone cannot replace the positive-genus
hypothesis; its divisors can lie on rational exceptional curves.

The only general geometric imports are resolution by point blowups and
Riemann--Hurwitz in characteristic zero. For the latter, see
[Stacks, Section 53.12, Lemmas 53.12.2 and 53.12.4](https://stacks.math.columbia.edu/tag/0C1B). The same source gives the
local tame differential calculation underlying (2). The source record
identifies the inspected statements; internal checking is separate from
external review.

## 2. The fixed integer BGV seed and its critical equation

For an odd prime p, now take the same formulas over Q:

    z=x+x^p*y,
    P=z+(xz)^p,
    Q=y+x^(p(p-1))*z^(p(p-1)+1),
    F0=(P,Q),                 J=Jac(P,Q).                (4)

The positive-characteristic origin and exact BGV theorem are recorded in
[the lifting proof](BGV_ODD_PRIME_LIFTING.md). Here all calculations,
critical curves and genera are in characteristic zero.

In the birational coordinates x,z, direct differentiation gives

    x*J=(1-p)*x+p*z+p*(p+1)*x^p*z^p
        -p*(p-1)*x^(p+1)*z^(p-1)
        -p*(p-1)*x^(p^2)*z^(p^2-p+1)
        +p*x^(p^2+p)*z^(p^2).                            (5)

Set u=xz. Substituting x=u/z and multiplying (5) by z^p gives

    u*z^(p-1)*J
       =(1-p)*u*z^(p-1)+p*z^(p+1)
        +p*(p+1)*u^p*z^p-p*(p-1)*u^(p+1)*z^(p-2)
        -p*(p-1)*u^(p^2)*z+p*u^(p^2+p).                 (6)

This display is an identity in the function field, with J evaluated after
x=u/z and y=(z-x)/x^p. Finally put

    N=2p-1,       v=z/u^p,
    H_p(v)=v^(p+1)+(p+1)*v^p-(p-1)*v+1.

Equation (6) becomes

    u*z^(p-1)*J
       =u^(p^2-p+1) *
          [p*u^N*H_p(v)-(p-1)*v^(p-2)*(v+p)].            (7)

Consequently the critical curve is birational to

    u^N=R_p(v),
    R_p(v)=(p-1)*v^(p-2)*(v+p)/(p*H_p(v)).               (8)

The inverse functions are z=v*u^p, x=u/z, y=(z-x)/x^p.
No component of J=0 is lost on x=0 or z=0: respectively J=1 and
J=1-p there. The birational correspondence therefore accounts for every
critical component once irreducibility of (8) is proved.

## 3. Squarefreeness and all branch points

**Lemma.** For every odd prime p, H_p is squarefree, H_p(0)=1, and
H_p(-p)=p^2-p+1-p^p is nonzero.

The two evaluations are immediate; p^p>p^2-p+1 for p>=3. To prove
squarefreeness, define

    T_p(v)=(p-1)*v^2+(p+1)*(p-2)*v-(p+1).

The exact polynomial identity

    (p+1)*(v+p)*H_p-v*(v+p+1)*H_p'=-p*T_p              (9)

shows that a repeated root of H_p must be a root of T_p. With v=-1+w,

    T_p(-1+w)=(p-1)*w^2+p*(p-3)*w-p*(p-1).             (10)

This is Eisenstein at p. In a p-adic algebraic closure every root w has
valuation 1/2, with v_p(p)=1. On the other hand,

    H_p(-1)=0,       H_p'(-1)=p*(p-1),
    H_p(-1+w)=w^(p+1) mod p.                            (11)

Thus the linear term of H_p(-1+w) has valuation 3/2 at such a root.
Every term of degree k>=2 and k<=p has valuation at least 1+k/2>=2,
and the top term has valuation (p+1)/2>=2. The linear term is the
unique term of least valuation, so H_p(-1+w) is nonzero. No common root
of H_p,H_p' exists. This proves the lemma uniformly, including p=3.

At v=-p the rational function R_p has valuation one. Hence the Kummer
polynomial u^N-R_p(v) is Eisenstein at that point over the algebraic
closure of Q, and is irreducible of degree N. This proves geometric
irreducibility. It is reduced, since its u-derivative N*u^(N-1) is
nonzero at its generic point. Thus J has one geometrically irreducible
critical component, with multiplicity one.

The divisor of R_p has precisely these places and orders:

| Place | Number of places | Order of R_p |
| --- | ---: | ---: |
| v=0 | 1 | p-2 |
| v=-p | 1 | 1 |
| Roots of H_p | p+1 | -1 |
| v=infinity | 1 | 2 |

There are no collisions among these places. For a local equation
u^N=t^a*unit in characteristic zero, let d=gcd(N,a). Taking roots of
the unit in the completed algebraically closed residue field gives d
points above the place, each of ramification index N/d. The total
ramification contribution is therefore N-d.

Put d0=gcd(2p-1,p-2). The contribution at zero is N-d0; all other
listed places have coprime order and contribute N-1 each. Riemann--Hurwitz
therefore gives

    2g-2=-2N+(N-d0)+(p+3)*(N-1),

or

    g=p^2+p-(3+d0)/2.                                   (12)

Here d0 divides 3, so g is p^2+p-2 or p^2+p-3. In particular it is
positive for every odd prime. The first values are g=10 at p=3 and
g=27 at p=5.

## 4. Noncontraction and rational corrections

On the critical function field, P=u^p*(v+1). Therefore

    P^N=R_p(v)^p*(v+1)^N.                               (13)

The right side has positive order p at v=-p, where v+1 is a unit.
It is not identically zero. Hence P is nonconstant on the critical
curve, which is not contracted by F0.

The multiplicity-one zero of J and formula (2) give ramification index
e=2 over its image curve. The theorem in Section 1 applies.

**Corollary.** For every odd prime p and every pair of dominant rational
maps alpha,beta of the plane over a characteristic-zero field, the rational
map beta composed F0 composed alpha does not have nonzero constant
Jacobian. In particular it cannot be a polynomial Keller map.

For nonclosed coefficient fields apply the argument after geometric base
change. Dominance and the differential identities persist. No bound is
placed on the degrees of either rational correction.

This conclusion concerns the fixed characteristic-zero seed (4). Adding
p-multiples to its coefficients changes the critical curve and falls
outside the statement. An algebraic correction whose coordinates do not
belong to the rational source function field also falls outside it.
The previous radial construction is one such algebraic correction;
its nonrational output obstruction requires its separate proof.
Thus there is no claim that every possible polynomial lift of the
positive-characteristic BGV map is excluded.

## 5. Exact controls and scope

The [checker](verify_bgv_critical_genus.py) reuses the standard-library
sparse arithmetic in the existing lifting checker. For p=3,5,7 it checks
the original Jacobian, both birational critical identities, (9)-(11),
the branch orders and the numerical Riemann--Hurwitz calculation.
It rejects a deliberately corrupted critical identity.

An independent Singular 4.3.2 normal.lib computation of the p=3 curve

    3*u^5*(v^4+4*v^3-2*v+1)-2*v*(v+3)=0

returns geometric genus 10. A rational-curve control returns zero.
The exact Singular input, version, normal.lib digest and actual output
are retained in the [receipt](bgv_critical_genus_verification.json).

These finite controls do not prove squarefreeness or the genus formula
for all p. Sections 1-4 supply the complete written argument, with the
two stated classical imports. The checker does not verify the general
surface-resolution theorem or establish a JC2 breakthrough.

Replay from the repository root:

    python3 research/jacobian/mixed-characteristic-lifting/verify_bgv_critical_genus.py

The replay invokes Singular for the separate p=3 genus calculation.
Use --write-receipt PATH only to record an intentional new receipt.
