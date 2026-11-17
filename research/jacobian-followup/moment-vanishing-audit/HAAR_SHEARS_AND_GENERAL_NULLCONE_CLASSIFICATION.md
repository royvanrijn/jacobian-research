# Haar shears and the general moment nullcone classification

15 November 2026. Internal written proof combining the repository's
[Gaussian shear mechanism](GAUSSIAN_SHEAR_NORMAL_FORM_AND_TRANSFERS.md)
with the invariant-functional setting of
[Long's spherical-pair draft](sources/long_spherical_pairs_e9edac2.tex).
The construction below re-proves moment preservation algebraically and
does not import a Gaussian contour theorem, Long's radial witness, or the
root-doublet branching argument. Classical compact forms, highest-weight
theory, complete reducibility, affineness for reductive stabilizers, and
the torus constant-term theorem remain inputs.
The [source audit](haar_shear_source_audit.json) identifies the selected
principal-subgroup source and the exact finite certificates.
The [checker](verify_haar_shear_witnesses.py) verifies concrete examples;
the general proof is not formally verified or independently specialist-reviewed.
No novelty or priority claim is made.

**Main result.** For connected complex reductive G and reductive H,
with no sphericity assumption, pure Reynolds moments detect the G-nullcone
in C[G/H] exactly when G_der is contained in H. Whenever they fail to
detect it, there is a pure-null element outside the nullcone whose every
fixed marked sequence nevertheless eventually vanishes. Thus these new
witnesses do not disprove the Mathieu--Zhao property.

The retained [Rodrigues survivor](../../archive/non-elliptic/extended-geometry/TWO_PAIR_SIC_BIDEGREE33_RODRIGUES_SURVIVOR.md)
already separates pure-nullity, nullcone membership and mixed safety for
a related contraction functional. The result here is the uniform
homogeneous-space construction and the complete Reynolds nullcone criterion.

The construction also removes the spherical assumption from the earlier
disconnected nullcone criterion. The remaining exceptional questions concern
the global Mathieu--Zhao property, not moment detection of the nullcone.

## 1. A functional-preserving replica of a unipotent derivation

Let A be a commutative complex algebra with a rational G-action, and let
h:A -> C be G-invariant. Let F be the locally nilpotent derivation arising
from a unipotent one-parameter subgroup. If q is in ker F, then

    D=qF,      D^j=q^j F^j,      Phi_t=exp(tD)              (1)

is locally nilpotent and Phi_t is a polynomial algebra automorphism for
each t in C, with inverse Phi_(-t). Moreover

    h(Da)=h(F(qa))=0,
    h(Phi_t(a))=h(a)                   for every a in A.   (2)

The exponential is finite on each element, so (2) needs no analytic change
of contour or claim that Phi_t maps a compact real form to itself.

For P_t=Phi_t(p),

    h(P_t^m)=h(p^m),
    h(Q P_t^m)=h(Phi_(-t)(Q) p^m).                         (3)

Therefore eventual marked vanishing for every fixed Q is transported
bijectively by Phi_t. This is the Haar analogue of the retained Gaussian
shear argument. Such a transformation cannot turn a safe element into
a Mathieu--Zhao failure witness.

## 2. A real self-dual representation is available on every active quotient

Let G be connected reductive, H reductive, X=G/H, and A=C[X].
Choose compatible compact forms L inside K of H inside G. Complex
conjugation on K/L extends to an antilinear involution star on A:
conjugate matrix coefficients on K are coefficients of the dual rational
representation, and right-L-invariance extends to H by Zariski density.
In particular A contains the dual of every irreducible type that occurs.

Suppose G_der is not contained in H. Since X is affine and regular
functions separate points, A has an irreducible type U_lambda on which
G_der acts nontrivially. Let lambda*=-w_0 lambda denote the dual highest
weight, and choose nonzero highest vectors of types lambda and lambda*
in A. Their product is nonzero, because A is a domain, and is a highest
vector of type

    mu=lambda+lambda*,           mu*=mu.                   (4)

This type is nontrivial on G_der and trivial on the connected center.

It has a compact real structure. Indeed the Cartan component of
U_lambda tensor U_lambda* is the unique irreducible component of highest
weight mu. Identify that tensor product with End(U_lambda), with K acting
by conjugation after choosing an invariant Hermitian product. Hermitian
adjoint is an antilinear involution commuting with K. It preserves the
Cartan component: the latter is self-dual and has multiplicity one.
Its restriction is therefore a real structure J on a model V_mu.

Occurrence of type mu in A implies V_mu^H is nonzero. Explicitly,
evaluation at eH on a nonzero function module is an H-fixed vector in
its dual and cannot vanish identically, by translation and transitivity.
Since mu is self-dual, this gives the asserted fixed space in V_mu.

The space V_mu^H equals V_mu^L and is J-stable. Choose a nonzero J-real
vector xi in it. Matrix coefficients of V_mu with this xi form an
irreducible, star-stable copy U of type mu in A. Hence a highest vector
p in U has

    q=p^star in U

as a lowest vector. Both are regular functions on X, and on K/L,

    q=conjugate(p),
    C_s=h((pq)^s)=integral_(K/L) |p|^(2s)>0       (s>=1).   (5)

Strict positivity follows from nonzero p, Zariski density, continuity,
and full support of compact Haar measure.

This argument avoids assuming that an arbitrary function-module copy
is star-stable, or that an arbitrary self-dual representation has real
rather than quaternionic type.

## 3. Principal weights and the separating invariant

Choose a principal sl_2 triple of derivations E,F,H_0 compatible with
the highest-weight convention. Its semisimple element is 2rho^vee.
One explicit construction takes

    E=sum_i e_i,        H_0=2rho^vee=sum_i c_i h_i,
    F=sum_i c_i f_i,

using simple Chevalley generators. The relations
[H_0,E]=2E, [H_0,F]=-2F and [E,F]=H_0 follow directly.
The cocharacter 2rho^vee is integral in the coroot lattice.
The chosen F belongs to the negative nilpotent algebra and generates
an algebraic unipotent action, so its derivation on A is locally nilpotent.

These principal conventions are recorded in
[Seppanen--Tsanov, Section 1.2](https://arxiv.org/abs/1503.07105v1).
Only this classical input is used from that paper, not its GIT quotient
or unstable-locus theorems. Products of simple factors are handled by the
diagonal triple, and a central cover does not change the argument.

Put

    n=mu(H_0)>0.

Because H_0 is strictly dominant, U has a unique highest H_0-weight n
and a unique lowest weight -n, at its G-extremal vectors p and q.
Its principal sl_2 restriction therefore has a unique top string of
highest weight n. In particular F^n p is a nonzero multiple of q.
Rescale F (and E inversely) so that

    Fq=0,             F^n p/n!=q.                          (6)

Apply Section 1 with this q. Then

    P_t=Phi_t(p)=sum_(j=0)^n (t^j/j!) q^j F^j p,          (7)

and its H_0-weights are

    n-j(n+2),                j=0,...,n.                   (8)

Every pure moment of p is zero by its nonzero weight, so (3) proves
all positive pure moments of P_t vanish for every t. The element p
also has eventual marked vanishing: its positive weight dominates the
finitely many weights of every fixed marker. Equation (3) transfers
that property to P_t, uniformly on any fixed finite-dimensional
marker space.

It remains to prove P_t is outside the G-nullcone for t!=0.
Let pi_mu and pi_((n+1)mu) be the G-isotypic projections on a finite
G-module containing P_t. They are well-defined by complete reducibility.

A copy of type mu has principal weights between -n and n. Thus (8)
gives

    pi_mu(P_t)=p+t v,                 H_0 v=-2v.          (9)

All terms with j>=2 have weight below -n. The j-th term of (7) is in
the product of j+1 copies of U. That tensor product has no principal
weight greater than (j+1)n. Consequently for j<n it cannot contain
type (n+1)mu, whose highest principal weight is n(n+1).

The last term of (7) is t^n q^(n+1). It is a nonzero G-lowest vector
of weight -(n+1)mu and is entirely of type (n+1)mu. Therefore

    pi_((n+1)mu)(P_t)=t^n q^(n+1).                       (10)

The invariant homogeneous polynomial, of degree n+2 in a module element f,

    I(f)=h((pi_mu f)^(n+1) pi_((n+1)mu) f)

now gives

    I(P_t)=t^n h((pq)^(n+1))=t^n C_(n+1).                (11)

Indeed each occurrence of v in (9) lowers the total weight by n+2,
so all such terms have zero Reynolds value. Equation (5) makes C_(n+1)
strictly positive. Thus I(P_t)!=0 whenever t!=0, whereas I vanishes
at zero. This proves P_t is outside the G-nullcone.

No converse from marked safety to instability is used. In fact (3) and
(11) show that converse is false, even for connected acting groups.

## 4. Complete nullcone criteria

For connected G, define

    Z_X={f:h(f^m)=0 for every m>=1},
    N_X={f:0 is in the closure of G.f}.

**Theorem 1.** Let G be connected complex reductive and H closed reductive,
possibly disconnected. Then

    Z_X=N_X  iff  G_der is contained in H.                 (12)

If G_der is not contained in H, (7) gives a one-parameter family with
vanishing pure moments and eventual marked vanishing for every parameter,
but it meets N_X exactly at t=0.

The preceding construction proves the negative direction. For the positive
direction G/H is an algebraic torus; Duistermaat--van der Kallen and strict
character separation identify pure-nullity with a contracting cocharacter.
The reverse inclusion N_X subset Z_X always holds because every pure
moment is an invariant homogeneous polynomial.

The full ordinary pure-moment ideal of the parameter family (7) is zero
in C[t], while its set-theoretic nullcone inverse image is {0}.
Equation (11) gives a concrete invariant separating every other parameter.
It is not asserted to generate the full scheme-theoretic nullcone ideal.

**Theorem 2.** Let G have reductive identity component K and let H have
reductive identity component. Put D=[K,K]. No componentwise spherical
hypothesis is required. If D is contained in H, form the component tori
and character map exactly as in the
[earlier component theorem](SPHERICAL_PAIRS_AND_COMPONENT_NULLCONES.md):

    H_i=K intersect g_i H g_i^(-1),
    theta: direct_sum_i X^*(K/H_i) -> X^*(K/D),
    theta((chi_i))=sum_i chi_i.

Then

    Z_X=N_X  iff  D is contained in H and theta is injective. (13)

If D is not contained in H, apply Theorem 1 on K/(K intersect H) and
extend the witness by zero to its open-and-closed component of X.
The normalized functional is the equal average of component functionals.
The G-nullcone equals the K-nullcone because G/K is finite, so the
extended witness stays outside it and remains eventually marked-null.

If D is contained in H, the retained finite-trace theorem GVC2FTS and
Duistermaat--van der Kallen give componentwise strict character separators.
A simultaneous separator exists for all pure-null functions precisely when
theta is injective. A nonzero relation in ker theta gives the monomial
counterexample from the earlier component proof. Those examples are
eventually marked-null as well.

Thus both obstructions in (13) have witnesses that do not themselves
violate the Mathieu--Zhao property. That global property remains a
different question.

## 5. Explicit witnesses on both exceptional quotients

Use coordinates x,y,z,w on SL_2 as in the earlier notes and the derivation

    F=z*d_x+w*d_y.

These substitution derivations form an sl_2 triple with
E=x*d_z+y*d_w and H_0=x*d_x+y*d_y-z*d_z-w*d_w.
They exponentiate to the same SL_2 orbit action after a Weyl change of
convention; reversing the usual function-action cocharacter changes none
of the conclusions.

Choose the following binary forms in a standard compact embedding:

    f_8(X,Y)=X^8+14X^4Y^4+Y^8,                  Gamma=2O,
    f_12(X,Y)=X*Y*(X^10+11X^5Y^5-Y^10),         Gamma=2I.

Put p=f_n(x,y), q=f_n(z,w). Both forms are invariant under the displayed
finite groups and obey q=conjugate(p) on SU_2. The explicit family is

    P_t=f_n(x+t q z, y+t q w).                              (14)

It is regular and right-Gamma-invariant. The inverse shear replaces t
by -t, preserves the determinant, and induces the automorphism (1).

The checker supplies exact generator matrices in Q(zeta_8) and Q(zeta_5),
respectively. It verifies determinants, compact conjugation relations,
group orders 48/120, projective element orders, invariance of the forms,
and fixed dimensions through n. It identifies the same conjugacy types as
the earlier binary-polyhedral models. A conjugate subgroup receives the
witness by the corresponding right translation; no explicit conjugator
to the earlier quaternion orientation is needed or claimed.

For t!=0 the exact minimum matrix-entry degree is

    deg(P_t)=n(n+1):                 72 for2O, 156 for2I.

The top term t^n q^(n+1) has substitution weight -n(n+1), which cannot
occur in a representative of lower matrix-entry degree.

For these concrete modules the checker also verifies a smaller stable
projection. With r=n/2, the V_n component is

    p+t kappa F^(r+1)p,

where

| n | kappa | Corresponding binary form |
|---:|---:|---|
| 8 | 7/4290 | X^8+t(1568/143)X^3Y^5 |
| 12 | -11/3907008 | X^12-t(3630/323)X^5Y^7 |

For t!=0 the largest root multiplicities are 3 and 5, strictly below
n/2=4 and6. This gives an independent non-nullcone test.

The calculation uses Omega=d_x*d_w-d_y*d_z and the identity

    Omega^r(q Fp)
       = r! (n+r+1)!/(n+1)! * kappa F^(r+1)p.

The scalar before normalization is 672 for n=8 and -39600 for n=12.
The checker independently applies the left-Casimir projector and obtains
the same component multiplied by det^r before imposing det=1. The
Fischer decomposition is valid because the finite fixed spaces of all
types 1,...,n-1 vanish; the weight -2 excludes the trivial type.

For completeness, any finite subgroup of SL_2 can be treated without its
classification. Conjugate it into SU_2, take the product of the translates
of a nonzero linear form, and multiply it by its compact-real conjugate
binary form. After a phase adjustment this gives a nonzero invariant real
form of even degree n=2|Gamma|. The same shear and invariant (11) work.
The two forms above give much smaller degrees for the exceptional cases.

## 6. Exact marked towers and sharp cutoff bounds

For (14) at t=1 and every s,m>=1, (3) gives

    h(q^s P_1^m)=delta_(m,s) C_s,
    C_s=h((pq)^s)>0.                                      (15)

If f_n(X,Y)^s=sum_i c_i X^(ns-i)Y^i, then

    C_s=(1/(ns+1)) sum_i |c_i|^2/binom(ns,i).               (16)

In particular C_1=8/15 for2O and25/1092 for2I.

Every marker Q of matrix-entry degree at most e satisfies

    h(Q P_1^m)=0 whenever n*m>(n+1)e.                     (17)

Under the inverse shear, each matrix coordinate has minimum substitution
weight at least -(n+1). Thus Phi_(-1)(Q) has minimum weight at least
-(n+1)e, while p^m has weight nm.

This coefficient is sharp. For Q=p^s of exact degree e=ns, the last
possible index m=(n+1)s has

    h(p^s P_1^((n+1)s))=C_((n+1)s)>0.                     (18)

Only the last lowering term in Phi_(-1)(p^s) has the weight needed to
pair with p^((n+1)s); its sign is positive because n is even. Hence the
sharp slopes are 9/8 and13/12. These are cutoff bounds for a fixed pure-null
function, not a failure of eventual vanishing: the marker changes with s.

The checker directly integrates the first weight-allowed pure moments,
m=5 for2O and m=7 for2I, by the determinant-one Haar monomial formula.
Both cancel exactly. This is independent finite support for the all-order
proof (2), not the basis for claiming infinitely many pure moments vanish.

## 7. What is closed and what remains

Equations (12) and (13) close the full moment-nullcone classification for
complex reductive homogeneous pairs, including disconnected stabilizers
and disconnected acting groups. The proof of the connected statement
does not require sphericity or the special radial identity.

The Mathieu--Zhao classification is not closed by these witnesses. Their
eventual marked vanishing is proved in (3), (15) and (17). The previous
exact reduction of the remaining global Mathieu question to failure on
SL_2/2O and SL_2/2I still applies, but only its Mathieu--Zhao half remains
undecided. A future failure witness needs an infinite nonzero marked tail;
being outside the group nullcone is now explicitly insufficient.

The construction is motivated by the repository's Gaussian-preserving
shears. It does not claim that a nonlinear automorphism preserves the
ordinary Laplacian or supplies a new Jacobian-conjecture implication.
Those transfer boundaries remain as recorded in the Gaussian shear note.
