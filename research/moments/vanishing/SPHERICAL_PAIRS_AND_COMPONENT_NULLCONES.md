# Mathieu properties and moment nullcones for disconnected spherical pairs

15 November 2026. Internal written deduction from
[Long's connected spherical classification](https://github.com/long-mathematics/mathieu-second-conjecture-spherical-pairs/blob/e9edac262b4b64eb68feb0e6e80496cefe4a1dee/mathieu_second_conjecture_spherical_pairs.tex)
and the repository's [finite-trace theorem](../../foundations/extended-geometry/BINARY_GVC_UNIFORM_FACE_TERMINATION.md).
The connected classification is an imported research-draft theorem. The
extension below has no independent specialist review or formal verification;
no novelty or priority claim is made. The [source audit](spherical_pairs_source_audit.json)
pins the imported source and records its proof dependencies.
The [exact checker](verify_spherical_components.py) supplies bounded circuit
and lattice controls, rather than an all-group verification.

The Mathieu--Zhao property extends to disconnected acting groups using the
same condition on the connected semisimple action. Equality between the
pure-moment zero set and the group nullcone needs an additional condition:
the character spaces of the different components must be jointly independent.
Thus the two properties in Long's connected theorem separate.

**Later result.** The [Haar-shear theorem](HAAR_SHEARS_AND_GENERAL_NULLCONE_CLASSIFICATION.md)
now closes moment-nullcone detection for all reductive homogeneous pairs,
without sphericity. Its witnesses are eventually marked-null. The two
exceptional global Mathieu--Zhao questions remain open.

## 1. Groups and the component character map

Work over C. Let G be a linear algebraic group with reductive identity
component K=G^0, and let H be a closed subgroup with reductive identity
component. Set X=G/H and D=[K,K]. Assume the connected homogeneous space

    X_0 = K/(K intersect H)

is spherical for K. We call this the componentwise spherical hypothesis.
Conjugation makes every other K-orbit spherical as well. There need not be
one Borel orbit dense in the disconnected space X.

Let g_1=1,...,g_r represent G/(KH). The K-orbits

    X_i = K g_i H/H = K/H_i,
    H_i = K intersect (g_i H g_i^(-1))

are precisely the connected components of X. They are affine, open and
closed. Indeed H_i has reductive identity component and finite component
group, so the usual reductive homogeneous-space affineness applies.
The orbit decomposition is finite, and each orbit is connected, since K is
connected. Consequently

    A = C[X] = product_i C[X_i].

Let h:A -> C be the normalized G-invariant Reynolds functional. Define

    Z = {f : h(f^m)=0 for every m>=1},
    N = {f : 0 belongs to the Zariski closure of G.f}.

The closure is taken in a finite-dimensional rational submodule containing f.
When D is contained in H, it is contained in every H_i: D is characteristic
in K and hence normal in G. Put

    T = K/D,                  T_i = K/H_i,
    M = X^*(T),               M_i = X^*(T_i) subset M.

Each H_i is normal in K, each T_i is a connected algebraic torus, and X_i
is identified with T_i by the chosen base point g_iH. Disconnected H_i may
make M_i nonsaturated; this causes no problem. The maps K -> T_i induce

    theta : direct_sum_i M_i -> M,    (chi_i)_i |-> sum_i chi_i,
    Theta : T -> product_i T_i.

Write L_i=M_i tensor Q. The following conditions are equivalent:

* theta is injective;
* the subspaces L_i form an internal direct sum in M tensor Q;
* dim_Q(sum_i L_i)=sum_i dim_Q(L_i);
* Theta is surjective.

The equivalence uses character duality for tori. The character groups are
free, so rational injectivity is equivalent to integral injectivity.
Surjectivity of a torus morphism does not require surjectivity on the
integral cocharacter lattices; a finite-index obstruction disappears after
multiplying a rational cocharacter.

## 2. Classification theorem

**Theorem.** Under the componentwise spherical hypothesis:

1. ker h is Mathieu--Zhao if and only if D is contained in H.
2. Z=N if and only if D is contained in H and theta is injective.
3. If D is not contained in H, fixed P,Q in A satisfy

       h(P^m)=0,       h(Q P^m)>0       for every m>=1.

   Thus P is outside N and the Mathieu property fails.
4. If D is contained in H but theta is not injective, there is an explicit
   P in Z outside N, constructed from any nonzero relation in ker theta.
   Nevertheless, every fixed Q has h(Q P^m)=0 for all sufficiently large m.

For connected X there is just one M_i, and its inclusion in M is injective.
Thus the two criteria agree whenever X is connected, even if G is
disconnected. For connected G the theorem recovers Long's classification.

## 3. Reynolds integration and finite component traces

Let h_i be the normalized K-invariant functional on C[X_i]. The restriction
of h to a component summand is a scalar multiple of h_i by uniqueness.
G permutes the components transitively and transports these functionals, so
the scalars agree. Normalization gives the exact formula

    h(f) = (1/r) sum_i h_i(f_i).                         (1)

In particular every component idempotent has Reynolds value 1/r. A regular
function on a component extends by zero to X, because the component is
open and closed; this is an algebraic zero extension.

Suppose D is contained in H. The functionals h_i are constant term on C[M_i].
Embed all C[M_i] in C[M]. These embeddings preserve constant terms, since
M_i -> M is injective. The finite-trace theorem GVC2FTS says

    sum_i CT(f_i^m)=0 for every m>=1
       implies CT(f_i^m)=0 for every i,m.                (2)

For clarity, its arithmetic input is all-order, rather than positivity of
the complex values CT(f_i^m). After finite-type specialization preserving
a chosen nonzero moment, work over a number field. Fix m and 1<=k<=r.
At arbitrarily large completely split good primes p, bounded signed digits
and Frobenius give

    CT(f_i^[m(1+p+...+p^(k-1))])
       = CT(f_i^m)^k modulo a prime over p.

The bound p>max absolute coordinate in every support of f_i^m eliminates
carries. The trace identity supplies the first r power sums of the component
m-th moments at infinitely many primes. They therefore vanish in
characteristic zero. Newton identities force every component moment zero.
This is the precise use of the retained finite-trace proof.

Duistermaat--van der Kallen now gives, for each nonzero f_i,

    0 is outside conv(Supp(f_i)) in M_i tensor R.        (3)

Choose a strictly positive rational linear functional on each such support.
For any fixed Q_i, the weights of Q_i f_i^m become strictly positive for
large m. Its constant term is then zero. Taking the maximum of the finitely
many component thresholds and using (1) proves the Mathieu property.
This direction needs no sphericity once D is contained in H.

If D is not contained in H, the connected pair (K,K intersect H) satisfies
Long's hypotheses and has nontrivial semisimple action. Its fixed witnesses
p,q have h_0(p^m)=0 and h_0(q p^m)>0 for every m>=1. Zero extend them from
X_0 to X. Formula (1) divides both sequences by r and proves assertion 3.

## 4. Nullcones require one simultaneous separator

First N is contained in Z for every such G: f |-> h(f^m) is an invariant
homogeneous polynomial on each finite-dimensional submodule.

Moreover the G-nullcone equals the K-nullcone. Choose representatives a_j
of G/K. Then

    closure(G.f) = union_j a_j closure(K.f).

This is a finite union of closed sets. Since each a_j fixes zero, zero
belongs to this union exactly when it belongs to closure(K.f).
Hilbert--Mumford for the connected reductive group K supplies a
one-parameter subgroup in every unstable case. Positive weights of f
eventually dominate the finite weights of any fixed Q. Hence instability
implies eventual marked vanishing, just as in Long's connected proof.
In particular the witnesses of assertion 3 are outside N.

Assume D is contained in H. Write the actual torus weights of all nonzero
component terms of f in the single ambient lattice M. With translation
action on functions the weight of the monomial chi is -chi; negating every
weight together changes none of the following separation conclusions.
The torus nullcone criterion is

    f in N iff there exists one rational linear functional on M tensor Q
       strictly positive on the union of all nonzero component supports. (4)

It is essential to retain a single functional in (4). Moment vanishing only
gives separate functionals in (3).

If theta is injective, extend the separate rational functionals from the
direct sum of the L_i to M tensor Q. It is strictly positive on every
component support at once. Scaling gives an integral cocharacter of T.
A further positive multiple lifts to a cocharacter of K: the connected
central torus Z(K)^0 maps onto K/D with finite kernel. Reversing the
cocharacter for the translation convention sends f to zero. Thus Z subset N.

Conversely, if theta has nonzero kernel, choose a nonzero integral relation

    chi_i in M_i,              sum_i chi_i=0.

Clear denominators first if the relation was computed over Q. Set P_i to
the character monomial chi_i when chi_i is nonzero, and set P_i=0 otherwise.
Each nonzero monomial has every positive constant term zero, so P is in Z.
A common strictly positive functional would give

    0 = ell(sum_i chi_i) = sum_(chi_i != 0) ell(chi_i) > 0,

a contradiction. Thus P is outside N by (4). This construction proves
necessity and assertion 4. The already proved Mathieu property supplies
eventual vanishing for every fixed marker Q.

Combining both cases proves every assertion of the theorem.

## 5. Examples that distinguish the criteria

**Shared character line.** Take G=C^* times C_2 and H=1. Then X has two
torus components, M_1=M_2=Z inside M=Z, and theta(a,b)=a+b.
The pair P=(z,z^(-1)) is pure-null and stable away from zero. Indeed on
its two-dimensional coefficient module the product of the coefficients
is G-invariant and equals 1. The functional is

    h(f_1,f_2) = (CT(f_1)+CT(f_2))/2.

Its kernel is Mathieu--Zhao by the theorem. This directly disproves
equivalence of Mathieu--Zhao and moment-nullcone equality for disconnected
acting groups, even with toral identity component.

**Independent component directions.** Let G=(C^*)^2 semidirect C_2 with C_2
swapping the two coordinates, and take H=ker(e_1) inside the identity torus.
The two character lattices are Z e_1 and Z e_2, so theta is injective.
X is disconnected and positive-dimensional, but Z=N and ker h is
Mathieu--Zhao. Neither connectedness nor a zero-dimensional quotient is
necessary.

**Pairwise independence is insufficient.** Let T=(C^*)^2 and let C_3 act on
its character lattice by the order-three matrix with columns

    A = [[0,-1],[1,-1]].

Take G=T semidirect C_3 and H=ker(e_1) inside T. The three component
lattices, in some order, are

    Z e_1,        Z e_2,        Z(-e_1-e_2).

Every pair of their rational spans intersects trivially, but their sum is
not direct. The three monomials in these directions give a pure-null P
outside N. The criterion must test the whole direct-sum map.

**Group quotients.** If H=1, the spherical hypothesis restricts the connected
identity group to a torus. For a virtually toral G, ker h is always
Mathieu--Zhao; Z=N holds exactly when G is connected or finite. In positive
torus dimension every extra component repeats the full character space.

## 6. What the linked source adds and what remains open

Long's connected theorem allows disconnected stabilizers and all central
forms. Its negative proof has three parts: every active simple factor of
rank at least two exposes an actual defining root doublet; a line--doublet
projection transfers the radial moment witness; the residual SL_2 factors
use an explicit Weyl-even SL_2/N(T) circuit or diagonal blocks. Sphericity
is used only to exclude finite SL_2 stabilizer projections in that residual
argument. The source's numerical weight 1 warning is valid: occurrence of
that weight alone is not a defining two-dimensional summand.

The imported normalizer circuit is

    U=u^2, V=v^2, S=tau^2,       UV=(1-S)^2,
    R=S(U+2)(2-S(U^2+2U+2)),
    P_N=(1+U)(V-R),              Q_N=U.

The all-order coefficient formula is

    h_N(B(U) P_N^m) = c_(2m) [U^m] B(U)(1+U)^(m-1),
    c_j=4^j (j!)^2/(2j+1)!.

The checker independently expands the regular polynomial circuit and
integrates balanced monomials exactly on the sphere. It also verifies the
component lattice relations, the order-three action, finite support
separators and the shared-line invariant. These are transcription and
boundary controls; the general Lie and orbit-closure steps use the written
proof and imported theorem.

For nonspherical homogeneous spaces the finite SL_2 projection sector is
still outside the connected classification used here. The source already
settles the active higher-rank sector without sphericity. This extension
does not settle the remaining nonspherical sector, minimal witness degrees,
positive characteristic, Gaussian GMC in higher degree, or any Jacobian
or ordinary-Laplacian conjecture. In particular, an eventual mixed cutoff
and a cocharacter contracting the whole tuple are different conclusions.

The [later nonspherical reduction](NONSPHERICAL_PAIRS_AND_EXCEPTIONAL_SL2.md)
extends the classification to every connected stabilizer and to component
groups without A4/S4/A5 quotients. The unrestricted connected remainder
reduces to two exceptional models; their moment decisions remain open.
