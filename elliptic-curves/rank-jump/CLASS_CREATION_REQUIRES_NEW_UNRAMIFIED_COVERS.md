# Class creation: new unramified covers beyond the global pool

Follow-up: [other generic fibres do not supply the retained jump blocks](OTHER_GENERIC_FIBRES_DO_NOT_SUPPLY_THE_JUMP_BLOCK.md)
tests 112 equation-only fibre correspondences. Rational matches are global
duplicates; a quadratic sign obstruction kills the entire rational trace
span on every nonrational component.

The large-jump problem now has a more specific target than an enlarged
collision statistic. Several successful fibres necessarily have an
additional **unramified quadratic-cover block of their cubic number field**
outside the entire inherited global pool. On fresh 103b2 at 3726/881 its
dimension is at least **seven**; on historic ICARM356 it is at least **nine**.
The extensions split at every bad place. Extra ramification at rational
primes is therefore not what these strict dimensions introduce.

These are new deductions from existing rank lower bounds and equation-only
capacity certificates. They are not independently constructed extensions,
newly measured class groups, or a condition predicting which t succeeds.
The class-creation mechanism remains open. The bounded point-independent
principalization experiment below supplies no new class.

## The two capacity bounds locate the missing structure

For a smooth rational specialization t, let H=H1(Q,E_t[2]), W⊂H be the
rational Kummer image, and G⊂W the marked generic subgroup, of dimension m.
Let F⊂H be the specialization of the global pool L from the
[root-curve capacity theorem](ROOT_CURVE_TORSION_AND_REAL_CAPACITY.md).
Thus G⊂F and dim F≤m+q with **q=2** for every panel family. No assertion
that all of F lies in the arithmetic Selmer group is needed.

Let U be the strict Selmer space: the cubic Kummer classes are locally
square over 2, infinity and every bad rational prime, and unramified
elsewhere. Put k=dim(G∩U). The
[strict-boundary certificate](FRESH_STRICT_BLOCK_NECESSITIES.md) bounds
the full localized Selmer image by h=m−k+a. Consequently, if rank E_t≥R,

\[
 \dim(W\cap U)\ge\max(k,R-h),\qquad
 \dim(F\cap U)\le k+q.
\]

The second inequality follows from the injection
(F∩U)/(G∩U)→F/G. Subtracting gives the sharper joint necessity

\[
\boxed{
 \dim\frac{W\cap U}{W\cap U\cap F}
 \ge\max(0,R-m-a-q)=\max(0,J-a-2).
}
\tag{1}
\]

This isolates strict classes outside **all of L**, not just outside G.
The argument remains valid if specialization on L has a kernel.

| Successful fibre | Retained gain J | Boundary excess cap a | All rational directions outside F ≥ | Strict rational directions outside F ≥ |
|---|---:|---:|---:|---:|
| 103b2, 3726/881 | +10 | 1 | 8 | **7** |
| 11952, −2448/11 | +10 | 5 | 8 | **3** |
| 11952, 110314/102227 | +10 | 4 | 8 | **4** |
| 11952, 2828/2015 | +10 | 2 | 8 | **6** |
| ICARM356, historic R17 | +12 | 1 | 10 | **9** |
| ICARM385, historic R17 | +12 | 8 | 10 | **2** |
| ICARM398, MW16 | +14 | 9 | 12 | **3** |

The [sixteen-row artifact](../../artifacts/generated-results/elliptic-curves/rank_jump_strict_cover_creation_capacity_v1.json)
also retains every matched low and incomplete high. All completed lows
have lower bound zero in (1), not a proof of absence. Five strict-boundary
rows remain UNKNOWN, including the fresh MW16 +11; its total outside-pool
lower bound remains nine. No missing factorization was retried.

## What the additional arithmetic objects are

Write K_t=Q(theta) for the irreducible cubic field and S_K for primes
above S={2, infinity, bad primes}. Class field theory identifies

\[
 U_t\simeq\operatorname{Hom}
 \left(\operatorname{Cl}(\mathcal O_{K_t,S_K}),\mathbb F_2\right).
\tag{2}
\]

A nonzero class is a quadratic extension K_t(sqrt(alpha))/K_t unramified
at all finite places and split at S_K and at the real places. A block of
r independent classes has compositum of degree 2^r over K_t. Its extra
directions cannot be accounted for by further splitting or a CT switch
inside F. In (1) one may choose representatives independent modulo F;
they are also independent as cubic squareclasses.

There is an equivalent quartic description. A nontrivial norm-square
class of a non-Galois cubic K gives an S4-quartic field Q_alpha with cubic
resolvent K. Their maximal-order discriminants satisfy

\[
 \operatorname{Disc}(Q_\alpha)
 =\operatorname{Disc}(K)\,
 N_{K/\mathbb Q}\mathfrak d(K(\sqrt\alpha)/K).
\]

Thus strict classes give quartic fields of **the same field discriminant**
as K. At each p∈S their local quartic algebra is Q_p×(K⊗Q_p).
This uses the norm-square/quartic correspondence and discriminant identity
of Cohen–Thorne,
[Theorem 2.2](https://link.springer.com/article/10.1007/s40993-015-0001-y).
The local assertion follows by making the squareclass trivial over K⊗Q_p:
the four-element torsor becomes the origin together with the three roots
of the cubic. It also applies at the real place.

Accordingly the required event can be stated without exceptional points:
**new independent S4 extensions with the same cubic resolvent and no added
relative ramification, satisfying the bad-place splitting conditions**.
Polynomial discriminants do not suffice; an index in a nonmaximal order
can introduce spurious discriminant factors.

This identifies the required object, not its cause. There is no canonical
map comparing the class groups of two different K_t in the matched panel.
Here “new” means outside F_t. It does not assert a class-group
specialization map from the neighbouring number field or a measured
increase in absolute class-group dimension. Nor does the capacity theorem
prove that every extra direction is strict; the last two table columns
deliberately distinguish those statements.

## A concrete sufficient incidence certificate at t

An independent class-creation construction can be certified by the following
data, using only the equation and generic classes:

1. Elements alpha_i∈K_t and fractional ideals I_i with the exact identities
   (alpha_i)=I_i^2 and N(alpha_i) a positive rational square.
2. Local square certificates for alpha_i at every prime above S and every
   real place. The ideal identity already removes odd-prime ramification
   outside S; the dyadic square tests are indispensable.
3. An Artin-character matrix at auxiliary good prime ideals, of enough rank
   to prove independence. To prove independence modulo G, append all the
   generic classes in the same squareclass-character calculation.

Obtaining k+3 independent strict classes proves that at least one lies
outside F, because dim(F∩U)≤k+2. More generally k+2+r prove at least r
dimensions of class creation outside F. These are sufficient **incidence**
certificates, irrespective of the retained rank labels. None proves the
classes rational: local solubility is built in, and rationality versus Sha
still requires the second stage.

The missing condition on t is an equation-defined construction of the
I_i and their simultaneous square principalizations, with the local and
independence certificates. Restating “the class quotient has enough
2-rank” would not supply that condition. No such specialization criterion
has been obtained here. A general norm-square class or a principal ideal
square alone is also insufficient: it can be trivial, inherited, or fail
the strict local tests.

## Bounded constructor test: no additional class

The frozen [protocol](PRIME_SQUARE_CLASS_CONSTRUCTOR_PROTOCOL.json) tests
the first 24 residue-degree-one prime ideals in the existing equation-only
early relation pool of the supplemental MW16-05, t=3/17 reference. This
is a calibration on a retained +6 example, not a new matched panel or a
prospective parameter search. Its additional strict block is known to
exist from rank accounting, but neither those points nor their derived
classes enter this constructor.

For each P, form I=P^2, reduce its positive trace lattice, and test all 49
primitive directions in [-2,2]^3 up to overall sign. An element alpha∈I
with |N(alpha)|=N(I) proves (alpha)=I by equal ideal index. The same test
runs on the principal ideals of the first three generic Kummer classes
as controls. The entire worker has a 30-second cap and a 256 MiB maximum
PARI stack. There are no extra ideals, enlarged boxes, adaptive weights,
class-group computations or elliptic point searches.

The [completed result](../../artifacts/generated-results/elliptic-curves/rank_jump_prime_square_class_constructor_v2.json)
finds **zero** candidate generators and one generator for each of the
three controls. Each control generator equals its original generic
Kummer representative coefficient by coefficient. The
[independent replay](../../artifacts/generated-results/elliptic-curves/rank_jump_prime_square_class_constructor_verification_v1.json)
checks all **1323 norms** with rational multiplication matrices, and checks
the exact ideal bases and target norms. The initial v1 UNKNOWN artifact
is preserved: its worker stopped before enumeration because PARI's
nfbasistoalg requires a column vector. The completion fixes only that
interface error; bounds and selection are unchanged.

The null result excludes just these lattice vectors. It proves neither
that P^2 is nonprincipal nor that the strict class group is small. Also,
a class of order two need not have a representative among these 24 prime
ideals. Recovering all controls checks this limited generator **visibility**;
it does not validate a detector of class-group rank. No independent new
Kummer class, additional CT entry, or solubility theorem results.

Replay the narrow calculations:

```sh
timeout 30 sage -python elliptic-curves/rank-jump/verify_prime_square_class_constructor.py check
python3 elliptic-curves/rank-jump/strict_cover_creation_capacity.py check
```

## Priorities after this test

1. **Incidence:** construct the unramified, S-split excess in (1). The most
   decisive matched test is still 103b2 at 3726/881 versus −1049/2296:
   both have k=0, yet the successful +10 needs seven strict dimensions
   beyond F, despite a smaller extra boundary cap (one versus two).
   A point-independent construction of three independent strict classes
   is already enough to cross the global-pool capacity on either fibre.
   It is a falsifiable first milestone, not an assertion that the low
   control lacks those classes.
2. **Incidence, alternative construction:** ramified covers over Q(t)
   whose ramification disappears upon specialization. Their geometric
   obstruction rank must meet the existing 8/9/10/12 capacity deficits;
   an explanation restricted to the inherited two-dimensional remainder
   cannot succeed. Base-place ramification and ramification at rational
   primes must not be conflated.
3. **Solubility, after construction:** compute the additional block's CT
   form and remaining obstruction to rationality. The six-direction
   fixed-incidence twist remains the control for this stage. Vanishing
   CT is necessary, not sufficient in general.

The individual-prime-square extractor is weak in this frozen test. Total
collision support, generic cover counts and inherited CT switches are
already inadequate explanations; expanding their correlation panels is
not the next action. Agent1 receives no new scoring feature. The usable
output is a precise necessary class supply and a proof gate for recognizing
its construction. The implication from a condition on t to that supply,
and then to simultaneous rational solubility, is still missing.
