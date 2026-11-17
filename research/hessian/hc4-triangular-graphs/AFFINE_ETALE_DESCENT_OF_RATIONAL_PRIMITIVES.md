# Affine etale descent removes the polynomial-representative obstruction

HC4MYMW0. 16 November 2026. Characteristic zero. Written all-degree proof.
Full HC4/JC2 and the original finite trace-dual/marked-energy implication
remain OPEN. No novelty, independent replay, formal verification or external
review is asserted.

This corrects the interpretation of HC4MYMV0, while retaining its proved
two-condition equivalence and its graph-family exclusion. Under MV0's ACTUAL
smooth affine, SURJECTIVE ETALE source-image hypotheses, its polynomial
representative condition follows from its rational primitive condition.
It is therefore not an additional open algebra obstruction in that setting.
The frozen MV0 proof, checker and receipt are preserved as historical evidence.

## 1. The general additive algebra statement

Let k be a field of characteristic zero. Let

    q:X=Spec(R) -> V=Spec(C)

be a SURJECTIVE ETALE morphism of smooth integral affine k-varieties, with
function fields L/E. This map is generically finite and separable; it need
not be finite. Regard all the following differentials as subgroups of
Omega^1_(L/k). Then

    {g in L : dg in Omega^1_(R/k) + Omega^1_(E/k)} = R + E.       (1)

The sums in (1) are ADDITIVE sums. No claim about integral closure, finite
modules, multiplicative orders or trace-dual stability is part of (1).

### A regular differential cannot hide a function pole

On a smooth integral characteristic-zero variety T, a rational function f
whose differential df is regular is itself regular. At the generic point
of any divisor, write f=t^(-a)*b, a>0, with b a unit. The divisor is
generically smooth over k, and there is a regular derivation D near that
generic point with D(t) a unit. Hence

    D(f)=t^(-a-1)*(-a*b*D(t)+t*D(b))

has order -a-1. The nonzero coefficient a cannot disappear in characteristic
zero. This contradicts regularity of df. Thus f has no divisorial poles;
normality extends it over all points. The argument applies separately to
the finitely many disjoint integral components of a smooth scheme.

### Regular differences on the actual source relation

Take g with dg=alpha+eta, alpha regular on X and eta rational over E.
Set T=X times_V X, with projections p1,p2. Both projections are etale.
Thus T is smooth, and every nonempty component dominates X: its image is
open and dense. In particular g has well-defined rational pullbacks on
every component. Put

    z=p2^*g-p1^*g.

The rational E-form eta has equal pullbacks, so

    dz=p2^*alpha-p1^*alpha

is REGULAR on the WHOLE relation T. By the pole argument z is a regular
function on T, not just an element at its generic points. On the triple
fibre product it satisfies the exact additive cocycle identity

    z23-z13+z12=0.                                            (2)

All these statements concern the actual fibre products of the specified
source morphism. No arbitrary generic-sheet identifications or missing
finite source points have been substituted.

### Affine faithfully flat descent supplies a regular representative

Since C -> R is faithfully flat, the augmented Amitsur complex

    0 -> C -> R -> R tensor_C R -> R tensor_C R tensor_C R -> ...

is exact. With d0(r)=1 tensor r-r tensor 1, its degree-one differential is
d1(z)=z23-z13+z12. Consequently (2) gives an r in R with

    z=p2^*r-p1^*r.                                            (3)

For completeness, exactness follows without assuming module finiteness.
After the faithfully flat base change C -> R the complex has one extra
first tensor factor. A contracting homotopy multiplies that factor by
the first ordinary tensor factor. In the alternating unit-insertion
differential, its first term gives the identity and the remaining terms
cancel those from the opposite composition. This proves dh+hd=1 in
every degree, including the augmented degree. Faithful flatness reflects
vanishing of the original cohomology modules. This is the standard proof
of [Stacks Lemma35.3.6, tag023M](https://stacks.math.columbia.edu/tag/023M),
whose statement and proof were read; [tag03FJ](https://stacks.math.columbia.edu/tag/03FJ)
identifies the same affine-cover cocycle complex. The argument above
reproduces the needed algebra rather than appealing to a classification.

Equations (3) imply that g-r has equal rational pullbacks. At the generic
point the relation ring is L tensor_E L, so faithful flat equalizer
descent for the field extension gives g-r in E. Thus g belongs to R+E.
The reverse inclusion in (1) follows just by differentiation.

AFFINENESS is used to solve (2) globally by a single r in R. SURJECTIVITY
supplies faithful flatness, and ETALE smoothness supplies the regular pole
test on the relation. No finite-map trace averaging is used. In positive
characteristic the pole test fails, even for finite etale covers: on the
Artin--Schreier cover k[y]/k[y^p-y], g=y^(-p) has dg=0 but is not in
k[y]+k(y^p-y), because its pole at y=0 cannot be the principal part of a
base function without matching poles at the other roots of y^p-y.

## 2. Complete exactness criterion for every permitted rational slice

Retain MR0/MV0's actual assumptions and notation

    R=C[x,y], u=h(y), deg(h)=m>=2,
    E=C(u,s), s=x/h'(y)+beta(y), beta in C(y),
    C=R intersect E smooth affine, X -> V=Spec(C) SURJECTIVE ETALE,
    Pi=Tr_(L/E)/m, gamma=Pi(beta), beta0=beta-gamma,
    Omega=ds wedge du.

Then the following are equivalent, in EVERY degree:

    (a) Omega has a global regular Kahler primitive on V;
    (b) beta0*du=dG for some rational G in C(y);
    (c) beta0(y)*h'(y)*dy has zero residue at every finite point. (4)

Trace G can always be normalized to zero. Indeed dPi(G)=Pi(beta0)du=0,
so Pi(G) is constant. Rational partial fractions prove (b)<->(c).
The infinity residue is then zero as well. Condition (c) is a criterion
for the supplied rational beta, not an assertion that all such beta pass.

For (a)->(b), the actual source primitive is x*dy+dF0 with F0 in R by
polynomial Poincare. Taking normalized field trace, as in MV0 and KO2,
gives G=F0-Pi(F0) in C(y) with dG=beta0 du and trace zero.

For (b)->(a), observe the identity

    dG=(s-gamma)*du-x*dy.                                    (5)

Here x*dy is a regular source form and (s-gamma)du is a rational E-form.
Apply (1): there is F0 in R with F0-G=H in E. In particular the formerly
separate POLYNOMIAL SOURCE REPRESENTATIVE EXISTS AUTOMATICALLY. The form

    theta=(s-gamma)*du+dH

pulls back to x*dy+dF0. Regularity descends along the surjective etale
map, at EVERY point of V, and dtheta=Omega. This proves (4).

One can see the signs directly on the relation. Write alpha_i=x_i*dy_i.
Then d(G2-G1)=alpha1-alpha2. Splitting its regular cocycle by
G2-G1=F02-F01 makes alpha_i+dF0_i equal on both sheets. Thus these
regular forms glue; neither a meromorphic extension assumption nor an
unproved polynomial degree bound is left in the argument.

## 3. Consequences and the remaining HC4 frontier

For an actual intermediate Keller factor the supplied target pair already
forces exact Omega. Thus ALL its permitted beta satisfy (4). Conversely,
under the actual smooth affine surjective etale hypotheses, (4) accounts
for the WHOLE obstruction to a regular primitive; there is no further
polynomial-representative test or increasing-degree F0 search to perform.
PF138 established the analogous completeness on whole translation atlases;
the relation argument here requires no such atlas and permits non-Galois
relative covers, including the first 2+3 / S5 model.

After writing beta=gamma+G'/h', the field slice becomes

    s-gamma=(x+G'(y))/h'(y).                                 (6)

G and G' can have poles. Equation (6) is a RATIONAL symplectic shear,
not a polynomial source automorphism. It does not remove the permitted
rational slice or change the actual source image to a plane. Base slice
shifts beta->beta+q(h) leave beta0 unchanged.

MV0's explicit mixed graph cover and all-degree multi-root graph exclusion
remain valid: that family fails (b) itself. Its source-image realization
does not have a plane target. PF130's different exact mixed control and
PF138's exact moving-label surfaces are not excluded. No general mixed
source image or actual target pair is constructed here.

The next missing input must use the actual target COFRAME or the FULL
finite B-order / trace-dual multiplication, beyond a regular primitive.
Equation (1) controls additive polar parts of a rational function once
its differential splits in the stated way. It does not show an arbitrary
source polynomial integral over B, stabilize the finite trace lattice,
or imply Tr(Lambda_L*Lambda_L) subset B. MU0's conditional marked quadratic
trace pole is retained; its marker cannot be dropped. No Mathieu ordinary
marked-period comparison or full target-prime regularity follows.

## 4. Correction and assurance boundary

The MV0 iff with BOTH G and F0 remains true, but its descriptions of F0 as
an independent remaining obstruction were too strong. The prior final
response made the same overstatement. MW0 discharges that condition by
global descent under exactly those hypotheses. Frozen mathematical files
are unchanged; the live ledger, methods, work items and navigation record
this correction and the completed residue-only criterion. No registered
mathematical checker is replayed. Supporting new controls check descent
signs, rational differentials and characteristic boundaries; they are not
independent verification of the unbounded proof.
