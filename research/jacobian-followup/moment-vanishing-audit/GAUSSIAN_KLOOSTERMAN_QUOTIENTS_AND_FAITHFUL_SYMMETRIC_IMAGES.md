# Gaussian Kloosterman quotients and faithful symmetric-power images

16 November2026. Internal all-degree construction and exact finite controls.
Independent replay, formal verification and external review remain false;
priority remains unassessed. This is a general actual-connection realization
bridge, not a regular-transcendental Gaussian example or a complete safe
normal-form classification. See the [constructor](gaussian_kloosterman_transfer.py),
[checker](verify_gaussian_kloosterman_transfer.py),
[receipt](gaussian_kloosterman_verification.json) and
[source audit](gaussian_kloosterman_source_audit.json).

## 1. Full marked image from a rational pole packet

Use circular covariance E[ZW]=1, real covariance E[T^2]=1, and the
retained pure phase

    B=W(1+Z)-(1+3Z/2+Z^2/2)T^2.

Its FULL original polynomial-Weyl image is K[t,(1-t)^-1], with
Phi_B[Z]=t/(1-t), as sealed in the
[finite-multiplier proof](FINITE_MULTIPLIER_MIXED_DECISION.md).
Add k independent real Gaussians Y_i, k>=2, and put

    C_k=product_i Y_i^2, P_k=C_k B,
    a_m=((2m-1)!!)^k, h_lambda=sum_(m>=0) a_m*(lambda*t)^m,
    theta=t*D_t, A_lambda=1-lambda*t*(2theta+1)^k.       (1)

Scaling B by lambda realizes h_lambda. The seed is1, and the ORIGINAL
linear marker Z has period h_lambda-1. No factorial is missing: a_m
multiplies the ORDINARY coefficient of a marked B exponential period.
For EVERY B whose full image is a finite rational pole localization
with distinct nonzero pole scalars S, the same construction gives

    M_(C_k B)=sum_(lambda in S) D_t*h_lambda.           (2)

For S empty this is K[t]. Polynomial parts stay polynomial. For every
pole order p and even spectator marker product Y_i^(2r_i), the exact
Euler transfer is

    E[product Y_i^(2r_i)*(1-lambda*t*C_k)^(-p)]
      = [product_(j=1..p-1)(theta+j)/(p-1)!
         *product_i product_(j=0..r_i-1)(2theta+2j+1)] h_lambda.

Odd spectator markers vanish. This proves the reverse inclusion in(2)
for ALL original amplitudes, not just the unit or one selected period.
The recurrence a_(m+1)=(2m+1)^k*a_m proves A_lambda*h_lambda=1.

The Newton weight k*ord_t(coefficient)-theta_degree gives initial
symbol 1-lambda*2^k*T*xi^k, prime in K[T,T^-1,xi]. Its xi degree is
the actual operator order. The same leading-symbol argument as in the
[sealed square transfer](SQUARED_SPECTATOR_INVARIANCE_AND_UNBOUNDED_GAUSSIAN_IMAGES.md)
proves A_lambda irreducible over K((t))<D_t>, and distinct lambda give
pairwise nonisomorphic modules. Each h_lambda is nonrational: its
coefficients have exact Gevrey order k. Consequently

    rank M=1+k*|S|,
    K(t)M/K(t)1=direct_sum_lambda W_lambda,
    W_lambda=K(t)<D_t>/(A_lambda).                       (3)

All W_lambda have pure slope1/k at0, and total irregularity1 each;
the unit adds slope0. At infinity they are regular singular. No
actual-image kernel or guessed minimal annihilator is used.

## 2. Exact Kloosterman identification and classical Galois import

Work over C, and first retain ONE nonzero pole lambda. Set

    z=(-1)^k/(2^k*lambda*t), y=z^(1/2)*u.

Inversion gives theta_t=-theta_z. After the displayed Kummer gauge,
the scalar identity is EXACT:

    -z * z^(-1/2) A_lambda z^(1/2)
       =theta_z^k-z.                                  (4)

The right side is the rank-k Kloosterman connection. The primary
published source is [Fresan--Sabbah--Yu, Duke Math.J.171(2022),1649-1747,
final arXivv5](https://arxiv.org/html/1810.06454v5): Proposition2.4 and
its complete proof identify this operator, irreducibility and local
types; the proof of Proposition2.7 explicitly restates Katz's differential
Galois theorem. Its group is SL_k for odd k and Sp_k for even k.
The latter notation means a symplectic group in dimension k.
The original Katz proof is imported through this explicit primary
restatement, not claimed independently verified here. We import no
motivic automorphy, finite-field monodromy or Hodge-number result.

The half-Kummer twist leaves the CONNECTED group unchanged. The full
quotient group H_k is Sp_k for even k; for odd k>=3 it is SL_k times
central mu_2, acting by scalar sign. Indeed the untwisted determinant
is trivial, while the twisted determinant is z^(k/2). The connected
Kloosterman Picard-Vessiot group has no nontrivial finite quotient, so
its field is disjoint from the nontrivial quadratic Kummer extension.
For even k the scalar -I already belongs to Sp_k; for odd k it lies
outside SL_k. At time infinity the quotient has monodromy -J_k,
where J_k is one unipotent Jordan block. These statements concern the
actual quotient(3), not a larger exponential source.

For a multi-pole packet we claim (2),(3), but do NOT infer an independent
product of Galois groups: odd-k opposite poles give dual modules and
can have shared Picard-Vessiot fields.

## 3. Full affine Galois group and the nonsplit unit

For one pole the exact generic image V has

    0 -> K(t)1 -> V -> W_lambda -> 0.                  (5)

The extension is nonsplit. Every nonzero Laurent series f has the same
initial t-term as A_lambda*f, so A_lambda*f=0 has no nonzero Laurent
solution. A section of(5), composed with the ACTUAL scalar-series
embedding of V, would give a nonzero map W_lambda->C((t)), impossible.
The unit is essential by the same argument and irreducibility.
The global image has rational core C[t] and has exactly one normalized
formal period map: A_lambda is bijective on C[[t]].

The differential Galois group G of V fixes the unit and surjects onto
H_k. Its kernel is a closed subgroup U of the additive k-space of
extension rows, stable under H_k. Such a subgroup is a vector subspace
in characteristic zero. Irreducibility of the H_k representation forces
U=0 or U=Ga^k. If U=0 the group is reductive and complete reducibility
would split(5). Therefore U=Ga^k, and filling every translation fiber
over H_k gives the FULL group

    G=Ga^k semidirect H_k.                             (6)

Module and solution conventions can dualize the standard representation;
the displayed classical affine group is unchanged up to isomorphism.
In particular G^0 is NONSOLVABLE for every k>=2. Gaussian pure phases
therefore realize arbitrarily large classical reductive groups and their
full affine unit extensions. A proposed exhaustive description by finite
algebraic-exponential induction and Liouvillian quadratures is false.
The Kloosterman classification is classical credit, not a new theorem
about those differential equations.

The Liouvillian implication uses the classical Kolchin criterion, as
restated by [Singer0712.4124v1](https://arxiv.org/html/0712.4124v1),
Section5 Definition6 and Theorem7 (printed1.5.6--1.5.7), together with
the standard fixed-field consequence of Section3 Corollary10(3): a
connected Picard-Vessiot group leaves no nontrivial finite algebraic
subextension. We do not import Singer's general Galois algorithm or
family constructibility results. An arbitrary annihilating equation
can have extra Liouvillian solutions; the full ACTUAL module(3),(5)
is essential for the conclusion here.

## 4. Independent sums realize FAITHFUL symmetric powers

Take s>=1 independent copies of the SAME phase P_k and add them. Every
original marked period is a product of s marked P_k periods, with unit
markers filling unused factors. Thus multiplication gives a surjection
from Sym^s_(K[t]) M_(P_k) to the FULL actual image of this new Gaussian
phase. Here it is injective, not an arbitrary tensor-product assumption.

To prove this, work in a Picard-Vessiot field of V. The k coordinates
h,theta*h,...,theta^(k-1)*h of the actual normalized period form a
particular inhomogeneous solution vector. Over the Picard-Vessiot field
of W_lambda, the full translation group Ga^k in(6) translates that
vector by an invertible homogeneous fundamental matrix times arbitrary
constant k-vectors. Therefore these k functions are algebraically
independent over that field, and hence over K(t). A polynomial relation
would vanish on an entire affine k-space, impossible. All monomials
of total degree<=s are consequently independent over K(t).

It follows that the generic multiplication map is injective. The global
module M is torsion-free over the PID K[t], hence flat. Its symmetric
power is a direct summand of its tensor power by characteristic-zero
averaging, hence torsion-free as well. A generic kernel would detect
every global kernel. Thus the FULL global image equality is

    M_(sum s copies P_k)=Sym^s_(K[t]) M_(P_k),
    rank_actual=binomial(k+s,s).                        (7)

Its Galois group is still(6): the invariant unit makes Sym^s V faithful,
since the subspace 1^(s-1)*V recovers the original representation.
The filtration by polynomial degree has graded pieces Sym^j W_lambda,
j=0,...,s; these are irreducible for SL/Sp as restated in the same
Proposition2.7. No finite rank table proves this all-s assertion.

## 5. Identical ordinary moment geometry, different connection complexity

Use the polynomial coefficient family

    B_e=W(1+eZ)-e(1+3eZ/2+e^2 Z^2/2)T^2.

For e!=0, the Gaussian change Z'=eZ,W'=W/e gives B_e=eB. For e=0 the
phase is W. Its square-product extensions and their independent s-copy
sums are pure over QQ[e] and every Q-algebra base change. All moments
of an original monomial marker have the single e-exponent
m+total(W degree)-total(Z degree), when nonzero; independence preserves
this scaling in sums. Negative exponents force zero polynomial moments.
Conversely one-copy markers Z^d, or Z^m when d>=m, attain the lower
exponents, with nonzero rational spectator factors. Hence the ENTIRE
ordinary profile is independent of k and s:

    V_(d,0)=(1), V_(0,m)=(0) for m>0,
    V_(d,m)=(e^max(m-d,0)), d>=1,m>0.                  (8)

This is equality of ordinary ideals, including all nilpotent and embedded
base-change structure, not a radical statement. It gives the sealed
sharp forbidden ideals for every REAL C,E>=0:

    J_(C,E)=(1) if C<1;
    J_(C,E)=(e^floor(C+E)) if C>=1.                    (9)

Over QQ[e]/e^L the full sharp region is C>=1,E>=0,C+E>=L. At e=0 the
rank is1, while every e!=0 has rank binomial(k+s,s) and group(6).
Thus the complete ordinary moment/cutoff geometry cannot recover even
the connected reductive group, its nonsplit unit, or symmetric-power rank.

A sharper SAME-RANK comparison uses two five-variable degree-eight pure
phases: P_old=Y^4 B (with one unused Gaussian) and
P_new=4Y_1^2 Y_2^2 B. Both have actual rank3, exact Gevrey order2,
identical ordinary profiles, and the same ramified formal type at0.
For the quotient operators

    A_old=1-t(4theta+1)(4theta+3),
    A_new=1-4t(2theta+1)^2,

put w=1/(2sqrt(t)). The old homogeneous equation has EXACT solutions
w^(1/2)exp(+-w); the new one, after y=w*u, is
u''+w^-1*u'-u=0, with formal types w^(1/2)exp(+-w) times formal units
for y. Both full unit extensions
split formally at0 because their inhomogeneous equations have unique
formal series solutions. Yet the old group has solvable identity
component (algebraic exponentials and variation-of-constants quadratures),
whereas the new group is Ga^2 semidirect SL_2. We do NOT assert that
their formal types at infinity agree, or compute their Stokes matrices.

This rules out classifying actual Gaussian connections from ordinary
cutoffs, rank, Gevrey order and ramified origin formal type alone.

## 6. Completion and source boundaries

The constructor implements arbitrary k and faithful symmetric connection
matrices; a finite checker retains original Wick moments, full momentum
marker transfers, recurrence/connection/gauge identities, the rank-three
comparison and polynomial e-scaling. Classical Galois and translation-
independence arguments remain written proofs, not computed Galois groups
or independent/formal/external verification.

This opens a genuinely non-solvable core beyond the earlier filtered
algebraic-exponential architectures. It does NOT produce a regular
transcendental Gaussian image: the positive slope at0 survives. The next
structural gate is a transfer realizing regular period connections inside
the ORIGINAL full Gaussian image, or a valid general obstruction. An
inverse Fourier/Borel transform or a toric period alone does not supply
that transfer. No classification of every Gaussian phase is claimed.
All earlier Oaku operators,torsion witnesses,specialization thresholds,
finite-primary backend and terminal UNKNOWN jobs remain sealed.
Sabbah/Maisonobe retain convergent rational-slope credit; nonconvergent
formal/all-ring priority remains separately unassessed. Direct Mathieu/SIC
work stays separately owned; the broad general goal remains active.
