# Intrinsic marked Stein action and full rank-one Gaussian lattices

16 November2026. Internal all-degree proofs without a weight-frame premise.
Independent replay,formal verification and external review remain false;
priority unassessed. The [constructor](gaussian_stein_rank_one.py),
[checker](verify_gaussian_stein_rank_one.py),
[certificate](gaussian_stein_rank_one_certificate.json) and
[receipt](gaussian_stein_rank_one_verification.json) retain finite controls.
Older filtered,marker-filtration,Oaku and actual-image packets stay sealed.

## 1. A full ORIGINAL-marker identity over every Q-algebra

Let R be any commutative Q-algebra,C a fixed nonsingular rational symmetric
Gaussian covariance,and P in R[X]. Put Delta_C=sum C_ij partial_i partial_j
and N=X.grad-Delta_C. Its inverse on nonconstant Hermite components is
the finite polynomial operator

    N^-1=exp(-Delta_C/2) Euler^-1 exp(Delta_C/2),

where Euler^-1 kills the constant and divides each positive ordinary
homogeneous degree d by d. Define on ORIGINAL polynomial amplitudes

    Gamma_P(A)=grad(P)^T C grad(N^-1 A).                 (1)

There are no formal coefficient queries or infinite operations in(1).
Let Phi[A]=E[A exp(tP)],f=Phi[1]. Gaussian divergence of
C grad(N^-1 A) is exactly A-E[A]. Integration by parts gives

    Phi[A]=E[A] f+t Phi[Gamma_P(A)].                    (2)

This holds coefficientwise over EVERY R,with no positivity,convergence,
purity or safety assumption. It is the polynomial Malliavin--Stein
identity,credited to [Azmoodeh--Gasbarra--Gaunt1912.04605v4,
Section2,Definition2.7/Lemma2.3](https://arxiv.org/html/1912.04605v4).
The algebraic divergence computation supplies the all-ring extension.

Let V=Phi(K[X]) over a characteristic-zero field K,and let K0 be the
kernel of this ORIGINAL spatial-amplitude map. Equation(2) shows Gamma_P
preserves K0,so it acts on the exact quotient K[X]/K0=V. On actual periods
its action is exactly

    B_f(h)=(h-h(0)f)/t.                                (3)

This does NOT mean that Gamma preserves an arbitrary ambient annihilator
kernel or a guessed cyclic scalar presentation. The period kernel is the
one in this statement. No kernel computation is needed to prove(2)-(3).

The full polynomial-Weyl image M is K[t]V;it is also D_t-stable because
partial_t Phi[A]=Phi[PA]. Usually V is NOT already a K[t]-module. For
P=X^2/2,V is spanned by(1-t)^(-j-1/2),j>=0,whereas t*f is not in V.
All later full-image statements explicitly use M,not that false shortcut.

The origin-saturation theorem ensures(3) also preserves M:its numerator
belongs to M and has zero initial coefficient,so division by t remains
in M. Moreover B_f(t h)=h and ker B_f=Kf. Normalization gives

    B_f(h)/f=B_1(h/f).                                 (4)

## 2. An intrinsic criterion for ACTUAL rank one

For a formal series q over a field,the backward-shift orbit
span{B_1^j q:j>=0} is finite dimensional iff q is RATIONAL. This is the
elementary linear-recurrence characterization of rational generating
functions. Indeed a polynomial relation in B_1 yields a constant-coefficient
recurrence for its coefficients;conversely a rational denominator yields
such a recurrence after its finite initial polynomial part. The orbit is
eventually zero iff q is a POLYNOMIAL.

Consequently,in EVERY dimension and phase degree,

    rank_actual M=1
      iff Gamma_P is LOCALLY FINITE on K[X]/K0,         (5)
    Gamma_P locally NILPOTENT on K[X]/K0
      iff M=f K[t].                                   (6)

Local finiteness/nilpotence means every vector has its own finite orbit
dimension/nilpotence index;there is no one global exponent. In(5),all
ratios Phi[A]/f are rational iff the generic span is one dimensional.
In(6),all ratios are polynomial iff K[t]V=fK[t],because f itself occurs.
These are criteria on the FULL actual marked quotient,not bounded tests
on selected amplitudes or iterates. Higher-rank Kloosterman/Airy images
therefore necessarily have some infinite-dimensional marked Stein orbit.

In(6),D-stability forces r=f'/f to be a polynomial. Hence

    f=exp(psi(t)), psi in tK[t], M=exp(psi)K[t].        (7)

Conversely that full image makes the action locally nilpotent. This
classifies the actual images with locally nilpotent marked Stein action
without claiming an exhaustive PHASE normal form. It includes the known
Gaussian realization P=Z+W^d,whose seed is exp(t^(d+1)). That phase is
UNSAFE despite locally nilpotent Stein action. Thus local nilpotence
alone is not a general safety test. Adding polynomiality of f makes
(7) force psi=0 and gives pure safety.

Under purity f=1,(3) is the usual shift. Then

    P safe iff B_1 locally nilpotent on M,
    M intersect K(t) = locally finite part of B_1 on M. (8)

This identifies the old rational core intrinsically. For a nonpure seed,
the locally finite part is instead M intersect fK(t),which need not be
the ordinary rational core. The two must not be interchanged.

## 3. Full rank-one GLOBAL lattices, including nonpure seeds

Assume the actual image has rank one,and retain faithful finite normalized
generators h_i/f in K(t),including1. Then r=f'/f belongs to K(t),is regular
at the normalization origin,and M=fN for a K[t]-submodule N of K(t)
containing K[t],stable under partial_t+r. It is generated by the supplied
rational functions under these operations. Work over a geometric splitting
field for the following local description;the resulting divisors descend.

At a finite point a!=0,let L_a be the supremum of pole orders allowed
in N. Each L_a is a nonnegative integer or infinity. Polynomial partial
fractions and the Chinese remainder theorem isolate every point and every
lower principal part,so these bounds specify N COMPLETELY. There are only
finitely many nonzero bounds:derivatives introduce no new poles outside
the denominators of r and of the finite generators.

The exact rules are:

* At a pole of r of order>=2,L_a=infinity.
* At a simple pole with residue k NOT a positive integer,L_a=infinity.
* At a simple pole with residue k a positive integer,L_a=k or infinity.
  It is infinity exactly when an original supplied normalized generator
  has pole order greater than k.
* At a regular point of r,L_a=0 or infinity,with infinity exactly when
  a supplied normalized generator has a pole there.

To prove the rules,take a principal part c(t-a)^(-j). The next derivative
under partial_t+r has leading coefficient(k-j)c at order j+1 when r
has a simple pole of residue k. A higher-order pole of r produces an
even higher nonzero pole. Starting with1 generates the mandatory ladder
up to k when k is a positive integer;it never stops otherwise. A generator
above k restarts an unbounded ladder. Conversely a finite bound k is
preserved when that leading coefficient cancels. Rational local partial
fractions prove both inclusions globally;no unknown source-rank equality
is inferred from a scalar operator.

Let S be the infinity-bound points. Put

    d=product_(a in S)(1-t/a),
    p=product_(a notin S, L_a>0)(1-t/a)^L_a.

Then p(0)=d(0)=1,p and d are coprime,and

    N=p^-1 K[t,d^-1], M=(f/p)K[t,d^-1].                (9)

This is the full GLOBAL image,not just its generic connection. The finite
apparent-zero polynomial p is essential for nonpure seeds. The old pure
rational-localization theorem is exactly p=1,f=1. The routine implements
(9) over Q by irreducible denominator factors and exact residue arithmetic;
it consumes faithful rank-one data and does NOT certify Gaussian realization
or actual rank from finite jets.

Safety and purity now have exact independent tests:

    safe iff d=1 AND r=p'/p,
    pure iff r=0.                                     (10)

Indeed M=K[t] iff the normalized gauge f/p is constant1 and no actual
pole is allowed. If M=K[t],multiplication by a nonconstant polynomial
gauge would fail D-stability,while a nonconstant rational denominator
would make poles;these are excluded by the local rules. The normalization
at zero fixes the constant. A rank-one connection is regular on P1 iff
r has at most simple finite poles and r=O(1/t) at infinity. Integer
gauge p does not alter that criterion.

No general quasi-unipotence theorem for arbitrary exponential Gaussian
direct images is imported here. An ABSTRACT regular rank-one connection
such as r=1/(1+t^2) has a transcendental solution exp(arctan t). Its
Gaussian realizability is NOT established. The rigid-connection theorem
requiring quasi-unipotence cannot be used to infer it for arbitrary input.

## 4. The exact remaining nonpure-safe obstruction

For a safe phase,M=K[t] and its unit period f is polynomial. If m=deg f>0,
the marked Stein action has the canonical invariant decomposition

    K[t]= fK[t] direct_sum K[t]_(degree<m).              (11)

On the first summand it is the ordinary nilpotent shift after division
by f. On the second it is an INVERTIBLE m-dimensional companion action,
with characteristic polynomial

    xi^m f(1/xi).                                     (12)

Its nonzero eigenvalues are reciprocal roots of f,and the Jordan length
at one root is its multiplicity. By contrast,an actual normalized rational
pole in S gives an UNBOUNDED Stein Jordan tower at the reciprocal point.
Formula(9) distinguishes these finite and infinite chains exactly.

Thus the general field safety=>purity problem is reduced to whether an
ORIGINAL Gaussian safe image can have the nonzero finite spectral summand
in(11). The abstract model K[t] with seed(1-t)^2 has such a summand and
is consistent as a D-module with a chosen seed;it is NOT a Gaussian
counterexample. The actual polynomial-seed compensator with that seed
instead has a pole generator of normalized order3,so S contains1 and
its Stein chain is unbounded:it is unsafe. The previous filtered theorem
excludes the finite-summand possibility in its whole supplied architecture.

This is an intrinsic all-degree obstruction and a full rank-one image
classification. It does not prove the remaining Gaussian realization
exclusion,discard the original marker filtration,or restore any universal
algebraic-exponential claim in higher rank.

## 5. Source audit: full Stein ideals and indexing boundaries

[Azmoodeh--Gasbarra--Gaunt2212.07321v2,29April2025](https://arxiv.org/html/2212.07321v2),
Theorem3.1(b) and its proof,gives the Gaussian principal right Stein ideal.
Theorem5.1 and its proof gives common nonzero operators for finite families
with holonomic densities/characteristic functions. The latter is a direct
warning that ONE scalar Stein operator need not identify an actual law
or image. We use the full original-marker identity(2) instead. This is
an arXiv source;no unverified final journal citation is invented.

Two indexing boundaries in the selected v2 displays are retained. The
family H_k partial^m-H_(k+m) contains zero at m=0,so a basis version must
omit those entries. The displayed telescoping RHS uses an index one too
high;the polynomial-Weyl identity actually used is

    H_k partial^m-H_(k+m)
      =(partial-x) sum_(r=0)^(m-1) H_(k+r) partial^(m-1-r), m>=1.  (13)

It follows at once from the source's preceding correct Lemma3.2 and is
checked coefficientwise here. At k=0,m=1,the printed RHS with H_1 already
differs from partial-x. These presentation corrections do not refute the
principal-right-ideal result. We import no probabilistic characterisation
theorem into arbitrary complex or nonreduced rings.

Derksen--van den Essen--Zhao1506.05192v1 Proposition4.6 already supplies
the classical necessary condition P(0)=0 for eventual pure vanishing.
Its complete proof was reread;it does not force E[P]=0. No new constant-
term theorem is claimed. General quasi-unipotence/regular-transcendental
realization and the finite-spectral-summand exclusion remain separate gates.
Sabbah/Maisonobe retain convergent rational-slope credit;formal nonconvergent
and all-ring extension priority remain separately unassessed.

The focused controls include15original Stein images,75original IBP moment
identities,an ordinary nilpotent control,nine faithful rank-one models,
three finite companion matrices,three Jordan ladders,three corrected
source identities and five rejects. Universal statements use the written
arguments. Broad goal ACTIVE;no deeper finite campaign,delegation,outreach,
terminal UNKNOWN restart or assurance promotion.
