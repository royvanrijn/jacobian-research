# Computing a given finite local image by partial saturation

Given a certified finite ambient quotient for the elliptic slice,
its exact rational image can be determined without finding a full
Mordell--Weil basis. Only primes dividing that finite quotient's order
need to be saturated. Every division test is finite, and every successful
adjoining step decreases a positive integer index.

## 1. The supplied quotient and the initial subgroup

Retain the elliptic curve `E`, the group `Gamma=E(Q)` and the
four-class map `alpha:Gamma->A` from
[the elliptic arithmetic note](ELLIPTIC_EXCEPTION_ARITHMETIC.md).
That note proves rank two and the full rational torsion group
`Tcal={O,T_0}`. The finite certificate
[ETC3a](ELLIPTIC_EXCEPTION_COMPLETE_THREE_ADIC_CLASSIFICATION.md)
supplies two independent points `P_1,P_2`. Thus

\[
 \Gamma_0=\langle P_1,P_2\rangle\oplus\mathcal T
                  \subseteq\Gamma                  \tag{1.1}
\]

has finite index. Its displayed points are a basis of the free
quotient of `Gamma_0`; they are not asserted to be a basis of
`Gamma`.

Fix a finite set of places `S` and put

\[
 G_S=A\times\prod_{p\in S_f}E(\mathbb Q_p),\qquad
 H_S=\overline{\{(\alpha(P),(P)_p):P\in\Gamma\}}.
\]

The input to the algorithm is a **certified** open product subgroup
and its finite ambient quotient

\[
 W=\{1_A\}\times\prod_{p\in S_f}E_{d_p}(p)
       \subseteq H_S,\qquad F=G_S/W.               \tag{1.2}
\]

Here the containment has been proved, as in
[EQA1](ELLIPTIC_EXCEPTION_QUANTITATIVE_FINITE_LOCAL_ASYMPTOTICS.md).
The input includes an exact finite-group description of `F`, the
local coset interpretation of its coordinates, and an effective exact
reduction map

\[
 \varphi:\Gamma\longrightarrow F
\]

on any explicitly given rational point. This is the scope of the
given-quotient hypothesis. The construction of this input for an
arbitrary new set of primes is not part of the result.

## 2. Every rational division test is finite

**Lemma FQS1 (decidable rational division).** For a positive prime
`ell` and a nonidentity point `R in E(Q)`, the complete set of
`Q in E(Q)` with `[ell]Q=R` is finite and can be determined by
exact rational arithmetic and finite polynomial factorization.

**Proof.** The chord-and-tangent group law computes multiplication
as rational functions on the nonsingular short Weierstrass curve.
The function `X([ell]Q)` is invariant under `Q -> -Q`. The fixed
subfield of `Q(X,y)` under `y -> -y` is `Q(X)`, so it has the form

\[
 X([\ell]Q)=R_\ell(X(Q)),\qquad R_\ell\in\mathbb Q(X).
\]

This rational function is nonconstant. Indeed, at the identity the
formal parameter satisfies `[ell]t=ell t+O(t^2)`, and
`X=t^(-2)+O(1)`. Thus `X([ell]Q)` has a pole there. This also
follows from the formal-depth calculation of EQA1. The precise
retained local input is Milne, *Elliptic Curves* (2006), Chapter II,
Theorem 4.1 and Lemma 5.9, printed pp.62--63 and 67:
[primary book](https://www.jmilne.org/math/Books/ectext6.pdf),
[retained copy](sources/milne-elliptic-curves-2006.pdf).

Write the computed rational function in reduced form `U(X)/V(X)`.
Every rational preimage has its `X` coordinate among the rational
roots of the nonzero polynomial

\[
 U(X)-X(R)V(X).                                     \tag{2.1}
\]

After clearing denominators these roots can be found by integer
factorization and the rational-root theorem, including a possible
zero root. For each root with `V(X)!=0`, test whether the curve's
right side is a rational square, take its at most two rational
square roots, and check `[ell]Q=R` by the exact group law. No
identity preimage is possible since `R!=O`. Reduction of `U/V`
ensures that a zero of `V` maps to the identity rather than to
the finite point `R`. Thus the procedure includes every rational
preimage and rejects any extraneous coordinate candidate.
\(\square\)

This proof does not require an imported numerical bound for the
degree or number of division points. Nonconstancy and the explicit
rational function suffice for a finite decision procedure.

## 3. Saturating just the primes visible in the quotient

Suppose an intermediate subgroup has exact group data

\[
 \Gamma'=\mathbb Z B_1\oplus\mathbb Z B_2\oplus\mathcal T.
                                                               \tag{3.1}
\]

It contains the full rational torsion and has finite index in
`Gamma`. Call it `ell`-saturated if every `Q in Gamma` with
`[ell]Q in Gamma'` already belongs to `Gamma'`.

**Theorem FQS2 (termination and the exact finite image).** Starting
from (1.1), the following procedure terminates:

1. Factor `|F|`.
2. For each prime `ell` dividing `|F|`, enumerate the nonzero
   classes in `Gamma'/ell Gamma'` using (3.1).
3. Apply FQS1 to a representative of every such class. If a
   rational division point `Q` is found, replace `Gamma'` by
   `Gamma'+ZQ`, update its exact basis data, and restart the tests
   at that prime. If none is found, move to the next prime.

At termination,

\[
 \gcd([\Gamma:\Gamma'],|F|)=1,\qquad
                 \boxed{\varphi(\Gamma')=\varphi(\Gamma).}      \tag{3.2}
\]

The unknown initial index need not be computed or bounded in advance.
No saturation is required at primes not dividing `|F|`.

**Why the finite tests certify saturation.** If a tested nonzero
representative `R` has `[ell]Q=R`, then `Q` cannot be in
`Gamma'`: otherwise `R` would be zero in `Gamma'/ell Gamma'`.
The new coset `Q+Gamma'` has exact order `ell`. Hence

\[
 [\Gamma'+\mathbb Z Q:\Gamma']=\ell.               \tag{3.3}
\]

Conversely, if `Q notin Gamma'` but `[ell]Q in Gamma'`, write
`[ell]Q=R+[ell]B` for a chosen representative `R` and `B in Gamma'`.
Then `Q-B` is a rational division point of `R`. The representative
cannot be the identity: in that case `Q-B` would be rational
`ell`-torsion, all of which is already in `Gamma'`. Thus the
nonzero representatives suffice precisely because the full rational
torsion was included at the start.

Equivalently, saturation says that the finite abelian quotient
`Gamma/Gamma'` has no element of order `ell`, hence that its order
is prime to `ell`. Each successful step divides its positive integer
order by `ell`, by (3.3), so there can be only finitely many such
steps. Each individual round consists of finitely many terminating
division tests. Later adjoining at a different prime only divides
the index by that other prime and preserves the saturation already
proved. This establishes termination and the first assertion of
(3.2). If `|F|=1`, the prime list is empty and this assertion is
immediate.

**Maintaining exact group data.** Representatives can be written

\[
 R=aB_1+bB_2+\tau,\qquad 0\le a,b<\ell,
       \quad \tau\in\mathcal T/\ell\mathcal T.
\]

For this curve there are `ell^2-1` nonzero classes for odd `ell`
and seven for `ell=2`. For any successful adjoining step at least
one of `a,b` is nonzero: otherwise `[ell]Q` is torsion, making
`Q` rational torsion as well, contrary to (3.3).

In the old free coordinates, the new generator has vector
`(a/ell,b/ell)`. The new free lattice is therefore

\[
 \mathbb Z^2+\mathbb Z(a/\ell,b/\ell)
       =\ell^{-1}\operatorname{im}_{\mathbb Z}
          \begin{pmatrix}\ell&0&a\\0&\ell&b\end{pmatrix}.
                                                               \tag{3.4}
\]

Integer Hermite or Smith reduction of this matrix provides a basis
and its explicit integer combinations of the three columns. Apply
the same combinations to `B_1,B_2,Q` to obtain two actual rational
points forming a basis of the new free quotient. Together with
the unchanged full torsion group they generate the entire new
subgroup: any remaining difference has zero free coordinate and
is already in `Tcal`. This maintains (3.1) after every step using
only known integer relations and rational point addition.

**Why the image is exact.** Put `C=varphi(Gamma)` and
`C'=varphi(Gamma')`. The quotient `C/C'` is an image of
`Gamma/Gamma'`, so its order divides `[Gamma:Gamma']`. It also
divides `|F|`. Coprimality in (3.2) therefore makes it trivial.
\(\square\)

## 4. Computing the Haar coefficient from that image

**Corollary FQS3 (exact coefficient from the supplied quotient).**
Under the given-quotient hypotheses of section 1, the procedure
determines `H_S/W` exactly. With the local weight calculation of
EQA2 it computes each rational finite-place coefficient `mu_(S,a)`
exactly. Consequently it decides whether any admissible full field
target of class `a` passes the specified finite set `S`.

**Proof.** Enumerate the finite subgroup generated in `F` by the
final explicit basis and torsion images. By FQS2 this is
`C=varphi(Gamma)`. The image of the closure in a finite discrete
group is the image of its dense rational subgroup, so `H_S/W=C`.
The supplied containment `W subset H_S` is essential: it gives

\[
 H_S=\pi^{-1}(C),\qquad \pi:G_S\longrightarrow F.
\]

Finite-group enumeration expresses every `f in C` as a combination
of the known generators and hence gives a rational representative
`P_f`. Let `a_f=alpha(P_f)` and let `r_(f,p)` be the soluble
conditional Haar measure in `P_f+E_(d_p)(p)`. EQA2 calculates these
rational numbers by finite local arithmetic and rational geometric
series. Thus

\[
 \boxed{\mu_{S,a}=\frac1{|C|}
          \sum_{f\in C:\ a_f=a}\ \prod_{p\in S_f}r_{f,p}}       \tag{4.1}
\]

is an exact finite calculation. It retains the actual finite joint
image and all correlations. Comparing its rational value to zero
and applying the finite-local dichotomy EQA2a proves the final
assertion. \(\square\)

The input (1.2), including its reduction operation and local coset
data, remains a required certificate. This result does not construct
that input for an arbitrary `S`, certify a full basis of `E(Q)`,
bound the running time, or address the all-primes limit. It identifies
the finite saturation obligations needed after that local quotient
has been supplied.

The proof is the rational-division decision, exact lattice basis
update, finite-index descent and finite-group argument above.
Its rank, torsion, initial independence, local-quotient and local-weight
inputs are the retained sources linked in sections 1--4. No new
finite experiment is used as an unrestricted proof, and no formal
verification or external review is claimed.
