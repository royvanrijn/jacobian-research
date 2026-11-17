# Full-source Kummer lifts and a single original marked resolvent

17 November 2026. Written all-degree realization lemma. Both FULL exceptional
Mathieu--Zhao statements remain **OPEN**. Independent replay, formal verification
and external review are not claimed. Review remains deferred.

A nonconstant rational Kummer system with its unit retained in the FULL
complement source yields ONE fixed original marker and an exact ORDINARY
Cauchy formula. No purity or fibre semisimplicity is needed. A character
of the reduced boundary 1-motive alone does not supply the required lift
to that full source; that distinction must remain explicit.

## 1. Exact input and conclusion

Let X=SL2/Gamma for either exceptional group (in fact any finite Gamma),
with normalized original Haar period h and original volume Omega. Fix a
nonconstant original regular phase f, and write over K=C(s)

    B_f=H^3_dR(X_K minus {f=s}),
    u=[Omega], period(u)=1.

Suppose an ACTUAL differential-module morphism from the rational Kummer
connection K_g to B_f is supplied, preserving the unit, where

    g in C(s)^* is NONCONSTANT,
    K_g=K*u direct_sum K*v,
    partial_s u=0, partial_s v=(g'/g)u.              (1)

This is a lift into B_f, BEFORE the compact-middle quotient. It may come,
for example, from a genuine geometric pullback of a full rank-two target
source. No equality of weights or abstract torus coordinates supplies it.

Then there are ONE fixed original regular marker Q in C[X], a nonzero
polynomial D in C[s], and a polynomial P in C[s] such that the actual
continued ordinary Cauchy germ satisfies

    C_Q(s)=h(Q/(s-f))=D(s) log g(s)+P(s).             (2)

The branch and constant of the logarithm in (2) are chosen by the actual
physical period near infinity; a constant change is absorbed into P.
The data follow by finite algebraic de Rham linear algebra from (1).
The marker does not depend on the moment order.

In particular Q has an infinite nonzero marked tail. At every finite a
with ord_a g!=0 its Cauchy continuation has a genuine integer-exponent
logarithm, even if D(a)=0. If f is pure, (2) is already an actual pure
failure with this fixed original marker; this lemma does not produce
such a pure f.

## 2. Why the Kummer logarithm survives the actual period kernel

The Kummer connection (1) has no rational flat complement to its unit.
Indeed such a complement would give a rational primitive of g'/g.
A nonconstant rational g has a zero or pole on P1, where g'/g has a
nonzero residue; derivatives of rational functions have residue zero.
Every rank-one differential submodule avoiding the unit would project
isomorphically to its trivial quotient and would be such a complement.
Thus its only rank-one differential submodule is the unit.

The actual period map from B_f to the space of continued scalar functions
is a differential-module morphism and maps u to literal 1. Its composite
with (1) therefore has zero kernel: a nonzero proper kernel would be the
unit, which is not killed. This uses the elementary rank-two Kummer fact;
it is NOT a universal essential-unit assertion for arbitrary period sources.

The image of v has derivative g'/g, so is log g plus a constant. Subtract
that constant times u. Hence the literal logarithm belongs to the actual
full marked-function module. No arbitrary solution of a differential
equation is substituted for the continued physical period.

## 3. Simple Cauchy forms span the full generic source

For every original regular R put

    alpha_R=R Omega/(s-f), C_R=h(R/(s-f)).

These classes span B_f over C(s). Here is the full-source argument.
The generic fibre F_s is smooth affine of dimension 2. Its relative
canonical module is trivialized by Omega/df, and its regular top forms
are all R*(Omega/df), with coefficients in K and original R in C[X].
Because the fibre is affine, its top de Rham cohomology is the quotient
of regular top forms by exact forms. Thus these forms span H^2_dR(F_s).
They are precisely the residues of the alpha_R, up to the fixed common
orientation. The Gysin sequence and H^4(X)=0 show that alpha_R span
B_f modulo the ambient line. That line is already in their span, by the
identity of actual forms

    Omega=s*alpha_1-alpha_f.                         (3)

This is the same residue spanning argument retained in the
[original-residue theorem](ALGEBRAIC_HAAR_LOG_COEFFICIENTS_AND_ORIGINAL_RESIDUES.md),
now used before projecting to the graph. It retains genus, compact-middle
classes and every original marker relation. No hypothesis of semisimple
local fibre monodromy is needed.

Finite dimensionality therefore writes the supplied image of v as a
finite combination sum_i a_i(s)[alpha_(R_i)], a_i in C(s). Taking ACTUAL
periods gives

    log g=sum_i a_i(s) C_(R_i)(s).                   (4)

Choose nonzero D clearing these rational denominators and put b_i=D a_i
in C[s]. Define the SINGLE original regular marker

    Q=sum_i b_i(f) R_i.                              (5)

For a polynomial b, the exact ordinary identity is

    b(s) C_R(s)=C_(b(f)R)(s)
             +h(R*(b(s)-b(f))/(s-f)).               (6)

The last term is a polynomial in s. Equations (4)--(6) prove (2).
This construction uses supplied algebraic de Rham data; it is not an
implemented universal resolution or marker algorithm, and it does not
assert that every intrinsic local class admits input (1).

## 4. The original infinity condition forces a balanced divisor

Every actual C_R is holomorphic in 1/s near infinity and is O(1/s), by
uniform convergence on the compact original Haar core. Equation (2)
therefore forces

    ord_infinity g=0.                               (7)

Otherwise log g has a nonzero log s term at infinity. Multiplication by
nonzero D and addition of P cannot remove its monodromy, whereas C_Q
has none there. Thus any Kummer lift of the form (1) necessarily has
at least two distinct finite zeros/poles. This is an ordinary-period
constraint, not a conclusion from an exponential generating function.

Scale g by a constant so g(infinity)=1 and choose log g(infinity)=0.
Factor its finite divisor as

    g(s)=product_a (s-a)^(n_a), sum_a n_a=0,
    D(s)=sum_(j=0)^d D_j s^j.

Expanding (2) at infinity gives the ALL-ORDER exact marked moments

    h(Q f^m)=-sum_(j=0)^d D_j
                  [sum_a n_a a^(m+j+1)]/(m+j+1), m>=0.  (8)

The polynomial P cancels precisely the polynomial part of D log g and
contributes no negative powers. Formula (8) uses the same fixed Q at
every order. Its tail cannot terminate: at a finite zero/pole of g,
D log g+P has a nonzero logarithmic germ. A polynomial zero D(a) merely
multiplies that logarithm by a finite power; it does not make it rational.

For g=s/(s+1), D=1, formula (8) gives (-1)^(m+1)/(m+1), the negative of
the retained uniform-interval marker tail. Reversing the Kummer coordinate
reverses this sign. This checks the ordinary normalization against the
existing actual Long/Bezout source without introducing a new phase.

## 5. The reduced boundary source does not automatically give (1)

The genuine boundary 1-motive and its corrected character valuations are
constructed from

    E_f=B_f/C_f,

where C_f is the canonical globally nontrivial compact-middle summand.
A global rational character of its semiabelian group can give a rank-two
Kummer subsource in E_f. Pulling it back to B_f gives an extension by
C_f, not an automatic Kummer embedding (1). The actual marked-function
module likewise maps to its boundary quotient; that quotient need not
embed in scalar Cauchy functions.

If C_f=0, the distinction disappears and an actual unit-preserving
Kummer subsource of E_f supplies (1). More generally the preimage
extension must be split compatibly with the ambient unit. If a proposed
rational Kummer character in E_f has nonzero order at infinity, (7)
is an obstruction to such a lift, not a contradiction in E_f: the
compact-middle quotient can create local nonsplitting at infinity.

For the regular pullbacks of Long's FULL source the lift is actually
supplied by ordinary cohomological pullback before any quotient, so (1)
holds and the earlier explicit marker is recovered in scope. A generic
character of the torus alone need not even extend to a homomorphism of
the whole semiabelian group; neither step is assumed here.

These lifting conditions belong ONLY to this literal-log realization
mechanism. They are not additional gates in the completed all-genus
intrinsic local safety criterion under purity. That criterion already
detects every nonzero reduced local component by some fixed original
marker, even when its full ordinary period is more complicated than (2).

No new exceptional pure phase, full pure-locus splitting, or ordinary
algebraicity theorem is proved. The conormal amplification route remains
unproved; the McKay/finite-model exploration supplied no replacement for
its missing identities and no primary theorem was imported from that
exploration. Capacity, automorphism, primitive-branch and completed family
transfer results and the full original goal are preserved. Both exceptional
cases remain OPEN; external review is deferred.
