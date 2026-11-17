# Hamiltonian symbol constraints on algebraic right actions

14 November 2026. **QHI1.** A normalized symbol ideal of a rank-one
Weyl module is preserved by the Hamiltonian vector field of every
symbol in its right idealizer. Consequently any remaining symbol
defect must lie where those Hamiltonian fields fail to span the
ambient tangent space. A symplectic polynomial frame in the full
right-symbol algebra therefore forces the left module to be free.

For the explicit charge-one module in
[QPF1](PROJECTIVE_FRAME_BOUNDARY.md), the full graded idealizer is
K plus all polynomials of degree at least two. More strongly, it
admits no scalar Weyl pair with either member of Bernstein degree
at most two, regardless of the degree of the other member.

These are algebraic right-action constraints. They do not infer a
right-symbol frame from a classical h=0 map, exclude all higher-degree
right actions, or establish a bounded lift of Long's map.

## Normalize the symbol ideal

Let A=A_n(K) in characteristic zero, with Bernstein filtration,
division ring D and polynomial symbol ring
\[
              S=K[x_1,\ldots,x_n,p_1,\ldots,p_n].
\]
Degree and principal symbol extend to D by fractions. Products have
additive degree, and
\[
 \deg[a,b]\le\deg a+\deg b-2,\qquad
 \sigma_{\deg a+\deg b-2}([a,b])=\{\sigma a,\sigma b\}.
 \tag{1}
\]
The formula follows first for Weyl polynomials, then for inverses
from [a,b^{-1}]=-b^{-1}[a,b]b^{-1} and the quotient rule.

Any finitely generated torsion-free rank-one left A-module has a
fractional-ideal realization. Clear a common right denominator to
obtain a nonzero left ideal I in A. Let g be the homogeneous gcd of
gr I and choose G in A with principal symbol g. The fractional module
\[
                            M=I G^{-1}
\]
has induced symbol ideal
\[
                 J=\operatorname{gr}M=g^{-1}\operatorname{gr}I
                            \subset S.                 \tag{2}
\]
It has gcd one. In particular J localized at every height-one prime
of S equals that localization of S. This normalization is described
in [Kanel-Belov--Fedorov, Section 7.1](https://arxiv.org/html/2610.33919v1);
only PBW and Ore symbol calculus are needed below.

Write
\[
                         E=\{r\in D:Mr\subset M\}       \tag{3}
\]
for the right multiplier algebra, which is the opposite of the
left endomorphism ring. Since J has gcd one, (J:J)=S: a rational
multiplier of J lies in every height-one localization of the normal
polynomial ring S. Therefore
\[
                             \operatorname{gr}E\subset S.
 \tag{4}
\]
Indeed J sigma(r) is contained in J for every r in E.

A nonzero negative-degree element of E would contradict (4).
A degree-zero element has scalar principal symbol; subtracting
that scalar gives a negative-degree element and hence zero. Thus
E at degree at most zero is exactly K.

The module M is free if and only if J=S. If J=S, lift its degree-zero
generator 1 to b in M. Leading-symbol subtraction, using the
bounded-below filtration, proves M=A b. Conversely, a principal
ideal has principal associated graded and becomes S after gcd
normalization.

## Hamiltonian preservation is stronger than multiplication

For every h in gr E,
\[
                              \{J,h\}\subset J.        \tag{5}
\]
It suffices to prove this for a homogeneous principal symbol h of
a multiplier r in E. Let d=deg r>0 and choose H in A with symbol h.

For m in M with symbol j, the expression
\[
                 mr-Hm=m(r-H)+[m,H]                    \tag{6}
\]
belongs to M. If r-H has degree d-1, its leading symbol rho_1
therefore satisfies J rho_1 subset J, since the commutator in (6)
has degree at most deg(m)+d-2. By (J:J)=S, rho_1 is a polynomial.
Subtract a polynomial lift of rho_1 from r-H. Equivalently adjust H
by a polynomial of degree d-1. We may now assume
\[
                              \deg(r-H)\le d-2.
\]
Put rho equal to the degree-(d-2) symbol of r-H, or zero if its
degree is smaller. The degree-(deg(m)+d-2) part of (6) gives
\[
                              j\rho+\{j,h\}\in J.      \tag{7}
\]
For each height-one prime P choose j in J outside P. Since
{j,h} is polynomial, equation (7) makes rho regular at P.
Normality of S gives rho in S. Thus j rho already belongs to J,
and (7) proves (5). The argument applies to every homogeneous j
in J; linearity completes the proof. For d=1 the same argument
forces the possible degree-minus-one rho to vanish.

This is a necessary condition on the full graded right-action
algebra, not just on the chosen generator symbols.

## A transversality criterion for a free left frame

Let h_1,...,h_r belong to gr E, and let H_i={-,h_i}. Over the
algebraic closure, every point of V(J) satisfies
\[
                \operatorname{rank}(H_1,\ldots,H_r)<2n.
 \tag{8}
\]
If the fields had full rank at a point of V(J), some maximal minor
would be invertible in its local ring. The coordinate derivatives
would then be local linear combinations of the H_i. Equation (5)
would make the nonzero local ideal J stable under every coordinate
derivative. Choose a nonzero element and differentiate a lowest
nonzero Taylor term sufficiently many times. This produces a unit
in the ideal, a contradiction.

It follows that if those fields span everywhere, V(J) is empty,
J=S, and M is free. In particular it suffices that gr E, or the
full associated graded of the proposed right-image algebra,
contain polynomials f_1,...,f_(2n) with
\[
                           (\{f_i,f_j\})=J_{\rm symp}.  \tag{9}
\]
Their Hamiltonian fields span all coordinate derivations by the
polynomial symplectic Jacobian identity. Polynomial degrees are
unrestricted in this criterion.

For rank one, a proper normalized ideal has finite colength and
contains a power of every linear coordinate. If gr E contained a
nonzero linear form ell, equation (5) would preserve J under the
constant nonzero vector field H_ell. Repeatedly differentiating a
power of a linear coordinate transverse to that field would give
1 in J. Thus every nonfree rank-one module has no degree-one
element in its right idealizer.

An algebraic isomorphism of fractional realizations multiplies the
module on the right by an element of D^*. It conjugates E, preserving
degree and symbols. Gcd normalization also removes that fractional
factor from J. Thus (5)--(8) do not depend on this choice of algebraic
realization. Different algebraizations of one formal class need not
be isomorphic; this does not make the criterion an invariant of
all formal gauges.

## Include cancellations in the full right-symbol algebra

One must use all elements of the right-image algebra, not merely
the principal symbols of the displayed Weyl generators. For every
integer N>=2, the exact automorphism
\[
                         X=x,\qquad P=p+x^N
\]
has displayed principal symbols x and x^N, whose Hamiltonian fields
have rank one. But P-X^N=p also belongs to its right-image algebra.
The full graded image is S and contains the required transverse
linear frame.

Conversely, a classical h=0 symplectic frame is not necessarily a
frame of generic Bernstein symbols. [QEF1](ALGEBRAIC_DUAL_AND_POLARIZATION.md)
already shows why classical agreement alone cannot supply algebraic
normalization. Criterion (9) must be proved for the actual graded
right-image algebra.

## The exact charge-one idealizer

Work now in rank one with [v,u]=1. The generic QPF1 ideal is
A u^2+A(uv+2). Normalize it by right multiplication by u^{-1}:
\[
                      M=A u+A(v+u^{-1}),\qquad
                      J=(u,v).                         \tag{10}
\]
Its idealizer E lies in K[u,u^{-1}][v]. In fact u is in M, so
u r is in M for every r in E; since M is contained in this Laurent
differential-operator ring, so is r=u^{-1}(u r).

Set
\[
                    H=uv,\qquad L=v^2-2u^{-2}.
\]
Right multiplication by u^2,u^3,H,L preserves M. With m_1=u and
m_2=v+u^{-1}, explicit polynomial left coefficients are
\[
 \begin{array}{c|cc}
 r&m_1r&m_2r\\ \hline
 u^2&u^2m_1&(vu+1)m_1\\
 u^3&u^3m_1&(vu^2+u)m_1\\
 H&-m_1+u^2m_2&v^2m_1\\
 L&-v^2m_1+2uvm_2&v^2m_2 .
 \end{array}                                           \tag{11}
\]
Their principal symbols include all three quadratics u^2,uv,v^2.
Successive commutators with L applied to u^3 provide all cubic
symbols. Products of quadratics and cubics provide every homogeneous
polynomial of degree at least two.

Equation (5) applied to J=(u,v) excludes every linear symbol, while
(4) excludes negative-degree symbols. Hence
\[
                 \operatorname{gr}E=K\oplus S_{\ge2}.  \tag{12}
\]
This is an exact all-degree description, not a truncated idealizer
calculation. It also computes codim_S(gr E)=2, consistent with the
charge-one identity lambda=2 nu discussed in the cited paper.
No arithmetic nonlifting theorem is used.

## No quadratic member of a Weyl pair in this idealizer

Suppose P,Q in E satisfy [P,Q]=kappa in K^*. Both have degree at
least two by (12). We show neither can have degree two.

Assume deg P=2 and let f be its quadratic principal symbol. If
f has two distinct linear factors over the algebraic closure,
the leading commutator equation forces
\[
                        \{f,\sigma(Q)\}=0.
\]
A linear coordinate change makes f a nonzero multiple of uv.
A homogeneous polynomial commuting with uv is a scalar power of uv,
so deg Q is even and sigma(Q)=a f^(deg(Q)/2).
Subtracting a P^(deg(Q)/2) from Q decreases its degree and preserves
the nonzero constant commutator. Repeating reaches degree at most
one, where E consists only of scalars, a contradiction.

It remains to treat f=a ell^2 for a nonzero linear form ell.
The three elements
\[
                 e=u^2/2,\qquad f_0=-L/2,\qquad h=H+1/2
\]
satisfy the sl_2 relations
\[
 [h,e]=2e,\qquad [h,f_0]=-2f_0,\qquad [e,f_0]=h.         \tag{13}
\]
Their adjoint actions preserve each finite-dimensional Bernstein
piece of E. Since E consists of Laurent differential polynomials,
ad h is diagonalizable with integer monomial weights. The raising
and lowering actions ad e and ad f_0 are consequently locally
nilpotent on each such piece. Their finite exponentials are
filtered automorphisms of E. On principal symbols they give the
two elementary symplectic linear shears, so they can send ell
to a nonzero multiple of u.

After such an automorphism, P=a u^2+b exactly: its difference
from a u^2 has degree at most one and hence is scalar by (12).
Regard Q as a Laurent differential operator in v. The leading
commutator with u^2 lowers positive differential order by exactly
one. A nonzero scalar commutator therefore forces
\[
                  Q=-\frac{\kappa}{2a u}v+r(u),
                     \qquad r(u)\in K[u,u^{-1}].        \tag{14}
\]
Every multiplication operator u^j for j>=2 belongs to E, by
(11) and products. Subtract the terms of r(u) with exponent at
least two. The remaining Q still belongs to E, has Bernstein
degree at most one, and still has commutator kappa with P.
Equation (12) again makes it scalar, a contradiction.

Thus every possible scalar Weyl pair in this explicit E must have
both Bernstein degrees at least three. The other member's degree
was never bounded. This statement concerns the specified normalized
Bernstein filtration; it does not exclude all embeddings A_1->E
and it is not a theorem about differential order two on arbitrary
coefficient curves.

## Use in the remaining Dixmier program

For a candidate finite projective graph model, normalize its symbol
ideal and study the full right-symbol algebra. Transversality of its
Hamiltonian fields removes the left-frame obstruction by (8). Failure
of transversality identifies where a boundary defect could remain;
it is a necessary condition, not evidence that such a defect exists.

This supplements [QPF1's parameter-compatible filtered-frame criterion](PROJECTIVE_FRAME_BOUNDARY.md).
The generic symbol criterion gives generic left freeness. Obtaining
a frame over the parameter ring with the prescribed classical
reduction still requires parameter descent and compatibility. Neither
criterion constructs the right actions in the first place.

The positive target for Long remains an actual algebraic model,
its full right Weyl action, and a compatible frame giving one finite
spatial bound as required by [QFL2](BOUNDED_REALIZATION.md). The
classical frame alone does not verify (9), and the explicit
charge-one calculation does not exclude every possible Long model.

[Exact controls](check_hamiltonian_idealizer.py) verify the rational
Ore identities in (11), the sl_2 commutators, the cubic-symbol ladder,
the quadratic-symbol centralizer identity and the scalar commutator
in (14). Normalization, Hamiltonian preservation, transversality,
the all-degree idealizer description, local nilpotence of the
adjoint actions and the unrestricted-partner exclusion use the
written proofs. No independent theorem replay, formal verification,
external review, novelty priority, bounded Long lift or DC2 solution
is claimed.
