# The Wright connection and the multiple fibre gate

14 November 2026. This is a literature lead and an explicit geometric
control for the actual-normalization compatibility question. The later
[period obstruction](VOLUME_PERIOD_OBSTRUCTION.md) and
[marked-ring comparison](WRIGHT_MARKED_MODEL.md) now identify the sole
surviving compatible boundary pair with Wright's index-three ring inside
the actual source algebra. They do not establish a compatible embedding
for an arbitrary actual normalization or exclude the remaining ring problem.

## The relevant older problem

Rodriguez Diaz, [On the integrality of etale extensions of polynomial rings,
arXiv 2403.02219v2](https://arxiv.org/html/2403.02219v2), discusses an
integral closure that is a primary submodule of the source polynomial ring.
Theorem 3.5 claims an affine-line bundle over P1 preserving the source open;
Proposition 3.6 identifies index three. Its discussion of Wright's 1997
work gives affine-line-bundle coordinate rings and identifies an unresolved
index-three case. The final Theorem 3.12 assumes the primary-submodule
property for every etale extension; it is not a theorem excluding one
arbitrary specified quartic packet.

This places our forced N=F1\B in the same abstract index-three surface
class. That literature comparison alone does not identify the marked
embedding U inside N or establish
that a proposed affine-line ruling preserves the actual omitted prime
with multiplicity one.

Wright's original paper is *Affine surfaces fibered by affine lines over
the projective line*, Illinois J. Math.41(1997),589-605,
[DOI 10.1215/ijm/1256068982](https://doi.org/10.1215/ijm/1256068982).
The publisher PDF still gives a service check. An author-uploaded primary
copy was readable in the later period/marking pass; the
[access record](sources/wright-original-access.json) preserves that
provenance and its unpinned-content limit. No theorem from the original
paper is adopted as a dependency of the new algebraic identities.
The pinned 2024 paper is a source to audit, not an additional dependency
of PF2RULED2.

## A smooth omitted affine line does not determine fibre multiplicity

Consider two affine-plane charts with transition

\[
 u=x^{-1},\qquad v=x^3y+x,
 \qquad x=u^{-1},\quad y=u^3v-u^2.
\]

They define an affine-line torsor over P1. In the second chart, let
R be u=0, with coordinate v. The first chart U=A2 is its complement.
For example the four functions

\[
 y,\quad xy,\quad x^2y,\quad x^3y+x
\]

extend across R with expressions

\[
 u^3v-u^2,\quad u^2v-u,\quad uv-1,\quad v.
\]

Now use the polynomial coordinates on U

\[
 A=x^2+y,\qquad D=x+A^2.
\]

Their Jacobian is -1 and the inverse is x=D-A^2, y=A-(D-A^2)^2. Along R,

\[
 A=u^{-2}(1-u^4+u^5v),\qquad
 D=u^{-4}\big((1-u^4+u^5v)^2+u^3\big).
\]

Both bracketed numerators are units along the entire R. Thus both
coordinate projections extend to maps into P1, but their infinity
fibres have multiplicities two and four along R. In particular the
reduced curve underlying either fibre is the same smooth A1.

These A,D are source coordinates on U; they are not regular functions
on N and do not supply a map N -> affine target A2. The example is a
control for a geometric inference, not a Keller counterexample.

Consequently a passage from a reduced omitted curve to an affine-line
bundle must prove the multiplicity of the chosen fibre. A surface
isomorphism or a different coordinate projection does not preserve that
multiplicity automatically. This is the precise gate to audit before
using a Wright-style coordinate ring as the actual finite normalization.

## Next use of this connection

Read the original Wright construction and check which hypotheses preserve
the given open embedding U, the omitted divisorial valuation and the
scheme-theoretic fibre. Then compare its boundary valuation with the B
used in the [pole reduction](POLE_REDUCTION.md). Pole height depends on
that completion; an abstract isomorphism of index-three surfaces does
not transport the bound without a valuation comparison.

The [checker](verify_pole_reduction.py) verifies the displayed transition,
global function expressions, coordinate inverse, Jacobian and pole
multiplicities exactly. It does not certify the existence of a compatible
bundle ruling for an actual quartic normalization.

