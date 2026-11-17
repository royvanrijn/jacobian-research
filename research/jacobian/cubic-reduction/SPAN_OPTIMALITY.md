# The saved 89-dimensional scalar-cubic span cannot give fewer cubes

The later [two-new-form obstruction](TWO_CUBE_OBSTRUCTION.md) also excludes
the unrestricted two-form escape route identified below.

Let P be the fixed padded eighteen-output cubic tensor in
[KCUB5](POWER89.md), and let S be the scalar cubic span of the saved89
selected cubes. Over every characteristic-zero extension field:

**Every shared-cube representation of P whose scalar cubes all lie
in S spans all of S. In particular it uses at least89 cubes.**

This allows all forms and all output coefficient columns to change;
it is stronger than the [single-new-form obstruction](SINGLE_CUBE_OBSTRUCTION.md).
It does not prove89 minimal for P among arbitrary representations, or
for cubic-linear Keller maps. A smaller representation must use cubes
outside S. The unrestricted fixed-tensor interval remains
[34,89](FIXED_TENSOR_BOUND.md). No border-rank statement is made.

## All cubes in S

Write a=z2, b=z4, t=z18 and wi=z(12+i), for 1<=i<=5. The projective
linear forms whose cubes lie in S are exactly:

- 65 isolated points, each already among the selected forms;
- the eight projective lines
  span(b,t), span(a,b+t), span(a,b-t), and span(wi,t) for i=1,...,5.

The [complete classification](power-span-split/classification.json)
contains rational bases for every point and line, and the indices of
all original forms. Of the89 selected points,65 are isolated and24 lie
on these lines. The isolated cubes and the24 selected line cubes give
a direct-sum basis of S. All cubes on any of the eight lines lie in the
24-dimensional span of those24 cubes.

The classification is certified by408 nodes of rational linear-factor
splitting trees, including236 factor identities. On each source-support
case, a nonzero homogeneous membership equation factors into linear
forms; every zero must lie in at least one corresponding hyperplane.
At each leaf every membership equation vanishes identically on the
retained linear subspace. Thus the tree proves equality of the zero
set with its leaf subspaces over any characteristic-zero extension.
No primary-decomposition routine is trusted by the final checker.

The15 coordinate-support cases cover every possible form: its nonzero
source coordinates must be a clique of the co-occurrence graph. The
independent checker exhausts all2^18 coordinate subsets, reconstructs
the membership-equation row spaces directly from the89 cubes, checks
every factorization and hyperplane branch, and merges the leaf spaces.

## An annihilator contradiction

Suppose a representation of P by cubes in S spans a proper subspace.
Choose a nonzero linear functional Lambda on S annihilating every cube
in that representation. It also annihilates every output of P.

Every isolated point is compulsory: its cube is independent of all
other isolated cubes and all line cubes, and its original coefficient
column in P is nonzero. Hence Lambda vanishes on all65 isolated cubes.
The padding output of P is t^3, so Lambda(t^3)=0.

For each wi, the only terms of P involving wi are -t^2 wi in output i.
Only the line span(wi,t) can contribute such terms. At least three
distinct forms wi+r t with nonzero wi coefficient are necessary: with
at most two, canceling the wi^3 and wi^2 t coefficients also cancels
the wi t^2 coefficient. Since Lambda(t^3)=0, the polynomial

\[
r\longmapsto\Lambda((w_i+r t)^3)
\]

has degree at most two. It vanishes at at least three distinct values,
so it is identically zero. Thus Lambda vanishes on the full cubic span
of each of these five binary lines.

The remaining annihilator conditions have rank87 in the89-dimensional
dual-coordinate calculation. Their two-dimensional solution space is
Lambda=u D+v E, where D and E vanish on the65 isolated cubes and, on
the24-dimensional line-cube space, extract the coefficients of bt^2
and at^2 respectively. The checker verifies both the rank and these
explicit functionals. Their restrictions to the three remaining lines
are

\[
\begin{aligned}
\Lambda((x b+y t)^3)&=3u x y^2,\\
\Lambda((x a+y(b+t))^3)&=3y^2(u y+v x),\\
\Lambda((x a+y(b-t))^3)&=3y^2(u y+v x).
\end{aligned}
\]

If u is nonzero, the zero cubes on these lines are contained in the
five projective points

\[
[a],\ [b],\ [t],\ [u a-v(b+t)],\ [u a-v(b-t)].
\]

Their scalar cubic span has dimension at most five. If u=0, then v is
nonzero: the zero cubes are the whole line span(b,t), plus [a], again
with scalar cubic span at most five.

Project S onto its nine-dimensional special-line summand, using the
nine original cubes on the three special lines as a basis. The isolated
cubes project to zero. Cubes on the five wi lines project to multiples
of t^3, already included in the preceding five-dimensional bounds.
But the projected eighteen-output tensor has scalar output rank seven,
as checked on the original coefficient matrix. It cannot be spanned by
the projected Lambda-zero cubes. This contradiction proves the claim.

## Consequence for two new forms

Suppose a shorter representation retains some original forms and adds
two new ones. In the quotient of the scalar cubic space by S, their
output coefficient columns multiply their two residual cubes to zero.
If those two columns were independent, both residual cubes would be
zero, contradicting the theorem. Thus a possible improvement with two
new forms requires **proportional output columns** and **two cubes
outside S whose nonzero residuals cancel**. Every deleted original
coefficient column must have that same projective direction.

This gives a concrete next construction problem, rather than another
basis change inside S.

## Exact finite-pool two-form test

The [two-form exchange calculation](power_two_exchange.py) combines the
previous weighted and quadratic-pencil pools into930 distinct forms.
Their cubes have rank154. There are151 cubes inside S and445 projective
directions among the nonzero residual cubes. Only1,774 outside-form
pairs have proportional residuals. For each pair, the unique residual
cancellation is tested with all retained coefficients freely re-solved.
At most two old forms cancel, so no shortening occurs in this finite
pool. [Primary receipt](power-two-exchange/checkpoint.json).

A [separate Fraction implementation](verify_power_two_exchange.py)
reconstructs the quotient projection and verifies every admissible pair.
[Independent receipt](power-two-exchange/validation.json). The151
inside-span candidates use the theorem above; pairs with independent
residuals or one zero residual cannot improve the representation.
This is not an unrestricted two-new-form exclusion.

The same independent replay checks214 proportional-column raw triples.
In209 the three source forms are independent, giving scalar derivative
rank three; the remaining five are the exact t^2 wi binary monomials,
also requiring three cubes. Thus these raw triples do not give a
three-to-two replacement with all other coefficient columns fixed.
Re-solving those columns in a larger, unrestricted form search remains
open.

## Reproduction and status

[Primary factor-tree construction](power_span_split.py) uses SymPy only
for exact factor selection. [Independent verification](verify_power_span_optimality.py)
uses standard-library Fraction arithmetic and the retained sparse
polynomial engine. It checks the full classification and every linear
algebra premise of the written annihilator argument. Its
[validation receipt](power-span-split/validation.json) pins the trees,
classification, sources and input witness. The separate
[primary decomposition](power_span_components.py) is retained as discovery
evidence, not a premise of the proof.

The integer89-variable cubic-linear witness is unchanged. Ordinary
cubic12/36, normalized three-variable support[10,48], and global
simplicity minima are separate, still-open targets.
