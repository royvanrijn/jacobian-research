# Complete all-degree closure of positive p dependence in the y graph

Let `y=R(x,p)` in the scaled Meng--Yang potential, with `L*N != 0`.
After taking this genuine polynomial graph, allow every actual polynomial
automorphism of the retained three duals over K[x], with polynomial inverse.
There is no bound on either degree of R or the automorphism word length.

**HC4MYGY2. If deg_p R>0 and the resulting actual four-variable Hessian
has nonzero constant determinant, its gradient is a polynomial automorphism.**
The remaining no-simple-boundary-root case in fact forbids unit Hessian.
This closes the whole positive-p y chart, including repeated Newton faces,
exceptional leading cancellation and all finite-base escapes. It is not
a degree-by-degree continuation of the quadratic argument.

**HC4MYGY3. In the entire family y=R(x,p), including p-independent R,
existence of a nonautomorphic unit-Hessian member is exactly existence
of a unit-ideal plane Keller pair.** Thus this whole chart has only the
previous base-only plane-surjectivity endpoint. This is a statement about
**nonautomorphic** members, not nonexistence of every unit-Hessian member.
The plane endpoint remains open, as do the separate nonlinear-p graph
singleton-zero endpoint and other omitted-base/mixing classes.

## 1. The remaining hypothetical case

Work over C. Any characteristic-zero candidate and its retained-dual
inverse use finitely many coefficients and may be embedded into C.
Invertibility, when obtained, descends by uniqueness of the polynomial
inverse. Set

\[
 F=3+2xR,\quad
 D=-\frac{Np}{4x^2}-\frac{(M-13N/(2x^2))^2}{4L},\quad
 \tau=\frac{dD}{dx},\quad E=D-x\tau.
\]

The previous [boundary-root theorem](CRITICAL_COMPONENTS.md) proves
gradient invertibility if R0=R(0,p) has a simple zero. An identically zero
R0 is excluded by [HC4MYGU0](RELATIVE_COORDINATE_ENDPOINT.md), and a
multiple zero at `p=-M/(2L)` is excluded by its critical two-plane.
It therefore suffices to rule out unit Hessian when

\[
 R_0\ne0,\qquad R_0\text{ has no simple zero},\qquad
                         R_0(-M/(2L))\ne0.                     \tag{1}
\]

Assume such a unit-Hessian retained-dual pullback exists. Let S be the
product of the p-dependent multiplicity-one factors of F. By
[HC4MYGM0/1/3](CRITICAL_NORMALIZATION.md), the entire relative critical
curve C=S=0 is nonempty, smooth, closed in the (x,p) plane, disjoint
from all other F components, and etale under x to C*. In particular
x and R_p are units on C, and

\[
            K=-32LN^2x^6R_p^2\frac{d\tau}{dx}\ne0              \tag{2}
\]

is one constant on the whole curve. The retained Tang rank-one
realization identifies the entire C with a plane Keller fibre B=0,
carrying tau to its Keller mate and E to a plane polynomial.
Its generic dual value fibres are those of `E+s*tau+V*B`, with
scalar s independent of the function x on C. These are exactly the
already established hypotheses, not a new plane classification.

## 2. Compare every Euler jump (HC4MYGY0)

Let a=deg_p S. Let b count the other distinct p roots of F, called
inactive roots, and let v count these roots with their multiplicities.
Thus `r=a+b`, `n=deg_p F=a+v`, and `v>=2b`.

For generic fixed x=s nonzero, the native three-dual value fibre has

\[
 \chi(f_s^{-1}(\lambda))=1+r-\#\{p_i:F(s,p_i)=0,\ D(s,p_i)=\lambda\}, \tag{3}
\]

where roots are counted distinctly. Indeed the open p stratum contributes
1-r; over each root the quadratic equation in T has two distinct roots
except at its value D, when it has one. The q coordinate is free there.
This formula remains valid for inactive roots of any multiplicity.

In the rank-one realization the same fibre has Euler
`1-chi(C)+#(E+s*tau=lambda on C)`. On C,

\[
                  d(E+s\tau)=(s-x)\,d\tau.
\]

Both dx and d tau are nowhere zero on the affine C, by (2). Thus its
only finite critical points have x=s, and all are simple ramification
points. They are precisely the a active roots, with the same values D
as in (3). For generic s, all distinct p roots have distinct values D,
because D(s,p) is affine with nonzero coefficient in p.

Take projective normalizations of the C components. The remaining
defects of the value-map fibre count are the boundary points q where
both E and tau are finite. Write their values alpha_q,beta_q, and let
k_q be the generic local degree of E+s*tau at q. Comparing (3) with
the compact value-map fibre count, after canceling the simple active
ramification contributions, gives the exact identity

\[
 \prod_{i\ {m inactive}}(Y-D_i(x))
     =\prod_q(Y-\alpha_q-\beta_qx)^{k_q}.                       \tag{4}
\]

This equality first holds for generic scalar x=s and hence as a monic
polynomial identity over C(x). The left side is squarefree, since the
inactive roots are distinct and D is affine in p. Consequently every
k_q=1, the affine functions on the right are distinct, and there is
exactly one such boundary point for each inactive root. In particular
every inactive root is rational and has

\[
 D_i=\alpha_i+\beta_ix,\qquad
 p_i=p_*(x)-\frac4N x^2(\alpha_i+\beta_ix),\quad
 p_*=-\frac{169N}{4L}x^{-2}+\frac{13M}{L}-\frac{M^2}{LN}x^2.    \tag{5}
\]

These boundary points cannot be finite-base escapes: there E,tau have
poles. They cannot lie over x=infinity either. For if x=t^-e and
`tau-beta` has order d, the identity `dE=-x*d tau` and regularity of
both functions require d>e, with `ord(E-alpha)=d-e`. Since k_q=1,
d=e+1. Hence `ord(d tau/dx)=2e+1`. The valuation of (2) would then be
`-6e+2*ord R_p+2e+1`, an odd integer, never zero.

All boundary points in (4) therefore lie over x=0 and have D regular
there. In particular **b>0 forces at least one active zero-base branch
on which D is regular**. This is a statement about the entire value
fibre, not a theorem turning one rational special fibre into a generic one.

## 3. One root-product inequality eliminates every regular-D branch (HC4MYGY1)

At x=0 every p root of F has a pole, because F(0,p)=3. Count all
Puiseux roots with their ramification multiplicities. Among active roots:

- A roots have D regular at zero;
- ell roots have p-pole exponent lambda<2; their exponent sum is mu;
- P=a-A-ell are the remaining roots, all with D having a pole.

Here `mu>0` when ell>0. Every inactive root has D regular by (5), and
therefore the same leading term `p~c*x^-2`, where `c=-169N/(4L)`
uses the original potential coefficient L.
So do the A regular-D active roots.

The [global pole budget](GLOBAL_CRITICAL_BUDGET.md) already gives

\[
              a+b\ge5\ell+P,\qquad\text{or}\qquad4\ell\le A+b. \tag{6}
\]

The low-pole roots cost5 per sheet; every other pole-D branch costs
more than its x index. All other boundary contributions are nonnegative.

Suppose A>0 and choose a regular-D active branch with primitive
parameter x=t^e. Write D as a power series in t, and let d be the
least exponent with nonzero coefficient other than 0 and e. It exists
because D'' cannot vanish in (2). Then

\[
 \operatorname{ord}_t D''=d-2e,\quad
 j:=\operatorname{ord}_tR_p=-2e-d/2.
\]

Since j is an integer, **d is even and d>=2**. Under conjugation
`t -> zeta*t`, zeta^e=1, the constant and linear-in-x terms cancel,
so differences between conjugate D functions have order at least d.

Use the reciprocal polynomial `H(x,w)=w^nF(x,1/w)`. Its leading
coefficient in w is F(x,0), a unit at zero. All its n roots w_i=1/p_i
tend to zero. On the chosen simple active root, H_w is that unit times
the product of its differences from **all** other w roots, counting
inactive multiplicity. The identity

\[
 w_i-w_j=\frac{4x^2}{N p_i p_j}(D_i-D_j)                        \tag{7}
\]

gives these lower bounds in the chosen t valuation:

| Other roots | Number | Difference order at least |
|---|---:|---:|
| Conjugates of the chosen branch | e-1 | 6e+d |
| Other regular-D active roots | A-e | 6e |
| Inactive roots, with multiplicities | v | 6e |
| Remaining non-low active roots | P | 2e |
| Low-pole active roots | ell | total e*mu |

For the regular-D rows, p_i p_j has pole order4e and D_i-D_j is
regular. The conjugate row also uses its extra order d. For a low-pole
root, w_j dominates w_i and has order e*lambda. The remaining roots
have w order at least2e. These rows exhaust n-1 roots.

On the other hand `F_p=-w^(2-n)H_w` and `R_p=F_p/(2x)` give the
exact order

\[
                   \operatorname{ord}_tH_w=(2n-5)e-d/2.
\]

Summing the table and using n=a+v gives

\[
 (2n-5)e-d/2\ge(e-1)(6e+d)+6e(A-e)+6ev+2eP+e\mu,
\]

or equivalently

\[
 2\ell-\mu\ge4(A+v)-1+(1-1/(2e))d\ge4(A+v).                  \tag{8}
\]

Twice (8) and (6) would force `0>=2mu+7A+8v-b`. Since v>=2b,
the right side is at least `2mu+7A+15b`, strictly positive when A>0.
This contradiction proves **A=0**, in all ramification indices and
degrees at once.

Section 2 now forces b=0. All active zero-base branches have D with a
pole. Each contributes strictly more than its x index to the global
budget. Their indices sum to a, so the budget's right side is strictly
greater than a=r, again a contradiction. This proves that no unit
Hessian exists under the remaining boundary assumptions (1).

The proof uses no finite-base escape assumption: all escape contributions
were retained in the global budget and only discarded as nonnegative
when deriving (6). It uses no bound on repeated leading faces, genus,
Puiseux characteristic exponents, or inactive multiplicities.

## 4. The full chart and its exact remaining endpoint (HC4MYGY2/3)

Combining Section 3 with the boundary alternatives in Section 1 proves
HC4MYGY2 for every positive p degree. More specifically, if F has any
active p-dependent factor, unit Hessian is impossible: a simple boundary
root would give a separate A1 critical component and hence reducibility,
while every other boundary case is now excluded. If all p-dependent
factors of F are repeated, a unit member can only survive through the
previous simple-boundary-root invertibility case. No sufficiency or
unit-Hessian construction for that repeated case is claimed.

If deg_p R=0, the graph is base-only. R(0)=0 forbids unit Hessian;
R(0) nonzero gives precisely all relative polynomial coordinates by
[HC4MYGU1](RELATIVE_COORDINATE_ENDPOINT.md). Such a coordinate has an
empty relative critical scheme, so any unit-Hessian gradient in that
class is nonautomorphic. HC4MYGU2/3 identify its existence exactly with
a unit-ideal plane Keller pair, equivalently a nonsurjective plane
Keller map. Conversely each such pair is realized in every nonzero
boundary base-only relative-coordinate orbit. Thus the **whole**
p-only y graph has this exact nonautomorphic endpoint; positive p
dependence adds none.

## Verification and remaining scope

[verify_y_graph_closure.py](verify_y_graph_closure.py) and its
[receipt](receipt-y-graph-closure.json) check exact identities and controls for
the jump defect count, the odd infinity valuation, the root-product
coefficient bookkeeping and inequality, and controls distinguishing
active root multiplicity from inactive degree. The every-value Euler
comparison and all-degree Puiseux-root product argument are written
proofs, not conclusions of a degree search. Inputs are the retained
Tang alternative, the earlier whole-curve normalization and Euler
proofs, and the explicit classical line imports of the simple-boundary
case. No additional one-fibre rationality theorem is imported.
No formal verification, independent external review or novelty priority
is asserted.

Arbitrary y=R(x,p) is closed for nonautomorphic unit-Hessian gradients.
General x=R(y,p) beyond the affine-p class, q/r-dependent omitted bases,
genuine free-base/dual mixing and uncovered embeddings remain. The
unit-ideal plane-surjectivity and separate singleton-zero properness
endpoints remain unresolved. No unrestricted HC4 or JC2 solution is claimed.
