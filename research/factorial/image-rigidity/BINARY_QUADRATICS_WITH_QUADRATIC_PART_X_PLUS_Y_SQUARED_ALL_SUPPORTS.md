# Binary quadratics with quadratic part a(x+y)^2 and all supports

Let K be a characteristic-zero field and let
Lambda_2(x^i y^j)=i!j!, Lambda_1(X^i)=i!. For a nonzero polynomial

 f=a(x+y)^2+d x+e y+k

its actual support size is
s=3*1_(a!=0)+1_(d!=0)+1_(e!=0)+1_(k!=0).
The three quadratic coefficients are a,2a,a and are simultaneously nonzero.

**Theorem 1 (the full displayed family).** For all a,d,e,k in K, provided
f!=0, and every n>=1, at least one of
Lambda_2(f^n),...,Lambda_2(f^(n+s-1)) is nonzero.

The new zero-constant and affine arguments are proved below. The remaining
cases import exactly the three written sector theorems specified in Section3.
No general rank-one orientation A(ux+vy)^2 or general binary quadratic theorem
is asserted. Independent rescaling of variables does not preserve Lambda_2.

## 1. Unequal linear coefficients and zero constant

Write M_r=Lambda_2(f^r), now with k=0, a!=0 and d!=e. Define
H_j(t)=Lambda_1((aX^2+tX)^j), H_0=1. The monomial Stein identity
Lambda_2(partial_x P)=Lambda_2(P)-Lambda_1(P(0,y)), and its y version,
give the polynomial identity

 (r+1)(d-e)M_r=H_(r+1)(d)-H_(r+1)(e), r>=0.              (1)

The native one-variable factorial recurrence, obtained from
Lambda_1(P')=Lambda_1(P)-P(0) applied to powers of aX^2+tX, is

 H_j(t)=2a j(2j-1)H_(j-1)(t)+j(j-1)t^2 H_(j-2)(t),
                                                        j>=2. (2)

Its initial values are H_0=1,H_1=2a+t. These are polynomial identities;
there is no division by a or a Hessian.

**Lemma 2 (zero-constant windows).** If a!=0,d!=e and k=0, the displayed
family has no four consecutive zero positive moments. If d+e=0, it has
the stronger exclusion of two consecutive zero moments.

First suppose S=d+e!=0. Four consecutive zero M values starting at n>=1
give, by (1), common values
R_l=H_l(d)=H_l(e) for l=n+1,...,n+4.
Subtract (2) at index l+2 for l=n+1,n+2. Since H_(l+2) and H_(l+1)
also agree, the difference is

 (l+2)(l+1)(d-e)(d+e)R_l=0.

The nonzero integer factors, d-e and S show R_(n+1)=R_(n+2)=0.
At least one of d,e is nonzero. For that boundary, the coefficient
j(j-1)t^2 in (2) is nonzero for every j>=2. A pair of adjacent H zeros
therefore propagates backward to H_0=0, contradicting H_0=1.

If S=0, then e=-d and d!=0. Put D_j=H_j(d)-H_j(-d).
The two H recurrences have the same coefficients, so D_j satisfies (2)
with t=d, for j>=2. Its initial values are D_0=0,D_1=2d!=0.
A pair of adjacent D zeros propagates backward to D_1=0, a contradiction.
By (1), two adjacent positive M zeros would give just such a pair.
This proves Lemma2.

If exactly one of d,e is zero, f has four terms, and the four-window is
its required window. If both are nonzero, it has five terms, and Lemma2
is stronger than its required window. No shorter window is inferred by
discarding a coefficient from another theorem.

## 2. All affine supports

**Lemma 3 (affine windows).** Every nonzero f=dx+ey+k satisfies the
Strong Factorial window of its actual support.

Put A_r=Lambda_2(f^r)/r!, A_0=1. The formal generating function is

 sum_(r>=0) A_r T^r=exp(kT)/[(1-dT)(1-eT)].

Indeed this is the product of the two geometric factorial integrals and
the constant exponential; equivalently it follows by literal multinomials.
Consequently, with negative-index terms omitted,

 A_r=(d+e)A_(r-1)-de A_(r-2)+k^r/r!.                    (3)

If d,e,k are all nonzero, a triple of consecutive zero coefficients would
give k^(n+2)=0 by (3), which is impossible. This is the actual three-window.

If k=0 and de!=0, (3) is homogeneous from r>=1, and its nonzero backward
coefficient de forces an adjacent zero pair back to A_0=0. Thus the actual
two-window holds. If exactly one linear coefficient and k are nonzero,
the generating function has a single denominator and the first-order identity

 A_r=d A_(r-1)+k^r/r!

(or its e version) excludes adjacent zeros. A sole nonzero monomial or
constant has every positive factorial moment nonzero. These cases exhaust
all nonzero affine supports and prove Lemma3.

## 3. Exhaustive coefficient partition and exact written imports

The following three written theorems are the only imported sector results.

(A) BINARY_QUADRATIC_OF_X_PLUS_Y_STRONG_FACTORIAL.md, in the same directory,
SHA256 49781e8ec35c9c579dabd488b23c6ecd464b9803d331513c8b151fc16f48ecb9,
proves the actual-support Strong Factorial window for every nonzero
g(x+y), deg(g)<=2, over every characteristic-zero field. Its exact theorem
handle is FACTORIAL-QUADRATIC-OF-X-PLUS-Y-STRONG-FACTORIAL.

(B) The complete written six-window proof
RANK_ONE_UNEQUAL_LINEAR_BINARY_QUADRATIC_STRONG_FACTORIAL.md (same directory),
SHA256 7ee19142bbf522a8a032597789cf5ca97efdfd65baa8e0b969917721ae133cdc,
proves no six consecutive zero positive moments for
a(x+y)^2+dx+ey+k when adek!=0 and d!=e. It includes the
normalized-coefficient number-field reduction and both the odd-prime
divisor and power-of-two index exclusions. Its coefficients and field
quantifiers are precisely those used here.

(C) The complete written weighted constant-square proof
WEIGHTED_CONSTANT_SQUARE_MOMENTS_AND_FIVE_TERM_QUADRATICS_IN_X_PLUS_Y.md (same directory),
SHA256 ca8ca48f911d87723299b84c43fec2da4ab6290d40b53f42c74093e6b920fe68,
proves adjacent coprimality for
K_j(A)=sum_(i=0)^j binom(j,i)(i+1)(2i)! A^i over every
characteristic-zero field, and deduces no five consecutive zero positive
moments for a(x+y)^2+dx+k with adk!=0. This is the exact (i+1)-weight.
The same theorem with x,y swapped proves the e-only case, because
Lambda_2 is invariant under this swap.

To prove Theorem1, partition all nonzero displayed polynomials as follows.

* If a=0, Lemma3 gives exactly the actual-support window.
* If a!=0 and d=e, write f=g(x+y) with g(T)=aT^2+dT+k.
  Import(A) covers every zero/nonzero pattern of d and k with its actual
  support in x,y, not the number of terms in g.
* If a!=0, d!=e and k=0, Lemma2 gives a four-window. The actual
  support is four when exactly one of d,e is nonzero, and five when
  both are nonzero. (They cannot both be zero when d!=e.) The four-window
  implies the required window in both cases.
* If a!=0, d!=e, k!=0 and de!=0, all six displayed monomials occur,
  so import(B) is exactly the required six-window.
* If a!=0, d!=e, k!=0 and de=0, exactly one linear coefficient is
  nonzero. The actual support is five; import(C), or its variable swap,
  gives exactly that five-window.

These cases are mutually disjoint and exhaustive. Every proof applies
over the original characteristic-zero field, using only the normalized
finite number field when a local valuation is required. No limiting
argument, coefficient deletion or arbitrary variable rescaling is used.
This proves Theorem1.

## 4. Scope and provenance

This is a written mathematical assembly and contains two new elementary
support lemmas, not a formal proof-assistant certificate or a priority
claim. The full binary quadratic Strong Factorial conjecture, SFC(2)
and SFC(3) do not follow: a general quadratic has independent coefficients
of x^2,xy,y^2, and even a general rank-one form A(ux+vy)^2 need not have
the fixed orientation x+y treated here. The stronger windows in Lemma2
are sufficient statements; no optimality is claimed.

The three source bytes and hashes above are pinned as written imports.
By restriction, Theorem1 implies the actual-support statements of all
three imported sector theorems and Lemma3. The converse assembly uses
all five cases. No reduction to one proper sector, or reverse implication
from the actual-support theorem to the stronger four- and two-windows of
Lemma2 on five-term supports or to adjacent K coprimality, is established
here. Those stronger conclusions are proved independently. The dependency graph uses the independent sector
proofs as inputs to the assembly; these restriction observations do not
supply a circular proof route.
