# Polynomial key chains exclude actual omitted primes at every degree

PF2RULED100 extends the [PF99 minimal-degree monic-key lemma](MONIC_KEY_JACOBIAN_OBSTRUCTION_AND_WHOLE_CHART_NECESSITY.md)
to an arbitrary finite chain of polynomial key forms. In particular, an
actual omitted prime of a normal affine surface containing the marked whole
plane must have a genuinely Laurent terminal key form. The same obstruction
applies directly to a plane Keller pair regular at a valuation, without
assuming an affine extension. No field degree, discrepancy, smoothness,
number of other primes or finite-generation hypothesis enters the algebraic
argument. This is a necessary obstruction, not a classification of the
remaining Laurent valuations or a proof of the full Jacobian conjecture.

Published [Domrina2000 / HC4MYKR0](../../hessian/hc4-triangular-graphs/QUADRATIC_SOURCE_RECONSTRUCTION_AND_QUARTIC_RETIREMENT.md)
already excludes every generic-degree-four plane Keller realization. The
new result concerns all degrees; it does not certify the historical internal
quartic proof chain. The broader Wright objective remains active.

## 1. Exact hypotheses and imported valuation input

Work over C. Let delta=-v be a normalized divisorial semidegree on C(x,y)
with delta(x)=omega_0>0. Use the **full** key-form chain
g_0=x,g_1=y,...,g_(n+1)=G in C[x,x^-1,y], with weights omega_i.
Put d_i=gcd(omega_0,...,omega_i), alpha_i=d_(i-1)/d_i. For 1<=j<=n,

    g_(j+1)=g_j^alpha_j - theta_j x^beta_(j,0)
                              product_(1<=i<j) g_i^beta_(j,i),
    alpha_j omega_j=beta_(j,0) omega_0 + sum_(1<=i<j) beta_(j,i) omega_i,
    0<=beta_(j,i)<alpha_i (i>=1), theta_j in C*,
    omega_(j+1)<alpha_j omega_j, d_(n+1)=1.

Here beta_(j,0) may be negative. No primitive-compactification or positive
terminal-weight assumption is made. The empty recurrence case n=0 is
included. Key forms and their generic descending Puiseux description are
imported from Mondal's [part II, Definition2.8 through Algorithm2.12 and
AppendixA](https://arxiv.org/html/1307.5577v4). In the generic expansion,
the leading coefficients of x,g_1,...,g_n are nonzero constants, while
the leading coefficient of G is a nonconstant polynomial h(xi). Thus h(xi)
is transcendental over C. These facts distinguish a terminal key from an
intermediate polynomial whose leading coefficient is still constant.

The separately published Proposition2.11 of [Mondal, MMJ66(2017),227–244](https://pinakimondal.org/papers/mmj.pdf),
PDF10/printed236, identifies “terminal key polynomial” with “all key forms
polynomial”, for delta != ordinary degree and delta(x)>0. Proposition2.10
immediately preceding it supplies the general two-boundary model; the
equivalence has no positive-terminal-weight hypothesis. Ordinary degree
is handled directly below. These classification statements are imported,
not independently reproved. The polynomial-chain algebra below is written
out independently, including all signs. Source pins and the exact fresh
page/section audit are in [the manifest](sources/polynomial-key-chain/manifest.json).

## 2. Polynomial recurrences, signs and a basis

Assume every g_i is polynomial. Inductively each is monic in y, with

    D_1=1, D_(i+1)=alpha_i D_i.

Indeed the subtracted monomial has y-degree
sum_(k<i) beta_(i,k)D_k <= D_i-1, by mixed-radix expansion; this also
works when some alpha_i=1. It follows that beta_(j,0)>=0: otherwise the
product of earlier monic polynomials, none divisible by x, would give a
genuine pole at x=0 that g_j^alpha_j cannot cancel. Conversely, if every
beta_(j,0)>=0, induction makes every key polynomial. This is precisely the
polynomial (also called algebraic) condition on the signed numerical key
sequence, without asserting that it defines an algebraic compactification.

All nonterminal weights omega_0,...,omega_n are nonnegative. If omega_j
were the first negative nonterminal weight, its recurrence would express
alpha_j omega_j as a nonnegative combination of earlier nonnegative
weights, a contradiction. A zero nonterminal weight has alpha_j=1,
since adjoining zero does not change the gcd. Its allowed basis exponent
is therefore zero.

The monomials

    x^b_0 g_1^b_1 ... g_n^b_n G^k,
    b_0,k>=0, 0<=b_i<alpha_i for 1<=i<=n                         (B)

form a C-basis of C[x,y]. To prove existence, divide a y-degree uniquely
into bounded mixed-radix digits b_1,...,b_n and the unrestricted last
digit k. The corresponding product is monic of that degree; subtract its
leading C[x]-coefficient times this product and induct on y-degree. For
uniqueness, a nonzero relation over C[x] has a highest y-degree, belonging
to a unique digit tuple, so its leading coefficient cannot vanish.
Expanding the C[x] coefficients gives (B).

For this basis, delta of a polynomial is **exactly the maximum** of its
nonzero monomial weights. To verify no cancellation, fix a weight and k.
If two bounded digit vectors have that weight, take their largest
differing index i. Their weight equation modulo d_(i-1) implies alpha_i
divides the difference of the i-th digits, so they are equal. Descend
through the digits; then omega_0>0 forces equality of b_0 too. Thus for
fixed k there is at most one term at that weight, with a nonzero constant
multiple of h(xi)^k as leading coefficient. Distinct k cannot cancel
because h(xi) is transcendental. This reproduces the needed polynomial
case of part II LemmaA.2 rather than assuming that values of sums always
equal maxima.

## 3. Regular algebra and the exact Jacobian obstruction

Let A_v=C[x,y] intersect O_v, equivalently {F:delta(F)<=0}.

* If omega_(n+1)>0, strict descent and the preceding sign argument force
  every earlier weight to be positive. Hence A_v=C.
* If omega_(n+1)=0, every earlier weight is positive. Hence A_v=C[G].
* If omega_(n+1)<0, any nonconstant monomial of (B) having k=0 has
  positive weight: zero-weight earlier keys have allowed exponent zero.
  Consequently every F in A_v has the form c+GH with H polynomial.
  The image of A_v in C[x,y]/(G) is exactly C. This is an inclusion
  A_v subset C+G*C[x,y], not an equality of these algebras.

In the first two cases any pair P,Q in A_v has zero Jacobian. In the
third, write P=p_0+GH and Q=q_0+GK. Direct differentiation gives

    Jac(P,Q)=G*( H Jac(G,K)+K Jac(H,G)+G Jac(H,K) ).                (J)

G is nonconstant because it is monic of positive y-degree. Therefore
Jac(P,Q) cannot be a nonzero constant. This is an unbounded identity,
not an inference from the finite controls.

Thus a plane Keller pair P,Q cannot both be regular at a divisorial
valuation with an entirely polynomial key chain. In particular every
actual omitted prime in its finite normalization must have a Laurent
key form, regardless of field degree or all other omitted primes.
By the imported terminal/all-polynomial equivalence, its **terminal**
key is Laurent. Ordinary degree has A_v=C and is excluded directly.

There is a stronger geometric statement. Let N be any normal affine
surface containing the whole marked U=A2 as a dense open, and let D
be an actual omitted prime. Some source polynomial coordinate has
positive delta: otherwise x,y regular at v_D give a second centre in U
for the same valuation, contradicting separatedness of N. Swap coordinates
if necessary. Restriction embeds O(N) into A_v. If A_v=C or C[G], the
affine coordinate algebra cannot have transcendence degree two. In the
negative case (J), every global function is constant on the curve
G=0 contained in U. That curve has distinct complex points, which affine
global functions must separate. Contradiction. This rules out polynomial
terminal keys even without a Keller map or any smoothness assumption.
It applies in every polynomial coordinate system with delta(x)>0.

## 4. Integer filter and the first necessary Laurent step

The [integer helper](key_sequence_obstruction.py) computes the bounded
beta digits by descending gcd congruences. At step i solve

    beta_i = (remaining/d_i)*(omega_i/d_i)^(-1) mod alpha_i,

then subtract beta_i omega_i; modulus one has digit zero. The last
remainder divided by omega_0 is beta_0. This uses integers and no
expanded polynomial. A negative beta_0 is necessary for an actual
affine/Keller omitted prime. Passing this filter certifies neither
affine realization, finiteness, a target algebra nor a Keller map.

If a polynomial prefix first reaches negative weight at a nonterminal
g_j, its next beta_0 must be negative: the earlier weights are
nonnegative and the target alpha_j omega_j is negative. Its leading
coefficient is still constant; treating g_j as terminal would give a
false exclusion. A polynomial zero-weight prefix has alpha_j=1, so
the next form is g_j-theta, with negative weight; the following step,
if present, is Laurent. A first Laurent step can instead occur while
all preceding weights are positive, through a numerical semigroup gap.
For example weights (3,2,1,-2) have g_2=y^3-x^2 and
g_3=g_2-c*x^-1*y^2, c!=0. Here beta_0=-1 although omega_2=1.
This is also part II Example2.13; no negative-prefix-only rule is used.

The excluded family of every finite depth L>=1 has weights

    (2^L, 3*2^(L-1), 3*2^(L-2), ..., 3, -1),
    g_2=y^2-x^3, g_(j+1)=g_j^2-g_(j-1) for 2<=j<=L.

All recurrence alpha_j=2 and all beta_0>=0. Its terminal y-degree is
2^L; at L>=2 this exceeds the PF99 minimal key degree two. Its canonical
index, from the retained general key formula, is 2^(L+1)-2. The new
theorem excludes the entire family at all depths, not just the controls.
The checker evaluates actual generic Laurent charts only for L<=3;
the L=64 and L=1000 controls use integer circuits without expansion.

The retained actual smooth affine control f^2 V=W(W-1),
W=f^2 x, V=x(f^2 x-1), is preserved. On its omitted prime W=1,f=0,
x=u^-2+t, f=u. In coordinates X=x,Y=x+f the full chain is

    X,Y,Y-X,(Y-X)^2-X^-1, with weights (2,2,-1,-4).

The intermediate polynomial Y-X has positive valuation but constant
leading coefficient 1; the terminal key is Laurent with nonconstant
leading coefficient t. Thus the example passes the new necessary
filter and is not overexcluded. In coordinates (x,f) its weights are
(2,-1,-4). No Keller realization is asserted for this affine example.

## 5. Failed cone route: do not import an unrestricted criterion

The retained published MMJ PDF literally states in Corollary1.10(1)
that the real cone of (ordinary degree,delta) is closed iff all/last
key forms are polynomial. Its Theorem1.8(2) similarly states an
unqualified infimum-attainment equivalence. The unrestricted negative
terminal-weight inference fails an exact control; it is not used here.
This concerns the scope of those specific clauses, not the separately
used Proposition2.11 or nonnegative-sign criterion Theorem1.4(1).

For the actual affine example x=u^-2+t,f=u and any polynomial F of
ordinary total degree d, every term of x^i f^j expands with exponents
j-2i+2k, 0<=k<=i, i+j<=d. All exponents lie in [-2d,d]. The substitution
is injective (u=f and t=x-f^-2), so cancellation leaves a nonzero
Laurent polynomial and cannot create exponents outside that interval.
Consequently

    -d <= delta(F) <= 2d.

Both edges occur: (deg f,delta f)=(1,-1), (deg x,delta x)=(1,2).
The real convex cone is therefore exactly the **closed** wedge generated
by these vectors, despite Laurent terminal key f^2-x^-1 of weight -4.
The degree-normalized infimum is attained by f as well.

The paper's own Example3.1, PDF15–16, also explicitly computes an
attained normalized infimum for keys x,y,y-x^-1 with Laurent terminal.
Thus no inference “Laurent terminal implies nonclosed cone” or
“non-finitely-generated degree/value semigroup” is registered. No
corrected universal version of those clauses or claim about an erratum
is asserted. The example and source-scope issue are retained so this
failed route is not restarted.

## 6. Verification and remaining work

Run `ulimit -v 524288; .venv/bin/python research/jacobian/affine-complements/verify_polynomial_key_chain.py`.
The [receipt](polynomial_key_chain_receipt.json) pins this proof, the helper,
the primary sources and PF99 predecessor. Controls include independent
monic divisions, actual generic Laurent leading coefficients, three
regular-algebra cases, zero-weight digits, exact Jacobian factorization,
the positive semigroup gap, the actual affine Laurent escape and the
closed-cone countercertificate. Controls check implementation and signs;
the written arguments establish the universal scope. No formal proof
assistant or external review is claimed.

The remaining all-degree frontier is genuinely Laurent terminal keys,
their cancellations, canonical budget and compatibility with the actual
finite target algebra; singular and higher-index boundary packets remain.
PF99 conditional whole-chart results and PF98 geometric peripheral data
are retained. The recursive marked-linking/dihedral comparison direction
suggested by the user remains available if the algebraic route stalls;
the sealed PF92 recursion is not extended beyond its longitude-one and
integral-core hypotheses here. Do not restart deeper coefficient sweeps
or quartic realization campaigns.
