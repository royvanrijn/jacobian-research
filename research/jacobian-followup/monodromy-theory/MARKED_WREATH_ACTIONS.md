# Marked-root wreath actions, arithmetic constants and atomicity

Work in characteristic zero. The explicit family and geometric group
below reproduce Harish's construction and Fonollosa's
[Theorem1.2 and Proposition2.2](https://arxiv.org/html/2609.33491v1).
The [source pin](sources/marked-wreath-manifest.json) reuses the retained
primary HTML; no novelty of that family or geometric computation is claimed.
The arithmetic calculation, complete block interval and consequences for
atomicity and deck groups are derived here.

[Generator](marked_wreath.py) · [checker](verify_marked_wreath.py) ·
[model](marked_wreath_model.json) · [receipt](marked_wreath_receipt.json).

## Statements

**KMFI69 (actual marked-root action family, reproduction).** For every
integer n>=3, put m=n-2. There is an explicit polynomial Keller self-map
of ordinary affine n-space, defined over Q, with determinant m/2,
generic degree N=nm and geometric permutation group on pairs
(i,j) in {1,...,n} times Z/m:
\[
 G_n=\begin{cases}
 (C_m)^n\rtimes S_n,&m\text{ odd},\\
 \{(v,\sigma):\sum_i v_i\equiv \operatorname{par}(\sigma)\pmod2\},
       &m\text{ even},
 \end{cases}                                                   \tag{1}
\]
where the second group is a subgroup of the same wreath product.
Its order is m^n n!/gcd(2,m). The even case is stated as a subgroup,
without an unjustified splitting assertion.

**KMFI70 (arithmetic group and exact constants).** The algebraic
constant field of its normal closure over Q is exactly Q(zeta_m).
The arithmetic group acts by
\[
 (i,j)\longmapsto(\sigma(i),u j+v_{\sigma(i)}),\qquad
 u\in(\mathbb Z/m)^\times.                                     \tag{2}
\]
For odd m there is no further condition. For even m it satisfies
\[
 \sum_i v_i\equiv\operatorname{par}(\sigma)+\chi_n(u)\pmod2,
 \quad
 \chi_n(u)=\begin{cases}
 0,&n\equiv0\pmod4,\\
 (u-1)/2\pmod2,&n\equiv2\pmod4.
 \end{cases}                                                   \tag{3}
\]
In the second case 4 divides m, so the character is well-defined.
The order is |G_n|*phi(m). Over an arbitrary characteristic-zero
coefficient field k, replace the unit group by the cyclotomic image
Gal(k(mu_m)/k); the constants are exactly k(mu_m).
For m=1 all marking groups and cyclotomic quotients are trivial.

**KMFI71 (complete blocks and stable atomicity).** Over algebraically
closed constants, the entire proper intermediate-field interval is
\[
 E_d=K(t,a^d),\qquad d\mid m,\qquad [L:E_d]=d,                  \tag{4}
\]
where L=K(t,a), a^m=P'(t). Each L/E_d is cyclic. The proper
blocks containing one chosen sheet have exactly the sizes d|m;
they lie inside its root block. Together with K these are all
intermediates. Every displayed Keller map is absolutely and stably
atomic, despite this divisor lattice of blocks.

**KMFI72 (generic deck groups of arbitrary finite abelian type).**
The geometric generic deck group of the map in KMFI69 is C_m,
while its polynomial source deck group is trivial by KMFI47.
Over k the generic deck group is mu_m(k); over Q its order is gcd(m,2).
Consequently every finite cyclic group is the geometric generic deck
group of an actual absolutely stably atomic Keller map. Every finite
abelian group occurs as geometric generic deck group of an actual
Keller map by taking independent products. These are deck groups,
not the corresponding monodromy groups.

None of these assertions decides exceptional degree11 PSL2(11), a
nonlinear specialization, or the full action spectrum.

## 1. An ordinary polynomial source, with its entire Jacobian

Use coordinates x,y,z3,...,zn and write p=1+xy. Define
\[
 \Phi(T)=(1-yT)^{n-1}(pT-x)+(pT-x)^2
 \left[ny(yT+\tfrac12)+z_3(pT+\tfrac x2)
                   +\sum_{j=4}^n z_jT^{j-2}\right].
                                                               \tag{5}
\]
Its linear coefficient is1. Write
\[
 \Phi(T)=\sum_{j=3}^n C_jT^j+A_2T^2+T-A_0.
\]
The polynomial map is F=(C3,...,Cn,A2,A0). In particular no
rational chart is being substituted for its source.

Here is a determinant calculation, also replayed on complete expanded
coordinates in the finite controls. On the dense open xp!=0 let
t=x/p and a=1/p. The high-coefficient derivative block with respect to
z3,...,zn is upper triangular, with diagonal p^3,p^2,...,p^2.
Also det d(t,a)/d(x,y)=-t*a^2. Holding the high coefficients fixed,
\[
 A_2=\frac{a^m-1-\sum_{j=3}^n jC_jt^{j-1}}{2t},
 \qquad A_0=\sum_{j=3}^n C_jt^j+A_2t^2+t.
\]
It follows directly that det d(A2,A0)/d(t,a)
=-m*a^(2m-1)/(2t). Their product is m/2 because ap=1
and 2m+1=2n-3. Polynomiality extends the identity across the omitted
dense-open complement. This proves the global Keller assertion,
including x=0 and p=0, rather than only a rational volume identity.

The same substitutions give Phi(t)=0, Phi'(t)=a^m.
Conversely x=t/a, y=(1-a)/t and upper-triangular back substitution
recover the source coordinates from t,a and the target coefficients.
Thus the generic source field is exactly L=K(t,a).

We use the original n-dimensional map. Fonollosa additionally proves
polynomial coordinate descents to dimension3; no reproduction of
that separate recursive descent is needed for the dimension-varying
occurrence claim registered here.

## 2. Root and marking fields

The target coefficients are independent by the constant nonzero
Jacobian. Set q=Cn. If t1,...,tn are the roots, the coefficient of T
imposes
\[
 q(-1)^{n-1}e_{n-1}(t_1,\ldots,t_n)=1.
\]
The generic root splitting field can therefore be written as the
purely transcendental field k(t1,...,tn), with q given by this equation.
The permutation-invariant subfield is the target field: its elementary
symmetric generators recover every coefficient. Hence its root group
is S_n over both closed and arbitrary characteristic-zero constants.

Let d_i=P'(t_i)=q*product_(j!=i)(t_i-t_j).
At the divisor t_i=t_j, q and all other root differences are units
generically, and d_i,d_j have valuation1, the others0. If
product d_i^(r_i) is an mth power in the splitting field, then
r_i+r_j=0 modm for every pair. Comparing three indices gives all
r_i equal and 2r_i=0. Thus the only possible nonzero relation is
(m/2,...,m/2), when m is even. With
delta=product_(i<j)(t_i-t_j),
\[
 \prod_i d_i=(-1)^{n(n-1)/2}q^n\delta^2.                       \tag{6}
\]
Over algebraically closed constants the putative even relation
indeed is an mth power. Kummer theory gives the base rotation group
of size m^n/gcd(2,m). For even m, choosing a_i^m=d_i shows
\[
             \prod_i a_i^{m/2}=c\,q^{n/2}\delta,\qquad
             c^2=(-1)^{n(n-1)/2}.                            \tag{7}
\]
The fixed scalar c imposes exactly the parity condition in(1);
the order calculation proves equality. Independent rotations, or
opposite rotations in two root blocks in the even case, make the
action transitive. The conjugates of (t,a) generate every a_i and
every root, so this is the normal closure and N=nm.

## 3. Arithmetic: the cyclotomic twist must be retained

The normal closure contains mu_m: the quotient of two markings
of the same root is a primitive mth root of unity. This is a
lower bound on its constants, not yet an equality.

For odd m all permutations in(2) are an upper bound for arithmetic
monodromy. For even m, c in(7) can be chosen as1 when n=0 mod4
and i when n=2 mod4. In the latter case m=0 mod4, so i already
belongs to k(mu_m). Applying a constant automorphism u to(7)
multiplies c by c^(u-1); permuting roots multiplies delta by
sgn(sigma), and rotating markings multiplies the left side by
(-1)^(sum v_i). This gives exactly(3).

For each cyclotomic automorphism there are precisely |G_n| possible
elements satisfying this upper bound. On the other hand the
geometric group is G_n and the constants already contain k(mu_m),
so the arithmetic order is at least
|G_n|*[k(mu_m):k]. The bounds coincide. This proves both the full
arithmetic group and the absence of any additional constants.

In particular no extra quadratic constant is left over from(7).
If c=i, the necessary fourth roots of unity are already present.
For n=6,m=4, pure conjugation of markings with roots fixed violates(3);
a lift with no color translation must exchange two root labels.
Ignoring this twist gives a wrong embedded arithmetic action.

For reference:
\[
\begin{array}{c|r|r|r|c}
 n&N&|G_{\rm geom}|&|G_{\rm arith}/\mathbb Q|&
                    \text{constants}\\\hline
 4&8&192&192&\mathbb Q\\
 5&15&29160&58320&\mathbb Q(\zeta_3)\\
 6&24&1474560&2949120&\mathbb Q(i)\\
 7&35&393750000&1575000000&\mathbb Q(\zeta_5)
\end{array}
\]

## 4. All blocks, with no implicit polynomial quotient

Fix alpha=(1,0), and let H be its stabilizer. H fixes every point
of root block1, and is transitive on all sheets outside that block.
Indeed its S_(n-1) quotient moves the other root labels freely.
It can rotate any one such block by any desired color; in the
even case use a different nonfirst root block to correct the
required parity. There are enough blocks since n>=3.

A block containing alpha is H-invariant. If it meets another root
block, it therefore contains all m(n-1) outside sheets. Its size
exceeds N/2, so it is the whole action. Every proper block is thus
contained in root block1.

The setwise stabilizer of root block1 induces its full regular
C_m translation group. Its blocks are the cosets of the unique
subgroups C_d for d|m, and these partitions are preserved by every
root permutation and translation in G_n. This proves completeness,
including the lack of blocks mixing different root labels.
Galois correspondence gives exactly(4), and K(t) is the unique
minimal proper intermediate above K when m>1.

Now suppose F=Q o R is a polynomial factorization into equal-dimensional
affine spaces over closed constants. The chain rule forces both factors
Keller. If both generic degrees exceed1, the intermediate field must
be some E_d with d>1, but its upper extension L/E_d is cyclic.
KMFI1--2 exclude every nontrivial Galois Keller extension. This is a
contradiction. A degree-one Keller factor is an automorphism.
The group action and field interval persist under independent
stabilization, proving absolute and stable atomicity.

There is no assertion that an intermediate E_d/K itself has, or lacks,
a separate Keller model. The obstruction here concerns the upper
extension L/E_d, which is the one a polynomial inner factor would
have to realize.

## 5. Deck symmetry and independent products

Any permutation commuting with a transitive group is determined by
its value at alpha, which must be fixed by H. There are exactly m
such fixed points here. Diagonal color translation supplies m
commuting permutations, so the geometric generic deck group is
exactly C_m.

There is also a direct source-coordinate explanation for the missing
polynomial symmetry. A generic deck rotation a->zeta*a fixes t and gives
\[
       x\longmapsto x/\zeta,\qquad
       y\longmapsto y+(1-\zeta)/x.
\]
For zeta!=1 the second coordinate has a genuine pole along the
ordinary source divisor x=0. Thus the full geometric generic deck
group has only its identity extending regularly. No boundary
compactification or conjecture is needed for this family-specific check.

Under arithmetic monodromy, diagonal translation by b commutes with
the cyclotomic part precisely when (u-1)b=0 modm for all cyclotomic u.
These are exactly mu_m(k). For Q, conjugation u=-1 gives at most
gcd(2,m) choices, and those choices are fixed by every odd unit
when m is even; for odd m only the identity remains.

Given a finite abelian group A=product C_(m_i), take independent
products of the maps with n_i=m_i+2, omitting trivial factors.
The geometric closure fields with independent target variables
are linearly disjoint and have product group and product sheet
stabilizer. Consequently
\[
 N_{\prod G_i}(\prod H_i)/\prod H_i
       =\prod N_{G_i}(H_i)/H_i=A.
\]
The product is a polynomial Keller self-map of ordinary affine
space, with constant determinant product(m_i/2).
One can work over a number field containing all mu_(m_i) to make
the entire generic deck group defined over the coefficient field.
KMFI47 still makes every polynomial fixed-target source deck
endomorphism the identity. A nontrivial product need not be atomic;
the cyclic individual factors are the atomic witnesses.

Finally, G_n contains an N-cycle: take an n-cycle on the root
blocks with total color translation1. For even n both its block
parity and total color parity are odd, so it satisfies(1).
Hence the cycle-gcd lcm invariant of KMFI10 is N. This is a
positive boundary control only: it does not prove that the exact
field pair is independently curve-induced.

The checker gives bounded exact polynomial Jacobians, permutation
orders, every principal block, arithmetic parity controls and pair-
valuation Smith forms. The arguments above supply the unbounded
statements. Formal verification, independent whole-proof review and
completion of OP-KMON are not claimed.
