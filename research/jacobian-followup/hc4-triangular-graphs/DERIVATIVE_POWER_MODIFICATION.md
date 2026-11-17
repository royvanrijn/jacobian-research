# Leading poles, a cyclic derivative surface, and normal rescaling

This note proves three further restrictions for the output pair
`(T,xy-P(T,z))`. The arguments allow every coefficient, source and embedding
degree. They do not solve either plane endpoint or unrestricted HC4/JC2.
The mathematical claims are `HC4MYIE0/1/2` in
[MATH_STATUS.json](../../MATH_STATUS.json).

Let `B=C[s,t]`, `L=Frac(B)`, let `K=(f,g)` have constant nonzero Jacobian
`j`, and write

```
P(T,z)=sum_(i=0)^m b_i(T) z^i,   a=b_m != 0.
```

The retained [normalized-label proof](LEADING_VALUATIONS_AND_CRITICAL_CURVES.md)
and [square-free pullback proof](ABELIAN_CLASSGROUP_OUTPUT_OBSTRUCTION.md)
give, for any root `h in L` of `P(f,h)+g=0`,

```
u=a(f)h in B,
Q(T,U,V)=U^m+sum_(i=1)^(m-1) b_i(T)a(T)^(m-1-i)U^i
                    +a(T)^(m-1)(b_0(T)+V),
R=Q_U(f,u),
J(f,u) R=-a(f)^(m-1) j.                         (1)
```

The polynomials `f,u` are algebraically independent: otherwise the root
equation would make `g` algebraic over `C(f)`. Every `f-c` is nonconstant
and square-free. Each of its prime divisors maps dominantly to the target
line `T=c`; in particular `b_0(f)+g` has valuation zero there. The rational
root theorem also gives that the denominator of `h` divides `a(f)`.

## 1. Every rational root needs an intermediate leading pole (HC4MYIE0)

For `m>=3`, a rational root must have a source prime divisor over some
leading zero `c` with

```
q=ord_c(a)>=2,      0 < -ord_delta(h) < q.       (2)
```

Thus a root with only regular branches and full-order leading poles is
impossible. This is a statement about actual source valuations, without
assuming irreducibility of any special `f`-fibre.

### 1.1 The derivative power forced by absence of intermediate poles

Suppose (2) never occurs, and set `k=m-2`. At a prime `delta` over `c`,
the valuation of `a(f)` is exactly `q`, since `f-c` is square-free.

If `h` is regular at `delta`, both `P_z(f,h)` and `J(f,h)` belong to the
local ring `B_(delta)`. Applying `J(f,-)` to the rational root equation
gives

```
P_z(f,h) J(f,h)=-j.
```

Both factors consequently have valuation zero. The identity
`R=a(f)^k P_z(f,h)` gives `ord_delta(R)=kq`. Also
`ord_delta(u)>=q`.

If `h` has a pole, its order is at most `q`, hence is exactly `q` by
assumption. The only possible cancellation of the leading term's valuation
`-(m-1)q` is with `b_(m-1)(f)h^(m-1)`. Therefore
`b_(m-1)(c)!=0`, and `u` specializes to `-b_(m-1)(c)` in the residue field.
At `T=c` the normalized polynomial is
`U^(m-1)(U+b_(m-1)(c))`, whose derivative is nonzero at that simple root.
Thus `ord_delta(R)=0`.

There are no prime factors of `R` outside `a(f)` by (1). Since `B` is a
UFD and every complex nonzero constant has a `k`th root, there is a
polynomial `D` with

```
D^k=R,       p=a(f)/D in B,       w=u/D in B.   (3)
```

Indeed the prime exponents of `D` are `q` on regular branches, zero on
full-pole branches. The numerator of each quotient in (3) has at least
these exponents. The case `k=1` is included. If `a` is constant, (1)
already contradicts algebraic independence, since `R` has positive degree
in `u`. Hence assume `a` nonconstant.

### 1.2 A polynomial map from the original plane

In target coordinates `(T,D,W,p)` define the surface equations

```
a(T)=D p,
m D W^(m-1)+sum_(i=1)^(m-1)
       i b_i(T) p^(m-1-i) W^(i-1)=1.           (4)
```

Equations (1) and (3), after cancellation in the source domain, give the
polynomial map

```
Phi:A2 -> A4,       (s,t) |-> (f,D,w,p),
```

with image in (4). The original root equation gives a further identity
on this image:

```
D w^m+sum_(i=1)^(m-1) b_i(f)p^(m-1-i)w^i
              +p^(m-1)(b_0(f)+g)=0.           (5)
```

There is no extension or replacement of the source plane here.

The open part `D!=0` of (4), using `u=DW` and `p=a(T)/D`, is isomorphic to

```
D^k=Q_U(T,u),       D!=0.                     (6)
```

This is integral. Over `C(T,u)`, the right side has `u` degree `m-1=k+1`.
At `u=infinity`, every root of `X^k-Q_U(T,u)` has valuation `-(k+1)/k`.
For a monic factor of degree `e` over `C(T,u)`, the valuation of its
nonzero constant coefficient would be `-e(k+1)/k`, an integer. Coprimality
of `k,k+1` forces `k|e`, proving irreducibility (with `k=1` immediate).
Gauss lemma for the monic polynomial proves the claimed domain statement.

Let `V_0` be the closure of this integral open surface. The image of `Phi`
has dimension two because `f,u=Dw` are algebraically independent; its
generic point lies in `D!=0`. Hence its closure is exactly `V_0`, and
`Phi` is generically finite onto `V_0`.

### 1.3 The omitted curve is C*

Choose any zero `c` of `a` and put `b=b_(m-1)(c)`. In (4) the locus
`T=c,p=0` is the curve

```
Gamma:  m D W^(m-1)+(m-1)b W^k=1,
D=[1-(m-1)b W^k]/[m W^(m-1)].                (7)
```

It is a closed curve isomorphic to `C*`, with coordinate `W` invertible.
Its dense open subset `D!=0` belongs to (6), so the entire curve belongs
to `V_0`. If a source point maps to this curve, (5) gives `D W=-b`.
Combining this with (7) gives

```
-b W^k=1.                                    (8)
```

If `b=0` no point is hit; otherwise at most `k` target points are hit.
Thus `Phi` omits a dense open of `Gamma`.

We explicitly reuse Jelonek--Lason, arXiv:1411.5011v2, Theorem 3.2:
the nonproperness set of a generically finite polynomial map from
complex affine space to affine space is covered by nonconstant polynomial
images of `A1`. The theorem applies to `Phi:A2->A4`; normality of `V_0`
is not required. Only its qualitative conclusion is used.

For completeness, its use really isolates `Gamma` as a component. A dense
principal open of `V_0` has finite preimage: clear denominators in monic
equations for the source coordinates over `Frac(O(V_0))`. This shows the
nonproperness set is contained in a proper closed subset of the surface,
of dimension at most one. Omitted points of `Gamma` are nonproper: every
neighborhood meets the dense image; a sequence of their preimages could
not stay bounded without producing a preimage of the omitted point.
Closedness of the nonproperness set therefore puts `Gamma` inside it.
Choose a point of `Gamma` outside its intersections with any other curve
components. A polynomial curve through this point supplied by the
imported theorem must lie in `Gamma`.

But along any polynomial `A1` map to (7), `W` is an invertible polynomial
and hence a nonzero constant, and then `D` is constant too. Such a map
cannot be nonconstant. This contradiction proves (2).

### 1.4 A coefficient-only pole test

At a leading zero `c`, put `q=ord_c(a)` and `n_i=ord_c(b_i)` for
`1<=i<m`, with `n_i=+infinity` if `b_i=0`, and `n_m=q`.
Define the finite pole spectrum

```
E_c={r in Z: 1<=r<=q and
          min({0} union {n_i-i*r:1<=i<=m}) occurs at least twice}.
```

The zero entry is the valuation of `b_0(f)+g`. A source pole order must
lie in `E_c`, since a unique term of least valuation cannot cancel. This
is a direct valuation argument, without a Newton polygon theorem import.

If `E_c` is contained in `{q}` for every leading zero, (2) is impossible,
so there is no rational root. A convenient sufficient condition is

```
gcd(a,a',b_(m-1))=1.                         (9)
```

At a simple leading zero the denominator bound already leaves only order
one. At a multiple zero, (9) gives `b_(m-1)(c)!=0`; if `0<r<q`, the
`(m-1)`st term has unique minimum, excluding that order. This covers
every square-free leading coefficient of arbitrary parameter degree,
and every leading coefficient coprime to its next coefficient, in every
`m>=3`. A further example is `T^3 z^3+z`: its only candidate nonzero
slope is `3/2`, so its integer pole spectrum is empty.

Absence of rational roots does **not** imply irreducibility for `m>=4`.
Factors of degrees `2+2`, for example, are outside this conclusion.

## 2. Normal monic rescaling gives full irreducibility (HC4MYIE1)

Let `m>=2`, `A(T)!=0`, and let `c_1,...,c_(m-1),d` be arbitrary
polynomials. Consider

```
P(T,z)=A(T)^(m-1) z^m
       +sum_(i=1)^(m-1) A(T)^(i-1)c_i(T)z^i+d(T).
F(T,U)=U^m+sum_(i=1)^(m-1)c_i(T)U^i.
```

Assume that at every multiple zero `c` of `A`,

```
gcd_U(F(c,U),F_U(c,U),F_T(c,U))=1.            (10)
```

This is vacuous when `A` is square-free. Then `P(f,Z)+g` is irreducible
over `L`, with no factor-degree qualification.

To prove this, substitute `U=A(T)z` and multiply by `A(T)`. The equation
becomes the monic cover

```
S(T,U,V)=F(T,U)+A(T)(V+d(T))=0.              (11)
```

It is an integral surface. As a polynomial linear in `V`, it is primitive
over `C[T,U]`: `gcd(A,F)=1` because `F` is monic in `U`. It is irreducible
over `C(T,U)`, and Gauss lemma proves integrality. It is finite monogenic
over `C[T,V]`, since it is monic of degree `m` in `U`.

Its singular points satisfy

```
A=0, F=0, F_U=0, F_T+A'(V+d)=0.             (12)
```

At a simple zero of `A`, there are finitely many choices of `U` (as `F`
is monic), and the last equation determines at most one `V` for each.
At a multiple zero, (10) excludes all singular points. Thus the singular
locus is finite. This hypersurface is Cohen--Macaulay by Stacks `00SB`,
hence satisfies `S_2`. The hypersurface derivation criterion, Stacks
`07PF`, and the finite common-zero locus of the partial derivatives give
regularity in codimension one. Serre's
criterion, Stacks `031S`, proves normality. These are explicit proof
imports, recorded in the [source manifest](sources/derivative-modification-import.json).

Apply the retained [normal monogenic theorem HC4MYHY0](NORMAL_MONOGENIC_COVERS.md)
to (11). It remains irreducible after Keller base change. The invertible
linear change `U=A(f)Z` over `L` carries it to
`A(f)(P(f,Z)+g)`. Therefore the latter equation is irreducible too.
No arbitrary nonnormal primitive model is treated as normal, and no
normalization is asserted globally monogenic without (11).

For example, every

```
P=T^2 z^3+T B(T)z^2+C(T)z+D(T)
```

closes, with completely arbitrary `B,C,D`. This includes `T^2 z^3+z`,
whose integer intermediate pole order one is permitted by the first
test. The two methods address complementary local behavior.

The repeated-root condition (10) matters. With `A=T^2` and `F=U^3`,
`U^3+T^2 V=0` is singular along `T=U=0`. This chosen model is nonnormal
and cannot invoke HC4MYHY0. Another retained gate can still close a
particular output; failure of (10) does not assert a counterexample.

## 3. Cubic and embedded-output consequences (HC4MYIE2)

For `P=a z^3+b z^2+c z+d`, reducibility over `L` would give a rational
root. Consequently `P(f,Z)+g` is irreducible whenever every integer pole
spectrum has no intermediate order, in particular whenever
`gcd(a,a',b)=1`. This closes **every square-free leading coefficient**,
with arbitrarily many leading zeros and no condition on `b^2-3ac`.
It also closes arbitrary leading multiplicities when `gcd(a,b)=1`.

Together with Section 2 and the retained constant-leading result, this
closes every cubic with `deg_T(a)<=2`, with all other coefficient degrees
unrestricted. The only extra case is
`a=alpha(T-c_0)^2`. If `b(c_0)!=0`, (9) applies. If `b(c_0)=0`, put
`A=sqrt(alpha)(T-c_0)`; then `b/A` is polynomial and Section 2 applies
with square-free `A`.

For every such irreducible cubic, and for every degree in the fully
irreducible family of Section 2, each closed stable embedding

```
A^(r+2) -> A^(r+3),    r>=1,
```

carrying a plane Keller pair through `(T,xy-P)` forces that pair to be
an automorphism. The retained [HC4MYHZ2](LEADING_FIBRE_UNIT_GATE.md)
constructs a common filtered basis, proves geometric integrality of every
affine height-one fibre, identifies the generic class lattice
`Z^m/diagonal`, and uses root transitivity to make every component
principal. Its global pencil argument then gives the inverse. All its
conditions apply to arbitrary polynomial `P` of positive `z` degree.
These output consequences use that existing proof, without a new
source Hessian invariance assertion.

The automorphic embedding control remains available for every `P`:
`(s,t,w)->(s,1,t+P(s,w),w)`. With the unused fourth source variable `v`,
the actual pullback of `Tv+yz+x^2` is
`1+sv+tw+wP(s,w)`, of Hessian determinant one. Its gradient inverse is
`(s,t,w,v)=(d,c-P(d,b)-bP_z(d,b),b,a-bP_T(d,b))` for gradient outputs
`(a,b,c,d)`.

## 4. Precise remaining boundary and replay

A surviving cubic rational root must have an actual intermediate source
pole over a multiple leading zero, necessarily a common zero of
`a,a',b`. It must also escape the normal rescaling criterion and all
earlier valuation, monodromy and critical-curve gates. The square-free
linear nonsquare-critical branch and the many-leading-level square-free
branch formerly left open are retired by this note. No classification
of all repeated-leading cubics is claimed. For example
`T^3 z^3+T z^2+z` passes the simple coefficient pole test with order one;
this observation alone asserts neither survival of all retained gates
nor any Keller realization.

The non-Keller control `f=s(1+st), g=-t(1+st), h=1/s` solves
`f^2 h^3-f h^2+g=0` and exhibits an intermediate pole. Its Jacobian
is nonconstant. It checks why local pole behavior must be analyzed;
it is not a Keller counterexample.

Higher-degree factorizations without linear factors, uncovered output
geometries, the singleton-zero endpoint, the unit-ideal endpoint and
unrestricted HC4/JC2 remain open. In particular no rational-root gate
here is a full arbitrary-degree irreducibility theorem.

Run [the exact replay](verify_derivative_modification.py) for polynomial
identities, pole-spectrum controls, the excluded singular-line model,
and retained byte pins:

```sh
timeout 60s .venv/bin/python research/jacobian-followup/hc4-triangular-graphs/verify_derivative_modification.py
```

The [receipt](receipt-derivative-modification.json) checks finite controls;
the all-degree proofs are the arguments above. The
[validation](validation-derivative-modification.json) distinguishes scoped
validation from unrelated missing-artifact failures in repository checks.
No formal verification, independent replay, external review or novelty
priority is claimed.
