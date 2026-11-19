# Two-sided protected rank-one defects

A protected copy of the four-strand braid `U` forces a dependence among
its four input roots and, separately, a dependence among its four input
covectors. Thus `k` disjoint protected blocks force `ar_C<=b-k`, with
no preserved-form hypothesis. The earlier disjunction and
`b-ceil(k/2)` bound in [Multiple protected obstructions](MULTIPLE_PROTECTED_OBSTRUCTIONS.md)
remain valid weaker consequences. The actual signed relations are essential:
the abstract opposite-side linear-algebra control there does not satisfy them.

This result supplies stronger augmentation upper bounds. A compatible
meridional lower bound, prime realization and hyperbolicity are separate
obligations; no knot of augmentation gap two is established here.

## 1. Statement and crossing trace

Use the left-to-right Hurwitz convention `T(a,c)=(c,c^-1*a*c)` and set

```
U=[1,2,-3,-3,2,-3,-3,2,1,1,2,2,2,-3].
```

**Theorem 1 (two-sided local defect).** Let four invertible matrices act
on any finite-dimensional space over any field, with nonzero displacement
rank one and one common determinant `mu!=0`. Propagate them through `U`.
If the first and fourth output matrices equal their respective inputs,
then the four input root lines span a space of dimension at most three,
and their four covector lines separately span a space of dimension at
most three. The assertion includes `mu=1` and characteristic two.

Write every matrix as `M=I-u*f`; then `f(u)=eta=1-mu` and
`M^-1=I+mu^-1*u*f`. A pair for a nonzero rank-one displacement is unique
up to reciprocal nonzero scaling. Conjugation by `J` sends its pair to
`(J*u,f*J^-1)`.

Here is the full crossing trace, also used in the
[current released proof](../../../results/prime-knot-augmentation-gaps/representations.tex),
Lemmas `lem:spans` and `lem:seven`. A row `(o,a,c,epsilon)` means
`M_c=M_o^epsilon*M_a*M_o^-epsilon`. The inputs are `(1,0,8,13)`;
the outputs are `(1,6,16,13)`. Repeating the two outer labels imposes
exactly the two assumed return equalities.

| Crossing | `o` | `a` | `c` | `epsilon` |
|---:|---:|---:|---:|:---:|
|0|0|1|2|-|
|1|8|2|3|-|
|2|3|13|14|+|
|3|14|3|4|+|
|4|4|8|9|-|
|5|9|14|15|+|
|6|15|9|10|+|
|7|10|4|5|-|
|8|10|0|1|-|
|9|1|10|11|-|
|10|5|11|12|-|
|11|12|5|6|-|
|12|6|12|13|-|
|13|13|15|16|+|

We first prove that four independent input covectors are impossible;
contragredient duality will then give the root assertion.

## 2. The two folds require only covector independence

Suppose the four input covectors are independent. Put
`A=M_2, B=M_8, C=M_13, D=M_0`. Replacing `M_1` by
`A=D^-1*M_1*D` adds a multiple of the covector of `D` to its
covector, preserving independence. Let `F` be the three-dimensional
covector span of `A,B,C`, and write `(v,nu)` for the pair of `D`.
Then `nu` is outside `F`. The internal crossings propagate covectors
inside `F`, regardless of possible dependencies among their roots.

Rows 0 and 8 give `D*A*D^-1=X^-1*D*X`, where `X=M_10`.
Write `(a,f),(x,xi)` for pairs of `A,X`; both covectors lie in `F`.
For some `lambda!=0`, equality of the rank-one displacements gives

```
a-nu(a)v = lambda*(v+mu^-1*xi(v)x),
f+mu^-1*f(v)nu = lambda^-1*(nu-nu(x)xi).
```

Comparing the coefficient of `nu` in the second equality gives
`f(v)=mu/lambda`; the remaining part gives `f=c*xi`, where
`c=-nu(x)/lambda!=0`. Apply `f` to the first equality. Its right
side evaluates to `mu+eta=1`, so `nu(a)=-lambda`.
The first equality therefore reduces to
`a=(lambda/mu)*xi(v)*x`. Since `c*xi(v)=mu/lambda`, this is
`a=x/c`. Thus `a*f=x*xi`, proving `X=A`.
The relation becomes `D*A*D=A*D*A`; row 9 then gives `M_11=D`.

Rows 10--12 take `(D,P)` to `(P,N)`, then `(N,Q)`, then `(Q,C)`,
where `P=M_5, N=M_12, Q=M_6`. For a pair `(p,pi)` of `P`, set
`x_1=pi(v), y_1=nu(p), h_1=1+x_1*y_1/mu`. Direct conjugation gives

| Matrix | Root | Covector |
|---|---|---|
|`N`|`v+(x_1/mu)p`|`nu-y_1*pi`|
|`Q`|`y_1*v+h_1*p`|`-(x_1/mu)nu+h_1*pi`|
|`C`|`h_1*v+(x_1/mu)(1+h_1)p`|`h_1*nu-y_1(1+h_1)pi`|

The last covector is in `F`, so `h_1=0` and `x_1*y_1=-mu`.
The remaining root and covector coefficients multiply to one;
therefore `P=C`. No root-independence hypothesis has been used.

It remains to exclude three independent covectors in the internal
seven-relation core. This stronger core lemma is the new step.

## 3. The seven-relation core has both spans at most two

**Lemma 2 (covector-only core obstruction).** Suppose the signed
relations in rows 1--7 hold, with `M_10=M_2` and `M_5=M_13`, and
all matrices have nonzero rank-one displacement and common nonzero
determinant. Their covector span has dimension at most two.

**Proof.** Reorder the relations into five propagations followed by
two remaining equalities:

```
M_15=M_9*M_14*M_9^-1,
M_2=M_15*M_9*M_15^-1,
M_3=M_8^-1*M_2*M_8,
M_4=M_14*M_3*M_14^-1,
M_13=M_3^-1*M_14*M_3;
M_4*M_9*M_4^-1=M_8,
M_2*M_13*M_2^-1=M_4.                                  (1)
```

All covectors are propagated from those of `M_8,M_14,M_9`.
If the span had dimension three, choose their pairs as
`(u,f),(v,g),(w,j)`, with `f,g,j` independent. The roots may
be dependent. Define

```
x=f(v), X=g(u), y=f(w), Y=j(u), z=g(w), Z=j(v),
eta=1-mu, h=1+zZ/mu, t=-zx/mu+h*y, ell=ZX+hY,
kappa=tX/mu+z*(h-eta/mu), Lambda=-ell*x+(eta+h)Z.         (2)
```

Using formal coefficients on `(u,v,w)` and on `(f,g,j)`, the first
five equalities in (1) give the following actual pairs. The root
coefficients need not be unique, and the calculation never compares
them using independence.

| Matrix | Root coefficients | Covector coefficients |
|---|---|---|
|`M_2`|`(0,-z/mu,h)`|`(0,Z,h)`|
|`M_3`|`(t/mu,-z/mu,h)`|`(-ell,Z,h)`|
|`M_4`|`(t/mu,-z/mu-kappa,h)`|`(-ell,Z+Lambda/mu,h)`|
|`M_13`|`(Lambda*t/mu^2,1-Lambda*z/mu^2,Lambda*h/mu)`|`(kappa*ell,1-kappa*Z,-kappa*h)`|

Let `P=f_4(w), Q=j(u_4)`. The sixth equality in (1) says, for a
nonzero `gamma`,

```
w-P*u_4=gamma*u,       j+(Q/mu)f_4=gamma^-1*f.
```

Only the second equality uses independent coordinates. Its `j`, `g`
and `f` coefficients give

```
Qh=-mu,     Lambda=-mu*Z,     gamma=h/ell.               (3)
```

Thus `Q,h,ell` are nonzero. For the last equality in (1), put
`R=f_2(u_13), S=f_13(u_2)` and take `beta!=0` with

```
u_13-R*u_2=beta*u_4,    f_13+(S/mu)f_2=beta^-1*f_4.
```

Its covector coefficients successively give

```
kappa=-beta^-1,     S=0,     kappa*Lambda=-mu.
```

Together with (3), these imply `kappa*Z=1`. Evaluating the displayed
pairs, with this identity and `zZ=mu(h-1)`, gives

```
S=kappa*(ell*t-h+mu*h^2),
Q=tY/mu-mu*h.
```

Therefore

```
ell*t=h-mu*h^2,                       tYh=mu^2(h^2-1).   (4)
```

Multiplying the definition of `kappa` by `mu*Z` and simplifying gives

```
tXZ=mu^2-mu^2*h^2+mu*h.                                (5)
```

Add (5) to the second equality of (4), and use `ell=ZX+hY`.
It follows that `t*ell=mu*h`. Compare with the first equality of
(4). Since `h!=0`,

```
mu*h=eta,        t*ell=mu*h!=0.                         (6)
```

The definition of `kappa` now simplifies to `kappa=tX/mu`;
hence `tXZ=mu` and `XZ/ell=1/h`.

Evaluate the root equality `w-P*u_4=(h/ell)u` using `g`.
The definition of `kappa` gives `g(u_4)=mu*kappa`.
Multiplying the resulting equality by `Z`, and using
`zZ=mu(h-1), kappa*Z=1, XZ/ell=1/h`, yields

```
mu*(h-1-P)=1.
```

Together with (6) this forces `P=-2`.

But (3) and `kappa*Z=1` also give `f_13=-kappa*f_4`.
Since `S=0`, we have `f_4(u_2)=0`. Directly
`f_4(v)=Lambda/mu=-Z`, so

```
0=f_4(u_2)=zZ/mu+hP=h-1+hP.
```

Consequently `h(P+2)=1+h=1/mu`, by (6). Its left side is zero
because `P=-2`, a contradiction. All divisions were by quantities
already shown nonzero. There is no division by two; if `mu=1`,
(6) already contradicts `h!=0`. This proves the lemma. QED.

The contragredient operation `M -> M^-T` preserves every signed
conjugation relation and has determinant `mu^-1`. Indeed its
nonzero displacement is `-mu^-1*f^T*u^T`, interchanging the root
and covector lines. Applying Lemma 2 to the contragredient proves
that the core root span also has dimension at most two.

To finish Theorem 1, the folds give the core relations. Its covector
span includes those of `M_2,M_8,M_13`, so under the assumed four
input covectors it has dimension three. This contradicts Lemma 2.
Thus the four input covectors are dependent. Apply the same whole
argument to the contragredient tuple to obtain dependence of the
four input roots as well. QED.

## 4. Disjoint blocks and the revised coloring target

**Theorem 3 (independent defects on both sides).** Under the braid,
tail and knot-closure hypotheses of Theorem 1 in
[Multiple protected obstructions](MULTIPLE_PROTECTED_OBSTRUCTIONS.md),
with `k` disjoint copies on `b>=4k` strands, every rank-one meridional
representation has root and covector dimensions

```
r<=b-k,       q<=b-k.                                  (7)
```

In particular an irreducible representation has dimension at most
`b-k`, and in the conventional complex `U=1, mu!=1` setting,

```
ar_C(K)<=b-k.                                         (8)
```

**Proof.** The tail fixes the two outer free meridians of each block.
Closure therefore supplies the two local return equalities.
Theorem 1 gives a nonzero relation among input roots supported on
each block, and another among input covectors supported on that
same block. The `k` relations on either side are independent because
their nonempty supports lie in disjoint strand sets. The maps from
`F^b` to the two spans consequently have kernels of dimension at
least `k`, proving (7). The input spans equal the global meridian
spans, since input meridians generate the group and conjugation
preserves both spans. In an irreducible representation the root
span is all the space and the covector common kernel is zero.
Cornwell's correspondence then gives (8). QED.

The explicit knot closures `K_(k,b,l)` of the linked note satisfy
these hypotheses. Two blocks now give an augmentation bound two
below their displayed strand count. A meridional lower bound is
still needed to turn that into a gap.

The [complete double-transposition coloring companion](PROTECTED_U_DOUBLE_TRANSPOSITION_CLASSIFICATION.md)
now proves the stronger support bound v<=2b-2k+2 for every exact-word
protected tail with a knot closure. Its support lower bound
mr>=ceil((v-1)/2), combined with (8), can certify a gap of at most one.
Thus changing the number of blocks or their double-transposition colors
cannot complete this particular lower-bound route. The augmentation
bound (8) remains valid; other lower invariants, finite meridian classes,
overlapping obstructions and local braids remain separate targets.

## 5. Primary input, exact checks and assurance

The crossing trace and preliminary fold formulas are taken from the
linked current proof, and the covector-only argument above establishes
the stronger local result directly. The augmentation correspondence is
C. R. Cornwell, *KCH representations, augmentations, and A-polynomials*,
J. Symplectic Geom. 15 (2017), 983--1017, Theorem 1.2 and Definition 1.8;
its [primary preprint](https://arxiv.org/pdf/1310.7526) explicitly
restricts augmentation rank to `epsilon(mu)!=1` and identifies it with
the corresponding irreducible KCH dimension. The local algebra here
has the broader stated scope. No historical priority is asserted.

The [exact symbolic checker](verify_two_sided_protected_obstructions.py)
reconstructs the five conjugations using a generic pairing matrix and
checks every coefficient of the propagation table. It verifies the
scalar reductions, rejects reversed-crossing and sign mutations,
and checks an actual rank-three reflection representation of `U`
with both outer returns. The latter attains both local bounds and
has determinant `-1`. It also checks the contragredient identity.
The written argument supplies the unrestricted theorem; its proof
uses integer identities and invertible scalars, so the rational
symbolic replay does not itself supply the characteristic-two scope.
That scope follows from the displayed proof without division by two.

```
.venv/bin/python research/meridional-rank/representations/verify_two_sided_protected_obstructions.py \
  --write-receipt research/meridional-rank/representations/two_sided_protected_obstructions_verification.json
```

The [receipt](two_sided_protected_obstructions_verification.json) records
source pins and exact replay scope. This is internal written and
computational evidence. No external review, formal verification,
prime realization, hyperbolicity or augmentation-gap-two example is
inferred from the checker.
