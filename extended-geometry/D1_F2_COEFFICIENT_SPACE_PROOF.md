# D1/F2 from coefficient spaces

This note isolates the reconstruction theorem from the weighted-map notation.
It proves a slightly stronger statement than the generic formulation of D1/F2:
over characteristic zero, the affine-marked Hessian map is an isomorphism on
an explicit open, and forgetting the mark is the universal simple-root cover
of degree `N-2`.

The only repository-specific input still needed for applications to stable
left-right equivalence is the assertion that the Hessian divisor and the
distinguished affine sheet are intrinsic and are preserved.  Once those two
objects have been identified, the reconstruction below is purely a statement
about one-variable polynomials.

## 1. The source and the twice-primitive operator

Fix `N>=4`, put `n=N-2`, and work first over a field `k` of characteristic
zero.  Every normalized seed has a unique expression

```text
H(W) = W^2(1-W)Q(W),   deg(Q)<=N-3,   Q(1)=1.          (1.1)
```

Indeed, the first two factors encode `H(0)=H'(0)=H(1)=0`, while
`H'(1)=-Q(1)`.  Thus the normalized coefficient space is the affine
hyperplane

```text
A_N = {Q in k[W]_{<=N-3} : Q(1)=1} ~= A^(N-3).         (1.2)
```

For a degree-`n` polynomial `K(W)=sum(k_i W^i)`, define its twice primitive
with zero initial conditions by

```text
J_K(W) = sum_{i=0}^n k_i W^(i+2)/((i+1)(i+2)).         (1.3)
```

Then `J_K''=K` and `J_K(0)=J_K'(0)=0`.  Write

```text
P_K(W)=J_K(W)/W^2.                                    (1.4)
```

Let `E_n` be the open in the projective coefficient space `P(k[W]_{<=n})`
on which

```text
k_0 k_n Disc(K) Disc(P_K) != 0.                        (1.5)
```

These four factors respectively impose an exact double zero at the origin,
exact degree `N`, a squarefree Hessian divisor away from `0,infinity`, and
`n` distinct nonzero primitive roots.

The coordinate torus acts by

```text
a.[K(W)] = [K(aW)].                                   (1.6)
```

The unmarked target is the quotient stack

```text
X_n = [E_n/G_m].                                      (1.7)
```

Now form the universal nonzero-root incidence

```text
Z_n = {([K],rho) in E_n x G_m : P_K(rho)=0}.           (1.8)
```

Thus `rho` marks a root sheet of the twice primitive `J_K`, not a point of
the Hessian divisor itself.  This is the coefficient-space meaning of the
affine root sheet.

The action is

```text
a.([K(W)],rho)=([K(aW)],rho/a),                        (1.9)
```

and the affine-marked target is

```text
Y_n=[Z_n/G_m].                                        (1.10)
```

Because `rho` is nonzero, the action on `Z_n` is free.  Every orbit has a
unique representative with `rho=1`, so `Y_n` is a scheme, not merely a
stack.

## 2. Exact reconstruction

Let `U_N` be the open in `A_N` on which `H` has exact degree `N`, `H/W^2`
is squarefree with nonzero constant term, and `H''` is squarefree of degree
`n` with nonzero constant term.

### Theorem

The map

```text
Phi: U_N -> Y_n,       H |-> ([H''],1)                 (2.1)
```

is an isomorphism of schemes.  Under this identification, the unmarked
Hessian-divisor map is the forgetful morphism

```text
pi: Y_n -> X_n.                                       (2.2)
```

The morphism `pi` is representable, finite etale, and of constant rank
`n=N-2`.

### Proof

Given `([K],rho)` in `Z_n`, set

```text
H_(K,rho)(W) = -J_K(rho W)/(rho J_K'(rho)).            (2.3)
```

The denominator is nonzero because `rho` is a simple root of `P_K`.  Formula
(2.3) is unchanged if `K` is multiplied by a scalar.  It is also unchanged
under (1.9), since

```text
J_(K(a .))(W)=a^(-2)J_K(aW).                           (2.4)
```

Thus (2.3) descends to `Y_n`.  Directly,

```text
H(0)=H'(0)=H(1)=0,   H'(1)=-1,                         (2.5)
```

and its Hessian divisor is the divisor represented by `K`.  Conversely,
for a normalized seed, the zero initial conditions give

```text
J_(H'')=H.                                             (2.6)
```

Equations (2.3) and (2.6) show that the two maps are inverse as morphisms,
not only on geometric points.

Before quotienting, `Z_n->E_n` is the zero scheme of the squarefree
degree-`n` polynomial `P_K`.  It is therefore finite etale of rank `n`.
Finite etale morphisms descend through the torus quotient, proving the final
claim.

This simultaneously proves:

* `dim(X_n)=n-1=N-3`;
* exact generic rerooting degree `n=N-2`;
* scheme-theoretic, rather than merely set-theoretic, recovery after marking;
* generic unramifiedness without a separate tangent calculation.

The `n` points over a Hessian divisor are exactly the choices of a nonzero
root `rho` of `J_K`.  Formula (2.3) is precisely normalized rerooting at that
root.

## 3. Nonemptiness and the excluded loci

The clean open is nonempty in every degree.  One simultaneous witness for
the two squarefreeness conditions is

```text
H_0(W)=(W^2-W^N)/(N-2).                               (3.1)
```

Here

```text
H_0/W^2=(1-W^n)/n,
H_0''=(2-N(N-1)W^n)/n,                                (3.2)
```

and both polynomials are squarefree with nonzero constant and leading
coefficients in characteristic zero.  The source is irreducible, so this
also shows that all further nonempty open conditions may be imposed
simultaneously.

The precise coefficient-space exclusions are:

| locus | equation or condition | effect if retained |
|---|---|---|
| degree drop | `lc(H)=0` | fewer than `n` nonzero roots and lower Hessian degree |
| enlarged zero cluster | `H''(0)=0` | the root at zero has multiplicity at least three |
| primitive-root discriminant | `Disc(H/W^2)=0` | the root cover ramifies or becomes nonreduced |
| Hessian discriminant | `Disc(H'')=0` | leaves the squarefree-divisor target |
| scaling symmetry | nontrivial `Stab(div(H''))` | stack degree remains `n`, but coarse rerootings can coincide |

The weighted exceptional and reconstruction-boundary equations used by the
larger construction are additional open restrictions.  They are not needed
for the coefficient-space theorem; one simply intersects `U_N` with them
when transferring the result back to that construction.

## 4. Stabilizers and symmetry fibers

Write `K=sum k_iW^i` with `k_0!=0`.  Over an algebraically closed field, its
scaling stabilizer is

```text
Stab([K]) = mu_g,
g = gcd{i>=1 : k_i!=0}.                               (4.1)
```

The generic stabilizer is trivial.  The symmetry locus is the finite union,
over divisors `m>1` of `n`, of the linear loci

```text
K(W) in k[W^m].                                       (4.2)
```

On the stack, `Y_n->X_n` remains finite etale of rank `n` along these loci.
On the coarse space, a divisor with stabilizer `mu_g` has `n/g` distinct
normalized rerootings: the `g` selected roots in one stabilizer orbit give
the same normalized seed.  Thus symmetry produces collapsed coarse fibers,
not additional generic fiber points.  Removing (4.2) is exactly what is
needed for the literal statement that a geometric coarse fiber consists of
`n` distinct seeds.

The mark kills every stabilizer: an element fixing a nonzero selected root
`rho` must be the identity.  This is the geometric reason that the marked
target is a scheme and why (2.1) stays injective even above the unmarked
symmetry locus.

## 5. Characteristic

The proof works verbatim over `Z[1/N!]`, hence over every field of
characteristic zero or characteristic `p>N`.  This is a clean uniform range,
not an artifact of the proof: reductions in smaller characteristic can
degenerate severely.

* In characteristic `2`, every second derivative is zero.
* In characteristic `3`, `(H'')'=H'''=0`; a positive-degree Hessian is
  inseparable, so the squarefree Hessian-clean open is empty.
* If `p` divides `N(N-1)`, the Hessian of every exact degree-`N` seed has
  degree less than `N-2`.
* Frobenius terms can destroy marked faithfulness even when the expected
  Hessian degree and squarefreeness survive.  In characteristic `p`, put

```text
L_p(W)=W^p-W^(p+1),
H_c=cG+(1-c)L_p.                                      (5.1)
```

  For every normalized `G` and every nonzero `c`, one has

```text
H_c(0)=H_c'(0)=H_c(1)=0,   H_c'(1)=-1,
H_c''=cG''.                                           (5.2)
```

  Hence all `H_c` have the same projective Hessian and the same selected
  root `1`.  Whenever two members remain in the clean open, F2 fails even
  set-theoretically; over a function field the fiber is positive
  dimensional.  Degree `8` in characteristic `5` supplies explicit clean
  examples below.

Thus reductions at primes dividing `N!` are useful collision detectors but
are not models of the characteristic-zero cover.  The uniform theorem
should be stated in characteristic zero (or over `Z[1/N!]`), not over an
unspecified base.

The coefficient-recovery failure itself has a sharp positive-characteristic
repair: the
[`p`-typical Hasse reconstruction theorem](HASSE_TYPICAL_SEED_RECOVERY.md)
recovers every degree-at-most-`N` polynomial modulo constants from the
`floor(log_p(N))+1` channels
`D^[1],D^[p],...,D^[p^floor(log_p(N))]`, and proves that every one of those
orders is necessary among families of raw Hasse derivatives.  This does not
repair the separate map-intrinsic problem.  In fact, the same note proves a
clean degree-twelve no-go theorem: distinct normalized `F_5` seeds can have
the same ordinary derivative and therefore define the identical weighted
polynomial map, while their fifth Hasse channels differ.  It then constructs
a dimension-preserving Frobenius correction of that map and proves exact
stable seed recovery on the odd-characteristic tame clean locus using the
full second-boundary factor data.

## 6. Exhaustive low-degree audit

[`verify_d1_f2_coefficient_reconstruction.py`](../scripts/verify_d1_f2_coefficient_reconstruction.py)
checks (2.3)--(2.6) over the symbolic function fields in degrees `4,5,6,8`,
including invariance under both coefficient scaling and coordinate scaling.

[`search_d1_f2_coefficient_collisions.py`](../scripts/search_d1_f2_coefficient_collisions.py)
exhausts the rational normalized seeds for the requested degrees.  It tests
both primitive and Hessian discriminants, groups by the marked projective
Hessian, and separately groups modulo coordinate scaling.

The main results are:

| `N` | field | clean seeds | marked collisions | largest unmarked rational fiber |
|---:|---:|---:|---:|---:|
| 4 | `F_3` | 0 | — | — |
| 4 | `F_5` | 3 | 0 | 2 = `N-2` |
| 5 | `F_3` | 0 | — | — |
| 5 | `F_7` | 27 | 0 | 3 = `N-2` |
| 6 | `F_3` | 0 | — | — |
| 6 | `F_5` | 0 | — | — |
| 6 | `F_7` | 201 | 0 | 4 = `N-2` |
| 8 | `F_3` | 0 | — | — |
| 8 | `F_5` | 1460 | 418 marked collision classes | 16 (failure) |
| 8 | `F_7` | 0 | — | — |
| 8 | `F_11` | 112760 | 0 | 6 = `N-2` |

The empty rows agree with the characteristic obstructions above.  Over
`F_5`, the degree-eight scan finds the clean marked collision

```text
H_1 = W^2+4W^3+W^6+3W^7+W^8,
H_2 = 2W^2+3W^3+4W^5+3W^6+W^7+2W^8,                 (6.1)
```

with

```text
[H_1'']=[H_2'']=[1,2,0,0,0,3,3].                     (6.2)
```

Both seeds pass the exact-degree, primitive-squarefree, and
Hessian-squarefree tests.  In the safe primes `p>N`, no marked collision was
found and the largest rational unmarked fiber has exactly the predicted
size.  Rational fiber counts are only diagnostics—non-rational root sheets
are invisible over `F_p`—whereas the incidence proof in Section 2 proves the
geometric degree.

## 7. Verdict for D1 and F2

On the explicit characteristic-zero open `U_N`, D1 and the coefficient
content of F2 are valid with the following precise strength:

```text
unmarked: finite etale of stack degree N-2;
free coarse locus: N-2 distinct geometric rerootings;
affine-marked: scheme-theoretic isomorphism, hence exact recovery.          (7.1)
```

The birational/normalization statement in the existing formulation follows
immediately, but is weaker than (7.1).  What remains external to this proof
is not polynomial reconstruction; it is the functorial identification of
the selected affine sheet inside the weighted stable-equivalence problem.
