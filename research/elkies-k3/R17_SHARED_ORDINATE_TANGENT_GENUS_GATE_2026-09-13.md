# The fixed generic-basis tangent construction has large genus

One tangent to a shared-ordinate cubic gives a candidate pair of sections over
one quadratic extension, starting from two generic sections. For every one of
the **136 unordered pairs** in the published R17 basis, the resulting covering
base has genus at least **86**. Thus this fixed construction cannot meet the
[infinite-rational-base goal](CORRELATED_QUADRATIC_GAINS_2026-09-12.md).

This does not exclude other generic words or iterated cubic constructions.
No independence assertion is needed or made for these candidate pairs.

## The elementary tangent identity

Write `f(x)=x^3+A*x+B`, and take generic sections
`R=(r,y_r)`, `S=(s,y_s)`, with nonzero ordinates. Put

```
a=f(r)=y_r^2,     b=f(s)=y_s^2,
F=3*r^2+A,        G=3*s^2+A,
N=a^2*G^3-b^2*F^3,
M=r*a*G^2-s*b*F^2,
U=r*N-3*M*a*G,
V=s*N-3*M*b*F.
```

The plane cubic `b*f(x1)=a*f(x2)` contains `(r,s)`. Where `a`, `G` and
`N` are nonzero, its tangent has slope `lambda=b*F/(a*G)` in the
`(x1,x2)` coordinates. Substituting `(x1,x2)=(r+h,s+lambda*h)` cancels
the constant and linear terms. The third intersection has

```
h=-3*(r-(a/b)*s*lambda^2)/(1-(a/b)*lambda^3),
x1=U/N,          x2=V/N.
```

Clearing denominators extends the polynomial identity through the other
charts. Define

```
H=U^3+A*U*N^2+B*N^3,
D=H*N.
```

Then `b*H=a*(V^3+A*V*N^2+B*N^3)`, and on `w^2=D` there are the exact
section formulas

```
P=(U/N, w/N^2),
Q=(V/N, (y_s/y_r)*w/N^2).
```

Both are anti-invariant. This supplies candidate sections, without proving
that they are independent or new. Reversing either inherited ordinate sign
does not change `a,b,N,M,D`, so all sign choices give the same branch class.
Swapping the two inherited sections also preserves the branch squareclass.

## A branch lower bound before point construction

The frozen packet consists only of the published equation and its 17 generic
polynomial sections. Take each unordered pair once and perform exactly one
tangent operation. The prime panel `1009,1013,1019` and a 40 CPU-second,
4 GiB cap were fixed before the calculation. No exceptional point or
specialization parameter is read.

For a nonzero polynomial `D` integral at an odd prime `p`, with nonzero
reduction, let `o_p` be the monic product of the irreducible factors of
`D mod p` having odd multiplicity. The reduced quadratic cover has

```
b_p=degree(o_p)+(degree(o_p) mod 2)
```

geometric branch points, including infinity. This is a lower bound on the
characteristic-zero branch degree. To see this directly, factor
`D=c*g*h^2` with `g` squarefree over Q. By Gauss's lemma, choose `g,h`
primitive and integral at `p`; nonzero reduction of `D` makes `c` a unit.
After reduction, the odd part comes from `g mod p`, hence its degree is
at most `degree(g)`. Rounding each degree up to an even integer proves
`b_p<=degree(g)+(degree(g) mod 2)`. Collisions and degree loss may weaken
the bound; they cannot create a false exclusion. Constant nonsquares alone
are not geometric branch points.

All 136 pairs have nonzero `N` and `D` reductions at the first prime, 1009.
Repeated polynomial gcds compute the odd parts without irreducible
factorization. Their degrees are below the characteristic, so no
inseparable remainder is omitted.

| Genus lower bound | Basis pairs |
|---:|---:|
| 86 | 1 |
| 87 | 135 |

These are lower bounds on the characteristic-zero normalized covers,
not claims of exact genus. In particular every candidate fails the genus
at-most-one condition, and no independence or rational-point search follows.

## Retained evidence and replay

The [frozen input](../artifacts/generated-results/elkies-k3-r17-shared-ordinate-tangents-v1/input.json),
[complete finite witnesses](../artifacts/generated-results/elkies-k3-r17-shared-ordinate-tangents-v1/result.json)
and [independent replay](../artifacts/generated-results/elkies-k3-r17-shared-ordinate-tangents-v1/independent-replay.json)
are retained. The checker uses only fractions and integer polynomial
arithmetic. It verifies all 17 inherited section equations over Q,
reconstructs every finite tangent identity and odd part, and requires exact
coverage of the 136 pairs. Four regressions reject missing pairs, altered
branch witnesses and nonintegral inputs, and retain the loss of branch
points under a colliding reduction.

```
python3 research/elkies-k3/scripts/verify_r17_shared_ordinate_tangents.py
python3 -m unittest discover -s research/tests -p 'test_r17_shared_ordinate_tangents.py' -q
```

The [execution receipt](../artifacts/generated-results/elkies-k3-r17-shared-ordinate-tangents-v1/execution.json)
records the completed bounded run. This excludes one tangent from each pair
of this fixed basis only. Other sections, further tangent iterations and
other correlated-cover identities retain their own unresolved gates.
