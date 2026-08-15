# R20 Brumer--Kramer audit

Status: exact arithmetic theorem and an unfinished class-group upper-bound
computation. No new rational point is claimed.

The canonical Fermigier--Mestre specialization
`fermigier-mestre-v1:u=28917/20` has an exact finite-reduction certificate for
twenty independent rational points and global minimal model

```text
y^2 + x*y + y = x^3 + x^2
  - 4437412060110743641525245114305*x
  + 3586842216822165612930264910099076801587288127.
```

The purpose of this audit is to move from undirected point search to the cubic
class-group obstruction controlling its 2-Selmer group.

## The cubic field

For a generalized Weierstrass model `[a1,a2,a3,a4,a6]`, nonzero 2-torsion
abscissas satisfy

\[
4x^3+b_2x^2+2b_4x+b_6=0.
\]

For R20 this primitive polynomial has ascending coefficients

```text
[14347368867288662451721059640396307206349152509,
 -17749648240442974566100980457218,
 5,
 4].
```

Putting `z=4*x` gives the monic defining polynomial

\[
f(z)=z^3+5z^2
-70998592961771898264403921828872z
+229557901876618599227536954246340915301586440144.
\]

Modulo 37 it is `z^3+5z^2+7z+33`, with no root in `F_37`; hence it is
irreducible. In particular `E(Q)[2]=0` and this defines the cubic subfield
`K` of the 2-division field.

Exact maximal-order computation gives

```text
Disc(K) =
17207612547621358265560224336784329653572551167050221201938192360

[O_K : Z[z]] = 712863540480000.
```

The field is totally real. The large index is important: using only the
polynomial discriminant would give incorrect ramification and class-group
input.

## Brumer--Kramer local term

Klagsbrun--Sherman--Weigandt restate the Brumer--Kramer bound as

\[
\dim_{\mathbf F_2}\operatorname{Sel}_2(E/\mathbf Q)
\leq
\dim_{\mathbf F_2}\operatorname{Cl}(K)[2]+u(E)+n(E),
\]

where `u(E)=2` for positive minimal discriminant,

\[
n(E)=\#\Phi_m+\sum_{p\in\Phi_a}(n_p-1),
\]

`Phi_m` contains the multiplicative primes with even minimal-discriminant
valuation, and `n_p` is the number of primes of `K` over an additive prime.
The primary source is
[*The Elkies Curve has Rank 28 Subject only to GRH*](https://arxiv.org/abs/1606.07178),
Proposition 3.1.

For R20, exact factorization of the conductor and minimal discriminant gives

```text
Phi_m = {3, 7, 13, 31, 79}.
```

The only additive prime is 17. Exact maximal-order prime decomposition gives
one prime above 17, with ramification index 3 and residue degree 1. Thus its
contribution is `1-1=0`, and

```text
u(E) = 2
n(E) = 5.
```

Combining this with the exact rank lower bound gives

\[
20
\leq \operatorname{rank}E(\mathbf Q)
\leq \dim \operatorname{Sel}_2(E/\mathbf Q)
\leq \dim \operatorname{Cl}(K)[2]+2+5.
\]

Therefore

\[
\boxed{\dim_{\mathbf F_2}\operatorname{Cl}(K)[2]\geq 13.}
\]

This lower bound is unconditional. It does not use GRH, BSD, numerical heights
or a class-group computation.

## The decisive upper-bound gate

The target value is now exact and unusually sharp:

- if an unconditional class-group quotient certificate proves
  `dim Cl(K)[2] <= 13`, then `rank E(Q)=20` and this fixed fibre is closed;
- a GRH-conditional upper bound 13 gives the same conclusion conditional on
  the relevant GRH assumptions;
- if the 2-rank is greater than 13, the surplus identifies the correct relative
  2-Selmer/cover lane rather than authorizing more blind point search.

PARI's ordinary `bnfinit` result is conditional on GRH. The useful rigorous
shortcut is `bnfcertify(bnf,1)`: it only needs to certify that the true class
group is a quotient of the computed group. If the computed 2-rank is 13, this
upper quotient certificate and the independent lower bound above prove exact
2-rank 13 without expanding the fundamental units.

## Replay

```bash
PYTHONPATH=elliptic-curves/cas .venv/bin/python \
  elliptic-curves/cas/audit_r20_brumer_kramer.py

PYTHONPATH=elliptic-curves/cas .venv/bin/python \
  elliptic-curves/cas/audit_r20_brumer_kramer.py --check

PYTHONPATH=elliptic-curves/cas .venv/bin/python -m unittest \
  elliptic-curves/tests/test_r20_brumer_kramer.py
```

The pinned machine record is
[`elliptic_r20_brumer_kramer.json`](../../artifacts/generated-results/elliptic_r20_brumer_kramer.json).
