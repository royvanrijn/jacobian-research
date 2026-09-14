# The two common-quartic control surfaces cannot supply an MW17 parent

Neither retained control K3 admits a rational Jacobian fibration of
arithmetic rank16 or17. The genus-one control has rational Picard rank
at most12, and the genus-zero control at most10. Therefore **every**
Q-defined Jacobian fibration on those surfaces has arithmetic MW rank at
most10 and8, respectively. Changing the fibration cannot close the transfer
gate on either surface.

The [two-gain identities and infinite rational bases](COMMON_QUARTIC_SINGULARITY_LOCUS_2026-09-14.md#3-two-explicit-controls-with-all-base-and-independence-gates-closed)
remain valid. Exact inherited ranks and geometric Picard ranks remain
unknown. The [requested MW17 construction](CORRELATED_QUADRATIC_GAINS_2026-09-12.md)
remains open on other parent surfaces.

## 1. Two complete small-prime counts

Use exactly the previous controls, with `s=t^4+1`,
`A=D*(2*s+1)-1`, and `B=D*s^2`.

| Control base | D | Good prime p | Counts over t=0,...,p-1,infinity | #X(F_p) |
|---|---|---:|---|---:|
| Genus1 | `t^4+t+1` |5|7, 3, 8, 8, 8, 7|41|
| Genus0 | `t^2*(t^2+1)` |7|8, 5, 7, 6, 6, 7, 5, 5|49|

For each prime, the certificate verifies that `4*A^3+27*B^2` retains
degree24 and is squarefree. In characteristics5 and7 this gives a smooth
proper K3 model with24I1 fibres; the leading discriminant coefficient also
certifies smooth infinity. At a nodal fibre, its discriminant has a simple
zero, so the total Weierstrass surface is smooth there. Count the nodal cubic
itself, including its unique point at infinity, as a fibre of that surface.
Replacing it by its normalization would change the surface count.

The K3 infinity chart is `v=1/t`, `X=v^4*x`, `Y=v^6*y`. Its fibre is
`Y^2=X^3+A_8*X+B_12`, which explains the final entry in each count row.

## 2. An unconditional rational Picard bound

Let `rho_Q` be the rank of the Galois-invariant Neron-Severi group of a
K3 surface over Q. At good reduction, its divisor classes inject into
the Frobenius-p eigenspace of H2. Smooth proper base change identifies
the22-dimensional cohomology, and the trace formula gives

```
T_p = #X(F_p)-1-p^2.
```

Every remaining eigenvalue has complex absolute value p by
[Deligne's Weil theorem](https://numdam.org/item/PMIHES_1974__43__273_0/).
Consequently

```
T_p >= rho_Q*p-(22-rho_Q)*p = (2*rho_Q-22)*p,
rho_Q <= floor((22*p+T_p)/(2*p)).                      (1)
```

Only the injection of divisor classes is needed. Equality with a
Frobenius eigenspace, the Tate conjecture and a full zeta polynomial are
not used. The specialization and elliptic-surface interpretation are
standard; see [Schutt--Shioda, sections6 and14.6](https://arxiv.org/pdf/0907.0298).

The two traces are15 and-1. Substitution in (1) gives

```
genus1 control: rho_Q <= floor((110+15)/10) = 12,
genus0 control: rho_Q <= floor((154-1)/14)  = 10.
```

For any Q-defined Jacobian fibration, the fibre, zero and independent
rational section classes give at least `2+rank MW(Q(t))` independent
rational divisor classes. Reducible fibres can only add to this
requirement. Thus the rank bounds10 and8 apply to every such fibration,
without enumerating neighbours or constructing another equation.
These are arithmetic bounds. Divisor classes defined only over an
extension need not contribute the eigenvalue p, so (1) is not a geometric
Picard bound.

## 3. The same obstruction holds in whole coefficient residue classes

Let another degree-(8,12) integral short K3 model reduce coefficientwise
to the first control modulo5, or to the second modulo7. Its discriminant
has the same good squarefree reduction and its finite surface has the same
point count. Equation (1) gives the same rational Picard bound.

In particular, within the reverse construction

```
A=D*(2*s+1)-1,   B=D*s^2,
```

every p-integral coefficient choice with `D,s` congruent to the corresponding
control modulo p is excluded as an MW16/MW17 parent, provided the displayed
degree bounds hold. This is an entire arithmetic residue class, without a
coefficient-height bound. Rational coefficients with denominators prime to
p are included. Other residue classes and other parent constructions remain
untested. A further finite cover changes the surface and is outside the
same-surface conclusion.

## 4. Evidence and scope of the computation

The [input](../artifacts/generated-results/elkies-k3-common-quartic-control-rank-gate-v1/input.json)
binds the two original control definitions and both implementations. The
[producer](scripts/certify_common_quartic_control_rank_gate.sage) uses Sage
elliptic-curve cardinalities on smooth fibres and square characters on nodal
fibres. The [independent checker](scripts/verify_common_quartic_control_rank_gate.py)
uses only integer arithmetic: it reconstructs A,B and the discriminants,
checks polynomial gcds, and enumerates every affine `(x,y)` pair at every
base value, including infinity. Both control checkpoints are retained.

```
python3 research/elkies-k3/scripts/verify_common_quartic_control_rank_gate.py
python3 -m unittest discover -s research/tests -p 'test_common_quartic_control_rank_gate.py' -q
```

The [independent receipt](../artifacts/generated-results/elkies-k3-common-quartic-control-rank-gate-v1/independent-replay.json)
records successful replay. The frozen producer cap was10 CPU seconds and
4GiB. The [preview](../artifacts/generated-results/elkies-k3-common-quartic-control-rank-gate-v1/preview.json)
records selection of the first good prime from the fixed panel5 through43;
this was discovery of a proof on two existing models, not a prospective
high-rank search. Its first invocation lacked SymPy in the system Python;
the installed Sage runtime ran the same script successfully. No failed run
was used as evidence. The written cohomological argument is not formally
verified.

The unresolved construction still needs an actual arithmetic MW17 parent,
one quadratic cover with infinitely many rational points and two certified
new independent directions. The controls establish the low-genus and
independence mechanisms; the present result proves that their parent
surfaces cannot supply the required inherited group.
