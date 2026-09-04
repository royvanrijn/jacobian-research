# Determinantal closure of the planar `(72,108)` Laurent systems

## Result and claim boundary

The two Laurent systems in Guccione--Guccione--Horruitiner--Valqui (GGHV),
Proposition 4.3, are inconsistent over characteristic zero.

The repository already contained exact unit certificates for both systems. The
new contribution here is a structural reconstruction that removes the only
large opaque step in Case 1 and completes the case-specific interface audit:

1. the final Newton tail and the two polygons are derived by an exact finite
   census;
2. the first Wronskian block descends to its intrinsic `mu_7` quotient before
   elimination;
3. Case 2 is organized as a binary-octic complete intersection;
4. the Case-1 hard ideal is proved to be exactly the two-generator ideal
   `(h,N)` by a small determinantal computation.

Consequently, **assuming the published general GGV/GGHV reduction and degree
classification**, the degree pair `(72,108)` (and its swap) is excluded and the
published frontier becomes `125`. This is not a proof of the planar Jacobian
conjecture, and it is not a from-scratch reproof of the general theorem ladder
or of the earlier admissible-chain census.

The case-specific front-end audit is recorded in
[`PROP43_EXHAUSTIVENESS_AUDIT.md`](PROP43_EXHAUSTIVENESS_AUDIT.md).

## 1. Correct Laurent systems

Put

\[
 t=xy^2,\qquad z=y^{-1}.
\]

Then `[t,z]_{x,y}=-1` and `x^2=t^2z^4`. The upper bands have the form

\[
 P=A(t)z^2+B(t)z+C(t),\qquad
 Q=D(t)z^3+E(t)z^2+F(t)z+G(t),
\]

and `[P,Q]=x^2` is equivalent to

\[
2AD'-3A'D=t^2, \tag{J4}
\]

\[
2(AE'-A'E)+(BD'-3B'D)=0, \tag{J3}
\]

\[
(2AF'-A'F)+(BE'-2B'E)-3C'D=0, \tag{J2}
\]

\[
2AG'+(BF'-B'F)-2C'E=0, \tag{J1}
\]

\[
BG'-C'F=0. \tag{J0}
\]

The two Newton-polygon alternatives are

\[
\begin{aligned}
N(P)&=\operatorname{conv}\{(0,0),(1,0),(8,14),(8,16),(0,8)\},\\
N(Q)&=\operatorname{conv}\{(0,0),(2,1),(12,21),(12,24),(0,12)\},
\end{aligned}
\]

and

\[
\begin{aligned}
N(P)&=\operatorname{conv}\{(0,0),(1,0),(8,14),(8,16)\},\\
N(Q)&=\operatorname{conv}\{(0,0),(2,1),(12,21),(12,24)\}.
\end{aligned}
\]

The `(0,1)` printed in the displayed Case-1 statement of Proposition 4.3 is an
apparent typo for `(8,14)`: the proof's final monomial map gives `(8,14)`, and
`(0,1)` is interior to the edge `(0,0)--(0,8)`.

## 2. Finite Newton-tail closure

For a current corner `(A,B)` and proper proportional continuation `(a,b)`, set

\[
g=\gcd(A-a,B-b),\qquad
(p,q)=\frac1g(A-a,B-b).
\]

The auxiliary-element divisibility used in the cited GGV/GGHV framework
reduces to the scale-free condition

\[
\boxed{Bp-Aq\mid A-B}.
\]

At `(A,B)=(24,7)`, the defect is the prime `17`. The complete lattice census
leaves only

\[
(17,5),\quad(10,3),\quad(3,1),
\]

plus the forbidden diagonal point `(1,1)`. Every proper survivor lies on the
primitive direction `(7,2)`, so the successor normal is `(-2,7)`.

All candidates satisfy `2a=7b-1`. For the terminal endpoint set
`{(-k,0),(k+1,1)}`, the two parallelism determinants reduce to

\[
D_-=(k-1)(5b-1),\qquad
2D_+=-(2k-3)(5b-1).
\]

Thus `k=1` is forced, with

\[
\operatorname{en}(P)=(-1,0),\qquad
\operatorname{en}(Q)=(2,1).
\]

The three proper proportional candidates are collinear interior points of the
same `(7,2)` edge, not new Newton vertices. Hence there is no third Laurent
polygon.

## 3. Quotient-first first block

Let `u=a_7`. The first Wronskian block has the residual action

\[
P(t,z)\longmapsto \zeta^{-1}P(\zeta t,z),\qquad
Q(t,z)\longmapsto \zeta^{-2}Q(\zeta t,z),\qquad \zeta^7=1.
\]

On coefficients, `a_i -> zeta^(i-1) a_i`. Introduce the invariants

\[
q=u^7,\qquad x_i=a_i u^{i-1}\quad(2\le i\le6).
\]

After the invertible rescaling `t=u tau`, the complete `(J4)` block retains its
normalized form and all eleven triangular solves and six obstruction equations
belong to

\[
\mathbf Q[q,x_2,\ldots,x_6].
\]

The degree-35 presentation is therefore a rank-seven pullback of an intrinsic
degree-five quotient. The exact quotient graph is checked independently in
`Q[q]/(H_5(q))`; the degree-35 field is not intrinsic to the first block.

## 4. Case 2: binary-octic complete intersection

After solving `(J3)` and `(J2)`, the two cubic residuals are affine-linear in
`h`. Eliminating `h` homogenously against either quartic residual produces two
binary octics `N_7(r,s),N_9(r,s)`. The archived exact chart calculations imply
that they are coprime, hence form a complete intersection of type `(8,8)`.

Its Hilbert vector is

```text
1,2,3,4,5,6,7,8,7,6,5,4,3,2,1
```

so the quotient has length `64`, socle degree `14`, and

\[
(r,s)^{15}\subseteq(N_7,N_9).
\]

Combining this with the exact origin Bezout identity gives the unit ideal. This
repackages the existing Case-2 certificate globally and makes the separate
singular-Cramer subcase redundant.

## 5. Case 1: the hard ideal is `(h,N)`

Work over the exact quintic field

\[
L=\mathbf Q[w]/(w^5-w^4+3w^3+3w^2+26)
\]

and set `S=L[h,u_1,u_2]`. On either sign branch, the quotient-first equations
have the form

\[
G_0=b_h(hu_3-N),\qquad G_i=a_i u_3+b_i.
\]

Eliminating `u_3` gives four hard residuals

\[
F_i=h b_i+a_iN,\qquad 1\le i\le4.
\]

Let

\[
M=\begin{pmatrix}
a_1&a_2&a_3&a_4\\
b_1&b_2&b_3&b_4
\end{pmatrix},
\qquad
m_{ij}=a_i b_j-a_j b_i.
\]

The basic syzygies are

\[
a_iF_j-a_jF_i=h m_{ij}, \tag{5.1}
\]

\[
b_iF_j-b_jF_i=-N m_{ij}. \tag{5.2}
\]

### 5.1 Determinantal two-generator lemma

The reduction used here is a general ring-theoretic lemma. Let `R` be a
commutative ring, let

\[
F_i=h b_i+a_iN,\qquad I=(F_1,\ldots,F_n),\qquad
m_{ij}=a_i b_j-a_j b_i.
\]

Suppose that, for some index `k` and some set of pairs `E`,

\[
1=dF_k+\sum_{(i,j)\in E}c_{ij}m_{ij}. \tag{5.3a}
\]

Multiplication by `h`, followed by
`h m_{ij}=a_iF_j-a_jF_i`, gives

\[
h=hdF_k+\sum_{(i,j)\in E}c_{ij}(a_iF_j-a_jF_i)\in I. \tag{5.3b}
\]

If, in addition, the special-fibre coefficients generate the unit ideal,
equivalently

\[
1=sh+\sum_iq_i a_i, \tag{5.3c}
\]

then

\[
\begin{aligned}
N
 &=shN+\sum_iq_i a_iN\\
 &=\sum_iq_iF_i+h\left(sN-\sum_iq_i b_i\right)\in I.
\end{aligned} \tag{5.3d}
\]

Since every `F_i` already lies in `(h,N)`, these two hypotheses imply

\[
\boxed{I=(h,N).} \tag{5.3e}
\]

In particular, an exact unit-ideal computation for (5.3a) is enough: the
standard-basis transformation coefficients need not be tracked or printed.
This distinction is useful here because coefficient tracking is far more
expensive than deciding the characteristic-zero unit ideal.

### 5.2 Adjacent-minor unit identity

An exact characteristic-zero standard-basis calculation proves

\[
\boxed{(m_{12},m_{23},m_{34},F_1)=S.} \tag{5.3}
\]

Thus there exist `c_12,c_23,c_34,d in S` with

\[
1=c_{12}m_{12}+c_{23}m_{23}+c_{34}m_{34}+dF_1.
\]

Multiplying by `h` and using (5.1) gives a four-generator identity

\[
\begin{aligned}
h={}&(hd-c_{12}a_2)F_1
 +(c_{12}a_1-c_{23}a_3)F_2\\
 &+(c_{23}a_2-c_{34}a_4)F_3
 +(c_{34}a_3)F_4.
\end{aligned} \tag{5.4}
\]

Therefore

\[
\boxed{h\in(F_1,F_2,F_3,F_4).} \tag{5.5}
\]

The exact characteristic-zero replay decides the unit ideal (5.3). By the
lemma above this already proves (5.5); an explicit characteristic-zero dump of
the multipliers is not a proof obligation.

The six minors alone are **not** the unit ideal: at the pinned good prime `71`
their rank-drop scheme has length `72`. The compatibility equation `F_1` is
essential. This is why the naive unit-minor shortcut fails while (5.3) works.

### 5.3 Special-fibre Bezout identity

Modulo `h`, the first two coefficients are affine-linear in `u_2`:

\[
a_1\equiv\alpha_1u_2+\beta_1,\qquad
a_2\equiv\alpha_2u_2+\beta_2.
\]

Exact arithmetic in `L` gives

\[
\Delta=\alpha_1\beta_2-\alpha_2\beta_1\ne0.
\]

Hence, with

\[
c_1=-\alpha_2/\Delta,\qquad c_2=\alpha_1/\Delta,
\]

we have

\[
c_1a_1+c_2a_2\equiv1\pmod h.
\]

It follows that

\[
c_1F_1+c_2F_2\equiv N\pmod h.
\]

Together with (5.5), this proves `N` belongs to the hard ideal. The reverse
inclusion is immediate from `F_i=h b_i+a_iN`. Therefore

\[
\boxed{(F_1,F_2,F_3,F_4)=(h,N).} \tag{5.6}
\]

The 89 MB membership certificate is thus replaced by a codimension-two
complete-intersection statement and a small adjacent-minor standard-basis
calculation.

### 5.4 Exclusion of the two sign branches

The first sign branch is reconstructed directly from
`case1_branch1_after_w.pkl` and checked coefficient-by-coefficient against the
archived degree-five residuals. The second branch is transported by the exact
sign involution fixing `h` and negating `(u_1,u_2)`, with the recorded row
signs.

Once (5.6) is known, every solution has `h=N=0`. The existing exact `h=0`
certificates then show that the remaining Case-1 equations generate the unit
ideal on this special fibre. Both exhaustive sign branches are therefore
empty.

## 6. Consequence

The repository now has the following audited chain:

```text
published GGV/GGHV structural inputs
        -> exact finite Proposition-4.3 tail audit
        -> exactly two corrected Laurent systems
        -> quotient-first exact coefficient systems
        -> Case 1 and Case 2 characteristic-zero contradictions
```

Accordingly:

- there is no surviving standard pair among the two Proposition-4.3 systems;
- `(72,108)` and `(108,72)` are closed as exact computer-assisted
  consequences of the published reduction;
- conditional on the published degree classification, every planar Jacobian
  counterexample has maximum degree at least `125`.

The older `F_2` connector/braid programme remains interesting as an independent
geometric derivation, but it no longer blocks closure of this degree pair.

## 7. Replay

The top-level verifier composes the exact adjacent-minor decision, both
serialized special-fibre unit certificates, and the sign-branch transport:

```bash
python scripts/verify_jc72_108_case1_determinantal_closure.py
```

For cleanup-only verification, check both archive manifests, the three
otherwise unmanifested transport inputs, the compact four-row reconstruction,
and the special-fibre Bezout row without starting Singular or multiplying the
89 MB certificate:

```bash
python scripts/verify_jc72_108_case1_determinantal_closure.py \
  --audit-existing-only
```

The standard-library front-end audit is:

```bash
python plane-jc/cas/verify_prop43_exhaustiveness.py
```

The quotient-first checks are:

```bash
python plane-jc/cas/verify_mu7_quotient_firstblock.py
python plane-jc/cas/research_72_108_gap_attack.py
python plane-jc/cas/verify_firstblock_quotient_graph.py
```

Generate the lower-level exact Case-1 rank-drop input directly:

```bash
python scripts/research_jc72_108_case1_rankdrop.py /tmp/case1-rankdrop.sing
Singular -q /tmp/case1-rankdrop.sing
```

Expected markers include:

```text
CASE1_RECONSTRUCTION_PASS
CASE1_SPECIAL_FIBRE_PYTHON_PASS
CASE1_SPECIAL_FIBRE_BEZOUT_PASS
CASE1_N_MOD_H_PASS
RANKDROP_EXACT_UNIT_PASS
```

For a compact coefficient-tracking regression at the pinned good prime `71`,
reconstruct and verify the specialization of (5.4):

```bash
python scripts/research_jc72_108_case1_rankdrop.py \
  /tmp/case1-rankdrop-p71.sing --prime 71 --lift \
  --write-certificate /tmp/case1-h-certificate-p71.txt
Singular -q /tmp/case1-rankdrop-p71.sing
```

The final markers are:

```text
RANKDROP_UNIT_IDENTITY_PASS
CASE1_H_DETERMINANTAL_PASS
```

This finite-field identity is an explicit regression for the syzygy assembly,
not the characteristic-zero proof. The exact proof is the unit-ideal marker
`RANKDROP_EXACT_UNIT_PASS` together with the ring-theoretic lemma in Section
5.1. An exact coefficient-tracking `liftstd` run is optional and substantially
more expensive; no such coefficient dump is claimed or required here.
