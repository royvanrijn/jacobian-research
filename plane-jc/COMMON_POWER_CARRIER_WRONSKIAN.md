# Primitive common-power carriers and universal Hurwitz passports

## Result and scope

The exact checker
[`cas/verify_common_power_carrier_wronskian.py`](cas/verify_common_power_carrier_wronskian.py)
extracts the carrier-Wronskian argument from the F2 `(75,125)` calculation and
proves it for every primitive common-power edge in characteristic zero.

Let

\[
q=y,\qquad v=xy^k,\qquad x=vq^{-k},\qquad k\geq 3,
\]

so that

\[
dx\wedge dy=-q^{-k}\,dq\wedge dv. \tag{1}
\]

Suppose a constant-Jacobian pair has leading carrier edge

\[
P=a q^{-km}c(v)^m+\cdots,
\qquad
R=b q^{-kn}c(v)^n+\cdots, \tag{2}
\]

where `ab != 0`, `gcd(m,n)=1`, `deg(c)=k`, and the roots of `c` have
multiplicities

\[
\mu_1+\cdots+\mu_r=k,
\qquad \gcd(\mu_1,\ldots,\mu_r)=1. \tag{3}
\]

The last condition is the **primitive carrier hypothesis**.  Then ordinary
target shears put the first nonconstant, nonshear target coefficient at the
forced descent

\[
\boxed{\delta=k(m+n-1)+1}. \tag{4}
\]

If this coefficient is `H(v)`, then necessarily

\[
\boxed{
H(v)=\frac{\operatorname{rad}(c)(v)N(v)}{c(v)^{m+n}},
\qquad \deg N=k-r-1.
} \tag{5}
\]

Writing `D=rad(c)N`, the entire first Jacobian row is the low-degree equation

\[
\boxed{
kD'-(k-1)\frac{c'}cD=\kappa,
\qquad \kappa\ne0.
} \tag{6}
\]

It is independent of the common powers `(m,n)`.  Those powers determine only
the target monomial, the number of removable shears, and the descent (4).

The carrier residue map is, up to nonzero target scaling,

\[
\boxed{g(v)=\frac{D(v)^k}{c(v)^{k-1}}}. \tag{7}
\]

Equation (6) gives

\[
D\frac{g'}g=\kappa. \tag{8}
\]

Thus `g` has no critical points away from its zeros, poles, and `v=infinity`.
It is an explicit three-point map.  This converts every primitive carrier
solution into a finite Hurwitz problem before any lower Laurent bands are
considered.

This is a necessary reduction, not an existence theorem for every passport
below and not a realization of a global Keller map.

The pinned exact output is
[`../artifacts/generated-results/jc2_common_power_carrier_wronskian.json`](../artifacts/generated-results/jc2_common_power_carrier_wronskian.json).

## 1. Unimodular target coordinates

Choose nonnegative integers `A,B` satisfying

\[
Bn-Am=1. \tag{9}
\]

They exist because `gcd(m,n)=1`.  Set

\[
\pi=P^A/R^B,
\qquad
h=P^n/R^m. \tag{10}
\]

The logarithmic exponent matrix is

\[
\begin{pmatrix}A&-B\\ n&-m\end{pmatrix},
\qquad
\det=Bn-Am=1. \tag{11}
\]

Consequently

\[
d\pi\wedge dh=\frac{\pi h}{PR}\,dP\wedge dR. \tag{12}
\]

At the carrier edge,

\[
\pi=q^kU(v)+\cdots,
\qquad U(v)=u/c(v),\quad u\ne0,
\qquad h=h_0+\cdots. \tag{13}
\]

The right side of (12) has `q`-order `k(m+n)`.  If the first normalized
coefficient of `w=h-h_0-shears` is `q^delta H(v)`, then

\[
d(q^kU)\wedge d(q^\delta H)
=q^{k+\delta-1}(kUH'-\delta U'H)\,dq\wedge dv. \tag{14}
\]

Equating orders proves (4).  In particular

\[
\gcd(k,\delta)=1, \tag{15}
\]

so the selected target ray `(k,delta)` is primitive and the transverse
ramification index of the generic carrier is one.

## 2. Why every earlier coefficient is a target shear

Suppose the first surviving coefficient before (4) occurs at descent `d`.
The corresponding homogeneous two-form row is

\[
kUH_d'-dU'H_d=0, \tag{16}
\]

hence

\[
H_d^k/U^d\in K^\times. \tag{17}
\]

At a root of `c` of multiplicity `mu_i`, rationality of `H_d` requires
`k | d*mu_i`.  The primitive condition (3) then implies `k | d`.  If
`d=jk`, equation (16) gives `H_d=constant*U^j`, exactly the coefficient
removed by the target shear `h -> h-constant*pi^j`.

Therefore the complete pre-Jacobian list is

\[
k,2k,\ldots,k(m+n-1). \tag{18}
\]

There are `m+n-1` such shears, and the next descent is (4), which is not
divisible by `k`.

## 3. The forced divisor of the Jacobian coefficient

At descent `delta`, equation (14) has a nonzero right side proportional to
`c^{-(m+n+1)}`.  Dividing by `U=u/c` gives

\[
kH'+\delta\frac{c'}cH=Kc^{-(m+n)},
\qquad K\ne0. \tag{19}
\]

At a root `alpha` of `c` of multiplicity `mu`, let `ell=ord_alpha(H)`.
The leading order and coefficient in (19) are

\[
\ell=1-\mu(m+n),
\qquad
k\ell+\delta\mu=k-\mu(k-1). \tag{20}
\]

For `k>=3` the last integer is never zero.  Thus every finite order in (20)
is exact.  At a finite point outside `c=0`, a pole of `H` would give an
uncancelled pole in `H'`, so no other finite poles occur.  It follows that

\[
H=\operatorname{rad}(c)N/c^{m+n} \tag{21}
\]

for a polynomial `N` coprime to `c`.

At infinity the two possible local orders are

\[
\delta
\quad\hbox{and}\quad
k(m+n)-1. \tag{22}
\]

Primitivity implies `r>=2`.  The second order in (22) would require
`deg(N)=1-r<0`, so it cannot occur.  The first gives

\[
\deg(N)=k-r-1. \tag{23}
\]

In particular `r<=k-1`: a squarefree degree-`k` carrier polynomial is
impossible.

Substitution of (21) into (19), using (4), cancels every occurrence of
`m+n` and gives (6).  If `N` had a root of multiplicity `nu>1` away from
`c=0`, the left side of (6) would vanish to order `nu-1`, contradicting
`kappa!=0`.  Hence `N` is squarefree as well as coprime to `c`.

## 4. The fixed-carrier problem is linear

For a reconstructed candidate `c`, multiply (6) by `c`:

\[
\boxed{
\Phi_c(D,\kappa)
=kcD'-(k-1)c'D-\kappa c=0,
\qquad \deg D\le k-1.
} \tag{24}
\]

This is a homogeneous linear coefficient system in the `k` coefficients of
`D` and the scalar `kappa`.  It is the reusable kernel/cokernel interface for
carrier reconstruction.

Moreover the projective solution is unique whenever it exists.  Indeed, a
nonzero element of the homogeneous kernel

\[
kcD'-(k-1)c'D=0 \tag{25}
\]

would make `D^k/c^(k-1)` constant.  At every root of multiplicity `mu_i`
this requires

\[
k\operatorname{ord}_{\alpha_i}(D)=(k-1)\mu_i. \tag{26}
\]

Since `gcd(k,k-1)=1`, every `mu_i` would be divisible by `k`.  Their sum is
`k`, forcing the one-part partition `(k)`, contrary to primitivity.  Thus
(25) has zero kernel, and any two nonzero solutions of (24) are proportional.

The checker builds this matrix exactly.  For

\[
c=v(v-1)^2(v^2-3v+3), \tag{27}
\]

the `9 x 6` matrix has rank five and its one-dimensional kernel gives

\[
D=v(v-1)(v^2-3v+3) \tag{28}
\]

up to scale.  Replacing the last two roots by the generic pair `2,3` raises
the rank to six and kills the kernel.  Thus admissibility is already a
determinantal condition on `c`, while reconstruction of `D` is linear.

## 5. Universal three-point passport

Let

\[
s=\#\{i:\mu_i=1\},
\qquad L=k-r-1,
\qquad f=kL+s. \tag{29}
\]

At a simple root of `c`, equation (7) has a simple zero.  At a root of `N`,
it has a zero of order `k`.  At a root of `c` with `mu_i>=2`, it has a pole
of order

\[
p_i=(k-1)\mu_i-k. \tag{30}
\]

These positive integers sum to `f`.  Since `deg(D)=k-1`, numerator and
denominator in (7) have the same degree, and `g(infinity)` is finite and
nonzero.  Equation (8) gives

\[
g(v)-g(\infty)\sim \text{constant}\,v^{-(k-2)}, \tag{31}
\]

so the ramification index at infinity is `k-2`.  Equation (8) also proves
that every other point in that fiber is unramified.

After scaling the third value to one, the universal passport is therefore

\[
\boxed{
(k^L,1^s)
\;\bigm|\;
(p_i:\mu_i\ge2)
\;\bigm|\;
(k-2,1^{f-k+2}).
} \tag{32}
\]

The Riemann-Hurwitz identity is

\[
L(k-1)+\sum_{\mu_i\ge2}(p_i-1)+(k-3)=2f-2. \tag{33}
\]

Thus every solution of (24) determines a genus-zero three-point Hurwitz
class with no unidentified branch values.

### Geometric monodromy parity

The sign of a branch permutation is `(-1)` to the ramification contribution
of its fiber.  If `k` is odd, then `k-1` and `k-3` are even.  The first and
third permutations in (32) are therefore even, and (33) forces the pole
permutation to be even as well.  Hence

\[
\boxed{k\text{ odd}\quad\Longrightarrow\quad G_{\rm geom}\subseteq A_f.}
\tag{34}
\]

If `k` is even, the distinguished cycle of length `k-2` has odd sign, so

\[
\boxed{k\text{ even}\quad\Longrightarrow\quad G_{\rm geom}\not\subseteq A_f.}
\tag{35}
\]

This is a geometric statement.  Galois action on the constant field can
still enlarge the arithmetic group, as happens for the F2 degree-six map
with geometric group `A_6` and arithmetic group `S_6`.

### Simple/double corollary

If the partition has `t` double roots and `s=k-2t>=1` simple roots, then it
is primitive and

\[
L=t-1,
\qquad f=t(k-2), \tag{36}
\]

while (32) becomes

\[
\boxed{
(k^{t-1},1^{k-2t})
\mid ((k-2)^t)
\mid (k-2,1^{(t-1)(k-2)}).
} \tag{37}
\]

For `t=1`, the degree is `k-2` and the two nontrivial fibers are totally
ramified.  Hence the map is cyclic after source and target normalization.

## 6. Recovery of the F2 carrier rows

The F2 calculation has `(k,m,n)=(5,3,5)` and therefore

\[
\delta=5(3+5-1)+1=36. \tag{38}
\]

Its two root partitions are

\[
(2,1,1,1)
\quad\hbox{and}\quad
(2,2,1). \tag{39}
\]

Formula (37) gives respectively

\[
1^3\mid3\mid3,
\qquad
(5,1)\mid(3,3)\mid(3,1,1,1). \tag{40}
\]

These are exactly the cyclic cubic and degree-six Belyi maps certified in
[`F2_75_125_CARRIER_WRONSKIAN_CLASSIFIER.md`](F2_75_125_CARRIER_WRONSKIAN_CLASSIFIER.md).
The new theorem explains why those passports were forced before their
coefficients were solved.

## 7. Explicit exceptional loci

Two loci are deliberately outside the theorem.

1. **`k=2`.**  A double root makes the coefficient in (20) vanish, and the
   infinity equation is resonant.  It needs a separate logarithmic analysis.
2. **Imprimitive root partition.**  If
   `a=gcd(mu_1,...,mu_r)>1`, then `c=c_0^a` over an algebraic closure and
   rational homogeneous coefficients may occur already at descent `k/a`.
   They are not ordinary integral-power shears in `pi`, whose descent is
   `k`.  This is a genuinely ramified transverse problem, not a harmless
   omission.

The exact checker enumerates every multiplicity partition through `k=24` as
a regression audit, verifies all primitive passports and Riemann-Hurwitz
identities, and labels every imprimitive partition as deferred.  The theorem
itself is the symbolic argument above and is not bounded by that census.

## 8. Implementation consequence

For any common-power Laurent branch satisfying (1)--(3), the appropriate
pipeline is now:

1. enumerate the finite passport (32), or its relevant constrained
   subpassport;
2. reconstruct candidate three-point maps over their number fields;
3. recover `c` from the prescribed pole multiplicities;
4. solve the single linear system (24) for `(D,kappa)`;
5. impose distinguished carrier points and field-of-moduli compatibility;
6. only then solve the successive lower Laurent deformation maps.

No simultaneous Groebner elimination of the carrier and lower bands is
needed.  Failure of the kernel in (24), of the marked-point normalization,
or of Galois descent kills a dessin before the lower deformation calculation.
