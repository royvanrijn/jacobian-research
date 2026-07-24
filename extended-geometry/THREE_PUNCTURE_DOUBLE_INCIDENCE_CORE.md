# Three-puncture double-incidence core

This note isolates a polynomial determinant mechanism with two interacting
inverse equations.  It is the first core in the repository whose selected
critical complete intersection has normalization

\[
\mathbb P^1\setminus\{0,1,\infty\}.
\]

It is not yet a Keller map.  The result is a sharp pre-algebraization:
the determinant ledger and the three-puncture normalization are both exact,
and the remaining affine-space problem is reduced to a coupled Jacobian
completion which cannot retain all three base coordinates.

Work over a characteristic-zero field \(k\).  Fix integers \(a,b\geq1\).

## 1. Two inverse equations and their diagonal primitive

On \(\mathbb A^4_{s,t,u,v}\), put

\[
D_0=1-su,\qquad D_1=1-tv,\qquad \ell=t-s.
\]

The vector field

\[
\partial_\Delta=\partial_s+\partial_t
\]

fixes \(\ell,u,v\).  Define

\[
\begin{aligned}
R_{a,b}(s,t,u,v)
&=\int_0^s
 (1-\xi u)^a(1-(\xi+\ell)v)^b\,d\xi\\
&=\sum_{i=0}^a\sum_{j=0}^b\sum_{q=0}^j
 \binom ai\binom bj\binom jq
 \frac{(-u)^i(-v)^j\ell^{j-q}s^{i+q+1}}
 {i+q+1}.
\end{aligned}                                      \tag{1}
\]

This is a polynomial over \(k\).  The fundamental theorem of calculus,
with \(\ell\) fixed, gives

\[
\boxed{\partial_\Delta R_{a,b}=D_0^aD_1^b.}         \tag{2}
\]

Consequently the polynomial map

\[
\boxed{
\Phi_{a,b}:\mathbb A^4\longrightarrow\mathbb A^4,
\quad
(s,t,u,v)\longmapsto
(t-s+1,u,v,R_{a,b})
}                                                    \tag{3}
\]

has

\[
\boxed{\operatorname {Jac}(\Phi_{a,b})=-D_0^aD_1^b} \tag{4}
\]

for the displayed source and target order.  Thus the complete determinant
ledger has two independently selected prime columns, with multiplicities
\(a\) and \(b\).  Unlike a one-root suspension, neither factor is a power
of the other.

## 2. The selected critical curve has three punctures

Let

\[
C=V(D_0,D_1,t-s+1)\subset\mathbb A^4.               \tag{5}
\]

On \(C\),

\[
t=s-1,\qquad u=s^{-1},\qquad v=(s-1)^{-1}.          \tag{6}
\]

Therefore

\[
\boxed{
k[C]\simeq k[s,s^{-1},(s-1)^{-1}]
}                                                    \tag{7}
\]

and

\[
\boxed{
C\simeq\mathbb P^1\setminus\{0,1,\infty\}.
}                                                    \tag{8}
\]

The gradients of \(D_0,D_1,t-s+1\) have rank three along \(C\), so this is
a smooth complete intersection, not a puncture count produced by
normalizing a singular plane model.  Its geometric unit lattice is

\[
\mathcal O(C)^*/\overline k^*
\simeq\mathbb Z\,[s]\oplus\mathbb Z\,[s-1],          \tag{9}
\]

of rank two.  The two determinant components select the two primitive
finite-puncture characters separately.

This distinguishes the core from the established marked
\(\mathbb A^1\) and \(\mathbb G_m\) cores before any compactification or
conductor decoration is added.

## 3. Exact reciprocal suspension and why it is not affine space

On

\[
U=\mathbb A^4\setminus V(D_0D_1)
\]

adjoin \(z\) and set

\[
Z=\frac{z}{D_0^aD_1^b}.
\]

Then

\[
(s,t,u,v,z)\longmapsto
(t-s+1,u,v,R_{a,b},Z)                               \tag{10}
\]

has constant determinant \(-1\).  This proves that the ledger is balanced
on the boundary complement.

It is not the requested polynomial algebraization.  The source
\(U\times\mathbb A^1\) has the two independent units \(D_0,D_1\), whereas
affine space has no nonconstant units.  Polynomial stabilization does not
remove this obstruction.

## 4. The one-coordinate affine completion gate

Suppose an attempted five-variable completion retains

\[
L=t-s+1,\qquad u,\qquad v
\]

as three target coordinates.  Use \(r=s\), so \(t=r+L-1\).  The remaining
two outputs are polynomials \(A(r,z),B(r,z)\) over

\[
K=k(L,u,v).
\]

The full Keller equation is, up to sign,

\[
\det\frac{\partial(A,B)}{\partial(r,z)}\in K^*.      \tag{11}
\]

Thus every such completion is already a plane Keller problem over \(K\).
The extra determinant coordinate has not used the three-puncture geometry;
it has merely moved the unresolved step into dimension two.

There is an unconditional conclusion for the first bounded ansatz.  If

\[
A=A_0(r)+a(r)z,\qquad B=B_0(r)+b(r)z,                \tag{12}
\]

then the coefficient of \(z\) in (11) is

\[
ba'-ab'.                                            \tag{13}
\]

If (11) is a nonzero constant, (13) vanishes.  Over \(K(r)\), \(a/b\) is
constant.  After a constant target row operation, (11) becomes

\[
b(r)H'(r)\in K^*.
\]

Both factors are therefore units of \(K[r]\).  A second affine target
change recovers \(r\), and then \(z\).  Hence:

\[
\boxed{\text{Every affine-linear one-coordinate completion retaining
\(L,u,v\) is a polynomial automorphism.}}            \tag{14}
\]

The same argument includes the cases \(a=0\) or \(b=0\).

## 5. Consequence for the affine-space search

Equations (1)--(9) answer the determinant and normalization halves of the
research question positively.  Equations (10)--(14) leave the
affine-space half open and impose a sharp gate:

1. at least one of \(L,u,v\) must be replaced, not retained;
2. one extra coordinate with affine-linear coupling is insufficient;
3. a surviving map must use a genuinely nonlinear three-variable block,
   or two coupled modification coordinates;
4. on the slice defining \(C\), the two boundary characters must remain
   independently selected, or the construction descends back to a
   \(\mathbb G_m\) quotient.

This is structurally parallel to the translated Davenport determinant
frontier, but it does not use the Davenport polynomial, a Gassmann pair, a
marked root, or conductor gluing.  It supplies a smaller universal test
core for the proposed dimension-four surgery.

The next finite search should replace \(u\) (or \(v\)) by an output
nonlinear in \((r,u,z)\), impose constant determinant coefficient by
coefficient, and reject every solution whose target algebra recovers
\((r,z)\).  Retaining all of \(L,u,v\) should not be searched further.

## 6. Verification

Run

```bash
.venv/bin/python scripts/verify_three_puncture_double_incidence.py
```

The checker verifies (1)--(4) for \(1\leq a,b\leq4\), the complete
intersection rank, the normalization parametrization, the two independent
puncture units, and the affine-linear completion identity.
