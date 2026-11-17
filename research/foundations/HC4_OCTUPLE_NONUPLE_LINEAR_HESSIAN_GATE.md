# The octuple/nonuple linear Hessian gate

## Status

This note proves `HC4NHM3`: on the generic-corank-one ternary-quintic
Schur stratum, no packet with an essential line and

\[
\det\operatorname{Hess}(h_5)=x^8\ell
\quad\text{or}\quad
\det\operatorname{Hess}(h_5)=x^9
\]

exists. Here \(\ell\) is a line not proportional to \(x\). Together with
`HC4NHM2`, this closes the one-component partition \(P=x^4\) of the clean
quartic-denominator packet from `HC4NHM1`. It does **not** treat the other
quartic partitions, generic lower-Smith components, or unrestricted
`HC(4)`.

Replay the exact identities with

~~~bash
.venv/bin/python scripts/verify_hc4_octuple_nonuple_linear_hessian_gate.py
~~~

The normalization input is
[`HC4_NONREDUCED_HESSIAN_SCHUR_MODULE.md`](HC4_NONREDUCED_HESSIAN_SCHUR_MODULE.md).
The moving-line, conic-kernel, and constant-kernel coefficient gates used
below are proved in
[`HC4_DIRECT_SEPTUPLE_LINEAR_HESSIAN_GATE.md`](HC4_DIRECT_SEPTUPLE_LINEAR_HESSIAN_GATE.md).

Work over an algebraically closed field of characteristic zero. Put

\[
C=\operatorname{Hess}(h_5),\qquad d=\nabla s_3,
\tag{0.1}
\]

assume

\[
\det C=x^mR_{9-m},\qquad m\in\{8,9\},
\tag{0.2}
\]

where \(R_1=\ell\) is coprime to \(x\) for \(m=8\), and \(R_0=1\)
for \(m=9\). Assume that \(C\bmod x\) has generic rank two, that \(x\)
is essential in the minimal denominator clearing \(C^{-1}d\), and that

\[
\det C\mid d^{\mathsf T}\operatorname{adj}(C)d.
\tag{0.3}
\]

## 1. The ten-row ladder

At the generic point of \(x=0\), the Smith form is
\(\operatorname{diag}(1,1,x^m)\). The Schur condition bounds the pole order
\(b\) of \(C^{-1}d\) by

\[
1\le b\le\left\lfloor\frac m2\right\rfloor=4.
\tag{1.1}
\]

For \(m=8\), the simple residual line cannot be an essential denominator
component; for \(m=9\), there is no residual component. Thus the minimal
denominator is \(P=x^b\), and

\[
e=x^bC^{-1}d\in K[x,y,z]_{b-1}^3.
\tag{1.2}
\]

Let the saturated kernel on \(L=(x=0)\) be
\(\mathcal K_L=\mathcal O_L(-\kappa)\), and let \(B_L\) be its normalized
corank-two defect. The self-dual determinant formula gives

\[
\deg B_L=6-2\kappa.
\tag{1.3}
\]

The nonzero residue of \(e\) is a section of
\(\mathcal O_L(b-1-\kappa)\), while \(\deg B_L\ge0\). Hence

\[
0\le\kappa\le\min(b-1,3).
\tag{1.4}
\]

For either value of \(m\), the complete numerical ladder is therefore

| pole \(b\) | primitive kernel degrees \(\kappa\) | defect lengths |
|---:|---:|---:|
| \(1\) | \(0\) | \(6\) |
| \(2\) | \(0,1\) | \(6,4\) |
| \(3\) | \(0,1,2\) | \(6,4,2\) |
| \(4\) | \(0,1,2,3\) | \(6,4,2,0\) |

There are ten rows. The only row not already present in the septuple ladder
has

\[
(b,\kappa,\deg B_L)=(4,3,0).
\tag{1.5}
\]

## 2. The rows of kernel degree one and two

Every primitive degree-one kernel has projective-line image. A primitive
degree-two kernel has either projective-line image or a nondegenerate conic
image. Section 3 of `HC4NHM2` proves that a moving line image is incompatible
with a repeated Hessian line. Sections 4--6 prove that every conic image
fails the first-normal divisibility needed even for \(x^2\mid\det C\).
Those proofs use only repeatedness of the line, not exact multiplicity seven.
Consequently every row with \(\kappa=1\) or \(2\) is empty for (0.2).

## 3. The defect-free cubic-kernel row

Let \(e|_L=(a,b,c)\) be the primitive kernel on (1.5). If its coordinates
span a plane, its projective image is a line and the preceding moving-line
argument applies. Otherwise they span all of \(K^3\). The basepoint-free
map \(\mathbf P^1\to\mathbf P^2\) has degree three and nondegenerate image,
so it is a rational plane cubic. In particular, \(a\ne0\).

Use the boundary expansion

\[
h_5=F_5+xG_4+\frac{x^2}{2}H_3+\cdots.
\tag{3.1}
\]

Since the defect is zero, symmetry gives a nonzero scalar \(\lambda\) with

\[
\operatorname{adj}(C|_L)=\lambda ee^{\mathsf T}.
\tag{3.2}
\]

The \(xx\)-cofactor is therefore

\[
\boxed{\det\operatorname{Hess}(F)=\lambda a^2.}
\tag{3.3}
\]

Thus the binary quintic Hessian is a perfect square of a cubic. The complete
square-divisor classification in `HC4NHM2` immediately makes this finite.
The squarefree orbit is represented by \(u^5+v^5\), whose Hessian
\(400u^3v^3\) is not a square. The exceptional \(2+1+1+1\) orbit has

\[
-600v^2(u-2v)^2(u^2+4uv+6v^2),
\tag{3.4}
\]

again not a square. The normalized \(3+1+1\) and \(2+2+1\) forms have
Hessians

\[
\begin{aligned}
\operatorname{HessDet}(u^3v(u-v))
 &=-8u^4(2u^2-3uv+3v^2),\\
\operatorname{HessDet}(u^2v^2(u-v))
 &=-8u^2v^2(3u^2-4uv+3v^2),
\end{aligned}
\tag{3.5}
\]

whose residual quadratics are squarefree. The pure fifth power has zero
binary Hessian. Hence only the root types \(4+1\) and \(3+2\) remain.

### 3.1 Root type \(4+1\)

Normalize \(F=uv^4\), so \(a=v^3\). Solving the last two rows of
\((C|_L)e=0\) first forces the \(u^4\)-coefficient of \(G\) to vanish.
Write

\[
G=g_0v^4+g_1uv^3+g_2u^2v^2+g_3u^3v.
\tag{3.6}
\]

Then the complete solution is

\[
\begin{aligned}
b&=-g_0v^3+g_2u^2v+2g_3u^3,\\
c&=-\frac v4(g_1v^2+2g_2uv+3g_3u^2).
\end{aligned}
\tag{3.7}
\]

The first row requires
\(H=-(G_ub+G_vc)/v^3\) to be polynomial. Modulo \(v^3\), its numerator is

\[
G_ub+G_vc
\equiv \frac{21}{4}g_3^2u^5v+5g_2g_3u^4v^2.
\tag{3.8}
\]

Characteristic zero forces \(g_3=0\). Equations (3.7) then show that
\(a,b,c\) share the factor \(v\), contradicting primitivity.

### 3.2 Root type \(3+2\)

Normalize \(F=u^2v^3\), so \(a=uv^2\). Polynomiality in the last two
kernel equations forces

\[
G=g_1uv^3+g_2u^2v^2+g_3u^3v.
\tag{3.9}
\]

They then give

\[
\begin{aligned}
b&=\frac{-g_1uv^2+g_3u^3}{2},\\
c&=-\frac{g_2uv^2+2g_3u^2v}{3}.
\end{aligned}
\tag{3.10}
\]

All of \(a,b,c\) share the factor \(u\). This again contradicts the
saturated primitive kernel. The defect-free cubic row is empty.

## 4. The constant-kernel rows

Only \(\kappa=0\) remains. A constant kernel transverse to \(L\) is
normalized to \(\partial_x\). The boundary equations give \(G=H=0\), and
successive determinant faces kill \(J,K,c\); the complete determinant is
then zero. This was (7.1)--(7.3) of `HC4NHM2` and already follows from
divisibility by \(x^7\).

For a tangent kernel, normalize to \(\partial_z\) and write

\[
\begin{aligned}
F&=\alpha y^5,&G&=\beta y^4,\\
H&=h_0y^3+h_1y^2z+h_2yz^2+h_3z^3,\\
J&=j_0y^2+j_1yz+j_2z^2,&K&=k_0y+k_1z.
\end{aligned}
\tag{4.1}
\]

If \(\alpha\ne0\), divisibility by \(x^7\) successively forces

\[
h_1=h_2=h_3=j_1=j_2=k_1=0.
\tag{4.2}
\]

If \(\alpha=0\), generic boundary rank gives \(\beta\ne0\), and
divisibility by \(x^8\) forces the same vanishing (the last step uses
\([x^7]\det C=-\beta k_1^2y^2/3\)). In either case every surviving term of
\(h_5\) is independent of \(z\). Hence \(C\) has a zero row and column and

\[
\det C=0.
\tag{4.3}
\]

All four constant-kernel rows are empty.

## 5. Conclusion and historical next packet

Sections 2--4 eliminate all ten rows for both multiplicities.

> **Theorem `HC4NHM3` -- Octuple/nonuple linear Hessian--Schur exclusion.**
> Under (0.1)--(0.3), with \(C\bmod x\) of generic rank two and \(x\)
> essential for the minimal Schur denominator, no packet with
> \(\det C=x^8\ell\), \(x\nmid\ell\), or with \(\det C=x^9\), exists.

As in `HC4NHM2`, the Hessian boundary is eliminated before curl-freeness or
the scalar cleared Schur equation is needed. In the clean
\(\det C=P^2\ell\), \(\deg P=4\), branch of `HC4NHM1`, this closes the
partition

\[
P=x^4.
\tag{5.1}
\]

The next partition is \(P=x^3y\). Since
\(\det C=P^2\ell\), the residual line has three incidence types:

\[
\det C=x^7y^2,\qquad x^6y^3,\qquad x^6y^2z,
\tag{5.2}
\]

after normalizing a distinct residual line to \(z=0\). The first is exactly
the nonsquarefree-cofactor septuple line left outside `HC4NHM2`. In all three,
the denominator has support on both repeated lines, so the residue degree on
either normalization is global rather than determined by one local pole;
that two-component compatibility was the next finite gate. The first incidence
type is classified explicitly in
[`HC4_TWO_LINE_QUARTIC_DENOMINATOR_PACKET.md`](HC4_TWO_LINE_QUARTIC_DENOMINATOR_PACKET.md).
`HC4NHM5--6` subsequently close every `3+1` incidence, `HC4NHM7--8` close
the `2+2` and `2+1+1` partitions, and `HC4NHM9--12` close the squarefree
partition. Thus the split-linear clean denominator branch is closed. Clean
nonlinear components, positive-defect packets, and lower-Smith components
remain separate.
