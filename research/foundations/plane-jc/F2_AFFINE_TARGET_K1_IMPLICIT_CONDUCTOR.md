# F2 `k=1` implicit-quintic and conductor-gradient theorem

> **Status.**  The normalized `k=1` parametrization
> `p=t^3+a*t`, `q=t^5+b*t^4+c*t^2+d*t` has an exact implicit quintic with
> only twelve `(P,Q)` support positions and top homogeneous form `P^5`.
> Pulling its two partial derivatives back to the normalization gives
> `F_P(p,q)=q'(t)C(t)` and `F_Q(p,q)=-p'(t)C(t)`, where the degree-eight
> polynomial `C` is the resultant of the collision quartic with the quadratic
> defining each collision pair.  On the generic four-node locus, `C` is
> exactly the eight-point conductor preimage.  This supplies the explicit
> target equation and conductor ideal needed by the source pullback
> compiler; it does not factor that pullback or locate the missing affine
> boundary divisor.

The resultant, support, gradient factorization, and exact generic witness are
replayed by
[`verify_f2_affine_target_k1_implicit_conductor.py`](../scripts/verify_f2_affine_target_k1_implicit_conductor.py).

## 1. Exact implicit equation

Use the normal form from
[`F2_AFFINE_TARGET_K1_COLLISION.md`](F2_AFFINE_TARGET_K1_COLLISION.md):

\[
 p=t^3+at,\qquad q=t^5+bt^4+ct^2+dt.             \tag{1.1}
\]

Define

\[
 F(P,Q)=-\operatorname{Res}_t(P-p(t),Q-q(t)).     \tag{1.2}
\]

Direct elimination gives

\[
\begin{aligned}
F={}&P^5+(ab+b^3+3c)P^4+3bP^3Q
 +(-abc+4ad+3b^2d+3c^2)P^3\\
&+(-5a^2-4ab^2+3bc+3d)P^2Q\\
&+(a^3c+a^2b^2c-a^2bd-2abc^2+5acd+3bd^2+c^3)P^2\\
&-5aPQ^2+(a^3b-3a^2c-5abd+3cd)PQ\\
&+(a^4d+a^3b^2d-2a^2bcd+2a^2d^2+ac^2d+d^3)P\\
&-Q^3+(2a^2b-2ac)Q^2\\
&+(-a^5-a^4b^2+2a^3bc-2a^3d-a^2c^2-ad^2)Q.
                                                               \tag{1.3}
\end{aligned}
\]

Its support is

\[
\begin{split}
&(5,0),(4,0),(3,1),(3,0),(2,1),(2,0),\\
&(1,2),(1,1),(1,0),(0,3),(0,2),(0,1).            \tag{1.4}
\end{split}
\]

Thus `deg_P F=5`, `deg_Q F=3`, `deg F=5`, and the degree-five homogeneous
part is exactly `P^5`, recovering the atlas constraint without a generic
implicitization black box.

## 2. Collision quartic becomes the conductor polynomial

Recall the collision quartic

\[
 R(u)=u^4+bu^3+au^2+(2ab-c)u-(a^2+d).            \tag{2.1}
\]

For a root `u`, the two normalization parameters are the roots of

\[
 z^2-uz+(u^2+a)=0.                               \tag{2.2}
\]

Consequently the polynomial cutting out every normalization preimage of a
collision is

\[
 \boxed{
 C(t)=\operatorname{Res}_u
 \left(R(u),t^2-ut+(u^2+a)\right),\qquad\deg C=8.} \tag{2.3}
\]

The elimination identity (2.3) packages the four unordered collision pairs
without adjoining the quartic roots.

## 3. Exact gradient/conductor identity

Differentiating (1.3), substituting (1.1), and factoring gives

\[
 \boxed{
 F_P(p(t),q(t))=q'(t)C(t),\qquad
 F_Q(p(t),q(t))=-p'(t)C(t).}                     \tag{3.1}
\]

The chain-rule identity `F_Pp'+F_Qq'=0` is then automatic.  On the immersed
normalization locus, `(p',q')` is the unit ideal, so the pullback of the
target Jacobian ideal `(F_P,F_Q)` is the principal conductor divisor `(C)`.
On the generic four-node packet, `C` has the eight distinct normalization
points lying over the nodes.

For the exact witness `(a,b,c,d)=(1,0,0,0)`,

\[
 C=t^8+3t^6+4t^4+2t^2+1,qquad
 \operatorname{Disc}(C)=4{,}000{,}000,            \tag{3.2}
\]

so the eight-point conductor is visibly reduced on a nonempty open set.

## 4. Compiler consequence

The first target-side pullback object is now explicit.  A source compiler can
apply the inverse affine target normalization to (1.3), substitute the F2
Laurent map, and factor one twelve-support expression.  Equation (3.1)
simultaneously transports the four nodal conductor pairs and identifies the
point set where the target normalization ceases to be locally one-to-one.

The subsequent
[`fixed-coordinate Keller-pullback theorem`](F2_AFFINE_K1_KELLER_PULLBACK.md)
restores the target normalization, computes the carrier jet, proves that the
affine pullback is reduced, and identifies its entire affine conductor with
four explicit node fibers.

<!-- status-consumer: PF2K1PB1 6f837229017243c4 -->

The later
[`conductor-conservation theorem`](F2_AFFINE_TARGET_K1_CONDUCTOR_CONSERVATION.md)
uses the adjoint differential `-dt/C` to prove that this same polynomial is
the exact conductor divisor on every degenerate `k=1` stratum, not just on
the immersed four-node open set.

<!-- status-consumer: PF2K1CC1 f152c82ef2d54c32 -->

What remains is genuinely boundary-side: supply the complete fixed-coordinate
F2 Laurent pair, factor or value `F(P(x,y),Q(x,y))` at the unresolved source
boundary, identify the boundary divisor with `e>1`, and compute its
`(e,f,E^2,b)` and boundary-local nonunit-`Fitt_1` data.  This theorem does not
perform that factorization, exclude `(75,125)`, or prove `JC(2)`.

## Reproduction

```bash
.venv/bin/python scripts/verify_f2_affine_target_k1_implicit_conductor.py
```
