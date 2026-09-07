# Incomplete commuting polynomial flows

Write \(F=(F_1,F_2,F_3)\) for the map in FACTS.md, and let

\[
X_i=(DF)^{-1}e_i.
\]

Since \(\det DF=-2\), the entries of \((DF)^{-1}\) are polynomial. Moreover
\(X_i(F_j)=\delta_{ij}\). The \(X_i\) frame the tangent bundle. Applying a
commutator to every \(F_k\) gives \([X_i,X_j](F_k)=0\); invertibility of
\(DF\) then gives \([X_i,X_j]=0\).

## An explicit incomplete integral curve

The square roots in this formula are essential. Starting at
\(p=(1,-3/2,13/2)\), choose the branch with \(\sqrt{1}=1\) and put

\[
\rho=\sqrt{1-4s},\qquad
q(s)=\left(\rho^{-1},-\frac32\rho,\frac{13}{2}\rho^2\right).
\]

Exact substitution gives \(F(q(s))=(-1/4+s,0,0)\). Differentiating and using
the invertibility of \(DF\) gives \(\dot q=X_1(q)\). The first coordinate has
a pole of Puiseux order \(1/2\) at \(s=1/4\). Hence \(X_1\) is not complete.
Putting \(1-4s\), rather than its square root, in the first two coordinates
does not give this identity.

## Classification of constant linear combinations

For \(v=(\alpha,\beta,\gamma)\), set
\(X_v=\alpha X_1+\beta X_2+\gamma X_3\). Along every local integral curve,

\[
F(q(s))=F(q(0))+sv. \tag{1}
\]

The exact image calculation gives

\[
F(\mathbb C^3)=\mathbb C^3\setminus\Gamma,\quad
\Gamma=\{(4/(27c^2),4/(3c),c):c\ne0\}.
\]

If \(X_v\) were complete, (1) would imply that translation by \(sv\) preserves
the image for every \(s\), and therefore preserves \(\Gamma\). But \(\Gamma\)
has no nonzero constant translation direction: its tangent is

\[
g'(c)=(-8/(27c^3),-4/(3c^2),1),
\]

and no fixed nonzero vector is proportional to this for all \(c\). Consequently

\[
\boxed{X_v\text{ is complete if and only if }v=0.}
\]

This also produces an incomplete orbit for every nonzero \(v\): choose
\(g\in\Gamma\) for which the line \(g+\mathbb Cv\) is not contained in
\(\Gamma\), take a nearby \(a=g-s_0v\notin\Gamma\), and lift \(a\). If all
lifts existed through \(s_0\), one would obtain a preimage of \(g\), impossible.

## Escape divisor and generic exponents

In target coordinates \((a,b,c)\), put

\[
Q=27a^2c^2-18abc+16a+b^3c-b^2.
\]

For a direction \(v\), the parameter-space escape divisor is

\[
E_v:\quad Q(a+s\alpha,b+s\beta,c+s\gamma)=0. \tag{2}
\]

This is the pullback of the nonproperness hypersurface. Over a smooth point
of \(V(Q)\), two of the three sheets escape and the remaining sheet has a
finite limit. Use the cubic inverse model

\[
P(T)=cT^3-2T^2+bT-2a,\qquad x=2/P'(t).
\]

At a transverse encounter, a double root splits as
\(t-t_0\asymp(s-s_0)^{1/2}\). Hence
\(P'(t)\asymp(s-s_0)^{1/2}\) and
\(x\asymp(s-s_0)^{-1/2}\). Thus the generic blow-up exponent is \(1/2\),
including the explicit orbit. If the line has contact order \(m\) with the
smooth discriminant and the local unfolding is otherwise nondegenerate, the
exponent becomes \(m/2\).

The singular locus of \(V(Q)\) is \(\Gamma\). There \(P\) has a triple root.
For a generic transverse one-parameter deformation,
\(t-t_0\asymp(s-s_0)^{1/3}\), so
\(P'(t)\asymp(s-s_0)^{2/3}\) and
\(x\asymp(s-s_0)^{-2/3}\). Accordingly the generic exponent at \(\Gamma\)
is \(2/3\).

### Complete classification at the missing curve

Fix \(g(c_0)\in\Gamma\), put \(t_0=2/(3c_0)\), and write \(h=s-s_0\).
After replacing \(t\) by \(t_0+T\), the inverse cubic along the target line is

\[
c_0T^3+h(A_0+A_1T+A_2T^2+\gamma T^3),
\]

where

\[
A_0=\gamma t_0^3+\beta t_0-2\alpha,\quad
A_1=3\gamma t_0^2+\beta,\quad A_2=3\gamma t_0.
\]

This gives all nonzero directions:

| condition | root scale | scale of \(P'\) | blow-up exponent |
|---|---:|---:|---:|
| \(A_0\ne0\) | \(T\asymp h^{1/3}\) | \(h^{2/3}\) | \(2/3\) |
| \(A_0=0,\ A_1\ne0\) | \(T=0\) or \(T\asymp h^{1/2}\) | \(h\) | \(1\) |
| \(A_0=A_1=0\) | \(T\asymp h\) for the simple sheet | \(h^2\) | \(2\) |

The first row is transverse to the tangent plane of the discriminant. The
second lies in that tangent plane but is not tangent to \(\Gamma\). In the last
row, \(v\) is proportional to

\[
g'(c_0)=(-8/(27c_0^3),-4/(3c_0^2),1);
\]

the target line is tangent to \(\Gamma\). The cubic has a persistent double
root, while its unique simple sheet escapes with exponent \(2\). Requiring
\(A_0=A_1=A_2=0\) forces \(v=0\), so the table is exhaustive.

These exponents refer to \(x=2/P'(t)\), hence to the largest affine coordinate
after reconstruction. The table upgrades the earlier generic \(2/3\) statement
to a full directional classification at every point of \(\Gamma\).

Run .venv/bin/python archive/tooling/commuting_flows.py for the exact certificates.

## A nonconstant complete combination

Although no nonzero *constant* combination is complete, the frame does contain
a natural complete combination with polynomial coefficients. The weighted
Euler field on the target is

\[
H=-2a\,\partial_a-b\,\partial_b+c\,\partial_c.
\]

It satisfies \(H(Q)=-2Q\), preserves \(\Gamma\), and has the complete flow

\[
(a,b,c)\longmapsto(e^{-2t}a,e^{-t}b,e^tc).
\]

Its lift through \(F\) is

\[
Y=(-2F_1)X_1-F_2X_2+F_3X_3
  =x\,\partial_x-y\,\partial_y-2z\,\partial_z.
\]

Thus \(Y\) is polynomial and complete, with flow

\[
(x,y,z)\longmapsto(e^tx,e^{-t}y,e^{-2t}z).
\]

This cleanly separates the two phenomena: translational directions in the flat
frame are all incomplete, while a multiplicative symmetry adapted to the
nonproperness geometry is complete.

There is also a small rigidity result. If \(h(q)=Mq+d\) is any affine target
vector field tangent to \(\Gamma\), direct coefficient comparison on

\[
g(c)=(4/(27c^2),4/(3c),c)
\]

shows that \(d=0\) and

\[
M=\mu\,\operatorname{diag}(-2,-1,1).
\]

Hence \(H\) spans all affine target vector fields preserving \(\Gamma\). In
particular, the weighted complete flow is the unique affine candidate, up to
constant rescaling.

## Polynomial target symmetries through degree two

An exact linear calculation gives a 13-dimensional space of polynomial target
vector fields of total degree at most two tangent to \(\Gamma\). Nine
dimensions vanish identically on \(\Gamma\). They are obtained by putting the
three quadratic relations

\[
R_1=12a-b^2,\qquad R_2=9ac-b,\qquad R_3=3bc-4
\]

independently in any of the three vector-field components.

Modulo this vanishing subspace, the tangent space is four-dimensional.
Representatives and their induced equations for the parameter \(c\) are

\[
\begin{array}{c|c}
\text{field representative}&\dot c\\ \hline
(-27a^2/2,-27ab/4,3b/4)&c^{-1}\\
(-3ab/2,-9a,1)&1\\
(-2a,-b,c)&c\\
(-2b/9,-4/3,c^2)&c^2.
\end{array}
\]

Only \(\dot c=\lambda c\) is complete on
\(\Gamma\simeq\mathbb C^\ast\). Indeed the other three flows reach \(c=0\),
develop a pole, or acquire branching in finite complex time. It follows that
every complete degree-at-most-two ambient field tangent to \(\Gamma\) must have
the form

\[
\lambda H+W,\qquad W|_\Gamma=0.
\]

This is a necessary condition, not yet a classification of completeness:
the nine-dimensional correction \(W\) can change the ambient dynamics even
though it is invisible on \(\Gamma\). Likewise, completeness of a target field
does not automatically imply completeness of its inverse-Jacobian lift,
because a lifted sheet may still escape over \(V(Q)\).

For example, \((3bc-4)\partial_a\) is a complete triangular target field
vanishing on \(\Gamma\), but its lift is \((3F_2F_3-4)X_1\). On the explicit
orbit with \(F_2=F_3=0\), this is a nonzero constant rescaling of the incomplete
\(X_1\), so the lift is incomplete. This supplies a concrete warning that
preservation of the omitted curve is necessary but not sufficient.

There are exactly two coordinate-shear fields among the nine displayed kernel
generators:

\[
U=(3bc-4)\partial_a,\qquad V=(12a-b^2)\partial_c.
\]

Their coefficients are independent of the differentiated coordinate, so both
are locally nilpotent and complete. More strongly, every
\(\lambda U+\mu V\) is complete: \(b\) is constant, and \((a,c)\) satisfies an
affine linear constant-coefficient system.

Nevertheless, both basic lifts are incomplete. The preceding \(X_1\) orbit
handles \(U\). For \(V\), put

\[
p(r)=\left(-\frac{r+1}{r},1,-r(r-3)\right).
\]

Exact substitution gives

\[
F(p(r))=(0,1,1-r^2).
\]

Here \(12a-b^2=-1\), so after the time change \(r=\sqrt{1+s}\), this is an
integral curve of the lift \(-X_3\). It starts at
\(p(1)=(-2,1,2)\) and escapes as \(r\to0\), at finite time \(s=-1\), with
exponent \(1/2\).

Thus there is already a two-parameter family of complete polynomial target
flows preserving the omitted curve whose basic inverse-Jacobian lifts fail to
be complete. The weighted Euler symmetry remains special because its lift is
compatible with the source scaling itself.

Run .venv/bin/python archive/tooling/low_degree_symmetries.py for the enumeration.
