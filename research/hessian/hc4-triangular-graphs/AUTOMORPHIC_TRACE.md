# All-degree obstruction for every base-preserving dual automorphism

> The [rational-curve extension](CURVE_TRACE.md), `HC4MYGC1/2/3`, now
> permits every polynomial base automorphism as well. It also closes
> the standard graph eliminating p, with the stated relative affine
> dual extension, and reduces its nonlinear-dual case to a globally
> affine-dual Keller endpoint. The Hamiltonian and trace lemmas below
> remain part of that proof.

14 November 2026. **HC4MYGV2** closes every polynomial graph after an
arbitrary dual automorphism over K[x,y], with unrestricted degrees and
word lengths. It extends the [first-slice result](FIRST_SLICE_INVERTIBILITY.md)
to genuinely y-dependent first traces. The proof uses two separately
identified external inputs: the reproduced [Tang alternative](../tang-hc4-audit/STRUCTURAL_REVIEW.md)
and the classical Rentschler criterion in [Furter, Lemma 4.2](https://www.math.u-bordeaux.fr/~jpfurter/planePolynomialAutomorphismsOfFixedMultidegree.pdf).
The latter is retained with its [source receipt](sources/automorphic-trace-import.json).

## 1. The trace criterion and the resulting graph theorem

Write h(y,P,Q)=f(0,y,P,Q). First work over C. Assume

1. h is a coordinate of C(y)[P,Q]; and
2. h(0,P,Q)=Ls(P,Q)^2+Ms(P,Q), where L is nonzero and s is nonconstant.

**HC4MYGV1.** If f has a jointly affine constant vector two-plane and
its Hessian determinant is a nonzero constant, then its gradient is
polynomially invertible. The proof below exhausts planes both contained
in and transverse to x=0. It uses Rentschler, but no Tang classification
or assumption of JC2.

**HC4MYGV2.** Under the same trace assumptions, every f with nonzero
constant Hessian determinant has a polynomially invertible gradient.
Otherwise HCTANG2 supplies the plane contradicted by HC4MYGV1.
This conclusion descends to any characteristic-zero field by the
finite-coefficient-field and unique-formal-inverse argument used in
the standard graph theorem.

In particular take the scaled Meng--Yang Psi=LA^2+MA+NB, with LN nonzero.
Allow arbitrary polynomial p(x,y,P,Q), q(x,y,P,Q), R(x,y,P,Q), requiring
only that their x=0 dual trace (p0(y,P,Q),q0(y,P,Q)) be an automorphism
of K[y][P,Q]. Then

\[
 h=Lp_0^2+(M+4Ny^2)p_0+Nyq_0.
\]

Over K(y), the pair (p0,h/(Ny)) is an automorphism, by a triangular
shear of the trace pair. Its y=0 specialization has
h(0)=Ls^2+Ms with s=p0(0,P,Q) a coordinate, hence nonconstant.
Both trace assumptions follow. No condition on the second divisor
or full off-slice dual Jacobian is needed.

Every dual automorphism over K[x,y] specializes to such an automorphic
trace over K[y]. Consequently **every** base-dependent polynomial dual
automorphism, including arbitrary alternating words, is covered, with
arbitrary polynomial r=R. This proves gradient invertibility for every
unit-Hessian member of that entire class. It does not assert that such
members or affine fibre planes are absent.

## 2. The Hamiltonian pencil lemma

Use {a,b}=a_U b_V-a_V b_U and D_a={a,-} on C[U,V].

**HC4MYGV0.** If b+t a is a coordinate of C(t)[U,V], then D_a is
locally nilpotent. If a has a polynomial Jacobian mate A with
{A,a} a nonzero constant, then (A,a) is a polynomial automorphism.

For a coordinate c, choose a mate d with C(t)[U,V]=C(t)[c,d].
Its Hamiltonian derivation is a nonzero C(t)-scalar multiple of
partial_d, so it is locally nilpotent. Apply this to c=b+t a.
There are finite exponents n_U,n_V killing U,V under (D_b+t D_a).
These are polynomial identities in t, even though the mate may have
rational t coefficients. The coefficient of t^n in

\[
 (D_b+tD_a)^n(g)
\]

is exactly D_a^n(g). Thus D_a kills U and V after sufficiently many
iterations, proving local nilpotence on the whole polynomial ring.
This argument does not suppose the two derivations commute, and it
does not infer nilpotence from finitely many tested parameter values.

The external Rentschler criterion says that such an a is a polynomial
H(z) in a coordinate z; see Furter Lemma 4.2. The Jacobian mate gives

\[
 \{A,a\}=H'(z)\{A,z\}\in C^\times.
\]

Both polynomial factors must be nonzero constants, so H is linear and
a itself is a coordinate. In coordinates (a,w), the constant Jacobian
then makes A a nonzero constant multiple of w plus a polynomial in a.
This is a triangular automorphism. This proves HC4MYGV0 without JC2.

Only the classical implication D_a locally nilpotent => a=H(z) is
imported here. We do not use a finite nilpotence bound or Furter's
separate closedness theorem.

## 3. A purely dual affine direction has a scalar trace coefficient

If a nonzero constant pure-dual direction belongs to the affine plane,
choose constant linear dual coordinates T,S making it partial_S.
Then f_SS=0 implies

\[
 h(y,T,S)=C(y,T)S+B(y,T).
\]

Over C(y), h is a coordinate. Hence D_h is locally nilpotent. On the
subring C(y)[T], D_h acts as -C(y,T) partial_T. Such a derivation is
locally nilpotent only if C has T degree zero (or is zero). Indeed, a
positive degree d gives successive nonzero leading terms: degree
1+n(d-1), with nonzero leading coefficient in characteristic zero;
degree one gives the same nonzero linear factor at every iteration.
It follows that C=C(y) is independent of T.

At y=0, affine S dependence of Ls^2+Ms forces s_S=0 by polynomial
degree along S. Thus C(0)=0. If C vanishes identically, coordinate h
depends only on T and must be linear in T over C(y). Polynomiality
would then give T degree at most one at y=0, whereas Ls(T)^2+Ms(T)
has degree at least two. Consequently C is not identically zero.

This scalar-coefficient argument uses coordinate, not merely Keller,
information. It therefore does not apply to the unrestricted Keller
trace encoding in HC4MYGF3.

## 4. No entire affine fibre plane lies in x=0

Suppose W is contained in ker(dx). It has a nonzero pure-dual vector
because its intersection with ker(dy) has positive dimension. Use
Section 3 to write h=C(y)S+B(y,T), with C(0)=0 and C not identically
zero. If W were the pure-dual plane, the slice at y=0 would be affine
in P,Q, contradicting the nonconstant square Ls^2+Ms. Thus W has
another direction that moves y.

The mixed derivative between that direction and partial_S is a
nonzero constant multiple of C'(y), so joint affineness forces C'=0.
Together with C(0)=0 this gives C=0, contradicting Section 3. This
excludes the whole tangent-plane case; transversality is not assumed.

## 5. A transverse plane with a purely dual intersection is triangular

Now W is transverse to x=0. Let its intersection with ker(dx) be
purely dual, giving partial_S as above. Normalize another generator
to (1,eta,tau,sigma) in x,y,T,S coordinates. Whole-plane transport
of f_S from x=0 uniquely gives

\[
 f_S=C(y-\eta x).
\]

In constant linear coordinates X=x, Y=y-eta x, T'=T-tau x,
S'=S-sigma x, the plane is span(partial_X,partial_S') and

\[
 f=C(Y)S'+XU(Y,T')+V(Y,T').
\]

The Hessian determinant is (C'(Y)U_T')^2. A nonzero constant determinant
forces C' and U_T' to be nonzero constants. Thus the base pair
(U,C) is triangular and invertible. The explicit inverse from
HC4MYGF1 applies with the nonzero constant C' replacing Nk, and with
the harmless constant term of C retained. Here C(0)=0 already removes
that term. This proves invertibility in this case.

## 6. A transverse plane whose intersection moves y has an LND endpoint

Finally normalize the x-invisible plane direction to (0,1,c,d).
Choose the other generator with x component one, then subtract a
multiple of this direction to make its y component zero. Constant
linear coordinates

\[
 X=x,\quad Y=y,\quad U=P-cy-ex,\quad V=Q-dy-fx
\]

make W=span(partial_X,partial_Y). Joint affineness gives exactly

\[
 f=b(U,V)+X A(U,V)+Y a(U,V).
\]

On X=0 the dual change is only a translation depending on y, so the
trace-coordinate property persists: b+y a is a coordinate over C(y).
HC4MYGV0 makes D_a locally nilpotent. Unit Hessian gives

\[
 \det\operatorname{Hess}f=\{A,a\}^2\in C^\times,
\]

so (A,a) is an automorphism by that same lemma. Recover U,V from the
two fibre gradient outputs A,a. The other two outputs then recover
X,Y through the polynomial inverse of D(A,a)^T after subtracting
grad(b). This proves full gradient invertibility and completes the
plane exhaustion in HC4MYGV1.

## Scope, remaining mechanisms and verification

The former y-dependent dual-automorphism continuation is now closed
in all degrees and word lengths. The weaker automorphic-first-trace
criterion also allows arbitrary off-slice corrections whose full
dual map is not invertible. No preservation of a Hessian determinant
under nonlinear reparameterization was assumed: the hypothesis is
always on the actual pulled-back potential. The marked v2 collision
still transfers only when its two pulled-back p covectors agree.

This does not close arbitrary nonlinear changes of the base x,y,
other eliminated-coordinate charts, or genuinely mixed/coisotropic
reductions. Nor does it close arbitrary polynomial substitutions with
only a Keller first trace: HC4MYGF3 identifies that larger class with
the plane conjecture. General HC4 and JC2 remain open.

The [checker](verify_automorphic_trace.py) checks the normal forms,
block determinants, coefficient transport, noncommuting Hamiltonian
pencil identities and exact representative inverses. It includes a
genuinely y-dependent nontriangular trace automorphism, with explicit
inverse, and a degenerating coordinate pencil whose leading member
is a square. The latter prevents replacing the Jacobian-mate condition
by coordinatehood of the generic pencil alone. The unit-Hessian control
f=(x+1)Q+yP has a tangent fibre plane and a generic-coordinate trace,
but its special trace Q is linear; it checks why the nonlinear square
at y=0 is required when excluding tangent planes. The all-degree and
all-plane arguments are written proofs. [Actual receipt](receipt-automorphic-trace.json).

```sh
timeout 60s .venv/bin/python research/hessian/hc4-triangular-graphs/verify_automorphic_trace.py
```

The classical Rentschler import and the separate Tang reproduction
are external theorem inputs. Their scope is explicit; no formal
verification, independent external review or novelty priority is claimed.
