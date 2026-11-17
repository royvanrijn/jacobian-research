# All-degree graph obstruction with matching dual traces on two divisors

> The [first-slice theorem](FIRST_SLICE_INVERTIBILITY.md) (`HC4MYGF1/2`)
> now closes every unit-Hessian member with a base-independent Keller
> first trace, regardless of its second trace. A surviving plane has an
> explicit triangular inverse. The no-plane theorem below remains stronger
> on its matching-trace subclass; its second-trace equation is no longer
> an open counterexample gate under that first-trace hypothesis.

14 November 2026. **HC4MYGD1** is the elementary no-plane theorem below;
**HC4MYGD2** is its consequence using the reproduced structural alternative
[HCTANG2](../tang-hc4-audit/STRUCTURAL_REVIEW.md). This extends the
[triangular theorem](PROOF.md) to a class containing arbitrarily long
base-dependent alternating words. It does not cover every such word.

## Statement

Let the coefficient field have characteristic zero, let LN be nonzero,
and retain the scaled Meng--Yang formulas for A, B and Psi from
[the standard graph proof](../hc4-polynomial-graphs/PROOF.md). Put u=1+xy.
Let p0(P,Q), q0(P,Q) be any plane Keller pair: its Jacobian determinant
is a nonzero constant. No inverse map is assumed. For arbitrary polynomials
E, H, R in x,y,P,Q, set

\[
 p=p_0+xuE,\qquad q=q_0+xuH,\qquad r=R,
 \qquad f=\Psi(x,y,p,q,r).
\]

**HC4MYGD1.** The polynomial f has no constant vector two-plane along
which it is jointly affine. This conclusion requires no Jacobian or
invertibility hypothesis on the full base-dependent substitution (p,q).

**HC4MYGD2.** If det Hess(f) is a nonzero constant, grad(f) is a
polynomial automorphism. This follows from HC4MYGD1 and HCTANG2, with the
same characteristic-zero descent as in the standard graph proof. It uses
no plane degree bound or assumption of JC2. Consequently this class
cannot produce any gradient counterexample, in any polynomial degree.

The condition says exactly that both traces, on x=0 and u=0, equal the
same base-independent Keller pair. The ideals (x) and (u) are comaximal
since u-xy=1, so their intersection is (xu).

## 1. Every affine direction invisible to x is purely dual

On x=0 the potential is

\[
 f_0=Lp_0^2+Mp_0+4Ny^2p_0+Nyq_0.
\]

Suppose D_v^2 f=0 for a constant vector v with zero x component. If
its y component is nonzero, normalize v=(0,1,c,d) and write
delta=c partial_P+d partial_Q. Restricting D_v^2 f to x=0 and comparing
coefficients of y gives successively

\[
 \delta^2p_0=0,\qquad \delta^2q_0=-16\delta p_0,
 \qquad 2L(\delta p_0)^2+8Np_0+2N\delta q_0=0.
\]

Applying delta to the last equation gives -24N delta p0=0. Hence
delta p0=0 and delta q0=-4p0. If delta=0 this forces p0=0, contradicting
the Keller hypothesis. Otherwise the nonzero constant vector (c,d)
generates the unit ideal. Multiplication by the polynomially invertible
Jacobian matrix of (p0,q0) preserves its component ideal. Its image here
has component ideal (p0), which is proper since p0 is nonconstant.
This is impossible. Thus v has zero y component.

## 2. The sole possible direction determines both trace coordinates

Take a nonzero purely dual affine direction delta. Restriction to x=0
and then comparison at y=0 gives delta^2(Lp0^2+Mp0)=0. Along any
translation line a nonconstant p0 of degree m would give degree 2m,
so delta p0=0. The invertible-Jacobian ideal argument now gives
delta q0=k for a nonzero constant k.

Choose constant linear dual coordinates T,S with delta=partial_S.
Then p0 depends only on T and q0=kS+h(T). The Keller determinant is
k p0'(T), up to a nonzero constant coordinate factor, so

\[
 p_0=aT+d_0,\qquad q_0=kS+h(T),\qquad ak\ne0.
\]

There is at most one line of possible x-invisible affine directions:
two independent purely dual directions annihilating p0 would make it
constant. If no such direction exists, there is already no jointly
affine two-plane, since any two-plane meets the kernel of dx nontrivially.

## 3. Transport from x=0 disagrees with the second divisor

Suppose such a two-plane W exists. It contains the S direction and
has a second generator (1,eta,tau,sigma) in x,y,T,S coordinates.
Joint affineness gives f_SS=0 and D_v f_S=0. Thus f_S=C(x,y,T), with

\[
 C_x+\eta C_y+\tau C_T=0,\qquad C(0,y,T)=Nky.
\]

Polynomial translation along the constant flow uniquely determines

\[
 C(x,y,T)=Nk(y-\eta x).
\]

Now restrict to u=0, namely y=-1/x, in K[x,x^{-1},T,S]. The corrections
xuE and xuH vanish, so the dual traces remain the same p0,q0. The
universal elimination identity gives, with A_*=-x^3 R(x,-1/x,T,S),

\[
 f_*=L A_*^2+(M-5N/x^2)A_*+\frac{2N}{x}(kS+h(T)).
\]

Since f_SS=0, f_* has S degree at most one. If A_* had positive
S degree n, the displayed polynomial would instead have degree 2n:
the other two summands have degrees at most n and one. Hence A_*
is independent of S, and C(x,-1/x,T)=2Nk/x. The transported formula
requires C(x,-1/x,T)=Nk(-1/x-eta x). Equating them gives

\[
 \eta x^2+3=0,
\]

an impossible identity in the free variable x. This proves HC4MYGD1.
The argument excludes the whole jointly affine plane; individual affine
directions need not be absent.

## Alternating words, collisions and the remaining boundary

There is a necessary gate when only the x=0 trace is prescribed. Any
hypothetical plane must first produce the adapted coordinates of Step 2.
Write qbar=q(x,-1/x,T,S), Kstar=M-5N/x^2 and
Vstar=-x^3 R(x,-1/x,T,S)+Kstar/(2L). Integrating the transported
S derivative on the second divisor forces

\[
 \bar q=-\frac{Lx}{2N}V_*^2
         -\frac{k}{2}(1+\eta x^2)S+H_*(x,T),
 \qquad H_*\in K[x,x^{-1},T].
\]

This follows from the same elimination identity even if the second p
trace differs or the second q trace is nonlinear. In particular, if
qbar is affine in S, the square forces Vstar to be S-independent and
its S coefficient must be exactly -k(1+eta x^2)/2, independent of T.
Matching qbar=kS+h(T) fails this gate. A different trace satisfying the
gate is only a necessary survivor, not a plane or a constant-Hessian
potential. It identifies the next boundary equation for unrestricted
words without treating the matching-trace exclusion as their closure.

An exact control passes this boundary gate but fails in the interior.
Compose the three shears
Q1=Q+xP, P1=P-(3/2)xy^2 Q1, Q2=Q1+xP1, and set R=0.
Their determinant is one and their x=0 trace is the identity. On u=0,
Q2=-Q/2+xP/2, so the potential's Q derivative is -N/x, exactly the
transport of Ny for eta=0. Nevertheless A_Q is nonzero, and the full
potential satisfies f_QQ=2L A_Q^2 nonzero. The first divisor confines
any hypothetical plane's x-invisible direction to Q, so this control
has no such plane despite passing the second-trace equation.

Over K[x,y], take any finite composition of elementary shears

\[
 (P,Q)\longmapsto(P+xuF(x,y,Q),Q),\qquad
 (P,Q)\longmapsto(P,Q+xuG(x,y,P)).
\]

Every shear is an automorphism with Jacobian one and is the identity
modulo xu. Their composition is still the identity modulo xu, regardless
of length or degree. Composing on either side with any base-independent
Keller substitution retains the required two traces. Thus both the
arbitrary corrections and this automorphism subclass are covered.
This is an ideal-preservation argument for all lengths, not a finite
word enumeration. A shear coefficient divisible only by x need not
vanish on u=0 and is outside this sufficient condition.

For the marked v2 ambient collision, the pullback covector is always
-grad(p)/2. Therefore point containment and equality of the two p
covectors are still both required to transfer that particular collision.
The corrections xuE vanish on the two proof divisors, but neither marked
point lies on either divisor; those corrections do not automatically
preserve the marked jets. The invertibility conclusion applies whether
or not a marked collision transfers.

Unrestricted base-dependent words with different divisor traces,
nonlinear base changes, other graph charts, and mixed or coisotropic
reductions remain open. Unmarked unit-Hessian existence is not excluded.
No formal verification, external review or novelty priority is claimed.

## Reproduction

```sh
timeout 60s .venv/bin/python research/hessian/hc4-triangular-graphs/verify_two_divisor_extension.py
```

The checker derives the universal directional coefficients and trace
transport, checks both divisor identities with arbitrary correction
functions, and verifies a nontriangular base-dependent word and a
negative boundary control. The unbounded degree and word-length
arguments are the written proof above. [Receipt](receipt-two-divisor.json).
