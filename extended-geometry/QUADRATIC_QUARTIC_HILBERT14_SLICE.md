# Additive invariants of the normalized quadratic--quartic slice

## Result

Work over a field \(k\) of characteristic zero.  Write

\[
 A=aT^2+cT+z,\qquad
 B=pT^4+dT^3+yT^2+xT+w
\]

and put

\[
 m=ad+cp,\qquad \rho=\operatorname{Res}(A,B).
\]

The normalized tangent slice is the affine sixfold

\[
 X_{2,4}=\{m=1,\ \rho=1\}\subset\mathbb A^8,\qquad
 R=k[X_{2,4}].                                      \tag{1}
\]

Two natural commuting locally nilpotent derivations on \(R\) have an exact
finitely generated common kernel.  The calculation is the next
saturation-ledger experiment after the normalized \((2,3)\) slice.

## 1. Three gauge derivations

The low-coefficient factor translation \(B\mapsto B+tA\) gives

\[
 D_0=(0,0,0,0,0,a,c,z).                             \tag{2}
\]

The tangent-preserving Euclidean shear

\[
 B\longmapsto B+t(aT^2-2cT)A
\]

gives

\[
 D_2=(0,0,0,a^2,-ac,az-2c^2,-2cz,0).               \tag{3}
\]

Finally, take minus \(a\) times the upper binary-variable shear and correct
the quartic by \(B\mapsto B+6ptTA\).  Its infinitesimal generator is

\[
 D_{10}=(0,-2a^2,-ac,0,2ap,-3ad+6cp,-2ay+6pz,-ax).
                                                               \tag{4}
\]

All three derivations annihilate \(m\) and the resultant as ambient
polynomials.  Their displayed triangular forms prove local nilpotence.  The
brackets are

\[
 [D_{10},D_0]=[D_2,D_0]=0,\qquad
 [D_{10},D_2]=-2acD_0.                              \tag{5}
\]

The commuting pair used below is therefore \(D_{10},D_0\).  It has no
artificial rank drop along \(a=0\): there \(cp=1\),

\[
 D_{10}(y)=6,\qquad D_0(x)=c,
\]

so the two tangent vectors are independent everywhere on that divisor.

## 2. Generic quotient

Define the five common invariants

\[
 \begin{aligned}
 U={}&4az-c^2,\\
 N={}&-4a^3x+4a^2cy-3ac^2d-12acpz+4az+5c^3p-c^2,\\
 M={}&-8(-2a^4w+a^3cx+2a^3yz-a^2c^2y-3a^2cdz\\
      &\hspace{34mm}-2a^2pz^2+ac^3d+4ac^2pz-c^4p).
 \end{aligned}                                      \tag{6}
\]

On \(X_{2,4}\) they satisfy

\[
 \boxed{M^2+4UN^2=256a^4.}                          \tag{7}
\]

The exact ambient error in (7) is

\[
\begin{aligned}
M^2+4UN^2-256a^4\rho
={}&4U^2(m-1)
 (8a^3x-8a^2cy-4a^2dz+7ac^2d\\
 &\quad+20acpz-4az-9c^3p+c^2).
\end{aligned}
\]

This is the complete common quotient after \(a\) is inverted.  Indeed put

\[
 s=-{c\over2a^2},\qquad
 Y={4a^2y-6acd+3c^2p\over4a^2},\qquad t={Y\over a},
                                                               \tag{8}
\]

and

\[
\begin{aligned}
 X={}&{4a^3x-4a^2cy+3ac^2d+12acpz-5c^3p\over4a^3},\\
 W={}&-{ -8a^4w+4a^3cx-2a^2c^2y+ac^3d
                 +6ac^2pz-2c^4p\over8a^4},\\
 L={}&W-{UY\over4a^2}.
\end{aligned}                                      \tag{9}
\]

Then

\[
 D_{10}(s)=1,\quad D_0(s)=0,\qquad
 D_{10}(t)=0,\quad D_0(t)=1,
\]

while

\[
 N=U-4a^3X,\qquad M=16a^4L+U^2p.                   \tag{10}
\]

Equations (8)--(10), together with \(m=1\), reconstruct
\((c,z,d,y,x,w)\) from
\((a,p,U,N,M,s,t)\) after \(a\) is inverted.  Consequently, if

\[
 B_0=
 k[a,p,U,N,M]\big/(M^2+4UN^2-256a^4),
\]

then

\[
 R_a=(B_0)_a[s,t],\qquad
 (\ker D_{10}\cap\ker D_0)_a=(B_0)_a.               \tag{11}
\]

## 3. The boundary saturation ladder

The first quotient model is not saturated at \(a=0\).  On that boundary

\[
 p=c^{-1},\qquad U=-c^2,\qquad N=4c^2,\qquad M=8c^3,
                                                               \tag{12}
\]

so \(B_0/(a)\) remembers only \(c\), whereas the boundary quotient by the
two actions has coordinates \(c,z,d\).

Three regular common invariants repair the defect:

\[
\begin{aligned}
C={}&-4(a^2x-acy+2c^2d+3cpz-5z),\\
Q={}&4p^2z+2d-ad^2,\\
S={}&-8(-2a^3p^3w+a^2cp^3x+2a^2p^3yz-ac^2p^3y\\
 &\quad-3acdp^3z-2ap^4z^2+2c^3dp^3+c^2dp^2
       +4c^2p^4z+cdp+d).
\end{aligned}                                      \tag{13}
\]

The letter \(Q\) is used here for the invariant called `R` in the checker,
to avoid confusing it with the coordinate ring.  In \(R\),

\[
 aC=N+4U,\qquad
 aQ=p^2U+1,\qquad
 aS=p^3M-8.                                        \tag{14}
\]

Thus these are precisely the first three regularized boundary fractions.
They are not a bounded guess: each equality in (14) has an explicit ambient
error divisible by \(m-1\).

The possible Hensel ladder terminates at \(S\).  Substituting (14) into (7)
and cancelling one factor of \(a\) gives

\[
\begin{aligned}
0=\Phi={}&
4C^2Qa^2p^4-4C^2ap^4-32CQ^2a^2p^2+64CQap^2-32Cp^2\\
&+64Q^3a^2-192Q^2a+192Q+S^2a+16S-256a^3p^6.
\end{aligned}                                      \tag{15}
\]

Modulo \(a\), this is the linear equation

\[
 S=2Cp^2-12Q.                                      \tag{16}
\]

Hence division of any further Hensel remainder by \(a\) produces a
polynomial in the already adjoined generators; there is no infinite
coefficient ladder.

## 4. Exact global kernel

Set

\[
 B=k[a,p,U,N,M,C,Q,S]\subset R.                     \tag{17}
\]

By (14), \(B_a=(B_0)_a\).  It remains to exclude invariants supported only
at \(a=0\).

The boundary is explicit.  Since \(cp=1\) and the resultant of the linear
form \(cT+z\) with \(B\) is linear in \(w\),

\[
 R/(a)=k[c^{\pm1},z,d,y,x].
\]

The two derivations act freely in the \(y,x\) directions, and the invariant
boundary coordinates from (13) are

\[
 p=c^{-1},\qquad
 C=8(z-c^2d),\qquad
 Q=2d+{4z\over c^2}.                                \tag{18}
\]

Their Jacobian with respect to \((c,z,d)\) is

\[
 -{48\over c^2}\ne0.
\]

More explicitly,

\[
 d={2Q-p^2C\over12},\qquad
 z={C\over24}+{Q\over6p^2}.                         \tag{19}
\]

Relations (12), (14), and (16) show that \(B/(a)\) is generated by
\(p^{\pm1},C,Q\).  Equations (18)--(19) prove that

\[
 B/(a)\longrightarrow R/(a)
\]

is injective.

Now let \(f\in R\cap B_a\), and write \(f=b/a^n\) with \(n\) minimal.
Reducing \(a^nf=b\) modulo \(a\), boundary injectivity forces
\(b\in aB\), contradicting minimality unless \(n=0\).  Combining this with
(11) proves the exact result

\[
\boxed{\ker D_{10}\cap\ker D_0
 =k[a,p,U,N,M,C,Q,S].}                              \tag{20}
\]

In particular the common invariant ring of this natural
\(\mathbb G_a^2\)-action is finitely generated.  Formula (20) is an equality
of the displayed subring with the kernel; it does not claim that the eight
generators are minimal, nor does it classify all LNDs on \(X_{2,4}\).

## 5. The noncommuting triple intersection

The third action does not commute with \(D_{10}\), but (5) implies that it
preserves their common kernel: if \(D_{10}(f)=D_0(f)=0\), then

\[
 D_{10}(D_2f)=D_2(D_{10}f)-2acD_0(f)=0,
 \qquad D_0(D_2f)=0.
\]

On the generators in (17), the induced LND is

\[
\begin{aligned}
D_2(a)&=D_2(U)=D_2(N)=D_2(M)=D_2(C)=0,\\
D_2(p)&=a^2,\qquad
D_2(Q)=2apU,\qquad
D_2(S)=3ap^2M.                                     \tag{21}
\end{aligned}
\]

After inverting \(a\), \(p/a^2\) is a slice.  Relations (14) recover \(Q,S\)
from \(p\), so the localized kernel is generated by
\(a,U,N,M,C\).  Put

\[
 A_3=k[a,U,N,M,C]
 \subset R.
\]

The two relations are

\[
 aC=N+4U,\qquad M^2+4UN^2=256a^4.                  \tag{22}
\]

There is again no hidden boundary saturation.  Modulo \(a\),

\[
 A_3/(a)
 =k[U,M,C]/(M^2+64U^3),
\]

and its map to the boundary quotient is the cusp parametrization

\[
 U=-c^2,\qquad M=8c^3,\qquad C=8(z-c^2d).           \tag{23}
\]

The map \(k[U,M]/(M^2+64U^3)\to k[c]\) in (23) is injective: every class has
a unique form \(f(U)+Mg(U)\), whose image has disjoint even and odd powers
of \(c\).  The coordinate \(C\) is algebraically independent from \(c\).
Thus \(A_3/(a)\to R/(a)\) is injective, and the same minimal-pole argument
used for (20) proves

\[
\boxed{\ker D_{10}\cap\ker D_0\cap\ker D_2
 =k[a,U,N,M,C].}                                   \tag{24}
\]

Equivalently, eliminating \(N=aC-4U\) gives the hypersurface presentation

\[
 k[a,U,M,C]\big/
 \left(M^2+4U(aC-4U)^2-256a^4\right).               \tag{25}
\]

The noncommuting triple therefore also has a finitely generated invariant
ring.  Its cusp boundary is geometrically nonnormal as a quotient model,
but it does not produce a Hilbert--14 escape ladder.

## 6. Hilbert--14 assessment

The \((2,4)\) slice is a stronger positive control than the \((2,3)\)
example:

- the primitive commuting pair remains rank two on the leading boundary;
- the generic quotient is again a small trinomial;
- the boundary initially loses two quotient coordinates;
- three explicit saturation classes restore them;
- the last defining relation is boundary-linear, forcing termination.
- the additional noncommuting Euclidean action cuts the result down to the
  finitely generated cusp-boundary hypersurface (25).

Thus neither of the first two unequal-degree normalized factorization slices
produces a non-finitely-generated invariant intersection.  The next useful
fork is no longer simply “increase one degree.”  One should move to a slice
with at least two independent leading boundary divisors, where a multigraded
conductor quotient can support a genuine escape ladder.

All polynomial identities, local nilpotence, brackets, generic slices,
boundary formulas, and the terminating relation are checked by
[`verify_quadratic_quartic_additive_invariants.py`](../scripts/verify_quadratic_quartic_additive_invariants.py).
