# Smooth value families and the motivic homotopy barrier

Theorems HC4MYJZ0/1/2, over C. The complete unstable A1-homotopy type
over a dense open value line supplies no additional exclusion of the
remaining normalized models. This strengthens the earlier passing Euler
and Hodge controls to the relative motivic homotopy category itself.
It does not construct a Keller pair, an action on the native family,
a source equivalence, or a unit-Hessian member.

Write H_*(S) for the pointed Morel--Voevodsky A1-homotopy category over
a smooth affine complex base S. Every equivalence below is in that
category, allowing zigzags. Define

    M_N = wedge_(i=1)^N (P1 smash Gm)
        = wedge_(i=1)^N (S_s^1 smash Gm smash Gm).

Here S_s^1 is the simplicial circle; N=0 means the point object. A
constant type over S means the base change of M_N, not a trivialization
as a family of schemes. Base points and coordinate trivializations are
chosen; no canonical marked equivalence is asserted.

## 1. The two imported homotopy inputs

We use [Asok--Dubouloz--Ostvaer, 2112.08241v2](https://arxiv.org/html/2112.08241v2),
read as a full primary paper, not only the earlier retained abstract:

- Theorem 2.3.1, homotopy purity: a smooth closed pair has cofiber the
  Thom space of its normal bundle. We use the usual ambient/complement
  notation Y/(Y-X) for X closed in Y; the displayed HTML formula reverses
  those variable names, while the surrounding statement and proof identify
  the standard theorem.
- Theorems 3.2.1 and 3.2.2: a smooth morphism to affine space with a
  section, an A1-equivalence away from a split degeneration locus, has
  homotopy type the indicated P1 suspension of its pointed special
  fibre, or the wedge of those suspensions for several split points.

All applications have trivial affine-space bundles off the stated locus,
explicit sections, and smooth closed fibres. Vector-bundle projections
and Zariski locally trivial affine-line bundles are A1-equivalences by
local triviality and Nisnevich descent. The homotopy theorems are proved
external inputs, not independently verified by the symbolic checker.

## 2. Every smooth cotangent value family has type M_N (HC4MYJZ0)

Let (F,G):A2_b -> A2 be a plane Keller map with exactly N points
Z in its zero fibre, let d(b) be arbitrary, and set

    H(b,P,Q)=d(b)+P F(b)+Q G(b).

The zero fibre is reduced by the Keller condition. Choose a dense open
value base S contained in A1_s minus {d(b_i): b_i in Z}. On H=s there
are no points over Z. Projection to U=A2 minus Z is the entire family

    X_S -> U times S.

It is a Zariski locally trivial affine-line torsor: over F!=0 solve
P=(s-d-QG)/F, and over G!=0 solve Q=(s-d-PF)/G. This is an actual
projection and these are actual local trivializations, independent of
properness of the plane Keller map. Hence X_S is A1-equivalent over S
to the constant family U times S.

We now identify that type for arbitrary finite N, without assuming a
Keller map exists for any particular N. Move the N distinct complex
points by an actual polynomial plane automorphism to

    Z_N={(X,z)=(0,i): 1<=i<=N}, C_N(z)=product_(i=1)^N(z-i).

The finite-point interpolation and its inverse are retained in JV2.
For N>=1 the smooth affine threefold

    W_N={X P+C_N(z) Q=1}

projects onto U_N=A2 minus Z_N as an affine-line torsor. Its charts
are the same two solved charts as above, so W_N and U_N are
A1-equivalent. Projection W_N -> A1_z is smooth: if C_N!=0 the
Q derivative is a unit, and if C_N=0 then XP=1. It has the section
X=P=1,Q=0. Away from the N simple roots it is a trivial A2-bundle;
at every root its fibre is Gm_X times A1_Q.

Theorem 3.2.2 therefore gives W_N equivalent to M_N. Its hypotheses
also hold after base change to S. For N=0, U=A2 is contractible and
the same conclusion with M_0 follows directly. Thus every smooth
cotangent value family above has constant relative type M_N, for
every N, every d and every Keller degree. We have neither enumerated
N nor equated one zero-fibre cardinality with generic Keller degree.

## 3. Hyperbolic suspensions of punctured lines (HC4MYJZ1)

Let B_S=A1_t times S, and let r_1,...,r_N be pairwise disjoint
sections. Suppose a chosen section t_0 avoids all of them. Put

    C_S=B_S minus union_i{t=r_i}, pointed by t_0.

Point A1_S by t_0. It is A1-contractible over S. Homotopy purity for
the N disjoint sections, with their coordinate normal trivializations,
therefore gives

    Sigma_s(C_S,t_0) equivalent to wedge_N P1_S.         (1)

Indeed the cofiber of C_S -> A1_S is both this simplicial suspension
and the wedge of the N trivial normal Thom spaces. Smashing (1) with
Gm gives

    P1 smash (C_S,t_0) equivalent to (M_N)_S.             (2)

Now let p(z,w) be a polynomial over O(S) whose whole zero scheme
is smooth over S and is isomorphic over S to C_S. Form

    Y_S={x y=p(z,w)}.

Choose the point (z,w) corresponding to t_0. Projection Y_S -> A1_x
over S has section y=0,(z,w)=(z_0,w_0). It is smooth: for x!=0
solve y, and at x=0 smoothness of p=0 makes one relative derivative
of p nonzero. Over x!=0 it is the trivial A2-bundle in z,w. Its
zero fibre is A1_y times C_S. Theorem 3.2.1 and (2) yield

    Y_S equivalent to P1 smash C_S equivalent to (M_N)_S. (3)

The claim is about the whole zero scheme and actual relative smoothness.
An open parametrization of just one component would not justify it.
It makes no assertion that Y_S is isomorphic to W_N or has its torsor
projection, and no actual Ga action is transported by (3).

## 4. Every normalized test has exactly the required relative type (HC4MYJZ2)

Use a value coordinate s. Shrink to the following explicit dense opens;
removing s=0 is convenient for a common section, not a new critical value:

| Native model | Base S | N |
|---|---|---:|
| coordinate endpoint | any dense open | 0 |
| nondegenerate quadratic endpoint | s!=0 | 1 |
| h_D=2xy+(z^2-1)w+D(z) | s!=0,D(1),D(-1) | 2 |
| K2=w(2xy+z^2-1)+z | s!=0,1,-1 | 2 |
| S_rho=2PV+B_rho(u,v) | s!=0,1,-1,rho | 3 |

Repeated excluded values are removed once. Every row has relative
A1-homotopy type the constant (M_N)_S. The argument is over the
whole indicated open value base, not merely point by point or after
an algebraic closure of C(s).

For h_D use the suspension x*(2y)=s-B(z,w), B=(z^2-1)w+D(z).
Its whole zero curve is

    z in A1 minus {1,-1}, w=(s-D(z))/(z^2-1).

At z=+-1 the excluded value condition prevents additional points.
The inverse is z. Choose t_0=z=0,w=D(0)-s, a regular section over S.
Section 3 gives N=2 for arbitrary D. In particular h0 and h1 pass.

For K2 retain the exact C[s]-algebra presentation from JS1:

    2xy+(s+uv)^2+v-1=0,
    s=K2,u=w,v=1-2xy-z^2,z=s+uv.

Both inverse identities hold globally. This is an isomorphism of the
total hypersurface family over the value line, not an asserted four-
coordinate source automorphism using s as an independent source
coordinate. Its whole exceptional plane curve is

    z in A1 minus {1,-1}, v=1-z^2,u=(z-s)/(1-z^2).

If v=0 then s^2=1, excluded on S; there are no extra points. The
inverse is z=s+uv and the section z=0 gives u=-s,v=1. Section 3
again gives N=2.

For S_rho use x=2P,y=V and p=s-B_rho. JQ0's whole curve formula is

    t in A1 minus {-1,1,s},
    u=(s-t)/(t^2-1), v=(t-rho)(t^2-1)/(s-t),
    t=uv+rho.

The conditions s!=+-1,rho remove every exceptional component.
These are inverse isomorphisms over S. The three omitted sections
-1,1,s are pairwise disjoint there. As s!=0, t_0=0 avoids them,
giving the regular section u=-s,v=rho/s. All hypotheses of section 3
hold for every fixed complex rho, including rho=+-1, and give N=3.
The moving fourth puncture on the projective completion therefore
does not obstruct the relative P1-suspended homotopy type.

The coordinate row has fibre A3 and is trivial over the value line.
For the quadratic row, split its nondegenerate quadratic into two
hyperbolic pairs. On s!=0, rescaling the two coordinates of one
isotropic plane identifies its family with {X P+z Q=1} times S;
this is W_1 after translating z. Hence it has type M_1. These are
actual linear/scalar family identifications, not arbitrary nonlinear
Hessian invariance.

Any hypothetical nonautomorphic unit-Hessian member with N finite
Morse critical points has the cotangent description of section 2,
and exactly this zero-fibre count by HG0. Its critical values agree
with the native ones under the actual source automorphism. Thus the
same open S can be used, and its relative type is the same constant
(M_N)_S. For both endpoint alternatives the respective counts are N=0,1.

Consequently no invariant that factors only through this unmarked
relative unstable A1-homotopy type over S can distinguish any remaining
native family from the required cotangent family. This includes its
stable motivic consequences. The earlier whole-value Hodge controls
and complete local volume-period controls remain consistent with this
stronger passing test.

This does not trivialize the torsor class, identify actual schemes or
actions, supply globally Keller coefficient generators, or identify
source orbits. Gluing across critical values, marked volume/quotient
data, geometry at infinity retaining more than this homotopy object,
and actual graph embedding recovery remain separate. The exact
contracted universal frontier, E2, the individual full-source orbits,
both old endpoints and unrestricted HC4/JC2 remain open. There is no
new generic-fibre homotopy calculation to deepen within these tests.

## Reproducibility and assurance

[The context](sources/motivic-family-barrier-context.json) pins the full
versioned primary paper, the exact proved homotopy imports, their checked
section/smoothness/splitting hypotheses, and retained structural proofs.
[The checker](verify_motivic_family_barrier.py) checks the actual torsor
charts, family inverse identities, full punctured-curve formulas,
sections and exceptional-value controls.
[Receipt](receipt-motivic-family-barrier.json).
The symbolic replay does not implement motivic homotopy theory or
independently prove the imported theorems. Independent replay, formal
verification and external review are not claimed.

