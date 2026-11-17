# Leading valuations and critical curves close further nonmonic outputs

15 November 2026. Two independent mechanisms reach beyond the nilpotent
and finite-fibre gates. Leading-degree drops constrain every source
factor degree, in arbitrary z degree. Separately, a rational root would
force the normalized critical curve to be polynomially parametrizable.
This closes large genuinely nonabelian cubic families with no special
Keller-fibre hypothesis and with unrestricted parameter degrees.

Throughout, K=(f,g) is complex plane Keller, B=C[s,t], L=C(s,t),
and P(T,z)=sum_(i=0)^m b_i(T)z^i, m>=2, has a=b_m!=0.
For a closed stable embedding j:A^(r+2)->A^(r+3), r>=1, carrying
K through (T,xy-P(T,z)), the retained
[global filtration](LEADING_FIBRE_UNIT_GATE.md#3-all-degree-nonmonic-output-equations-and-their-remaining-escape-hc4myhz2)
works for every P here. Whenever P(f,Z)+g is irreducible over L,
root transitivity kills the invariants of Z^m/Z(1,...,1), and the
retained component-pencil and global primitive argument forces K
invertible. This uses the transitivity proof in HC4MYHY1, not an
assertion of full splitting-group preservation.

## 1. Leading-degree drops constrain all factor degrees (HC4MYIC0)

Choose distinct complex levels c1,...,ck with k>=1. Write
q_j=ord_(T=c_j)(a)>0. Suppose

\[
        ord_{c_j}(b_i)\ \geq\ i q_j/m,
        \qquad 1\leq i<m,\quad 1\leq j\leq k,            \tag{1}
\]

where a zero coefficient has infinite order. These are integer-order
inequalities, equivalently at least ceil(i*q_j/m). Put
delta=gcd(m,q1,...,qk). Every nonconstant factor of P(f,Z)+g
over L has degree divisible by m/delta. In particular

\[
            delta=1\quad\Longrightarrow\quad
                P(f,Z)+g\text{ irreducible over }L.        \tag{2}
\]

Consequently every carried output embedding above has invertible K
under (1)-(2), with all degrees unrestricted and no irreducibility
assumption on any f=c fibre. More generally the number of irreducible
source factors is at most delta.

To prove this, take any source prime h dividing f-c_j. The
[square-free pullback argument](ABELIAN_CLASSGROUP_OUTPUT_OBSTRUCTION.md#1-no-abelian-intermediate-extension-hc4myhw0)
gives v_h(f-c_j)=1 and dominant image of h=0 in T=c_j. The target
polynomial V+b0(T) is irreducible and distinct from T-c_j, so their
pullbacks are coprime by that same argument. Thus

\[
 v_h(b0(f)+g)=0,\qquad v_h(b_i(f))=ord_{c_j}(b_i).       \tag{3}
\]

Extend v_h to an algebraic closure of L, normalized to have value
group Z on L. Every root alpha of P(f,Z)+g has value exactly
-q_j/m. If its value were larger, the constant term in (3) would
be the unique term of minimum value. If smaller, the leading term
would be the unique minimum: for i<m its value difference is at
least (m-i)*(-q_j/m-v_h(alpha))>0. Both contradict a sum equal
to zero. This elementary argument proves the required single-slope
conclusion directly; a Newton-polygon irreducibility theorem is not
silently imported.

A monic factor of degree e over L has constant term, up to sign,
the product of e roots, of value -e*q_j/m. This value is integral
on L, hence m/gcd(m,q_j) divides e. Taking all j gives the
least common multiple m/delta, proving the factor-degree statement
and (2). No source factor coefficients were assumed integral at the
leading-coefficient divisor.

For example, all

\[
             P=a(T)R(T,z)+b(T),\qquad R\text{ monic in }z,
                 \deg_z R=m,
\]

close if gcd(m, all root multiplicities of a)=1. The coefficients
of R and b are otherwise arbitrary. One simple zero of a suffices
for every m, even with nonnilpotent root group. Multiple levels can
cooperate: for m=6 and a=T^2*(T-1)^3, the two places force factor
degrees divisible by three and by two, hence by six. Neither place
alone has its multiplicity coprime to six.

## 2. Rational labels force affine-line critical curves (HC4MYIC1)

Suppose h in L is a rational root of P(f,Z)+g. Set u=a(f)*h and

\[
 Q(T,U,V)=a(T)^{m-1}\big(P(T,U/a(T))+V\big).             \tag{4}
\]

This is a polynomial, monic of degree m in U. Its derivative Q_U
is independent of V. Since Q(f,u,g)=0, u is integral over C[f,g]
and hence over the normal source ring B; therefore u belongs to B.
Differentiating along f gives the exact polynomial identity

\[
 J(f,u)\,Q_U(f,u)=-a(f)^{m-1}J(f,g).                   \tag{5}
\]

The right side is nonzero. Thus Phi=(f,u):A2->A2_(T,U) is a
dominant generically finite polynomial map, and the image of Phi
over a(T)!=0 omits every point of Q_U=0. Both polynomial factors
on the left of (5) have only vertical source divisors over zeros
of a. If J(f,u) vanishes on a smooth irreducible component of such
an f-fibre, u is constant along that component: the nowhere-zero
Hamiltonian tangent derivation J(f,.) kills its restriction. This
does not make g constant or declare K critical there.

Every irreducible component of Q_U=0 must be polynomially
parametrizable by A1. To justify this global conclusion, take a
nonzero target polynomial d for which Phi is finite over d!=0;
it exists by clearing the denominators of the algebraic equations
of s,t over C(f,u). Hence the nonproperness set S_Phi is contained
in the plane curve d=0. Every omitted target point belongs to
S_Phi: properness on a neighbourhood would make its image closed
there and, by dominance, the whole neighbourhood.

Now apply the explicit primary input
[Jelonek--Lason, Theorem 3.2](https://arxiv.org/html/1411.5011v2):
S_Phi is covered by nonconstant polynomial images of A1. A component
C of Q_U=0 cannot be vertical, because Q_U has constant nonzero
leading U coefficient m. Choose a point on C outside a=0 and the
other components of d=0. Its polynomial curve in S_Phi must lie
in C, since an irreducible image in d=0 through that point lies in
its unique component. This supplies the claimed parametrization.
Only the qualitative parametric-curve conclusion of that theorem
is imported, not its degree bound; its hypotheses hold for the
constructed Phi, which need not itself be Keller.

Equivalently each such critical component has normalization A1.
Indeed a polynomial parametrization lifts to the affine normalization
by integrality and normality of C[t]. It extends to a map of smooth
projective curves from P1. Riemann--Hurwitz forces the target genus
to be zero, and surjectivity permits at most one omitted point at
infinity; an affine curve has at least one. Thus its normalization
is A1. Conversely such a normalization provides a parametrization.
These curve consequences use the retained characteristic-zero
Riemann--Hurwitz input from HC4MYHZ1.

There is a quantitative unit consequence as well. If Q_U has ell
distinct roots beta_i in C(T), each beta_i belongs to C[T] by
monicity and normality of C[T]. Each polynomial u-beta_i(f) divides
Q_U(f,u), so (5) makes it a unit of B[1/a(f)]. Their classes are
independent even modulo C(f)^*: a multiplicative relation among
them would be a rational-function identity in the transcendental
u over C(f), and its orders at u=beta_i(f) force every exponent
to vanish. The exact unit lattice of HC4MYHZ0 therefore gives

\[
          \sum_{a(c)=0}(r_c-1)\ \geq\ ell,               \tag{6}
\]

where r_c is the number of full irreducible f=c components. This
counts rational critical branches, not all branches after a parameter
field extension. In particular two distinct rational critical branches
cannot be supported by a rank-one localized-unit defect.

There is also an independent no-pole criterion. If a is square-free
and a divides b_(m-1), then P(f,Z)+g has **no rational root** in L.
In lowest terms h=A/B, the UFD rational-root argument gives B|a(f).
Any pole therefore has integer order k>=1 at a source component
over a simple zero of a. At this valuation the leading term has
value 1-m*k. The next term has value at least 1-(m-1)*k, and
every lower positive-degree term has value at least -i*k for i<=m-2;
the constant term has value zero by (3). The leading term is the
unique minimum, a contradiction. Thus h is polynomial. But then

\[
                  J(f,g)=-P_z(f,h)J(f,h)
\]

is impossible: f,h are algebraically independent because f,g are,
so P_z(f,h) is a nonconstant polynomial for m>=2 and cannot divide
a nonzero constant. Constant nonzero a is included, with no poles
to consider. Rational-root exclusion alone is not claimed to imply
irreducibility in degrees four and higher.

## 3. Cubic outputs beyond the nilpotent and finite-fibre gates (HC4MYIC2)

Let

\[
         P=a(T)z^3+b(T)z^2+c(T)z+d(T),\qquad a\ne0,
              D(T)=b(T)^2-3a(T)c(T).                     \tag{7}
\]

Every closed stable embedding carrying K through (T,xy-P) forces
K invertible if either of these independent conditions holds:

- a is square-free and a divides b;
- D!=0 and its square-free part has degree at least two.

All coefficient and embedding degrees are unrestricted, and no
special K1-fibre hypothesis is required. These conditions do not
assume that the generic root group is abelian or nilpotent.

A reducible cubic over L has a rational root. The first condition
is excluded by Section 2. For the second, (4) and its derivative are

\[
 Q=U^3+bU^2+acU+a^2(d+V),\qquad Q_U=3U^2+2bU+ac.
\]

Under Y=3U+b, the critical curve is Y^2=D(T). Write
D=H(T)^2*D0(T), with D0 square-free and nonconstant of degree n.
Its affine normalization is y^2=D0(T), with y=Y/H integral over
the original curve. If n=2 this is a rational curve with two distinct
points at infinity; no polynomial A1 can dominate it, because a
nonconstant P1 map cannot put the preimages of two distinct points
inside its single omitted point. If n>=3, Riemann--Hurwitz for the
double cover gives genus floor((n-1)/2)>=1, likewise excluding an
A1 parametrization. This contradicts Section 2. Hence the cubic
stays irreducible over L, and the retained root-transitivity/global
filtration argument gives the inverse for K.

Here D is the **critical-point discriminant**, not the discriminant
of P+V in z. For example P=T*z^3+z satisfies the first condition
and has S3 generic root group: its root discriminant is
-4T-27T^2*V^2, nonsquare over C(T,V). This closes a partial leading
degree drop for which the single-slope criterion in Section 1 fails.

For a complementary family take, for every integer e>=1,

\[
           P=T*z^3+z^2-\tfrac13 T^{2e-1}z.
\]

Then D=1+T^(2e) is square-free of degree 2e, so the second condition
closes every member. The first condition fails since T does not
divide b=1, and T=0 drops degree from three only to two, so Section 1
does not close this example either. The generic group is S3: P+V
is irreducible and its root discriminant has two simple V zeros,
as its V-discriminant is 16*D^3!=0. Thus this is a genuine new
nonnilpotent family, with unbounded T degree and no fibre hypothesis.

The genuine embedding controls j=(s,1,t+P(s,w),w), with an unused
v, still have inverse (T,y-P(T,z),z,v). The independent quadratic
T*v+y*z+x^2 pulls back to 1+s*v+t*w+w*P(s,w), with actual unit
Hessian and the retained explicit gradient inverse. These examples
have invertible gradients and do not construct nonautomorphic K.

The cubic residual is now sharply constrained. Apart from cases
closed by another retained gate, a nonautomorphic realization would
require a rational root label, positive localized-unit defect, S3
root monodromy, and D with square-free part of degree at most one;
it must also fail the square-free-leading/divisible-quadratic condition.
If D is a nonzero square, (6) strengthens the required unit defect
to at least two, since Q_U has two distinct rational critical roots.
The case D=0 is already cyclic over C(T,V) after rational translation,
and is closed by HC4MYIA2. The remaining conditions are necessary,
not a construction or a proof that the residual is nonempty.
General higher-degree intransitive factorizations not closed by (2),
other output geometries, both plane endpoints and unrestricted HC4/JC2
remain open.

The [source record](sources/valuation-critical-import.json) identifies
the new nonproperness input and retained Hurwitz input explicitly.
The [checker](verify_valuation_critical.py) and
[receipt](receipt-valuation-critical.json) give exact finite controls
separate from these unrestricted proofs. No formal verification,
independent replay, external review or novelty priority is asserted.
