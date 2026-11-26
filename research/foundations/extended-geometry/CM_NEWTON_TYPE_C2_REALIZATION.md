# A geometric realization of the rank-four type \(C_2\) case

**Theorem.** Let \(F\) be a totally real field of degree \(e\ge2\)
and put \(E=F(i)\). There is a principally polarized, absolutely
simple abelian variety \(A/K\) over a number field \(K\subset\mathbb R\)
of dimension \(4e\) with full geometric rational endomorphism
algebra \(E\). It has one conjugate Hodge pair of
signature \((2,2)\) and \(e-1\) definite pairs. For an auxiliary
prime \(\ell\), its connected derived monodromy over
\(\overline{\mathbb Q}_\ell\) is a product
of \(e\) copies of \(\mathrm{Sp}_4\). Each factor acts by its
standard representation on the two conjugate coefficient blocks
and trivially on the others. The connected projective spectrum
field is \(F\) and the minimal endomorphism field is the quadratic
extension \(KE/K\). Consequently \(A\) has \(\mu\)-ordinary
reduction at natural density one over \(K\) and ordinary
reduction at natural density \(1/2\).

For \(F=\mathbb Q(\sqrt2)\), the construction below uses the
quaternion algebra \((-1,\sqrt2)_F\). It produces an eightfold by
real specialization on a three-dimensional Shimura variety. The
argument establishes existence of a number-field fiber and does
not compute equations for that fiber. Its derived block
images have type \(C_2\), so the example falls under the
[linked rank-four theorem](CM_NEWTON_LINKED_RANK_FOUR_PFAFFIAN.md)
with totally real projective spectrum field.

## Quaternionic and polarized representations

Choose \(a\in F\) positive at one real embedding \(v_*\) and
negative at every other embedding, using weak approximation.
Let \(D=(-1,a)_F\) with generators \(i,j\) satisfying
\(i^2=-1\), \(j^2=a\), and \(jz=\bar z j\) for \(z\in E=F(i)\).
Let \(x\mapsto x^*\) denote canonical quaternion conjugation.
At \(v_*\) the algebra \(D\) is split, while at the other real
places it is Hamilton's quaternion algebra. On the right
\(D\)-module \(V=D^2\) take the hermitian form
\[
 q(x,y)=x_1^*y_1+x_2^*y_2 .
\]
Write each coordinate as \(u+jv\) with \(u,v\in E\). Projection
to the \(E\)-summand of \(D=E\oplus jE\) gives the hermitian
form \(H=\operatorname{pr}_E q\) on the rank-four \(E\)-space
\(V\). In these coordinates its matrix is
\(\operatorname{diag}(1,-a,1,-a)\). It is positive definite
at the other real places and has signature \((2,2)\) at \(v_*\).
The rational alternating form
\[
          \psi(x,y)=\operatorname{Tr}_{E/\mathbb Q}(-i H(x,y))
\]
has the property that right multiplication by \(e\in E\)
is adjoint to right multiplication by \(\bar e\).

Let \(G_s=\operatorname{Res}_{F/\mathbb Q}\mathrm U(q)\), where
\(\mathrm U(q)\) consists of right \(D\)-linear isometries.
Its action on \(V\) preserves \(H\) and \(\psi\) and commutes
with the right \(E\)-action. It has absolute type \(C_2\)
over each embedding of \(F\). To see the representations,
split \(D\) as \(M_2\) over an algebraic closure. The
involution on \(M_2(D)\) adjoint to \(q\) becomes the
symplectic involution with matrix the direct sum of
two invertible skew matrices of size two. Its isometry
group is \(\mathrm{Sp}_4\) acting on the column space
of \(\operatorname{Hom}(\mathbb C^2,\mathbb C^4)\).
Each right \(E\)-embedding picks one of the two columns.
Over \(F_{v_*}=\mathbb R\) the same calculation gives
\(\mathrm{Sp}_4(\mathbb R)\). At each other real place
it gives the compact quaternionic isometry group
\(\mathrm{Sp}(2)\) since \(q\) is positive definite there.
[Oda–Schwermer, §1](https://www.jstage.jst.go.jp/article/tmj/61/1/61_1_83/_pdf).

Let
\[
 T=\{z\in\operatorname{Res}_{E/\mathbb Q}\mathbb G_m:
                 z\bar z\in\mathbb G_m\}
\]
act by right scalars, with \(\mathbb G_m\) on the right
embedded diagonally over all real places. Set
\(G=T G_s\subset\mathrm{GSp}(V,\psi)\).
The derived group of \(G\) is \(G_s\).
Define \(h:\mathbb S\to G_{\mathbb R}\) on real points
as follows. For \(z=b+di\in\mathbb C^\times\), let \(h(z)\)
act by left multiplication by \(b+di\) at \(v_*\), and by
right multiplication by \(b+di\) at the other real places.
Each of these maps is algebraic over \(\mathbb R\).
It lies in \(G(\mathbb R)\): at \(v_*\) the isometry part is
left multiplication by \(z/|z|\) after taking out the real
scalar \(|z|\), and all scalar parts have norm \(|z|^2\).
The Hodge structure defined by \(h\) is polarized by \(\psi\).
At \(v_*\), \(h(i)\) acts on each \((u,v)\) by
\((iu,-iv)\), so \(\psi(x,h(i)x)>0\) uses the two opposite
signs of \(H\). At other real places \(h(i)\) is right
multiplication by \(i\) and positivity follows from the
positive definite \(H\). This gives \((2,2)\) at \(v_*\)
and the definite types elsewhere.

Let \(X\) be the \(G(\mathbb R)\)-conjugacy class of \(h\).
Then \((G,X)\) is a Shimura datum of Hodge type. On the
noncompact factor at \(v_*\), conjugation by \(h(i)\)
is the Cartan involution of the Siegel datum for
\(\mathrm{Sp}_4(\mathbb R)\). Indeed left \(i\) on the
column representation is the standard compatible complex
structure for the symplectic matrix just described. On
the compact factors the adjoint action of \(h\) is trivial.
The corresponding adjoint Hodge weights are \(0,(1,-1),(-1,1)\),
where \(0\) denotes \((0,0)\). The rational adjoint group is
\(\mathbb Q\)-simple by restriction of scalars from \(F\),
and \(h\) is nontrivial on it. Its weight acts through
rational scalar matrices. These are Deligne's Shimura
conditions. The complex domain is the Siegel domain of
degree two, of dimension three.
[Milne, *Introduction to Shimura varieties*](https://www.jmilne.org/math/xnotes/svi.html).

## Real descent with the labeled coefficient action

On \(D^2\) let \(r\) act by quaternion conjugation
by \(j\) on each coordinate. In the right \(E\)-basis
\(1,j\), it acts by \((u,v)\mapsto(\bar u,\bar v)\).
Since \(j^2\in F^\times\), \(r^2=1\). The map \(r\)
is \(E\)-antilinear and \(\psi\)-anti-symplectic.
It normalizes \(G_s\), because the coordinatewise
quaternion automorphism preserves the canonical
involution and transports \(q\) to its conjugate.
It normalizes \(T\) through \(c\), and
\(rh(z)r^{-1}=h(\bar z)\) at every real place.
Thus \(\vartheta=\operatorname{Int}(r)\) on \(G\)
is an involution taking \(X\) to its conjugate.

The two rational eigenspaces of \(r\) are complementary
Lagrangians for \(\psi\). Choose dual rational bases and
let \(\Lambda\) be their integral span. It is \(r\)-stable
and self-dual. Choose a positive integer \(m\) such
that the \(c\)-stable order \(\mathbb Z+m\mathcal O_E\)
acts integrally on \(\Lambda\). Use full level \(N\ge3\)
prime to \(m\). The enclosing PEL group \(G_{\mathrm{PEL}}\)
of \(E\)-linear symplectic similitudes contains \(G\),
and the PEL datum contains \(h\). The \(G\)-level obtained
by intersecting with the PEL lattice stabilizer and its
principal level \(N\) subgroup is \(\vartheta\)-stable.
In characteristic zero this gives a canonical
Hodge-type morphism to a fine PEL moduli space with
the indicated order action, principal polarization
and level.
[Kottwitz, §5](https://webusers.imj-prg.fr/~christophe.cornut/M2/Ref/Kottwitz_PointsShiModp.pdf).

Work over a finite CM field \(k\subset\mathbb C\)
containing \(E\), both reflex fields and \(\zeta_N\).
One can take the normal closure of \(E\) times
\(\mathbb Q(\zeta_N)\), since both cocharacters are
defined over that normal closure. Indeed, over it the
real embedding factors and the \(E\)-orientations split.
In the split factor the Hodge cocharacter can be
written using the projectors for left \(i\) with
eigenvalues \(i\) and \(-i\), and at other factors
it is scalar on each oriented \(E\)-block.
The reflex fields are nonreal, since complex
conjugation swaps the weights on a definite oriented
block and neither group permutes the \(E\)-labels by
conjugation. Put \(k^+=k\cap\mathbb R\).

Apply the complex conjugation isomorphism for
canonical models followed by \(\vartheta^{-1}\).
This gives an effective descent datum to \(k^+\) on
the Hodge-type space. In fact the conjugation rule on
complex double quotients is
\([h',g]\mapsto[\vartheta^{-1}(\overline{h'}),
 \vartheta^{-1}(g)]\), where
\(\overline{h'}=h'\circ c_{\mathbb S}\).
Its square is the identity and the point \([h,1]\)
is fixed. This is the involution criterion of
Blasius–Guerberoff, applied to the explicit
\(\vartheta\). The canonical conjugation maps and
the involutions commute with the morphism into the
PEL space. Equivalently, this compatibility is checked
on the double quotients using the same \(r\) in both
groups, and then holds for the canonical morphisms.
[Blasius–Guerberoff, §§2.2–2.3 and Theorem 2.3.1](https://msp.org/ant/2017/11-10/ant-v11-n10-p02-p.pdf).

On the fine PEL functor the same descent replaces the
action \(\iota\) by \(\iota\circ c\) and the level
\(\eta:\Lambda/N\Lambda\to A[N]\) by \(\eta\circ r\)
on a conjugate object. This restores the determinant
condition, since \(\overline{P_e(T)}=P_{\bar e}(T)\)
for the signature polynomial, and restores the pairing
condition, since \(r\) reverses the alternating form
as complex conjugation reverses \(\zeta_N\).
It keeps the underlying polarized abelian scheme,
and both relabelings square to the identity.
Thus the PEL universal abelian scheme, forgetting the
labels and the level, descends and pulls back to
the descended Hodge-type space. Its labeled \(E\)-action
over \(k\) is conjugated by \(c\) in this descent.
The geometric component containing \([h,1]\) is
fixed by complex conjugation, so its stabilizer
in the Galois group defines a field \(k_0\supset k^+\)
inside \(\mathbb R\). The component over \(k_0\)
is smooth and geometrically connected and has this
real point.

## A number-field fiber with connected type \(C_2\)

The arithmetic subgroup of the complex uniformization
intersects \(G_s\) in an arithmetic subgroup. It is
Zariski dense in \(G_s\). The real group \(G_s(\mathbb R)\)
has one noncompact factor, \(\mathrm{Sp}_4(\mathbb R)\),
and compact factors at the other real places. The
arithmetic subgroup projects to a lattice in the
noncompact factor, since the remaining factors are
compact. Borel density gives a surjective Lie
projection of its rational Zariski closure onto the
noncompact Lie algebra. This closure contains the
factor: the kernel of its projection to the compact
factors must project onto the noncompact simple Lie
algebra, since otherwise that algebra would be a
quotient of a compact Lie algebra. The rational
closure then contains the conjugate factors as well.
Betti–étale comparison gives the full \(G_s\) in the
geometric algebraic monodromy of the universal family.
[Borel–Harish-Chandra, Theorem 7.8](https://www.mathi.uni-heidelberg.de/~wienhard/retreat10/references/borel_harishchandra_2.pdf),
[Borel, density theorem](https://www.jstor.org/stable/1970150),
[Noot, §1.8](https://www.numdam.org/article/CM_1995__97_1-2_161_0.pdf).

To retain this group on a real number-field fiber, let
\(\mathcal G\) be the compact arithmetic image of the
family on the auxiliary Tate module. Take a finitely
generated open normal pro-\(\ell\) subgroup \(U\).
Its Frattini subgroup \(\Phi(U)\) is open and
characteristic. A closed subgroup of \(\mathcal G\)
surjecting onto \(\mathcal G/\Phi(U)\) intersects \(U\)
in a subgroup surjecting onto \(U/\Phi(U)\), so it
contains \(U\) by the pro-\(\ell\) Burnside basis theorem
and therefore equals \(\mathcal G\).

Take the connected finite étale cover for
\(\mathcal G/\Phi(U)\) and compose it with an étale
chart of the descended component over affine space.
Shrink the chart and target so that this composite
is finite étale. A neighborhood of the given smooth
real point meets the chosen Zariski open. The real
inverse function theorem gives a local real section
on an open set in the target. Hilbert irreducibility
with weak approximation gives a \(k_0\)-rational
parameter in this open set whose fiber in the
composed cover is connected. Its image in the
chart fiber is a single closed point and has
residue field \(K\subset\mathbb R\) by the local
section. Connectedness of the cover fiber gives
surjectivity onto the Frattini quotient and hence
onto the full arithmetic image at this fiber.
[Noot, Proposition 1.3](https://www.numdam.org/article/CM_1995__97_1-2_161_0.pdf),
[Serre, *Topics in Galois Theory*, Theorem 3.5.3](https://indico.math.cnrs.fr/event/11410/attachments/4760/7351/SerreTopicsGaloisTheory.pdf).

At this fiber the connected \(\ell\)-adic group
contains \(G_s\) and is contained in \(G_{\mathbb Q_\ell}\).
The latter inclusion follows because the Hodge
cocharacter of each fiber lies in \(G\) and Deligne's
absolute Hodge theorem puts the connected \(\ell\)-adic
group inside the corresponding Mumford–Tate group.
The derived group is therefore exactly \(G_s\).
[Deligne, Proposition 2.9 and Theorem 2.11](https://www.jmilne.org/math/Documents/Deligne82.pdf),
[Pink, Theorem 5.6](https://people.math.ethz.ch/~pink/ftp/AMT-v3.pdf).

Over \(\overline{\mathbb Q}_\ell\) the only possible
extra linear intertwiners under the connected group
are between conjugate standard blocks over the same
embedding of \(F\). Indeed different embeddings of
\(F\) see independent factors and cannot intertwine,
and each block is irreducible under its own factor.
If a conjugate pair were intertwined, it would give
isomorphic local representations after a finite
extension defining the actions and connected
monodromy. Choose the embedding of this extension
into \(\mathbb C_\ell\) so that the chosen pair
corresponds to a definite pair of the original
signature. This is possible by varying the local
embedding of a normal closure of \(E\). The
two Hodge–Tate weight lists are then all zero and
all one, so the representations cannot be isomorphic.
Faltings's theorem now gives geometric endomorphism
algebra of dimension \(2e\), which is exactly \(E\)
by the given \(E\)-action. A field has no nontrivial
idempotent, so the fiber is absolutely simple.
[Faltings, Theorems 3–4](https://math.uchicago.edu/~drinfeld/Deligne%27s_conjecture_Manin_conf/Faltings_argument/Faltings.pdf),
[Faltings, *\(p\)-adic Hodge theory*](https://public.websites.umich.edu/~bhattb/almost_purity_2011/Faltings_-_p-adic-Hodge-theory.pdf).

Conjugate blocks have identical projective spectra
on connected monodromy because they see the same
standard factor and the center contributes only
scalars. Blocks over different embeddings of \(F\)
have distinct projective spectra, as one may vary
their independent derived factors. Thus the connected
projective spectrum field is \(F\). The labeled
action is defined over \(Kk\), which is quadratic
over \(K\), and \(Kk=KE\) since \(k^+\subset K\)
and \(E\) is nonreal. The descent acts nontrivially
on \(E\) by \(c\), so this is exactly the
endomorphism field.

Now the restriction-kernel hypothesis of the
[linked rank-four theorem](CM_NEWTON_LINKED_RANK_FOUR_PFAFFIAN.md)
holds with \(\Gamma=\{1,c\}\), and each derived
block image is simple. There is one balanced class
since \(E_0=F\) and there is one balanced pair.
Thus \(\mu\)-ordinarity has density one. At the
identity class of \(\operatorname{Gal}(KE/K)\) the
reference slopes are ordinary. At the nontrivial
class the definite coefficient places have
reference slope \(1/2\), by the half-power comparison
for \(c\). Since \(e\ge2\) there is at least one
such place. Therefore the ordinary primes have
density \(1/2\) by finite Chebotarev.
