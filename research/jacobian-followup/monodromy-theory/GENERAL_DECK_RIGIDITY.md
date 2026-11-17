# Regular source deck rigidity for every polynomial Keller map

**KMFI47.** Let F:A^m->A^m be a polynomial Keller map over any
characteristic-zero field. If a polynomial endomorphism sigma of
the source satisfies F composed with sigma=F, then sigma is the
identity. In particular the polynomial source deck group is
trivial in every dimension, after every characteristic-zero
constant extension and every identity stabilization.

This does **not** say that the generic deck group
Aut(k(source)/k(target))=N_G(H)/H is trivial.
Nonidentity generic deck transformations, when present, are
rational maps which cannot extend regularly to the whole source.
It does not restrict generic monodromy to the identity and does
not settle the inverse action problem. KMFI30 already supplies
actual imprimitive Keller examples with generic deck C2.

[Exact checker](verify_regular_deck.py) · [receipt](regular_deck_receipt.json) ·
[primary-source pin](sources/regular-deck-manifest.json).

## 1. A polynomial endomorphism over F has finite order

F is dominant and generically finite. F sigma=F makes sigma dominant.
On rational function fields it induces an injective endomorphism
of the finite extension L/K, fixing K=k(F).
As an injective endomorphism of a finite-dimensional K-vector
space it is surjective. It is an automorphism of the finite
separable extension L/K and has finite order.
Thus a power of sigma is the identity on the function field,
and hence as a polynomial morphism.
In particular sigma is an automorphism with polynomial inverse.

The group of such automorphisms injects into Aut_K(L), a finite
group. It is enough to rule out one nontrivial finite subgroup.

## 2. Etaleness makes the deck action free

All coefficients, the identities F sigma=F and sigma^r=id, and a
nonzero coefficient witnessing sigma!=id descend to a finitely
generated characteristic-zero field. Embed it into C.
No assumption is lost; a contradiction over C suffices.

At a fixed point of a deck automorphism, the complex inverse
function theorem makes F biholomorphic on a neighborhood.
The equation F sigma=F then forces sigma to be the identity on
a neighborhood, hence identically as a polynomial map.
Consequently every nonidentity deck automorphism has no fixed
point. Any nontrivial finite subgroup H acts freely on A^m_C.

Its affine invariant quotient Q=A^m_C/H exists, and the free
finite action gives a finite etale covering of degree |H|.
In particular Q is a finite-type complex variety.

## 3. Euler characteristic rules out the free quotient

The compactly supported Euler characteristic of complex affine
space is1, and that of a finite-type complex variety is an integer.
For a finite etale cover of degree d, it is multiplied by d.
The latter standard fact is explicitly recorded in
[Tony Feng's notes on Ravi Vakil's lectures, Section12.2](https://math.berkeley.edu/~fengt/EquivAG.pdf),
with a retained versioned source pin.

Apply it to the free quotient:

\[
 1=\chi_c(\mathbf A^m)=|H|\,\chi_c(Q).
\]

An integer chi_c(Q) makes this impossible for |H|>1.
Hence no nontrivial sigma exists. This proves KMFI47.

The argument also applies to a smooth affine source X with
chi_c(X)=1 and a generically finite etale map to any smooth target.
More generally a finite regular deck subgroup must have order
dividing chi_c(X), if that integer is nonzero.
No assertion is made that the rational generic deck group has
that divisibility property.

## 4. Necessary hypothesis controls

In characteristic p, F(x)=x^p-x is Keller with derivative -1
and admits every translation x->x+a,a in F_p as a regular deck
automorphism. Characteristic zero is essential.

For the non-etale map x->x^n on A1 in characteristic zero,
multiplication by an n-th root of unity is a regular deck
automorphism with a fixed ramification point.

On G_m, x->x^n is finite etale and has regular deck group mu_n.
Here chi_c(G_m)=0, so the Euler equation is consistent.
Affine-space source geometry is also essential.

These exact controls test the boundaries of the proof; they do
not replace its all-dimensional argument.
No novelty priority, formal verification or independent whole-proof
review is asserted. The full Keller monodromy classification is open.
