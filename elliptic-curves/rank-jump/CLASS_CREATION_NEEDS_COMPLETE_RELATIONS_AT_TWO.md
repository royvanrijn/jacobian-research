# Class creation needs complete relations at two

The generic-half-ideal result gives a sufficient incidence constructor,
but its observed Artin rank is **not the number of quadratic characters**.
A change from ten cyclic factors of order four to ten factors of order two
changes that detection rank from zero to ten without increasing the number
of independent quadratic characters. These are abstract comparison groups,
not a demonstrated transition between two number fields in the panel.

This separates two possible arithmetic events that the next matched test
must distinguish: a larger space of unramified covers, and a change in how
existing characters detect order-two ideals. Neither event establishes
rational solubility.

## Three dimensions, with different meanings

For the cubic field K_t and the complete retained bad set S, put
C_t=Cl(O_K,S). Write its two-primary part as

    C_t(2) = direct sum_i Z/2^(a_i),  a_i >= 1.

Define c=#i, rho4=#{i:a_i>=2}, and e1=#{i:a_i=1}. Then

    dim U_t = c,
    dim(2C_t/4C_t) = rho4,
    rank(C_t[2] -> C_t/2C_t) = e1 = c-rho4.

The first equality uses the established strict class-field identification.
For the last equality, the order-two element in a cyclic factor is
2^(a_i-1): it survives modulo two exactly when a_i=1. Consequently the
full Artin evaluation of quadratic characters on order-two ideals has
rank e1. Restricting its columns to the generic half ideals gives

    d = dim image(Phi_S(G) -> C_t/2C_t) <= e1.

The earlier certificates give d>=10,8,6 using oracle characters, hence
e1>=10,8,6. They do **not** determine c or rho4. For example, if c=10+epsilon
on the first control, they imply rho4<=epsilon. The prior sufficient
criterion remains valid: with b=dim Phi_S(G) and g0=dim(G intersect U),
generic classes plus norm-square S-units contain at least max(0,b-g0)
extra strict classes. Since b>=d, the certified lower bounds still apply.
No earlier ideal-independence certificate is withdrawn.

This terminology concerns arithmetic detection by Artin characters. It
has no connection to point-search visibility or the half-lattice chart budget.

## An exact relation condition at a specialization

Suppose n chosen ideals give a certified presentation of C_t at the prime
two. Let the columns of the integer matrix R_t be the complete relations
at two, including the classes of primes inverted in S. Precisely, the
natural map from Z^n/im R_t to C_t must induce an isomorphism on their
two-primary parts. Odd kernel and cokernel are harmless here. A list of
verified principal relations alone does not prove this hypothesis.

Then

    c(t) = n - rank_F2(R_t mod 2).

If A_t records the generic half ideals in this presentation, with their
order-two identities checked, then

    d(t) = rank_F2([R_t | A_t]) - rank_F2(R_t).

For a full-rank finite presentation, let N_q count vectors x modulo q
with R_t^T x=0. The two counts suffice:

    log2 N_2 = c,
    log2 N_4 = c+rho4,
    e1 = 2 log2 N_2 - log2 N_4.

Thus the required endpoint can be extracted from a presentation complete
at two; computing odd class-group invariants and expanding a full unit
basis is not logically necessary. Producing that certified presentation
for the research fields remains a computational problem.

For a global pool F_t containing G_t with dim(F_t/G_t)<=q, the equation-only
condition

    n-rank_F2(R_t) >= g0(t)+q+r

forces at least r strict classes outside F_t. Indeed,
dim(U_t intersect F_t)<=g0(t)+q, while dim U_t=c(t). The retained uniform
bound is q=2; the published original R17 family has the stronger q=0 from
[the completed global calculation](THE_LAST_GLOBAL_CLASS_IS_ABSENT.md).
These are Selmer-incidence classes; rational points still require a
separate solubility argument.

This is an exact certificate condition on the arithmetic at t, **not** a
derived polynomial condition on t or an explanation for why its matrix
loses rank. Matrices from different number fields have no canonical common
ideal basis. An actual specialization mechanism must supply the relations
and prove their change, rather than interpret an unfinished computation as
a vanishing relation.

## Why incomplete relations cannot certify creation

If R_partial contains only some true relations and the ideal generators
surject onto the class group, its cokernel surjects onto C_t. Missing
relations can therefore inflate the apparent two-rank. The single matrix
[2] predicts C2, while adding the missing relation [1] gives the trivial
group. Without certified generation even this upper-bound interpretation
is unavailable. Exact verification of each retained relation does not
close the gap. Nor is ordinary lattice saturation the right operation:
saturating a full-rank relation lattice inside Z^n destroys its torsion.

The relevant completeness is equality with the *true* relation lattice
after tensoring with Z_2. A full certified class group is sufficient, but
not necessary. PARI distinguishes conditional BNF output from certification;
see its [class-group documentation](https://pari.math.u-bordeaux.fr/dochtml/html-stable/General_number_fields.html#bnfcertify).
The previous two timeouts returned no such presentation.

## Bounded test and the next arithmetic endpoint

The [protocol](CLASS_CREATION_RELATION_GATE_PROTOCOL.json) tested 49 abstract
two-factor groups with factors drawn from 1,2,3,4,6,8,12. Independent exhaustive
enumeration of characters modulo two and four agrees with Smith-form
calculations. Unimodular changes of generators and relations preserve every
reported invariant. Explicit ten-factor examples verify the detection
counterexample above; a missing-relation example verifies the false-positive
failure. These computations construct no new arithmetic classes.

The [artifact](../../artifacts/generated-results/elliptic-curves/rank_jump_class_creation_relation_gate_v1.json)
also audits all sixteen frozen panel rows and the three earlier ideal
certificates. Equation-only absolute c and rho4 remain UNKNOWN in these
inputs. The rank-derived strict lower bound for fresh case-02 is nine;
case-03's lower bound is zero, not an upper bound. Thus a larger absolute
class rank on the high fibre has not yet been measured. The established
capacity theorem instead proves new classes relative to its inherited
global pool, which is a different and valid conclusion.

The next arithmetic test remains the frozen 103b2 pair, with no point inputs:
obtain a presentation certified at two, then compute (c,rho4,d,g0) before
unmasking outcomes. A certified small c on the low fibre and large c on the
high would establish a paired incidence difference. Equal c with different
d would direct attention to the elementary factors and generic correction
map, without establishing greater absolute cover capacity. Equal incidence
with different rational gain would instead require the solubility stage.
No larger BNF retry or parameter sweep is launched by this note.

The strongest unresolved implication is still: an explicit condition on t
forces enough independent unramified, S-split covers outside the global
pool. The finite-presentation formula specifies what must be certified;
it does not close that implication. Neither the new point-derived
determinant-1092 parent candidate nor a larger Artin rank supplies the
missing point-independent proof. No candidate-selection feature is ready
for Agent 1 from this calculation.

Replay:

```sh
timeout 30 sage -python elliptic-curves/rank-jump/class_creation_relation_gate.py check
```
