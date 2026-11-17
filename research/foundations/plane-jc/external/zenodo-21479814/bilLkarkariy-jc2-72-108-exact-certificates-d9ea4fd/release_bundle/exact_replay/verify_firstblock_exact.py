#!/usr/bin/env python3
import sympy as sp
import exact_core as ec
# 1. Exact triangular J4 reconstruction.
expr=sp.expand(2*ec.A0*sp.diff(ec.D0,ec.t)-3*sp.diff(ec.A0,ec.t)*ec.D0-ec.t**2)
sub=sp.cancel(expr.subs(ec.dsubs))
coeff=[sp.factor(sp.Poly(sub,ec.t).nth(k)) for k in range(sp.Poly(sub,ec.t).degree()+1)]
nonzero=[q for q in coeff if q!=0]
assert len(nonzero)==6
# 2. Lex basis gives exact zero to every original residual.
for q in ec.res:
    qq=sp.together(q.subs(ec.LEX)).as_numer_denom()[0]
    rem=sp.rem(sp.Poly(qq,ec.U,domain=sp.QQ),sp.Poly(ec.Hprim,ec.U,domain=sp.QQ)).as_expr()
    assert sp.expand(rem)==0
# 3. H sparse shape and irreducibility over Q.
Hp=sp.Poly(ec.Hprim,ec.U,domain=sp.QQ)
fac=sp.factor_list(Hp.as_expr(),ec.U)
assert len(fac[1])==1 and sp.Poly(fac[1][0][0],ec.U).degree()==35 and fac[1][0][1]==1
print('J4_D_solves',len(ec.dsubs),'residuals',len(nonzero))
print('LEX_relations',len(ec.rels),'H_degree',Hp.degree(),'H_terms',len(Hp.terms()))
print('H_irreducible_Q',True)
print('all_firstblock_residuals_vanish_mod_lex',True)
