#!/usr/bin/env sage-python
"""Exact T-to-minus-T action on the supplied generic section span."""
import sys,argparse
from pathlib import Path
from importlib.machinery import SourceFileLoader
from sage.all import QQ,ZZ,matrix,vector,lcm
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';sys.path.insert(0,str(CAS))
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
source=SourceFileLoader('mestre_generic_height',str(CAS/'mestre_generic_height_audit.sage')).load_module()
D=ROOT/'artifacts/local/elliptic-curves/mestre-generic-height-v1'
OUT=ROOT/'artifacts/generated-results/elliptic-curves/mestre_base_involution_v1.json'

def compute():
    paths=[Path(__file__).resolve(),CAS/'mestre_generic_height_audit.sage'];rows=[]
    for u in cert.read(D/'protocol.json')['outer_parameters']:
        path=D/('u'+u)/'height.json';d=cert.read(path);paths.append(path)
        E,A,B,P=source.construct(u);geometry=source.setup_height(A,B);T=A.parent().gen()
        indices=d['seed_indices'];basis=[P[i] for i in indices]
        G=matrix(QQ,d['covariant_height_gram']).matrix_from_rows_and_columns(indices,indices)
        action=matrix(QQ,11,11);relations=[]
        for j,point in enumerate(basis):
            conjugate=E([c(-T) for c in point.xy()]);h=source.height(conjugate,A,B,geometry)[0]
            if h!=G[j,j]:raise ArithmeticError('base involution changed section height')
            pairings=vector(QQ,[(source.height(P+conjugate,A,B,geometry)[0]-G[i,i]-h)/2 for i,P in enumerate(basis)])
            coeff=G.solve_right(pairings);den=lcm([c.denominator() for c in coeff]);word=[ZZ(den*c) for c in coeff]
            if den*conjugate!=sum((a*P for a,P in zip(word,basis)),E(0)):raise ArithmeticError('exact involution action relation failed')
            action.set_column(j,coeff);relations.append({'target_index':j,'multiplier':int(den),'coefficients':list(map(int,word))})
        if action*action!=matrix.identity(QQ,11) or action.transpose()*G*action!=G:raise ArithmeticError('involutive isometry check failed')
        plus=11-(action-matrix.identity(QQ,11)).rank();minus=11-plus
        rows.append({'outer_u':u,'raw_covariant_basis_indices':indices,'action_matrix':[[str(c) for c in r] for r in action.rows()],
          'exact_group_relations':relations,'invariant_rank':int(plus),'anti_invariant_rank':int(minus),
          'quotient_geometric_MW_rank':6,'geometric_K3_MW_lower_bound':int(6+minus),
          'geometric_NS_lower_bound':int(2+5+6+minus)})
        print('u'+u,'SIGMA + / -',plus,minus,'GEOMETRIC NS >=',2+5+6+minus,flush=True)
    return {'schema':'elliptic-curves.mestre-base-involution.v1','status':'PASS','rows':rows,
      'sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in paths},
      'scope':'Every column of the T->-T action is checked by exact rational-function group law, and the matrix is an involutive height isometry. The invariant geometric part is the rank6 rational elliptic quotient over Qbar(s), s=T^2. Combining it with the certified anti-invariant span gives a geometric MW/NS lower bound, not additional rational sections, a full NS marking, or a new fibration.'}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args();result=compute()
    if a.check:assert result==cert.read(OUT)
    else:
        if OUT.exists():raise FileExistsError('preserve involution proof')
        checkpoint(OUT,result)
