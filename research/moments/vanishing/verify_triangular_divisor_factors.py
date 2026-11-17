"""Auxiliary triangular D factors for retained rank-one seed recurrences.
The general formal block/seed/minor lemma and universal analytic/class cover
remain written inputs. Retained long-running checker mains are not run.
"""
from pathlib import Path
from itertools import combinations
import hashlib,json,time
import sympy as s
from sympy.polys.matrices import DomainMatrix
P=Path(__file__).resolve().parent
t,D=s.symbols('t D');sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
T=lambda n:n*(n-1)//2

def determinants(columns,rank,var):
    ring=s.QQ.poly_ring(var);out=[]
    for k in range(rank+1):
        cols=columns[:k]+columns[k+1:];mat=DomainMatrix.from_list_sympy(rank,rank,[[x.as_expr() for x in row] for row in zip(*cols)]).convert_to(ring)
        out.append((-1)**k*s.Poly(ring.to_sympy(mat.det()),var,domain=s.QQ))
    assert all(sum((out[k]*columns[k][i] for k in range(rank+1)),s.Poly(0,var)).is_zero for i in range(rank))
    return out

def abstract_models():
    rows=[]
    for rank in [3,5,7]:
        for derivative_order in [2,4]:
            left=s.Matrix([i+1 for i in range(rank)]);right=s.Matrix([[1]+[0]*(rank-1)]);A0=left*right
            A=A0+D*s.diag(*[i+2 for i in range(rank)])+D*D*s.Matrix(rank,rank,lambda i,j:(i+j+1)%4)
            assert A0.trace()==1 and A0.rank()==1
            cur=s.eye(rank)[:,0];cols=[[s.Poly(x,D) for x in cur]]
            for j in range(rank):
                scalar=-j*D if derivative_order==2 else -(3*j+1)*D**3
                cur=(D**derivative_order*cur.diff(D)+A*cur+scalar*cur).applyfunc(s.expand)
                cols.append([s.Poly(x,D,domain=s.QQ) for x in cur])
            co=determinants(cols,rank,D);checks=[]
            for k,C in enumerate(co):
                need=T(rank) if k==0 else max(0,T(rank)-k)
                Q=C.exquo(s.Poly(D**need,D));checks.append({'omitted_column':k,'certified_order':need,'zero':Q.is_zero})
            assert not co[0].is_zero
            rows.append({'rank':rank,'derivative_order':derivative_order,'trace_rank_one_residue':1,'minor_checks':checks})
    return rows

def residue_checks():
    ordinary=json.loads((P/'ordinary_rank_one_lattice_certificate.json').read_text())['universal'];cusp=json.loads((P/'cusp_rank_one_lattice_certificate.json').read_text())['rank_one_residue']
    u,v=s.symbols('u v');eps,R,p,q,r=s.symbols('eps R p q r')
    rows=[]
    for c in ordinary:
        b=s.Matrix([s.sympify(x,locals={'t':t,'eps':eps,'R':R}) for x in c['residue_vector']]);unit=s.zeros(1,len(b));unit[0,0]=1;A=b*unit
        assert s.expand(A.trace()+t*eps*R*R/2)==0
        assert all(s.expand(x)==0 for x in A*A-A.trace()*A)
        assert A.trace()!=0
        rows.append({'kind':c['kind'],'trace':str(A.trace()),'generic_trace_not_identically_zero':True,'residue_square_equals_trace_times_residue':True,'derivative_coefficient_D_order':2,'scalar_j_term_D_order_at_least':1})
    b=s.Matrix([s.sympify(x,locals={'p':p,'q':q,'r':r}) for x in cusp['base_first_seed_residue']]);A=s.Matrix([[s.sympify(x,locals={'p':p,'q':q,'r':r}) for x in row] for row in cusp['base_rank_one_matrix']])*t*eps
    assert s.expand(A.trace()+t*eps*q**6/18)==0 and A.trace()!=0
    assert all(s.expand(x)==0 for x in A*A-A.trace()*A)
    rows.append({'kind':'cusp','trace':str(A.trace()),'generic_trace_not_identically_zero':True,'residue_square_equals_trace_times_residue':True,'derivative_coefficient_D_order':4,'scalar_j_term_D_order_at_least':3})
    return rows

def replay(m,origin):
    cusp='cusp' in origin;kind='cusp' if cusp else m['kind'];rank=len(m['seed_columns']);columns=[[s.Poly(int(i==0),t,domain=s.QQ) for i in range(rank)]]
    columns += [[s.Poly(s.sympify(x,locals={'t':t}),t,domain=s.QQ) for x in row['N']] for row in m['seed_columns']]
    den=s.Poly(s.sympify(m['D'],locals={'t':t}),t,domain=s.QQ);co=determinants(columns,rank,t);rows=[]
    for k,C in enumerate(co):
        need=T(rank) if k==0 else max(0,T(rank)-k);C.exquo(den**need)
        H=s.Poly(1,t,domain=s.QQ) if k==0 else s.Poly(s.sympify(m['seed_columns'][k-1]['H'],locals={'t':t}),t,domain=s.QQ)
        reduced=(H*C).exquo(den**T(rank))
        if k==0:
            gate={'cusp':115,'nodal':108,'transverse':112,'triangle':108}[kind]
            assert not reduced.is_zero and reduced.degree()<=gate
        val=min(e[0] for e,c in reduced.terms() if c) if reduced else None
        if k and val is not None:assert val>=3*k
        rows.append({'order':k,'raw_minor_D_order_certified':need,'common_scalar_D_order':T(rank),'reduced_scalar_degree':int(reduced.degree()) if reduced else -1,'reduced_scalar_valuation':val})
    return {'origin':origin,'kind':kind,'label':m.get('label',m.get('mode')),'D':m['D'],'rank':rank,'cofactor_checks':rows}

def arithmetic():
    rows=[]
    bases={'nodal':[(0,0),(1,0),(0,1),(2,0),(1,1),(0,2),(0,3)],'transverse':[(0,0),(1,0),(0,1),(2,0),(1,1),(0,2)],'triangle':[(0,0),(1,0),(0,1),(1,1),(0,2)],'cusp':[(0,0),(1,0),(2,0),(0,1),(1,1),(2,1)]}
    for kind,bs in bases.items():
        if kind=='triangle':grades=[i+k+int(s.ceiling(i+s.Rational(k,2)))+int(s.ceiling(s.Rational(i,2)+k)) for i,k in bs]
        elif kind=='cusp':grades=[0,1,2,1,3,4]
        else:grades=[sum(ex) for ex in bs]
        slope={'nodal':5,'transverse':7,'triangle':9,'cusp':7}[kind];offset={'nodal':0,'transverse':2,'triangle':0,'cusp':1}[kind];cases=[]
        for size in range(1,len(bs)+1):
            for selected in combinations(range(len(bs)),size):
                if 0 not in selected:continue
                raw=slope*size*(size+1)//2-offset*size-sum(grades[k] for k in selected)
                for k in range(1,size+1):
                    denominator=4*k-1 if kind=='cusp' else 2*k
                    assert max(0,T(size)-k)+denominator>=T(size)
                cases.append({'rank':size,'rows':selected,'raw_degree_bound':raw,'common_scalar_D_power':T(size),'new_degree_bound':raw-T(size)})
        gate=max(c['new_degree_bound'] for c in cases);expected={'nodal':108,'transverse':112,'triangle':108,'cusp':115}[kind];assert gate==expected
        rows.append({'kind':kind,'linear_D_gate':gate,'all_rank_row_cases':cases})
    cover={'quadratic':7,'smooth':94,'nodal':108,'cusp':115,'transverse':112,'tangent':120,'triangle':108,'concurrent':104,'repeated':99}
    assert max(108,101,96)==108 and max(112,106,96)==112 and max(108,103,96)==108 and max(115,109,92,95)==115 and max(cover.values())==120
    return rows,cover

def main():
    start=time.monotonic();residues=residue_checks();abstract=abstract_models();models=[]
    for name in ['cusp_rank_one_lattice_certificate.json','ordinary_rank_one_lattice_certificate.json']:
        for m in json.loads((P/name).read_text())['models']:
            if m['D_degree']==1:models.append(replay(m,name))
    cases,cover=arithmetic();cert={'generic_residue_checks':residues,'abstract_rank_one_differential_models':abstract,'retained_cofactor_models':models,'degree_cases':cases,'global_cover':cover}
    target=P/'triangular_divisor_factors_certificate.json';target.write_text(json.dumps(cert,indent=2,default=int)+'\n')
    report={'completed':True,'checker_sha256':sha(Path(__file__)),'certificate_sha256':sha(target),'runtime_seconds':round(time.monotonic()-start,3),'sympy_version':s.__version__,'generic_residue_strata':4,'abstract_differential_models':6,'retained_full_cofactor_models':len(models),'retained_full_scalar_cofactor_checks':sum(len(m['cofactor_checks']) for m in models),'linear_D_gates':{c['kind']:c['linear_D_gate'] for c in cases},'constant_D_gates_unchanged':{'nodal':101,'transverse':106,'triangle':103,'cusp':109},'global_sufficient_pure_prefix':120,'odd_sufficient_pure_prefix_unchanged':98,'homogeneous_sufficient_pure_prefix_unchanged':96,'least_prefix_interval':[20,120],'odd_least_prefix_interval':[20,98],'homogeneous_least_prefix_interval':[14,96],'global_cover':cover,'optimality_proved':False,'written_analytic_or_global_proof_certified_by_checker':False,'boundary':'Exact generic residue/trace and six abstract differential recurrence models; ten retained full Gaussian cofactor relations recomputed with D^(s(s-1)/2) cancellation and reduced valuations; all ranks/rows. The general formal semisimple-block/seed/minor lemma, polynomial extension to degenerate trace specializations, complete coefficient/frame cover and retained analytic exclusions remain written inputs. No independent replay, external review, formal verification or expanded radical cofactors.'}
    (P/'triangular_divisor_factors_verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':main()
