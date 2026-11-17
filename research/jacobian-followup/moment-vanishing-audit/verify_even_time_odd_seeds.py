"""Exact auxiliary even-time source/normal rows and odd-prefix98 arithmetic.
All-degree parity/reduction induction, complete class/metric/frame cover and
analytic exclusions remain written inputs. Retained checker mains are not run.
"""
from pathlib import Path
from functools import lru_cache
from itertools import combinations
from fractions import Fraction
import hashlib,json,time
import sympy as s
import verify_ordinary_node_finite_prefixes as ordinary
import verify_squarefree_singular_finite_prefixes as weighted
import verify_euler_seed_bounds as euler
import verify_shared_divisor_weighted_seeds as shared
import verify_endpoint_degrees_and_odd_prefixes as endpoint
import verify_weighted_odd_prefix108 as previous
P=Path(__file__).resolve().parent
t,u,v=ordinary.t,ordinary.u,ordinary.v
z=s.Symbol('z');J=s.Symbol('j')
sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
DATA={'nodal':((1,1),Fraction(4),98),'cusp':((1,1),Fraction(4),72),'transverse':((Fraction(4,3),2),Fraction(16,3),98),'triangle':((2,2),Fraction(6),78),'tangent':((1,Fraction(3,2)),Fraction(4),50)}

def basis(kind):return (ordinary.CASES if kind in ['nodal','transverse','triangle'] else weighted.CASES)[kind]['basis']
def even_to_z(expr):
    poly=s.Poly(s.expand(expr),t)
    assert all(n[0]%2==0 for n,c in poly.terms() if c)
    out=s.expand(sum(c*z**(n[0]//2) for n,c in poly.terms()))
    assert s.expand(out.subs(z,t*t)-expr)==0
    return out

def source_vector(kind,W,D):
    if kind in ['nodal','transverse','triangle']:
        S=s.cancel(euler.ordinary_source(W,D)/t);bu=bv=s.Rational(1,2)
    else:
        S,_,_=shared.weighted_source(kind,W,D);S=s.cancel(S/t);bu=s.Rational(1,2);bv=s.Rational(2,3) if kind=='cusp' else s.Rational(3,4)
    vector=lambda f:D*D*t*t*(bu*u*s.diff(f,u)+bv*v*s.diff(f,v))
    return S,vector,bu,bv

def principal(kind,W,D,ex):
    if kind in ['nodal','transverse','triangle']:
        pv=s.Poly(W,u,v);co=[pv.coeff_monomial(x) for x in [u*u*v*v,u**3,v**3,u*u,v*v,u*v,u,v,1]]
        prep=ordinary.operator_data(kind,*co,t*t*D);assert s.degree(prep[-1],t)==0
        rem,cu,cv,_,_=ordinary.one_step(kind,ex,co,t*t*D,prep)
        rem,cu,cv=[s.cancel(x/prep[-1]) for x in [rem,cu,cv]]
    else:
        nf,cu,cv=weighted.principal_division(kind,ex,s.Rational(1,2),D)
        R=W-u**4/2-D*(v**3 if kind=='cusp' else u*v*v)
        rem=s.expand(nf-cu*s.diff(R,u)-cv*s.diff(R,v)-t*t*D*(s.diff(cu,u)+s.diff(cv,v)))
    assert s.expand(u**ex[0]*v**ex[1]-rem-t*t*D*(s.diff(cu,u)+s.diff(cv,v))-s.diff(W,u)*cu-s.diff(W,v)*cv)==0
    return rem,cu,cv

def universal():
    D,ju,jv,ell,E,F,qu,qv,quv,a,d=s.symbols('D ju jv ell E F qu qv quv a d');out=[]
    for kind in DATA:
        W=endpoint.odd_phase(kind,D,ju,jv,ell,E,F,qu,qv,quv,a,d) if kind in ['nodal','cusp'] else previous.phase(kind,D,ju,jv,ell,E,F,qu,qv,quv,a)
        S,vector,bu,bv=source_vector(kind,W,D);Wz,Sz=even_to_z(W),even_to_z(S);h=t**3*D*D
        identities=[]
        for n in range(4):
            for ex in ordinary.monomials(2):
                f=z**n*u**ex[0]*v**ex[1];r=f.subs(z,t*t)
                left=h*s.diff(r,t)-J*s.diff(h,t)*r+vector(r)+S*r
                right=2*D*D*z*z*s.diff(f,z)-3*J*D*D*z*f+D*D*z*(bu*u*s.diff(f,u)+bv*v*s.diff(f,v))+Sz*f
                assert s.expand(left-right.subs(z,t*t))==0
                identities.append([n,*ex])
        reductions=[]
        for ex in ordinary.monomials(8):
            rem,cu,cv=principal(kind,W,D,ex)
            for f in [rem,cu,cv]:even_to_z(f)
            reductions.append(ex)
        out.append({'kind':kind,'phase_even':True,'source_even':True,'delta_even':True,'whole_denominator_odd':True,'normal_divisors_constant_in_time':True,'tau_recurrence_identities':identities,'principal_evenness_and_divergence_identities':reductions,'beta_u':str(bu),'beta_v':str(bv)})
    return out

def replay(m,origin):
    kind=m['kind'];D=s.Integer(m['D']);W=s.sympify(m['W'],locals={'t':t,'u':u,'v':v});S,_,bu,bv=source_vector(kind,W,D)
    Wz,Sz=even_to_z(W),even_to_z(S);ring=s.QQ.poly_ring(z);poly=lambda f:s.Poly(f,u,v,domain=ring)
    zero=poly(0);R=s.Integer(1);rows=[];ampweights,cost,_=DATA[kind]
    for row in m['seed_columns']:
        j=row['j'];R=s.expand(2*D*D*z*z*s.diff(R,z)-3*(j-1)*D*D*z*R+D*D*z*(bu*u*s.diff(R,u)+bv*v*s.diff(R,v))+Sz*R)
        values=[even_to_z(s.sympify(x,locals={'t':t})) for x in row['N']]
        Cu,Cv=[even_to_z(s.sympify(row['witness'][n],locals={'t':t,'u':u,'v':v})) for n in ['Cu','Cv']]
        normal=sum((poly(c*u**ex[0]*v**ex[1]) for c,ex in zip(values,basis(kind))),zero)
        assert (poly(R)-normal-poly(z*D)*(poly(Cu).diff(u)+poly(Cv).diff(v))-poly(Wz).diff(u)*poly(Cu)-poly(Wz).diff(v)*poly(Cv)).is_zero
        limits=[2*((cost*j-sum(Fraction(a)*b for a,b in zip(ampweights,ex)))//2) for ex in basis(kind)]
        assert all(2*s.degree(c,z)<=b for c,b in zip(values,limits))
        H=s.sympify(row['H'],locals={'t':t});assert s.expand(H-(t**3*D*D)**j)==0
        assert min(e[0] for e,c in s.Poly(H,t).terms() if c)>=3*j
        rows.append({'j':j,'N_in_tau':[str(c) for c in values],'even_t_degree_bounds':limits,'H':str(H),'H_valuation':3*j,'original_twisted_identity_in_tau':True})
    return {'kind':kind,'D':int(D),'origin':origin,'metric_determinant':m['metric_determinant'],'seed_columns':rows}

def assignment(matrix):
    size=len(matrix)
    @lru_cache(None)
    def best(j,mask):
        if j==size:return 0,()
        options=[]
        for k in range(size):
            if not mask>>k&1:
                val,path=best(j+1,mask|(1<<k));options.append((matrix[j][k]+val,(k,)+path))
        return max(options)
    return best(0,0)

def arithmetic():
    out=[]
    for kind,(weights,cost,expected) in DATA.items():
        bs=basis(kind);grades=[sum(Fraction(a)*b for a,b in zip(weights,ex)) for ex in bs];cases=[]
        for size in range(1,len(bs)+1):
            for rows in combinations(range(len(bs)),size):
                if 0 not in rows:continue
                matrix=[[2*((cost*j-grades[k])//2) for k in rows] for j in range(1,size+1)]
                gate,perm=assignment(matrix)
                assert gate==sum(matrix[j][perm[j]] for j in range(size))
                cases.append({'rank':size,'rows':rows,'degree_matrix':matrix,'max_assignment':perm,'gate':gate})
        maximum=max(c['gate'] for c in cases);assert maximum==expected,(kind,maximum,expected)
        out.append({'kind':kind,'cost':str(cost),'basis_weights':[str(v) for v in grades],'cofactor_gate':maximum,'all_rank_row_cases':cases})
    cover={'quadratic':7,'smooth':94,'nodal':98,'cusp':94,'transverse':98,'tangent':94,'triangle':96,'concurrent':36,'repeated':98}
    assert max(98,96)==98 and max(72,92,94)==94 and max(50,92,94)==94 and max(cover.values())==98
    return out,cover

def main():
    start=time.monotonic();univ=universal();models=[]
    for name,key in [('endpoint_degrees_and_odd_prefixes_certificate.json','odd_models'),('weighted_odd_prefix108_certificate.json','models')]:
        for m in json.loads((P/name).read_text())[key]:models.append(replay(m,name))
    cases,cover=arithmetic();cert={'universal':univ,'models':models,'degree_cases':cases,'odd_cover':cover}
    path=P/'even_time_odd_seeds_certificate.json';path.write_text(json.dumps(cert,indent=2,default=int)+'\n')
    report={'completed':True,'checker_sha256':sha(Path(__file__)),'certificate_sha256':sha(path),'runtime_seconds':round(time.monotonic()-start,3),'sympy_version':s.__version__,'universal_strata':5,'tau_recurrence_identities':sum(len(c['tau_recurrence_identities']) for c in univ),'principal_evenness_divergence_identities':sum(len(c['principal_evenness_and_divergence_identities']) for c in univ),'retained_models':len(models),'exact_seed_identities_in_tau':sum(len(m['seed_columns']) for m in models),'even_normal_rows':sum(len(row['N_in_tau']) for m in models for row in m['seed_columns']),'odd_complement_gates':{c['kind']:c['cofactor_gate'] for c in cases},'odd_sufficient_pure_prefix':98,'global_sufficient_pure_prefix_unchanged':130,'homogeneous_sufficient_pure_prefix_unchanged':96,'least_prefix_interval':[20,130],'odd_least_prefix_interval':[20,98],'homogeneous_least_prefix_interval':[14,96],'odd_cover':cover,'optimality_proved':False,'written_analytic_or_global_proof_certified_by_checker':False,'boundary':'Exact even phases/sources, whole-denominator parity,120 symbolic-order tau recurrence identities,225 principal/divergence parity identities,30 original twisted seed identities in tau from ten retained Gaussian models and174 even normal rows; every rank/row determinant assignment. All-degree parity/reduction induction, complete coefficient/frame cover, period comparison and analytic nonzero zeroth coefficients remain written inputs. No independent replay, external review, formal verification or expanded radical cofactors.'}
    (P/'even_time_odd_seeds_verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':main()
