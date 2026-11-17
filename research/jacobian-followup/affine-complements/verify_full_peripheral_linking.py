"""PF2RULED101: exact full peripheral controls and universal Z[t] certificates."""
from pathlib import Path
import sys,json,hashlib,platform,math,itertools
import sympy as S
from sympy.matrices.normalforms import smith_normal_form
P=Path(__file__).resolve().parent;sys.path.insert(0,str(P))
import full_peripheral_linking as N
import composed_peripheral_cover as Old
from longitude_joining_transfer_certificate import certificate
import verify_cyclic_loop_first_arrival as K


def natural_seed(q=3,regular=False):
    A=(0,2,1,0) if q%12==3 else (0,1,2,0);B=(0,2,1,1)
    if regular:
        G=[a for a in itertools.product(range(3),repeat=4) if (a[0]*a[3]-a[1]*a[2])%3==1]
        def mul(a,b):return tuple(sum(a[2*i+h]*b[2*h+j] for h in range(2))%3 for i in range(2) for j in range(2))
        return [tuple(G.index(mul(a,b)) for b in G) for a in [A,B]]
    V=[v for v in itertools.product(range(3),repeat=2) if v!=(0,0)]
    return [tuple(V.index(tuple(sum(a[2*i+h]*v[h] for h in range(2))%3 for i in range(2))) for v in V) for a in [A,B]]


def scale(matrix,factor):return [[str(S.Rational(x)*factor) for x in row] for row in matrix]


def main():
    checks=[]
    def check(name,value):
        assert value,name;checks.append(name)
    def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
    predecessors=['polynomial_key_chain_receipt.json','proper_finite_marked_end_receipt.json',
                  'recursive_linking_endpoints_receipt.json']
    for name in predecessors:
        old=json.loads((P/name).read_text());check('predecessor PASS '+name,old['outcome']=='PASS')
        for f,h in old['source_sha256'].items():check('sealed predecessor '+name+' '+f,sha(P/f)==h)
    manifest=P/'sources/full-peripheral-linking/manifest.json';m=json.loads(manifest.read_text())
    for row in m['primary_sources']:
        f=P/row['path'];check('exact primary pin '+f.name,f.stat().st_size==row['bytes'] and sha(f)==row['sha256'])
    author=json.loads((P/'dihedral_author_comparison.json').read_text())
    pin=json.loads((P/'sources/recursive-linking-endpoints/author_code_pin.json').read_text())
    check('unmodified original author source and completed comparison retained',
          sha(P/pin['path'])==pin['sha256'] and author['outcome']=='PASS')
    # No parameter sampling: two fixed unimodular reductions over Z[t].
    universal=[certificate(e) for e in [1,3]]
    for c in universal:
        check('universal Z[t] unit certificate e='+str(c['normal_index']),
              c['basis_determinant'] in (-1,1) and c['H1_free_rank']==2 and
              c['H1_torsion']==[] and c['parameter_samples_used_as_proof'] is False)
    # Independent literal path tracker, restricted to SMALL word controls.
    def wp(w,n):return w*n if n>=0 else [-g for g in reversed(w)]*(-n)
    def literal(w,seed,rho,keep):
        d=len(rho[0]);row=[0]*(d*len(rho));s=seed
        for g in w:
            if g>0:
                row[(g-1)*d+s]+=1;s=N.inverse(rho[g-1])[s]
            else:
                h=-g-1;s=rho[h][s];row[h*d+s]-=1
        return s,[row[g*d+j] for g,j in keep]
    s4rho=[(1,0,2,3),(0,2,3,1)]
    for rho,cables in [(s4rho,[]),([(1,0,2),(1,2,0)],[]),
                       ([(1,0,2),(1,2,0),(0,1,2)],[(3,2)])]:
        data=N.peripheral_cover(2,3,cables,rho);rels,mu,lam=K.words(2,3,cables)
        check('independent literal relator rows '+str(cables)+str(len(rho[0])),
              [literal(w,s,rho,data['keep'])[1] for w in rels for s in range(len(rho[0]))]==data['relation_rows'])
        for b in data['boundaries']:
            s=b['seed'];e=b['normal_index'];a=b['shear'];f=b['longitudinal_degree']
            check('independent literal full peripheral lift '+str((len(rho[0]),cables,s)),
                  literal(wp(mu,e),s,rho,data['keep'])==(s,b['meridian']) and
                  literal(wp(mu,a)+wp(lam,f),s,rho,data['keep'])==(s,b['longitude']))
    s4=[]
    for reverse in [False,True]:
        data=N.peripheral_cover(2,3,[],s4rho,reverse);out=N.rational_core_linking(data)
        b=data['boundaries'][0];rows=data['relation_rows']
        check('nonzero shear and rational core order exact '+str(reverse),
              (b['normal_index'],b['longitudinal_degree'],b['shear'])==(4,1,2) and
              out['matrix']==[['1/2']] and out['core_orders']==[4] and
              out['rational_row_denominators']==[2] and out['cleared_relation_torsion_orders']==[2] and
              out['closed_H1']==dict(free=0,torsion=[8]) and N.reduced_h1(rows)==dict(free=1,torsion=[2]))
        direct=smith_normal_form(S.Matrix(rows+[b['meridian']]),domain=S.ZZ)
        diag=[abs(int(direct[i,i])) for i in range(min(direct.shape)) if direct[i,i]!=0]
        check('independent unreduced Smith normal filling '+str(reverse),[x for x in diag if x>1]==[8])
        fills={}
        for slope in [-2,-1,0,1,2,3]:
            row=[x+slope*y for x,y in zip(b['meridian'],b['longitude'])]
            fills[str(slope)]=N.reduced_h1(rows,[row])
        s4.append(dict(data=out,arbitrary_integral_fillings_m_plus_sL=fills))
    check('both spanning trees retain rational matrix torsion defects and all filling groups',s4[0]==s4[1])
    abstract=dict(relation_rows=[[0]],boundaries=[dict(meridian=[2],longitude=[1],normal_index=2,
                  longitudinal_degree=1,shear=1)])
    other=N.rational_core_linking(abstract)
    check('same rational matrix does not determine actual integral core order or filled group',
          other['matrix']==[['1/2']] and other['core_orders']==[2] and other['closed_H1']==dict(free=0,torsion=[2]))
    guard=dict(relation_rows=[[0,0]],boundaries=[dict(meridian=[1,0],longitude=[0,1],normal_index=1,
               longitudinal_degree=1,shear=0)])
    try:N.rational_core_linking(guard)
    except ValueError:check('non-rational-null core rejected rather than assigning a linking number',True)
    else:check('non-rational-null core rejected rather than assigning a linking number',False)
    # Dihedral overlap: matrices from the sealed PF92 method and the original
    # authors-code receipt agree OFF the diagonal. No fresh author execution claimed.
    overlap=[]
    for p in [3,5,7]:
        rho=[tuple((-i)%p for i in range(p)),tuple((i+1)%p for i in range(p))]
        d=N.peripheral_cover(2,p,[],rho);now=N.rational_core_linking(d)
        old=Old.integral_null_core_linking(Old.peripheral_cover(2,p,[],rho))
        check('full periphery reduces to old integral preferred framing degree '+str(p),
              now['matrix']==[[str(x) for x in row] for row in old['matrix']] and
              all(b['longitudinal_degree']==1 and b['shear']==0 for b in d['boundaries']))
        off=[S.Rational(now['matrix'][i][j]) for i in range(len(d['boundaries'])) for j in range(i)]
        want=[S.Integer(2)] if p==3 else list(map(S.Rational,next(x for x in author['higher_degree_controls'] if x['p']==p)['offdiagonal']))
        check('unmodified author off-diagonal comparison retained degree '+str(p),off==want)
        overlap.append(dict(degree=p,matrix=now['matrix'],author_offdiagonal=[str(x) for x in want],
                            author_diagonal_framing_identified=False,new_author_execution_claimed=False))
    seeds=[];family=[]
    for regular in [False,True]:
        rho=natural_seed(regular=regular);data=N.peripheral_cover(2,3,[],rho)
        seed=N.rational_core_linking(data);size=24 if regular else 8
        want=([['-3' if i==j else '1' for j in range(4)] for i in range(4)] if regular else [['-9','3'],['3','-1']])
        check('actual longitude joins meridian cycles seed degree '+str(size),
              seed['matrix']==want and seed['closed_H1']==dict(free=0,torsion=[]) and
              seed['longitudinal_degrees']==[2]*len(want) and seed['shears']==[0]*len(want))
        reverse=N.rational_core_linking(N.peripheral_cover(2,3,[],rho,True))
        check('both trees retain full longitude-joining seed degree '+str(size),seed==reverse)
        seeds.append(dict(degree=size,rho=[list(x) for x in rho],result=seed))
        M=data['meridian_monodromy'];L=data['longitude_monodromy'];z=tuple(M[L[i]] for i in range(size))
        for params in ([[1],[1,2,1],[1]*12,[10**30]] if not regular else [[1]]):
            cables=[(36*t+1,36*t-5) for t in params];R=math.prod(r for r,k in cables)
            result=N.rational_core_linking(N.peripheral_cover(2,3,cables,rho+[z]*len(cables)))
            check('all-length exact family control '+str((size,len(params),params[0])),
                  result['matrix']==scale(seed['matrix'],R*R) and
                  result['closed_H1']==seed['closed_H1'] and result['core_orders']==seed['core_orders'] and
                  result['longitudinal_degrees']==seed['longitudinal_degrees'])
            family.append(dict(degree=size,parameters=params,winding_product=str(R),result=result,
                               expanded_words_or_diagrams=False))
    # Added torsion is retained rather than assumed absent from every companion.
    rho=natural_seed(q=9);data=N.peripheral_cover(2,9,[],rho);before=N.rational_core_linking(data)
    M=data['meridian_monodromy'];L=data['longitude_monodromy'];z=tuple(M[L[i]] for i in range(8))
    after=N.rational_core_linking(N.peripheral_cover(2,9,[(37,31)],rho+[z]))
    check('exact transfer retains unrelated closed torsion summand',
          before['closed_H1']==after['closed_H1']==dict(free=0,torsion=[3]) and
          after['matrix']==scale(before['matrix'],37**2))
    receipt=dict(as_of='2026-11-16',outcome='PASS',check_count=len(checks),checks=checks,
        general=dict(arbitrary_longitude_monodromy=True,joint_orbit_boundary_tori=True,
                     integral_surface_balance_eL_minus_am=True,rational_linking_symmetric_Fe_equals_a=True,
                     exact_core_orders_and_minimal_denominator_torsion_defect=True,
                     arbitrary_filling_integer_rows_retained=True,all_finite_iterated_lengths_composed=True),
        universal_cable_theorem=dict(orbit_types=[[1,2,0],[3,2,0]],parameters='r=36t+1,k=36t-5,t>=1',
             exact_Zt_unimodular_certificates=universal,arbitrary_cover_degree=True,arbitrary_finite_length=True,
             F_new='r^2 F_old',closed_H1_preserved=True,core_class_multiplied_by_winding=True,
             core_order_new='q/gcd(q,r)',finite_samples_used_as_universal_proof=False),
        shear_and_torsion_control=s4[0],same_matrix_different_integral_control=other,
        dihedral_overlap=overlap,nontrivial_longitude_seeds=seeds,all_length_family_controls=family,
        added_torsion_control=dict(before=before,after=after),
        assurance=dict(written_full_peripheral_and_linking_proofs=True,unbounded_parameter_unit_certificate=True,
                       finite_controls_not_universal_proof=True,published_quartic_retirement_preserved=True,
                       original_author_comparison_reused_truthfully=True,
                       complete_filled_torsion_linking_form_computed=False,new_covers_identified_as_dihedral=False,
                       arbitrary_Laurent_Keller_boundaries_excluded=False,new_higher_degree_Keller_exclusion_proved=False,
                       arbitrary_affine_or_Keller_realization_asserted=False,formal_verification=False,
                       external_review=False,research_goal_completed=False),
        source_sha256={name:sha(P/name) for name in [Path(__file__).name,'full_peripheral_linking.py',
          'longitude_joining_transfer_certificate.py','FULL_PERIPHERAL_RATIONAL_LINKING_AND_LONGITUDE_JOINING_RECURSION.md',
          'sources/full-peripheral-linking/manifest.json','composed_peripheral_cover.py','peripheral_boundary_lattices.py',
          'dihedral_author_comparison.json','sources/recursive-linking-endpoints/author_code_pin.json',pin['path']]+predecessors+
          [r['path'] for r in m['primary_sources']]},runtime=dict(python=platform.python_version(),sympy=S.__version__))
    (P/'full_peripheral_linking_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print('PASS',len(checks),'controls; universal Z[t] cable theorem and full-peripheral rational linking; actual Keller gate open')


if __name__=='__main__':main()
