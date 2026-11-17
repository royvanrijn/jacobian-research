"""Finite exact controls for the universal written first-marker/real-region proof."""
from pathlib import Path
import copy,hashlib,json,sys
from fractions import Fraction
import sympy as s
import gaussian_real_cutoff_producer as G
p=Path(__file__).resolve().parent
sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()

def main():
    examples=[];files={};root=p/'gaussian_real_cutoff_examples'
    for name,rho,decision in [('isotropic_linear','1','DECIDED'),('real_linear','0','UNSAFE'),('isotropic_quadratic','2','DECIDED')]:
        d=json.loads((root/name/'input.json').read_text());v=json.loads((root/name/'result.json').read_text())
        assert v['status']==decision and G.verify(d,v['certificate'])
        law=G.order_law(d,v['certificate']);assert law==v['order_law'] and law['marker_slope']==rho
        if decision=='DECIDED':assert v['region']['vertices']==[[law['asymptotic_coefficient'],'0']]
        examples.append(dict(name=name,status=decision,order_law=law,certificate=v['certificate']))
    # Direct original heat orders,including genuine nonmonotonicity.
    orders=[]
    for name,expected in [('isotropic_linear',lambda m:m),('real_linear',lambda m:m%2),('isotropic_quadratic',lambda m:2*m)]:
        d=G.field_family(json.loads((root/name/'input.json').read_text()))
        actual=[G.heat_order(d,m) for m in range(13)]
        assert actual==[expected(m) for m in range(13)]
        orders.append(dict(name=name,original_orders=actual,monotonicity_claimed=name!='real_linear'))
    # Exact all-index nilpotent-jet family;the proof contains its binomial identity.
    jet=dict(dimension=3,parameter_count=1,phase_terms=[[[2,0,0,0],'1'],[[0,0,2,1],'1']],
             parameter_ideal=[[[[2],'1']]],precision_matrix=[[0,1,0],[1,0,0],[0,0,1]])
    assert [G.heat_order(jet,m) for m in range(1,7)]==[2*m-2 for m in range(1,7)]
    visits=[0]
    def jet_oracle(C,E):
        if E>=1:return dict(status='DECIDED',zero=True,certificate='written exact d_m=2m-2 and C>=1/2')
        visits[0]+=1;m=4 if visits[0]==1 else 1
        marker=[0,2*m-2,0]
        value=G.O.original_moment(jet,marker,m)
        cs=s.symbols('c:1');gb=s.groebner([cs[0]**2],*cs,domain=s.QQ)
        assert gb.reduce(value)[1]==m*s.factorial(2*m-2)*cs[0]
        return dict(status='DECIDED',zero=False,witness=dict(marker=marker,power=m,polynomial=G.O.encode_polynomial(value,cs)))
    refined=G.refine_region(Fraction(1,2),lambda m:G.heat_order(jet,m),jet_oracle)
    assert refined['vertices']==[['1/2','1']]
    assert [v['power'] for v in refined['supporting_original_points']]==[4,1]
    nil=[]
    for kind,phase,H,L,points,expected_vertices in [
        ('pure nilpotent',[[[1,0,1],'1']],[[0,1],[1,0]],2,[(1,1)],[['0','1'],['1','0']]),
        ('nonpure nilpotent',[[[1,1],'1']],[[1]],3,[(1,1),(2,0)],[['0','2']])]:
        data=dict(dimension=len(H),parameter_count=1,phase_terms=phase,precision_matrix=H,parameter_ideal=[[[[L],'1']]])
        actual=[G.heat_order(data,m) for m in range(1,L+1)]
        assert actual==[q for m,q in points]+[None]
        vertices=[[str(C),str(E)] for C,E in G.vertices(0,points)]
        assert vertices==expected_vertices
        nil.append(dict(name=kind,orders=actual,vertices=vertices,nilpotence_index=L))
    oracle=[]
    data=dict(dimension=1,parameter_count=1,phase_terms=[[[1,1],'1']],parameter_ideal=[[[[3],'1']]])
    for E,zero in [('2',True),('1',False)]:
        result=G.rational_zero(data,'1/2',E,root/'unused_exact_nilpotent',timeout=10)
        assert result['status']=='DECIDED' and result['zero']==zero
        oracle.append(result)
    # The completed ordinary-ideal adapter attempt is UNKNOWN;never rerun it.
    capped=root/'rational_ordinary_ideal/result.json';unknown=json.loads(capped.read_text())
    assert unknown['status']=='UNKNOWN' and unknown['underlying_result']['status']=='UNKNOWN'
    rejected=[];original=examples[0]['certificate'];data=json.loads((root/'isotropic_linear/input.json').read_text())
    def reject(name,fn):
        try:fn()
        except (ValueError,TypeError,KeyError):rejected.append(name)
        else:raise AssertionError(name)
    def tamper(field):
        c=copy.deepcopy(original)
        if field=='scalar':c['Euler']['scalar_terms'].append([[0,0],'1'])
        elif field=='raising':c['Euler']['operator_terms'].append([[0]*5,'1'])
        elif field=='start':c['specialization_start']=1
        elif field=='witness':c['ordinary']['left_witnesses'][0].append([[0]*6,'1'])
        G.verify(data,c)
    for field in ['scalar','raising','start','witness']:reject('damaged '+field,lambda field=field:tamper(field))
    altered=copy.deepcopy(data);altered['phase_terms'][0][1]='2'
    reject('changed original phase',lambda:G.verify(altered,original))
    reject('negative rational cutoff',lambda:G.rational_zero(data,'-1','0',root/'unused_reject'))
    reject('invalid origin source precision',lambda:G.source(dict(dimension=1,phase_terms=[[[1],'1']],precision_matrix=[[0]])))
    for f in root.rglob('*'):
        if f.is_file():files[str(f.relative_to(p))]=sha(f)
    packet=dict(generated_QQ_examples=examples,original_heat_order_controls=orders,nonreduced_jet_region=refined,
                nilpotent_regions=nil,rational_zero_oracles=oracle,rejections=rejected,
                capped_rational_ordinary_ideal=dict(status='UNKNOWN',file=str(capped.relative_to(p)),sha256=sha(capped),restarted=False),
                universal_claims_use_written_proofs=True,general_presented_ring_jet_certificate_backend_implemented=False,
                independent_replay=False,formal_verification=False,external_review=False)
    out=p/'gaussian_real_cutoff_certificate.json';out.write_text(json.dumps(packet,indent=2)+'\n')
    runtime=dict(Python=sys.version,Sympy=s.__version__,sealed_imports_sha256={str(Path(mod.__file__).relative_to(p)):sha(Path(mod.__file__)) for mod in [G.W,G.S,G.O,G.O.F,G.O.G]})
    run=p/'gaussian_real_cutoff_runtime.json';run.write_text(json.dumps(runtime,indent=2)+'\n')
    receipt=dict(completed=True,canonical_proof_sha256=sha(p/'GAUSSIAN_FIRST_MARKER_SLOPES_AND_SHARP_REAL_CUTOFF_POLYHEDRA.md'),
        producer_sha256=sha(p/'gaussian_real_cutoff_producer.py'),checker_sha256=sha(Path(__file__)),certificate_sha256=sha(out),runtime_sha256=sha(run),
        generated_QQ_examples=3,exact_graph_identities=6,original_heat_orders=39,nonreduced_jet_orders=6,
        nonreduced_region_refinements=2,nilpotent_regions=2,rational_zero_oracles=2,rejections=len(rejected),retained_example_hashes=files,
        ordinary_rational_ideal_run='terminal UNKNOWN at10second cap;not restarted',historical_cusp_and_cone_restarted=False,
        independent_replay=False,formal_verification=False,external_review=False)
    (p/'gaussian_real_cutoff_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:v for k,v in receipt.items() if 'sha256' not in k and k!='retained_example_hashes'}))

if __name__=='__main__':main()
