#!/usr/bin/env python3
"""Separate equation-only 302 parent inputs from a retrospective core bridge."""
from pathlib import Path
import retrospective as r

SOURCE=r.OUT/'curve302_recovered_mw17_parent_v1.json'
PROOF=r.OUT/'curve302_recovered_mw17_parent_proof_v1.json'
OLD=r.OUT/'curve302_parent_blocks_inputs_v1.json'
INPUT=r.OUT/'rank_jump_curve302_new_parent_inputs_v1.json'
OUTPUT=r.OUT/'rank_jump_curve302_new_parent_intake_v1.json'

def compute():
    from sage.all import QQ,ZZ,PolynomialRing,matrix
    source=r.read(SOURCE);proof=r.read(PROOF)
    for name,sha in proof['input_sha256'].items():assert r.digest((r.ROOT/name).read_bytes())==sha
    assert proof['status']=='PASS_FULL_ARITHMETIC_MW17_PARENT'
    assert proof['generic_arithmetic_MW_rank']==17 and proof['full_generic_saturation']
    assert proof['specialization_parameter']=='0'
    R=PolynomialRing(QQ,'t');F=R.fraction_field()
    def decode(x):return F(R(x['numerator'])/R(x['denominator']))
    aa=list(map(decode,source['a_invariants']));assert aa[:3]==[1,1,1]
    model=[a(0) for a in aa];points=[]
    for pair in source['basis_weierstrass_coordinates']:
        x,y=map(decode,pair);x0,y0=x(0),y(0)
        assert y0*y0+x0*y0+y0==x0**3+x0*x0+model[3]*x0+model[4]
        points.append([str(x0),str(y0)])
    assert len(points)==17
    b2=aa[0]**2+4*aa[1];b4=aa[0]*aa[2]+2*aa[3];b6=aa[2]**2+4*aa[4]
    b8=aa[0]**2*aa[4]+4*aa[1]*aa[4]-aa[0]*aa[2]*aa[3]+aa[1]*aa[2]**2-aa[3]**2
    delta=-b2*b2*b8-8*b4**3-27*b6**2+9*b2*b4*b6
    assert delta.denominator().degree()==0 and delta.numerator().degree()==24 and delta(0)!=0
    assert delta.numerator().gcd(delta.numerator().derivative()).degree()==0
    # Arithmetic constructor projection: deliberately omit embeddings, height
    # reconstruction, the public D basis and all exceptional-point coordinates.
    projection={'schema':'rank-jump.curve302-new-parent-inputs.v1',
        'a_invariants':source['a_invariants'],'generic_basis':source['basis_weierstrass_coordinates'],
        'specialization_parameter':'0','specialized_model':list(map(str,model)),
        'specialized_generic_points':points,
        'provenance_boundary':'Retrospectively reconstructed parent. Generic sections are authorized inputs; no fourteen-direction quotient data is included. This projection cannot erase the parent construction\'s retrospective provenance.'}
    r.write_new(INPUT,projection)
    # This separate retrospective calculation never enters the input above.
    C=matrix(ZZ,r.read(OLD)['parents'][0]['columns']);B=matrix(ZZ,source['basis_embedding_in_public_D'])
    assert C.dimensions()==B.dimensions()==(31,17)
    U=C.change_ring(QQ).solve_right(B);assert C*U==B and all(c.denominator()==1 for c in U.list())
    U=matrix(ZZ,U);assert abs(U.det())==1
    smith=B.smith_form()[0];assert [smith[i,i] for i in range(17)]==[1]*17
    out={'schema':'rank-jump.curve302-new-parent-intake.v1','status':'PASS',
        'source_proof_bindings_verified':len(proof['input_sha256']),
        'generic_rank':17,'generic_basis_saturated':True,'height_determinant':source['height_determinant'],
        'specialization_parameter':'0','specialization_is_smooth':True,
        'discriminant_degree':24,'geometric_singular_fibres':'24 I1',
        'specialized_generic_points_verified':17,
        'old_core_to_new_basis_matrix':[list(map(int,row)) for row in U.rows()],
        'basis_change_determinant':int(U.det()),'old_core_equals_new_generic_image':True,
        'displayed_quotient':'Z^14','full_fibre_rank':proof['E302_exact_rank'],
        'construction_endpoint':'No additional strict class constructed by this intake',
        'interpretation':'The old core-relative retrospective arithmetic now has a certified generic parent. Its existing point-derived classes remain evaluation data, not arithmetic constructor inputs.',
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in [Path(__file__),SOURCE,PROOF,OLD,INPUT]}}
    r.write_new(OUTPUT,out);print('PASS: new generic image equals prior primitive core; separate equation-only projection',flush=True)

if __name__=='__main__':compute()
