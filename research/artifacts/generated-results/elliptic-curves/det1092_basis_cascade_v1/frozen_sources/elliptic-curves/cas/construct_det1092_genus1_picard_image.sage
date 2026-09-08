#!/usr/bin/env sage-python
"""Full Picard-image restriction inputs for two already frozen genus1 curves.

One old parity coset, <=200000 exact nodes; all19 NS generators accounted
for. No marked seed, point search, new curve selector or V3 input.25s cap.
"""
import hashlib,json,signal,sys
from pathlib import Path
from sage.all import QQ,ZZ,PolynomialRing,EllipticCurve,matrix,vector
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
CAS=ROOT/'elliptic-curves/cas';SOURCE=ART/'det1092_norm8_seed_cover_v2'
DIR=ART/'det1092_genus1_picard_image_v1'
sys.path.insert(0,str(CAS));from visibility_lattice_v2 import ExactParity
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def retain(p,d):
    payload=json.dumps(d,indent=2,sort_keys=True)+'\n'
    if p.exists():assert p.read_text()==payload
    else:p.write_text(payload)
def rec(f):return {'numerator':list(map(str,f.numerator().list())),
                  'denominator':list(map(str,f.denominator().list()))}
def main():
    paths=[ART/'curve302_recovered_mw17_parent_v1.json',SOURCE/'generic.json',
       SOURCE/'equation-only-cover-00.json',SOURCE/'equation-only-cover-01.json',
       ART/'det1092_rr_generic_point_controls_v2/protocol.json',CAS/'visibility_lattice_v2.py',Path(__file__)]
    parent,generic=map(read,paths[:2]);roster=read(paths[4])
    protocol={'classification':'generic Picard-image construction on two frozen members',
      'rule':'Use only the completed norm8 pencil and the two fixed equation payloads. Account for O,F and all17 original sections, normalized at the inherited P14 intersection. Freeze the old64-prime pool before any marked first-point evaluation. Do not substitute a smaller inherited subgroup.',
      'limits':{'wall_seconds':25,'fixed_members':2,'NS_generators':19,
        'parity_cosets':1,'CVP_node_limit':200000,'primes':64,'new_addresses':0,
        'marked_point_inputs':0,'point_searches':0,'V3_inputs':0,'pilot_changes':0,
        'Selmer_dimensions':0,'class_groups':0},
      'primes':roster['primes'],'inputs':{str(p.relative_to(ROOT)):sha(p) for p in paths}}
    retain(DIR/'protocol.json',protocol)
    G=matrix(ZZ,parent['generic_height_gram']);w=vector(ZZ,[0]*14+[1,-1,0])
    U=G.LLL_gram();wr=U.inverse()*w;Gr=U.transpose()*G*U
    cvp=ExactParity(Gr.rows()).solve(wr,wr,node_limit=200000);assert cvp['norm']==8
    vertical=[];kernel=[]
    D=vector(QQ,[2,4,*w]);zero=vector(QQ,[1,QQ(G[14,14])/2,*[int(i==14) for i in range(17)]])
    kernel=[D,zero]
    for r in cvp['minima']:
        v=U*vector(ZZ,r);x=(w-v)/2;assert all(c in ZZ for c in x)
        assert x*G*x==w*G*x
        section=vector(QQ,[1,(x*G*x)/2,*x]);kernel.append(section)
        vertical.append({'word':list(map(int,x)),'height':int(x*G*x),
                         'NS_coordinates':list(map(str,section))})
    kr=matrix(QQ,kernel).rank();upper=19-kr
    frame={'classification':'exact kernel lower bound, hence full-image rank upper bound',
      'generic_Gram':[[int(v) for v in row] for row in G.rows()],
      'D':list(map(str,D)),'new_zero_section':list(map(str,zero)),
      'vertical_old_sections':vertical,'kernel_rank':int(kr),'Picard_image_rank_upper_bound':int(upper),
      'CVP_nodes':cvp['nodes'],'minimum_vectors':len(cvp['minima']),
      'argument':'O and every listed vertical old section have intersection zero with each smooth member, hence restrict trivially. The normalized P14 restriction is zero. D is a fibre and restricts trivially. These classes span a kernel of the stated rank in the full rank19 geometric NS group.',
      'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [paths[0],paths[1],DIR/'protocol.json']}}
    retain(DIR/'frame.json',frame)
    print('CHECKPOINT_PICARD_KERNEL',kr,'full image rank at most',upper,'vertical sections',len(vertical),flush=True)
    R=PolynomialRing(QQ,'t');K=R.fraction_field();t=R.gen()
    def dec(d):return K(R(d['numerator']))/R(d['denominator'])
    ai=list(map(dec,parent['a_invariants']));E=EllipticCurve(K,ai)
    basis=[E(list(map(dec,P))) for P in parent['basis_weierstrass_coordinates']]
    C=basis[14]-basis[15];cx=C[0]+E.b2()/12;cy=C[1]+(E.a1()*C[0]+E.a3())/2
    h=R(generic['pole_h']);shift=R(generic['shift'])
    section_functions=[]
    for i,P in enumerate(basis):
        xx=P[0]+E.b2()/12;yy=P[1]+(E.a1()*P[0]+E.a3())/2
        m=(yy+cy)/(xx-cx);z=(h*m+shift)/(h*h);W=(2*xx+cx-m*m)/h
        degree=ZZ(G[i,i]-(G*w)[i])
        assert max(z.numerator().degree(),z.denominator().degree())==degree
        section_functions.append({'index':i,'intersection_degree':int(degree),'z':rec(z),'W':rec(W)})
    retain(DIR/'generic-restrictions.json',{'classification':'generic-only restriction functions',
       'sections':section_functions,'inputs':{str(p.relative_to(ROOT)):sha(p) for p in paths[:2]}})
    for index,payload_path in enumerate(paths[2:4]):
        payload=read(payload_path);q=R(payload['quartic_coefficients']);zvalue=QQ(payload['cover_parameter_z'])
        source=payload['infinite_base_point_source'];t0=QQ(source['t0']);s=QQ(source['s'])
        shifted=q(t+t0);q0,q1,q2,q3,q4=shifted.list();assert q0==s*s
        cubic=[q2,q1*q3-4*s*s*q4,s*s*q3*q3+q1*q1*q4-4*s*s*q2*q4]
        Pbar=list(map(QQ,source['generator_XY']))
        J=EllipticCurve(QQ,[0,cubic[0],0,cubic[1],cubic[2]]);assert J(Pbar)
        divisors=[]
        for row in section_functions:
            i=row['index'];degree=row['intersection_degree']
            if i==14:
                assert degree==1
                divisors.append({'index':i,'status':'ZERO_NORMALIZED_BASEPOINT','degree':1});continue
            if degree==0:
                assert dec(row['z'])!=zvalue
                divisors.append({'index':i,'status':'ZERO_VERTICAL_SECTION','degree':0});continue
            z=dec(row['z']);g=R(z.numerator()-zvalue*z.denominator()).monic()
            assert g.degree()==degree and g.gcd(g.derivative())==1
            W=dec(row['W']);assert g.gcd(W.denominator())==1 and g(t0)
            Wmod=(W.numerator()*W.denominator().inverse_mod(g))%g
            assert (Wmod*Wmod-q)%g==0
            du=t-t0;Xmod=((2*s*(Wmod+s)+q1*du)*(du*du).inverse_mod(g))%g
            Ymod=(((Xmod*Xmod-4*s*s*q4)*du-q1*Xmod-2*s*s*q3)/(2*s))%g
            assert (Ymod*Ymod-Xmod**3-cubic[0]*Xmod*Xmod-cubic[1]*Xmod-cubic[2])%g==0
            divisors.append({'index':i,'status':'EXACT_RESTRICTION_DIVISOR','degree':int(degree),
              'g':list(map(str,g.list())),'Wmod':list(map(str,Wmod.list())),
              'Jacobian_Xmod':list(map(str,Xmod.list())),'Jacobian_Ymod':list(map(str,Ymod.list()))})
        result={'classification':'all inherited NS generators accounted for',
          'index':index,'z':str(zvalue),'quartic':list(map(str,q.list())),
          'basepoint':[str(t0),str(s)],'Jacobian_cubic_coefficients':list(map(str,cubic)),
          'old_O_restriction':'ZERO_D_DOT_O_IS_ZERO',
          'old_F_restriction':{'point':list(map(str,Pbar)),'identity':'F|C-2P0 = Pbar-P0'},
          'old_section_restrictions':divisors,'rank_upper_bound':int(upper),
          'marked_point_evaluated':False,
          'inputs':{str(p.relative_to(ROOT)):sha(p) for p in
                    [payload_path,DIR/'frame.json',DIR/'generic-restrictions.json',DIR/'protocol.json']}}
        retain(DIR/f'case-{index:02d}-generic.json',result)
        print('CHECKPOINT_COMPLETE_PICARD_RESTRICTIONS',index,[r['degree'] for r in divisors],flush=True)
if __name__=='__main__':
    signal.alarm(25);DIR.mkdir(parents=True,exist_ok=True);main()
