#!/usr/bin/env sage-python
"""One genus-one cover through the permitted first seed, with calibration explicit.

The generic pencil is frozen and checkpointed first. Two fixed members only:
z=0 and the unique member through the historical first point. No later point,
point search, new original-fibre address, or pilot input. 25-second cap.
"""
import hashlib,json,signal
from pathlib import Path
from sage.all import QQ,ZZ,GF,PolynomialRing,EllipticCurve,gcd,matrix,vector
ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
DIR=ART/'det1092_norm8_seed_cover_v2'
PRIMES=[17,47,53,61,67,71,79,83,89,101,107,113,127,137,149,179,191,197]
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def retain(p,d):
    payload=json.dumps(d,indent=2,sort_keys=True)+'\n'
    if p.exists():assert p.read_text()==payload
    else:p.write_text(payload)
def rec(f):return {'numerator':list(map(str,f.numerator().list())),
                  'denominator':list(map(str,f.denominator().list()))}
def jacobian_point(f,t0,w0):
    ee,dd,cc,bb,aa=f.list()
    I=12*aa*ee-3*bb*dd+cc*cc
    J=72*aa*cc*ee+9*bb*cc*dd-27*aa*dd*dd-27*bb*bb*ee-2*cc**3
    E=EllipticCurve(QQ,[-27*I,-27*J])
    R=f.parent();t=R.gen()
    g=R([dd*dd/16-cc*ee/6,cc*dd/12-bb*ee/2,
        cc*cc/12-bb*dd/8-aa*ee,bb*cc/12-aa*dd/2,bb*bb/16-aa*cc/6])
    # Homogeneous quartic Jacobian covariants evaluated at(t0,1).
    hv=(f.derivative()(t0)*(4*g(t0)-t0*g.derivative()(t0))
        -(4*f(t0)-t0*f.derivative()(t0))*g.derivative()(t0))/8
    P=E([36*g(t0)/w0**2,108*hv/w0**3])
    exposures=[];bound=ZZ(0);good=[]
    for p in PRIMES:
        if any(v.denominator()%p==0 for v in [*E.a_invariants(),*P[:2]]):
            exposures.append({'prime':p,'status':'SKIP_DENOMINATOR'});continue
        if E.discriminant().valuation(p)!=0:
            exposures.append({'prime':p,'status':'SKIP_BAD_REDUCTION'});continue
        Ep=EllipticCurve(GF(p),E.a_invariants());order=ZZ(Ep.cardinality())
        bound=gcd(bound,order);good.append((p,Ep))
        exposures.append({'prime':p,'status':'GOOD','order':int(order)})
    witness=None
    for p,Ep in good:
        multiple=bound*Ep(list(P[:2]))
        if not multiple.is_zero():
            witness={'prime':p,'bound':int(bound),'multiple':list(map(int,multiple))};break
    assert witness is not None,'FROZEN_PRIMES_DO_NOT_CERTIFY_NONTORSION'
    return {'Jacobian_a_invariants':list(map(str,E.a_invariants())),
       'rational_nontorsion_point':list(map(str,P[:2])),
       'finite_exposures':exposures,'nontorsion_witness':witness,
       'source':'Quartic covariant image of inherited section14 intersection, not the unlock point.'}
def main():
    paths=[DIR/'generic.json',ART/'curve302_recovered_mw17_parent_v1.json',
       ART/'det1092_first_centre_rr_net_replay_v1.json',
       ART/'det1092_rr_generic_point_controls_v2/protocol.json',Path(__file__)]
    protocol={'classification':'explicit retrospective first-seed calibration after generic construction',
       'rule':'Keep the completed generic norm8 orbit20124 pencil. Compute the unique parameter through the historical first302 witness. Compare precisely z=0 and that calibrated member, keeping the old prime pool and nine original-fibre addresses. Do not replace a failed member or import any later point.',
       'limits':{'wall_seconds':25,'generic_pencils':1,'fixed_members':2,
          'historical_first_point_inputs':1,'later_point_inputs':0,'point_searches':0,
          'old_original_addresses':9,'new_original_addresses':0,'V3_inputs':0,
          'pilot_changes':0,'class_groups':0,'Selmer_dimensions':0},
       'proof_primes':PRIMES,'inputs':{str(p.relative_to(ROOT)):sha(p) for p in paths}}
    retain(DIR/'calibration-protocol.json',protocol)
    generic,parent,first,roster=map(read,paths[:4])
    for name,digest in generic['inputs'].items():assert sha(ROOT/name)==digest
    T=PolynomialRing(QQ,'t');K=T.fraction_field();Z=PolynomialRing(QQ,'z')
    def dec(d):return K(T(d['numerator']))/T(d['denominator'])
    ai=list(map(dec,parent['a_invariants']));E=EllipticCurve(K,ai)
    G=matrix(QQ,parent['generic_height_gram']);word=vector(ZZ,generic['selection']['word'])
    basis=[E(list(map(dec,P))) for P in parent['basis_weierstrass_coordinates']]
    C=sum((n*P for n,P in zip(word,basis)),E(0));assert word*G*word==8
    cx=C[0]+E.b2()/12;cy=C[1]+(E.a1()*C[0]+E.a3())/2
    h=T(generic['pole_h']);shift=T(generic['shift'])
    nx=T(generic['nx']);ny=T(generic['ny'])
    assert nx==cx*h*h and ny==cy*h**3
    coefficients=list(map(Z,generic['quartic_t_coefficients_in_z']))
    E0=EllipticCurve(QQ,[f(0) for f in ai])
    # Only this previously certified first point is read, never the combined
    # historical/V3 construction audit mentioned in its provenance.
    P=E0(list(map(QQ,first['reconstructed_point_literal302'])))
    X=P[0]+E.b2()(0)/12;Y=P[1]+(E.a1()(0)*P[0]+E.a3()(0))/2
    slope=(Y+cy(0))/(X-cx(0));zstar=(h(0)*slope+shift(0))/h(0)**2
    Wstar=(2*X+cx(0)-slope*slope)/h(0)
    retain(DIR/'calibration.json',{'classification':'retrospective derivation; not blind selection',
       'z_star':str(zstar),'t_star':'0','W_star':str(Wstar),
       'historical_point':list(map(str,P[:2])),
       'inputs':{str(p.relative_to(ROOT)):sha(p) for p in paths[:3]}})
    covers=[]
    for index,zvalue in enumerate([QQ(0),zstar]):
        f=T([v(zvalue) for v in coefficients])
        assert f.degree()==4 and f.gcd(f.derivative()).degree()==0
        assert f.gcd(T(E.discriminant())).degree()==0
        M=h*h*zvalue-shift;m=K(M)/h
        x0=(-cx+m*m)/2-E.b2()/12;x1=K(h)/2
        y0=m*(x0+E.b2()/12-cx)-cy-(E.a1()*x0+E.a3())/2
        y1=(m-E.a1()/2)*x1
        line=[K(h)*E.a3()/2-M*E.b2()/12+(M*nx+ny)/(h*h),
              K(h)*E.a1()/2-M,K(h)]
        assert all(v.denominator().degree()==0 for v in line)
        line=list(map(T,line));assert gcd(line).degree()==0
        assert all(v.degree()<=bound for v,bound in zip(line,[8,4,2])) and line[2].degree()==2
        section=generic['sections'][0]
        def at(d):return Z(d['numerator'])(zvalue)/Z(d['denominator'])(zvalue)
        t0,W0=at(section['t_of_z']),at(section['W_of_z'])
        assert W0 and W0*W0==f(t0)
        proof=jacobian_point(f,t0,W0)
        cross=G*word
        H=(2*G).augment(matrix(QQ,17,1,list(cross))).stack(matrix(QQ,1,18,list(cross)+[8]))
        assert H.is_positive_definite() and H.det()==4*2**17*1092
        outcomes=[]
        for case in roster['cases']:
            tau=QQ(case['parameter']);value=f(tau)
            split=bool(value>=0 and value.numerator().is_square() and value.denominator().is_square())
            row={'label':case['label'],'parameter':case['parameter'],'value':str(value),'split':split}
            if value>=0:
                row['square_root_bounds']=list(map(str,[value.numerator().isqrt(),value.denominator().isqrt()]))
            if split:
                root=value.sqrt();maps=[x0(tau)+x1(tau)*root,y0(tau)+y1(tau)*root]
                row['literal_point']=list(map(str,maps))
                Es=EllipticCurve(QQ,[v(tau) for v in ai]);Ps=Es(maps)
                if index==1 and tau==0:
                    Cs=Es([C[0](0),C[1](0)])
                    assert Ps==P or Cs-Ps==P
                    row['first_point_or_centre_complement']=True
            outcomes.append(row)
        if index==1:assert Wstar*Wstar==f(0)
        result={'classification':'generic member' if index==0 else 'retrospectively calibrated member',
          'index':index,'z':str(zvalue),'quartic_coefficients':list(map(str,f.list())),
          'maps':{k:rec(v) for k,v in zip(['x0','x1','y0','y1'],[x0,x1,y0,y1])},
          'RR_line':[[str(c) for c in v.list()] for v in line],
          'inherited_basepoint':[str(t0),str(W0)],'positive_rank_base':proof,
          'pullback_height_gram':[[str(v) for v in row] for row in H.rows()],
          'height_Schur_complement':'4','pullback_rank_lower_bound':18,
          'anti_invariant_height':'16','control_outcomes':outcomes}
        retain(DIR/f'cover-{index:02d}.json',result);covers.append(result)
        print('CHECKPOINT_GENUS_ONE_COVER',index,'positive rank base',
              'splits',sum(r['split'] for r in outcomes),flush=True)
    retain(DIR/'covers.json',{'status':'CANDIDATE_POSITIVE_RANK_GENUS_ONE_COVER_THROUGH_FIRST_UNLOCK',
       'classification':'new constructive deduction; independent replay required',
       'covers':covers,'limits':protocol['limits'],
       'boundary':'The first witness calibrates z_star. This is a lower-genus carrier through that point, not a target-blind selection of z_star or a new302 seed. No later point or V3 input.',
       'inputs':{str(p.relative_to(ROOT)):sha(p) for p in paths+[DIR/'calibration-protocol.json',DIR/'calibration.json']}})

if __name__=='__main__':
    signal.alarm(25);main()
