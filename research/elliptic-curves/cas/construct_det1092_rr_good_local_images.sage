#!/usr/bin/env sage-python
"""Class-group-free local 2-descent on the frozen 64 odd primes.

Use generic section restrictions only. Normalize each divisor separately;
a nonunit generator does not discard the other generators at that prime.
"""
import hashlib,json,signal
from collections import Counter
from pathlib import Path
from sage.all import QQ,ZZ,GF,PolynomialRing,matrix,vector,power_mod
from sage.version import version as sage_version
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
OUT=ART/'det1092_rr_good_local_images_v1';LOCAL=ROOT/'artifacts/local/elliptic-curves/det1092-rr-good-local-images-v1'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def retain(p,d):
    payload=json.dumps(d,indent=2,sort_keys=True)+'\n'
    if p.exists():assert p.read_text()==payload
    else:p.write_text(payload)
def build():
    generic_path=ART/'det1092_rr_full_inherited_jacobian_v1.json';generic=json.loads(generic_path.read_text())['generic_r_functions']
    proto_path=ART/'det1092_rr_generic_point_controls_v2/protocol.json';old=json.loads(proto_path.read_text());primes=old['primes']
    sources=[ART/'det1092_rr_generic_point_controls_v2'/('case-%02d.json'%i) for i in range(9)]+[
        ART/'det1092_rr_residual_jacobian_class_v1.json']
    protocol={'classification':'frozen local-descent calculation using generic divisor inputs only',
              'primes':primes,'source_only_cases':9,'retrospective_calibration_cases':1,
              'limits':{'wall_seconds':25,'prime_case_pairs':640,'generic_divisors_per_case':17,
                        'new_prime_searches':0,'point_searches':0,'global_Selmer_runs':0,'class_group_runs':0,
                        'external_computations':0,'pilot_changes':0},
              'inputs':{str(p.relative_to(ROOT)):sha(p) for p in sources+[generic_path,proto_path,Path(__file__)]}}
    OUT.mkdir(parents=True,exist_ok=True);LOCAL.mkdir(parents=True,exist_ok=True);retain(OUT/'protocol.json',protocol)
    R=PolynomialRing(QQ,'T');T=R.gen();cases=[];totals=Counter()
    for i,path in enumerate(sources):
        source=json.loads(path.read_text());curve=source['curve'] if i==9 else source
        q=R(curve['q']);scale=QQ(curve['scale']);u=QQ(curve['u'])
        x0=QQ(source['inherited_x_values'][source['base_pair_index']] if i==9 else source['base_x'])
        gs=[(R(row['numerator'])-u*R(row['denominator'])).monic() for row in generic]
        assert all(g.degree()==row['degree'] for g,row in zip(gs,generic))
        trials=[]
        for p in primes:
            F=GF(p);S=PolynomialRing(F,'X');X=S.gen();trial={'p':p}
            qv=min(a.valuation(p) for a in q if a);qm=QQ(p)**(-qv);fp=S((qm*q).list())
            adjusted_scale=scale/qm
            trial['q_multiplier']=str(qm);trial['adjusted_scale_valuation']=int(adjusted_scale.valuation(p))
            if fp.degree()!=6 or fp.gcd(fp.derivative())!=1:
                trial['status']='UNRESOLVED_NOT_SQUAREFREE_DEGREE_SIX_REDUCTION'
            elif adjusted_scale.valuation(p)%2:
                trial['status']='UNRESOLVED_ODD_SCALAR_TWIST'
            else:
                factors=sorted([g.monic() for g,e in fp.factor()],key=lambda g:(g.degree(),list(map(int,g))))
                degrees=[int(g.degree()) for g in factors];k=len(degrees);odd=any(d%2 for d in degrees)
                dimension=k-1-int(odd)
                anchor=R([1,-1/x0]) if x0 and x0.valuation(p)<0 else R([x0,-1])
                anchor_p=S(anchor.list())
                trial.update({'good_reduction_verified':True,'degrees':degrees,
                              'factors':[list(map(int,g.list())) for g in factors],
                              'local_Kummer_dimension':dimension,'anchor':list(map(str,anchor.list()))})
                if anchor_p.gcd(fp)!=1:
                    trial['status']='UNRESOLVED_BASEPOINT_MEETS_BRANCH'
                else:
                    raw=[];kept=[];omitted=[];normalizations=[]
                    for j,g in enumerate(gs):
                        gv=min(a.valuation(p) for a in g if a);mult=QQ(p)**(-gv);gp=S((mult*g).list())
                        normalizations.append(str(mult))
                        if gp.gcd(fp)!=1:
                            omitted.append(j);continue
                        bits=[]
                        for h in factors:
                            value=((-1)**g.degree()*gp*power_mod(anchor_p,-g.degree(),h))%h
                            sign=power_mod(value,(ZZ(p)**h.degree()-1)//2,h)
                            assert sign in [S(1),S(-1)];bits.append(int(sign==S(-1)))
                        assert sum(bits)%2==0
                        kept.append(j);raw.append(bits)
                    raw_matrix=matrix(GF(2),raw).transpose() if raw else matrix(GF(2),k,0)
                    parity=vector(GF(2),[d%2 for d in degrees])
                    projection=matrix(GF(2),1,k,list(parity)).right_kernel().basis_matrix()
                    block=projection*raw_matrix;rank=block.rank();assert rank<=dimension
                    trial.update({'status':'COMPLETE_LOCAL_KUMMER_IMAGE' if rank==dimension else 'INCOMPLETE_GENERIC_LOCAL_SPAN',
                                  'generic_divisor_indices':kept,'omitted_nonunit_divisors':omitted,
                                  'divisor_p_normalizations':normalizations,
                                  'raw_character_rows':[list(map(int,row)) for row in raw_matrix.rows()],
                                  'scalar_quotient_projection':[list(map(int,row)) for row in projection.rows()],
                                  'image_rows':[list(map(int,row)) for row in block.rows()],
                                  'generic_image_rank':int(rank),
                                  'basis_generic_divisors':[kept[j] for j in block.pivots()]})
            trials.append(trial);retain(LOCAL/('case-%02d-prime-%04d.json'%(i,p)),trial)
        counts=Counter(row['status'] for row in trials);totals.update(counts)
        result={'classification':'verified application; independent replay required',
                'case_index':i,'arm':'RETROSPECTIVE_CALIBRATION' if i==9 else 'GENERIC_SOURCE_ONLY_COHORT',
                'q':list(map(str,q.list())),'scale':str(scale),'u':str(u),'base_x':str(x0),
                'generic_divisors':[list(map(str,g.list())) for g in gs],
                'trials':trials,'counts':dict(counts),'source_path':str(path.relative_to(ROOT)),
                'full_global_Selmer_group':'NOT_COMPUTED'}
        retain(OUT/('case-%02d.json'%i),result);cases.append({'case_index':i,'counts':dict(counts)})
    result={'classification':'verified application and new deduction; local conditions only',
            'status':'PASS_FROZEN_GOOD_PRIME_LOCAL_IMAGE_CONSTRUCTION',
            'sage_version':sage_version,'cases':cases,'totals':dict(totals),
            'dimension_formula':'For k irreducible factors, dim J(Qp)/2J(Qp)=k-1 if all degrees are even, otherwise k-2, at the verified good odd primes.',
            'proof':'Good odd reduction identifies the local Kummer image with the unramified finite-field image. Its dimension is the Frobenius-fixed dimension of even branch subsets modulo complement. Unit squareclasses have even norm-character sum and are quotiented by scalar degree parity. The recorded generic classes attain that dimension in every row marked complete.',
            'scope':'Only the declared64 odd primes. Dyadic, real, bad-place and untested-prime conditions, and the global unramified class group, are not computed. A complete local image is not a full global Selmer group.',
            'limits':protocol['limits'],'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [OUT/'protocol.json',Path(__file__)]}}
    retain(OUT/'manifest.json',result);print(result['status'],dict(totals),flush=True)
if __name__=='__main__':signal.alarm(25);build()
