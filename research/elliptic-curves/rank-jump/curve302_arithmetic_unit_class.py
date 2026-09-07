#!/usr/bin/env python3
"""Bounded arithmetic units from repeated reduced ideals on302."""
import argparse
from pathlib import Path
import subprocess
import sys
import retrospective as r
from early_relation_pool import elimination

PROTOCOL=Path(__file__).with_name('CURVE302_ARITHMETIC_UNIT_CLASS_PROTOCOL.json')
ARITH=r.OUT/'rank_jump_curve302_strict_constructor_arithmetic_v1.json'
ANCHORS=r.OUT/'rank_jump_generic_only_class_anchors_v1.json'
OUTPUT=r.OUT/'rank_jump_curve302_arithmetic_unit_class_v1.json'
WORK=r.ROOT/'artifacts/local/rank-jump-curve302-arithmetic-unit-class-v1'

def compute():
    from sage.all import QQ,ZZ,AA,GF,PolynomialRing,pari,prime_range
    sys.path.insert(0,str(r.ROOT/'elliptic-curves/cas'))
    from research_runtime.local_kummer import LocalSquareclasses
    d=r.read(ARITH);a=r.read(ANCHORS)['cases'][1];policy=r.read(PROTOCOL)
    assert d['generic_strict_dimension']==0 and a['generic_dimension']==17
    pari.allocatemem(64000000,policy['bounds']['pari_stack_bytes'],silent=True)
    R=PolynomialRing(QQ,'z');f=R(d['cubic_ascending']);nf=pari.nfinit([pari(f),d['S_finite']])
    gammas=[pari.Mod(pari(R(x['beta_ascending'])),pari(f)) for x in d['generic_classes']]
    bad=[P for p in d['S_finite'] for P in pari.idealprimedec(nf,p)]
    gcds=[pari.idealadd(nf,ZZ(row['norm']).sqrt(),g) for row,g in zip(d['generic_classes'],gammas)]
    starts=[{'role':'unit_ideal','ideal':pari.idealhnf(nf,1)}]
    for mask0 in a['generic_coefficient_masks']:
        mask=int(mask0);indices=[i for i in range(17) if mask>>i&1];I=pari.idealhnf(nf,1);beta=pari.Mod(1,pari(f))
        for i in indices:I=pari.idealmul(nf,I,gcds[i]);beta*=gammas[i]
        for P in bad:
            v=int(pari.idealval(nf,beta,P));assert v%2==0
            correction=v//2-int(pari.idealval(nf,I,P))
            if correction:I=pari.idealmul(nf,I,pari.idealpow(nf,P,correction))
        assert pari.idealpow(nf,I,2)==pari.idealhnf(nf,beta)
        J,multiplier=pari.idealred(nf,[I,1]);assert pari.idealmul(nf,J,multiplier)==I
        starts.append({'role':'generic_unramified_half_ideal','mask':mask0,'ideal':J})
    local=[LocalSquareclasses(nf,p) for p in d['S_finite']];roots=f.roots(AA,multiplicities=False)
    def polynomial(v):return R([QQ(pari.lift(v).polcoef(i)) for i in range(3)])
    blocks=[]
    for p0 in prime_range(policy['bounds']['proof_prime_start'],policy['bounds']['proof_prime_end']+1):
        p=int(p0)
        if f.discriminant()%p==0:continue
        rr=f.change_ring(GF(p)).roots(multiplicities=False)
        if len(rr)==3 and all(polynomial(g).change_ring(GF(p))(x)!=0 for g in gammas for x in rr):blocks.append((p,list(map(int,rr))))
    def character(g):
        b=polynomial(g);bits=[]
        for p,rr in blocks:
            ff=b.change_ring(GF(p))
            for x in rr:
                val=int(ff(x));assert val;bits.append(int(pow(val,(p-1)//2,p)==p-1))
        return r.pack(bits)
    def signature(g):
        b=polynomial(g)
        return r.pack([int(bit) for L in local for bit in L.signature(g)]+[int(b(x)<0) for x in roots])
    chars=[character(g) for g in gammas];signatures=[signature(g) for g in gammas]
    assert r.rank(chars)==17 and r.rank(signatures)==17
    retained=[];seen_units=set();trials=[];candidates=[]
    def strict_test():
        rank,kernel=elimination(signatures)
        result=[]
        for mask in kernel:
            cc=0
            for i,c in enumerate(chars):
                if mask>>i&1:cc^=c
            if r.reduce(cc,r.basis(chars[:17])):
                result.append({'coefficient_mask':str(mask),'character':str(cc),
                    'generic_indices':[i for i in range(17) if mask>>i&1],
                    'unit_indices':[i-17 for i in range(17,len(chars)) if mask>>i&1]})
        return result
    directions=[[0,i,j] for i in range(0,129,16) for j in range(0,129,16)]
    for start_index,start in enumerate(starts):
        images={};I=start['ideal'];rows=[]
        for direction in directions:
            J,multiplier=pari.idealred(nf,[I,1],direction);assert pari.idealmul(nf,J,multiplier)==I
            alpha=pari.nfbasistoalg(nf,multiplier);key=str(J)
            row={'direction':direction,'reduced_ideal_hnf':key,'multiplier_GP':str(multiplier)}
            if key in images:
                other,other_direction=images[key];unit=alpha/other
                assert pari.idealhnf(nf,unit)==pari.idealhnf(nf,1)
                norm=QQ(pari.nfeltnorm(nf,unit));assert norm in [-1,1]
                if norm<0:unit=-unit
                encoded=str(pari.nfalgtobasis(nf,unit))
                if encoded not in seen_units:
                    seen_units.add(encoded);cc=character(unit)
                    if r.rank(chars+[cc])>r.rank(chars):
                        sig=signature(unit);retained.append({'unit_ascending':list(map(str,polynomial(unit).list())),
                            'start_index':start_index,'directions':[other_direction,direction],
                            'norm_sign_correction':-1 if norm<0 else 1,'character':str(cc),'local_signature':str(sig)})
                        chars.append(cc);signatures.append(sig)
                        print('UNIT',len(retained),'GLOBAL_RANK',r.rank(chars),'LOCAL_RANK',r.rank(signatures),flush=True)
                        candidates=strict_test()
                        if candidates:rows.append(row);break
            else:images[key]=(alpha,direction)
            rows.append(row)
        path=WORK/('start_%02d.json'%start_index)
        r.write_new(path,{'role':start['role'],'generic_mask':start.get('mask'),'starting_ideal_hnf':str(I),'trials':rows})
        trials.append({'path':str(path.relative_to(r.ROOT)),'sha256':r.digest(path.read_bytes())})
        print('START',start_index,'IMAGES',len(images),'UNITS',len(retained),'EXTRA',len(candidates),flush=True)
        if candidates:break
    return {'schema':'rank-jump.curve302-arithmetic-unit-class.v1','status':'PASS',
        'positive_endpoint':'CANDIDATE_PENDING_INDEPENDENT_REPLAY' if candidates else 'NOT_REACHED',
        'retained_units':retained,'strict_candidates':candidates,'generic_plus_units_character_rank':r.rank(chars),
        'generic_plus_units_local_rank':r.rank(signatures),'distinct_unit_values':len(seen_units),'trials':trials,
        'character_blocks':[{'p':p,'roots':rr} for p,rr in blocks],
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in [Path(__file__),PROTOCOL,ARITH,ANCHORS,
            r.ROOT/'elliptic-curves/cas/research_runtime/local_kummer.py']},
        'boundary':'Equation/generic-only arithmetic unit relations and exact local tests. A negative result closes only this bounded attempt. No class-group completeness, rank prediction, rational solubility or large-jump explanation.'}

def capture():
    WORK.mkdir(parents=True,exist_ok=True)
    with (WORK/'worker.log').open('x') as log:
        try:
            p=subprocess.run([sys.executable,__file__,'worker'],stdout=log,stderr=log,timeout=r.read(PROTOCOL)['bounds']['worker_seconds'])
            error='worker failure' if p.returncode else None
        except subprocess.TimeoutExpired:error='bounded timeout'
    if error and not OUTPUT.exists():r.write_new(OUTPUT,{'status':'UNKNOWN','reason':error})
    print(r.read(OUTPUT)['status'],flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['capture','worker']);a=p.parse_args()
    if a.mode=='worker':r.write_new(OUTPUT,compute())
    else:capture()
