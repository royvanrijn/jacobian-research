#!/usr/bin/env python3
"""Replay polynomial hypotheses, model transports and generic Kummer witnesses."""
import argparse
from fractions import Fraction as Q
from math import isqrt,lcm
from pathlib import Path
import retrospective as r
import generic_selmer_capacity as source
import complete_generic_selmer_geometry as completion
import verify_fixed_field_transfer_geometry as vg
from verify_fresh_symbolic_discriminant import value


def compute():
    gcd_witnesses=[]
    def coprime_mod_prime(a,b):
        # A degree-preserving coprime reduction certifies coprimality over Q.
        # Avoid coefficient swell from rational Euclid on historical models.
        a=list(map(Q,a));b=list(map(Q,b))
        for p in (1009,1013,1019,1021,1031,1033,1039,1049):
            assert all(p%d for d in range(2,isqrt(p)+1))
            if any(x.denominator%p==0 for x in a+b):continue
            aa=[x.numerator*pow(x.denominator,-1,p)%p for x in a]
            bb=[x.numerator*pow(x.denominator,-1,p)%p for x in b]
            if not aa[-1] or not bb[-1]:continue
            while bb!=[0]:
                while len(aa)>=len(bb) and aa!=[0]:
                    k=len(aa)-len(bb);c=aa[-1]*pow(bb[-1],-1,p)%p
                    for j,v in enumerate(bb):aa[j+k]=(aa[j+k]-c*v)%p
                    while len(aa)>1 and aa[-1]==0:aa.pop()
                aa,bb=bb,aa
            if len(aa)==1:
                gcd_witnesses.append({'degrees':[len(a)-1,len(b)-1],'prime':p});return [Q(1)]
        raise AssertionError('No bounded coprimality witness')
    inp=r.read(source.INPUT);raw=r.read(source.OUTPUT)
    completed=r.read(completion.OUTPUT);cinp=r.read(completion.INPUT)
    for obj in (completed,cinp):
        for p,sha in obj['bindings'].items():assert r.digest((r.ROOT/p).read_bytes())==sha
    old=next(x for x in inp['families'] if x['family']==completion.FAMILY)
    new=cinp['families'][0]
    assert new['A']==old['A'][::-1] and new['B']==old['B'][::-1]
    assert Q(new['irreducibility_witness_parameter'])*Q(old['irreducibility_witness_parameter'])==1
    assert next(x for x in raw['rows'] if x['family']==completion.FAMILY)['status']=='UNKNOWN'
    # In-memory adapters let the immutable polynomial verifier check all seven
    # rows, including the explicitly verified reversed coordinate presentation.
    replay_input={**inp,'families':[new if x['family']==completion.FAMILY else x for x in inp['families']]}
    replay_output={**raw,'rows':[completed['rows'][0] if x['family']==completion.FAMILY else x for x in raw['rows']]}
    saved=vg.source.INPUT,vg.source.OUTPUT,vg.source.PROVENANCE
    original_gcd=vg.gcd
    read=r.read
    try:
        vg.source.INPUT=source.INPUT;vg.source.OUTPUT=source.OUTPUT;vg.source.PROVENANCE=source.PROVENANCE
        def replay_read(path):
            if path==source.INPUT:return replay_input
            if path==source.OUTPUT:return replay_output
            return read(path)
        r.read=replay_read
        vg.gcd=coprime_mod_prime
        geometric_replay=vg.compute()
    finally:
        vg.source.INPUT,vg.source.OUTPUT,vg.source.PROVENANCE=saved
        r.read=read
        vg.gcd=original_gcd
    families=[]
    for x in replay_output['rows']:
        # I_n has n components at each geometric zero of multiplicity n.
        root_rank=sum(q['degree']*(q['multiplicity']-1) for q in x['squarefree_factors'])
        euler=sum(q['degree']*q['multiplicity'] for q in x['squarefree_factors'])
        assert euler==24 and x['geometric_monodromy']=='S3'
        capacity=22-2-root_rank
        # Enumerate every possible characteristic-zero K3 Picard rank.
        picard_checks=[]
        for rho in range(2+root_rank,21):
            mw=rho-2-root_rank;sha=22-rho
            assert mw+sha==capacity
            picard_checks.append({'rho':rho,'geometric_MW_rank':mw,'geometric_Sha_two_dimension':sha})
        families.append({'family':x['family'],'euler_number':euler,'root_lattice_rank':root_rank,
            'geometric_two_torsion_dimension':0,'geometric_selmer_dimension':capacity,
            'arithmetic_global_pool_dimension_upper_bound':capacity,'picard_rank_checks':picard_checks,
            'theorem_inputs':'Kummer exact sequence, Shioda-Tate, geometric Sha=Br and dim Br[2]=22-rho; see proof note.'})
    by_family={x['family']:x for x in families};equations={x['family']:x for x in inp['families']};cases=[]
    for x in inp['cases']:
        f=equations[x['family']];t=Q(x['parameter']);Af=value(list(map(Q,f['A'])),t);Bf=value(list(map(Q,f['B'])),t)
        model,points=r.short(x['model'],x['generic_sections']);As,Bs=map(Q,model[3:])
        assert Af*Bf*As*Bs and 4*Af**3+27*Bf**2
        q=(Bs/Bf)/(As/Af);assert q>0
        un,ud=isqrt(q.numerator),isqrt(q.denominator)
        assert un*un==q.numerator and ud*ud==q.denominator
        u=Q(un,ud);assert As==Af*u**4 and Bs==Bf*u**6
        scale=lcm(*(Q(a).denominator for a in model));Am=As*scale**4;Bm=Bs*scale**6
        scaled_model=['0','0','0',str(Am),str(Bm)]
        scaled_points=[[str(Q(a)*scale**2),str(Q(b)*scale**3)] for a,b in points]
        blocks=[]
        for p,roots in x['independence_blocks']:
            assert p>2 and all(p%d for d in range(2,isqrt(p)+1))
            actual=r.roots_at(str(Am),str(Bm),p)
            assert actual is not None and list(actual)==roots
            blocks.append((p,tuple(roots)))
        sigs=[r.point_signature(scaled_model,P,blocks) for P in scaled_points]
        assert sigs==x['independence_signatures'] and r.rank(sigs)==len(points)
        d=by_family[x['family']]['geometric_selmer_dimension'];assert d>=len(points)
        cases.append({'token':x['token'],'family':x['family'],'parameter':x['parameter'],
            'short_model_scale_u':str(u),'generic_mod_two_dimension':len(points),
            'geometric_selmer_dimension':d,'extra_global_pool_capacity_upper_bound':d-len(points),
            'generic_good_prime_witness_count':len(blocks)})
    paths=[Path(__file__),Path(source.__file__),source.PROTOCOL,source.INPUT,source.OUTPUT,source.PROVENANCE,
        Path(vg.__file__),Path(r.__file__),Path(__file__).with_name('verify_fresh_symbolic_discriminant.py'),
        Path(completion.__file__),completion.INPUT,completion.OUTPUT]
    return {'schema':'rank-jump.generic-selmer-capacity-verification.v1','status':'PASS',
        'families':families,'cases':cases,'geometry_replay_rows':geometric_replay['rows'],
        'modular_coprimality_witnesses':gcd_witnesses,
        'bindings':source.bind(paths),
        'boundary':'Portable verification of hypotheses, exact model isomorphisms, generic mod-two independence and dimension arithmetic. The cohomological theorem is proved in the accompanying note, not mechanically formalized. No exceptional points or rank labels read; no arithmetic Selmer basis computed.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);args=p.parse_args();result=compute()
    if args.mode=='build':r.write_new(source.VERIFICATION,result)
    else:assert result==r.read(source.VERIFICATION)
    print('PASS',[(x['family'],x['geometric_selmer_dimension']) for x in result['families']],len(result['cases']),'specialization/generic checks',flush=True)
