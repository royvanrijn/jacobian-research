#!/usr/bin/env python3
"""Three bounded trace-coset checks for a global carrier-solubility criterion."""
import argparse
import ast
from pathlib import Path
import subprocess
import retrospective as r

HERE=Path(__file__).resolve().parent
PROTOCOL=HERE/'NORM_SIX_CARRIER_SOLUBILITY_PROTOCOL.json'
ATLAS=r.ROOT/'artifacts/generated-results/elkies-2026-equation-bisections-full.json'
GRAM=r.ROOT/'elkies-k3/scripts/verify_elkies_2026_published_r17_target.sage'
PAIR=r.OUT/'rank_jump_native_pair_collapse_locus_inputs_v1.json'
PRIOR=r.OUT/'rank_jump_paired_character_moments_verification_v1.json'
INPUT=r.OUT/'rank_jump_norm_six_carrier_solubility_inputs_v1.json'
OUTPUT=r.OUT/'rank_jump_norm_six_carrier_solubility_v1.json'
WORK=r.ROOT/'artifacts/local/rank-jump-norm-six-carrier-solubility-v1'


def capture():
    tree=ast.parse(GRAM.read_text())
    node=next(n for n in tree.body if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='Gpub' for t in n.targets))
    G=ast.literal_eval(node.value.args[1])
    atlas=r.read(ATLAS);labels=['orbit-0911e','orbit-1795d','orbit-11278','orbit-030cb']
    covers={label:next(c for c in atlas['bisections'] if c['label']==label) for label in labels}
    rows=[]
    for name,ls,outcome in [('FG',['orbit-0911e','orbit-1795d'],'YES'),('FD',['orbit-1795d','orbit-11278'],'YES'),('AD',['orbit-030cb','orbit-11278'],'NO')]:
        rows.append({'id':name,'labels':ls,'traces':[covers[x]['published_basis_w'] for x in ls],
                     'quadratics':[covers[x]['residual_chord']['q_coefficients'] for x in ls],
                     'previously_certified_global_solubility':outcome})
    r.write_new(INPUT,{'schema':'rank-jump.norm-six-carrier-solubility-inputs.v1','gram':G,'pairs':rows,
        'known_FG_translate':r.read(PAIR)['generic_word'],
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in (PROTOCOL,ATLAS,GRAM,PAIR,PRIOR)}})


def compute():
    from sage.all import ZZ,QQ,matrix,vector,identity_matrix,pari,PolynomialRing
    inp=r.read(INPUT)
    for path,sha in inp['bindings'].items():assert r.digest((r.ROOT/path).read_bytes())==sha
    G=matrix(ZZ,inp['gram']);assert G.det()==948 and G.is_positive_definite()
    R=PolynomialRing(QQ,'t');rows=[]
    for old in inp['pairs']:
        traces=[vector(ZZ,x) for x in old['traces']];assert all(v*G*v==10 for v in traces)
        q1,q2=map(R,old['quadratics']);assert q1.is_squarefree() and q2.is_squarefree() and q1.gcd(q2)==1
        w=sum(traces);H=(2*identity_matrix(ZZ,17)).stack(matrix(ZZ,[w])).row_module(ZZ).basis_matrix()
        assert H.nrows()==17 and abs(H.det())==2**16
        g=H*G*H.transpose();U=matrix(ZZ,pari(g).qflllgram());assert abs(U.det())==1
        V=U.transpose()*H;gg=V*G*V.transpose()
        raw=pari(gg).qfminim(6);signed_count=int(raw[0]);vecs=matrix(ZZ,raw[2]).columns()
        found=[]
        for v in vecs:
            z=v*V;norm=int(z*G*z)
            assert 0<norm<=6 and all((z[i]-w[i])%2==0 for i in range(17))
            found.append({'vector':list(map(int,z)),'norm':norm})
        assert signed_count==2*len(found)
        has6=any(v['norm']==6 for v in found)
        if old['previously_certified_global_solubility']=='NO':assert not has6
        row=old|{'coset_lattice_basis':[[int(x) for x in row] for row in H],
                 'reduced_basis':[[int(x) for x in row] for row in V],
                 'signed_short_vector_count':signed_count,'unoriented_vectors':found,
                 'has_norm_six':has6,'norm_residue_mod_four':int(w*G*w)%4,
                 'criterion_global_solubility':'YES' if has6 else 'UNKNOWN'}
        if old['id']=='FG':
            S=vector(ZZ,[ZZ(x) for x in inp['known_FG_translate']]);z=2*S-w
            assert z*G*z==6 and any(z==vector(ZZ,v['vector']) or -z==vector(ZZ,v['vector']) for v in found)
            row['known_intersection_vector']=list(map(int,z));row['intersection_number']=int(z*G*z/2-2)
        rows.append(row)
    return {'schema':'rank-jump.norm-six-carrier-solubility.v1','status':'PASS','layer':'solubility','rows':rows,
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in (INPUT,PROTOCOL,Path(__file__),HERE/'retrospective.py')},
        'boundary':r.read(PROTOCOL)['boundary']}


def run():
    WORK.mkdir(parents=True,exist_ok=True);log=WORK/'worker.log';execution=WORK/'execution.json'
    if not log.exists():
        with log.open('x') as out:
            try:
                p=subprocess.run(['/home/royvanrijn/.local/bin/sage','-python',str(Path(__file__).resolve()),'worker'],stdout=out,stderr=out,timeout=60)
                result={'status':'COMPLETE' if p.returncode==0 else 'FAILED','returncode':p.returncode}
            except subprocess.TimeoutExpired:result={'status':'TIMEOUT'}
        r.write_new(execution,result)
    print(r.read(execution));print(log.read_text())


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['capture','run','worker','check']);a=p.parse_args()
    if a.mode=='capture':capture()
    elif a.mode=='run':run()
    else:
        data=compute()
        if a.mode=='worker':r.write_new(OUTPUT,data)
        else:assert data==r.read(OUTPUT)
        for row in data['rows']:print(row['id'],row['signed_short_vector_count'],'signed vectors; norm6',row['has_norm_six'],flush=True)
