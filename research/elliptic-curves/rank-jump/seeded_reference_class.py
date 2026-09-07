#!/usr/bin/env python3
"""Inject equation-only principal relations into an isolated PARI class computation."""
import argparse
from pathlib import Path
import re
import subprocess
import time
import retrospective as r

HERE=Path(__file__).resolve().parent
PROTOCOL=HERE/'SEEDED_REFERENCE_CLASS_PROTOCOL.json'
REFERENCE=r.OUT/'rank_jump_reference_strict_class_construction_inputs_v1.json'
POOL=r.OUT/'rank_jump_early_relation_pool_inputs_v1.json'
INPUT=r.OUT/'rank_jump_seeded_reference_class_inputs_v1.json'
OUTPUT=r.OUT/'rank_jump_seeded_reference_class_v1.json'
WORK=r.ROOT/'artifacts/local/rank-jump-seeded-reference-class-v1'
UPSTREAM=r.ROOT/'artifacts/local/rank-jump-reference-strict-class-construction-v1/source-audit/buch2.c'
RUNTIME=Path('/home/royvanrijn/.local/share/jacobian-sage-10.9')

INJECT='''
  if (RJ_SEEDS && TRIES == 1)
  {
    long j;
    pre_allocate(&cache, lg(RJ_SEEDS)+RELSUP+RU);
    for (j = 1; j < lg(RJ_SEEDS); j++)
      try_elt(&cache, &F, nf, gel(RJ_SEEDS,j), fact);
    err_printf("RJ_SEEDED %ld MISSING %ld\\n", cache.last-cache.base, cache.missing);
    need = 0;
  }
'''
MAIN='''
int main(int argc, char **argv)
{
  GEN inp, nf, answer;
  if (argc != 2) return 2;
  pari_init(64000000, 500000);
  paristack_setsize(64000000, 268435456);
  DEBUGLEVEL = 1;
  setrand(stoi(20260908));
  inp = gp_read_file(argv[1]);
  nf = nfinit0(gel(inp,1), 0, DEFAULTPREC);
  RJ_SEEDS = gclone(gel(inp,2));
  err_printf("RJ_BASIS %Ps\\n", nf_get_zk(nf));
  answer = Buchall_param(nf, 0., 0., 4, nf_FORCE, DEFAULTPREC);
  pari_printf("%Ps\\n", answer);
  gunclone(RJ_SEEDS);
  pari_close();
  return 0;
}
'''

def bindings(paths):
    return {str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in paths}

def prepare():
    from sage.all import QQ,PolynomialRing,pari
    WORK.mkdir(parents=True,exist_ok=True)
    ref=r.read(REFERENCE);pool=r.read(POOL)
    assert list(map(str,ref['cubic_ascending']))==list(map(str,pool['cubic_ascending']))
    R=PolynomialRing(QQ,'z');f=R(ref['cubic_ascending'])
    nf=pari.nfinit([pari(f),ref['S_finite']]);assert str(nf.disc())==ref['field_discriminant']
    seeds=[]
    for row in pool['relations']:
        c=pari.nfalgtobasis(nf,pari.Mod(pari(R(row['alpha_ascending'])),pari(f)))
        den=pari.denominator(c);c*=den
        content=pari.content(c);c/=content
        assert all(x.type()=='t_INT' for x in c)
        seeds.append(list(map(str,c)))
    assert len(seeds)==4134
    data={'schema':'rank-jump.seeded-reference-class-inputs.v1','cubic_ascending':list(map(str,f.list())),
          'S_finite':ref['S_finite'],'field_discriminant':ref['field_discriminant'],
          'integral_basis_GP':list(map(str,nf.nf_get_zk())), 'primitive_seed_columns':seeds,
          'bindings':bindings([REFERENCE,POOL,PROTOCOL,Path(__file__)])}
    r.write_new(INPUT,data)
    cols='['+','.join('['+','.join(s)+']~' for s in seeds)+']'
    (WORK/'reference.gp').write_text('[['+str(f)+','+str(ref['S_finite'])+'],'+cols+']\n')
    (WORK/'control.gp').write_text('[[z^3-11,[3,11]],[[1,1,0]~,[2,1,0]~]]\n')
    print('Prepared',len(seeds),'primitive equation-only seed elements',flush=True)

def compile_worker():
    policy=r.read(PROTOCOL);raw=UPSTREAM.read_bytes()
    assert r.digest(raw)==policy['upstream']['buch2_c_sha256']
    src=raw.decode();names=re.findall(r'\n(?:GEN|void|long|int)\n(\w+)\(',src)
    prefix=''.join('#define '+n+' rank_jump_seeded_'+n+'\n' for n in names)
    prefix+='#include "pari.h"\nstatic GEN RJ_SEEDS = NULL;\n'
    marker='  if (computed)\n  {';assert src.count(marker)==1
    src=src.replace(marker,INJECT+'\n'+marker)
    marker='  k = add_rel_i(cache, R, nz, m, 0, 0, &rel, in_rnd_rel);';assert src.count(marker)==1
    src=src.replace(marker,marker+'\n  if (k > 0 && m && typ(m) != t_INT) err_printf("RJ_ELEMENT %Ps\\n", m);')
    modified=WORK/'seeded_buch2.c';binary=WORK/'seeded_buch2'
    with modified.open('x') as stream:stream.write(prefix+src+MAIN)
    cmd=['cc','-O2','-I'+str(RUNTIME/'include/pari'),str(modified),'-L'+str(RUNTIME/'lib'),
         '-Wl,-rpath,'+str(RUNTIME/'lib'),'-lpari','-lm','-o',str(binary)]
    proc=subprocess.run(cmd,capture_output=True,text=True,timeout=policy['bounds']['compile_seconds'])
    result={'returncode':proc.returncode,'stdout':proc.stdout,'stderr':proc.stderr,'command':cmd,
            'modified_source_sha256':r.digest(modified.read_bytes()),'renamed_exports':names,
            'runtime_sha256':r.digest((RUNTIME/'lib/libpari-gmp-tls.so.2.17.3').read_bytes())}
    r.write_new(WORK/'compile.json',result);assert proc.returncode==0,result

def run_one(token,seconds):
    start=time.monotonic()
    with (WORK/(token+'.log')).open('x') as log, (WORK/(token+'-bnf.gp')).open('x') as out:
        try:
            p=subprocess.run([str(WORK/'seeded_buch2'),str(WORK/(token+'.gp'))],
                             stdout=out,stderr=log,timeout=seconds)
            status='RETURNED_UNCERTIFIED' if p.returncode==0 else 'UNKNOWN'
            reason=None if p.returncode==0 else 'worker failure'
        except subprocess.TimeoutExpired:status='UNKNOWN';reason='bounded timeout'
    data={'status':status,'reason':reason,'elapsed_seconds':time.monotonic()-start}
    r.write_new(WORK/(token+'-terminal.json'),data);return data

def capture():
    from sage.all import pari
    policy=r.read(PROTOCOL);assert r.read(INPUT)['bindings']==bindings([REFERENCE,POOL,PROTOCOL,Path(__file__)])
    compile_worker()
    control=run_one('control',policy['bounds']['control_seconds'])
    assert control['status']=='RETURNED_UNCERTIFIED',control
    bnf=pari((WORK/'control-bnf.gp').read_text());assert list(map(int,bnf.bnf_get_cyc()))==[2]
    control['ordinary_class_group']=[2];control['bnfcertify']=int(pari.bnfcertify(bnf));assert control['bnfcertify']==1
    print('Control PASS; starting bounded seeded reference',flush=True)
    reference=run_one('reference',policy['bounds']['reference_seconds'])
    log=(WORK/'reference.log').read_text();elements=[x[len('RJ_ELEMENT '):] for x in log.splitlines() if x.startswith('RJ_ELEMENT ')]
    accepted=WORK/'reference-accepted-elements.json';r.write_new(accepted,{'basis_GP':r.read(INPUT)['integral_basis_GP'],'elements_GP':elements})
    result={'schema':'rank-jump.seeded-reference-class.v1','control':control,'reference':reference,
            'seed_checkpoint_lines':[x for x in log.splitlines() if x.startswith('RJ_SEEDED ')],
            'accepted_principal_element_occurrences':len(elements),
            'explicit_additional_strict_class':'NOT_YET_CONSTRUCTED',
            'bindings':bindings([Path(__file__),PROTOCOL,INPUT,WORK/'compile.json',WORK/'reference.log',accepted,WORK/'reference-bnf.gp']),
            'boundary':'Principal relations are construction inputs, not strict classes. No class-group completeness, rank, absence, rationality or Sha claim.'}
    r.write_new(OUTPUT,result);print(reference,flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['prepare','capture']);a=p.parse_args();globals()[a.mode]()
