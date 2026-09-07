#!/usr/bin/env python3
"""Memory-bounded relation-only variant; no compact transformation matrix."""
import argparse
from pathlib import Path
import re
import subprocess
import time
import retrospective as r
import seeded_reference_class as previous

PROTOCOL=Path(__file__).with_name('SEEDED_REFERENCE_RELATIONS_PROTOCOL.json')
WORK=r.ROOT/'artifacts/local/rank-jump-seeded-reference-relations-v1'
OUTPUT=r.OUT/'rank_jump_seeded_reference_relations_v1.json'

def capture():
    from sage.all import pari
    p=r.read(PROTOCOL);WORK.mkdir(parents=True,exist_ok=True)
    data=r.read(previous.INPUT)
    assert data['bindings']==previous.bindings([previous.REFERENCE,previous.POOL,previous.PROTOCOL,Path(previous.__file__)])
    src=previous.UPSTREAM.read_text()
    assert r.digest(src.encode())==r.read(previous.PROTOCOL)['upstream']['buch2_c_sha256']
    names=re.findall(r'\n(?:GEN|void|long|int)\n(\w+)\(',src)
    prefix=''.join('#define '+n+' rank_jump_relations_'+n+'\n' for n in names)
    prefix+='#include "pari.h"\nstatic GEN RJ_SEEDS = NULL;\n'
    marker='  if (computed)\n  {';assert src.count(marker)==1
    inject=previous.INJECT.replace('      try_elt(&cache, &F, nf, gel(RJ_SEEDS,j), fact);',
        '    { pari_sp seed_av = avma; try_elt(&cache, &F, nf, gel(RJ_SEEDS,j), fact); set_avma(seed_av); }')
    src=src.replace(marker,inject+'\n'+marker)
    marker='  k = add_rel_i(cache, R, nz, m, 0, 0, &rel, in_rnd_rel);';assert src.count(marker)==1
    src=src.replace(marker,marker+'\n  if (k > 0 && m && typ(m) != t_INT) err_printf("\\nRJ_ELEMENT %Ps\\n", m);')
    main=previous.MAIN.replace('nf_FORCE, DEFAULTPREC','0, DEFAULTPREC')
    modified=WORK/'relations_buch2.c';binary=WORK/'relations_buch2'
    with modified.open('x') as stream:stream.write(prefix+src+main)
    rt=previous.RUNTIME
    cmd=['cc','-O2','-I'+str(rt/'include/pari'),str(modified),'-L'+str(rt/'lib'),
         '-Wl,-rpath,'+str(rt/'lib'),'-lpari','-lm','-o',str(binary)]
    comp=subprocess.run(cmd,capture_output=True,text=True,timeout=30)
    compiled={'returncode':comp.returncode,'stdout':comp.stdout,'stderr':comp.stderr,'command':cmd,
        'modified_source_sha256':r.digest(modified.read_bytes()),'renamed_exports':names,
        'runtime_sha256':r.digest((rt/'lib/libpari-gmp-tls.so.2.17.3').read_bytes())}
    r.write_new(WORK/'compile.json',compiled);assert comp.returncode==0,compiled
    terminals={}
    for token in ['control','reference']:
        start=time.monotonic()
        with (WORK/(token+'.log')).open('x') as log,(WORK/(token+'-bnf.gp')).open('x') as out:
            try:
                run=subprocess.run([str(binary),str(previous.WORK/(token+'.gp'))],
                    stdout=out,stderr=log,timeout=p['bounds'][token+'_seconds'])
                status='RETURNED_UNCERTIFIED' if run.returncode==0 else 'UNKNOWN'
                reason=None if run.returncode==0 else 'worker failure'
            except subprocess.TimeoutExpired:status='UNKNOWN';reason='bounded timeout'
        terminal={'status':status,'reason':reason,'elapsed_seconds':time.monotonic()-start}
        r.write_new(WORK/(token+'-terminal.json'),terminal);terminals[token]=terminal
        if token=='control':
            assert status=='RETURNED_UNCERTIFIED',terminal
            b=pari((WORK/(token+'-bnf.gp')).read_text());assert list(map(int,b.bnf_get_cyc()))==[2]
            assert pari.bnfcertify(b)==1
            print('Control PASS; starting relation-only reference',flush=True)
    log=(WORK/'reference.log').read_text()
    elements=re.findall(r'RJ_ELEMENT (\[[^\n]+\]~)',log)
    accepted=WORK/'accepted-elements.json'
    r.write_new(accepted,{'basis_GP':data['integral_basis_GP'],'elements_GP':elements})
    result={'schema':'rank-jump.seeded-reference-relations.v1','terminals':terminals,
        'seed_checkpoint_lines':[x for x in log.splitlines() if x.startswith('RJ_SEEDED ')],
        'accepted_principal_element_occurrences':len(elements),
        'explicit_additional_strict_class':'NOT_YET_CONSTRUCTED',
        'bindings':previous.bindings([Path(__file__),PROTOCOL,Path(previous.__file__),previous.INPUT,
            WORK/'compile.json',WORK/'reference.log',accepted,WORK/'reference-bnf.gp']),
        'boundary':'Exact accepted elements retained for independent parity-kernel extraction. Floating class-group output gives no unverified completeness claim.'}
    r.write_new(OUTPUT,result);print(terminals['reference'],len(elements),'accepted elements',flush=True)

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('mode',choices=['capture']);a.parse_args();capture()
