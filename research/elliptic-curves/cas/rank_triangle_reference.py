#!/usr/bin/env python3
"""Reuse already certified equation-support primes for the separate302 reference."""
from pathlib import Path
import json
import rank_triangle_panel as panel
from research_runtime.supervisor import run,Limits

def main():
    source=panel.ART/'curve302_recovered_quotient_local_filtration_v1.json'
    primes=[r['place'] for r in panel.read(source)['local_places'] if isinstance(r['place'],int)]
    row=next(r for r in panel.read(panel.OUT/'panel.json')['rows'] if r['id']=='302')
    out=panel.OUT/'arithmetic/302-cached-support';out.mkdir(parents=True,exist_ok=True)
    assert not (out/'protocol.json').exists()
    gp=out/'field.gp'
    gp.write_text(f'P={primes};for(i=1,#P,if(!isprime(P[i]),error("prime hint")));addprimes(P);\n'+panel.program(row,out))
    panel.save(out/'protocol.json',{'script_sha256':panel.digest(Path(__file__)),
        'prime_hint_source_sha256':panel.digest(source),'equation_program_sha256':panel.digest(gp),
        'primes':primes,'wall_seconds':15,'rss_bytes':3*2**30,'point_inputs':False,
        'purpose':'The unhinted reference attempt timed out in factorization. Reuse retained equation-support primes, not point-generated class information. No repeated BNF probe.'})
    sup=run([str(panel.GP),'-q','-f','-s',str(2*2**30),str(gp)],limits=Limits(15,3*2**30,pari_stack_bytes=2*2**30),
            log_path=out/'field.log',checkpoint_path=out/'field.supervisor.json',cwd=panel.ROOT)
    log=(out/'field.log').read_text();result={'id':'302','status':'UNKNOWN_FIELD','supervision':sup,
         'log_sha256':panel.digest(out/'field.log'),'g_upper':None,'grh_g':None,'rank_upper':None}
    if sup['outcome']=='completed' and 'DONE_FIELD|1' in log and 'FAIL|' not in log:
        f=dict(line.split('|',1) for line in log.splitlines() if '|' in line)
        local=[json.loads(line.split('|')[1]) for line in log.splitlines() if line.startswith('LOCAL|')]
        result.update(status='PASS_EQUATION_LOCAL_ARITHMETIC',bk_offset=int(f['BK_OFFSET']),
          field_discriminant=f['FIELD_DISC'],field_index=f['FIELD_INDEX'],signature=json.loads(f['SIGNATURE']),
          root_number=int(f['ROOT_NUMBER']),conductor=f['CONDUCTOR'],minimal_model=json.loads(f['MINIMAL_MODEL']),
          cubic=json.loads(f['CUBIC']),local=local,
          local_2_descent_dimensions={str(r[0]):int(f['LOCAL2_DIM']) if r[0]==2 else r[4]-1 for r in local},
          dyadic_local_dimension=int(f['LOCAL2_DIM']),real_local_dimension=int(f['REAL_LOCAL_DIM']),
          split=[line.split('|')[1:] for line in log.splitlines() if line.startswith('SPLIT|')])
    panel.save(out/'result.json',result);print(json.dumps(result,indent=2))
if __name__=='__main__':main()
