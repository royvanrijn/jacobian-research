#!/usr/bin/env python3
"""Read-only audit: reject Laurent input to the polynomial substitution helper.
All checker writes are intercepted; only this audit's own receipt is saved.
"""
from pathlib import Path
import contextlib,hashlib,importlib,io,json,time,traceback
import weighted_verify_utils as utils
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 original=utils.substitution;write=Path.write_text;calls=0;writes=[];bad=[]
 def guard(poly,values,n):
  nonlocal calls
  calls+=1
  if any(any(k<0 for k in e) for e in poly):
   stack=traceback.extract_stack()[:-1];bad.append(dict(caller=str(stack[-1]),stack=[str(f) for f in stack[-5:]],negative_exponents=[list(e) for e in poly if any(k<0 for k in e)][:8]));raise AssertionError('Laurent polynomial passed to polynomial-only substitution')
  return original(poly,values,n)
 def intercept(path,data,*args,**kwargs):writes.append(str(path));return len(data)
 utils.substitution=guard;Path.write_text=intercept;results=[]
 try:
  for name in ['verify_weighted_newfirst7','verify_weighted_unmarked7','verify_weighted_zero7','verify_weighted_123','verify_weighted_213','verify_weighted_113','verify_carrier_reparam12']:
   before=calls;count=len(writes);started=time.monotonic();output=io.StringIO()
   try:
    with contextlib.redirect_stdout(output):importlib.import_module(name).main()
    status='pass'
   except Exception as e:status='error';error=traceback.format_exc()
   row=dict(module=name,source_sha256=sha(HERE/(name+'.py')),status=status,substitution_calls=calls-before,intercepted_write_count=len(writes)-count,seconds=round(time.monotonic()-started,3))
   if status=='error':row['error']=error
   results.append(row)
 finally:Path.write_text=write;utils.substitution=original
 receipt=dict(status='pass' if all(r['status']=='pass' for r in results) else 'domain-failures-found',source_sha256=sha(Path(__file__)),results=results,domain_failures=bad,checker_files_actually_written=[],intercepted_checker_writes=writes,scope='Reject every polynomial substitution call receiving a negative input exponent while replaying the listed current arithmetic checkers. Original receipts and sources are untouched. This checks the domain boundary, not formal proof completeness.')
 (HERE/'substitution-domain-audit.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n')
 for r in results:print(r['module'],r['status'],r['substitution_calls'],r.get('error','').splitlines()[-1] if r.get('error') else '')
 print('DOMAIN FAILURES',len(bad),'INTERCEPTED WRITES',len(writes))
if __name__=='__main__':main()
