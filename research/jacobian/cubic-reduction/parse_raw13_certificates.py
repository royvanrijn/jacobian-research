#!/usr/bin/env python3
"""Parse retained complete Singular certificates directly into sparse rational polynomials."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json,re,time
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def parse(text,variables):
 names={v:i for i,v in enumerate(variables)};out={}
 if text=='0':return []
 assert '(' not in text and ')' not in text
 terms=re.findall(r'[+-]?[^+-]+',text);assert ''.join(terms)==text
 for term in terms:
  sign=Q(1)
  if term[0] in '+-':
   if term[0]=='-':sign=-sign
   term=term[1:]
  exponent=[0]*len(variables);coefficient=sign
  for factor in term.split('*'):
   name,power=(factor.split('^')+[None])[:2] if '^' in factor else (factor,None)
   if name in names:exponent[names[name]]+=int(power) if power is not None else 1
   else:assert power is None;coefficient*=Q(name)
  e=tuple(exponent);out[e]=out.get(e,Q(0))+coefficient
 return [[list(e),str(c)] for e,c in sorted(out.items(),reverse=True) if c]
def main():
 started=time.monotonic();out=HERE/'raw13-rigidity';rawpath=HERE/'raw13-specialization/input.json';normalpath=HERE/'raw13-normalized/input.json';raw=json.loads(rawpath.read_text());normal=json.loads(normalpath.read_text());manifest={}
 for tag,data,count in [('raw-unit-0',raw,1),('raw-unit-1',raw,1),('normalized-squares',normal,13)]:
  log=out/(tag+'.log');lines=log.read_text().splitlines();position=0;certs=[]
  for i in range(count):
   header=lines[position].split(':');assert header[0]=='CERT' and int(header[1])==i;position+=1;ngens=int(header[2]);mult=[parse(line,data['variables']) for line in lines[position:position+ngens]];position+=ngens
   if tag.startswith('raw-unit'):target=[[[0]*len(data['variables']),'1']]
   else:
    base=Q(raw['atoms'][i+2]['coefficient']) if i<11 else Q(1);e2=[0]*13;e2[i]=2;e1=[0]*13;e1[i]=1;target=[[e2,'1'],[e1,str(-2*base)],[[0]*13,str(base*base)]]
   certs.append(dict(target=target,multipliers=mult))
  assert position==len(lines)
  payload=dict(source_sha256=sha(Path(__file__)),exporter_sha256=sha(HERE/'certify_raw13_rigidity.py'),singular_source_sha256=sha(out/(tag+'.sing')),singular_log_sha256=sha(log),variables=data['variables'],certificates=certs)
  path=out/(tag+'-parsed.json');path.write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n');manifest[tag]=dict(path=path.name,sha256=sha(path));print('PARSED',tag,flush=True)
 receipt=dict(status='exact-primary-identities',source_sha256=sha(Path(__file__)),exporter_sha256=sha(HERE/'certify_raw13_rigidity.py'),raw_input_sha256=sha(rawpath),normalized_input_sha256=sha(normalpath),certificates=manifest,identities=15,seconds=round(time.monotonic()-started,3),scope='Direct sparse parsing of completed Singular logs. Independent polynomial replay is still required.')
 (out/'certification-direct.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');print(json.dumps(receipt),flush=True)
if __name__=='__main__':main()
