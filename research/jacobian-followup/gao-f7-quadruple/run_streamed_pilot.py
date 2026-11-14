"""One bounded modular local-DAG pilot, with no characteristic-zero inference."""
from pathlib import Path
import gzip,hashlib,json,subprocess,resource,time
HERE=Path(__file__).resolve().parent;START=time.monotonic();prime=1000000007
command=['/usr/local/bin/julia','--startup-file=no','--project='+str(HERE.parent/'gao-f7-three-point-certificate/julia'),str(HERE/'generate_streamed_dag.jl'),str(prime)]
def limits():
 resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3))
 resource.setrlimit(resource.RLIMIT_FSIZE,(30*1024**2,30*1024**2))
with (HERE/'streamed-dag-pilot.log').open('wb') as fd:
 try:
  proc=subprocess.run(command,stdout=fd,stderr=subprocess.STDOUT,timeout=90,preexec_fn=limits)
  outcome='completed';code=proc.returncode
 except subprocess.TimeoutExpired:outcome='timeout';code=None
artifact=HERE/('raw-dag-streamed-'+str(prime)+'.txt.gz');stats={'nodes':0,'polynomial_terms':0,'relation_terms':0}
if artifact.exists():
 with gzip.open(artifact,'rt') as fd:
  for line in fd:
   if line.startswith('NODE '):stats['nodes']+=1
   elif line.startswith('P '):stats['polynomial_terms']+=1
   elif line.startswith('R '):stats['relation_terms']+=1
receipt={'outcome':outcome,'returncode':code,'prime':prime,'seconds':time.monotonic()-START,'cap_seconds':90,'address_cap_GiB':4,'per_file_output_cap_MiB':30,'artifact_bytes':artifact.stat().st_size if artifact.exists() else 0,**stats,'theorem_exclusion':False,'boundary':'Untrusted modular DAG only. Full characteristic-zero identities and a separate checker are required. Streaming gzip addresses the measured raw-output cap. No automatic prime campaign follows this pilot.'}
(HERE/'streamed-dag-pilot.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt));print((HERE/'streamed-dag-pilot.log').read_text()[-1700:])
