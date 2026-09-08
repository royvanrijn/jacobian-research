#!/usr/bin/env python3
"""Landlock allowlist + new network namespace, enforced before starting Sage.

Usage: unshare -Unpf python3 isolate.py WORKER FIXTURE OUTPUT SCRATCH [selftest]
No search imports or mathematical choices here. ABI>=5 required, fail closed.
"""
import ctypes,os,sys
from pathlib import Path
worker,fixture,output,scratch=map(lambda s:str(Path(s).resolve()),sys.argv[1:5])
selftest=len(sys.argv)>5 and sys.argv[5]=='selftest'
runtime='/home/royvanrijn/.local/share/jacobian-sage-10.9'
lib=ctypes.CDLL(None,use_errno=True)
def call(n,*args):
    ans=lib.syscall(n,*args)
    if ans<0:raise OSError(ctypes.get_errno(),os.strerror(ctypes.get_errno()))
    return ans
abi=call(444,0,0,1)
assert abi>=5,abi
class Ruleset(ctypes.Structure):_fields_=[('fs',ctypes.c_uint64),('net',ctypes.c_uint64)]
class Rule(ctypes.Structure):
    _pack_=1
    _fields_=[('access',ctypes.c_uint64),('fd',ctypes.c_int32)]
fs=(1<<16)-1
ruleset=call(444,ctypes.byref(Ruleset(fs,3)),ctypes.sizeof(Ruleset),0)
read=1|4|8
allowed=[('/usr',read),('/lib',read),('/lib64',read),('/bin',read),(runtime,read),
         ('/etc',read),('/dev/null',2|4|32768),('/dev/urandom',4),('/dev/random',4),
         ('/proc/self',read),('/proc/meminfo',4),('/proc/cpuinfo',4),('/sys/devices/system/cpu',read),
         (worker,4),(output,fs),(scratch,fs)]
if not selftest:allowed.append((fixture,4))
for path,rights in allowed:
    if not Path(path).exists():continue
    fd=os.open(path,os.O_PATH|os.O_CLOEXEC)
    call(445,ruleset,1,ctypes.byref(Rule(rights,fd)),0);os.close(fd)
if lib.prctl(38,1,0,0,0)!=0:raise OSError(ctypes.get_errno(),'no_new_privs')
call(446,ruleset,0);os.close(ruleset)
os.chdir(scratch)
env={'PATH':runtime+'/bin:/usr/bin:/bin','HOME':scratch,'TMPDIR':scratch,'DOT_SAGE':scratch+'/sage',
     'PYTHONNOUSERSITE':'1','OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1','MKL_NUM_THREADS':'1',
     'BLIND_OUTPUT':output,'BLIND_FIXTURE':fixture}
os.execve(runtime+'/bin/python',[runtime+'/bin/python',worker]+(['--selftest'] if selftest else []),env)
