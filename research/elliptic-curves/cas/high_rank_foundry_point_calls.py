"""Durable point-call intents and exact-result reuse after process interruption."""
import hashlib
import json
from pathlib import Path
import time
from v3_warm_support import read, atomic


def install(job, backend):
    original=backend.execute
    def execute(search,mapping,height,seconds,gp_hash):
        request={'state_key':search.state.key,'mapping':mapping,'height':height,
                 'seconds':seconds,'gp_sha256':gp_hash,'sources':backend.sources()}
        request=json.loads(json.dumps(request))
        key=hashlib.sha256(json.dumps(request,sort_keys=True,separators=(',',':')).encode()).hexdigest()
        folder=Path(job)/'point-invocations'/key
        if (folder/'result.json').exists():
            record=read(folder/'result.json')
            return record,backend.replay(search,mapping,record)
        folder.mkdir(parents=True,exist_ok=True)
        atomic(folder/f'attempt-{time.time_ns()}.json',{'request':request,'started_at':time.time()},immutable=True)
        record,points=original(search,mapping,height,seconds,gp_hash)
        # Publish the complete raw point return before admission or chart handling.
        # If the process dies next, replay this result without another GP call.
        atomic(folder/'result.json',record,immutable=True)
        if backend.replay(search,mapping,record)!=points:
            raise ArithmeticError('journalled point return differs from exact replay')
        return record,points
    backend.execute=execute
