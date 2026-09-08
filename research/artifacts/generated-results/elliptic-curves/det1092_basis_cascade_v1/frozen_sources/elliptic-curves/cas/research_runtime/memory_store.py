"""Per-worker immutable facts with portable snapshots, without per-fact files.

Use for cheap finite-field/point signatures in a bounded worker. Durable
checkpoints must retain snapshot(); expensive shared number-field discoveries
continue to use FactStore. Keys and witness hashes match that disk backend.
"""
from copy import deepcopy
from hashlib import sha256
import threading

from .store import FactStore, canonical, digest, encoded


class MemoryFactStore:
    def __init__(self):
        self._records={}
        self._index={}
        self._accessed=set()
        self._lock=threading.RLock()

    _key=FactStore._key

    def get(self,namespace,inputs,*,version='1'):
        with self._lock:
            h=self._index.get(digest(self._key(namespace,inputs,version)))
            if h is None:return None
            self._accessed.add(h)
            return deepcopy(self._records[h]['value'])

    def require(self,namespace,inputs,*,version='1'):
        value=self.get(namespace,inputs,version=version)
        if value is None:raise KeyError(f'missing retained fact: {namespace}')
        return value

    def discover(self,namespace,inputs,build,*,version='1'):
        with self._lock:
            old=self.get(namespace,inputs,version=version)
            if old is not None:return old
            record={'key':self._key(namespace,inputs,version),'value':canonical(build())}
            h=sha256(encoded(record)).hexdigest()
            self._records[h]=record;self._index[digest(record['key'])]=h;self._accessed.add(h)
            return deepcopy(record['value'])

    def snapshot(self):
        with self._lock:
            return {'schema':'elliptic-curves.arithmetic-facts.v1',
                'facts':[{'sha256':h,'record':deepcopy(self._records[h])} for h in sorted(self._accessed)]}

    def import_snapshot(self,snapshot):
        if snapshot.get('schema')!='elliptic-curves.arithmetic-facts.v1':raise ValueError('unknown arithmetic snapshot schema')
        with self._lock:
            validated=[]
            for item in snapshot['facts']:
                record=deepcopy(item['record']);h=item['sha256'];key=record['key']
                if sha256(encoded(record)).hexdigest()!=h:raise ValueError('snapshot content hash mismatch')
                if key!=self._key(key['namespace'],key['inputs'],key['version']):raise ValueError('invalid snapshot key')
                address=digest(key)
                if address in self._index and self._index[address]!=h:raise ValueError('conflicting immutable fact')
                validated.append((address,h,record))
            # Reject conflicting duplicate keys within the snapshot as well.
            local={}
            for address,h,record in validated:
                if address in local and local[address]!=h:raise ValueError('conflicting duplicate snapshot key')
                local[address]=h
            for address,h,record in validated:
                self._index[address]=h;self._records[h]=record;self._accessed.add(h)
