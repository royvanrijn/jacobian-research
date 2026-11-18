#!/usr/bin/env python3
"""Independent finite membership authentication for primitive -8/-11 norms.

Python standard library only. Enumerate complete integral ellipses with gcd,
without prime-support sieving. Import each original checker solely to obtain
its actual sieve bytearray, then compare WHOLE BYTEARRAYS by exact equality at
the ACTUAL certificate bounds. Hashes identify inputs/results after equality;
they are not used as a substitute for byte comparison. No sumset masks,
spectral bounds or unrestricted classifications are certified here.
"""
import argparse
from hashlib import sha256
import json
import importlib.util
from math import gcd, isqrt
from pathlib import Path


def file_digest(path):
    return sha256(path.read_bytes()).hexdigest()


def find_root(script):
    for parent in script.parents:
        if (parent/'AGENTS.md').is_file() and (parent/'research/arithmetic').is_dir():
            return parent
    raise ValueError('Cannot locate repository; supply --repo-root.')


def coordinate_flags(d,bound):
    flags=bytearray(bound+1)
    visited=primitive=0
    if d==8:
        # Independent signs preserve Q and gcd; odd norms have odd x.
        for y in range(isqrt(bound//2)+1):
            twice_y_square=2*y*y
            for x in range(1,isqrt(bound-twice_y_square)+1,2):
                visited+=1
                if gcd(x,y)==1:
                    flags[x*x+twice_y_square]=1
                    primitive+=1
    elif d==11:
        # Simultaneous signs permit y>=0. Reflection x->-x-y preserves
        # Q and gcd and permits 2x+y>=0. 4Q=(2x+y)^2+11y^2.
        for y in range(isqrt(4*bound//11)+1):
            s=isqrt(4*bound-11*y*y)
            for x in range(-(y//2),(s-y)//2+1):
                visited+=1
                q=x*x+x*y+3*y*y
                if q>0 and q%2 and q%11 and gcd(x,y)==1:
                    flags[q]=1
                    primitive+=1
    else:
        raise ValueError('Unsupported discriminant.')
    if flags[0]!=0 or flags[1]!=1:
        raise ValueError('Zero/unit failure.')
    return flags,visited,primitive


def authenticate_flags(candidate,reference):
    if candidate!=reference:
        raise ValueError('Sieve and coordinate membership bytearrays differ.')


def load_sieve_checker(path,d):
    spec=importlib.util.spec_from_file_location('_primitive_norm_sieve_minus'+str(d),path)
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def negative_checks(candidate,reference,false_value):
    rejected=[]
    for label,position,original,corrupt in [('inert_prime_power_injected',false_value,0,1),
                                          ('unit_deleted',1,1,0)]:
        if candidate[position]!=original or reference[position]!=original:
            raise ValueError('Unexpected original mutation bit.')
        candidate[position]=corrupt
        try:
            authenticate_flags(candidate,reference)
        except ValueError:
            rejected.append(label)
        else:
            raise ValueError('Corrupted membership flags accepted: '+label)
        finally:
            candidate[position]=original
    authenticate_flags(candidate,reference)
    return rejected


def main():
    script=Path(__file__).resolve()
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo-root',type=Path)
    parser.add_argument('--output',type=Path,help='Optional JSON receipt; default is stdout only.')
    args=parser.parse_args()
    root=args.repo_root.resolve() if args.repo_root else find_root(script)
    base=root/'research/arithmetic'
    rows=[]
    for d,bound,false_value in [(8,1<<27,5**6),(11,1<<25,7**6)]:
        ledger_path=base/('primitive_norm_minus'+str(d)+'_source_ledger.json')
        ledger=json.loads(ledger_path.read_bytes())
        checker_path=base/ledger['checker']
        receipt_path=base/ledger['receipt']
        receipt=json.loads(receipt_path.read_bytes())
        checker_hash=file_digest(checker_path)
        if receipt['checker_sha256']!=checker_hash:
            raise ValueError('Retained receipt is stale against current checker bytes.')
        finite=receipt['finite_coverage' if d==8 else 'finite_core_coverage']
        declared=ledger['complete_scope' if d==8 else 'exact_classification']['finite_bound']
        if finite['bound']!=bound or declared!=bound:
            raise ValueError('Actual certificate bounds differ from the complete checked domain.')
        # Construct sieve flags first: its temporary prime table is released
        # before allocating the independent coordinate array.
        sieve_flags=load_sieve_checker(checker_path,d).core_membership_flags(bound)
        flags,visited,primitive=coordinate_flags(d,bound)
        if len(sieve_flags)!=bound+1 or len(flags)!=bound+1:
            raise ValueError('Membership flags have the wrong byte length.')
        authenticate_flags(sieve_flags,flags)
        rejected=negative_checks(sieve_flags,flags,false_value)
        # Exact whole-bytearray equality above is the verification. Hashes
        # below only identify matching data and detect a stale input receipt.
        expected=finite['core_flags_sha256']
        actual=sha256(flags).hexdigest()
        if actual!=expected:
            raise ValueError('Current membership data differ from retained receipt.')
        rows.append(dict(discriminant=-d,bound_inclusive=bound,ellipse_pairs_visited=visited,
                         primitive_filtered_pairs=primitive,distinct_core_values=sum(flags),bytes=len(flags),
                         comparison='Exact equality of whole independently generated bytearrays.',
                         whole_bytearrays_equal=True,
                         full_membership_sha256=actual,retained_core_flags_sha256=expected,
                         complete_membership_match=True,negative_mutations_rejected=rejected,
                         checker_path=str(checker_path.relative_to(root)),checker_sha256=checker_hash,
                         receipt_path=str(receipt_path.relative_to(root)),receipt_sha256=file_digest(receipt_path),
                         source_ledger_path=str(ledger_path.relative_to(root))))
        del flags,sieve_flags
    # Ledger hashes are deliberately absent: ledgers may pin this receipt,
    # while original checker/receipt bytes above remain explicit pinned inputs.
    result=dict(scope='Independent complete finite membership authentication; no spectral or unrestricted conclusion.',
                method='Independent integral ellipse/gcd enumeration versus the original prime-support sieve; exact whole-bytearray equality, no sumset mask comparison.',
                checker_sha256=file_digest(script),rows=rows)
    out=json.dumps(result,indent=2)+'\n'
    if args.output:
        args.output.write_text(out)
    print(out,end='')


if __name__=='__main__':
    main()
