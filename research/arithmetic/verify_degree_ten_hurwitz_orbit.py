#!/usr/bin/env python3
"""Complete degree-ten double-transposition Hurwitz-orbit certificate.

Default replay compiles the embedded exact C++17 algorithm with assertions
enabled, runs its finite enumeration and all Hurwitz edges, and checks the
complete counts. No sampling or group-classification oracle is used.
Requires g++ and about 0.8 GB RAM. It writes temporary build files only under
the ignored repository .scratch directory.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import platform
import re
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / "research/meridional-rank/representations/alternating_target_prime_families_inputs.json"
CPP_SOURCE = r"""#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <unordered_map>
#include <vector>
using U=uint64_t;
using Clock=std::chrono::steady_clock;
const U VIS=U(1)<<63;
U compose(U a,U b){U c=0;for(unsigned j=0;j<10;j++)c|=((a>>(4*((b>>(4*j))&15)))&15)<<(4*j);return c;}
U inv(U a){U b=0;for(unsigned j=0;j<10;j++)b|=U(j)<<(4*((a>>(4*j))&15));return b;}
U dt(unsigned a,unsigned b,unsigned c,unsigned d){
 U p=0;for(unsigned j=0;j<10;j++)p|=U(j)<<(4*j);
 p^=U(a^b)<<(4*a);p^=U(a^b)<<(4*b);
 p^=U(c^d)<<(4*c);p^=U(c^d)<<(4*d);return p;
}
U mix(U a){a^=a>>30;a*=0xbf58476d1ce4e5b9ULL;a^=a>>27;a*=0x94d049bb133111ebULL;return a^(a>>31);}
int main(){
 auto start=Clock::now();
 std::vector<U>C;
 for(unsigned a=0;a<10;a++)for(unsigned b=a+1;b<10;b++)for(unsigned c=b+1;c<10;c++)for(unsigned d=c+1;d<10;d++){
  C.push_back(dt(a,b,c,d));C.push_back(dt(a,c,b,d));C.push_back(dt(a,d,b,c));
 }
 std::sort(C.begin(),C.end());assert(C.size()==630&&std::unique(C.begin(),C.end())==C.end());
 std::unordered_map<U,unsigned>idx;
 for(unsigned i=0;i<C.size();i++)idx[C[i]]=i;
 std::array<U,5>seedp={dt(0,1,2,3),dt(1,4,2,3),dt(0,5,1,6),dt(1,2,6,7),dt(1,8,6,9)};
 U identity=0;for(unsigned j=0;j<10;j++)identity|=U(j)<<(4*j);
 U P=identity,seed=0;
 for(unsigned j=0;j<5;j++){P=compose(P,seedp[j]);seed|=U(idx[seedp[j]])<<(10*j);}
 assert(((P>>12)&15)==3);
 std::unordered_map<U,std::vector<unsigned>>pairs;
 for(unsigned a=0;a<630;a++)for(unsigned b=0;b<630;b++)pairs[compose(C[a],C[b])].push_back(a|(b<<10));
 for(U element:C){
  assert(compose(element,element)==identity);
  unsigned support=0;for(unsigned j=0;j<10;j++)support+=((element>>(4*j))&15)!=j;
  assert(support==4);
 }
 std::vector<unsigned char> moves(630);for(unsigned i=0;i<630;i++)moves[i]=((C[i]>>12)&15)!=3;
 constexpr size_t CAP=size_t(1)<<26;constexpr size_t MASK=CAP-1;
 std::vector<U>table(CAP,0);
 U total=0,full=0;
 auto insert=[&](U key){U stored=key+1;size_t s=mix(key)&MASK;while(table[s]&&table[s]!=stored)s=(s+1)&MASK;assert(!table[s]);table[s]=stored;};
 for(auto const& item:pairs){
  U target=compose(inv(item.first),P);
  for(unsigned c=0;c<630;c++){
   auto it=pairs.find(compose(C[c],target));if(it==pairs.end())continue;
   for(unsigned ab:item.second){
    unsigned a=ab&1023,b=ab>>10;bool firstmoves=moves[a]||moves[b]||moves[c];
    for(unsigned de:it->second){
     ++total;unsigned d=de&1023,e=de>>10;
     if(firstmoves||moves[d]||moves[e]){++full;insert(U(a)|(U(b)<<10)|(U(c)<<20)|(U(d)<<30)|(U(e)<<40));}
    }
   }
  }
 }
 assert(total==48892572&&full==17006112);
 std::cout<<"enumeration "<<total<<" full "<<full<<" pairproducts "<<pairs.size()<<" seconds "<<std::chrono::duration<double>(Clock::now()-start).count()<<std::endl;
 std::vector<uint16_t>conj(630*630);
 for(unsigned a=0;a<630;a++)for(unsigned b=0;b<630;b++)conj[a*630+b]=idx.at(compose(compose(C[a],C[b]),C[a]));
 auto visit=[&](U key){size_t s=mix(key)&MASK;U stored=key+1;while(table[s]&&((table[s]&~VIS)!=stored))s=(s+1)&MASK;assert(table[s]);if(table[s]&VIS)return false;table[s]|=VIS;return true;};
 // Invalid-product controls must not occur in the retained domain.
 auto contains=[&](U key){size_t s=mix(key)&MASK;U stored=key+1;while(table[s]&&((table[s]&~VIS)!=stored))s=(s+1)&MASK;return table[s]!=0;};
 U wrong=(seed&~U(1023))|U(idx[seedp[1]]);
 assert(!contains(wrong));
 U swapped=(seed&~U(1048575))|U(idx[seedp[1]])|(U(idx[seedp[0]])<<10);
 assert(!contains(swapped));
 U seen=0;unsigned compcount=0;U initial=seed;std::vector<U>queue;queue.reserve(full);
 while(seen<full){
  if(compcount){initial=0;for(U slot:table)if(slot&&!(slot&VIS)){initial=slot-1;break;}}
  assert(visit(initial));queue.clear();queue.push_back(initial);
  for(size_t h=0;h<queue.size();h++){
   U key=queue[h];
   for(unsigned j=0;j<4;j++){
    unsigned a=(key>>(10*j))&1023,b=(key>>(10*(j+1)))&1023;
    U erased=key&~(U(1048575)<<(10*j));
    U pos=erased|(U(b)<<(10*j))|(U(conj[b*630+a])<<(10*(j+1)));
    U neg=erased|(U(conj[a*630+b])<<(10*j))|(U(a)<<(10*(j+1)));
    if(visit(pos))queue.push_back(pos);if(visit(neg))queue.push_back(neg);
   }
  }
  seen+=queue.size();++compcount;
  std::cout<<"component "<<compcount<<" size "<<queue.size()<<" seen "<<seen<<" seconds "<<std::chrono::duration<double>(Clock::now()-start).count()<<std::endl;
 }
 assert(compcount==1&&seen==full);
 std::cout<<"DONE "<<compcount<<" components"<<std::endl;
}
"""

def require(ok, message):
    if not ok:
        raise AssertionError(message)

def permutation(*cycles):
    out=list(range(10))
    for cycle in cycles:
        for a,b in zip(cycle,cycle[1:]+cycle[:1]):
            out[a]=b
    return tuple(out)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-receipt",type=Path)
    parser.add_argument("--cxx",default="g++")
    args=parser.parse_args()
    data=json.loads(INPUT.read_text())
    seed=(permutation((0,1),(2,3)),permutation((1,4),(2,3)),
          permutation((0,5),(1,6)),permutation((1,2),(6,7)),
          permutation((1,8),(6,9)))
    require([list(v[:8]) for v in seed[:4]]==data["seed_quartet"],
            "retained eight-sheet seed binding")
    product=tuple(range(10))
    for value in seed:
        product=tuple(product[value[j]] for j in range(10))
    require(product==permutation((0,5,1,8,2,6,9,7,4)),
            "literal ten-sheet product")
    version=subprocess.run([args.cxx,"--version"],capture_output=True,text=True,check=True).stdout.splitlines()[0]
    scratch=ROOT/".scratch"
    scratch.mkdir(exist_ok=True)
    flags=["-O3","-std=c++17","-UNDEBUG"]
    with tempfile.TemporaryDirectory(prefix="hurwitz-degree10-",dir=scratch) as build:
        build=Path(build)
        source=build/"orbit.cpp"
        binary=build/"orbit"
        source.write_text(CPP_SOURCE)
        subprocess.run([args.cxx,*flags,str(source),"-o",str(binary)],
                       capture_output=True,text=True,check=True)
        run=subprocess.run([str(binary)],capture_output=True,text=True,check=True)
    enumeration=re.search(r"enumeration (\d+) full (\d+) pairproducts (\d+)",run.stdout)
    require(enumeration is not None,"enumeration report")
    raw,full,pairs=map(int,enumeration.groups())
    components=[int(x) for x in re.findall(r"component \d+ size (\d+)",run.stdout)]
    require((raw,full,pairs)==(48892572,17006112,64144),"exact domain census")
    require(components==[17006112],"single full orbit")
    require("DONE 1 components" in run.stdout,"complete traversal")
    require(full%9==0 and full//9==1889568,"absolute quotient")
    report={
        "schema_version":1,
        "kind":"complete_finite_Hurwitz_orbit_certificate",
        "python":platform.python_version(),
        "cxx":version,
        "compiler_flags":flags,
        "arithmetic":"exact 64-bit packed permutations; no floating-point proof steps",
        "checker_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "embedded_cpp_sha256":hashlib.sha256(CPP_SOURCE.encode()).hexdigest(),
        "seed_input_sha256":hashlib.sha256(INPUT.read_bytes()).hexdigest(),
        "double_transposition_class_size":630,
        "pair_product_count":pairs,
        "raw_fixed_product_tuple_count":raw,
        "nontransitive_tuple_count":raw-full,
        "transitive_tuple_count":full,
        "B5_orbit_sizes":components,
        "absolute_Nielsen_class_size":full//9,
        "seed_tuple":[list(v) for v in seed],
        "target_product":list(product),
        "controls":["all class elements are double transpositions",
                    "duplicate tuple insertion rejected",
                    "wrong-product tuple excluded",
                    "unconjugated adjacent swap excluded",
                    "every Hurwitz neighbor belongs to the exhaustive domain",
                    "both positive and negative moves traversed"],
        "conclusion_with_seed_generation_proof":"all 17006112 transitive tuples form one A10 orbit",
        "scope":"Exact finite n=10 classification only. Known generation of the literal seed and preservation of generated subgroup identify the unique orbit as full A10. Does not prove uniform braid transitivity, unirationality or a rational-point theorem in higher degree.",
        "independent_replay":False,
        "formal_verification":False,
        "external_review":False,
        "result":"PASS",
    }
    output=json.dumps(report,indent=2,sort_keys=True)+"\n"
    if args.write_receipt:
        args.write_receipt.write_text(output)
    print(output,end="")

if __name__=="__main__":
    main()
