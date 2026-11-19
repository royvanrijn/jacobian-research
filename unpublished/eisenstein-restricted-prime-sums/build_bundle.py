#!/usr/bin/env python3
"""Build the standalone source archive; optionally verify a clean PDF rebuild."""
from pathlib import Path
from hashlib import sha256
from zipfile import ZipFile, ZipInfo, ZIP_DEFLATED
import argparse,json,subprocess,tempfile

def digest(path):
    return sha256(path.read_bytes()).hexdigest()

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify-rebuild',action='store_true')
    args=parser.parse_args()
    root=Path(__file__).resolve().parent
    files=[root/name for name in ('main.tex','references.tex','README.md','metadata.json','verification.json','build_bundle.py')]
    files+=sorted((root/'sections').glob('*.tex'))
    files+=sorted(path for path in (root/'certificates').rglob('*')
                  if path.is_file() and '__pycache__' not in path.parts
                  and path.suffix in ('.py','.java','.json','.pdf'))
    for path in files:
        if not path.is_file() or not path.resolve().is_relative_to(root):
            raise RuntimeError('Missing or external source: '+str(path))
    target=root/'source.zip'
    with ZipFile(target,'w') as archive:
        for path in sorted(files):
            info=ZipInfo(path.relative_to(root).as_posix(),(1980,1,1,0,0,0))
            info.compress_type=ZIP_DEFLATED
            info.external_attr=0o100644<<16
            archive.writestr(info,path.read_bytes())
    result={'archive':'source.zip','members':len(files),'archive_sha256':digest(target)}
    if args.verify_rebuild:
        with tempfile.TemporaryDirectory(prefix='exact-norms-source-') as temp:
            dest=Path(temp)
            with ZipFile(target) as archive:
                for name in archive.namelist():
                    if not (dest/name).resolve().is_relative_to(dest):
                        raise RuntimeError('Unsafe archive path')
                archive.extractall(dest)
            build=subprocess.run(['latexmk','-pdf','-jobname=paper','-interaction=nonstopmode','-halt-on-error','main.tex'],cwd=dest,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
            if build.returncode:
                raise RuntimeError(build.stdout[-8000:])
            log=(dest/'paper.log').read_text(errors='replace')
            if any(x in log for x in ('undefined references','undefined citations','Overfull','LaTeX Warning:')):
                raise RuntimeError('Clean rebuild diagnostics: '+chr(10).join(line for line in log.splitlines() if any(t in line for t in ('Warning','undefined','Overfull'))))
            original=digest(root/'paper.pdf'); rebuilt=digest(dest/'paper.pdf')
            if original!=rebuilt:
                raise RuntimeError('PDF is not byte-identical: '+original+' != '+rebuilt)
            result['clean_source_rebuild']='passed; byte-identical PDF'
            result['pdf_sha256']=original
    print(json.dumps(result,indent=2))

if __name__=='__main__':
    main()
