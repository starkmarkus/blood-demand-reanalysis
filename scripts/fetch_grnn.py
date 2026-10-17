"""Fetch the compatible historical GRNN implementation, preserving its source."""
from pathlib import Path
import urllib.request,tarfile,io,hashlib,json
root=Path(__file__).resolve().parents[1]
meta=json.loads((root/'inputs/grnn-package.json').read_text())
with urllib.request.urlopen(meta['url'],timeout=60) as response:data=response.read()
assert hashlib.sha256(data).hexdigest()==meta['sha256'],'Package hash mismatch'
out=root/'inputs/grnn'
out.mkdir(exist_ok=True)
with tarfile.open(fileobj=io.BytesIO(data),mode='r:gz') as tar:
 for member in tar.getmembers():
  if member.isfile() and (member.name.startswith('tsfgrnn/R/') or member.name=='tsfgrnn/DESCRIPTION'):
   target=out/member.name
   assert target.resolve().is_relative_to(out.resolve())
   target.parent.mkdir(parents=True,exist_ok=True)
   target.write_bytes(tar.extractfile(member).read())
print('Verified tsfgrnn',meta['version'],'source files.')
