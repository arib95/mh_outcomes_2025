# External artefact archives

Large release artefacts are hosted on Zenodo record 21269864:

https://zenodo.org/records/21269864

The Git repository tracks these files through git-annex URLs plus portable checksum files.

## Files

| Archive | File | MD5 from Zenodo | Size | URL |
| --- | --- | --- | --- | --- |
| LA5c | LA5c_20260708-230818.tar.gz.partaa | 3f8ab546a0ed112566c0ac77605931e0 | 5.4 GB | https://zenodo.org/records/21269864/files/LA5c_20260708-230818.tar.gz.partaa?download=1 |
| LA5c | LA5c_20260708-230818.tar.gz.partab | f33d1bd6727fcf9ecb5dc3396bea45de | 5.4 GB | https://zenodo.org/records/21269864/files/LA5c_20260708-230818.tar.gz.partab?download=1 |
| LA5c | LA5c_20260708-230818.tar.gz.partac | e5dc426f65d8218d2169775ab93fbae4 | 5.4 GB | https://zenodo.org/records/21269864/files/LA5c_20260708-230818.tar.gz.partac?download=1 |
| LA5c | LA5c_20260708-230818.tar.gz.partad | 8dc6178dda0e7f090fdddcff5e480fc2 | 2.0 GB | https://zenodo.org/records/21269864/files/LA5c_20260708-230818.tar.gz.partad?download=1 |
| TCP | TCP_20260708-230205.tar.gz | a93eb40d463168d695eb2cef55e53a8a | 4.7 GB | https://zenodo.org/records/21269864/files/TCP_20260708-230205.tar.gz?download=1 |

## Restore

With git-annex installed, retrieve the hosted files:

```bash
git annex get Artefacts/LA5c_20260708-230818.tar.gz.partaa \
  Artefacts/LA5c_20260708-230818.tar.gz.partab \
  Artefacts/LA5c_20260708-230818.tar.gz.partac \
  Artefacts/LA5c_20260708-230818.tar.gz.partad \
  Artefacts/TCP_20260708-230205.tar.gz
```

Reassemble LA5c and verify both archives:

```bash
cd Artefacts
cat LA5c_20260708-230818.tar.gz.partaa \
  LA5c_20260708-230818.tar.gz.partab \
  LA5c_20260708-230818.tar.gz.partac \
  LA5c_20260708-230818.tar.gz.partad \
  > LA5c_20260708-230818.tar.gz
shasum -a 256 -c LA5c_20260708-230818.tar.gz.sha256
shasum -a 256 -c TCP_20260708-230205.tar.gz.sha256
```
