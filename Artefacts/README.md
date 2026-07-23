# External artefact archives

The LA5c and TCP reproducibility archives are hosted on Zenodo:

- Record: https://zenodo.org/records/21269864
- DOI: https://doi.org/10.5281/zenodo.21269864

The Git repository tracks the split archive parts with git-annex. Portable per-part and reassembled-archive checksums are retained alongside them.

## Hosted files

| Archive | File | Size | Zenodo MD5 |
| --- | --- | ---: | --- |
| LA5c | `LA5c_20260722-130855.tar.gz.partaa` | 5.369 GB | `4fd5615b72fada13314fd58c5d97e88a` |
| LA5c | `LA5c_20260722-130855.tar.gz.partab` | 5.369 GB | `ca3b2217d03c4e5795e461b469b095fe` |
| LA5c | `LA5c_20260722-130855.tar.gz.partac` | 5.369 GB | `31cae6cfaa9e8364850b959bc7eb2584` |
| LA5c | `LA5c_20260722-130855.tar.gz.partad` | 5.369 GB | `98dc692133e39be0b28e01e51a64b8df` |
| LA5c | `LA5c_20260722-130855.tar.gz.partae` | 3.668 GB | `46f8e22574e22d64b1eda2cff3b36aa7` |
| TCP | `TCP_20260722-130902.tar.gz.partaa` | 5.369 GB | `c82368db248d5487056b116cd18cf5a6` |
| TCP | `TCP_20260722-130902.tar.gz.partab` | 5.369 GB | `9fae3f0be1f2f5cad94096e3bbb00fba` |
| TCP | `TCP_20260722-130902.tar.gz.partac` | 5.369 GB | `89e9f755134aa31e3aa749b203f41885` |
| TCP | `TCP_20260722-130902.tar.gz.partad` | 5.369 GB | `59c5029b039004342c56d44aaaec5999` |
| TCP | `TCP_20260722-130902.tar.gz.partae` | 0.158 GB | `448d4f0736d56570af28be7530d76f08` |

Permanent file URLs and both MD5 and SHA-256 values are recorded in `zenodo_upload_manifest.tsv`.

## Restore with git-annex

From the repository root:

```bash
git annex get 'Artefacts/LA5c_20260722-130855.tar.gz.part??' \
  'Artefacts/TCP_20260722-130902.tar.gz.part??'

cd Artefacts
cat LA5c_20260722-130855.tar.gz.part?? > LA5c_20260722-130855.tar.gz
cat TCP_20260722-130902.tar.gz.part?? > TCP_20260722-130902.tar.gz
shasum -a 256 -c LA5c_20260722-130855.tar.gz.sha256
shasum -a 256 -c TCP_20260722-130902.tar.gz.sha256
```

## Restore without git-annex

From `Artefacts/`, run:

```bash
./fetch_from_zenodo.sh
```

The script resumes partial downloads, verifies every part against Zenodo's MD5, reassembles both archives, and verifies their SHA-256 checksums.
