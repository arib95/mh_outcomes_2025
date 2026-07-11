#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$script_dir"

fetch() {
  local url="$1"
  local out="$2"
  if [[ -s "$out" ]]; then
    echo "exists: $out"
  else
    echo "downloading: $out"
    curl -L --fail --continue-at - --output "$out" "$url"
  fi
}

md5_check() {
  local file="$1"
  local expected="$2"
  local observed
  if command -v md5sum >/dev/null 2>&1; then
    observed="$(md5sum "$file" | awk '{print $1}')"
  else
    observed="$(md5 -q "$file")"
  fi
  if [[ "$observed" != "$expected" ]]; then
    echo "MD5 mismatch for $file" >&2
    echo "expected: $expected" >&2
    echo "observed: $observed" >&2
    exit 1
  fi
}

sha256_check() {
  local manifest="$1"
  if command -v sha256sum >/dev/null 2>&1; then
    sha256sum -c "$manifest"
  else
    shasum -a 256 -c "$manifest"
  fi
}

fetch 'https://zenodo.org/records/21269864/files/LA5c_20260708-230818.tar.gz.partaa?download=1' 'LA5c_20260708-230818.tar.gz.partaa'
fetch 'https://zenodo.org/records/21269864/files/LA5c_20260708-230818.tar.gz.partab?download=1' 'LA5c_20260708-230818.tar.gz.partab'
fetch 'https://zenodo.org/records/21269864/files/LA5c_20260708-230818.tar.gz.partac?download=1' 'LA5c_20260708-230818.tar.gz.partac'
fetch 'https://zenodo.org/records/21269864/files/LA5c_20260708-230818.tar.gz.partad?download=1' 'LA5c_20260708-230818.tar.gz.partad'
fetch 'https://zenodo.org/records/21269864/files/TCP_20260708-230205.tar.gz?download=1' 'TCP_20260708-230205.tar.gz'

md5_check 'LA5c_20260708-230818.tar.gz.partaa' '3f8ab546a0ed112566c0ac77605931e0'
md5_check 'LA5c_20260708-230818.tar.gz.partab' 'f33d1bd6727fcf9ecb5dc3396bea45de'
md5_check 'LA5c_20260708-230818.tar.gz.partac' 'e5dc426f65d8218d2169775ab93fbae4'
md5_check 'LA5c_20260708-230818.tar.gz.partad' '8dc6178dda0e7f090fdddcff5e480fc2'
md5_check 'TCP_20260708-230205.tar.gz' 'a93eb40d463168d695eb2cef55e53a8a'

cat \
  LA5c_20260708-230818.tar.gz.partaa \
  LA5c_20260708-230818.tar.gz.partab \
  LA5c_20260708-230818.tar.gz.partac \
  LA5c_20260708-230818.tar.gz.partad \
  > LA5c_20260708-230818.tar.gz

sha256_check 'LA5c_20260708-230818.tar.gz.sha256'
sha256_check 'TCP_20260708-230205.tar.gz.sha256'
