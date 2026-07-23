#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
manifest="$script_dir/zenodo_upload_manifest.tsv"
cd "$script_dir"

[[ -f "$manifest" ]] || {
  echo "Missing manifest: $manifest" >&2
  exit 1
}

file_size() {
  stat -f '%z' "$1" 2>/dev/null || stat -c '%s' "$1"
}

md5_of() {
  if command -v md5sum >/dev/null 2>&1; then
    md5sum "$1" | awk '{print $1}'
  else
    md5 -q "$1"
  fi
}

fetch_and_verify() {
  local url="$1"
  local file="$2"
  local expected_bytes="$3"
  local expected_md5="$4"
  local actual_bytes observed_md5

  if [[ -f "$file" ]]; then
    actual_bytes="$(file_size "$file")"
    if [[ "$actual_bytes" == "$expected_bytes" ]]; then
      observed_md5="$(md5_of "$file")"
      if [[ "$observed_md5" == "$expected_md5" ]]; then
        echo "verified already: $file"
        return
      fi
      echo "Complete-sized file has the wrong MD5: $file" >&2
      exit 1
    fi
    if (( actual_bytes > expected_bytes )); then
      echo "Local file is larger than expected: $file" >&2
      exit 1
    fi
    echo "resuming: $file"
  else
    echo "downloading: $file"
  fi

  curl --location --fail --show-error \
    --retry 5 --retry-all-errors --retry-delay 5 \
    --continue-at - --output "$file" "$url"

  actual_bytes="$(file_size "$file")"
  [[ "$actual_bytes" == "$expected_bytes" ]] || {
    echo "Size mismatch for $file: expected $expected_bytes, found $actual_bytes" >&2
    exit 1
  }
  observed_md5="$(md5_of "$file")"
  [[ "$observed_md5" == "$expected_md5" ]] || {
    echo "MD5 mismatch for $file" >&2
    exit 1
  }
  echo "verified: $file"
}

while IFS=$'\t' read -r archive file bytes size_gb md5 sha256 url; do
  [[ "$file" == "file" || -z "$file" ]] && continue
  [[ -n "$url" ]] || {
    echo "Missing Zenodo URL for $file" >&2
    exit 1
  }
  case "$file" in
    LA5c_20260722-130855.tar.gz.part??|TCP_20260722-130902.tar.gz.part??) ;;
    *)
      echo "Unexpected filename in manifest: $file" >&2
      exit 1
      ;;
  esac
  fetch_and_verify "$url" "$file" "$bytes" "$md5"
done < "$manifest"

cat LA5c_20260722-130855.tar.gz.part?? > LA5c_20260722-130855.tar.gz
cat TCP_20260722-130902.tar.gz.part?? > TCP_20260722-130902.tar.gz

if command -v sha256sum >/dev/null 2>&1; then
  sha256sum -c LA5c_20260722-130855.tar.gz.sha256
  sha256sum -c TCP_20260722-130902.tar.gz.sha256
else
  shasum -a 256 -c LA5c_20260722-130855.tar.gz.sha256
  shasum -a 256 -c TCP_20260722-130902.tar.gz.sha256
fi
