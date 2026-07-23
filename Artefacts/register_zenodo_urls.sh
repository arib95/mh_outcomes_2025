#!/usr/bin/env bash
set -euo pipefail

export COPYFILE_DISABLE=1

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "$script_dir/.." && pwd)"
url_manifest="$script_dir/zenodo_urls.tsv"

cd "$repo_root"

if [[ -z "$(git config --get annex.uuid || true)" ]]; then
  git annex init "zenodo-artefacts" || true
fi

[[ -n "$(git config --get annex.uuid || true)" ]] || {
  echo "git-annex initialization did not establish a repository UUID." >&2
  exit 1
}

registered=0
while IFS=$'\t' read -r file url; do
  [[ "$file" == "file" || -z "$file" ]] && continue
  [[ -n "$url" ]] || {
    echo "Missing URL for $file in $url_manifest" >&2
    exit 1
  }
  case "$file" in
    LA5c_20260722-130855.tar.gz.part??|TCP_20260722-130902.tar.gz.part??) ;;
    *)
      echo "Unexpected filename in URL manifest: $file" >&2
      exit 1
      ;;
  esac

  path="Artefacts/$file"
  [[ -f "$path" ]] || {
    echo "Missing local upload part: $path" >&2
    exit 1
  }

  key=""
  for annex_attempt in 1 2 3; do
    key="$(git annex lookupkey "$path" 2>/dev/null || true)"
    [[ -n "$key" ]] && break
    # On macOS/exFAT, an AppleDouble journal error can occur either before
    # or after the useful operation. A retry consumes the transient entry.
    git annex add --force -- "$path" || true
  done
  [[ -n "$key" ]] || key="$(git annex lookupkey "$path" 2>/dev/null || true)"
  [[ -n "$key" ]] || {
    echo "Could not annex: $path" >&2
    exit 1
  }

  for register_attempt in 1 2 3; do
    git grep -F -- "$url" git-annex -- '*.log.web' >/dev/null 2>&1 && break
    git annex registerurl "$key" "$url" || true
  done
  git grep -F -- "$url" git-annex -- '*.log.web' >/dev/null 2>&1 || {
    echo "Zenodo URL was not recorded for $path" >&2
    exit 1
  }

  echo "registered: $path"
  registered=$((registered + 1))
done < "$url_manifest"

echo "Registered and verified $registered Zenodo URLs with git-annex."
