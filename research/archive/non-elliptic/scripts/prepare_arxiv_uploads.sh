#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
paper_dirs=(
  "papers/gaussian-moments-two-variables"
  "papers/sparse-minimality-gaussian-moments-dimension-three"
)

for relative_paper_dir in "${paper_dirs[@]}"; do
  paper_dir="$repo_root/$relative_paper_dir"
  upload_dir="$paper_dir/arxiv-upload"

  latexmk -cd -pdf -interaction=nonstopmode -halt-on-error \
    "$paper_dir/main.tex" >/dev/null

  mkdir -p "$upload_dir"
  unexpected_files="$(
    find "$upload_dir" -mindepth 1 -maxdepth 1 \
      ! -name main.tex ! -name main.bbl -print
  )"
  if [[ -n "$unexpected_files" ]]; then
    printf 'Refusing to overwrite unexpected arXiv upload contents:\n%s\n' \
      "$unexpected_files" >&2
    exit 1
  fi

  cp "$paper_dir/main.tex" "$upload_dir/main.tex"
  cp "$paper_dir/main.bbl" "$upload_dir/main.bbl"
  printf '\n\\typeout{get arXiv to do 4 passes: Label(s) may have changed. Rerun}\n' \
    >>"$upload_dir/main.tex"

  if grep -nE '(^|[^\\])%' "$upload_dir/main.tex"; then
    printf 'Unescaped LaTeX comments remain in %s\n' "$upload_dir/main.tex" >&2
    exit 1
  fi

  validation_dir="$(mktemp -d "$repo_root/tmp/pdfs/arxiv-upload.XXXXXX")"
  cp "$upload_dir/main.tex" "$validation_dir/main.tex"
  cp "$upload_dir/main.bbl" "$validation_dir/main.bbl"
  for _pass in 1 2 3 4; do
    (
      cd "$validation_dir"
      pdflatex -interaction=nonstopmode -halt-on-error main.tex >/dev/null
    )
  done

  if grep -Eq 'undefined (citations|references)|Citation .* undefined|Reference .* undefined' \
    "$validation_dir/main.log"; then
    printf 'Unresolved citations or references in %s\n' "$relative_paper_dir" >&2
    exit 1
  fi

  find "$validation_dir" -depth -delete
  printf 'Prepared and verified %s/arxiv-upload\n' "$relative_paper_dir"
done
