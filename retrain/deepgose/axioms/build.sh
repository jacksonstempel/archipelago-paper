#!/bin/bash
# DeepGO-SE's ontology axioms: GO-Plus of the GO release whose go-basic.obo labels the data, reduced to the OWL 2 EL
# profile (makeEl.groovy) and normalized (Normalizer.groovy).
#
#   bash retrain/deepgose/axioms/build.sh                (data root from $ARCHIPELAGO_DATA)
#
# Input:  $ARCHIPELAGO_DATA/external/deepgose/go-plus.owl, from
#         https://release.geneontology.org/2026-08-05/ontology/extensions/go-plus.owl (the GO release 2026-07-26, the
#         release of datasets/ontology/go-basic.obo).
# Output: $ARCHIPELAGO_DATA/inputs/deepgose/axioms/{go-plus-el.owl, go-plus.norm, normalize.log}; go-plus.norm is the
#         file train.py reads (467,171 statements).
# Every file is checked against the SHA-256 of the paper's run. The normalizer's numbering of auxiliary classes sets the
# size of DeepGO-SE's class table, so a different go-plus.norm is a different model, not a harmless variation.
# Needs Java 11 and Groovy 2.4 (the paper used 2.4.21; Grape downloads the Java libraries on first use), about 80 GB of
# memory (JAVA_OPTS=-Xmx80g); the paper's run fitted in a two-hour job on 4 CPUs.
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
root=${ARCHIPELAGO_DATA:?set ARCHIPELAGO_DATA to the data root}
in=$root/external/deepgose/go-plus.owl
out=$root/inputs/deepgose/axioms

check() {  # file, expected SHA-256, what it is
  local got
  got=$(sha256sum "$1" | cut -d' ' -f1)
  if [[ $got != "$2" ]]; then
    echo "SHA-256 mismatch for $3 ($1): got $got, the paper's run had $2" >&2
    exit 1
  fi
  echo "ok $3 $got"
}

check "$in" 52f51e70a94667f887f89205daeb442abe5f6f8b6cd8dda51ff5ad7fe7524bf5 "GO-Plus (release 2026-07-26)"
mkdir -p "$out"
export JAVA_OPTS=${JAVA_OPTS:--Xmx80g}
groovy "$here/makeEl.groovy" "$in" "$out/go-plus-el.owl.tmp"
mv "$out/go-plus-el.owl.tmp" "$out/go-plus-el.owl"
check "$out/go-plus-el.owl" 1daf207964edfd813bc57c93402b4ad4db437270732792410e3d20b8f9fe954a "GO-Plus, EL profile"
groovy "$here/Normalizer.groovy" -i "$out/go-plus-el.owl" -o "$out/go-plus.norm.tmp" | tee "$out/normalize.log"
mv "$out/go-plus.norm.tmp" "$out/go-plus.norm"
check "$out/go-plus.norm" 13c7b4d30f3abf69ae843737a2c4f6284934dd53c31220e626e6b534fb835e77 "normal forms"
