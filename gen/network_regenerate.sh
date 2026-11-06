#!/bin/sh
set -eu
: "${FORTSYM_BUILD:?Set FORTSYM_BUILD to a configured and built FortSym tree}"
network_sym_root=$(dirname "$FORTSYM_BUILD")
if [ "$(git -C "$network_sym_root" rev-parse HEAD)" != b9ea95c5a0475b016b1e149875a9bec9d8562533 ]; then
    echo 'Use FortSym revision b9ea95c5a0475b016b1e149875a9bec9d8562533' >&2
    exit 1
fi
network_task_dir=$(mktemp -d)
trap 'rm -rf "$network_task_dir"' EXIT
${FC:-gfortran} -ffree-line-length-none -I"$FORTSYM_BUILD/include" \
    gen/network/src/generate_network.f90 -L"$FORTSYM_BUILD/lib" -lfortsym \
    -Wl,-rpath,"$FORTSYM_BUILD/lib" -o "$network_task_dir/generate-network"
"$network_task_dir/generate-network" src/generated
