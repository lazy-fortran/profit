#!/bin/sh
set -eu
: "${FORTSYM_BUILD:?Set FORTSYM_BUILD to a configured and built FortSym tree}"
network_sym_root=$(dirname "$FORTSYM_BUILD")
if [ "$(git -C "$network_sym_root" rev-parse HEAD)" != 063008070086aaf2d78f6e9d9a8fddff44c74a46 ]; then
    echo 'Use FortSym revision 063008070086aaf2d78f6e9d9a8fddff44c74a46' >&2
    exit 1
fi
network_task_dir=$(mktemp -d)
trap 'rm -rf "$network_task_dir"' EXIT
${FC:-gfortran} -ffree-line-length-none -I"$FORTSYM_BUILD/include" \
    gen/network/src/generate_network.f90 -L"$FORTSYM_BUILD/lib" -lfortsym \
    -Wl,-rpath,"$FORTSYM_BUILD/lib" -o "$network_task_dir/generate-network"
"$network_task_dir/generate-network" src/generated
