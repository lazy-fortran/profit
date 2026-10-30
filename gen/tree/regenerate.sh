#!/bin/sh
set -eu
: "${FORTSYM_BUILD:?Set FORTSYM_BUILD to a configured and built FortSym tree}"
FORTSYM_ROOT=$(dirname "$FORTSYM_BUILD")
if [ "$(git -C "$FORTSYM_ROOT" rev-parse HEAD)" != 063008070086aaf2d78f6e9d9a8fddff44c74a46 ]; then
    echo 'Use FortSym revision 063008070086aaf2d78f6e9d9a8fddff44c74a46' >&2
    exit 1
fi
task_dir=$(mktemp -d)
trap 'rm -rf "$task_dir"' EXIT
${FC:-gfortran} -ffree-line-length-none -I"$FORTSYM_BUILD/include" \
    gen/tree/generate_tree.f90 -L"$FORTSYM_BUILD/lib" -lfortsym \
    -Wl,-rpath,"$FORTSYM_BUILD/lib" -o "$task_dir/generate-tree"
"$task_dir/generate-tree" src/generated
