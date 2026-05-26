#!/bin/bash

set -e

SRC_BASE="/pnfs/icarus/scratch/users/gputnam/Ar23+_iterE/ICARUSSpringMC"
DEST="/exp/icarus/data/users/nsommagg/PDF_FOLDER"

IDS=(
148 47 90 114 49 67 155 126 51 33 137 106 133 131 46 60 22 68 32 50 9 3 24 87 136 74 2 144 53 35 118 52 38 69 66 59 92 100 81 139 116 65
)

mkdir -p "$DEST"

for id in "${IDS[@]}"; do
    SRC="$SRC_BASE/84005943_${id}/out1.flat.caf.root"

    if [[ -f "$SRC" ]]; then
        RUN=$(basename "$(dirname "$SRC")")
        DEST_FILE="${DEST}/${RUN}_out1.flat.caf.root"

        echo "Copying $SRC -> $DEST_FILE"
        cp "$SRC" "$DEST_FILE"
    else
        echo "WARNING: missing file $SRC"
    fi
done

echo "Done."
