#!/bin/bash

set -e

SRC_BASE="/pnfs/icarus/scratch/users/gputnam/Ar23+_iterE/ICARUSSpringMC"
DEST="/exp/icarus/data/users/nsommagg/SELECTION_FOLDER"

IDS=(
157 99 58 17 117 132 163 11 5 104 167 26 123 86 153 80 168 19 128 109
63 124 62 98 15 82 16 25 122 121 140 14 48 103 57 149 96 159 36 43
12 55 170 88 34 40 39 127 72 101 111 41 73 44 143 91 29 107 158 13
45 112 102 1 146 94 6 42 0 97 165 93 18 77 129 115 56 110 27 145
31 71 162 70 147
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
