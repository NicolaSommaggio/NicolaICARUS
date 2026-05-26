#!/bin/bash

set -e

SRC_BASE="/pnfs/icarus/scratch/users/gputnam/Ar23+_iterE/ICARUSSpringMC"
DEST="/exp/icarus/data/users/nsommagg/TRAINING_FOLDER"

IDS=(
95
169
152
30
8
113
119
160
84
10
78
75
54
89
166
76
150
85
164
154
21
28
156
138
7
79
4
20
130
83
125
141
37
61
120
135
134
23
161
151
105
64
142
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
