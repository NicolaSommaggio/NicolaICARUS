#!/bin/bash

ffiles=(
  1d
  2d
  2d_DNN_NOPulseTrains
  2d_DNN_PulseTrains
  1d_YZ
  2d_YZ
  2d_DNN_NOPulseTrains_YZ
  2d_DNN_PulseTrains_YZ
)

particles=(
  muon
#  proton
#  pion
)

ROOT_MACRO="confronto1d2d.C"

cp "$ROOT_MACRO" "${ROOT_MACRO}.bak"

for i in "${!ffiles[@]}"; do
  file="${ffiles[i]}"
  for j in "${!particles[@]}"; do
  particle="${particles[j]}"

  echo -e "\033[1;32mEseguo particle='$particle' file='$file'\033[0m"

root -l << EOF
.L $ROOT_MACRO
confrontoDatiMClight("$particle", "$file");
.q
EOF

  if [ $? -ne 0 ]; then 
    echo -e "\033[1;31mErrore ROOT per particle='$particle' file='$file'\033[0m" 
  fi
done
done
