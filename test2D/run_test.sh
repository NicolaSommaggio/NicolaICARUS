#!/bin/bash

fdatas=(
  msotgia_v10_06_00_07_BNB_1d_respun2_caf_numu
  msotgia_v10_06_00_07_BNB_2d_respun2_caf_numu
  msotgia_v10_06_00_07_BNB_stage1_noPulseTrains_numu_caf
  msotgia_v10_06_00_07_BNB_stage1_pulseTrains_numu_caf
  msotgia_v10_06_00_07_BNB_numu_caf
  msotgia_v10_06_00_07_BNB_yzsim_wc_numu_caf
  msotgia_v10_06_00_07_BNB_yzsim_wcdnn_noPulseTrains_numu_caf
  msotgia_v10_06_00_07_BNB_yzsim_wcdnn_pulseTrains_numu_caf
)

ffiles=(
  /exp/icarus/data/users/nsommagg/test1d.root
  /exp/icarus/data/users/nsommagg/test2d.root
  /exp/icarus/data/users/nsommagg/test2d_DNN_NOPulseTrains.root
  /exp/icarus/data/users/nsommagg/test2d_DNN_PulseTrains.root
  /exp/icarus/data/users/nsommagg/test1d_YZ.root
  /exp/icarus/data/users/nsommagg/test2d_YZ.root
  /exp/icarus/data/users/nsommagg/test2d_DNN_NOPulseTrains_YZ.root
  /exp/icarus/data/users/nsommagg/test2d_DNN_PulseTrains_YZ.root
)

ROOT_MACRO="MacroDataLoader.C"

cp "$ROOT_MACRO" "${ROOT_MACRO}.bak"

for i in "${!fdatas[@]}"; do
  fdata="${fdatas[i]}"
  outfile="${ffiles[i]}"

  sed -i \
  -e "s|const std::string fdata = .*;|const std::string fdata = \"$fdata\";|" \
  -e "s|filename = .*;|filename = \"$outfile\";|" \
  "$ROOT_MACRO"

  #logfile="MacroDataLoader_${i}.log"
  #echo "Eseguo index=$i fdata='$fdata' outfile='$outfile'"

  echo -e "\033[1;32mEseguo index=$i fdata='$fdata' outfile='$outfile'\033[0m"

  #cafe -b -q "$ROOT_MACRO" > "$logfile" 2>&1 || \
  #echo "Errore ROOT per index $i (vedi $logfile)" >&2

  cafe -b -q "$ROOT_MACRO" || \
  echo "Errore ROOT per index $i" >&2

done
