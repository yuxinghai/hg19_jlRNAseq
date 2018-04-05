#!/usr/bin/bash
proj_name="JLRNAseq"
proj="JLRNAseq"
echo "track ${proj_name}"
echo "compositeTrack on"
# container multiWig
echo "configurable on"
echo "shortLabel ${proj_name}"
echo "longLabel ${proj_name}"
# aggregate transparentOverlay
echo "dragAndDrop subTracksq"
echo "type bigWig 0 100"
echo "autoScale off"
echo "viewLimits 0:100"
echo "configurable on"
echo "visibility full"
echo "maxHeightPixels 64:64:11"
echo -e "showSubtrackColorOnUi on \n"

NAMES=(F1 F2 G1 G2 N01 N02 NT1 NT2 PS1 PS2 V1 V2)
COLS=(50,205,50 128,128,128 70,130,180 153,50,204 255,0,255 255,69,0)
rep_num=2
reps=(_rep1 _rep2)
declare -a cols=($(for i in ${COLS[@]}; do for ((n=1; n<=$rep_num; n++)) do echo "$i";done; done))
T_names=(PSF PSF Control Control NONO NONO NEAT1_T NEAT1_T PSPC1 PSPC1 NEAT1_V2 NEAT1_V2)
numSamples=${#NAMES[@]}

for (( i=0; i<${numSamples}; i++ )); do
    name=${NAMES[$i]}
    T_name=${T_names[$i]}${reps[($i%2)]}
    col=${cols[$i]}
    #echo forward file
    echo "track ${T_name}_F"
    echo "bigDataUrl ftp://202.114.67.213/yuxinghai/${proj}/data/${name}.Forward.bw"
    echo "shortLabel ${T_name}_F"
    echo "longLabel ${T_name}_F"
    echo "type bigWig"
    echo "autoScale on"
    echo "color ${col}"    
    echo -e "parent ${proj_name} \n"
done



for (( i=0; i<${numSamples}; i++ )); do
    name=${NAMES[$i]}
    T_name=${T_names[$i]}${reps[($i%2)]}
    col=${cols[$i]}
    #echo reverse file
    echo "track ${T_name}_R"
    echo "bigDataUrl ftp://202.114.67.213/yuxinghai/${proj}/data/${name}.Reverse.bw"
    echo "shortLabel ${T_name}_R"
    echo "longLabel ${T_name}_R"
    echo "type bigWig"
    echo "autoScale on"
    echo "color ${col}"    
    echo -e "parent ${proj_name} \n"

done


