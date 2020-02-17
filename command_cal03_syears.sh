job=cal03_syears
rm ${job}.e
rm ${job}.o
bsub -W 168:00 -n 20 -q long -R span[hosts=1] -R rusage[mem=4000] -e ${job}.e -o ${job}.o -J ${job} -cwd "/project/uma_colin_gleason/yuta/DA/missouri/pyHRR/" /project/uma_colin_gleason/yuta/DA/missouri/pyHRR/batch_${job}.sh
