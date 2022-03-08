logFile=timing/log
timingFile=timing/timelog


rm $logFile
rm $timingFile


for iter in {1..20}
do
    ./Pele3d.llvm.ex inputs/inputs.0d_cvode_faster > $logFile
    python timing/getTiming.py $logFile $timingFile 
done

python timing/printElapsedTime.py $timingFile 
