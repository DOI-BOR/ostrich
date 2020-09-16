@echo off

copy ostIn_FParallel.txt ostIn.txt

REM initialize counter used by SaveBest.bat
echo 0 > Counter.txt

REM replace C:\bin with location of OSTRICH installation
START C:\bin\fmpirun.exe -n 5 -t 24h C:\bin\OstrichFMPI.exe

REM pause
