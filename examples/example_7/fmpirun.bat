@echo off

copy ostIn_parallel.txt ostIn.txt

REM replace C:\bin with location of OSTRICH installation
START C:\bin\fmpirun.exe -n 5 -t 24h C:\bin\OstrichFMPI.exe

REM pause
