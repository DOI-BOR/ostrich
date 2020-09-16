@echo off
REM replace C:\bin with location of OSTRICH installation
copy OstIn_parallel.txt OstIn.txt
START C:\bin\fmpirun.exe -n 5 -t 24h C:\bin\OstrichFMPI.exe
REM pause
