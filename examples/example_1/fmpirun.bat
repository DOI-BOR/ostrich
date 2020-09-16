@echo off
REM replace C:\bin with location of OSTRICH installation
START C:\bin\fmpirun.exe -n 4 -t 24h C:\bin\OstrichFMPI.exe
REM pause
