@echo off

copy ostIn_Parallel.txt ostIn.txt

REM replace C:\bin with location of OSTRICH installation
mpiexec -n 4 C:\bin\OstrichMPI.exe
pause
