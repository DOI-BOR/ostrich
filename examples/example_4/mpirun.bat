@echo off

copy ostIn_Parallel.txt ostIn.txt

REM initialize counter used by SaveBest.bat
echo 0 > Counter.txt

REM replace C:\bin with location of OSTRICH installation
mpiexec -n 5 C:\bin\OstrichMPI.exe

pause
