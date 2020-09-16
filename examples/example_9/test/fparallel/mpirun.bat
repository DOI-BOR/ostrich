@echo off
REM replace C:\bin with location of OSTRICH installation
copy OstIn_parallel.txt OstIn.txt
mpiexec -n 5 C:\bin\OstrichMPI.exe
pause
