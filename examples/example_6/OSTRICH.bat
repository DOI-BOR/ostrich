@echo off

copy ostIn_serial.txt ostIn.txt

DEL PostProcess.out

REM Replace C:\bin with location of OSTRICH installation
C:\bin\Ostrich.exe

pause
