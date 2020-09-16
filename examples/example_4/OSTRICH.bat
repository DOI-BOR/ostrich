@echo off

copy ostIn_Serial.txt ostIn.txt

REM initialize counter used by SaveBest.bat
echo 0 > Counter.txt

REM Replace C:\bin with location of OSTRICH installation
H:\bin\Ostrich.exe

pause
