@echo off

copy ostIn_Serial.txt ostIn.txt

REM Replace C:\bin with location of OSTRICH installation
DEL OstErrors*.txt
REM DEL OstModel*.txt
DEL OstOutput*.txt
C:\bin\Ostrich.exe
pause
