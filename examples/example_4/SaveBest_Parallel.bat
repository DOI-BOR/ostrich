@ECHO OFF

REM set a fixed path for the "best" folder
SET ROOT=H:\Matott\Software\Ostrich\WebReleases\v17.10.30\Windows\Demos\Demo4\test\parallel

REM read in the counter
SET /P A=<%ROOT%\Counter.txt

REM increment the counter
SET /A A=%A%+1

REM create the best directory
mkdir %ROOT%\best 2>NUL

REM Create a subdir for the latest best solution.
REM This is optional, could just overwrite previous best
REM if only interest in saving files associated with the
REM ulitmate best objective function.
mkdir %ROOT%\best\%A% 

REM copy model files to new subdir
copy * %ROOT%\best\%A% 
del %ROOT%\best\%A%\*.exe

REM save counter to file for next time
echo %A% > %ROOT%\Counter.txt
