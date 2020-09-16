
REM create directory
mkdir mod%1% 2>NUL
cd mod%1%
mkdir trial%2% 2>NUL
cd trial%2%
mkdir run%3% 2>NUL
cd run%3%

REM copy input and output files
copy ..\..\debug.out .\debug.out
copy ..\..\..\errors.csv .\errors.csv
copy ..\..\..\errors.dat .\errors.dat
copy ..\..\..\extract.csv .\extract.csv
copy ..\..\..\extract.dat .\extract.dat
copy ..\..\..\forces.csv .\forces.csv
copy ..\..\..\forces.dat .\forces.dat
copy ..\..\..\head.dat .\head.dat
copy ..\..\..\headerr.dat .\headerr.dat
copy ..\..\..\phreatic.out .\phreatic.out
copy ..\..\..\progress.out .\progress.out
copy ..\..\..\progress.tmp .\progress.tmp
copy ..\..\..\solution .\solution
copy ..\..\..\split.dat .\split.dat
copy ..\..\..\test.csv .\test.csv

REM save PostProcess args
ECHO "rank  = %1%"  > PostProcess.out
ECHO "trial = %2%" >> PostProcess.out
ECHO "run   = %3%" >> PostProcess.out
ECHO "type  = %4%" >> PostProcess.out
ECHO ************* >> PostProcess.out

REM return to workdir
cd ..\..\..
