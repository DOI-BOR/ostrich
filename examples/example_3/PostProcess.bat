
REM create directory
mkdir proc%1% 2>NUL
cd proc%1%
mkdir run%3% 2>NUL
cd run%3%

REM copy input and output files
copy ..\..\DCS2.in .\DCS2.in
copy ..\..\DCS2.out .\DCS2.out

REM save PostProcess args
ECHO "rank  = %1%"  > PostProcess.out
ECHO "trial = %2%" >> PostProcess.out
ECHO "run   = %3%" >> PostProcess.out
ECHO "type  = %4%" >> PostProcess.out
ECHO ************* >> PostProcess.out

REM return to workdir
cd ..\..



