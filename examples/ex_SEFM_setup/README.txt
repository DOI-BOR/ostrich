The following instructions describe setup of OSTRICH for SEFM on Windows.

Place configuration files in directory that contains your SEFM project folder and SEFM code folder. Configuration files include:
	-ostin.txt
	-run_sefm.bat
	-CommandLineParmInput.csv
	-CommandLineParmInput_temp.csv
	-SEFM_wrapper.py and functions.py
	-obj_in.txt
	-ostrich.exe
	-execute.bat

The following updates need to be made to the configuration files to run OSTRICH:

-ostin.txt:
instructions for updating ostin.txt are provided in OSTRICH documentation found on the github repo. Parameter codes in the
Parameter Specification section must match a parameter code in CommandLineParmInput_temp.csv

-run_sefm.bat:
update path on line 7 to point to python installation

-CommandLineParmInput.csv:
Update rows starting at row 36 to include storm template information

-CommandLineParmInput_temp.csv:
Update to indicate which parameters should be calibrated. Parameters with parameter code inside dollar signs (i.e. $max_infil1$)
will be calibrated during OSTRICH simulation. Parameters with number values are constant. Parameter code must match those
specified in ostin.txt

-SEFM_wrapper.py and functions.py:
No updates required

-obj_in.txt:
Update line 1 to indicate event names
Update Line 2 to indicate reservoir name. The name must match the reservoir name in SEFM
Update line 3 to include objective functions of interest. options include NSE, KGE, pkerr, and PBIAS
Update line 4 to determine number of decimals for rounding
Update line 5 to show folder name of SEFM project

-execute.bat:
Update line 5 for correct SEFM model folder name and SEFM project folder name and project file name




