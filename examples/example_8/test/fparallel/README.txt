# input file must be named MyZonesIn.txt
# output file will be called MyZonesOut.txt

# to run the program use the nZones.bat batch file. 
# If any of the included .exe files are moved to a different
# directory, be sure to modify the batch files accordingly.

# NOTE: the mapping of layer configurations is given at the bottom of this file.

#Annotated input file
9  #configuration of first layer (a value between 1 and 13)
10 #configuration of second layer (a value between 1 and 13)
10 #configuration of third layer (a value between 1 and 13)
10 #configuration of fourth layer (a value between 1 and 13)
14 #configuration of fifth layer (a value between 1 and 14)
14 #configuration of sixth layer (a value between 1 and 14)
1  # type of contaminant (1 = benzene, 2 = TCE, 3 = 1,2-DCB
# if not using pre-emption, leave out the CostThreshold and
# FluxThrehsold lines. 
CostThreshold 100.00  # if enabling pre-emption, specify the cost threshold
# if enabling pre-emption, specify the penalty method (APM or MPM)  followed
# by the penalty multiplier and the performance constraint.
FluxThreshold MPM 1E6 5E-6 

#Annotated output file
Param[01] = 09 #summary of layer 1 configuration
Param[02] = 10 #summary of layer 2 configuration
Param[03] = 10 #summary of layer 3 configuration
Param[04] = 10 #summary of layer 4 configuration
Param[05] = 14 #summary of layer 5 configuration
Param[06] = 14 #summary of layer 6 configuration
Pollutant = Benzene #summary of pollutant type
Cost Threshold : none specified #summary of cost-based pre-emption threshold (if using)
Height : 0.600000 # liner height (0.15 m per active layer)
Complexity : 3 # complexity factor, high values indicate barrier will likely be difficult to construct.
Material Costs # summary of layer costs
Layer : Cost
01 : 2.31     # layer 1 costs
02 : 10.65    # layer 2 costs
03 : 10.65    # layer 3 costs
04 : 10.65    # layer 4 costs
05 : 0.00     # layer 5 costs
06 : 0.00     # layer 6 costs
Total : 34.27 # overall cost of liner (not including performance penalties)
Flux Threshold : none specified #summary of performanc-based pre-emption threshold (if using)
Mass : 1.216289E-005 # amount of mass exiting after 100 years (>= 1.0 if pre-empted)

#######################################################################################################
#######################################################################################################
# Mapping of layer configrations
#######################################################################################################
#######################################################################################################
coded      sorptive      percentage          sand/bentonite
value      material      amendment            composition            cost     
1	     BTEA            3%          87% sand, 10% bentonite     19.33
2	     BTEA            6%          84% sand, 10% bentonite     36.50
3	     BTEA            9%          81% sand, 10% bentonite     53.46
4	    HDTMA            3%          87% sand, 10% bentonite     20.32
5	    HDTMA            6%          84% sand, 10% bentonite     37.68
6	    HDTMA            9%          81% sand, 10% bentonite     53.46
7	    Shale            3%          87% sand, 10% bentonite      2.08
8	    Shale            6%          84% sand, 10% bentonite      2.20
9	    Shale            9%          81% sand, 10% bentonite      2.31
10	    GAC              3%          87% sand, 10% bentonite     10.65
11	    GAC              6%          84% sand, 10% bentonite     19.34
12	    GAC              9%          81% sand, 10% bentonite     28.03
13	    none             0%          90 % sand, 10% bentonite     1.96
14          n/a (no layer)   0%          not applicable (no layer)    0.00




