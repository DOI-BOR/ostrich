/******************************************************************************
File      : Tokens.h
Author    : L. Shawn Matott
Copyright : 2006, L. Shawn Matott

Contains a list of parsing tokens.

Version History
05-16-06    lsm   created
******************************************************************************/
#ifndef TOKENS_H
#define TOKENS_H

#define BEG_GCFG_TOK "BeginGeneralConfiguration"
#define END_GCFG_TOK "EndGeneralConfiguration"
#define SOL_TOLR_TOK "Tolerance"
#define MAX_STEP_TOK "MaxTimeStep"
#define CNC_UNIT_TOK "ConcentrationUnits"
#define TIM_UNIT_TOK "TimeUnits"
#define SIM_TIME_TOK "SimulationTime"
#define MAX_RTIM_TOK "MaxRunTime"
#define ODE_SOLV_TOK "Solver"
#define ODE_SIMR_TOK "simr" //semi-implicit midpoint rule
#define ODE_ADRK_TOK "ark" //adaptive Runge-Kutta
#define ODE_CVDE_TOK "cvode" //CVODE dense solver 
#define CHG_BLNC_TOK "ChargeBalance" //enable/disable adjustment of pH to maintain charge balance
#define MSS_BLNC_TOK "MassBalance"   //enable/disable adjustment of mass to maintain non-negative solids mass
#define BEG_AQSP_TOK "BeginAqueousSpecies"
#define END_AQSP_TOK "EndAqueousSpecies"
#define RXN_SEPR_TOK "-->" //reaction separator
#define BEG_AQRX_TOK "BeginAqueousReactions"
#define END_AQRX_TOK "EndAqueousReactions"
#define OPT_ENBL_TOK "yes" //enable some option
#define OPT_DSBL_TOK "no" //disable some option
#define BEG_MNSP_TOK "BeginMineralSpecies"
#define END_MNSP_TOK "EndMineralSpecies"
#define BEG_IXSP_TOK "BeginIonExchangeSpecies"
#define END_IXSP_TOK "EndIonExchangeSpecies"
#define BEG_SFSP_TOK "BeginSurfaceComplexationSpecies"
#define END_SFSP_TOK "EndSurfaceComplexationSpecies"
#define BEG_MNRX_TOK "BeginMineralReactions"
#define END_MNRX_TOK "EndMineralReactions"
#define BEG_IXRX_TOK "BeginIonExchangeReactions"
#define END_IXRX_TOK "EndIonExchangeReactions"
#define BEG_SFCX_TOK "BeginSurfaceComplexationReactions"
#define END_SFCX_TOK "EndSurfaceComplexationReactions"
#endif /* TOKENS_H */


