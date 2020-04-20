/******************************************************************************
File      : SuperMuseUtility.cpp
Author    : L. Shawn Matott
Copyright : 2007, L. Shawn Matott

This file contains a c-style interface to the global SuperMUSE instance.

Version History
07-16-07    lsm   Created.
******************************************************************************/
#include "SuperMuseUtility.h"

#include "ModelABC.h"

#include "SuperMUSE.h"

/* -------------------------------
Global flag to indicate whether or 
not the SuperMUSE system has been 
'enabled' by the user.
------------------------------- */
bool gUseSuperMUSE = false;
SuperMUSE * gSuperMUSE = NULL;

/******************************************************************************
IsSuperMUSE()

Returns true if SuperMUSE is enabled, false otherwise.
******************************************************************************/
bool IsSuperMUSE(void)
{
   return gUseSuperMUSE;
}/* end IsSuperMUSE() */

/******************************************************************************
GetSuperMusePtr()

Returns a pointer to a global instance of the SuperMUSE class.
******************************************************************************/
SuperMUSE * GetSuperMusePtr(void)
{
   return gSuperMUSE;
}/* end GetSuperMusePtr() */

/******************************************************************************
EnableSuperMUSE()

Enables the use of the SuperMUSE system, which is assumed to be present based
on user-specification of the SuperMUSE field of the OstIn.txt input file.
******************************************************************************/
void EnableSuperMUSE(void)
{
   gUseSuperMUSE = true;
}/* end EnableSuperMUSE() */

/******************************************************************************
DisableSuperMUSE()

Disables the use of the SuperMUSE system. If SuperMUSE tasker reports an error
or times out, Ostrich will disable the SuperMUSE system and revert to serial
(single processor) execution.
******************************************************************************/
void DisableSuperMUSE(void)
{
   gUseSuperMUSE = false;
}/* end DisableSuperMUSE() */

/******************************************************************************
InitSuperMUSE()

Initialize SuperMUSE based on user input file.
******************************************************************************/
void InitSuperMUSE(FILE * pFile, ModelABC * pModel)
{
   gSuperMUSE = new SuperMUSE(pFile, pModel);
}/* end InitSuperMUSE() */

/******************************************************************************
CleanSuperMUSE()

Cleanup SuperMUSE environment variables.
******************************************************************************/
void CleanSuperMUSE(void)
{
   gSuperMUSE->EnvVarCleanup();
}/* end CleanSuperMUSE() */

/******************************************************************************
DestroySuperMUSE()

Perform garbage collection of the SuperMUSE class.
******************************************************************************/
void DestroySuperMUSE(void)
{   
   delete gSuperMUSE;
   gSuperMUSE = NULL;
}/* end DestroySuperMUSE() */

/******************************************************************************
WriteSuperMuseSetupToFile()

Write out SuperMUSE configuration.
******************************************************************************/
void WriteSuperMuseSetupToFile(FILE * pFile)
{
   if(gSuperMUSE != NULL)
   {
      gSuperMUSE->WriteSetup(pFile);
   }
}/* end WriteSuperMuseSetupToFile() */
