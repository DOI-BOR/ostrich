/******************************************************************************
File     : FilePair.cpp
Author   : L. Shawn Matott and Vijaykumar Raghavan
Copyright: 2003, L. Shawn Matott and Vijaykumar Raghavan

FilePair classes are used to associate a template Ostrich input file with a
Model input file. Keywords in the template input file are replaced with 
properly formatted model parameter values and a valid Model input file is 
generated. In combination with the FilePipe class, the FilePair class provides 
the optimization and gridding algorithms with a convenient interface for 
altering model parameters.

Version History
03-09-03    lsm   added copyright information and initial comments.
08-20-03    lsm   created version history field and updated comments.
08-11-04    lsm   added FilePipe class member, as part of RAM fragmentation fix
******************************************************************************/
#include <string.h>

#include "FilePair.h"
#include "FilePipe.h"

#include "Exception.h"

/******************************************************************************
CTOR

Creates a file pair.
******************************************************************************/
FilePair::FilePair(IroncladString in, IroncladString out)
{
   int len;

   len = (int)strlen(in) + 1;
   NEW_PRINT("char", len);
   m_pInFile = new char[len];

   len = (int)strlen(out) + 1;
   NEW_PRINT("char", len);
   m_pOutFile = new char[len];   
   MEM_CHECK(m_pOutFile);

   m_pNxt = NULL;

   //create file pipe and read template file into memory
   NEW_PRINT("FilePipe", 1);
   m_pPipe = new FilePipe(in, out);   
   MEM_CHECK(m_pPipe);
   m_pPipe->FileToString();

   strcpy(m_pInFile, in);
   strcpy(m_pOutFile, out);

   IncCtorCount();
} /* end CTOR */

/******************************************************************************
Destroy()

Frees up the file pair list.
******************************************************************************/
void FilePair::Destroy(void)
{   
   delete [] m_pInFile;
   delete [] m_pOutFile;
   delete m_pPipe;
   delete m_pNxt;
   IncDtorCount();
}/* end Destroy() */

/******************************************************************************
InsertPair()

Insert a FilePair into the file pair list.
******************************************************************************/
void FilePair::InsertPair(FilePair * pNxt)
{
   FilePair * pCur;

   pCur = this;

   while(pCur->GetNext() != NULL){ pCur = pCur->GetNext();}

   pCur->SetNext(pNxt);
}/* end InsertPair() */

/******************************************************************************
SetNext()

Sets the pointer to the next FilePair in the file pair list.
******************************************************************************/
void FilePair::SetNext(FilePair * pNxt)
{
   delete m_pNxt;
   m_pNxt = pNxt;
}/* end SetNext() */
