/******************************************************************************
File     : FilePair.h
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
******************************************************************************/
#ifndef FILE_PAIR_H
#define FILE_PAIR_H

#include "MyHeaderInc.h"

//forward declararation
class FilePipe;

/******************************************************************************
class FilePair

A container class for the names of two files. Implemented as a linked list.
******************************************************************************/
class FilePair
{
   public:
       FilePair(IroncladString in, IroncladString out);
      ~FilePair(void){ DBG_PRINT("FilePair::DTOR"); Destroy(); }
      void Destroy(void);
      void InsertPair(FilePair * pNxt);
      FilePair * GetNext(void) { return m_pNxt;}
      FilePipe * GetPipe(void) { return m_pPipe;}

   private:
      void SetNext(FilePair * pNxt);

      StringType m_pInFile;
      StringType m_pOutFile;
      FilePipe * m_pPipe;
      FilePair * m_pNxt;      
}; /* end class FilePair */

#endif /* FILE_PAIR_H */



