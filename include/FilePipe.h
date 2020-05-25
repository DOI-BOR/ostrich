/******************************************************************************
File     : FilePipe.h
Author   : L. Shawn Matott and Vijaykumar Raghavan
Copyright: 2003, L. Shawn Matott and Vijaykumar Raghavan

FilePipe classes are used to perform operations on in input/output file pair 
(i.e. FilePair class). Upon creation, a FilePipe reads the cntents of the input
file into RAM (as a character array). The FindAndReplace() routine then allows
the contents of the input file to be altered. Typically, this is used to 
replace keywords in the template input file with properly formatted model 
parameter values. Finally, when a FilePipe is destroyed, the (possibly 
modified) input file is written to the desired output file. In combination 
with the FilePair class, the FilePipe class provides the optimization and 
gridding algorithms with a convenient interface for altering model parameters.

Version History
03-09-03    lsm   added copyright information and initial comments.
08-20-03    lsm   created version history field and updated comments.
08-17-04    lsm   RAM fragmentation fixes
******************************************************************************/
#ifndef FILE_PIPE_H
#define FILE_PIPE_H

#include "MyHeaderInc.h"

/******************************************************************************
class FilePipe
   The class mainly deals with the I/O opearations of a file.
   
   Hence each object of this class is associated with a file and
   the whole file is first read and stored into a STRING.

   Operations are done on the STRING and the STRING is then 
   written into a file when the object of this class is destroyed.
********************************************************************************/
class FilePipe
{ 
   private:
      StringType m_pInFile;
      StringType m_pOutFile;

      /* template file string */
      StringType m_pDataStr;
      int m_DataSize;

      /* replacement file string */
      StringType m_pRepStr;
      int m_RepSize;

   public :
      FilePipe(IroncladString in, IroncladString out);
	  ~FilePipe(void){ DBG_PRINT("FilePipe::DTOR"); Destroy(); }
      void Destroy(void);

      int FindAndReplace(IroncladString find, IroncladString replace);
      void FileToString(void);
      void StringToFile(void);
      StringType GetTemplateFileName(void){ return m_pInFile;}
      StringType GetModelInputFileName(void){ return m_pOutFile;}
}; /* end class FilePipe */

#endif /* FILE_PIPE_H */




