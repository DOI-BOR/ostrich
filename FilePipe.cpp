/******************************************************************************
File     : FilePipe.cpp
Author   : L. Shawn Matott and Vijaykumar Raghavan
Copyright: 2003, L. Shawn Matott and Vijaykumar Raghavan

FilePipe classes are used to perform operations on an input/output file pair 
(i.e. FilePair class). Upon creation, a FilePipe reads the contents of the input
file into RAM (as a character array). The FindAndReplace() routine then allows
the contents of the input file to be altered. Typically, this is used to 
replace keywords in the template input file with properly formatted model 
parameter values. Finally, when a FilePipe is destroyed, the (possibly 
modified) input file is written to the desired output file. In combination 
with the FilePair class, the FilePipe class provides the optimization and 
gridding algorithms with a convenient interface for altering model parameters.

Version History
03-05-04    lsm   FindAndReplace() no longer aborts if search string is not found
03-09-03    lsm   added copyright information and initial comments.
08-20-03    lsm   created version history field and updated comments.
07-08-04    lsm   added more descriptive FindAndReplace() error message
08-17-04    lsm   RAM fragmentation fixes
12-10-04    lsm   Added check for existence of model output file.
******************************************************************************/
#include "mpi_stub.h"
#include <string.h>

#include "FilePipe.h"

#include "Exception.h"
#include "Utility.h"

/******************************************************************************
CTOR

Associates an input and output filename with the FilePipe class.
******************************************************************************/
FilePipe::FilePipe(IroncladString in, IroncladString out)
{
   int len;

   len = (int)strlen(in) + 1;
   NEW_PRINT("char", len);
   m_pInFile = new char[len];
   

   len = (int)strlen(out) + 1;
   NEW_PRINT("char", len);
   m_pOutFile = new char[len];
   
   MEM_CHECK(m_pOutFile);

   strcpy(m_pInFile, in);
   strcpy(m_pOutFile, out);
   
   m_pDataStr = NULL; 
   m_DataSize = 0;


   m_pRepStr = NULL; 
   m_RepSize = 0;

   IncCtorCount();
} /* end CTOR */

/******************************************************************************
Destroy()

Frees up memory associated with the filenames of the FilePipe.
******************************************************************************/
void FilePipe::Destroy(void)
{   
   delete [] m_pInFile;
   delete [] m_pOutFile;

   delete [] m_pDataStr;
   m_DataSize = 0;

   delete [] m_pRepStr;
   m_RepSize = 0;

   IncDtorCount();
} /* end Destroy() */

/******************************************************************************
FileToString()

Reads the input file of the FilePipe and stores it into a string.
******************************************************************************/
void FilePipe::FileToString(void)
{
   int fileSize;
   int i;
   FILE * pFile;

   pFile = fopen(m_pInFile, "r");  
   if(pFile == NULL)
   {
      FileOpenFailure("FilePipe::CTOR", m_pInFile);
   }/* end if() */

   /*
   count number of chars in file, 
   so that fileStr can be sized
   */
   fileSize = 0;
   while(feof(pFile) == 0) 
   {
      fileSize++;
      fgetc(pFile);
   }/* end while() */
   fileSize--;

   //resize fileStr, if needed
   if(m_DataSize < (fileSize+1))
   {
      delete [] m_pDataStr;
      m_DataSize = fileSize + 1;
      NEW_PRINT("char", m_DataSize);
      m_pDataStr = new char[m_DataSize];
      MEM_CHECK(m_pDataStr);
   }/* end if() */
      
   //fill fileStr
   rewind(pFile);   
   for(i = 0; i < fileSize; i++)
   {
      m_pDataStr[i] = (char)(fgetc(pFile));
   }/* end for() */
   m_pDataStr[i] = 0;

   fclose(pFile);

   //resize replace string, if needed
   if(m_RepSize < m_DataSize)
   {
      delete [] m_pRepStr;
      m_RepSize = m_DataSize;
      NEW_PRINT("char", m_RepSize);
      m_pRepStr = new char[m_RepSize];      
      MEM_CHECK(m_pRepStr);
   }/* end if() */

   //initialize replace string
   strcpy(m_pRepStr, m_pDataStr);

} /* end FileToString() */

/******************************************************************************
StringToFile()

Writes the string to the output file of the FilePipe. Then, the replacement
string is re-initialized (to prepare for next round of FindAndReplace().
******************************************************************************/
void FilePipe::StringToFile(void)
{   
   FILE * pFile;
   char msg[DEF_STR_SZ];
   
   pFile = fopen(m_pOutFile,"w");

   if(pFile == NULL)
   {
      sprintf(msg, "Couldn't open model output file: |%s|", m_pOutFile);
      LogError(ERR_FILE_IO, msg);
      ExitProgram(1);
   }

// ****************************************
// The following approach is MUCH slower
// than using fwrite() for moderate to large
// file sizes.
//   unsigned int i;
//   for(i = 0; i < strlen(m_pRepStr); i++)
//   {
//      fprintf(pFile, "%c", m_pRepStr[i]);
//   }
// *****************************************   
   int nChars = strlen(m_pRepStr);
   fwrite(m_pRepStr, 1, nChars, pFile);

   fclose(pFile);

   //reset the replace string
   strcpy(m_pRepStr, m_pDataStr);
} /* end StringToFile() */

/******************************************************************************
FindAndReplace()

Finds a search string and replaces it with the replace string.

Returns a 1 if a replacement is made, a 0 if no replacement is made.
******************************************************************************/
int FilePipe::FindAndReplace(IroncladString find, IroncladString replace)
{
   StringType pTmp;
   int len, count, diff;

   //count the number of replacements required
   count = MyStrOccur(m_pRepStr, find);
   if(count == 0){ return 0;}

   //compute the size difference between find and replace
   diff = count * ((int)strlen(replace) - (int)strlen(find));
   if(diff < 0){ diff = 0;}
   
   //compute the minimum size needed for the replacement string
   len = (int)strlen(m_pRepStr) + diff + 1;

   //resize the replacement string, if necessary
   if(m_RepSize < len)
   {
      //printf("***** resizing replacement string ****\n");
      //create buffer and copy the existing string into it.
      m_RepSize = len;
      NEW_PRINT("char", m_RepSize);
      pTmp = new char[m_RepSize];      
      MEM_CHECK(pTmp);
      strcpy(pTmp, m_pRepStr);

      //delete current m_pRepStr and swap in pTmp
      delete [] m_pRepStr;
      m_pRepStr = pTmp;
   }/* end if() */ 

   //perform the replacements
   MyStrRep(m_pRepStr, find, replace);

	return 1;
} /* end FindAndReplace() */

