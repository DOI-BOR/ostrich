/******************************************************************************
File      : ValueExtractor.cpp
Author    : L. Shawn Matott and Vijaykumar Raghavan
Copyright : 2003, L. Shawn Matott and Vijaykumar Raghavan

This class uses the instructions in the observation file to read the output file
of the model program.

Version History
06-12-03    lsm   added copyright information and initial comments.
08-20-03    lsm   created version history field and updated comments.
03-05-04    lsm   added PSO, fixed filename parse error in obs. parse
                  added support for alternative tokens
08-11-04    lsm   made into a linked list to reduce file I/O
03-09-05    lsm   added support for Fortran-style scientific number formation
                  (e.g. 1.000D-4 vs. 1.000E-04)
******************************************************************************/
#include <string.h>
#include <stdlib.h>

#include "ValueExtractor.h"

#include "Exception.h"
#include "Utility.h"

/******************************************************************************
CTOR

Reads model output file into a string.
******************************************************************************/
ValueExtractor::ValueExtractor(IroncladString file, bool bQuitOnErr, double errVal)
{     
   int len;

   len = (int)strlen(file) + 1;
   NEW_PRINT("char", len);
   m_FileName = new char[len];
   MEM_CHECK(m_FileName);
   strcpy(m_FileName, file);

   m_DataStr = NULL;
   m_DataSize = 0;
   m_bQuitOnError = bQuitOnErr;
   m_ErrorVal = errVal;

   m_pNxt = NULL;

   IncCtorCount();
} /* end CTOR */

/******************************************************************************
Destroy()

Frees up memory used by member variables.
******************************************************************************/
void ValueExtractor::Destroy(void)
{
   delete [] m_FileName;
   delete [] m_DataStr;
   m_DataSize = 0;
   delete m_pNxt;

   IncDtorCount();
}/* end Destroy() */

/******************************************************************************
Insert()

Inserts a ValueExtractor() into the linked list, if one with the same file
name has not already been inserted.
******************************************************************************/
void ValueExtractor::Insert(IroncladString name)
{   
   ValueExtractor * pCur, * pIns;

   pCur = GetByName(name);

   //already in the list? Then return....
   if(pCur != NULL){ return;}

   pCur = this;
   NEW_PRINT("ValueExtractor", 1);
   pIns = new ValueExtractor(name, m_bQuitOnError, m_ErrorVal);

   while(pCur->GetNext() != NULL){ pCur = pCur->GetNext();}
   pCur->SetNext(pIns);
} /* end CTOR */

/******************************************************************************
ReadOutputFiles()

Converts the output files (identified by m_pFileName) into strings.
******************************************************************************/
void ValueExtractor::ReadOutputFiles(void)
{   
   ValueExtractor * pCur;

   pCur = this;

   while(pCur != NULL)
   {
      pCur->FileToString();
      pCur = pCur->GetNext();
   }
} /* end ReadOutputFiles() */

/******************************************************************************
GetByName()

Returns a pointer to the ValueEctractor associated with the given file name.
Returns NULL if no ValueExtractor exists for the given file.
******************************************************************************/
ValueExtractor * ValueExtractor::GetByName(IroncladString name)
{   
   ValueExtractor * pCur;

   pCur = this;

   while(pCur != NULL)
   {
      if(strcmp(pCur->GetName(), name) == 0){ return pCur;}
      pCur = pCur->GetNext();
   }
   return NULL;
} /* end GetByName() */

/******************************************************************************
FileToString()

Reads a file into a string.
******************************************************************************/
void ValueExtractor::FileToString(void)
{
   int fileSize;
   int i;
   FILE * pFile;

   pFile = fopen(m_FileName, "r");
   if(pFile == NULL)
   {
      FileOpenFailure("ValueExtractor::CTOR", m_FileName);
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

   //size fileStr, if necessary
   if(m_DataSize < (fileSize+1))
   {
      //printf("**** Resizing value extractor file string ****\n");
      delete [] m_DataStr;
      m_DataSize = (fileSize+1);
      NEW_PRINT("char", m_DataSize);
      m_DataStr = new char[m_DataSize];
      MEM_CHECK(m_DataStr);
   }/* end if() */

   //fill fileStr
   rewind(pFile);
   for(i = 0; i < fileSize; i++)
   {
      m_DataStr[i] = (char)(fgetc(pFile));
   }/* end for() */
   m_DataStr[i] = 0;

   fclose(pFile);
} /* end FileToString() */

/******************************************************************************
ExtractValue()

Positions file (name) string at the line containing the search arg. Then uses 
the line and col values to locate the position of the desired value. This value 
is then extracted, converted to a double and stored in the val argument.

Returns true if extraction is successful, false if an error occurs.
******************************************************************************/
bool ValueExtractor::ExtractValue
(
   IroncladString name, 
   IroncladString search, 
   int line, 
   int col,
   char tok,
   double * val
)
{
   bool errVal;
   ValueExtractor * pExtractor;

   pExtractor = GetByName(name);
   errVal = pExtractor->ExtractValue(search, line, col, tok, val);
   return errVal;
} /* end extractValue() */

/******************************************************************************
ExtractValue()

Positions file string at the line containing the search arg. Then uses 
the line and col values to locate the position of the desired value. This value 
is then extracted, converted to a double and stored in the val argument.

Returns true if extraction is successful, false if an error occurs.
******************************************************************************/
bool ValueExtractor::ExtractValue
(
   IroncladString search,
   int line, 
   int col,
   char tok,
   double * val
)
{   
   UnchangeableString curPos; //current position within the file string
   char * msg;
   char tmp;   
   int i, j;

   int max_msg_size = GetMaxLineSizeInString(m_DataStr);
   msg = new char[max_msg_size];

   if(strcmp(search, "OST_NULL") == 0)
   {
      curPos = m_DataStr;
   }
   else
   {
      curPos = strstr(m_DataStr, search);
      if(curPos == NULL)
      {      
         //set error
         sprintf(msg, "extractValue(): strstr() failed : couldn't find |%s|", search);
         LogError(ERR_FILE_IO, msg);
		 *val = m_ErrorVal;
       delete [] msg;
       if(m_bQuitOnError == true)
		    return false;
       else
          return true;
      }/* end if() */
   }/* end else() */
   
   //advance to the desired line 
   i = 0;
   while(i < line)
   {
      tmp = *curPos;
      curPos++;

      if(tmp == '\n') {i++;}
      else if(tmp == 0)
      {
         LogError(ERR_FILE_IO, "extractValue(): could not locate line");
		   *val = m_ErrorVal;
         delete [] msg;
         if(m_bQuitOnError == true)
		      return false;
         else
            return true;
      }/* end else if() */
   } /* end while() */

   //advance to the desired column
   for(i = 0; i < col; i++)
   {
      if(tok == ' ')
      {
         j = ExtractString(curPos, msg);
      }
      else
      {
         j = ExtractColString(curPos, msg, tok);
      }
      j = CheckExtraction(j, i, col, "ExtractValue()");
	  if(j == -1)
	  {
		  *val = m_ErrorVal;
        delete [] msg;
        if(m_bQuitOnError == true)
		     return false;
        else
           return true;
	  }
      curPos += j;
   } /* end for() */

   /*------------------------------------------------------
   Convert value from string to double. Some programs (such 
   as those written in Fortran) output number formats that
   are incompatible with a given implementations of the 
   atof() function. Therefore, the string is manipulated
   slightly prior to calling atof().
      1. Replace occurences of 'd' and 'D' with 'E' so that
         1.000D-003 --> 1.000E-003
   ------------------------------------------------------*/
   MyStrRep(msg,"D", "E");
   MyStrRep(msg,"d", "E");
   *val = atof(msg);
   delete [] msg;

   return true;
} /* end extractValue() */

