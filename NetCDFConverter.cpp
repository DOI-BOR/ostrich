/******************************************************************************
File      : NetCDFConverter.cpp
Author    : Tsu-Wei Webber Chen
Copyright : 2010, Tsu-Wei Webber Chen

An implementation of NetCDF converter.

Version History
02-15-10    tws   created 
******************************************************************************/
#include <string.h>
#include <string>
using namespace std;

#include "NetCDFConverter.h"
#include "ADOConnection.h"

#include "Exception.h"
#include "Utility.h"

/******************************************************************************
DeleteASCIIFile()

Delete the ASCII file that contains converted responses.
******************************************************************************/
void NetCDFConverter::DeleteASCIIFile(void)
{
   if (strncmp(m_accessType, "Read", 4) == 0)
   {
      //deletes converted file, if it exists
      string s_fileName(m_fileName);
      s_fileName = s_fileName.substr(0, s_fileName.find_last_of('.'));
      s_fileName += ".txt";
		if (remove(s_fileName.c_str()) != 0)
		{
			int i = 0;	
		}
   } 
}/* end ReadResponses() */

/******************************************************************************
ReadResponse()

Read the requested response and append to an ASCII file.
******************************************************************************/
void NetCDFConverter::ReadResponse(void)
{
   if (strncmp(m_accessType, "Read", 4) == 0)
   {
	   string s_fileName(m_fileName);
	   string cmd(m_command);

      s_fileName = s_fileName.substr(0, s_fileName.find_last_of('.'));
      s_fileName += ".txt";
	   cmd = "nc2text " + string(m_fileName) + " " + string(m_arrayName) + "[" + string(m_itemPos) + "]";
		ExecuteCommandLine(cmd.c_str(), true, s_fileName.c_str(), m_name);
   } 
}/* end ReadResponses() */

/******************************************************************************
WriteParameter()

Write the requested parameter value to the database.

returns true if write was made, false otherwise (i.e. the database entry doesn't
match the requested parameter name).
******************************************************************************/
bool NetCDFConverter::WriteParameter(char * pName, char * pValue)
{
   if((strcmp(pName, m_param) == 0) && (strcmp(m_accessType, "Write") == 0))
   {
	   string s_fileName(m_fileName);
	   string cmd(m_command);

      s_fileName = s_fileName.substr(0, s_fileName.find_last_of('.'));
      s_fileName += ".txt";
	   cmd = "echo " + string(pValue) + " | text2nc "+ string(m_fileName) + " " + string(m_arrayName) + "[" + string(m_itemPos) + "]";
		ExecuteCommandLine(cmd.c_str(), false, s_fileName.c_str(), m_name);
      return true;
	}
   else
   {
      return false;
   }
}/* end WriteParameter() */

/******************************************************************************
InsertDbase()

Insert a database conversion at the end of the list.
******************************************************************************/
void NetCDFConverter::InsertDbase(DatabaseABC * pNxt)
{
   if(m_pNxt == NULL) m_pNxt = pNxt;
   else m_pNxt->InsertDbase(pNxt);
}/* end InsertDbase() */

/******************************************************************************
CTOR

Initializes member variables to defaults
******************************************************************************/
NetCDFConverter::NetCDFConverter(void)
{
   m_bIsEmpty = true;
   m_pNxt = NULL;
   strcpy(m_command, "");
   strcpy(m_accessType, "");
   strcpy(m_fileName, "");
   strcpy(m_arrayName, "");
   strcpy(m_itemPos, "");
   strcpy(m_param, "");
   strcpy(m_name, "");
}/* end default CTOR */

/******************************************************************************
Initialize()

Initializes info.
******************************************************************************/
void NetCDFConverter::Initialize(char * line)
{
   m_bIsEmpty = false;
   strcpy(m_command, "");
   strcpy(m_accessType, "");
   strcpy(m_fileName, "");
   strcpy(m_arrayName, "");
   strcpy(m_itemPos, "");
   strcpy(m_param, "");
   strcpy(m_name, "");

   char * info = line;
   int j;
   //info file name
   j = ExtractString(info, m_fileName);
   info += j;
   //info access type
   j = ExtractString(info, m_accessType);
   info += j;
   //info array name
   j = ExtractString(info, m_arrayName);
   info += j;
   //info item position
   j = ExtractString(info, m_itemPos);
   info += j;
   if (strncmp(m_accessType, "Read", 4) == 0)
   {
     //info column
	 j = ExtractString(info, m_name);
   }
   else if (strncmp(m_accessType, "Write", 5) == 0)
   {
	 //info column
	 j = ExtractString(info, m_param);
   }
}/* end default CTOR */

/******************************************************************************
Destroy()

Free up memory in the linked list.
******************************************************************************/
void NetCDFConverter::Destroy(void)
{
   DatabaseABC * pCur, * pNxt;
   for(pCur = this; pCur != NULL; pCur = pNxt)
   {
      pNxt = pCur->GetNext();
      delete pCur;
   }
}/* end Destroy() */

/******************************************************************************
ReadFromFile()

Read in the type conversion section and create a linked list of Access 
converters.

Returns false if the section does not exist or if the section exists but does
not contain any access conversions.
******************************************************************************/
bool NetCDFConverter::ReadFromFile(void)
{
   NetCDFConverter * pNew = NULL;
   int j;
   char tmpFileType[DEF_STR_SZ];
   char * lineStr;
   FILE * pFile;

   IroncladString pFileName = GetOstFileName();
   pFile = fopen(pFileName, "r");

   if(pFile == NULL) 
   {
      FileOpenFailure("NetCDFConverter::ReadFromFile()", pFileName);
      return false;
   }/* end if() */

   if(CheckToken(pFile, "BeginTypeConversion", pFileName) == false)
   {
      return false;
   }

   //make sure correct tokens are present
   rewind(pFile);
   FindToken(pFile, "BeginTypeConversion", pFileName);
   FindToken(pFile, "EndTypeConversion", pFileName);
   rewind(pFile);

   FindToken(pFile, "BeginTypeConversion", pFileName);
   lineStr = GetNxtDataLine(pFile, pFileName);

   rewind(pFile);
   FindToken(pFile, "BeginTypeConversion", pFileName);
   lineStr = GetNxtDataLine(pFile, pFileName);

   while(strstr(lineStr, "EndTypeConversion") == NULL)
   {
      j = ExtractString(lineStr, tmpFileType);
	   lineStr += j;
	   //determine conversion type
	   if (strncmp(tmpFileType, "NetCDF", 6) == 0)
	   {
         if(m_bIsEmpty == true)
         {
	         Initialize(lineStr);
         }
         else //add to the linked list
         {
            pNew = new NetCDFConverter();
            pNew->Initialize(lineStr);
            InsertDbase((DatabaseABC *)pNew);
         }
	   }
	   else //unsupported type
	   {
         LogError(ERR_FILE_IO, "Unsupported database type");
	   }

      lineStr = GetNxtDataLine(pFile, pFileName);
   }/* end while() */
   fclose(pFile);

   if(m_bIsEmpty == true) return false;
   return true;
}/* end ReadFromFile() */

