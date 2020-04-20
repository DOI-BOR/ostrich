/******************************************************************************
File      : AccessConverter.cpp
Author    : Tsu-Wei Webber Chen
Copyright : 2010, Tsu-Wei Webber Chen

An implementation of Access converter.

Version History
01-25-10    tws   created 
******************************************************************************/
#include <string.h>
#include <string>
using namespace std;

#include "AccessConverter.h"
#include "ADOConnection.h"

#include "Utility.h"
#include "Exception.h"

/******************************************************************************
DeleteASCIIFile()

Delete the ASCII file that contains converted responses.
******************************************************************************/
void AccessConverter::DeleteASCIIFile(void)
{
   if (strncmp(m_accessType, "Read", 4) == 0)
   {
      //deletes converted file, if it exists
      string s_fileName(m_fileName);
      s_fileName = s_fileName.substr(0, s_fileName.find_last_of('.'));
      s_fileName += ".txt";
      remove(s_fileName.c_str());
   } 
}/* end ReadResponses() */

/******************************************************************************
ReadResponse()

Read the requested response and append to an ASCII file.
******************************************************************************/
void AccessConverter::ReadResponse(void)
{
   CreateConnectionString();
   ADOConnection connection(m_connectionString);
   if (strncmp(m_accessType, "Read", 4) == 0)
   {
	   connection.Read(m_table, m_keyColumn, m_key, m_column, m_name, m_fileName);
   } 
}/* end ReadResponses() */

/******************************************************************************
WriteParameter()

Write the requested parameter value to the database.

returns true if write was made, false otherwise (i.e. the database entry doesn't
match the requested parameter name).
******************************************************************************/
bool AccessConverter::WriteParameter(char * pName, char * pValue)
{
   if((strcmp(pName, m_param) == 0) && (strcmp(m_accessType, "Write") == 0))
   {
	   CreateConnectionString();
	   ADOConnection connection(m_connectionString);
		connection.Write(m_table, m_keyColumn, m_key, m_column, pValue);
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
void AccessConverter::InsertDbase(DatabaseABC * pNxt)
{
   if(m_pNxt == NULL) m_pNxt = pNxt;
   else m_pNxt->InsertDbase(pNxt);
}/* end InsertDbase() */

/******************************************************************************
CTOR

Initializes member variables to defaults
******************************************************************************/
AccessConverter::AccessConverter(void)
{
   m_bIsEmpty = true;
   m_pNxt = NULL;
   strcpy(m_connectionString, "");
   strcpy(m_accessType, "");
   strcpy(m_fileName, "");
   strcpy(m_table, "");
   strcpy(m_keyColumn, "");
   strcpy(m_key, "");
   strcpy(m_column, "");
   strcpy(m_param, "");
   strcpy(m_name, "");
}/* end default CTOR */

/******************************************************************************
Initialize()

Initializes info.
******************************************************************************/
void AccessConverter::Initialize(char * line)
{
   m_bIsEmpty = false;
   strcpy(m_connectionString, "");
   strcpy(m_accessType, "");
   strcpy(m_fileName, "");
   strcpy(m_table, "");
   strcpy(m_keyColumn, "");
   strcpy(m_key, "");
   strcpy(m_column, "");
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
   //info table
   j = ExtractString(info, m_table);
   info += j;
   //info key column
   j = ExtractString(info, m_keyColumn);
   info += j;
   //info key
   j = ExtractString(info, m_key);
   info += j;
   //info column
   j = ExtractString(info, m_column);
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
void AccessConverter::Destroy(void)
{
   DatabaseABC * pCur, * pNxt;
   for(pCur = this; pCur != NULL; pCur = pNxt)
   {
      pNxt = pCur->GetNext();
      delete pCur;
   }
}/* end Destroy() */

/******************************************************************************
Convert

Peforms conversion.
******************************************************************************/
void AccessConverter::Convert(void)
{
	CreateConnectionString();
	ADOConnection connection(m_connectionString);
	if (strncmp(m_accessType, "Read", 4) == 0)
	{
	  connection.Read(m_table, m_keyColumn, m_key, m_column, m_name, m_fileName);
	} 
	else if (strncmp(m_accessType, "Write", 5) == 0)
	{
          //hard coded changes to test write function
	  connection.Write(m_table, m_keyColumn, m_key, m_column, "50.00");
	}
}/* end Convert() */

/******************************************************************************
CreateConnectionString

Creates the appropriate connection string.
******************************************************************************/
void AccessConverter::CreateConnectionString(void)
{
   strcpy(m_connectionString, "Provider=Microsoft.ACE.OLEDB.12.0;Data Source=");
   strcat(m_connectionString, m_fileName);
}/* end CreateConnectionString() */

/******************************************************************************
ReadFromFile()

Read in the type conversion section and create a linked list of Access 
converters.

Returns false if the section does not exist or if the section exists but does
not contain any access conversions.
******************************************************************************/
bool AccessConverter::ReadFromFile(void)
{
   AccessConverter * pNew = NULL;
   int j;
   char tmpFileType[DEF_STR_SZ];
   char * lineStr;
   FILE * pFile;

   IroncladString pFileName = GetOstFileName();
   pFile = fopen(pFileName, "r");

   if(pFile == NULL) 
   {
      FileOpenFailure("AccessConverter::ReadFromFile()", pFileName);
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
	   if (strncmp(tmpFileType, "Access", 6) == 0)
	   {
         if(m_bIsEmpty == true)
         {
	         Initialize(lineStr);
         }
         else //add to the linked list
         {
            pNew = new AccessConverter();
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
