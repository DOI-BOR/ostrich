/******************************************************************************
ile      : ADOConnection.cpp
Author    : Tsu-Wei Webber Chen
Copyright : 2010, Tsu-Wei Webber Chen

An implementation of ADO.NET connection.

Version History
08-01-07    tws   created 
******************************************************************************/
#include "ADOConnection.h"

#ifdef WIN32
#include <string>
#include <iostream>
using namespace std;

#include "Exception.h"

/******************************************************************************
CTOR

Initializes parameters with the appropriate connection string.
******************************************************************************/
ADOConnection::ADOConnection(char * connectionString)
{
	string cString(connectionString);
   bstrConnect = cString.c_str();
   hr = CoInitialize(NULL);
   if (FAILED(hr))
   {
	   LogError(ERR_FILE_IO, "ADOConnection(): Failed to CoInitialize() COM.");
   }
}/* end default CTOR */

/******************************************************************************
Read

Reads an ADO database.
******************************************************************************/
void ADOConnection::Read(char * table, char * keyColumn, char * key, char * column, char * name, char * fileName)
{
   try
   {
	   FILE * pFile;
	   string s_fileName(fileName);
      s_fileName = s_fileName.substr(0, s_fileName.find_last_of('.'));
      s_fileName += ".txt";
	   pFile = fopen(s_fileName.c_str(), "a"); 
	   string s_table(table);
	   string s_keyColumn(keyColumn);
	   string s_key(key);
	   string s_column(column);
	   string s_query = "SELECT " + s_column + " FROM " + s_table + " WHERE " + s_keyColumn + "=" + s_key;
       ADODB::_ConnectionPtr pConn("ADODB.Connection");
       hr = pConn->Open(bstrConnect, "admin", "", ADODB::adConnectUnspecified);
       if (SUCCEEDED(hr))
       {
            // Prepare SQL query.
			_bstr_t query = s_query.c_str();

            // Excecute the query and create a record set
            ADODB::_RecordsetPtr pRS("ADODB.Recordset");
            hr = pRS->Open(query, 
                    _variant_t((IDispatch *) pConn, true), 
                    ADODB::adOpenUnspecified,
                    ADODB::adLockUnspecified, 
                    ADODB::adCmdText);
            if (SUCCEEDED(hr))
            {
				ADODB::Fields* pFields = NULL;
                hr = pRS->get_Fields(&pFields);
                while (!pRS->AdoNSEOF)
                {
                    for (long nIndex=0; nIndex < pFields->GetCount(); nIndex++)
                    {
						string outputLine(name);
						string value(_bstr_t(pFields->GetItem(nIndex)->GetValue()));
						outputLine += (" " + value + "\n");
						fprintf(pFile, outputLine.c_str());
                    }
                    pRS->MoveNext();
                }
            }
            
            pRS->Close();
            pConn->Close();
        }
        else
        {
			LogError(ERR_FILE_IO, "ADOConnection(): Unable to connection to data source - " + bstrConnect);
        }
	    fclose(pFile);
    }
    catch(_com_error& e)
    {
		LogError(ERR_FILE_IO, "ADOConnection(): " + e.Description());
    }

    CoUninitialize();
}/* end Read */
/******************************************************************************
Write

Writes to an ADO database.
******************************************************************************/
void ADOConnection::Write(char * table, char * keyColumn, char * key, char * column, const char * param)
{
   try
   {
	   string s_table(table);
	   string s_keyColumn(keyColumn);
	   string s_key(key);
	   string s_column(column);
	   string s_param(param);
	   string s_query = "UPDATE " + s_table + " SET " + s_column + "='" + s_param + "' WHERE " + s_keyColumn + "=" + s_key;
       ADODB::_ConnectionPtr pConn("ADODB.Connection");
       hr = pConn->Open(bstrConnect, "admin", "", ADODB::adConnectUnspecified);
       if (SUCCEEDED(hr))
       {
            // Prepare SQL query.
			_bstr_t query = s_query.c_str();

            // Excecute the query and create a record set
            ADODB::_RecordsetPtr pRS("ADODB.Recordset");
            hr = pRS->Open(query, 
					pConn.GetInterfacePtr(), 
                    ADODB::adOpenUnspecified,
                    ADODB::adLockUnspecified, 
                    ADODB::adCmdText);
            pConn->Close();
        }
        else
        {
			LogError(ERR_FILE_IO, "ADOConnection(): Unable to connection to data source - " + bstrConnect);
        }
    }
    catch(_com_error& e)
    {
		LogError(ERR_FILE_IO, "ADOConnection(): " + e.Description());
    }

    CoUninitialize();
}/* end Write */

#else /* not WIN32, no support */
/******************************************************************************
CTOR

Initializes parameters with the appropriate connection string.
******************************************************************************/
ADOConnection::ADOConnection(char * connectionString)
{

}/* end default CTOR */
/******************************************************************************
Read

Reads an ADO database.
******************************************************************************/
void ADOConnection::Read(char * table, char * keyColumn, char * key, char * column, char * name, char * fileName)
{
}/* end Read */
/******************************************************************************
Write

Writes to an ADO database.
******************************************************************************/
void ADOConnection::Write(char * table, char * keyColumn, char * key, char * column, const char * param)
{
}/* end Write */
#endif /* not WIN32 */

