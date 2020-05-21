/******************************************************************************
File     : ADOConnection.h
Author   : Tsu-Wei Webber Chen
Copyright: 2010, Tsu-Wei Webber Chen

ADOConnection connects to an ADO.NET database.

Version History
01-25-10    twc   Created file.
******************************************************************************/
#ifndef ADO_CONNECTION_H
#define ADO_CONNECTION_H

#ifdef WIN32

#import "msado15.dll"  \
    rename( "EOF", "AdoNSEOF" )

#include "MyHeaderInc.h"

/******************************************************************************
class ADOConnection
******************************************************************************/
class ADOConnection
{
   public:
      ADOConnection(char * connectionString);
      void Read(char * table, char * keyColumn, char * key, char * column, char * name, char * fileName);
      void Write(char * table, char * keyColumn, char * key, char * column, const char * param);
   private:
	   HRESULT hr;
	   _bstr_t bstrConnect;
}; /* end class ADOConnection */
#else /* not WIN32 */
/******************************************************************************
class ADOConnection
******************************************************************************/
class ADOConnection
{
   public:
     ADOConnection(char * connectionString);
     void Read(char * table, char * keyColumn, char * key, char * column, char * name, char * fileName);
     void Write(char * table, char * keyColumn, char * key, char * column, const char * param);
}; /* end class ADOConnection */
#endif /* not WIN32 */
#endif /* ADO_CONNECTION_H */

