/******************************************************************************
File     : AccessConverter.h
Author   : Tsu-Wei Webber Chen
Copyright: 2010, Tsu-Wei Webber Chen

NetCDFConverter deals with Network Common Data Form files.

Version History
02-15-10    twc   Created file.
******************************************************************************/
#ifndef NETCDF_CONVERTER_H
#define NETCDF_CONVERTER_H

#include "MyHeaderInc.h"

//parent class
#include "DatabaseABC.h"

/******************************************************************************
class NetCDFConverter

An instance of the DatabaseABC abstract base class.
******************************************************************************/
class NetCDFConverter : public DatabaseABC
{
   public:
	  NetCDFConverter(void);
     void Initialize(char * line);
	  bool ReadFromFile(void);
	  ~NetCDFConverter(void){ DBG_PRINT("NetCDFConverter::DTOR"); Destroy(); }
	  void Destroy(void);

     void InsertDbase(DatabaseABC * pNxt);
     DatabaseABC * GetNext(void) { return m_pNxt;}
     bool WriteParameter(char * pName, char * pValue);
     void ReadResponse(void);
     void DeleteASCIIFile(void);

   private:
      DatabaseABC * m_pNxt;
      bool m_bIsEmpty; //true if the converter has not been initialized      
	  char m_command[DEF_STR_SZ];
	  char m_accessType[DEF_STR_SZ];
      char m_fileName[DEF_STR_SZ];
      char m_arrayName[DEF_STR_SZ];
      char m_itemPos[DEF_STR_SZ];
      char m_param[DEF_STR_SZ];
      char m_name[DEF_STR_SZ];
}; /* end class NetCDFConverter */

#endif /* NETCDF_CONVERTER_H */

