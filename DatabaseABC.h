/******************************************************************************
File     : DatabaseABC.h
Author   : L. Shawn Matott
Copyright: 2010, L. Shawn Matott

Encapsulates an interface to database conversion classes developed by Webber Chen.
The conversion classes allow Ostrich to link with models that use database I/O 
instead of text files.

Version History
02-09-10    lsm   created
******************************************************************************/
#ifndef DATABASE_ABC_H
#define DATABASE_ABC_H

#include "MyHeaderInc.h"

/******************************************************************************
class DatabaseABC

An interface class for database conversion classes. Implemented as a linked list.
******************************************************************************/
class DatabaseABC
{
   public:
      virtual ~DatabaseABC(void) { DBG_PRINT("DatabaseABC::DTOR"); }
      virtual void Destroy(void) = 0;
      virtual bool ReadFromFile(void) = 0;
      virtual void InsertDbase(DatabaseABC * pNxt) = 0;
      virtual DatabaseABC * GetNext(void) = 0;
      virtual bool WriteParameter(char * pName, char * pValue) = 0;
      virtual void ReadResponse(void) = 0;
      virtual void DeleteASCIIFile(void) = 0;
   private:
}; /* end class DatabaseABC */

#endif /* DATABASE_ABC_H */



