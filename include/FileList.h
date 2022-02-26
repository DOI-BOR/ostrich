/******************************************************************************
File     : FileList.h
Author   : L. Shawn Matott
Copyright: 2008, L. Shawn Matott

FileList classes are used to store a collection of files that Ostrich needs to
delete when its done running. These files are executables and and extra input
files.  Files are deleted to conserve disk space, which is required for large
parallel runs.

Version History
07-05-08    lsm   created
******************************************************************************/
#ifndef FILE_LIST_H
#define FILE_LIST_H

#include "MyHeaderInc.h"
#ifdef GCC5X
#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;
#else
#include <filesystem>
namespace fs = std::filesystem;
#endif

#include <string>
#include <vector>

/******************************************************************************
class FileList

A container class for the names of files. Implemented as a linked list.
******************************************************************************/
class FileList
{
   public:
       FileList(IroncladString name);
      ~FileList(void){ DBG_PRINT("FileList::DTOR"); Destroy(); }
      void Destroy(void);
      void Insert(IroncladString name);
      void Cleanup(fs::path dir);
      FileList * GetNext(void){ return m_pNxt;}
      IroncladString GetName(void){ return m_Name;}

   private:
      char m_Name[DEF_STR_SZ];
      FileList * m_pNxt;
}; /* end class FileList */

#endif /* FILE_LIST_H */
