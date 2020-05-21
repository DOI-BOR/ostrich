/******************************************************************************
File     : FileList.cpp
Author   : L. Shawn Matott
Copyright: 2008, L. Shawn Matott

FileList classes are used to store a collection of files that Ostrich needs to
delete when its done running. These files are executables and and extra input 
files.  Files are deleted to conserve disk space, which is required for large
parallel runs.

Version History
07-05-08    lsm   created
******************************************************************************/
#include <string.h>
#include <stdlib.h>

#include "FileList.h"

#include "Utility.h"
#include "Exception.h"

/******************************************************************************
CTOR

Creates a file list.
******************************************************************************/
FileList::FileList(IroncladString name)
{
   strcpy(m_Name, name); 
   m_pNxt = NULL; 
   IncCtorCount();
} /* end CTOR */

/******************************************************************************
Destroy()

Frees up the file list.
******************************************************************************/
void FileList::Destroy(void)
{   
   delete m_pNxt;
   IncDtorCount();
}/* end Destroy() */

/******************************************************************************
Insert()

Insert an item into the file list.
******************************************************************************/
void FileList::Insert(IroncladString name)
{
   if(m_pNxt == NULL)
   {
      m_pNxt = new FileList(name);
   }
   else
   {
      m_pNxt->Insert(name);
   }
}/* end Insert() */

/******************************************************************************
Cleanup()

Delete the files in the list.
******************************************************************************/
void FileList::Cleanup(IroncladString dir)
{
   static bool bLogged = false;  //only log once to reduce output file size
   char tmp[DEF_STR_SZ];
   FileList * pCur;
   MY_CHDIR(dir);
   for(pCur = this; pCur != NULL; pCur = pCur->GetNext())
   {
      strcpy(tmp, pCur->GetName());
      if(tmp[0] == '"')
      {
         tmp[0] = ' ';
         tmp[strlen(tmp)-1] = ' ';
         MyTrim(tmp);
      }
      if(MY_ACCESS(tmp, 0) != -1)
      {
         #ifdef WIN32
            sprintf(tmp, "del %s 1>> %s 2>>&1", pCur->GetName(), GetOstExeOut());
         #else
            sprintf(tmp, "rm %s 2>&1 | >> %s", pCur->GetName(), GetOstExeOut());
         #endif
         system(tmp);
         if(bLogged == false)
         {         
            sprintf(tmp, "Ostrich deleted %s/%s", dir, pCur->GetName());
            LogError(ERR_CLEANUP, tmp);
         }/* end if(logged) */
      }/* end if(file exists) */
   }/* end for(each file) */
   bLogged = true;
   MY_CHDIR("..");
}/* end Cleanup() */
