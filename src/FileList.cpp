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
void FileList::Cleanup(std::filesystem::path dir)
{

    // Cleanup the files
    FileList* pCur;
    for (pCur = this; pCur != NULL; pCur = pCur->GetNext())
    {
        // Get the filename to be checked
        std::string file = pCur->GetName();

        // Construct the path to the file
        std::filesystem::path filepath = dir;
        filepath /= file;

        // If the file is in the archive directory, delete it
        if (std::filesystem::exists(filepath)) {
            // Delete the file
            std::filesystem::remove(filepath);

            // Log any removal error
            IroncladString filepathString = &filepath.string()[0];
            LogError(ERR_CLEANUP, filepathString);
        }
    }

    // Remove any empty directories in the archived run folder
    // Create a list of entries in the archive folder
    std::vector<std::filesystem::path> archivePaths;
    for (auto& p : std::filesystem::recursive_directory_iterator(dir)) {
        archivePaths.push_back(p);
    }

    // Loop over the files backward, deleting if empty. The reverse loop is necessary to prevent deleting folders that are subsequently checked.
    for (std::vector<std::filesystem::path>::reverse_iterator rit = archivePaths.rbegin(); rit != archivePaths.rend(); ++rit) {
        if (std::filesystem::is_empty(*rit)) {
            std::filesystem::remove(*rit);
        }
    }
}/* end Cleanup() */
