/******************************************************************************
File      : ValueExtractor.h
Author    : L. Shawn Matott and Vijaykumar Raghavan
Copyright : 2003, L. Shawn Matott and Vijaykumar Raghavan

This class uses the instructions in the observation file to read an output file
of the model program.

Version History
06-12-03    lsm   added copyright information and initial comments.
08-20-03    lsm   created version history field and updated comments.
03-05-04    lsm   added PSO, fixed filename parse error in obs. parse
                  added support for alternative tokens
08-11-04    lsm   made into a linked list to reduce file I/O
******************************************************************************/
#ifndef VALUE_EXTRACTOR_H
#define VALUE_EXTRACTOR_H

#include "MyHeaderInc.h"

/******************************************************************************
class ValueExtractor
******************************************************************************/
class ValueExtractor
{
   public:
      ValueExtractor(IroncladString file, bool bQuitOnErr, double errVal);
      ~ValueExtractor(void){ DBG_PRINT("ValueExtractor::DTOR"); Destroy();}
      void Destroy(void);

      void   Insert(IroncladString name);
      void   ReadOutputFiles(void);
      bool ExtractValue(IroncladString name, IroncladString search, int line, 
                          int col, char tok, double * val);

   private:
      void FileToString(void);
      bool ExtractValue(IroncladString search, int line, int col, char tok, double * val);
      ValueExtractor * GetByName(IroncladString name);
      ValueExtractor * GetNext(void)     { return m_pNxt;}            
      IroncladString   GetName(void)     { return m_FileName;}
      void SetNext(ValueExtractor * pNxt){ delete m_pNxt; m_pNxt = pNxt;}

      StringType m_FileName;

      StringType m_DataStr;
      int m_DataSize;
      bool m_bQuitOnError;
      double m_ErrorVal;

      ValueExtractor * m_pNxt;
}; /* end class ValueExtractor */
  
#endif /* VALUE_EXTRACTOR_H */

