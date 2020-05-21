/******************************************************************************
File     : ConversionTool.h
Author   : Tsu-Wei Webber Chen
Copyright: 2010, Tsu-Wei Webber Chen

ConversionTool is a base class which defines the set of common interface used to
convert different file types into plain ASCII that is readable by the Ostrich
program.

Version History
01-25-10    twc   Created file.
******************************************************************************/
#ifndef CONVERSION_TOOL_H
#define CONVERSION_TOOL_H

/******************************************************************************
class ConversionTool
******************************************************************************/
class ConversionTool
{
   public:
	  virtual void Convert(void)=0;
}; /* end class ConversionTool */

#endif /* CONVERSION_TOOL_H */