/******************************************************************************
File     : MemoryTracker.h
Author   : L. Shawn Matott
Copyright: 2013, L. Shawn Matott

Track memory usage of sub-routines.

Version History
01-29-13    lsm   added copyright information and initial comments.
******************************************************************************/
#ifndef MEMORY_TRACKER_H
#define MEMORY_TRACKER_H

void LogMemUsage(char * file, const char * tag);
double GetMemUsage(void);

#endif /* MEMORY_TRACKER_H */