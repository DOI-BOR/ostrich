/******************************************************************************
File      : FortranSupportUtilities.cpp
Author    : L. Shawn Matott
Copyright : 2011, L. Shawn Matott

C-style routines to help with various pre- and post-processing tasks.

Version History
03-17-11    lsm   created
******************************************************************************/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

#include "Exception.h"
#include "FortranSupportUtilities.h"
#include "Utility.h"

/******************************************************************************
TestSupportUtilities()
  
A stub for testing the utlity code.
******************************************************************************/
int TestSupportUtilities(void)
{
   double val;

   ExtractParameter("____DUR____", "input.tpl", "input.txt", false, &val, NULL, 0);
   printf("____DUR____ = ");
   WriteFixedFormat(stdout, val, "F11.2");
   printf("\n");

   ExtractParameter("____VEL____", "input.tpl", "input.txt", false, &val, NULL, 0);
   printf("____VEL____ = ");
   WriteFixedFormat(stdout, val, "E11.4");
   printf("\n");

   ExtractParameter("___F__", "input.tpl", "input.txt", true, &val, NULL, 0);
   printf("___F__ = ");
   WriteFixedFormat(stdout, val, "F6.4");
   printf("\n");

   ExtractParameter("__A_", "input.tpl", "input.txt", true, &val, NULL, 0);
   printf("__A_ = ");
   WriteFixedFormat(stdout, val, "F4.1");
   printf("\n");

   ExtractParameter("__P_", "input.tpl", "input.txt", true, &val, NULL, 0);
   printf("__P_ = ");
   WriteFixedFormat(stdout, val, "F4.2");
   printf("\n");

   ExtractParameter("___H__", "input.tpl", "input.txt", true, &val, NULL, 0);
   printf("___H__ = ");
   WriteFixedFormat(stdout, val, "F6.4");
   printf("\n");

   return 0;
}/* end TestSupportUtilities() */

/******************************************************************************
ExtractParameter()
  
Read initial value of a parameter from a model input file, using the 
corresponding template file as a guide.

returns true if successful, in which case the value is stored in val
******************************************************************************/
bool ExtractParameter(const char * name, const char * tpl, const char * inp, bool bFixed, double * val, char ** name_list, int num_params)
{
   bool bFound = false, bAdvanced = false;
   char * pTpl, * pInp, * pTok;
   char * tpl_line;
   char * inp_line;
   char pVal[DEF_STR_SZ];
   int p;

   FILE * pTPL = fopen(tpl, "r");
   if(pTPL == NULL) return false;

   FILE * pINP = fopen(inp, "r");
   if(pINP == NULL){ fclose(pTPL); return false;}

   //determine maximum line size in each file and allocate space accordinagly
   int max_tpl_line = 2*(GetMaxLineSize(pTPL)+2*num_params);
   int max_inp_line = 2*GetMaxLineSize(pINP);
   tpl_line = new char[max_tpl_line];
   inp_line = new char[max_inp_line];
  
   //search line by line for parameter name
   while((!bFound) && (!feof(pINP)) && (!feof(pTPL)))
   {
      fgets(tpl_line, max_tpl_line, pTPL);
      fgets(inp_line, max_inp_line, pINP);

      if(strcmp(tpl_line, inp_line) != 0)
      {
         // protect parameter names
         if(num_params > 0)
         {
            for(p = 0; p < num_params; p++)
            {
               MyStrProtect(tpl_line, name_list[p]);
            }
         }
         else
         {
            MyStrProtect(tpl_line, name);
         }

         //if parameter name found, extract value from input file
         if(strstr(tpl_line, name) != NULL)
         {
            MyStrDiff(tpl_line, inp_line);

            // unprotect parameter names
            if(num_params > 0)
            {
               for(p = 0; p < num_params; p++)
               {
                  MyStrUnProtect(tpl_line, name_list[p]);
               }
            }
            else
            {
               MyStrUnProtect(tpl_line, name);
            }

            //work our way through the tpl and inp lines until next parameter matches the target
            pTpl = tpl_line;
            pInp = inp_line;
            while(strncmp(pTpl, name, strlen(name)) != 0)
            {
               bAdvanced = false;
               for(p = 0; p < num_params; p++)
               {
                  if(strncmp(pTpl, name_list[p], strlen(name_list[p])) == 0)
                  {
                     //advance tpl line past parameter
                     pTpl += strlen(name_list[p]);
                     //advance inp line past value
                     while(IsNumeric(*pInp) == true)
                     {
                        pInp++;
                        if(*pInp == NULLSTR)
                        {
                           break;
                        }
                     }
                     //advance inp and tpl past matching text to position at next parameter
                     while(*pInp == *pTpl)
                     {
                        pInp++;
                        pTpl++;
                        if((*pInp == NULLSTR) || (*pTpl == NULLSTR))
                        {
                           break;
                        }
                     }/* end while() */
                     MyTrim(pInp);
                     MyTrim(pTpl);
                     bAdvanced = true;
                     break; /* break out of for() loop */
                  }/* end if() */
               }/* end for()*/
               // sanity check, did we somehow make it through all the parameters without finding a match?
               if(bAdvanced == false)
               {
                  break;
               }
            }/* end while() */

            //did we align tpl and inp at desired parameter?
            if(strncmp(pTpl, name, strlen(name)) == 0)
            {
               //found match, next value in pInp is the value we want
               pTok = pInp;
               while(IsNumeric(*pTok) == true)
               {
                  pTok++;
                  if(*pTok == NULLSTR)
                  {
                     break;
                  }
               }
               *pTok = NULLSTR;
               strcpy(pVal, pInp);

               *val = atof(pVal);
               bFound = true;
            }/* end if() */
         }/* end if() */
      }/* end if(tpl_line != inp_line) */
   }/* end while() */

   fclose(pTPL);
   fclose(pINP);
   delete [] tpl_line;
   delete [] inp_line;
   return bFound;
}/* end ExtractParameter() */

/******************************************************************************
GetMaxLineSize()
  
Determine largest line in a file.
******************************************************************************/
int GetMaxLineSize(FILE * pFile)
{
   int cur;
   int max = 0;
   while(!feof(pFile))
   {
      cur = 0;
      while(fgetc(pFile) != '\n')
      {
         cur++;
         if(feof(pFile))
         {
            break;
         }
      }
      if(cur > max)
      {
         max = cur;
      }
   }/* end while() */
   rewind(pFile);
   return max;
}/* end GetMaxLineSize() */

/******************************************************************************
WriteFixedFormat()
  
Write the parameter value to a file according to the specified fixed format.
******************************************************************************/
bool WriteFixedFormat(FILE * pFile, double val, const char * fmt)
{
   char t = toupper(fmt[0]);
   int w;
   int d;
   if(pFile == NULL){ return false;}

   switch(t)
   {
      case('F') :
      {
         if(strlen(fmt) < 4) { fprintf(pFile, "%f", val); return false;}
         if(strstr(fmt, ".") == NULL) { fprintf(pFile, "%f", val); return false;}
         sscanf(&(fmt[1]), "%d.%d", &w, &d);
         fprintf(pFile, "%*.*f", w,d,val);
         break;
      }/* end case('F') */
      case('E') :
      case('D') :
      {
         if(strlen(fmt) < 4) { fprintf(pFile, "%E", val); return false;}
         if(strstr(fmt, ".") == NULL) { fprintf(pFile, "%E", val); return false;}
         sscanf(&(fmt[1]), "%d.%d", &w, &d);
         fprintf(pFile, "%*.*E", w,d,val);
         break;
      }/* end case('E/D') */
      case('I') :
      {
         if(strlen(fmt) < 2) { fprintf(pFile, "%d", (int)val); return false;}
         sscanf(&(fmt[1]), "%d", &w);
         fprintf(pFile, "%*d", w, (int)val);
         break;
      }/* end case('I') */
      default :
      {
         fprintf(pFile, "%E", val);
         return false;
      }
   }/* end switch() */

   return true;
}/* end WriteFixedFormat() */

/******************************************************************************
GetFixedFormatValAsStr()
  
Write the parameter value to a string according to the specified fixed format.
******************************************************************************/
bool GetFixedFormatValAsStr(char * valStr, double val, char * fmt)
{
   if(fmt == NULL){ return false;}
   char t = toupper(fmt[0]);
   int w;
   int d;
   int offset;
   if(valStr == NULL){ return false;}

   static bool bFmtIsSet = false;

   if(!bFmtIsSet)
   {
      bFmtIsSet = true;
   }

   switch(t)
   {
      case('F') :
      {
         if(strlen(fmt) < 4) { sprintf(valStr, "%f", val); return false;}
         if(strstr(fmt, ".") == NULL) { sprintf(valStr, "%f", val); return false;}
         sscanf(&(fmt[1]), "%d.%d", &w, &d);
         sprintf(valStr, "%*.*f", w,d,val);         
         break;
      }/* end case('F') */
      case('E') :
      case('D') :
      {
         if(strlen(fmt) < 4) { sprintf(valStr, "%E", val); return false;}
         if(strstr(fmt, ".") == NULL) { sprintf(valStr, "%E", val); return false;}
         sscanf(&(fmt[1]), "%d.%d", &w, &d);
         if(val < 0.00) offset = 8;
         else offset = 7;
         if((d+offset) > w) d = w - offset;
         if(d == -1) d = 0; //remove decimal place
         if(d < 0)
         {
            char msg[200];
            sprintf(msg, "Parameter value %E does not fit into fixed format of %s\n", val, fmt);
            LogError(ERR_FIXD_FMT, msg);
            ExitProgram(1);
         }
         sprintf(valStr, "%*.*E", w,d,val);
         break;
      }/* end case('E/D') */
      case('I') :
      {
         if(strlen(fmt) < 2) { sprintf(valStr, "%d", (int)val); return false;}
         sscanf(&(fmt[1]), "%d", &w);
         sprintf(valStr, "%*d", w, (int)val);
         break;
      }/* end case('I') */
      default :
      {
         sprintf(valStr, "%E", val);
         return false;
      }
   }/* end switch() */

   valStr[w]=NULLSTR; //enforce maximum width by truncating final value
   return true;
}/* end GetFixedFormatValAsStr() */
