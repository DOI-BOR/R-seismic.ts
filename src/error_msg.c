#include "common.h"

char msgbuf[2048];

#ifdef USING_R
# include <R.h>
# include <Rinternals.h>
/* OOPS - print a string and exit */
void oops(char *name, char *msg)
{
  error("%s: %s\n",name, msg);
}
void smsg(char *name, char *msg)
{
  char *fmt = isalnum(*name) ? "%s: %s\n" : "%s%s\n";
  Rprintf(fmt, name, msg);
}
#else
/* OOPS - print a string and exit */
void oops(char *name, char *msg)
{
  fprintf(stderr,"%s: %s\n",name, msg);
  exit(1);
}
void smsg(char *name, char *msg)
{
  char *fmt = isalnum(*name) ? "%s: %s\n" : "%s%s\n";
  fprintf(stderr, fmt, name, msg);
  fflush(stderr);
}
#endif
