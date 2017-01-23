#include "common.h"

#include <R.h>
#include <Rinternals.h>

#ifdef USING_R
/* OOPS - print a string and exit */
void oops(char *name, char *msg)
{
	error("%s: %s\n",name, msg);
}
void smsg(char *name, char *msg)
{
	Rprintf("%s: %s\n", name, msg);
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
	fprintf(stderr,"%s: %s\n", name, msg);
}
#endif


