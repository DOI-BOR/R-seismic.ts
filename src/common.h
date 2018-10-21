#pragma once

#if ! defined(__GNUC__) && (defined(_WIN32) || defined (_WIN64) || defined(_WIN))
/* defs for Windows (but not if MINGW, CYGWIN, or other GNUC-based compiler) */
# include <windows.h>
# include <io.h>
# include <direct.h>
# define strcasecmp _stricmp
# define strncasecmp _strnicmp
# define popen _popen
# define pclose _pclose
# define fstat _fstat
# define stat _stat
# define isatty _isatty	/* deprecated POSIX usage */
# define creat _creat
# define unlink _unlink
# define mkdir _mkdir
# define write _write
# define close _close
# define fileno _fileno
# define swab _swab
# define cabs _cabs
# define _USE_MATH_DEFINES
# ifndef DllExport
#  define DllExport	__declspec( dllexport )
# endif
# ifndef DllImport
#  define DllImport	__declspec( dllimport )
# endif
# ifndef WIN32
#  define WIN32
# endif
#else
# include <unistd.h>
# include <sys/file.h>
# include <sys/time.h>
# include <dirent.h>
# ifndef DllExport
#  define DllExport
# endif
# ifndef DllImport
#  define DllImport	extern
# endif
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>

/* define BOOL type */
typedef int BOOL;
#ifndef TRUE
# define TRUE (1)
#endif
#ifndef FALSE
# define FALSE (0)
#endif

#if defined(__i386) || defined(__x86_64) || defined(__amd64)
# define MACHINE_IS_BIG_ENDIAN	FALSE
#else
# define MACHINE_IS_BIG_ENDIAN	TRUE
#endif

#ifndef RAD2DEG
# define RAD2DEG (180. / M_PI)
#endif
#ifndef DEG2RAD
# define DEG2RAD (M_PI / 180.)
#endif
#ifndef DEG2KM
# define DEG2KM (2 * M_PI * 6378.163 / 360)
#endif

#if defined(_WIN32) || defined (_WIN64) || defined(_WIN)
# define PATH_SEP ('\\')
#else
# define PATH_SEP ('/')
#endif
#ifndef BASENAME
# define BASENAME(s)	(strrchr((s),PATH_SEP) == NULL ? (s) : strrchr((s),PATH_SEP)+1)
#endif
#ifndef LASTCHAR
# define LASTCHAR(s)	((s)[strlen(s)-1])
#endif

#ifndef MIN
# define MIN(a,b)	((a) <= (b) ? (a) : (b))
#endif
#ifndef MAX
# define MAX(a,b)	((a) >= (b) ? (a) : (b))
#endif
#ifndef ABS
# define ABS(x)	((x) < 0 ? -(x) : (x))
#endif
#ifndef SIGN
# define SIGN(x)	(((x) < 0) ? -1 : 1)
#endif
#ifndef CHS
# define CHS(x)	(x) = (-x)
#endif
#ifndef ROUND
# define ROUND(x)	((int)((x) + ((x) < 0 ? -.5 : .5)))
#endif

#ifdef SIGN
# undef SIGN
# define SIGN(a,b)	((b) < 0 ? -ABS(a) : ABS(a))
#endif
