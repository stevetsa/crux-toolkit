#ifndef WINCRUX_H
#define WINCRUX_H
#ifdef WIN32
#include <io.h>
#include <time.h>
#include <windows.h>

// Rename some functions to the windows version
#define pclose _pclose
#define popen _popen
#define random rand
#define srandom srand

#define R_OK 04
#define W_OK 02
#define F_OK 00

int gettimeofday(struct timeval *tv, struct timezone *tz);
#endif
#endif