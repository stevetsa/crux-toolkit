/**
 * \file WinCrux.cpp
 * \brief Support functions need to compile Crux under native windows.
 */

#include <stdlib.h>
#include "carp.h"
#include "utils.h"
#include "WinCrux.h"

char *realpath(const char * file_name, char * resolved_name) {

  char * full_path_buffer = (char *) mymalloc(MAX_PATH * sizeof(char));
  size_t needed_buff_size 
    = GetFullPathName(file_name, MAX_PATH, full_path_buffer, NULL);

  if (needed_buff_size == 0) {
    // An error occurred.
    LPSTR lpMsgBuf;
    FormatMessage(
      FORMAT_MESSAGE_ALLOCATE_BUFFER 
      | FORMAT_MESSAGE_FROM_SYSTEM 
      | FORMAT_MESSAGE_IGNORE_INSERTS,
      NULL,
      GetLastError(),
      MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT), // Default language
      (LPTSTR) &lpMsgBuf,
      0,
      NULL 
    );
    carp(CARP_FATAL, lpMsgBuf);

  }

  if (needed_buff_size > MAX_PATH) {
    full_path_buffer = (char *) myrealloc(full_path_buffer, needed_buff_size);
    needed_buff_size 
    = GetFullPathName(file_name, MAX_PATH, full_path_buffer, NULL);
  }

  return full_path_buffer;
}

/* 
This code was placed in the public domain by the author, 
Sean Barrett, in November 2007. Do with it as you will. 
(Seee the page for stb_vorbis or the mollyrocket source 
page for a longer description of the public domain non-license). 
*/ 

#define WIN32_LEAN_AND_MEAN 
#include <windows.h> 

// Public domain code from https://mollyrocket.com/forums/viewtopic.php?p=2529
// map 'filename' and return a pointer to it. fill out *length and *un if not-NULL 
void *stub_mmap(const char *filename, SIMPLE_UNMMAP *un) 
{ 
   HANDLE f = CreateFile(filename, GENERIC_READ, FILE_SHARE_READ,  NULL, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL); 
   HANDLE m; 
   void *p; 
   if (!f) return NULL; 
   m = CreateFileMapping(f, NULL, PAGE_READONLY, 0,0, NULL); 
   if (!m) { CloseHandle(f); return NULL; } 
   p = MapViewOfFile(m, FILE_MAP_READ, 0,0,0); 
   if (!p) { CloseHandle(m); CloseHandle(f); return NULL; } 
   if (un) { 
      un->f = f; 
      un->m = m; 
      un->p = p; 
   } 
   return p; 
} 

void stub_unmmap(SIMPLE_UNMMAP *un) 
{ 
   UnmapViewOfFile(un->p); 
   CloseHandle(un->m); 
   CloseHandle(un->f); 
} 