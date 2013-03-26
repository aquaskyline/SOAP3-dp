//
//    lookupBuilder.h
//
//    SOAP2 / 2BWT
//
//    Copyright (C) 2012, HKU
//
//    This program is free software; you can redistribute it and/or
//    modify it under the terms of the GNU General Public License
//    as published by the Free Software Foundation; either version 2
//    of the License, or (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program; if not, write to the Free Software
//    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
//
///////////////////////////////////////////////////////////////////////////////////////////////

#ifndef __LOOKUP_BUILDER_H__
#define __LOOKUP_BUILDER_H__

#include <assert.h>
#include <errno.h>
#include <fcntl.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>
#include "../2bwt-lib/BWT.h"
#include "../2bwt-lib/HSP.h"
#include "LT.h"
#define Kibi *1024
#define Mibi *1024*1024
#define Gibi *1024ll*1024*1024
#ifndef MAX
#define MAX(a, b) ((a)<(b)?(b):(a))
#endif
#ifndef MIN
#define MIN(a, b) ((a)<(b)?(a):(b))
#endif
#define TRY(a) (assert((a) == -1 ? (printf("Error: [%d] %s\n", errno, strerror(errno)), 0) : 1))
#define TRYEQ(a, b) (assert((a) != b ? (printf("Error: [%d] %s\n", errno, strerror(errno)), 0) : 1))
typedef unsigned long long ULL;

int BuildLookupTable(char * PackedDNAFileName, char * LookupTableFileName, int lookupTableSize);

#endif

