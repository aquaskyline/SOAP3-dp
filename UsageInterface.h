/*
 *
 *    UsageInterface.c
 *    Soap3(gpu)
 *
 *    Copyright (C) 2011, HKU
 *
 *    This program is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU General Public License
 *    as published by the Free Software Foundation; either version 2
 *    of the License, or (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 *
 */

#ifndef __USAGE_INTERFACE_H__
#define __USAGE_INTERFACE_H__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

#include "Release.h"

void UIPrintUsageOverview ( char * program_name );
void UIPrintAlignmentOptions ( int casenum );
void UIprintUsageSingle ( char * program_name );
void UIprintUsagePair ( char * program_name );
void UIprintUsageSingleList ( char * program_name );
void UIprintUsagePairList ( char * program_name );


#endif
