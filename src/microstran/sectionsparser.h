/*	FRAME3DD: Static and dynamic structural analysis of 2D & 3D frames and trusses
	Copyright (C) 2007-2008 John Pye

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*//** @FILE
	parser for microstran properties.txt files
*/

#ifndef SECTIONSPARSER_H
#define SECTIONSPARSER_H

#include "config.h"
#include "sections.h"
#include "parse.h"

#ifdef __cplusplus
extern "C"{
#endif

MSTRANP_API cbool parseSections(parse *p, section_library *a);

#ifdef __cplusplus
}
#endif

#endif

