/*                                                                
**  Copyright (C) 2004-2008  Smithsonian Astrophysical Observatory 
*/                                                                

/*                                                                          */
/*  This program is free software; you can redistribute it and/or modify    */
/*  it under the terms of the GNU General Public License as published by    */
/*  the Free Software Foundation; either version 3 of the License, or       */
/*  (at your option) any later version.                                     */
/*                                                                          */
/*  This program is distributed in the hope that it will be useful,         */
/*  but WITHOUT ANY WARRANTY; without even the implied warranty of          */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           */
/*  GNU General Public License for more details.                            */
/*                                                                          */
/*  You should have received a copy of the GNU General Public License along */
/*  with this program; if not, write to the Free Software Foundation, Inc., */
/*  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.             */
/*                                                                          */

#include <ascdm.h>
#include <cxcregion.h>
#include "hdrlib2.h"


/* Hold info for an input image */
typedef struct {
    void *data;        // pixel values
    dmDataType dt;     // pixel datatype
    long *lAxes;       // axis lenghts
    short *mask;        // mask of valid pixels
    dmDescriptor *xdesc;  // X (or sky) coordinate descriptor
    dmDescriptor *ydesc;  // Y coordinate descriptor
    dmBlock *block; // The block image came from
    Header_Type *hdr;   // Header keywords
} Image;


typedef union {
    unsigned char null_byte;
    short null_short;
    unsigned short null_ushort;
    int null_int;
    long null_long;
    unsigned long null_ulong;
    float null_float;
    double null_double;
} NullValue;


Image *load_image( char *infile );


short *get_image_mask( dmBlock *inBlock, void *data, dmDataType dt, 
                       long *lAxes, regRegion *dss, NullValue null, short has_null, 
                       dmDescriptor *xAxis, dmDescriptor *yAxis );

double get_image_value( void *data, dmDataType dt, 
                        long xx, long yy, long *lAxes, 
                        short *mask );

dmDataType get_image_data( dmBlock *inBlock, void **data, long **lAxes,
                           regRegion **dss, NullValue *nullval, short *nullset );

short  get_image_wcs( dmBlock *imgBlock, dmDescriptor **xAxis, 
                      dmDescriptor **yAxis );


