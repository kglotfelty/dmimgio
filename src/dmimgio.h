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

short *get_image_mask( dmBlock *inBlock, void *data, dmDataType dt, 
                       long *lAxes, regRegion *dss, long null, short has_null, 
                       dmDescriptor *xAxis, dmDescriptor *yAxis );

double get_image_value( void *data, dmDataType dt, 
                        long xx, long yy, long *lAxes, 
                        short *mask );

dmDataType get_image_data( dmBlock *inBlock, void **data, long **lAxes,
                           regRegion **dss, long *nullval, short *nullset );

short  get_image_wcs( dmBlock *imgBlock, dmDescriptor **xAxis, 
                      dmDescriptor **yAxis );


