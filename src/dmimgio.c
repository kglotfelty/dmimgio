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

#include "dmimgio.h"

#include <dsnan.h>
#include <stdlib.h>
#include <string.h>
#include <cxcregion.h>
#include <liberr.h>
#include <float.h>


/*
 *  Load images using dmimgio routines
 */
Image* load_image(char *infile)
{
    // Load image

    Image *image;
    if (NULL == (image = calloc(1,sizeof(Image)))) {
        err_msg("ERROR: Cannot allocate memory for image\n");
        return(NULL);
    }

    if (NULL == (image->block = dmImageOpen(infile))) {
        err_msg("ERROR: Cannot load infile '%s'\n",infile);
        return(NULL);
    }

    // dmimgio
    regRegion *dss = NULL;
    NullValue null_val;
    short has_null;
    image->dt = get_image_data(image->block, &(image->data),
                    &(image->lAxes), &dss, &null_val, &has_null);
    if (dmUNKNOWNTYPE == image->dt) {
        return NULL;
    }
    get_image_wcs(image->block, &(image->xdesc), &(image->ydesc));
    image->mask = get_image_mask(image->block, image->data,
                    image->dt, image->lAxes, dss, null_val,
                    has_null, image->xdesc, image->ydesc);

    if (dss != NULL){
        regFree(dss);
        dss=NULL;
    }

    image->hdr = getHdr(image->block, hdrDM_FILE);

    return(image);
}







short *get_image_mask( dmBlock *inBlock, void *data, dmDataType dt, 
                       long *lAxes, regRegion *dss, NullValue null, short has_null, 
                       dmDescriptor *xAxis, dmDescriptor *yAxis )
{
  long npix = lAxes[0] * lAxes[1];
  short *mask;
  long xx, yy;
  mask = (short*)calloc( npix, sizeof(short));

  double xmin, xmax, ymin, ymax;
  
  xmin = ymin = -DBL_MAX;
  xmax = ymax = DBL_MAX;
  
  if (xAxis) {
    if (yAxis) {
        dmDescriptorGetRange_d( xAxis, &xmin, &xmax);
        dmDescriptorGetRange_d( yAxis, &ymin, &ymax);        
    } else {
        dmDescriptorGetRange_d( dmGetCpt(xAxis, 1), &xmin, &xmax);
        dmDescriptorGetRange_d( dmGetCpt(xAxis, 2), &ymin, &ymax);        
    }
  } 


  dmDescriptor *x_dss, *y_dss;
  double *x_lo, *x_hi;
  long x_num;
  double *y_lo, *y_hi;
  long y_num;
  char xname[100];
  char yname[100];

  
  if (yAxis) {
        dmGetName(xAxis, xname, 99);
        dmGetName(yAxis, yname, 99);
  } else {
        dmGetName(dmGetCpt(xAxis,1), xname, 99);
        dmGetName(dmGetCpt(xAxis,2), yname, 99);
  }
  x_dss = dmSubspaceColOpen(inBlock,xname );
  y_dss = dmSubspaceColOpen(inBlock,yname );
  
  if (x_dss) {
      dmSubspaceColGet_d(x_dss, &x_lo, &x_hi, &x_num);
  } else {
      x_num = 0;
      x_lo = x_hi = NULL;
  }

  if (y_dss) {
      dmSubspaceColGet_d(y_dss, &y_lo, &y_hi, &y_num);
  } else {
      y_num = 0;
      y_lo = y_hi = NULL;
  }



  
  for ( xx=lAxes[0]; xx--; ) {
    for ( yy=lAxes[1]; yy--; ) {
      double dat;
      long idx;
      idx = xx + ( yy * lAxes[0] );
      
      dat = get_image_value( data, dt, xx, yy, lAxes, NULL );
      
      /*
       * Now check if there is a NULL value
       */ 

      if ds_dNAN(dat) {
          continue;
      }

      if (has_null) {
          switch ( dt ) {
              case dmBYTE: {
                if ( dat == null.null_byte) continue;
                break;
              }
                
              case dmSHORT: {
                if ( dat == null.null_short) continue;
                break;
              }
                
              case dmUSHORT: {
                if ( dat == null.null_ushort) continue;
                break;
              }
                
              case dmLONG: {
                if ( dat == null.null_long) continue;
                break;
              }
                
              case dmULONG: {
                if ( dat == null.null_ulong) continue;
                break;
              }
                
              case dmFLOAT: {
                if ( dat == null.null_float) continue;
                break;
              }
              case dmDOUBLE: {
                if ( dat == null.null_double) continue;
                break;
              }
              default:
                break;
          
        } // end switch
          
      } // end has_null


            
      /* If the image has a data sub space (aka a region filter applied)
         then need to convert coords to physical and check */
      if (xAxis ) {
        double pos[2];
        double loc[2];
        pos[0]=xx+1;
        pos[1]=yy+1;
        
        if (yAxis) {  /* If no y axis, then xAxis has 2 components */
          dmCoordCalc_d( xAxis, pos, loc );
          dmCoordCalc_d( yAxis, pos+1, loc+1 );
        } else {
          dmCoordCalc_d( xAxis, pos, loc );
        }


        // Check region subspace
        if ( dss && !regInsideRegion( dss, loc[0], loc[1] ) )
          continue;


        // Check range
        if ((loc[0] < xmin) || (loc[0]>xmax)) continue;
        if ((loc[1] < ymin) || (loc[1]>ymax)) continue;

        // Check subspace
        int ii;
        short good = x_num > 0 ? 0 : 1;  // if no DSS, then good
        for (ii=0;ii<x_num;ii++) {
            if ((x_lo[ii] <= loc[0]) && (loc[0] <= x_hi[ii])) good=1;            
        }
        if (!good) 
            continue;

        good = y_num > 0 ? 0 : 1;   // if no DSS then good
        for (ii=0;ii<y_num;ii++) {
            if ((y_lo[ii] <= loc[1]) && (loc[1] <= y_hi[ii])) good=1;            
        }            
        if (!good)
            continue;


       }   // end if xAxis
      

      
      mask[idx] = 1;  // Valid Pixel, inside mask
    }
  }


  return(mask );
}


double get_image_value( void *data, dmDataType dt, 
                        long xx, long yy, long *lAxes, 
                        short *mask )
{

  long npix = xx + (yy * lAxes[0] );
  double retval;

  /* Okay, first get all the data from the different data types.  
     Cast everything to doubles */


  if (( xx < 0 ) || ( xx >= lAxes[0] ) ||
      ( yy < 0 ) || ( yy >= lAxes[1] )) {
    return(0);
  }


  switch ( dt ) {
    
  case dmBYTE: {
    unsigned char *img = (unsigned char*)data;
    retval = img[npix];
    break;
  }
    
  case dmSHORT: {
    short *img = (short*)data;
    retval = img[npix];
    break;
  }
    
  case dmUSHORT: {
    unsigned short *img = (unsigned short*)data;
    retval = img[npix];
    break;
  }
    
  case dmLONG: {
    long *img = (long*)data;
    retval = img[npix];
    break;
  }
    
  case dmULONG: {
    unsigned long *img = (unsigned long*)data;
    retval = img[npix];
    break;
  }
    
  case dmFLOAT: {
    float *img = (float*)data;
    retval = img[npix];
    break;
  }
  case dmDOUBLE: {
    double *img = (double*)data;
    retval = img[npix];
    break;
  }
  default:
    ds_MAKE_DNAN( retval );

  }


  if ( mask ) {
    if ( !mask[npix] ) {
      ds_MAKE_DNAN( retval );
    }
  }


  return(retval);

}


/* Load the data into memory,  check for DSS, null values */
dmDataType get_image_data( dmBlock *inBlock, void **data, long **lAxes,
                           regRegion **dss, NullValue *nullval, short *nullset )
{

  dmDescriptor *imgDesc;
  dmDataType dt;
  dmDescriptor *grp;
  dmDescriptor *imgdss;

  long naxes;
  long npix;
  char ems[1000];


  *dss = NULL;
  *nullset = 0;
  
  imgDesc = dmImageGetDataDescriptor( inBlock );

  /* Sanity check, only 2D images */
  naxes = dmGetArrayDimensions( imgDesc, lAxes );
  if ( naxes != 2 ) {
    err_msg("ERROR: Only 2D images are supported");
    return( dmUNKNOWNTYPE );
  }
  npix = (*lAxes)[0] * (*lAxes)[1];
  dt = dmGetDataType( imgDesc );


  /* Okay, first lets get the image descriptor */
  grp = dmArrayGetAxisGroup( imgDesc, 1 );
  dmGetName( grp, ems, 1000);
  imgdss = dmSubspaceColOpen( inBlock, ems );
  if ( imgdss )
    *dss = dmSubspaceColGetRegion( imgdss);
  
  
  switch ( dt ) 
    {
    case dmBYTE:
      *data = ( void *)calloc( npix, sizeof(char ));
      dmGetArray_ub( imgDesc, (unsigned char*) *data, npix );
      if ( dmDescriptorGetNull_ub( imgDesc, &(nullval->null_byte)) == 0 ) {
        *nullset=0;
      } else
        *nullset=1;
      break;
      
    case dmSHORT:
      *data = ( void *)calloc( npix, sizeof(short ));
      dmGetArray_s( imgDesc, (short*) *data, npix );
      if ( dmDescriptorGetNull_s( imgDesc, &(nullval->null_short)) == 0 ) {
        *nullset=0;
      } else
        *nullset=1;
      break;
      
    case dmUSHORT:
      *data = ( void *)calloc( npix, sizeof(short ));
      dmGetArray_us( imgDesc, (unsigned short*) *data, npix );
      if ( dmDescriptorGetNull_us( imgDesc, &(nullval->null_ushort)) == 0 ) {
        *nullset=0;
      } else
        *nullset=1;
      break;
      
    case dmLONG:
      *data = ( void *)calloc( npix, sizeof(long ));
      dmGetArray_l( imgDesc, (long*) *data, npix );
      if ( dmDescriptorGetNull_l( imgDesc, &(nullval->null_long)) == 0 ) {
        *nullset=0;
      } else
        *nullset=1;
      break;
      
    case dmULONG:
      *data = ( void *)calloc( npix, sizeof(long ));
      dmGetArray_ul( imgDesc, (unsigned long*) *data, npix );
      if ( dmDescriptorGetNull_ul( imgDesc, &(nullval->null_ulong)) == 0 ) {
        *nullset=0;
      } else
        *nullset=1;
      break;
      
    case dmFLOAT:
      *data = ( void *)calloc( npix, sizeof(float ));
      dmGetArray_f( imgDesc, (float*) *data, npix );
      if ( dmDescriptorGetNull_f( imgDesc, &(nullval->null_float)) == 0 ) {
        *nullset=0;
      } else
        *nullset=1;
      break;
      
    case dmDOUBLE:
      *data = ( void *)calloc( npix, sizeof(double ));
      dmGetArray_d( imgDesc, (double*) *data, npix );
      if ( dmDescriptorGetNull_d( imgDesc, &(nullval->null_double)) == 0 ) {
        *nullset=0;
      } else
        *nullset=1;
      break;
      
    default:
      err_msg("ERROR: Unsupported datatype");
      return( dmUNKNOWNTYPE );
    }

  return(dt);

}




/* Get the WCS descriptor */
short  get_image_wcs( dmBlock *imgBlock, dmDescriptor **xAxis, 
                      dmDescriptor **yAxis )
{
  

  dmDescriptor *imgData;
  long n_axis_groups;

  imgData = dmImageGetDataDescriptor( imgBlock );
  n_axis_groups = dmArrayGetNoAxisGroups( imgData );
  

  /* This is the usual trick ... can have 1 axis group w/ 
     dimensionality 2 (eg a vector column) or can have
     2 axis groups w/ dimensionaity 1 (eg 2 disjoint columns)*/
  if ( n_axis_groups == 1 ) {
    dmDescriptor *pos = dmArrayGetAxisGroup( imgData, 1 );
    dmDescriptor *xcol;
    long n_components;
    
    n_components = dmGetElementDim( pos );
    if ( n_components != 2 ) {
      err_msg("ERROR: could not find 2D image\n");
      return(-1);
    }
    
    xcol = dmGetCpt( pos, 1 );
    
    *xAxis = pos;
    *yAxis = NULL;
    
  } else if ( n_axis_groups == 2 ) {
    dmDescriptor *xcol;
    dmDescriptor *ycol;
  
    xcol = dmArrayGetAxisGroup( imgData, 1 );
    ycol = dmArrayGetAxisGroup( imgData, 2 );

    *xAxis = xcol;
    *yAxis = ycol;
    
  } else {
    err_msg("Invalid number of axis groups\n");
    *xAxis = NULL;
    *yAxis = NULL;
    return(-1);
  }

  return(0);

}
