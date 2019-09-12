/*
! Copyright (C) 1991-2004  ; All Rights Reserved ; Colorado State University
! Colorado State University Research Foundation ; ATMET, LLC
! 
! This file is free software; you can redistribute it and/or modify it under the
! terms of the GNU General Public License as published by the Free Software 
! Foundation; either version 2 of the License, or (at your option) any later version.
! 
! This software is distributed in the hope that it will be useful, but WITHOUT ANY 
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
! PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along with this 
! code; if not, write to the Free Software Foundation, Inc., 
! 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
!======================================================================================
*/

/*
 * Using these C wrapper routines to be able to utilize the C API for HDF5.
 * Unfortunately, the fortran .mod files that came with our installation of 
 * HDF5 do not work with the pgf90 compiler.
 */

#include "utils_sub_names.h"
#include "revu_sub_names.h"
#include <string.h>
#include "hdf5.h"

/*
 * Limits for arrays, strings, etc. Need to keep these in sync with like
 * named parameters in rhdf5_utils.f90
*/
#define RHDF5_MAX_DIMS    10

/*
 * Integer coding for HDF5 types. Need to keep these in sync with like named
 * parameters in rhdf5_utils.f90
 * These codes need to start with zero and be contiguous (0..n with no gaps)
 * for the rhdf5_type_names array. Make sure the entries in rhdf5_type_names
 * are consistent with the type numbers.
 */
#define RHDF5_NUM_TYPES     4

#define RHDF5_TYPE_STRING   0
#define RHDF5_TYPE_INTEGER  1
#define RHDF5_TYPE_FLOAT    2
#define RHDF5_TYPE_CHAR     3

char *rhdf5_type_names[RHDF5_NUM_TYPES] = {"STRING","INTEGER","FLOAT","CHAR"};

/* Prototypes for internal routines */

 void rf2c_get_mem_type (hid_t typid,int *type,int *size);
hid_t rf2c_set_mem_type (int dtype,int ssize);

/* Routines for HDF5 REVU IO. */

/********************** FILE ROUTINES ****************************************/

/*****************************************************************************/
 void rh5f_open (char *fname,int *f_flgs,int *fileid,int *hdferr)
{
unsigned flags;

if(*f_flgs == 1) flags = H5F_ACC_RDONLY;
if(*f_flgs == 2) flags = H5F_ACC_RDWR;

*fileid = H5Fopen (fname, flags, H5P_DEFAULT);

if (*fileid < 0)
  {
  *hdferr = *fileid;
  }
else
  {
  *hdferr = 0;
  }

return;
}

/*****************************************************************************/
 void rh5f_create (char *fname,int *f_flgs,int *fileid,int *hdferr)
{
unsigned flags;

if(*f_flgs == 1) flags = H5F_ACC_TRUNC;
if(*f_flgs == 2) flags = H5F_ACC_EXCL ;

*fileid = H5Fcreate (fname,flags,H5P_DEFAULT,H5P_DEFAULT);

if (*fileid < 0)
  {
  *hdferr = *fileid;
  }
else
  {
  *hdferr = 0;
  }

return;
}

/*****************************************************************************/
 void rh5f_close (int *fileid,int *hdferr)
{

*hdferr = H5Fclose (*fileid);

return;
}

/**************************** GROUP ROUTINES *********************************/

/*****************************************************************************/
 void rh5g_create (int *id,char *name,int *gid,int *hdferr)
{

 /* HDF5 1.8 API
 *  3rd arg is link creation property list
 *  4th arg is group creation property list
 *  5th art is group access property list (must be H5P_DEFAULT as of 4/19/12)
 */
*gid = H5Gcreate (*id,name,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);

if (*gid < 0)
  {
  *hdferr = *gid;
  }
else
  {
  *hdferr = 0;
  }
 
return;
}

/*****************************************************************************/
 void rh5g_open (int *id,char *name,int *gid,int *hdferr)
{

 /* HDF5 1.8 API
  * 3rd arg is group access property list
  */
*gid = H5Gopen (*id,name,H5P_DEFAULT);

if (*gid < 0)
  {
  *hdferr = *gid;
  }
else
  {
  *hdferr = 0;
  }
 
return;
}

/*****************************************************************************/
 void rh5g_close (int *grpid,int *hdferr)
{

*hdferr = H5Gclose (*grpid);

return;
}

/***************** LINK (GROUP, DATASET, etc) ROUTINES ***********************/

/*****************************************************************************/
 void rh5l_exists (int *id,char *name,int *exists)
{

/*
 * Third arg is "link access property list id"
 */
*exists = H5Lexists (*id,name,H5P_DEFAULT);

return;
}

/************************* DATASET ROUTINES **********************************/

/*****************************************************************************/
 void rh5d_open (int *id,char *name,int *dsetid,int *hdferr)
{

/* HDF5 1.8 API */
*dsetid = H5Dopen (*id,name,H5P_DEFAULT);

if (*dsetid < 0)
  {
  *hdferr = *dsetid;
  }
else
  {
  *hdferr = 0;
  }

return;
}

/*****************************************************************************/
 void rh5d_setup_and_write (int *id,char *name,int *dtype,int *ssize
                         ,int *ndims,int *dims,int *ext_dims,int *chunk_dims
                         ,int *deflvl,int *dsetid,void *data,int *hdferr)
{
/*
 * This routine will either create or extend a dataset and then
 * perform the write into that dataset.
 */
int i;
int cur_ndims;
hsize_t cur_dims[RHDF5_MAX_DIMS]    = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
hsize_t max_dims[RHDF5_MAX_DIMS]    = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
hsize_t new_dims[RHDF5_MAX_DIMS]    = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
hsize_t mem_dims[RHDF5_MAX_DIMS]    = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
hsize_t chunk_sizes[RHDF5_MAX_DIMS] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
hsize_t hs_count[RHDF5_MAX_DIMS]    = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
hsize_t hs_start[RHDF5_MAX_DIMS]    = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
hid_t mtypid;
hid_t mspcid;
hid_t fspcid;
hid_t dcplid;

*hdferr = 0;

/*
 * Set the memory type (translate integer code in dtype to HDF5 type)
 */
mtypid = rf2c_set_mem_type (*dtype, *ssize);
if (mtypid < 0)
  {
  printf("ERROR: rh5d_setup_and_write: unrecognized dtype: %d\n", *dtype);
  *hdferr = -1;
  return;
  }

/*
 * Check to see if the dataset exists aleady. If so extend it and do
 * the write. If not, then create it and do the write.
 */
if ( H5Lexists (*id,name,H5P_DEFAULT) )
  {
  /*
   * dataset exists --> extend and write
   *
   * Open the existing dataset
   */

  /* HDF5 1.8 API */
  *dsetid = H5Dopen (*id,name,H5P_DEFAULT);

  /*
   * Get the currently existing dimension information from the dataset
   * Close the dataspce after grabbing the current dimension information
   * since we need to reopen the dataspace after extending the dataset.
   */
  fspcid = H5Dget_space (*dsetid);
  if (fspcid < 0)
    {
    printf("ERROR: rh5d_setup_and_write: Cannot get file space id\
           from dataset\n");
    *hdferr = -1;
    return;
    }
  cur_ndims = H5Sget_simple_extent_dims (fspcid, cur_dims, max_dims);
  H5Sclose (fspcid);
  
  if (cur_ndims != *ndims)
    {
    printf("ERROR: rh5d_setup_and_write: number of dimensions in dataset\
           (%d) does not match specified number of dimensions (ndims): \
           %d\n",cur_ndims, *ndims);
    *hdferr = -1;
    return;
    }
  
  /*
   * Get the new dimension information from the input arguments
   */
  for (i = 0; i<*ndims; i++)
    {
    new_dims[i] = dims[i];
    /* hs_* vars are for the hyperslab selection later on */
    hs_count[i] = cur_dims[i];
    hs_start[i] = 0;
    if (ext_dims[i] == 1)
      {
      /* this is the dimension we are extended and assume for now */
      /* that we are extending by just one element                */
      mem_dims[i] = 1;

      if (new_dims[i] < cur_dims[i])
        {
        printf("ERROR: rh5d_setup_and_write: Shrinking an extendible \
                dimension is not supported\n");
        printf("ERROR:   Current dimension number: %d, Current dimension \
                size: %d\n", i+1, (int) cur_dims[i]);
        printf("ERROR:   New dimension number:     %d, New dimension size: \
                %d\n", i+1, (int) new_dims[i]);
        *hdferr = -1;
        return;
        }
      }
    else
      {
      mem_dims[i] = dims[i];
      }
    }
  
  /*
   * Create a memory dataspace to control the selection from the data 
   * buffer (input arg "data"). We are assuming that the entries in the
   * data buffer are organized so that the dimension that we are extending
   * is a size of one so if we set up the memory space with that in mind
   * the selection will happen automatically without the need to apply any 
   * selection mechanism on the memory space.
   */
  mspcid = H5Screate_simple (cur_ndims,mem_dims,NULL);
  
  /*
   * Extend the dataset and re-open the dataspace
   */
  H5Dset_extent (*dsetid,new_dims);
  fspcid = H5Dget_space (*dsetid);
  
  if ((mspcid < 0) || (fspcid < 0))
    {
    *hdferr = -1;
    return;
    }
  
  /*
   * Create a hyperslab to select just the region where the new data will be
   * stored into the file. That is, we are appending the new data just after
   * the end of the existing data so select out the existing data. To do this
   * first select the entire dataspace (extended) and subtract out the 
   * dataspace prior to the extension. The select operator H5S_SELECT_NOTB
   * is used to do the subtraction.
   */
  H5Sselect_all(fspcid);
  *hdferr = H5Sselect_hyperslab (fspcid,H5S_SELECT_NOTB,hs_start
                                ,NULL,hs_count,NULL);
  
  if (*hdferr != 0)
    {
    return;
    }

  /*
   * write and free up resources
   */
  *hdferr = H5Dwrite (*dsetid,mtypid,mspcid,fspcid,H5P_DEFAULT,data);

  if (*dtype == RHDF5_TYPE_STRING)
    {
    H5Tclose (mtypid);
    }
  H5Sclose (mspcid);
  H5Sclose (fspcid);
  }
else
  {
  /*
   * dataset does not exist --> create and write
   *
   * H5Screate_simple wants all of the elements of the dims and maxdim 
   * arrays (2nd, 3rd args) to be set to a non-zero value, even the elements
   * that are not being used (since they are beyond what ndims specifies).
   * Copy the used contents of dims (the first ndims entries) into an internal 
   * array which is initialized to all 1's so that the caller doesn't need
   * to bother with this detail.
   *
   * H5Pset_chunk wants elements not used in the chunk size array (2nd arg) 
   * to be set to zero and the used elements set to non-zero values.
   */
  for (i=0 ; i<*ndims; i++)
    {
    cur_dims[i] = dims[i];
    if (ext_dims[i] == 0)
      {
      max_dims[i] = dims[i];
      }
    else
      {
      max_dims[i] = H5S_UNLIMITED;
      }
    chunk_sizes[i] = chunk_dims[i];
    }
  
  /*
   * create the file dataspace
   */
  mspcid = H5Screate_simple (*ndims,cur_dims,max_dims);
  if (mspcid < 0)
    {
    printf("ERROR: rh5d_create: Cannot create memory space for dataset\n");
    *hdferr = -1;
    return;
    }
  
  /*
   * create the dataset
   * dataset creation properties
   */
  dcplid = H5Pcreate (H5P_DATASET_CREATE);
  if (dcplid < 0)
    {
    printf("ERROR: rh5d_create: Cannot make create properties for dataset\n");
    *hdferr = -1;
    return;
    }
  
  /*
   * Set up for data compression
   *   enable chunking of data
   *   enable shuffling of data
   *   enable deflation (compression) of data
   */
  H5Pset_chunk (dcplid,*ndims,chunk_sizes);
  H5Pset_shuffle (dcplid);
  H5Pset_deflate (dcplid,*deflvl);
  
  /*
   * dataset
   */

  /* HDF5 1.8 API
   * 3rd arg is type of the data elements - always using float for now
   * 5th arg is the link creation property list (eg, automatically build 
   *   intermediate groups)
   * 7th arg is the data access property list
   */
  *dsetid = H5Dcreate (*id,name,mtypid,mspcid,H5P_DEFAULT,dcplid,H5P_DEFAULT);

  if (*dsetid < 0)
    {
    printf("ERROR: rh5d_create: Cannot make create dataset\n");
    *hdferr = -1;
    return;
    }
  
  /*
   * write and free up resources
   */
  *hdferr = H5Dwrite (*dsetid,mtypid,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);

  if (*dtype == RHDF5_TYPE_STRING)
    {
    H5Tclose (mtypid);
    }
  H5Pclose (dcplid);
  H5Sclose (mspcid);
  }

return;
}

/*****************************************************************************/
 void rh5d_write (int *dsetid,int *mtypid,int *mspcid,int *fspcid
                ,int *xfplid,void *data,int *hdferr)
  {
  *hdferr = H5Dwrite (*dsetid,*mtypid,*mspcid,*fspcid,*xfplid,data);

  return;
  }

/***********************************************************************
 * This routine will open the dataset given by id (file id) and name, 
 * read in and return the dimension information.
 *
 * The dataset is left open so that the caller can read the data and
 * retrieve attributes so it is up to the caller to close the dataset.
 */
 void rh5d_read_get_dims (int *dsetid,int *ndims,int *dims,int *hdferr)
  {
  hid_t fspcid;
  hid_t dtypid;
  hsize_t dset_dims[RHDF5_MAX_DIMS];
  int i;

  *hdferr = 0;
  
  /*
   * Get the information about the dimensions
   * Need to read dimension sizes, from H5Sget_simple_extent_dims, into
   * an (hsize_t *) type since it's length is different than (int *) type.
   */
  fspcid = H5Dget_space (*dsetid);
  *ndims = H5Sget_simple_extent_dims (fspcid,dset_dims,NULL);
  for (i=0; i<*ndims; i++)
    {
    dims[i] = dset_dims[i];
    }

  H5Sclose (fspcid);
  if (*ndims < 0)
    {
    *hdferr = -1;
    }

  return;
  }

/**********************************************************************
 * This routine will open the dataset given by id (file id) and name, and
 * read the dataset into the buffer "data".
 *
 * The dataset is left open so that the caller can retrieve attributes
 * so it is up to the caller to close the dataset.
 *
 */
 void rh5d_read (int *dsetid,int *dtype,int *ms_ndims,int *ms_dims
               ,int *fs_ndims,int *fs_offset,int *fs_counts
               ,int *dsize,void *data,int *hdferr)
  {
  hid_t dtypid;
  hsize_t hs_file_counts[RHDF5_MAX_DIMS] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
  hsize_t hs_file_offset[RHDF5_MAX_DIMS] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
  hsize_t hs_mem_dims[RHDF5_MAX_DIMS]   = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
  hsize_t hs_mem_counts[RHDF5_MAX_DIMS] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
  hsize_t hs_mem_offset[RHDF5_MAX_DIMS] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
  hid_t file_dspace;
  hid_t mem_dspace;
  int i;
  int mtype;
  
  *hdferr = 0;
  
  /*
   * Get the information about the data type
   */
  dtypid = H5Dget_type (*dsetid);
  if (dtypid < 0)
    {
    *hdferr = -1;
    return;
    }

  rf2c_get_mem_type (dtypid,&mtype,dsize);
  if (mtype < 0)
    {
    *hdferr = -1;
    return;
    }
  if (*dtype != mtype)
    {
    printf("ERROR: rh5a_setup_and_read: Data type requested by caller does");
    printf("       not match type in file\n");
    printf("ERROR: Requested attribute type: %s\n", rhdf5_type_names[*dtype]);
    printf("ERROR: Data type found in HDF5 file: %s\n", rhdf5_type_names[mtype]);
    *hdferr = -1;
    return;
    }

  /*
   * set up the hyperslab selection
   * copy counts and offsets - the hdf5 interface doesn't like the argument 
   * integer pointer types for some reason
   */
  for (i=0; i<*fs_ndims; i++)
    {
    hs_file_counts[i] = fs_counts[i];
    hs_file_offset[i] = fs_offset[i];
    }

  /*
   * get the file data space
   */
  file_dspace = H5Dget_space (*dsetid);

  /*
   * file data space selection (hyperslab)
   */
  *hdferr = H5Sselect_hyperslab (file_dspace,H5S_SELECT_SET,hs_file_offset
                                ,NULL,hs_file_counts,NULL);
  if (*hdferr != 0)
    {
    return;
    }

  /*
   * memory data space selection (all)
   */
  for (i=0; i<*ms_ndims; i++)
    {
    hs_mem_dims[i] = ms_dims[i];
    hs_mem_offset[i] = 0;
    hs_mem_counts[i] = ms_dims[i];
    }
  mem_dspace = H5Screate_simple (*ms_ndims,hs_mem_dims,NULL);
  *hdferr = H5Sselect_hyperslab (mem_dspace,H5S_SELECT_SET,hs_mem_offset
                                ,NULL,hs_mem_counts,NULL);
  if (*hdferr != 0)
    {
    return;
    }

  /*
   * read the dataset, clean up and return
   */
  *hdferr = H5Dread (*dsetid,dtypid,mem_dspace,file_dspace,H5P_DEFAULT,data);

  if (*dtype == RHDF5_TYPE_STRING)
    {
    H5Tclose (dtypid);
    }

  return;
  }

/**********************************************************************/
 void rh5d_close (int *dsetid,int *hdferr)
{
*hdferr = H5Dclose (*dsetid);

return;
}

/******************* DATASPACE ROUTINES *******************/

/**********************************************************/
 void rh5s_close (int *dspcid,int *hdferr)
{
*hdferr = H5Sclose (*dspcid);

return;
}

/******************* PROPERTY ROUTINES *******************/

/**********************************************************/
 void rh5p_close (int *propid,int *hdferr)
{
*hdferr = H5Pclose (*propid);

return;
}

/****************** DATATYPE ROUTINES ******************/

/**********************************************************/
 void rh5t_close (int *dtypid,int *hdferr)
{
*hdferr = H5Tclose (*dtypid);

return;
}

/******************* ATTRIBUTE ROUTINES *******************/

/**********************************************************/
/* This routine will write a scalar attribute into the HDF5 file, given the 
 * object id and attribute name.
 *
 * value is a VOID pointer so the caller can use any type (char, int, float 
 * for now). vtype is used to denote what type value is. The encoding for 
 * vtype is contained in the the defines: RHDF5_TYPE_* (see top of this file).
 */
 void rh5a_write_anyscalar (int *id,char *name,void *value,int *vtype,int *hdferr)
  {
  hid_t a_id, amem_id, atype;
  htri_t attr_exists;
  int ssize;
  
  *hdferr = 0;
  
  /*
   * Set up memory type
   */
  if (*vtype == RHDF5_TYPE_STRING)
    {
    ssize = strlen(value) + 1;
    }
  else
    {
    ssize = 0;
    }
  
  atype = rf2c_set_mem_type (*vtype, ssize);
  
  if (atype < 0)
    {
    printf("ERROR: rh5a_write_anyscalar: unrecognized vtype: %d\n", *vtype);
    *hdferr = -1;
    return;
    }
  
  if (*hdferr == 0)
    {
    /*
     * got a recognized type so keep going
     */
  
    attr_exists = H5Aexists (*id,name);
  
    if (attr_exists)
      {
      /*
       * attribute exists already --> update mode
       */
  
      a_id = H5Aopen (*id,name,H5P_DEFAULT);
      }
    else
      {
      /*
       * attribute does not exist --> create mode
       *
       *
       * scalar dataspace to hold the string for the file write
       */
      amem_id = H5Screate (H5S_SCALAR);
    
      /*
       * create the attribute and save the result for returning the status in hdferr
       * set hdferr here since we will be closing the a_id pointer at the end
       */

      /* HDF5 1.8 API
       * 5th arg is attribute creation property list (must be H5P_DEFAULT as of 4/19/12)
       * 6th arg is attribute access property list (must be H5P_DEFAULT as of 4/19/12)
       */
      a_id = H5Acreate (*id,name,atype,amem_id,H5P_DEFAULT,H5P_DEFAULT);
      }
  
    if (a_id < 0)
      {
      *hdferr = a_id;
      }
  
    /*
     * If we have an error don't attempt the write, but do continue on to
     * the close statments.
     */
  
    if (*hdferr == 0)
      {
      /*
       * The third arg to H5Awrite is a pointer to the data. value is that already,
       * but this does put the onus on the caller to make sure that value is really
       * pointing to the type that matches what's in vtype
       */
      H5Awrite (a_id,atype,value);
      }
  
    /*
     * clean up the attribute related pointers
     */
    if (*vtype == RHDF5_TYPE_STRING)
      {
      /*
       * Only call the close when using the string type (which is a complex structure)
       */
      H5Tclose (atype);
      }
    if (attr_exists == 0)
      {
      /*
       * if created the attribute, then close the memory space id
       */
      H5Sclose (amem_id);
      }
    H5Aclose (a_id);
    }
  
  return;
  }

/***********************************************************************
 * This routine will read in the attribute given by the object id and
 * name. The caller requests a certain type (arg: dtype) and that is checked
 * against the type found in the hdf5 file.
 */
 void rh5a_read_anyscalar (int *dsetid,char *name,void *value,int *vtype,int *hdferr)
  {
  hid_t attrid;
  hid_t atypid;
  int atype;
  int asize;
  
  *hdferr = 0;
  
  /*
   * First open the attribute
   */
  attrid = H5Aopen (*dsetid,name,H5P_DEFAULT);
  if (attrid < 0)
    {
    *hdferr = -1;
    return;
    }

  /*
   * Get the information about the data type
   */
  atypid = H5Aget_type (attrid);
  if (atypid < 0)
    {
    *hdferr = -1;
    return;
    }
  rf2c_get_mem_type (atypid,&atype,&asize);
  if (atype < 0)
    {
    *hdferr = -1;
    return;
    }
  if (atype != *vtype)
    {
    printf("ERROR: rh5a_read_anyscalar: Attribute type requested by caller");
    printf("       does not match type in file\n");
    printf("ERROR: Attribute naame: %s\n", name);
    printf("ERROR: Requested attribute type: %s\n",rhdf5_type_names[*vtype]);
    printf("ERROR: Attribute type found in HDF5 file: %s\n",rhdf5_type_names[atype]);
    *hdferr = -1;
    return;
    }

  /*
   * read the dataset, clean up and return
   */
  *hdferr = H5Aread (attrid,atypid,value);

  if (*vtype == RHDF5_TYPE_STRING)
    {
    H5Tclose (atypid);
    }
  H5Aclose (attrid);

  return;
  }

/*************** DIMENSION SCALE ROUTINES ***************/
 void rh5ds_set_scale (int *dsetid,char *dimname,int *hdferr)
  {
  *hdferr = H5DSset_scale (*dsetid,dimname);

  return;
  }

 void rh5ds_attach_scale (int *dsetid,int *dsclid,int *index,int *hdferr)
  {
  *hdferr = H5DSattach_scale (*dsetid,*dsclid,*index);

  return;
  }

/*************** INTERNAL ROUTINES ***************/

/*******************************************************************
 * This routine translates the integer coded memory type to the
 * corresponding HDF5 type.
 *
 * Note that ssize input argument is only used for string types.
 */

hid_t rf2c_set_mem_type (int dtype,int ssize)
{
hid_t mtype;

switch (dtype)
  { 
  case RHDF5_TYPE_STRING:
    /*
     * define type for the string
     * C style string, null byte terminated, size taken from input argument ssize
     */
    mtype = H5Tcopy (H5T_C_S1);
    H5Tset_size (mtype,ssize);
    H5Tset_strpad (mtype,H5T_STR_NULLTERM);
    break;

  case RHDF5_TYPE_INTEGER:
    /*
     * integer
     */
    mtype = H5T_NATIVE_INT;
    break;

  case RHDF5_TYPE_FLOAT:
    /*
     * float
     */
    mtype = H5T_NATIVE_FLOAT;
    break;

  case RHDF5_TYPE_CHAR:
    /*
     * char
     */
    mtype = H5T_NATIVE_UCHAR;
    break;

  default:
    /*
     * unrecognized type, send back -1 to inform caller
     */
    mtype = -1;
    break;
  }

return(mtype);
}

/*******************************************************************
 * This routine finds and translates the HDF5 type to
 * the integer coded memory type.
 */
 void rf2c_get_mem_type (hid_t typid,int *type,int *size)
  {
  /*
   * First check the class to see if this is a string type
   * since string type is not a "native" type. All the other
   * types (integer, float, char) are native types which can
   * checked using H5Tequal. H5Tequal returns positive value
   * for "true", zero for "false" and negative value for "error".
   */

  if ( H5Tget_class (typid) == H5T_STRING )
    {
    *type = RHDF5_TYPE_STRING;
    }
  else
    {
    /*
     * Not a string type so check the native types
     *
     * Check the character types first since they are a subset of
     * the integer types.
     */

    if (( H5Tequal (typid, H5T_NATIVE_UCHAR) > 0) || 
        ( H5Tequal (typid, H5T_NATIVE_SCHAR) > 0 ))
      {
      *type = RHDF5_TYPE_CHAR;
      }
    else if ( H5Tequal (typid,H5T_NATIVE_INT) > 0 )
      {
      *type = RHDF5_TYPE_INTEGER;
      }
    else if ( H5Tequal (typid,H5T_NATIVE_FLOAT) > 0 )
      {
      *type = RHDF5_TYPE_FLOAT;
      }
    else
      {
      /*
       * unrecognized type
       */
      *type = -1;
      }
    }

  /*
   * grab the size of this type
   */
  *size = H5Tget_size (typid);

  return;
  }
