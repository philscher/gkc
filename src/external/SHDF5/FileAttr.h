/*
 * =====================================================================================
 *
 *       Filename: FileAttr.h
 *
 *    Description: Wrapper for using HDF-5 with CArray
 *
 *         Author: Paul P. Hilscher (2010-), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */

#ifndef FILE_ATTR_H_
#define FILE_ATTR_H_

#include "hdf5.h"
#include "hdf5_hl.h"

/**
*  
*  @brief Class for simplifying HDF-5 dataset creation
*
*  Warning if we use std::string a.c_str(), the resulting pointer
*  is invalid once a gets out of scope !
*
**/
class FileAttr
{
 
  // Copy Arrays
  void copy(hsize_t dim_from[], hsize_t dim_to[]) 
  {
    for(int n=0; n < ndim; n++) dim_to[n] = dim_from[n];
  };

  int ndim;                          ///< number of dimensions

  hid_t dataset_hdf,                 ///< HDF-5 dataset      id
        property_hdf,                ///< HDF-5 property     id
        memory_space_hdf,            ///< HDF-5 memory space id 
        typeId_hdf,                  ///< HDF-5 data type    id
        file_hdf,
        group_hdf;                   ///< Group-Id for dataset

  std::string name      ,            ///< Name of array
              loc_string;            ///< Location string

  hsize_t  dim [7],                  ///< Dimension of Array
           mdim[7],                  ///< Maximum dimension of array 
           cdim[7],                  ///< Chunk size 
           off [7],                  ///< (Boundary) off 
           stride1[7];               ///< Set unit stride

  bool do_write;                     ///< Does process write to data

 public: 

  /**
  *
  *
  *    cDdim - Chunk domain dimension
  *    cBdim - Chunk boundary dimension
  *    cmoff - Chunk Memory off
  *    cdoff - Chunk Data off
  *    dsdim - Dataset dimensions
  *    dmdim - Dataset maximum dimension
  *   
  **/
  FileAttr(std::string _name, hid_t group, hid_t file_id, int _ndim, 
           hsize_t  _dim[] , hsize_t _mdim[], 
           hsize_t _cdim[], hsize_t  moff[], 
           hsize_t cBdim[], hsize_t  _off[],   
           bool _write, 
           hid_t _typeId_hdf = H5T_NATIVE_DOUBLE, 
           bool createFile=true) 
  : ndim(_ndim), do_write(_write), typeId_hdf(_typeId_hdf), name(_name)
  {
    copy(_dim, dim); copy(_mdim, mdim); copy(_cdim, cdim); copy(_off, off);

    file_hdf = file_id;

    for(int n =0 ; n < 7; n++) stride1[n] = 1;

    memory_space_hdf = check( H5Screate_simple(ndim, cBdim, NULL) , DMESG("HDF-5 Error"));

    check(H5Sselect_hyperslab(memory_space_hdf, H5S_SELECT_SET, moff, stride1, cdim, NULL), DMESG("HDF-5 Error"));
       
    if(do_write == false) H5Sselect_none(memory_space_hdf);

    property_hdf = H5P_DEFAULT;

#ifdef USE_MPI // if domain decomposition is enables
    property_hdf = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(property_hdf, H5FD_MPIO_COLLECTIVE);
#endif

    // Create new data space 
    if(createFile == true) {

      // note that for the initial state, data size is set to zero
      dim[ndim-1] = 0;
      off[ndim-1] = -1;
      hid_t dspace = check ( H5Screate_simple(ndim, dim, mdim) , DMESG("HDF-5 Error"));
      hid_t dcpl   = H5Pcreate (H5P_DATASET_CREATE);
      H5Pset_chunk(dcpl, ndim, cdim);

      // usually the dataset is access once, tell HDF-5 thus to not to use reduce buffers 
      hid_t dapl   = H5P_DEFAULT;

      dataset_hdf = check(H5Dcreate2(group, name.c_str(), typeId_hdf, dspace, H5P_DEFAULT, dcpl, dapl), DMESG(name + "H5Dcreate2"));
      
      H5Pclose(dcpl); 
      H5Pclose(dapl); 
      H5Sclose(dspace);
      
    } else {   // open data space
      dataset_hdf = check(H5Dopen(group, name.c_str(), H5P_DEFAULT), DMESG("H5Dopen"));
    }
  }; 

  /**
  *    @brief Write Data to HDF-5 file
  *
  *    Write data, by first extending last dimensions.
  *
  **/
  template<typename T> void write(T data, int increase=+1) 
  {
    // extend data-dimensions
    dim[ndim-1] += increase;
    off[ndim-1] += increase;
      
    check( H5Dset_extent (dataset_hdf, dim) , DMESG("Extending Dataset")); 

    hid_t dspace = H5Dget_space (dataset_hdf);

    // select hyperslab (ignore if not part of reading family) & write 
    check( H5Sselect_hyperslab (dspace, H5S_SELECT_SET, off, stride1, cdim, NULL), DMESG("HDF-5 Error"));
    if(do_write == false) H5Sselect_none(dspace); 
    check( H5Dwrite(dataset_hdf, typeId_hdf, memory_space_hdf, dspace, property_hdf, data)  , DMESG("HDF-5 Error"));
    H5Sclose(dspace);
  };
   
  /**
  *    @brief Reads data from HDF-5 file
  *
  *    off=-1 refers to the last written data file 
  *    in last dimension.
  *
  **/
  template<typename T> void read(T data, int time_off=-1) 
  {
    hid_t dspace = H5Dget_space (dataset_hdf);
        
    // get dimensions of stored array
    hsize_t off[7];
    H5Sget_simple_extent_dims (dspace, off, NULL);

    // shift last dimensions
    hsize_t off_2[7] = { 0, 0, 0, 0, 0, 0, off[6]+time_off };
    
    // select hyperslab (ignore if not part of reading family) & read 
    check(H5Sselect_hyperslab (dspace, H5S_SELECT_SET, off_2, stride1, cdim, NULL), DMESG("HDF-5 Error"));
    if(do_write == false) H5Sselect_none(dspace); 
    check(H5Dread(dataset_hdf, typeId_hdf, memory_space_hdf, dspace, property_hdf, data), DMESG("HDF-5 Error"));

    H5Sclose(dspace);
  }

 ~FileAttr() 
  {
    // close dataset, -space and -property
    check( H5Dclose(dataset_hdf)     ,  DMESG("HDF-5 Error"));
    check( H5Pclose(property_hdf)    ,  DMESG("HDF-5 Error"));
    check( H5Sclose(memory_space_hdf),  DMESG("HDF-5 Error"));
  };
};

#endif // FILE_ATTR_H_

