/*
 * =====================================================================================
 *
 *       Filename:  FileAttr.h
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


class FileAttr
{
 
  // Copy Arrays
  void copy(hsize_t dim[], hsize_t dim_2[]) 
  {

    for(int n=0; n < n_dim; n++) dim[n] = dim_2[n];

  };

  int n_dim;                               ///< number of dimensions

  hid_t      dataset_hdf,                  ///< HDF-5 Dataset identifiers
             property_hdf, 
             memory_space_hdf, 
             typeId_hdf;

  std::string name;                        ///< Name of array

  hsize_t  dim[7],                         ///< Dimension of Array
           maxdim[7],                      ///< Maximum dimnesion of array 
           chunkdim[7],                    ///< Chunk size of current process 
           offset[7];                      ///< (Boundary) Offset 
  bool do_write;                           ///< Does process write to data

  public: 

  FileAttr(std::string _name, hid_t group, int n_dim, hsize_t _dim[], hsize_t mdim[], 
           hsize_t cdim[], hsize_t moffset[], hsize_t chunkBdim[], hsize_t _offset[],   
           bool _write, hid_t _typeId_hdf = H5T_NATIVE_DOUBLE, bool createFile=true) : n_dim(n_dim), do_write(_write), typeId_hdf(_typeId_hdf), name(_name)
  {

     copy(dim, _dim); copy(maxdim, mdim); copy(chunkdim, cdim);
    

       memory_space_hdf = check( H5Screate_simple(n_dim, chunkBdim, NULL) , DMESG(name + " : Creating Mememory_space_hdf Failed"));
       hsize_t  stride1[7] = { 1, 1, 1, 1, 1, 1, 1};
       check(H5Sselect_hyperslab(memory_space_hdf, H5S_SELECT_SET, moffset, stride1, chunkdim, NULL), DMESG(name + " : Selecting Mememory_space_hdf Hyperslab"));
       
       if(do_write == false) H5Sselect_none(memory_space_hdf);
     

    
       property_hdf = H5P_DEFAULT; 
#ifdef USE_MPI
       property_hdf = H5Pcreate(H5P_DATASET_XFER);
       H5Pset_dxpl_mpio(property_hdf, H5FD_MPIO_COLLECTIVE);
#endif
       copy(offset, _offset);
       
       // Create new data space 
     if(createFile == true) {

       hid_t dspace = check ( H5Screate_simple(n_dim, dim, maxdim) , DMESG(name + " : Creating Dataspace for scalarValue X"));
       hid_t dcpl   = H5Pcreate (H5P_DATASET_CREATE);
       H5Pset_chunk(dcpl, n_dim, chunkdim);
       dataset_hdf = check(H5Dcreate2(group, _name.c_str(), typeId_hdf, dspace, H5P_DEFAULT, dcpl, H5P_DEFAULT), DMESG(name + " H5Dcreate2 : x_temp_dset"));
       H5Pclose(dcpl); H5Sclose(dspace);
     
     } else {   // open data space

        hid_t dspace = check(H5Dopen(group, name.c_str(), H5P_DEFAULT), DMESG(" H5Dopen : PSF"));
        hid_t dcpl   = H5Dget_space (dspace);
        dataset_hdf = H5Dopen(group, name.c_str(), dcpl );
        H5Pclose(dcpl); H5Sclose(dspace);
    }
       

   }; 



    template<typename T> void write(T data, int increase=+1) 
    {
          dim[n_dim-1] += increase;
          offset[n_dim-1]+= increase;
          
          check( H5Dextend (dataset_hdf, dim) , DMESG("Extending PSF")); 
          hid_t dspace = H5Dget_space (dataset_hdf);
          hsize_t  stride1[7] = { 1, 1, 1, 1, 1, 1, 1};
          check ( H5Sselect_hyperslab (dspace, H5S_SELECT_SET, offset, stride1, chunkdim, NULL), DMESG(name + " : Selecting Hyperslab"));
          if(do_write == false) H5Sselect_none(dspace); 
          check(  H5Dwrite(dataset_hdf, typeId_hdf, memory_space_hdf, dspace, property_hdf, data)  , DMESG(name + " : H5DWrite Dataset"));
          H5Sclose(dspace);
        
    };
   
    /**
    *    @brief reads data
    *
    *    Offset=-1 refers to the last written data file 
    *    in last dimension.
    *
    *
    **/
    template<typename T> void read(T data, int offset=-1) 
    {
          
          hid_t dspace = H5Dget_space (dataset_hdf);
          hsize_t  stride1[7] = { 1, 1, 1, 1, 1, 1, 1};

          check (H5Sselect_hyperslab (dspace, H5S_SELECT_SET, offset, stride1, chunkdim, NULL), DMESG(name + " : Selecting Hyperslab for PowerSpectrumY"));
          check(H5Dread(dataset_hdf, typeId_hdf, memory_space_hdf, dspace, property_hdf, data), DMESG("H5Dread : PSF"));

          H5Sclose(dspace);
    }

   ~FileAttr() 
   {
       check( H5Dclose(dataset_hdf),    DMESG(name + " : H5Dclose"));
       check( H5Pclose(property_hdf),   DMESG(name + " H5Pclose"));
       check( H5Sclose(memory_space_hdf),  DMESG(name + " H5Pclose"));
   };

};







#endif // FILE_ATTR_H_


