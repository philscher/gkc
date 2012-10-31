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

#include <iostream>
#include <string>
#include <sstream>
#include <stdlib.h>

inline int check( int status, std::string file, int line, std::string error_text, bool doAbort=false) {
        if(status == -1 ) {
            
       // check rank from MPI!
            std::stringstream ss;
            ss << std::endl; 
            ss << "\033[0;m";
            ss << "\033[1;m"  <<"!.... " << file << "(" << line << ") " << error_text; 
            ss << "\033[0m" << std::endl << std::flush;

            std::cout << ss.str();
              
            // exit through abort so we can get stack trace
       if(doAbort==true) abort();
            abort();
            exit(0);
        }
        return status;
    }

#define DMESG(mesg)  std::string(__FILE__), __LINE__, std::string(mesg)



class FileAttr
{
 
  // Copy Arrays
  void copy(hsize_t dim_from[], hsize_t dim_to[]) 
  {
    for(int n=0; n < n_dim; n++) dim_to[n] = dim_from[n];
  };

  int n_dim;                          ///< number of dimensions

  hid_t dataset_hdf,                  ///< HDF-5 Dataset      id
        property_hdf,                 ///< HDF-5 Property     id
        memory_space_hdf,             ///< HDF-5 Memory space id 
        typeId_hdf;                   ///< HDF-5 Datatype     id

  std::string name;                   ///< Name of array

  hsize_t  dim[7],                    ///< Dimension of Array
           maxdim[7],                 ///< Maximum dimension of array 
           chunkdim[7],               ///< Chunk size 
           offset[7],                 ///< (Boundary) Offset 
           stride1[7];                ///< Set unit stride

  bool do_write;                      ///< Does process write to data

  public: 

   FileAttr(std::string _name, hid_t group, int n_dim, hsize_t _dim[], hsize_t mdim[], 
            hsize_t cdim[], hsize_t moffset[], hsize_t chunkBdim[], hsize_t _offset[],   
            bool _write, hid_t _typeId_hdf = H5T_NATIVE_DOUBLE, bool createFile=true) 
       : n_dim(n_dim), do_write(_write), typeId_hdf(_typeId_hdf), name(_name)
    
   {

       copy(_dim, dim); copy(mdim, maxdim); copy(cdim, chunkdim); copy(_offset, offset);
       for(int n =0 ; n < 7; n++) stride1[n] = 1;

       memory_space_hdf = check( H5Screate_simple(n_dim, chunkBdim, NULL) , DMESG("H5Screate_simple"));
       check(H5Sselect_hyperslab(memory_space_hdf, H5S_SELECT_SET, moffset, stride1, chunkdim, NULL), DMESG("H5Sselect_hyperslab"));
       
       if(do_write == false) H5Sselect_none(memory_space_hdf);
    
       property_hdf = H5P_DEFAULT; 
#ifdef USE_MPI
       property_hdf = H5Pcreate(H5P_DATASET_XFER);
       H5Pset_dxpl_mpio(property_hdf, H5FD_MPIO_COLLECTIVE);
#endif
       
       // Create new data space 
       if(createFile == true) {

          hid_t dspace = check ( H5Screate_simple(n_dim, dim, maxdim) , DMESG(name + " : H5Screate_simple"));
          hid_t dcpl   = H5Pcreate (H5P_DATASET_CREATE);
          H5Pset_chunk(dcpl, n_dim, chunkdim);
          dataset_hdf = check(H5Dcreate2(group, _name.c_str(), typeId_hdf, dspace, H5P_DEFAULT, dcpl, H5P_DEFAULT), DMESG(name + "H5Dcreate2"));
          H5Pclose(dcpl); H5Sclose(dspace);
     
       } else {   // open data space

          dataset_hdf = check(H5Dopen(group, name.c_str(), H5P_DEFAULT), DMESG("H5Dopen"));
       }
       

   }; 

    /**
    *    @brief Write Data to HDF-5 file
    *
    *    Writed data, by first extending last dimensions.
    *
    **/
    template<typename T> void write(T data, int increase=+1) 
    {
          // extend data-dimensions
          dim[n_dim-1] += increase;
          offset[n_dim-1]+= increase;
          check( H5Dextend (dataset_hdf, dim) , DMESG("Extending PSF")); 

          // open corresponding data-space (why cannot be open only once ?)
          hid_t dspace = H5Dget_space (dataset_hdf);

          // select hyperslab (ignore if not part of reading family) & write 
          check ( H5Sselect_hyperslab (dspace, H5S_SELECT_SET, offset, stride1, chunkdim, NULL), DMESG(name + " : Selecting Hyperslab"));
          if(do_write == false) H5Sselect_none(dspace); 
          check(  H5Dwrite(dataset_hdf, typeId_hdf, memory_space_hdf, dspace, property_hdf, data)  , DMESG(name + " : H5DWrite Dataset"));
          H5Sclose(dspace);
        
    };
   
    /**
    *    @brief Reads data from HDF-5 file
    *
    *    Offset=-1 refers to the last written data file 
    *    in last dimension.
    *
    *
    **/
    template<typename T> void read(T data, int time_offset=-1) 
    {
          hid_t dspace = H5Dget_space (dataset_hdf);
        
          // get dimensions of stored array
          hsize_t off[7];
          H5Sget_simple_extent_dims (dspace, off, NULL);

          // shift last dimensions
          hsize_t off_2[7] = { 0, 0, 0, 0, 0, 0, off[6]+time_offset };
          
          // select hyperslab (ignore if not part of reading family) & read 
          check(H5Sselect_hyperslab (dspace, H5S_SELECT_SET, off_2, stride1, chunkdim, NULL), DMESG("H5Sselect_hyperslab"));
          if(do_write == false) H5Sselect_none(dspace); 
          check(H5Dread(dataset_hdf, typeId_hdf, memory_space_hdf, dspace, property_hdf, data), DMESG("H5Dread"));

          H5Sclose(dspace);
    }

   ~FileAttr() 
   {
      // close dataset, -space and -property
       check( H5Dclose(dataset_hdf)     ,  DMESG("H5Dclose"));
       check( H5Pclose(property_hdf)    ,  DMESG("H5Pclose"));
       check( H5Sclose(memory_space_hdf),  DMESG("H5Pclose"));
   };

};







#endif // FILE_ATTR_H_


