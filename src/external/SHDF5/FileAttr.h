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
        typeId_hdf,                  ///< HDF-5 Datatype     id
        file_hdf,
        group_hdf;                    ///< Group-Id for dataset

  std::string name,                   ///< Name of array
              loc_string;             ///< Location string

  hsize_t  dim[7],                    ///< Dimension of Array
           maxdim[7],                 ///< Maximum dimension of array 
           chunkdim[7],               ///< Chunk size 
           offset[7],                 ///< (Boundary) Offset 
           stride1[7];                ///< Set unit stride

  bool do_write;                      ///< Does process write to data

  public: 

   /**
   *
   *
   *    cDdim - Chunk domain dimension
   *    cBdim - Chunk boundary dimension
   *    cmoff - Chunk Memory offset
   *    cdoff - Chunk Data offset
   *    dsdim - Dataset dimensions
   *    dmdim - Dataset maximum dimension
   *   
   **/
   FileAttr(std::string _name, hid_t group, hid_t file_id, int n_dim, 
               hsize_t _dim[], hsize_t mdim[], 
               hsize_t cdim[], hsize_t moffset[], 
               hsize_t chunkBdim[], hsize_t _offset[],   
               bool _write, 
               hid_t _typeId_hdf = H5T_NATIVE_DOUBLE, 
               bool createFile=true) 

       : n_dim(n_dim), do_write(_write), typeId_hdf(_typeId_hdf), name(_name)
    
   {
       copy(_dim, dim); copy(mdim, maxdim); copy(cdim, chunkdim); copy(_offset, offset);
   
       file_hdf = file_id;

       for(int n =0 ; n < 7; n++) stride1[n] = 1;

       memory_space_hdf = check( H5Screate_simple(n_dim, chunkBdim, NULL) , DMESG("H5Screate_simple"));
       check(H5Sselect_hyperslab(memory_space_hdf, H5S_SELECT_SET, moffset, stride1, chunkdim, NULL), DMESG("H5Sselect_hyperslab"));
       
       if(do_write == false) H5Sselect_none(memory_space_hdf);

       property_hdf = H5P_DEFAULT; 
#ifdef USE_MPI
       property_hdf = H5Pcreate(H5P_DATASET_XFER);
       H5Pset_dxpl_mpio(property_hdf, H5FD_MPIO_COLLECTIVE);
#endif
      
       // set fixed buffer ?! H5Pset_buffer( property_hdf, 32 * 1024 *1024, NULL, NULL);

       // Create new data space 
       if(createFile == true) {

          hid_t dspace = check ( H5Screate_simple(n_dim, dim, maxdim) , DMESG(name + " : H5Screate_simple"));
          hid_t dcpl   = H5Pcreate (H5P_DATASET_CREATE);
          H5Pset_chunk(dcpl, n_dim, chunkdim);

          // usually the dataset is access once, tell HDF-5 thus to not to use reduce buffers 
          hid_t dapl   = H5Pcreate (H5P_DATASET_ACCESS);
          H5Pset_chunk_cache(dapl, H5D_CHUNK_CACHE_NSLOTS_DEFAULT, H5D_CHUNK_CACHE_NSLOTS_DEFAULT, 1.);

          dataset_hdf = check(H5Dcreate2(group, _name.c_str(), typeId_hdf, dspace, H5P_DEFAULT, dcpl, dapl), DMESG(name + "H5Dcreate2"));
          H5Pclose(dcpl); 
          H5Pclose(dapl); 
          H5Sclose(dspace);
     
       } else {   // open data space

          dataset_hdf = check(H5Dopen(group, name.c_str(), H5P_DEFAULT), DMESG("H5Dopen"));

       }

       
      char id_name[1024];

        //H5Iget_name(group, id_name, 255 );
        H5Iget_name(dataset_hdf, id_name, 1024 );
        loc_string = std::string(id_name);
//     if(1==1) {

        //H5O_info_t object_info;
        //H5Oget_info( dataset_hdf, &object_info );
        //H5Oget_info( group, &object_info );
        //haddr_t dst_addr = object_info.addr;
        //H5Oopen_by_addr( group, dst_addr );
  //      std::cout << "Id Name : " << std::string(id_name) << std::endl;
      //  std::string id_str = std::string(id_name) + "/" + name.c_str();
//          
//       group_hdf = group;
//       H5Iinc_ref( group );
        
        H5Dclose(dataset_hdf);
        // make copy of group
     //}

   }; 

    /**
    *    @brief Write Data to HDF-5 file
    *
    *    Writed data, by first extending last dimensions.
    *
    **/
    template<typename T> void write(T data, int increase=+1) 
    {

        // always open & close dataset, otherwise we get memory leak 
   //       dataset_hdf = check(H5Dopen(group_hdf, name.c_str(), H5P_DEFAULT), DMESG("H5Dopen"));
  ///dset = H5Dopen (file, DATASET, H5P_DEFAULT);
  //space = H5Dget_space (dset);
  //ndims = H5Sget_simple_extent_dims (space, dims, NULL);
  //
  // We close and release dataset after each write, otherwise temporary buffers may accumulate
     hid_t dataset_hdf = check(H5Dopen(file_hdf, loc_string.c_str(), H5P_DEFAULT), DMESG("H5Dopen"));

          // extend data-dimensions
          dim[n_dim-1]    += increase;
          offset[n_dim-1] += increase;
          check( H5Dset_extent (dataset_hdf, dim) , DMESG("Extending Dataset")); 


          // open corresponding data-space (why cannot be open only once ?)
          hid_t dspace = H5Dget_space (dataset_hdf);

          // select hyperslab (ignore if not part of reading family) & write 
          check ( H5Sselect_hyperslab (dspace, H5S_SELECT_SET, offset, stride1, chunkdim, NULL), DMESG(name + " : Selecting Hyperslab"));
          if(do_write == false) H5Sselect_none(dspace); 
          check(  H5Dwrite(dataset_hdf, typeId_hdf, memory_space_hdf, dspace, property_hdf, data)  , DMESG(name + " : H5DWrite Dataset"));
          H5Sclose(dspace);
   
          H5Dclose(dataset_hdf);

        
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
 //      check( H5Dclose(dataset_hdf)     ,  DMESG("H5Dclose"));
       check( H5Pclose(property_hdf)    ,  DMESG("H5Pclose"));
       check( H5Sclose(memory_space_hdf),  DMESG("H5Pclose"));
//       check( H5Gclose(group_hdf),  DMESG("H5Pclose"));
   };

};







#endif // FILE_ATTR_H_


