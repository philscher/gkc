/*
 * =====================================================================================
 *
 *       Filename:  FileAttr.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  05/08/2011 03:38:16 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */

#ifndef FILE_ATTR_H_
#define FILE_ATTR_H_

class FileAttr
{
  void copy(hsize_t dim[], hsize_t dim_2[]) {for(int n=0;n<n_dim; n++) dim[n] = dim_2[n];};
  int n_dim;
  bool do_write;
  hid_t      dset, plist, mspace, type_id;
  std::string name;
  hsize_t    dim[7], maxdim[7], chunkdim[7], offset[7]; 

  public: 

  FileAttr(std::string _name, hid_t group, int n_dim, hsize_t _dim[], hsize_t mdim[], hsize_t cdim[], hsize_t moffset[], hsize_t chunkBdim[], hsize_t _offset[],   bool _write, hid_t _type_id = H5T_NATIVE_DOUBLE) : n_dim(n_dim), do_write(_write), type_id(_type_id), name(_name)
  {
     copy(dim, _dim); copy(maxdim, mdim); copy(chunkdim, cdim);

     mspace = check( H5Screate_simple(n_dim, chunkBdim, NULL) , DMESG(name + " : Creating Memspace Failed"));
     hsize_t  stride1[7] = { 1, 1, 1, 1, 1, 1, 1};
     check(H5Sselect_hyperslab(mspace, H5S_SELECT_SET, moffset, stride1, chunkdim, NULL), DMESG(name + " : Selecting Memspace Hyperslab"));
     if(do_write == false) H5Sselect_none(mspace);
     
     plist = H5P_DEFAULT; 
#ifdef HELIOS_PARALLEL_MPI
     plist = H5Pcreate(H5P_DATASET_XFER);
     H5Pset_dxpl_mpio(plist, H5FD_MPIO_COLLECTIVE);
#endif
//   };

//   int create(hid_t group, const char *path_name, hsize_t _offset[], bool no_preallocate=true)
//   {
   //     if(no_preallocate) dim[n_dim-1] = 0;
        copy(offset, _offset);
        hid_t dspace = check ( H5Screate_simple(n_dim, dim, maxdim) , DMESG(name + " : Creating Dataspace for scalarValue X"));
        hid_t dcpl   = H5Pcreate (H5P_DATASET_CREATE);
        H5Pset_chunk(dcpl, n_dim, chunkdim);
        dset = check(H5Dcreate2(group, _name.c_str(), type_id, dspace, H5P_DEFAULT, dcpl, H5P_DEFAULT), DMESG(name + " H5Dcreate2 : x_temp_dset"));
        H5Pclose(dcpl); H5Sclose(dspace);

    }; 


   // ToDo : Can we directly pass the blitz array ? and then do A.data() ?
    template<typename T, int W> void write(Array<T,W> A, int increase=+1) {
		write(A.data());
    } 
    template<typename T> int write(T data, int increase=+1) 
    {
          dim[n_dim-1] += increase;
          offset[n_dim-1]+= increase;
          
          check( H5Dextend (dset, dim) , DMESG("Extending PSF")); 
          hid_t dspace = H5Dget_space (dset);
          hsize_t  stride1[7] = { 1, 1, 1, 1, 1, 1, 1};
          check ( H5Sselect_hyperslab (dspace, H5S_SELECT_SET, offset, stride1, chunkdim, NULL), DMESG(name + " : Selecting Hyperslab for PowerSpectrumY"));
          if(do_write == false) H5Sselect_none(dspace); 
          check(  H5Dwrite(dset, type_id, mspace, dspace, plist, data)  , DMESG(name + " : H5DWrite Dataset"));
          H5Sclose(dspace);
        
          return HELIOS_SUCCESS;
    };


   ~FileAttr() {
       check( H5Dclose(dset),   DMESG(name + " : H5Dclose"));
       check( H5Pclose(plist),   DMESG(name + " H5Pclose"));
       check( H5Sclose(mspace),   DMESG(name + " H5Pclose"));
   };

};







#endif // FILE_ATTR_H_


