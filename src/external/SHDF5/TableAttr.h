/*
 * =====================================================================================
 *
 *       Filename:  TableAttr.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  12/07/2011 02:35:22 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */


#ifndef TABLE_ATTR_H_
#define TABLE_ATTR_H_

#include "hdf5.h"
#include "hdf5_hl.h"


class TableAttr
{
  int numCol; 
  size_t    offsets[32];
  size_t    sizes[32];
  hid_t     types[32];
  hid_t     nodeID;
  template<typename T> void copy(T inValues[], T copiedValues[]) {for(int n=0;n<numCol; n++) copiedValues[n] = inValues[n];};
  std::string name;

  public: 
  template<typename T>
  TableAttr(hid_t _nodeID, std::string  _name, int _numCol, const char *field_names[], size_t _offsets[], hid_t _types[], size_t _sizes[], T *table) : numCol(_numCol), nodeID(_nodeID) 
  {
            check((numCol > 32) ? -1 : 1, DMESG("TableAttr : number of fields is limited to 32. Increase TABLE_LIM_MAX"));
         
            name = _name;

            // copy data which is later used by append
            copy(_offsets,offsets);
            copy(_types,types);
            copy(_sizes, sizes);;
           
            // probably bads solution, but we should close twice (not working anyway)
            //            node = _node;
            //            H5Oincr_refcount( node );
            //node = check(H5Gcreate(_node, name.c_str(),H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT), DMESG("Error creating group file for Spectrum : H5Gcreate"));

            check(H5TBmake_table(name.c_str(), nodeID, name.c_str(), (hsize_t) numCol, (hsize_t) 0, sizeof(T), (const char**) field_names,
                                           offsets, types, 100, NULL, 0, table ), DMESG("H5Tmake_table : scalarValue"));
  }
  
  template<class T> int append(T *table, int n=1) {
        
            check(H5TBappend_records (nodeID, name.c_str(), n, sizeof(T), offsets, sizes, table), DMESG("Append Table"));
            return 0;


  };


   ~TableAttr() {
//       H5Gclose(node); 
   };

};


#endif // TABLE_ATTR_H
