/*
 * =====================================================================================
 *
 *       Filename: TableAttr.h
 *
 *    Description: Wrapper for simplifying usage of HDF-5 tables.
 *
 *         Author: Paul P. Hilscher (2010-), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */

#ifndef TABLE_ATTR_H_
#define TABLE_ATTR_H_

#include "hdf5.h"
#include "hdf5_hl.h"

/**
*   @brief Wrapper for HDF-5 tables
*
*   This is a simple class to facilitate handling of HDF-5 tables.
*   It basically stores the array structures as they are needed
*   in both, H5TBmake_table  and H5TBappend_records.
*
*
*   @todo How to use ... docu
*
**/
class TableAttr
{
  int numCol; 
  size_t    offsets[32];
  size_t    sizes[32];
  hid_t     types[32];
  hid_t     nodeID;

  // copy values 
  template<typename T> void copy(T inValues[], T copiedValues[]) {for(int n=0;n<numCol; n++) copiedValues[n] = inValues[n];};
  std::string name;

 public: 

  /**
  *   @brief please document me ...
  *
  **/
  template<typename T>
  TableAttr(hid_t _nodeID, std::string  _name, int _numCol, const char *field_names[], 
            size_t _offsets[], hid_t _types[], size_t _sizes[], T *table) : numCol(_numCol), nodeID(_nodeID) 
  {
    
    check((numCol > 32) ? -1 : 1, DMESG("TableAttr : number of fields is limited to 32. Increase TABLE_LIM_MAX"));
         
    name = _name;

    // copy data which is later used by append
    copy(_offsets,offsets);
    copy(_types,types);
    copy(_sizes, sizes);;
           
    check(H5TBmake_table(name.c_str(), nodeID, name.c_str(), (hsize_t) numCol, (hsize_t) 0, sizeof(T), (const char**) field_names,
                         offsets, types, 100, NULL, 0, table ), DMESG("H5Tmake_table : scalarValue"));
  }

  /**
  *   @brief please document me ...
  *
  **/
  template<class T> void append(T *table, int n=1) 
  {
    check(H5TBappend_records (nodeID, name.c_str(), n, sizeof(T), offsets, sizes, table), DMESG("Append Table"));
  };


  /**
  *   @brief please document me ...
  *
  **/
 ~TableAttr() 
 {
//       H5Gclose(node); 
 };
};

#endif // TABLE_ATTR_H
