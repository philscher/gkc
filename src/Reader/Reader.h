/*
 * =====================================================================================
 *
 *       Filename:  Reader.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  07/15/2011 04:05:59 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */

#ifndef READER_H
#define READER_H


#include "Global.h"
#include "Setup.h"

class Reader {
public:
    Reader(Setup *setup) {};
    virtual ~Reader() {};

    virtual double getValue(const double x, const double y, const double z) = 0;
};


#endif // READER_H
