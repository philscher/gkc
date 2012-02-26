/*
 * =====================================================================================
 *
 *       Filename: IfaceHelios.h
 *
 *    Description: Defined common interface for gkc classes
 *
 *         Author: Paul P. Hilscher (2011), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */

class IfaceHelios {

    protected:
      virtual void printOn(ostream &output) const = 0;
      //virtual ~IfaceHelios() { closeData(); };
      virtual ~IfaceHelios() { };
    public:
    friend ostream& operator<<(ostream& output, const IfaceHelios& ih) { ih.printOn(output); return output; };

    // Data Output Operation
};
    
class IfaceDataIO {
    virtual ~IfaceDataIO() { };
    virtual void initDataOutput(FileIO *fileIO) = 0;
    virtual void writeData(Timing timing, double dt) = 0;
    virtual void closeData() = 0;
};
