/*
 * =====================================================================================
 *
 *       Filename: Event.cpp
 *
 *    Description: Initiates events started by triggers
 *
 *         Author: Paul P. Hilscher (2010-), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */

#include "Event.h"

 Event::Event(Setup *setup, Grid *grid, Parallel *parallel, FileIO *fileIO, Geometry<HELIOS_GEOMETRY> *geo) 
{
    useEvent = setup->get("Event.use", 0);

};


void Event::checkEvent(Timing timing, Vlasov *vlasov, Fields *fields) {
    if(useEvent == false) return;

    // Example turn of collisionality
    /*
    if(timing.time > 50) {
        vlasov->collisionBeta = 0.;
        useEvent = false;
        std::cout << "Trigger Event";
    }
    */


};
