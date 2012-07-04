/*
 * =====================================================================================
 *
 *       Filename: Plasma.cpp
 *
 *    Description: Properties of plasma
 *
 *         Author: Paul P. Hilscher (2009-), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */



#include "Plasma.h"


Plasma::Plasma(Setup *setup, FileIO *fileIO, Geometry<HELIOS_GEOMETRY> *geo, const int _nfields) : nfields(_nfields) {
      
      species.resize(Range(0, SPECIES_MAX));
      species(0).n0 = 0.;

      debye2 = setup->get("Plasma.Debye2", 0. ); 
      B0     = setup->get("Plasma.B0"    , 1. );
       
      beta   =  setup->get("Plasma.Beta"    , 0. );
      w_p    =  setup->get("Plasma.w_p"    , 0. );
      global  = setup->get("Plasma.Global",   0 );
      
      cs      = setup->get("Plasma.cs",   1. );
      
      nfields  = ((beta > 0.e0)  ? 2 : 1);
      nfields  = (setup->get("Plasma.Bp", 0 ) == 1) ? 3 : nfields;
      
      // adiabatic species
      std::string species_name = setup->get("Plasma.Species0.Name"  , "Unnamed") + " (adiab.)";
      snprintf(species(0).name, sizeof(species_name.c_str()), "%s", species_name.c_str());
      
      species(0).n0   = setup->get("Plasma.Species0.Density" , 0. );
      species(0).T0   = setup->get("Plasma.Species0.Temperature" , 1. );
      species(0).q    = setup->get("Plasma.Species0.Charge" , 1. );
      species(0).m = 0.;
      // this is dummy for flux average
      species(0).doGyro  = setup->get("Plasma.Species0.FluxAverage", 0 );
      species(0).w_n     = setup->get("Plasma.Species0.Phase"      , 0.0 );
      species(0).w_T     = 0.;
      
      // Parse Kinetic species
      for(int s = 1; s <= SPECIES_MAX; s++) { 

        std::string key          = "Plasma.Species" + Num2String(s); 
        std::string species_name = setup->get(key + ".Name"  , "Unnamed");
        snprintf(species(s).name, sizeof(species_name.c_str()), "%s", species_name.c_str());
        species(s).m   = setup->get(key + ".Mass"       , 1. );
        species(s).n0  = setup->get(key + ".Density"    , 0. );
        species(s).T0  = setup->get(key + ".Temperature", 1. );
        species(s).q   = setup->get(key + ".Charge"     , 1. );
        species(s).gyroModel  = setup->get(key + ".gyroModel", (Nm > 1) ? "Gyro" : "Gyro-1" );
        species(s).f0_str = setup->get(key + ".F0"      , "n/(pi*T)^1.5*exp(-v^2/T)*exp(-m*B/T)" );
        species(s).f1_str = setup->get(key + ".F1"      , "0.");
        
        species(s).doGyro = (Nm > 1) ? 1 : 0;//(species(s).gyroModel == "Drift") ? 0 : 1;

        if(species(s).m < 1.e-10) check(-1, DMESG(std::string("Mass for species ") + std::string(species(s).name) + std::string(" choosen too low")));
  
        
        if(global) { 
            snprintf(species(s).n_name, 64, "%s", setup->get(key + ".n", "0." ).c_str());
            snprintf(species(s).T_name, 64, "%s", setup->get(key + ".T", "1." ).c_str());
        
            // Set Temperature and Density Gradient
            FunctionParser n_parser = setup->getFParser(); 
            FunctionParser T_parser = setup->getFParser();  
        
            check(((n_parser.Parse(species(s).n_name, "x") == -1) ? 1 : -1), DMESG("Parsing error of Initial condition n(x)"));
            check(((T_parser.Parse(species(s).T_name, "x") == -1) ? 1 : -1), DMESG("Parsing error of Initial condition T(x)"));

            // we not to normalize N, so that total density is equal in gyro-simulations 
            for(int x = NxLlB; x <= NxLuB; x++) { 
                species(s).n(x) = n_parser.Eval(&X(x));///sum(exp(-M));
                species(s).T(x) = T_parser.Eval(&X(x));
                } 
            } else {
                species(s).w_T = setup->get(key + ".w_T", 0.0 );
                species(s).w_n = setup->get(key + ".w_n", 0.0 );
                species(s).n   = species(s).n0;
                species(s).T   = species(s).T0;
                snprintf(species(s).n_name, 64, "%f", species(s).n0);
                snprintf(species(s).T_name, 64, "%f", species(s).n0);
            }

        species(s).update(geo, cs);

      }   

        // make some simple checks
        
        //Total charge density
        double rho0_tot = 0.;
        for(int s = 0; s <= NsGuD; s++) rho0_tot += species(s).q * species(s).n0;
        if(rho0_tot != 0.) check(setup->get("Plasma.checkTotalCharge", -1), DMESG("VIOLATING charge neutrality, check species q * n ! Exciting...")); 

      initDataOutput(fileIO);

    };

     void Plasma::printOn(ostream &output) const {
         output << 
               "Type       | " << (global ? " Global" : "Local") << std::endl
            << " Cs   : " << cs
            << "Model      |  " << ((do_gyro) ? "Gyrokinetic Model" : "Driftkinetic Model") << std::endl
            << "Species    | ";
            if(species(0).n0 != 0.) { output << "0. " << species(0).name << "  Density : " << species(0).n0 << " Charge : " << species(0).q << 
                                           " Temp : " << species(0).T0 << " FluxAvrg : " << (species(0).doGyro ? "Yes" : "No") << 
                                             " Phase : " << (species(0).w_n) <<  " (adiabatic) " <<  std::endl;  
            } else output << "0. -- No adiabatic species -- " << std::endl;
            
            cout.precision(5);
            for(int s = 1; s <= Ns; s++) {
            output   
            << "          +| " << s << ". " << species(s).name << "  Charge : " << species(s).q << "     Mass : " << species(s).m;
            if(global) {
             output << " T : " << species(s).T_name << " n0 : " << species(s).n_name  <<  std::endl;
            } else {
             output << " T0 : " << species(s).T0 << " n0 : " << species(s).n0 <<
             "  w_n : " << species(s).w_n << "  w_T : " << species(s).w_T << " Model : " << species(s).gyroModel << " doGyro : " << species(s).doGyro << std::endl;
            }
 }
            output <<  "           |  Debye Length^2 : " << debye2 << "   B0 : " << B0 << "  beta : " << beta << "  w_p : " << w_p << std::endl; 
     }


    void Plasma::initDataOutput(FileIO *fileIO) {
          hid_t plasmaGroup = check(H5Gcreate(fileIO->getFileID(), "/Plasma",H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT), DMESG("Error creating group for Geometry : H5Gcreate"));

         check(H5LTset_attribute_double(plasmaGroup, ".", "Debye2",  &debye2, 1), DMESG("H5LTset_attribute"));
         check(H5LTset_attribute_double(plasmaGroup, ".", "beta",  &beta, 1), DMESG("H5LTset_attribute"));
         check(H5LTset_attribute_double(plasmaGroup, ".", "B0",  &B0, 1), DMESG("H5LTset_attribute"));
        // check(H5LTset_attribute_string(plasmaGroup, ".", "Physics", ((plasma->nfields > 1)   ? "Electromagnetic" : "Electrostatic")), DMESG("H5LTset_attribute"));
         check(H5LTset_attribute_string(plasmaGroup, ".", "Model", ((do_gyro) ? "Gyrokinetic Model" : "Driftkinetic Model")), DMESG("H5LTset_attribute"));

         
         
          //////////////////////// Set Table for species.
         
         size_t species_offset[]     = { HOFFSET( Species, name ), HOFFSET( Species, q ), HOFFSET( Species, m ), HOFFSET( Species, n ), HOFFSET( Species, w_T ), HOFFSET( Species, w_n ), HOFFSET( Species, collision ) };
         size_t species_sizes[]      = { sizeof(species(0).name ), sizeof(species(0).q ), sizeof(species(0).m ), sizeof(species(0).n ), sizeof(species(0).w_T ), sizeof(species(0).w_n ), sizeof(species(0).collision ) };
         hid_t species_type[]        = { fileIO->s256_tid        , H5T_NATIVE_DOUBLE    , H5T_NATIVE_DOUBLE    , H5T_NATIVE_DOUBLE    , H5T_NATIVE_DOUBLE      , H5T_NATIVE_DOUBLE      , H5T_NATIVE_DOUBLE             };
         const char *species_names[] = { "Name"                  , "Charge"             , "Mass"               , "Density"            , "w_T"                  , "w_n"                  , "Collision"                   };

         check(H5TBmake_table("SpeciesTable", fileIO->getFileID(), "Species", (hsize_t) 7, (hsize_t) 0, sizeof(Species), (const char**) species_names,
                               species_offset, species_type, 32, NULL, 0, &species(0) ), DMESG("H5Tmake_table : Species"));
        
         // create table for all included species
          for(int s=0; s <= NsGuD; s++)  H5TBappend_records (fileIO->getFileID(), "Species", 1, sizeof(Species), species_offset, species_sizes, &species(s)); 

     /////////////////////
         
         
         H5Gclose(plasmaGroup);

    };

