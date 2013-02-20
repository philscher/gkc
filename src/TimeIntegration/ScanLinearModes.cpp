/*
 * =====================================================================================
 *
 *       Filename: ScanLinearModes.h
 *
 *    Description:Time Integration Interface for gkc. 
 *
 *         Author: Paul P. Hilscher (2011), 
 *
 *        License: GPLv3+
 * =====================================================================================
 */

#include "ScanLinearModes.h"

#include "Tools/Tools.h"

ScanLinearModes::ScanLinearModes(Setup *setup, Grid *grid, Parallel *_parallel, Vlasov *_vlasov, Fields *_fields, FileIO *fileIO) :
  parallel(_parallel)
{

  std::string modes = setup->get("ScanLinearModes.N", "1");
  
  minAbsError_re = setup->get("ScanLinearModes.MinAbsErrorFreq", 1.e-3);
  minAbsError_im = setup->get("ScanLinearModes.MinAbsErrorGrow", 1.e-3);
   
  InnerIterMax   = setup->get("ScanLinearModes.InnerIterMax", 100);
      
  //ArrayField0 = nct::allocate(fields->Field0)(&Field0);
  ArrayField0 = nct::allocate(nct::Range(0,Nq), grid->RzLD, grid->RkyLD, grid->RxLD)(&Field0_t);

  initData(setup, fileIO); 

}

ScanLinearModes::~ScanLinearModes()
{
  H5Gclose(scanGroupID);
  delete freqTable;
}


void ScanLinearModes::solve(Vlasov *vlasov, Fields *fields, TimeIntegration *timeIntegration, Eigenvalue *eigenvalue, Init *init, Visualization *visual)
{

  // use 3/5 point average and caluclate rms

  // set modes
  //for(int n = 0; n < 65; n++) modes.push_back(n*0.05+0.1);

  modes = Tools::logspace(-1., 1., 10);

  for(auto ky = modes.begin(); ky != modes.end(); ++ky) {

    /// Currently we cannot iterate over ky list, thus change Ly (global variable)
    Ly = 2.*M_PI / *ky; 

    Timing timing(0, 0.);
      
    // Include ZF and Nyquist in calculations
    // Note : don't screen out Nyquist and ZF ! 
    timeIntegration->setMaxLinearTimeStep(eigenvalue, vlasov, fields);
  
    //init->PerturbationPSFNoise((A6zz) vlasov->f0, (A6zz) vlasov->f);
    init->PerturbationPSFExp((A6zz) vlasov->f0, (A6zz) vlasov->f);
    
    ErrorVals w_im(5), w_re(5);

    do { // Outer iteration 

      double dt = 0.;
     
      for(int n = 0; n < InnerIterMax; n++) {
      // Iterate over couple of timesteps
        dt += timeIntegration->solveTimeStep(vlasov, fields, NULL, timing);
      }

      // Calculate change in growthrates
      #pragma omp single
      //, copyprivate(w1_im, w0_im, w0_re, w1_re)
      [&] (CComplex Field0    [Nq][NzLD][Nky][NxLD],
           CComplex Field0_tm1[Nq][NzLD][Nky][NxLD]) {

        // actually can only be done for Field0
        // We use first order upwind to calculate $\partial A / \partial t$
        const double field_0 =  __sec_reduce_add(cabs(Field0    [0:Nq][NzLlD:NzLD][0:Nky][NxLlD:NxLD]));
        const double field_1 =  __sec_reduce_add(cabs(Field0_tm1[0:Nq][NzLlD:NzLD][0:Nky][NxLlD:NxLD]));

        // Calculte real frequency through phase shift
        const double field_0_ph =  carg(parallel->reduce(__sec_reduce_add(Field0    [0:Nq][NzLlD:NzLD][0:Nky][NxLlD:NxLD]), Op::sum, DIR_X));
        const double field_1_ph =  carg(parallel->reduce(__sec_reduce_add(Field0_tm1[0:Nq][NzLlD:NzLD][0:Nky][NxLlD:NxLD]), Op::sum, DIR_X));
        
        w_im.addValue(-1./dt * ( field_1 - field_0) / field_0);
        w_re.addValue((field_1_ph - field_0_ph) / dt);
         
        // update 
        Field0_tm1[0:Nq][NzLlD:NzLD][0:Nky][NxLlD:NxLD] = Field0[0:Nq][NzLlD:NzLD][0:Nky][NxLlD:NxLD];

      } ((A4zz) fields->Field0, (A4zz) Field0_t);

      // if accuracy is reached save
    } while(   ( (w_im.getStdDev() > minAbsError_im) 
            ||   (w_re.getStdDev() > minAbsError_re) )
            && (timing <= timeIntegration->maxTiming) );
        
    std::stringstream ss;
    ss << *ky << "  Growthrate : " << w_im.getMean()   << " " << w_re.getMean() 
              << "  Error : "      << w_im.getStdDev() << " " << w_re.getStdDev() << std::endl;
      
    parallel->print(ss.str());
        
    ///////////////   Save results /////////////////////
    Frequency freqValue;

    freqValue.ky        = *ky;
    freqValue.frequency = Complex(w_re.getMean()  , w_im.getMean()     );
    freqValue.error     = Complex(w_re.getStdDev(), w_im.getStdDev());

    freqTable->append(&freqValue);
      
    visual->writeData(Timing(0, *ky), 0., true);

  }
}

void ScanLinearModes::initData(Setup *setup, FileIO *fileIO) 
{
  scanGroupID = fileIO->newGroup("LinearScan");

  // ********************* setup Table for EigenValues *****************
  Frequency freq_table;
         
  size_t freq_offsets[] = { HOFFSET(Frequency, ky), HOFFSET(Frequency, frequency), HOFFSET(Frequency, error) };
  size_t freq_sizes  [] = { sizeof(freq_table.ky), sizeof(freq_table.frequency) ,  sizeof(freq_table.error) };
  hid_t  freq_types  [] = { H5T_NATIVE_DOUBLE, fileIO->complex_tid , fileIO->complex_tid };
  const char * freq_names  [] = {"Mode", "Frequency", "Error"};

  freqTable = new TableAttr(scanGroupID, "Frequencies", 3, freq_names, freq_offsets, freq_types, freq_sizes, &freq_table);

}

