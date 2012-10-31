#include "MHD.h"


#include <cstdarg>


extern "C" double conj (const cmplxd z);
extern "C" double creal(const cmplxd z);
extern "C" double cimag(const cmplxd z);

template<class T> inline T max (T a, T b) { return a > b ? a : b; };
template<class T> inline T pow2(T a)      { return a*a;          };


MHD::MHD(std::string setup_filename, std::string setup_Xoptions)  
{
    int num_threads = 0;
    #pragma omp parallel shared(num_threads)
    {
      num_threads =  omp_get_num_threads();
    }

	 std::cout << "Running with - " << num_threads << " - OpenMP Threads" << std::endl;
    if(Nky % omp_get_num_threads() != 0) std::cout << "WARNING : Number of Threads not multiple of number of Ky modes";
    
    
    input = new Input(setup_filename, setup_Xoptions);

    
    // some standard parameters
    double Lx   = input->get("Lx" , 5.0 ); 
    double Ly   = input->get("Ly" , 2. * M_PI );     
      
    Viscosity   = input->get("Viscosity"  , 1.e-55);   
    Resistivity = input->get("Resistivity", 1.e-4 );   
    doNonLinear = input->get("doNonLinear", 0 );   
      
    if(Nx <= 16) check(-1, DMESG("Nx < 16, choose Nx >= 16"));
    
    initVariables(Lx, Ly);

    
    data = new DataOutput(input, Nx, Nky, Lx, Ly, X, Psi0, Viscosity, Resistivity);
  
}

void MHD::startMainLoop() 
{

     double Time = 0., TimeMax=0.; 
     int    Step = 0 , StepMax=0 ;
    
     TimeMax = input->get("TimeMax", 5000);
     StepMax = input->get("StepMax", -1  );
     
     double dt    = input->get("dt", max(TimeMax/StepMax, 0.01));

     
     #pragma omp parallel, shared(Step, Time)
     {

       for(; ((Step <= StepMax) || (StepMax == -1)) && (Time <= TimeMax);) 
       {
           
            ///////////// First-half step ////////////////
           #pragma omp for 
           for(int y_k = 0; y_k < Nky; y_k++) { 
            
            swap_Old_New(y_k);
	         calculateRHSLinear(y_k);
           };
      
           if(doNonLinear == true) calculateRHSNonLinear();
         
           #pragma omp for
           for(int y_k = 0; y_k < Nky; y_k++) { 
	      
           calculateVor(y_k, 0.5*dt);
	        setBoundary(y_k);
	        calculateDerivatives(y_k);
          
            ////////////// Second Half-Step ///////////////
	        calculateRHSLinear(y_k);
         
         };

         if(doNonLinear == true) calculateRHSNonLinear();

         #pragma omp for
         for(int y_k = 0; y_k < Nky; y_k++) { 

	      calculateVor(y_k, dt);
	      setBoundary(y_k);
	      calculateDerivatives(y_k);
   
         }
         
         // Data Output and Analysize is only handles by single thread
         #pragma omp single
         {
            Time += dt; Step++;

            data->output(Step, Time, dt, Nky, Nx, Phi, Psi);
                
            if(Step % 100 == 0) 
               std::cout << "\r" << " Step : " << Step << "/" << StepMax << " Time : " << Time  << "/" << TimeMax;// << std::flush;
          }
        }

     } // parallel section
     std::cout << " Simulation Finished " << std::endl;
}


void MHD::setBoundary(const int y_k) 
{

      // Set Zero-Boundary 
      Phi[y_k][0:4] = 0.;   Phi[y_k][Nx-4:4] = 0.;
      Vor[y_k][0:4] = 0.;   Vor[y_k][Nx-4:4] = 0.;
      Cur[y_k][0:4] = 0.;   Cur[y_k][Nx-4:4] = 0.;
      Psi[y_k][0:4] = 0.;   Psi[y_k][Nx-4:4] = 0.;
      
     dPhi[y_k][0:4] = 0.;  dPhi[y_k][Nx-4:4] = 0.;
     dVor[y_k][0:4] = 0.;  dVor[y_k][Nx-4:4] = 0.;
     dCur[y_k][0:4] = 0.;  dCur[y_k][Nx-4:4] = 0.;
     dPsi[y_k][0:4] = 0.;  dPsi[y_k][Nx-4:4] = 0.;
	   
    ddPhi[y_k][0:4] = 0.; ddPhi[y_k][Nx-4:4] = 0.;
      

}

MHD::~MHD() 
{
    delete data;
};

void MHD::initVariables(double Lx, double Ly) 
{
	  
  const double dx = 2. * Lx / ((double) (Nx-1)); 
    
  _kw_2_dx = 1./(2.*dx); _kw_dx2 = 1./(dx*dx);
      
  // Parse equilibrium (psi0)
  FunctionParser Psi0_parser = input->getFParser();
  check(((Psi0_parser.Parse(input->get("Psi0", "1./cosh(x)^2"), "x") == -1) ? 1 : -1), DMESG("Parsing error"));
  // Parse equilibrium (Perturbation)
  FunctionParser Pert_parser = input->getFParser();
  check(((Pert_parser.Parse(input->get("Pert", "1.e-8*(exp(-(x-1.2)^2)+exp(-(x+1.2)^2))"), "x") == -1) ? 1 : -1), DMESG("Parsing error"));

  
  // Set X-Axis & Equilibrium      
  for(int x=0; x < Nx; x++) {

    X[x] = -Lx + x * dx;

    double x_0  = X[x]       ; double psi0_0  = Psi0_parser.Eval(&x_0 );
    double x_p1 = X[x]+dx    ; double psi0_p1 = Psi0_parser.Eval(&x_p1);
    double x_m1 = X[x]-dx    ; double psi0_m1 = Psi0_parser.Eval(&x_m1);
    double x_p2 = X[x]+dx+dx ; double psi0_p2 = Psi0_parser.Eval(&x_p2);
    double x_m2 = X[x]-dx-dx ; double psi0_m2 = Psi0_parser.Eval(&x_m2);
         
    
     Psi0[x] = psi0_0;
    dPsi0[x] = (psi0_p1 - psi0_m1) * _kw_2_dx;
     Cur0[x] = (psi0_p1 - 2. * psi0_0 +  psi0_m1) * _kw_dx2;
    
    dCur0[x] = (psi0_p2 - 2. * psi0_p1 + 2. * psi0_m1 - psi0_m2) * _kw_2_dx * _kw_dx2; 
			
    
    // Set inital perturbation
    
    Psi[:][x]= Pert_parser.Eval(&x_0);
    Phi[:][x]= Pert_parser.Eval(&x_0);
    Vor[:][x]= 0.;
    
  }
      
  for(int y_k = 0; y_k < Nky; y_k++) { 
  
    // Set ky-Axis
    ky[y_k]   = 2.*M_PI*((cmplxd) 0.+ 1.I) * ((double) (y_k))/Ly;
    setBoundary(y_k);
    calculateDerivatives(y_k);
  }
}


void MHD::calculateDerivatives(const int y_k) 
{
   
	     // central differences (first derivative, 2nd order)
        // central difference (2nd  derivatrive, 2nd order)
        dPhi [y_k][1:Nx-2] = (Phi[y_k][2:Nx-2] -     Phi[y_k][0:Nx-2]) * _kw_2_dx;
	     ddPhi[y_k][1:Nx-2] = (Phi[y_k][2:Nx-2] - 2.* Phi[y_k][1:Nx-2] + Phi[y_k][0:Nx-2]) * _kw_dx2;
	     
	     dPsi [y_k][1:Nx-2] = (Psi[y_k][2:Nx-2] -     Psi[y_k][0:Nx-2]) * _kw_2_dx;
	     ddPsi[y_k][1:Nx-2] = (Psi[y_k][2:Nx-2] - 2.* Psi[y_k][1:Nx-2] + Psi[y_k][0:Nx-2]) * _kw_dx2;
        
        dVor [y_k][1:Nx-2] = (Vor[y_k][2:Nx-2] - Vor[y_k][0:Nx-2]) * _kw_2_dx;
	     dCur [y_k][1:Nx-2] = (Cur[y_k][2:Nx-2] - Cur[y_k][0:Nx-2]) * _kw_2_dx;
		 
}


void MHD::calculateRHSLinear(const int y_k) 
{
     // what is caluclated here ?!
     VorRHS[y_k][:] = ky[y_k] * (dPsi0[:]*Cur[y_k][:] - dCur0[:]*Psi[y_k][:]);  
     PsiRHS[y_k][:] = ky[y_k] * (dPsi0[:]*Phi[y_k][:]);
};     

void MHD::calculateRHSNonLinear()
{
     // Take care of complex-conjugates
     // Zonal flow (y_k=0) should cancel cancel out, however, we include
     // it for parallel efficiency and confirmation. 
     #pragma omp for 
     for(int y_k = 0; y_k < Nky; y_k++) { 
     
     for(int y_k_ = (-Nky+1); y_k_ < Nky; y_k_++) {
        
       // Triad condition (Mode matching condition)
       const int y_k__ = y_k - y_k_;       
       // Ignore modes which are out of range
       if(abs(y_k__) >= Nky) continue;
       
       
       // if y_k_ < 0 use complex conjugate value for negative modes
       if     ((y_k_ >= 0) && (y_k__ >= 0)) 
       {
       
         VorRHS[y_k][:]  +=  - ky[y_k_]*Vor[y_k_][:]*dPhi[y_k__][:] + ky[y_k__]*Phi[y_k__][:]*dVor[y_k_][:] 
                             + ky[y_k_]*Cur[y_k_][:]*dPsi[y_k__][:] - ky[y_k__]*Psi[y_k__][:]*dCur[y_k_][:];
         PsiRHS[y_k][:]  +=    ky[y_k_]*Phi[y_k_][:]*dPsi[y_k__][:] - ky[y_k__]*Psi[y_k__][:]*dPhi[y_k_][:];
       } 
       
       else if((y_k_ >= 0) && (y_k__ <  0)) 
       {
         		     VorRHS[y_k][:]  +=  - ky[y_k_]*Vor[y_k_][:]*conj(dPhi[-y_k__][:]) - ky[-y_k__]*conj(Phi[-y_k__][:])*dVor[y_k_][:] 
                                        + ky[y_k_]*Cur[y_k_][:]*conj(dPsi[-y_k__][:]) + ky[-y_k__]*conj(Psi[-y_k__][:])*dCur[y_k_][:];
                    PsiRHS[y_k][:]  +=    ky[y_k_]*Phi[y_k_][:]*conj(dPsi[-y_k__][:]) + ky[-y_k__]*conj(Psi[-y_k__][:])*dPhi[y_k_][:];
       }
       
       else if((y_k_ < 0) && (y_k__ >= 0)) 
       {
         		     VorRHS[y_k][:]  +=  + ky[-y_k_]*conj(Vor[-y_k_][:])*dPhi[y_k__][:] + ky[y_k__]*Phi[y_k__][:]*conj(dVor[-y_k_][:]) 
                                        - ky[-y_k_]*conj(Cur[-y_k_][:])*dPsi[y_k__][:] - ky[y_k__]*Psi[y_k__][:]*conj(dCur[-y_k_][:]);
                    PsiRHS[y_k][:]  +=  - ky[-y_k_]*conj(Phi[-y_k_][:])*dPsi[y_k__][:] - ky[y_k__]*Psi[y_k__][:]*conj(dPhi[-y_k_][:]);
       } 
       else // ((y_k_ <  0) && (y_k__ < 0)) 
       {
         		     VorRHS[y_k][:]  +=  + ky[-y_k_]*conj(Vor[-y_k_][:])*conj(dPhi[-y_k__][:]) - ky[-y_k__]*conj(Phi[-y_k__][:])*conj(dVor[-y_k_][:]) 
                                        - ky[-y_k_]*conj(Cur[-y_k_][:])*conj(dPsi[-y_k__][:]) + ky[-y_k__]*conj(Psi[-y_k__][:])*conj(dCur[-y_k_][:]);
                    PsiRHS[y_k][:]  +=  - ky[-y_k_]*conj(Phi[-y_k_][:])*conj(dPsi[-y_k__][:]) + ky[-y_k__]*conj(Psi[-y_k__][:])*conj(dPhi[-y_k_][:]);
       }

    } } 
}


void MHD::calculateVor(const int y_k, const double dt) 
{
 
   const double dtV=0.5*dt*Viscosity;
   const double dtR=0.5*dt*Resistivity;

   // crashes if static is used
   __declspec(align(64)) cmplxd Sub[Nx], Diag[Nx], Sup[Nx], vB[Nx];
  
   ////////////////// Setup Matrix and Solve for Vorticity   //////////////////

   // 2nd order central differences
   Sub [:] =     - dtV *        _kw_dx2                 ;
   Diag[:] =  1. - dtV * (-2. * _kw_dx2 + pow2(ky[y_k]));
   Sup [:] =     - dtV *        _kw_dx2                 ;

   // Calc ? 
   vB  [1:Nx-2] =  dt*VorRHS[y_k][1:Nx-2]        - Sub[1:Nx-2]  * VorOLD[y_k][0:Nx-2] 
                                           + (2.0-Diag[1:Nx-2]) * VorOLD[y_k][1:Nx-2] 
                                                 - Sup[1:Nx-2]  * VorOLD[y_k][2:Nx-2];
  
   // Set Zero boundary conditions 
   vB[0   :4] = 0.;
   vB[Nx-4:4] = 0.;
           
   // solution is directly stored in Vor
   solveTriDiagonalMatrix(Sub, Diag, Sup, vB, Vor, y_k);
     

   ////////////////// Setup Matrix and Solve for Flow (phi) ////////////////////

   // 2nd order central differences
   Sub [:] =     _kw_dx2;
   Diag[:] = -2.*_kw_dx2 + pow2(ky[y_k]);
   Sup [:] =     _kw_dx2;

   // Move Vortized
   vB[:] = Vor[y_k][:];
                
   // Set Zero boundary conditions 
   vB[0   :4] = 0.;
   vB[Nx-4:4] = 0.;
            
   // we store solution vector directly in Flow
   solveTriDiagonalMatrix(Sub, Diag, Sup, vB, Phi, y_k);
  
   ////////////////// Setup Matrix and Solve for Flux (psi) ////////////////////
        
   // 2nd order central differences
   Sub [:] =    - dtR        * _kw_dx2;
   Diag[:] = 1. - dtR * (-2. * _kw_dx2 + pow2(ky[y_k]));
   Sup [:] =    - dtR        * _kw_dx2;

   vB[1:Nx-2] = dt*PsiRHS[y_k][1:Nx-2]    - Sub[1:Nx-2]  * PsiOLD[y_k][0:Nx-2]
                                      +(2.-Diag[1:Nx-2]) * PsiOLD[y_k][1:Nx-2] 
                                           -Sup[1:Nx-2]  * PsiOLD[y_k][2:Nx-2];
		
   // we have zero boundary conditions 
   vB[0   :4] = 0.;
   vB[Nx-4:4] = 0.;
      
   // we store solution vector directly in Psi
   solveTriDiagonalMatrix(Sub, Diag, Sup, vB, Psi, y_k);

   //////////////////// Set new current ////////////////////////

   Cur[y_k][1:Nx-2]=        (Psi[y_k][0:Nx-2] 
                      - 2. * Psi[y_k][1:Nx-2] 
                           + Psi[y_k][2:Nx-2]) * _kw_dx2 + pow2(ky[y_k])*Psi[y_k][1:Nx-2]; 

}

void MHD::swap_Old_New(const int y_k) 
{
	   VorOLD[y_k][:] = Vor[y_k][:];
	   PsiOLD[y_k][:] = Psi[y_k][:];
}

void MHD::solveTriDiagonalMatrix(cmplxd Sub[Nx], cmplxd Diag[Nx], cmplxd Sup[Nx], cmplxd b[Nx], cmplxd X[Nky][Nx], const int y_k)
{

  // Note that 0, 1, 2, 3, [ 4,  ...,  Nx-4] , Nx-3,  Nx-2, Nx-1 is boundary domain
  for (int x = 5; x < Nx-4; x++)
  {
                const cmplxd m = Sub[x]/Diag[x-1];
                Diag[x]       += - m * Sup[x-1];
                   b[x]       += - m *   b[x-1];
  }

  // back substitution
  X[y_k][Nx-5] = b[Nx-5]/Diag[Nx-5];
  for (int x = Nx-6; x >= 4; x--) X[y_k][x]=(b[x]-Sup[x]*X[y_k][x+1])/Diag[x];
}

/* 
void MHD::solveTriDiagonalMatrix(cmplxd Sub[Nx], cmplxd Diag[Nx], cmplxd Sup[Nx], cmplxd b[Nx], cmplxd X[Nky][Nx], const int y_k)
{
        // Note that 0, 1, [ 2, 3, 4,  ...,  Nx-4 , Nx-3, ] Nx-2, Nx-1 is boundary domain
        for (int x = 3; x < Nx-2; x++)
        {
                const cmplxd m = Sub[x]/Diag[x-1];
                Diag[x]       += - m*Sup[x-1];
                   b[x]       += - m*  b[x-1];
        }

        // back substitution
        X[y_k][Nx-3] = b[Nx-3]/Diag[Nx-3];
        for (int x = Nx-4; x >= 2; x--) X[y_k][x]=(b[x]-Sup[x]*X[y_k][x+1])/Diag[x];
}


 * */
