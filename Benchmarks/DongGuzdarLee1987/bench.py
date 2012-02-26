"""
    Benchmark for collisionless drift Alfven waves.

    Ref. : F. Jenko, Computer Physics Communication 125, (2000) 196-209


    Solved is the Dispersion Relation from Eq. (26)

Take care of nomralization

alpha^2 = \frac{2}{\hat{\epsilon}}
\hat{b} = \beta_e \hat{\epsilon}
\hat{mu} = \mu_e \hat{\epsilon}



"""


from Execute_Thin import *
from sets import Set
import string 

N = 5


def linspace(start, stop, count):
	x = []
	for v in range(count) : x.append(start+v/float(count-1)*(stop-start))
        return x

shear = 0.025



if sys.argv[1] == "Exe" : 

     decomposition = "1:1:1:1:4"
     helios_options = []
     


     if sys.argv[2]   == "Fig1" : helios_paramfile_name  = "Dong_Fig1.helios"
     elif sys.argv[2] == "Fig2" : helios_paramfile_name  = "Dong_Fig2.helios"
     elif sys.argv[2] == "Fig3" : helios_paramfile_name  = "Dong_Fig3.helios"
     elif sys.argv[2] == "Fig4" : helios_paramfile_name  = "Dong_Fig4.helios"
     elif sys.argv[2] == "Fig5" : helios_paramfile_name  = "Dong_Fig5.helios"
     elif sys.argv[2] == "Fig6" : helios_paramfile_name  = "Dong_Fig6.helios"
     elif sys.argv[2] == "Fig8" : helios_paramfile_name  = "Dong_Fig8.helios"
     elif sys.argv[2] == "Fig9" : helios_paramfile_name  = "Dong_Fig9.helios"
     elif sys.argv[2] == "Fig10" : helios_paramfile_name = "Dong_Fig10.helios"
     elif sys.argv[2] == "Fig12" : helios_paramfile_name = "Dong_Fig12.helios"
     else : print "ERROR"

     if sys.argv[2] in ["Fig1", "Fig2", "Fig3", "Fig9", "Fig10"]:
     
        i = 1
	
        eta_end = 4.
        if sys.argv[2] == "Fig2" : eta_end = 13
        if sys.argv[2] == "Fig1" : eta_end = 5
        for eta in linspace(1.0,eta_end,N) : 
            option = "\"DataOutput.OutputFileName=" + sys.argv[2] + "_A_q" + string.zfill(i, 3)+\
		    	"_Eta_scan.h5; Plasma.Species1.w_T=%f;\"" % (eta)
            helios_options.append(option)

            if sys.argv[2] in ["Fig9", "Fig10"]:
          	    option = "\"DataOutput.OutputFileName=" + sys.argv[2] + "_B_q" + string.zfill(i, 3)+\
			    "_Eta_scan.h5; Plasma.Species1.w_T=%f;Plasma.Beta=0.01\"" % (eta)
          	    helios_options.append(option)

            i += 1

            
        
     if sys.argv[2] in ["Fig4", "Fig5", "Fig6"]:
        i = 1
	
        for beta in linspace(0.0, 0.02,N) : 
          	
		 
            option = "\"DataOutput.OutputFileName=" + sys.argv[2] + "_A_q" + string.zfill(i, 3)+\
				"_Beta_scan.h5; Plasma.Beta=%f;\"" % (beta)
            helios_options.append(option)
		    
            option = "\"DataOutput.OutputFileName=" + sys.argv[2] + "_B_q" + string.zfill(i, 3)+\
				"_Beta_scan.h5; Plasma.Beta=%f;Plasma.Species0.Density=1.;Grid.Ns=1;\"" % (beta)
            helios_options.append(option)
            i += 1
     
     if sys.argv[2] in ["Fig8"]:
        i = 1
	
        for beta in linspace(0.0, 0.02,N) : 
          	
		 
            option = "\"DataOutput.OutputFileName=" + sys.argv[2] + "_A_q" + string.zfill(i, 3)+\
				"_Beta_scan.h5; Plasma.Beta=%f;\"" % (beta)
            helios_options.append(option)
		    
            # No Gyro-averaging
            option = "\"DataOutput.OutputFileName=" + sys.argv[2] + "_B_q" + string.zfill(i, 3)+\
				"_Beta_scan.h5; Plasma.Beta=%f;Grid.Nm=1;\"" % (beta)
            helios_options.append(option)
            i += 1
            
            decomposition = "1:1:1:1:1"

     if sys.argv[2] in ["Fig12"]:
        i = 1
	
        for L in linspace(0.025, 0.20,N) : 
          	
		  
            option = "\"DataOutput.OutputFileName=" + sys.argv[2] + "_A_q" + string.zfill(i, 3)+\
                "_Beta_scan.h5; Plasma.Beta=0.0;Geometry.Shear=%f;\"" % L
            helios_options.append(option)
            
            option = "\"DataOutput.OutputFileName=" + sys.argv[2] + "_B_q" + string.zfill(i, 3)+\
                  "_Beta_scan.h5; Plasma.Beta=0.01;Geometry.Shear=%f;\"" % L
            #helios_options.append(option)
            i += 1
     
    
     
    
     # Finally execute all 
     #executeThin(helios_paramfile_name, helios_options, decomposition, queue='lh10247')
     executeThin(helios_paramfile_name, helios_options, decomposition, queue='vh10140')
     #executeThin(helios_paramfile_name, helios_options, decomposition, queue='eh')


if sys.argv[1] == "Ana":
    from tables import *
    from FitFunction import *
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy import interpolate

    n_start = int(sys.argv[2])
    n_stop  = int(sys.argv[3])	

    #determine Type

    result_A = [] ; result_B = []

    for fileh_name in sys.argv[5:]:
          try:
            fileh5  = openFile(fileh_name)
            phi2    = fileh5.root.Analysis.scalarValues.cols.phiEnergy[:]
            #n_start = -int(len(phi2)*9./10.)
            #n_stop  = -1
            T       = fileh5.root.Analysis.scalarValues.cols.Time[:]
            
            # divide by 2 for phi^2 and not phi, and sqrt(2) due to thermal velocity
            gamma = fitExpGrowthOptimize(T[n_start:n_stop], phi2[n_start:n_stop], 0)/2.
            
            print "File : " , fileh_name , " Fit start : ",  T[n_start] , "  stop: ", T[n_stop]

            if  sys.argv[4] in ["Fig1", "Fig2", "Fig3", "Fig9", "Fig10" ] : xv = fileh5.root.Species.cols.w_T[1]/fileh5.root.Species.cols.w_n[1]
            if  sys.argv[4] in ["Fig4" , "Fig5", "Fig6", "Fig8"] : xv = fileh5.root.Plasma._v_attrs.beta[0]
            if  sys.argv[4] in ["Fig12" ]                : xv = fileh5.root.Geometry._v_attrs.Shear[0]

            if "_A_" in fileh_name : result_A.append((xv,  gamma))
            if "_B_" in fileh_name : result_B.append((xv,  gamma))
	
            fileh5.close()
          except:
                print "Cannot open : ", fileh_name, "   (skipping)"
    
  	# Create subarrays by sorting thorugh eta
    result_A = array(result_A) ; result_B = array(result_B)



    pl1 = plot(result_A[:,0], result_A[:,1], 'r^', markersize=12.)
    
    try:
    # Interpolate result_A A
        tck = interpolate.splrep(result_A[:,0],result_A[:,1],s=0)
        xnew = linspace(result_A[:,0].min(), result_A[:,0].max(), 201)
        ynew = interpolate.splev(xnew,tck,der=0)
        plot(xnew, ynew, 'r-')
    except : print "Cannot interpolate"

    if sys.argv[4] in ["Fig4", "Fig5", "Fig6", "Fig8", "Fig9", "Fig10", "Fig12"]:

            # Interpolate result_B B
    	    
            pl2 = plot(result_B[:,0], result_B[:,1], 'bs', markersize=12.)
            try:
                tck = interpolate.splrep(result_B[:,0],result_B[:,1],s=0)
                xnew = linspace(result_B[:,0].min(), result_B[:,0].max(), 201)
                ynew = interpolate.splev(xnew,tck,der=0)
                plot(xnew, ynew, 'b-')
            except: print "Cannot Interpolate B"
     
     # Set XLABEL
    if sys.argv[4] in ["Fig1", "Fig2", "Fig3", "Fig9"]: xlabel("$\\eta_i$")
    if sys.argv[4] in ["Fig4", "Fig5", "Fig6", "Fig7", "Fig8", "Fig13", "Fig16" ]: xlabel("$\\beta_e$")
    if sys.argv[4] in ["Fig10"                       ]: xlabel("$\\eta_e$")
    if sys.argv[4] in ["Fig12"                       ]: xlabel("L")

     
    ylabel("Growthrate of$\\phi_{rms} \\gamma$") 
    
    if sys.argv[4] in ["Fig12"                       ]: 
      # try to approximate
      L = linspace(0.025,0.2,201)
      plot(L, interpolate.interp1d([0.025, 0.05, 0.2], [0.039, 0.043, 0.011], kind=2)(L), 'r-')
      plot(L, interpolate.interp1d([0.025, 0.09, 0.2], [0.001, 0.019, 0.003], kind=2)(L), 'b-')
      ylim((0.0,0.05))
    
    if sys.argv[4] in ["Fig9", "Fig10", "Fig12"]: legend((pl1, pl2), ("$\\beta = 0.00$", "$\\beta=0.01$"))
    
    savefig(sys.argv[4]+".png")










