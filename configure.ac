AC_INIT([GKC],[0.20],[mail@gkc.special-combo.net])
AM_INIT_AUTOMAKE([-Wall -Werror foreign]) 

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_SRCDIR([src])

AC_CANONICAL_HOST

AC_PROG_CXX($CXX g++ sunCC mpiicc icc)


####################################
# Set code options

# profiling
AC_ARG_ENABLE([profile], AC_HELP_STRING([--enable-profil], [gkc is compiled with profiling support]), [profile=yes], [profile=no])
if test x$profile = xyes ; then
   AC_DEFINE([PROFILE], [], [enables profiling])
fi

# debugging
AC_ARG_ENABLE([debug], AC_HELP_STRING([--enable-debug], [enables debug output messages and information for debugging]), [debug=yes], [debug=no])
if test x$debug = xyes ; then
   AC_DEFINE([DEBUG], [], [enables debug output messages])
fi
 
case $host in 
       i?86-*-darwin*)  AC_DEFINE([OS_DARWIN],[],[A]) ;;
       i?86-*-Linux*)   AC_DEFINE([OS_LINUX],[],[B])   ;;
esac


#################### Check for various libraries ###############################

AC_ARG_ENABLE([petsc], AC_HELP_STRING([--enable-petsc], [compile with PETSc support]), [petsc=yes], [petsc=no])
if test x$petsc = xyes ; then
   AC_DEFINE([GKC_PETSC], [], [enables PETSc support])
fi

AC_ARG_ENABLE([slepc], AC_HELP_STRING([--enable-slepc], [compile with SLEPc support]), [slepc=yes], [slepc=no])
if test x$slepc = xyes ; then
   AC_DEFINE([GKC_SLEPC], [], [enables SLEPc support])
fi

AC_ARG_ENABLE([mkl], AC_HELP_STRING([--enable-mkl]  , [compile with Intel Math Kernel library]), [mkl=yes], [mkl=no])
if test x$mkl = xyes ; then
   AC_DEFINE([GKC_MKL], [], [enables MKL support])
fi

AC_ARG_ENABLE([papi], AC_HELP_STRING([--enable-papi],   [compile with PAPI support]), [papi=yes], [papi=no])
if test x$papi = xyes ; then
   AC_DEFINE([GKC_PAPI], [], [enables papi support])
fi

#################  Set Paths ###################################

AC_ARG_WITH(hdf5dir, AC_HELP_STRING([--with-hdf5dir=prefix],
  [use specified prefix install directory of hdf5 library]),
  DIR_HDF5="$with_hdf5dir", DIR_HDF5="$DIR_HDF5DIR")
AC_SUBST(DIR_HDF5)

AC_ARG_WITH(petscdir, AC_HELP_STRING([--with-petscdir=prefix],
  [use specified prefix install directory of PETSc library]),
  DIR_PETSC="$with_petscdir", DIR_PETSC="$DIR_PETSC")
AC_SUBST(DIR_PETSC)

AC_ARG_WITH(slepcdir,AC_HELP_STRING([--with-slepcdir=prefix],
    [use specified prefix install directory of SLEPc library]),
  DIR_SLEPC="$with_slepcdir", DIR_SLEPC="$DIR_SLEPC")
AC_SUBST(DIR_SLEPC)


AC_ARG_WITH(mkldir, AC_HELP_STRING( [--with-mkldir=prefix],
  [use specified prefix install directory of MKL library]),
  DIR_MKL="$with_mkldir", DIR_MKL="$DIR_MKLDIR" )
AC_SUBST(DIR_MKL)


AC_ARG_WITH(papidir, AC_HELP_STRING([--with-papidir=prefix],
    [use specified prefix install directory of papi library]),
  DIR_PAPI="$with_papidir", DIR_PAPI="$DIR_PAPI")
AC_SUBST(DIR_PAPI)


##########################3333
# Libary stuff

AC_ARG_ENABLE(static, AC_HELP_STRING([--enable-static], [Build with gkc using satic libraries.]), [static=yes], [static=no])           
if test x$static = xyes ; then
   AM_CONDITIONAL([STATIC], [test x$static = xyes])
fi


AC_ARG_ENABLE(shared, AC_HELP_STRING([--enable-shared], [Build with gkc using shared libraries.]), [shared=no], [static=yes])           
if test x$shared = xyes ; then
   AM_CONDITIONAL([SHARED], [test x$shared = xyes])
fi


##############################################################################
# Configure parallel/serial run mode
#
mpi=no
openmp=no
opencl=no

AC_ARG_ENABLE(mpi, AC_HELP_STRING([--enable-mpi], [Build with MPI support.]), [mpi=yes], [mpi=no])           
if test x$mpi = xyes ; then
   AC_DEFINE([GKC_PARALLEL], [], [Defined when compiled parallel version])
   AC_DEFINE([GKC_PARALLEL_MPI], [], [Defined when compiled parallel version with MPI])
   AM_CONDITIONAL([MPI], [test x$mpi = xyes])
fi


AC_ARG_ENABLE(openmp, AC_HELP_STRING([--enable-openmp], [Build with OpenMP support.]), [openmp=yes], [openmp=no]) 
if test x$openmp = xyes ; then
   AC_DEFINE([GKC_PARALLEL], [], [Defined when compiled parallel version])
   AC_DEFINE([GKC_PARALLEL_OPENMP], [], [Defined when compiled with OpenMP support])
   AM_CONDITIONAL([OPENMP],  [test x$openmp  = xyes])
fi


AC_ARG_ENABLE(opencl, AC_HELP_STRING([--enable-opencl], [Build with OpenCL support.]), [opencl=yes], [opencl=no]) 
 echo " opencl = $opencl "
if test x$opencl = xyes ; then
   AC_DEFINE([GKC_PARALLEL], [], [Defined when compiled parallel version])
   AC_DEFINE([GKC_PARALLEL_OPENCL], [], [Defined when compiled with OpenCL support])
fi

if [[ "$mpi" == "no" -a "$openmp" == "no" ]] ; then 
    echo -e "SERIAL: no mpi, no openmp" 
   AC_DEFINE([SERIAL], [], [Defined when compiled serial version]) 
fi 

#########################################################################################################
#  Switch between different FFT libraries

################### p3dfft ####################################################
#AC_ARG_ENABLE(p3dfft, AC_HELP_STRING([--enable-p3dfft], [Build with P3DFFT support.]), [p3dfft=yes], [p3dfft=no]) 


################### fftw3 ####################################################

AC_ARG_ENABLE(fftw3, AC_HELP_STRING([--enable-fftw3], [Build with FFTW support.]), [fftw3=yes], [fftw3=no]) 
AC_ARG_WITH(fftw3dir,AC_HELP_STRING([--with-fftw3dir=prefix],
    [use specified prefix install directory of fftw3 library]),
  DIR_FFTW3="$with_fftw3dir", )
AC_SUBST(DIR_FFTW3)

if test x$fftw3 = xyes ; then
   AC_DEFINE([FFTW3], [], [Defined when compiled parallel version])
fi

###############

# Do I need them again ?!
AM_CONDITIONAL([MPI],    [test x$mpi = xyes])
AM_CONDITIONAL([STATIC], [test x$static  = xyes])
AM_CONDITIONAL([SHARED], [test x$shared  = xyes])
AM_CONDITIONAL([OPENMP], [test x$openmp  = xyes])
AM_CONDITIONAL([OPENCL], [test x$opencl  = xyes])


#### Set automake conditionals
AM_CONDITIONAL([MPIP],  [test x$mpiP    = xyes])
AM_CONDITIONAL([PETSC], [test x$petsc   = xyes])
AM_CONDITIONAL([SLEPC], [test x$slepc   = xyes])
AM_CONDITIONAL([MKL],   [test x$mkl     = xyes])
AM_CONDITIONAL([PAPI],  [test x$papi    = xyes])


AM_CONDITIONAL([FFTW3],  [test x$fftw3  = xyes])

##############################################################################
# Set optimization flags for
# different compilers
#


# mpip profiling
AC_ARG_ENABLE(mpiP, AC_HELP_STRING([--enable-mpiP], [Link mpiP for profiling.]), [mpiP=yes], [mpiP=no]) 
AM_CONDITIONAL([MPIP], [test x$mpiP = xyes])

# static linking of intel mpi library, ...
static=no 

# store environment variable CXXFLAGS 
CXXFLAGS_ENV=$CXXFLAGS
if test x$NO_ENV_CXXFLAGS = xtrue ; then
  CXXFLAGS=""
fi

# GNU C++ compiler
if test x$CXX = xg++ || test x$FC = xmpic++.openmpi ; then
  AC_DEFINE([CXX_GCC], , [Defined when compiling with g++])
  
  if test x$NO_ENV_CXXFLAGS = xtrue ; then
    if test x$debug = xyes ; then
      CXXFLAGS+=" -g -ggdb -O3 -Wall" # -Wextra"
    else
      if test x$profile = xyes ; then
        CXXFLAGS+=" -pg -O3"
      else
        CXXFLAGS+=" -O3 -Wall"
      fi
    fi

    if test x$openmp = xyes ; then
      # echo " openmp = $openmp : add -fopenmp to CXXFLAGS"
   CXXFLAGS+=" -fopenmp"
    fi
      CXXFLAGS+=" -std=c++0x"
  fi
fi

## Set CLANG

# Intel C++ compiler
if test x$CXX = xicc || test x$CXX = icc.openmpi ; then
  AC_DEFINE([CXX_ICC], , [Defined when compiling with Intels C++ compiler])
  
  if test x$NO_ENV_CXXFLAGS = xtrue ; then
    if test x$debug = xyes ; then
      CXXFLAGS+=" -g -O0 -Wcheck "
    else
      if test x$profile = xyes ; then
        CXXFLAGS+=" -p "
      else
        CXXFLAGS+=" -fast -vec "
      fi
    fi

    if test x$openmp = xyes ; then
   CXXFLAGS+=" -openmp "
    fi
        CXXFLAGS+="-std=c++0x -Wall"
  fi

  # depending on the distribution, we need use this in order __attribute__
fi


echo " -------Compiler Settings -------"
if test x$static = xyes ; then 
  echo " static: $static (static linking of Intel-MPI library) "
fi
echo " CXX           = $CXX"
echo " CXXFLAGS      = $CXXFLAGS "
echo " CXXFLAGS_USER = $CXXFLAGS_USER "
echo " LDFLAGS       = $LDFLAGS "
echo " LDFLAGS_USER  = $LDFLAGS_USER "
echo " LIBS          = $LIBS "
echo " LIBS_USER     = $LIBS_USER "
echo " --------------------------------"

echo " HDF5  DIR    = $DIR_HDF5  "
echo " PETSC DIR    = $DIR_PETSC "
echo " SLEPC DIR    = $DIR_SLEPC "
echo " MKL   DIR    = $DIR_MKL   "
echo " PAPI  DIR    = $DIR_PAPI  "

AC_CONFIG_FILES([
  Makefile
  src/Makefile 
])
AC_OUTPUT 

# Create output
echo \
===============================================================================
echo -e "Configure Summary: " 
echo -e Compiler:  $FC
if test x$NO_ENV_CXXFLAGS = xtrue ; then
  echo -e Compiler Flags: $CXXFLAGS \(not included from env CXXFLAGS: $CXXFLAGS_ENV\)
else
  echo -e Compiler Flags: $CXXFLAGS \(included from env CXXFLAGS: $CXXFLAGS_ENV\)
fi
echo -e " additional user compiler/linker flags and libraries can be set with"
echo -e " the environment variables: CXXFLAGS_USER, LDFLAGS_USER, LIBS_USER: "
echo -e " e.g.: env CXXFLAGS_USER=\" -O0\" ./configure" 
echo -e "CXXFLAGS_USER = $CXXFLAGS_USER"
echo -e OpenMP: $openmp
echo -e OpenCL: $opencl
echo -e MPI: $mpi
if test x$mpi = xyes ; then 
    echo -e "With MPI " 
fi

echo -e " Type in \" make \" now ... "
