program main
  use variablesMod
  use ioUtilityMod
  use xyzPois2DMod, only : BiCGStabCrt2D_MPI
  implicit none
  include "mpif.h"
  integer :: i,j

  call MPI_Init     ( ierr )
  call MPI_Comm_Rank( MPI_COMM_WORLD, myRank, ierr )
  call MPI_Comm_Size( MPI_COMM_WORLD, PEtot , ierr )

  write(6,*)
  write(6,*) "[main.f90] myRank == ", myRank
  write(6,*) "[main.f90] PEtot  == ", PEtot
  write(6,*)
  
  call load__source_and_boundary  
  call BiCGSTABCrt2D_MPI( xvec, source, dx, dy, LI, LJ )
  call save__results
  call MPI_Finalize ( ierr )
  
end program main
  
