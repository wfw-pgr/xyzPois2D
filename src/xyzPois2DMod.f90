module xyzPois2DMod
  use PBiCGStabMod
  ! --------------------   USAGE   ------------------------- !
  ! source :: include boundary condition of x.           --- !
  ! --- e.g.) x=phi_b at (i=1) --> source(1,j) = phi_b   --- !
  ! -------------------------------------------------------- !
contains
  
  ! =================================================================== !
  ! ===  BiCGSTABCrt2D_MPI  :: BiCGSTAB solver for Cartesian 2D     === !
  ! =================================================================== !
  subroutine BiCGSTABCrt2D_MPI( x, source, dx, dy, LI, LJ )
    implicit none
    include 'mpif.h'
    integer         , intent(in)    :: LI, LJ
    double precision, intent(in)    :: dx, dy
    double precision, intent(in)    :: source(LI,LJ)
    double precision, intent(inout) :: x(LI,LJ)
    integer                         :: ierr, iPE, LIbuf, LJbuf, Nitem_buf
    integer                         :: LIloc, LJloc, LIxLJ
    logical, save                   :: Flag__MakeMatrix = .true.
    double precision, allocatable   :: xbuf(:,:)
    integer         , parameter     :: Fr_=1, To_=2
    
    ! ------------------------------------- !
    ! --- [1] Make Matrix if No Matrix  --- !
    ! ------------------------------------- !
    if ( Flag__MakeMatrix ) then
       ! -- [1-1] variables for comm.   --  !
       call DefineCommunication( LI, LJ )
       LIloc = Nloc
       LJloc = Mloc
       LIxLJ = NxM
       ! -- [1-2] Make Matrix           --  !
       allocate( xloc(LIloc,LJloc), sloc(LIloc,LJloc) )
       allocate( Amat(LIxLJ,Ncomm), pt(LIxLJ,Ncomm),  &
            &    Nelem(LIxLJ)     , Minv(LIxLJ)       )
       xloc(:,:)  = 0.d0
       sloc(:,:)  = 0.d0
       Amat(:,:)  = 0.d0
       pt  (:,:)  = 0.d0
       Nelem(:)   = 0.d0
       Minv (:)   = 0.d0
       call XYLaplacian( dx, dy, LIloc, LJloc, LIxLJ )
       Flag__MakeMatrix = .false.
    endif
    
    ! ------------------------------------- !
    ! --- [2] call PBiCGSTAB_Engine     --- !
    ! ------------------------------------- !
    sloc(:,:) = source( FromTo(myRank,Fr_):FromTo(myRank,To_), 1:LJloc )
    call PBiCGSTAB_Engine_MPI( sloc, xloc )
    ! call MPI_Finalize( ierr )
    ! stop
    ! ------------------------------------- !
    ! --- [3] BroadCast Solution        --- !
    ! ------------------------------------- !
    do iPE=0, PEtot-1
       ! -- [3-1] Prepare iPE's Solution -- !
       LJbuf      = LJloc
       LIbuf      = FromTo(iPE,2) - FromTo(iPE,1) + 1
       Nitem_buf  = LIbuf * LJbuf
       allocate( xbuf(LIbuf,LJbuf) )
       if ( myRank.eq.iPE ) xbuf = xloc
       ! -- [3-2] BroadCast              -- !
       call MPI_Bcast( xbuf, Nitem_buf     , MPI_DOUBLE_PRECISION, &
            &          iPE , MPI_COMM_WORLD, ierr )
       ! -- [3-3] Save in x(From:To)     -- !
       x( FromTo(iPE,1):FromTo(iPE,2), 1:LJloc ) = xbuf(1:LIbuf,1:LJbuf)
       deallocate( xbuf )
    enddo
    
    return
  end subroutine BiCGSTABCrt2D_MPI


  ! =================================================================== !
  ! ===  XYLaplacian  :: Calculate Laplacian ( cartesian )          === !
  ! =================================================================== !
  subroutine XYLaplacian( dx, dy, LI, LJ, LIxLJ )
    implicit none
    include 'mpif.h'
    integer         , intent(in)  :: LI, LJ, LIxLJ
    double precision, intent(in)  :: dx, dy
    integer                       :: i, j, k, kL, kR
    logical                       :: Flag_SizeERROR
    double precision              :: coef(5)

    ! ------------------------------------------------------ !
    ! --- [1] Preparation                                --- !
    ! ------------------------------------------------------ !
    !  -- [1-1] Display                 --  !
    if ( myRank.eq.0 ) then
       write(6,'(2x,a)'       ) '[ XYLaplacian  @ myPBiCGSTAB ]'
       write(6,'(5x,2(a,i12))') '* A-Matrix :: ', LIxLJ, ' x ', LIxLJ
    endif
    !  -- [1-2] Coefficient             --  !
    coef(1) = 1.d0 / dy**2
    coef(2) = 1.d0 / dx**2
    coef(3) = - 2.d0 * ( 1.d0 / dx**2 + 1.d0 / dy**2 )
    coef(4) = 1.d0 / dx**2
    coef(5) = 1.d0 / dy**2
    !  -- [1-3] Initialization          --  !
    Amat(:,:) = 0.d0
    pt(:,:)   = 0
    k         = 1
    Minv(:)   = 1.d0
    
    ! ------------------------------------------------------ !
    ! --- [2] Make Matrix                                --- !
    ! ------------------------------------------------------ !
    !  -- [2-1] j=1 Boundary            --  !
    do i=1, LI
       Amat(k,1)   = 1.d0
       pt(k,1)     = k
       Nelem(k)    = 1
       k = k+1
    enddo
    !  -- [2-2] Main Region (j=2~LJ-1)  --  !
    do j=2, LJ-1
       do i=1, LI
          ! - Laplacian - !
          Amat(k,1)   = coef(1)
          Amat(k,2)   = coef(2)
          Amat(k,3)   = coef(3)
          Amat(k,4)   = coef(4)
          Amat(k,5)   = coef(5)
          pt(k,1)     = k-LI
          pt(k,2)     = k-1
          pt(k,3)     = k
          pt(k,4)     = k+1
          pt(k,5)     = k+LI
          Nelem(k)    = 5
          k = k+1
       enddo
    enddo
    !  -- [2-3] j=LJ Boundary           --  !
    do i=1, LI
       Amat(k,1) = 1.d0
       pt(k,1)   = k
       Nelem(k)  = 1
       k         = k+1
    enddo

    ! ------------------------------------------------------ !
    ! --- [3] Boundary ( x=xMin, x=xMax ) ( parallel )   --- !
    ! ------------------------------------------------------ !
    if ( PEtot.gt.1 ) then
       if ( ( myRank.gt.0 ).and.( myRank.lt.PEtot-1 ) ) then
          !  -- (Left) Neibour PE       --  !
          kL       = 1
          do j=1, LJ
             pt(kL,2)     = k
             kL           = kL + LI
             k            = k  + 1
          enddo
          !  -- (Right) Neibour PE      --  !
          kR       = LI
          do j=1, LJ
             pt(kR,4)     = k
             kR           = kR + LI
             k            = k  + 1
          enddo
       endif
       if ( myRank.eq.0 ) then
          !  -- (Left) Dirichlet        --  !
          kL = 1
          do j=1, LJ
             Amat(kL,1)   = 1.d0
             pt(kL,1)     = kL
             Nelem(kL)    = 1
             Minv(kL)     = 1.d0 
             kL           = kL + LI
          enddo
          !  -- (Right) Neibour PE      --  !
          kR = LI
          do j=1, LJ
             pt(kR,4)     = k
             kR           = kR + LI
             k            = k  + 1
          enddo
       endif
       if ( myRank.eq.PEtot-1 ) then
          !  -- (Left) Neibour PE       --  !
          kL = 1
          do j=1, LJ
             pt(kL,2)     = k
             kL           = kL + LI
             k            = k  + 1
          enddo
          !  -- (Right) Dirichlet       --  !
          kR = LI
          do j=1, LJ
             Amat(kR,1)   = 1.d0
             pt(kR,1)     = kR
             Nelem(kR)    = 1
             Minv(kR)     = 1.d0 
             kR           = kR + LI
          enddo
       endif
    endif
    
    ! ------------------------------------------------------ !
    ! --- [4] Boundary ( x=xMin, x=xMax ) ( single )     --- !
    ! ------------------------------------------------------ !
    if ( PEtot.eq.1 ) then
       !  -- Left  B.C. ( Dirichlet )   --  !
       kL = 1
       do j=1, LJ
          Amat(kL,1)   = 1.d0
          pt(kL,1)     = kL
          Nelem(kL)    = 1
          Minv(kL)     = 1.d0 
          kL           = kL + LI
       enddo
       !  -- Right B.C. ( Dirichlet )   --  !
       kR = LI
       do j=1, LJ
          Amat(kR,1)   = 1.d0
          pt(kR,1)     = kR
          Nelem(kR)    = 1
          Minv(kR)     = 1.d0 
          kR           = kR + LI
       enddo
    endif
    
    ! ------------------------------------------------------ !
    ! --- [5] Size Check                                 --- !
    ! ------------------------------------------------------ !
    Flag_SizeERROR = .false.
    if ( PEtot.gt.1 ) then
       if ( ( myRank.eq.0 ).or.( myRank.eq.PEtot-1 ) ) then
          if ( (k-1).ne.(LIxLJ+  LJ) ) Flag_SizeERROR = .true.
       else
          if ( (k-1).ne.(LIxLJ+2*LJ) ) Flag_SizeERROR = .true.
       endif
    endif
    if ( PEtot.eq.1 ) then
       if    ( (k-1).ne.(LIxLJ    ) ) Flag_SizeERROR = .true.
    endif
    if ( Flag_SizeERROR ) then
       write(6,*) ' [WARNING] Size and k_count are Incompatable.'
       write(6,*) "      myRank / PEtot :: ", myRank, " / ", PEtot
       write(6,*) "           LIxLJ+2LJ :: ", LIxLJ+2*LJ
       write(6,*) "           LIxLJ+ LJ :: ", LIxLJ+  LJ
       write(6,*) "            LIxLJ    :: ", LIxLJ
       write(6,*) "                 k-1 :: ", k-1
       stop
    endif

    return
  end subroutine XYLaplacian
  
  
end module xyzPois2DMod
