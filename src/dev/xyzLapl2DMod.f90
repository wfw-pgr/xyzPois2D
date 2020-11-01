module xyzLapl2DMod
contains

  ! ====================================================== !
  ! === generate xyz Cartesian Laplacian (2D)          === !
  ! ====================================================== !
  
  subroutine generate__xyzLaplacian2D( Amat, dx, dy, LI, LJ, LIxLJ )
    implicit none
    integer         , intent(in)    :: LI, LJ, LIxLJ
    double precision, intent(in)    :: dx, dy
    double precision, intent(inout) :: Amat(LIxLJ,), pt(LIxLJ,5), Nelem(LIxLJ)
    double precision                :: coef(5)
    
    ! ------------------------------------------------------ !
    ! --- [1] preparation                                --- !
    ! ------------------------------------------------------ !
    !  -- [1-1] Initialization          --  !
    Amat(:,:) = 0.d0
    pt(:,:)   = 0
    k         = 1
    Minv(:)   = 1.d0
    
    !  -- [1-2] Coefficient             --  !
    coef(1)   = 1.d0 / dy**2
    coef(2)   = 1.d0 / dx**2
    coef(3)   = - 2.d0 * ( 1.d0 / dx**2 + 1.d0 / dy**2 )
    coef(4)   = 1.d0 / dx**2
    coef(5)   = 1.d0 / dy**2

    ! ------------------------------------------------------ !
    ! --- [2] Make Matrix                                --- !
    ! ------------------------------------------------------ !
    
    !  -- [2-1] j=1 Boundary            --  !
    do i=1, LI
       A(k,1)      = 1.d0
       pt(k,1)     = k
       Nelement(k) = 1
       k = k+1
    enddo
    !  -- [2-2] Main Region (j=2~LJ-1)  --  !
    do j=2, LJ-1
       do i=1, LI
          ! - Laplacian - !
          A(k,1)      = coef(1)
          A(k,2)      = coef(2)
          A(k,3)      = coef(3)
          A(k,4)      = coef(4)
          A(k,5)      = coef(5)
          pt(k,1)     = k-LI
          pt(k,2)     = k-1
          pt(k,3)     = k
          pt(k,4)     = k+1
          pt(k,5)     = k+LI
          Nelement(k) = 5
          k = k+1
       enddo
    enddo
    !  -- [2-3] j=LJ Boundary           --  !
    do i=1, LI
       A(k,1) = 1.d0
       pt(k,1) = k
       Nelement(k) = 1
       k = k+1
    enddo



    
    return
  end subroutine generate__xyzLaplacian2D
  
end module xyzLapl2DMod
