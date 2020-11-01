module ioUtilityMod
contains

  
  ! ====================================================== !
  ! === load__source_and_boundary                      === !
  ! ====================================================== !
  subroutine load__source_and_boundary
    use variablesMod
    implicit none
    integer :: i, j, LK
    character(cLen) :: cmt

    ! ------------------------------------------------------ !
    ! --- [1] load Amatrix                               --- !
    ! ------------------------------------------------------ !
    
    open(lun,file=trim(srcFile),status="old")
    read(lun,*)
    read(lun,*)
    read(lun,*) cmt, LK, LJ, LI, nCmp
    allocate( source(LI,LJ), xvec(LI,LJ), xg(LI,LJ), yg(LI,LJ) )
    xvec(:,:) = 0.d0
    do j=1, LJ
       do i=1, LI
          read(lun,*) xg(i,j), yg(i,j), source(i,j)
       enddo
    enddo
    close(lun)
    write(6,*) "[load__source_and_boundary] srcFile is loaded...."


    ! ------------------------------------------------------ !
    ! --- [2] source term clean up                       --- !
    ! ------------------------------------------------------ !
    ! do j=2, LJ-1
    !    do i=2, LI-1
    !       source(i,j) = 0.d0
    !    enddo
    ! enddo

    ! ------------------------------------------------------ !
    ! --- [3] dx, dy                                     --- !
    ! ------------------------------------------------------ !
    dx = xg(2,1) - xg(1,1)
    dy = yg(1,2) - yg(1,1)
    
    return
  end subroutine load__source_and_boundary
  
  

  ! ====================================================== !
  ! === save results in a file                         === !
  ! ====================================================== !
  subroutine save__results
    use variablesMod
    implicit none
    integer :: i, j

    ! ------------------------------------------------------ !
    ! --- [1] save results                               --- !
    ! ------------------------------------------------------ !
    open(lun,file=trim(rslFile),status="replace")
    do j=1, LJ
       do i=1, LI
          write(lun,*) xg(i,j), yg(i,j), xvec(i,j)
       enddo
    enddo
    close(lun)
    write(6,*) "[save__results] rslFile is saved...."
    
    return
  end subroutine save__results
  
end module ioUtilityMod
