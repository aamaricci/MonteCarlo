program xy2d
  USE SCIFOR
  implicit none

  integer :: Nx
  integer :: Nsweep
  integer :: Nwarm
  integer :: Nmeas
  integer :: Ntemp
  integer :: seed
  real(8) :: DTheta
  real(8) :: Temp
  !
  call parse_input_variable(Nx,"Nx","inputXY2d.conf",default=10)
  call parse_input_variable(Nsweep,"Nsweep","inputXY2d.conf",default=100000)
  call parse_input_variable(Nwarm,"Nwarm","inputXY2d.conf",default=1000)
  call parse_input_variable(Nmeas,"Nmeas","inputXY2d.conf",default=100)
  call parse_input_variable(Temp,"Temp","inputXY2d.conf",default=1d0)
  call parse_input_variable(dtheta,"dtheta","inputXY2d.conf",default=0.1d0,comment="Range of variation for theta in unit of 2pi")
  call parse_input_variable(seed,"SEED","inputXY2d.conf",default=2342161)
  call save_input("inputXY2d.conf")


  open(unit=100,file="xy_2d.dat",position='append')
  call mersenne_init(seed)
  call MC_xy2D(Nx,Nsweep,Nwarm,Nmeas,100)
  close(100)

contains


  subroutine MC_xy2D(Nx,Nsweep,Nwarm,Nmeas,unit)
    integer                  :: Nx
    integer                  :: Nsweep
    integer                  :: Nwarm
    integer                  :: Nmeas
    integer                  :: unit
    real(8),dimension(Nx,Nx) :: Lattice
    !
    real(8)                  :: theta
    real(8)                  :: theta_flip
    real(8)                  :: E0
    real(8)                  :: Eflip
    real(8)                  :: Ediff
    real(8)                  :: P
    !
    integer                  :: Nacc
    integer                  :: Nave
    !
    real(8)                  :: Ene
    real(8)                  :: E_sum
    real(8)                  :: Esq_sum
    real(8)                  :: E_mean
    real(8)                  :: Esq_mean
    !
    real(8)                  :: Mag(2)
    real(8)                  :: M_sum(2)
    real(8)                  :: Msq_sum(2)
    real(8)                  :: M_mean(2)
    real(8)                  :: Msq_mean(2)
    !
    real(8)                  :: CV,Chi(2)
    integer                  :: iter
    integer                  :: i,j,k,N,Nlat

    !
    !
    Nlat=Nx*Nx
    !
    call Init_Lattice(Lattice)
    !
    Ene = Lattice_Energy(Lattice)
    Mag = Lattice_Magnetization(Lattice)
    ! !
    Nacc = 0
    Nave = 0
    E_sum   = 0d0
    M_sum   = 0d0
    Esq_sum = 0d0
    Msq_sum = 0d0    
    ! !
    call start_timer()
    do iter=1,Nsweep
       !
       !Lattice Sweep
       do i=1,Nx
          do j=1,Nx
             !
             theta      = Lattice(i,j)
             theta_flip = theta + pi2*dtheta*(2d0*mersenne()-1d0)                          
             E0         = Lattice_Eloc(Lattice,i,j)
             Eflip      = Lattice_Eloc(Lattice,i,j,theta_flip)
             Ediff      = Eflip - E0
             !
             P = exp(-Ediff/Temp)
             !
             if( min(1d0,P) > mersenne()) then
                Ene  = Ene + Ediff
                Mag  = Mag - [cos(theta),sin(theta)] + [cos(theta_flip),sin(theta_flip)]
                !
                Lattice(i,j) = theta_flip
                !
                Nacc = Nacc + 1
             end if
             !
             if(iter>Nwarm.AND.mod(iter,Nmeas)==0)then
                Nave    = Nave + 1
                E_sum   = E_sum + Ene
                M_sum   = M_sum + Mag
                Esq_sum = Esq_sum + Ene*Ene
                Msq_Sum = Msq_Sum + Mag*Mag
             endif
             !
          enddo
       enddo
       call eta(iter,Nsweep)
       !
    enddo
    call stop_timer
    !
    E_mean = E_sum/Nave
    M_mean = M_sum/Nave
    Esq_mean = Esq_sum/Nave
    Msq_mean = Msq_sum/Nave
    !
    Mag = abs(M_mean)/Nlat
    Ene = E_mean/Nlat
    Chi = (Msq_Mean - M_mean**2)/Temp/Nlat
    Cv  = (Esq_mean - E_mean**2)/Temp**2/Nlat
    !
    write(*,*)temp,Mag,Ene,Cv,Chi,dble(Nacc)/Nlat/Nsweep,Nave
    write(unit,*)temp,Mag,Ene,Cv,Chi,dble(Nacc)/Nlat/Nsweep,Nave
    do i=1,Nx
       do j=1,Nx
          write(200,*)i,j,Lattice(i,j)
       enddo
    enddo
    !
  end subroutine MC_Xy2D




  subroutine Init_Lattice(lattice)
    real(8),dimension(:,:) :: lattice
    integer                :: N,i,j
    !
    N=size(lattice,1)
    call assert_shape(lattice,[N,N])
    !
    do j=1,N
       do i=1,N
          lattice(i,j) = pi2*mersenne()
       enddo
    enddo
  end subroutine Init_Lattice


  function Lattice_Neighbors(Lattice,i,j) result(neigh)
    real(8),dimension(:,:) :: Lattice
    integer,dimension(4,2) :: neigh
    integer                :: i,j,N,k
    integer                :: i_sx,i_dx
    integer                :: j_up,j_dw
    !
    N=size(Lattice,1)
    !
    !PBC:
    i_dx = i+1 ;if(i_dx>N)i_dx=1
    j_up = j+1 ;if(j_up>N)j_up=1
    i_sx = i-1 ;if(i_sx<1)i_sx=N
    j_dw = j-1 ;if(j_dw<1)j_dw=N
    !
    neigh(1,:) = [i_dx,j]
    neigh(2,:) = [i,j_up]
    neigh(3,:) = [i_sx,j]
    neigh(4,:) = [i,j_dw]
  end function Lattice_Neighbors


  function Lattice_Eloc(Lattice,i,j,theta_) result(E0)
    real(8),dimension(:,:) :: Lattice
    integer                :: i,j
    real(8),optional       :: theta_
    integer                :: k,N
    real(8)                :: E0,ThetaNN,theta
    integer                :: NstNbor(4,2)
    !
    NstNbor = Lattice_Neighbors(Lattice,i,j)
    !
    if(.not.present(theta_))then
       theta = Lattice(i,j)
    else
       theta = theta_
    endif
    E0    = 0d0
    do k=1,4
       thetaNN = Lattice(NstNbor(k,1),NstNbor(k,2))
       E0 = E0 - 1d0*cos(thetaNN - theta)
    enddo
  end function Lattice_Eloc


  function Lattice_Energy(Lattice) result(Ene)
    real(8),dimension(:,:) :: Lattice
    real(8)                :: Ene
    integer                :: i,j,N
    !
    N=size(Lattice,1)
    call assert_shape(Lattice,[N,N])
    !
    Ene=0d0
    do i=1,N
       do j=1,N
          Ene = Ene + Lattice_Eloc(Lattice,i,j)
       enddo
    enddo
  end function Lattice_Energy


  function Lattice_Magnetization(Lattice) result(Mag)
    real(8),dimension(:,:) :: Lattice
    real(8)                :: Mag(2),theta
    integer                :: i,j,N
    !
    N=size(Lattice,1)
    call assert_shape(Lattice,[N,N])
    !
    Mag=0d0
    do i=1,N
       do j=1,N
          Theta = Lattice(i,j)
          Mag   = Mag + [cos(theta),sin(theta)]
       enddo
    enddo
  end function Lattice_Magnetization


end program xy2d







! subroutine Init_Probability(Probability)
!   real(8),dimension(5) :: Probability
!   integer               :: N,i,j,k
!   Probability(5) = exp(-2*4/Temp)
!   Probability(4) = exp(-2*2/Temp)
!   Probability(3) = exp(-2*0/Temp)
!   Probability(2) = exp( 2*2/Temp)
!   Probability(1) = exp( 2*4/Temp)
!   if( any(isnan(Probability)) ) stop "P undefined: overflow in Exp(-DeltaE/Temp)"
! end subroutine Init_Probability
