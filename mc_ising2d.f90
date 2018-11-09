program ising2d
  USE SCIFOR
  implicit none

  integer :: Nx
  integer :: Nsweep
  integer :: Nwarm
  integer :: Nmeas
  integer :: Ntemp
  integer :: seed
  real(8) :: Temp
  !
  call parse_input_variable(Nx,"Nx","inputISING2d.conf",default=10)
  call parse_input_variable(Nsweep,"Nsweep","inputISING2d.conf",default=100000)
  call parse_input_variable(Nwarm,"Nwarm","inputISING2d.conf",default=1000)
  call parse_input_variable(Nmeas,"Nmeas","inputISING2d.conf",default=100)
  call parse_input_variable(Temp,"Temp","inputISING2d.conf",default=1d0)
  call parse_input_variable(seed,"SEED","inputISING2d.conf",default=2342161)
  call save_input("inputISING2d.conf")


  open(unit=100,file="ising_2d.dat",position='append')
  call mersenne_init(seed)
  call MC_Ising(Nx,Nsweep,Nwarm,Nmeas,100)
  close(100)

contains


  subroutine MC_Ising(Nx,Nsweep,Nwarm,Nmeas,unit)
    integer                  :: Nx
    integer                  :: Nsweep
    integer                  :: Nwarm
    integer                  :: Nmeas
    integer                  :: unit
    integer,dimension(Nx,Nx) :: Ising
    !
    real(8)                  :: Ene,Mag,CV,Chi
    real(8)                  :: P,Ediff
    integer                  :: Spin,WeissField
    integer                  :: iter
    integer                  :: i,j,k,N,Nlat
    integer                  :: Nacc,Nave
    !
    real(8)                  :: M_sum,Msq_sum
    real(8)                  :: E_sum,Esq_sum
    real(8)                  :: M_mean,Msq_mean
    real(8)                  :: E_mean,Esq_mean
    !
    !
    Nlat=Nx*Nx
    !
    call Init_Ising(Ising)
    !
    Ene = Ising_Energy(Ising)
    Mag = sum(Ising)
    !
    Nacc = 0
    Nave = 0
    E_sum   = 0d0
    M_sum   = 0d0
    Esq_sum = 0d0
    Msq_sum = 0d0    
    !
    do iter=1,Nsweep
       !
       !Lattice Sweep
       do i=1,Nx
          do j=1,Nx
             !
             Spin       = Ising(i,j)
             WeissField = Ising_WF(Ising,i,j)
             Ediff      = 2d0*Spin*WeissField
             P          = exp(-Ediff/Temp)
             !      
             if( min(1d0,P) > mersenne()) then  !flip-spin: ACCEPT
                Ising(i,j) = -Spin
                Ene  = Ene + 2d0*Spin*WeissField
                Mag  = Mag - 2d0*Spin
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
       !
    enddo
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
    write(unit,*)temp,Mag,Ene,Cv,Chi,dble(Nacc)/Nlat/Nsweep,Nave
    !
  end subroutine MC_Ising



  subroutine Init_Ising(Ising)
    integer,dimension(:,:) :: Ising
    integer                :: N,i,j
    !
    N=size(Ising,1)
    call assert_shape(Ising,[N,N])
    !
    do j=1,N
       do i=1,N
          Ising(i,j) = sgn(2*mersenne()-1)
       enddo
    enddo
  end subroutine Init_Ising


  function Ising_Neighbors(Ising,i,j) result(neigh)
    integer,dimension(:,:) :: Ising
    integer,dimension(4,2) :: neigh
    integer                :: i,j,N,k
    integer                :: i_sx,i_dx
    integer                :: j_up,j_dw
    N=size(Ising,1)
    call assert_shape(Ising,[N,N])
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
  end function Ising_Neighbors


  function Ising_WF(Ising,i,j) result(WF)
    integer,dimension(:,:) :: Ising
    integer                :: i,j,k,N
    integer                :: WF
    integer                :: NstNbor(4,2)
    NstNbor = Ising_Neighbors(Ising,i,j)
    WF = 0
    do k=1,4
       WF = WF + Ising(NstNbor(k,1),NstNbor(k,2))
    enddo
  end function Ising_WF


  !Get energy = -Sx*WF!/2
  function Ising_Energy(Ising) result(Ene)
    integer,dimension(:,:) :: Ising
    real(8)                :: Ene
    integer                :: Spin
    integer                :: WeissField
    integer                :: i,j,N
    !
    N=size(Ising,1)
    call assert_shape(Ising,[N,N])
    !
    Ene=0d0
    do i=1,N
       do j=1,N
          Spin       = Ising(i,j)
          WeissField = Ising_WF(Ising,i,j)
          Ene = Ene - Spin*WeissField/2d0
       enddo
    enddo
  end function Ising_Energy




end program ising2d







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
