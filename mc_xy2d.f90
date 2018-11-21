program xy2d
  USE SCIFOR
  USE MPI
  implicit none

  integer :: Nx
  integer :: Nsweep
  integer :: Nwarm
  integer :: Nmeas
  integer :: Ntemp
  integer :: seed
  real(8) :: Dtheta
  real(8) :: Temp,rnd
  !
  integer :: MpiComm
  integer :: irank
  integer :: MpiRank
  integer :: MpiSize
  logical :: MpiMaster


  call Init_MPI(MpiComm,.true.)
  MpiSize   = get_size_MPI(MpiComm)
  MpiRank   = get_rank_MPI(MpiComm)
  MpiMaster = get_master_MPI(MpiComm)
  call parse_input_variable(Nx,"Nx","inputXY2d.conf",default=10)
  call parse_input_variable(Nsweep,"Nsweep","inputXY2d.conf",default=100000)
  call parse_input_variable(Nwarm,"Nwarm","inputXY2d.conf",default=1000)
  call parse_input_variable(Nmeas,"Nmeas","inputXY2d.conf",default=100)
  call parse_input_variable(Temp,"Temp","inputXY2d.conf",default=1d0)
  call parse_input_variable(Dtheta,"Dtheta","inputXY2d.conf",default=0.1d0)
  call parse_input_variable(seed,"SEED","inputXY2d.conf",default=2342161)

  if(MpiMaster)call save_input("inputXY2d.conf")

  do irank=0,MpiSize-1
     call Barrier_MPI(MpiComm)
     if(irank==MpiRank)then
        call random_number(rnd)
        seed = int(1000000d0*rnd)
        open(100,file="Seed_Rank"//str(MpiRank,4)//".dat")
        write(100,*)seed
        close(100)
     endif
  end do

  if(MpiMaster)open(unit=100,file="xy_2d.dat",position='append')
  call mersenne_init(seed)
  call MC_xy2D(Nx,Nsweep,Nwarm,Nmeas,100)
  if(MpiMaster)close(100)

  call Finalize_MPI()

contains


  subroutine MC_xy2D(Nx,Nsweep,Nwarm,Nmeas,unit)
    integer                  :: Nx
    integer                  :: Nsweep
    integer                  :: Nwarm
    integer                  :: Nmeas
    integer                  :: unit
    real(8),dimension(Nx,Nx) :: Lattice
    !
    real(8)                  :: theta,rnd_theta
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
    real(8)                  :: E_mean,E_mean_tmp
    real(8)                  :: Esq_mean,Esq_mean_tmp
    !
    real(8)                  :: Mag(2),aMag
    real(8)                  :: M_sum
    real(8)                  :: Msq_sum
    real(8)                  :: M_mean,M_mean_tmp
    real(8)                  :: Msq_mean,Msq_mean_tmp
    !
    real(8)                  :: CV,Chi
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
    ! 
    if(MpiMaster)call start_timer()
    do iter=1,Nsweep
       !
       !Lattice Sweep
       do i=1,Nx
          do j=1,Nx
             !
             theta      = Lattice(i,j)
             rnd_theta  = Dtheta*pi*(2d0*mersenne()-1d0)
             theta_flip = theta + rnd_theta !mod(theta + rnd_theta,pi2)
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
                M_sum   = M_sum + sqrt(dot_product(Mag,Mag))
                Esq_sum = Esq_sum + Ene*Ene
                Msq_Sum = Msq_Sum + dot_product(Mag,Mag)
             endif
             !
          enddo
       enddo
       if(MpiMaster)call eta(iter,Nsweep)
       !
       if(MpiMaster.AND.mod(iter,Nmeas)==0)&
            call animate_Lattice(Lattice,"Lattice_Animated",.false.)
    enddo
    if(MpiMaster)call stop_timer
    !
    E_mean_tmp = E_sum/Nave
    M_mean_tmp = M_sum/Nave
    Esq_mean_tmp = Esq_sum/Nave
    Msq_mean_tmp = Msq_sum/Nave
    !
    call AllReduce_MPI(MpiComm,E_mean_tmp,E_mean);E_mean=E_mean/MpiSize
    call AllReduce_MPI(MpiComm,M_mean_tmp,M_mean);M_mean=M_Mean/MpiSize
    call AllReduce_MPI(MpiComm,Esq_mean_tmp,Esq_mean);Esq_mean=Esq_mean/MpiSize
    call AllReduce_MPI(MpiComm,Msq_mean_tmp,Msq_mean);Msq_mean=Msq_mean/MpiSize
    !
    aMag = abs(M_mean)/Nlat
    Ene = E_mean/Nlat
    Chi = (Msq_Mean - M_mean**2)/Temp/Nlat
    Cv  = (Esq_mean - E_mean**2)/Temp**2/Nlat
    !
    if(MpiMaster)then
       write(*,*)temp,aMag,Ene,Cv,Chi,dble(Nacc)/Nlat/Nsweep,Nave
       write(unit,*)temp,aMag,Ene,Cv,Chi,dble(Nacc)/Nlat/Nsweep,Nave
       !
       call print_Lattice(Lattice,"Lattice")
       call animate_Lattice(Lattice,"Lattice_Animated",.true.)
       !
    endif
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



  subroutine animate_lattice(lattice,pfile,last)
    real(8),dimension(:,:) :: lattice
    character(len=*)       :: pfile
    logical         :: last
    integer                :: N,i,j
    integer,save           :: iter=0
    !
    N=size(lattice,1)
    call assert_shape(lattice,[N,N])
    !
    if(.not.last)then
       iter=iter+1           
       open(200,file=str(pfile)//str(iter,12)//".dat")
       do i=1,N
          do j=1,N
             write(200,*)i,j,to_positive_angle(Lattice(i,j))
          enddo
          write(200,*)""
       enddo
       close(200)
       return
    endif
    open(200,file="plot_"//str(pfile)//".gp")
    write(200,"(A)")"set terminal gif size 450,450 nocrop animate delay 50 enhanced font 'Times-Roman'" 
    write(200,"(A)")"set output '"//str(pfile)//".gif'"
    write(200,"(A)")"set size square"
    write(200,"(A)")"unset key"
    write(200,"(A)")"set xrange [0.5:"//str(N+0.5d0)//"]"
    write(200,"(A)")"set yrange [0.5:"//str(N+0.5d0)//"]"
    write(200,"(A)")"set cbrange [0:2*pi]"
    write(200,"(A)")"set xtics 5"
    write(200,"(A)")"set ytics 5"
    write(200,"(A)")"set cbtics ('0.00'0, '' pi/2, '3.14' pi, '' 3*pi/2, '6.28' 2*pi)"
    write(200,"(A)")"xf(r,phi) = r*cos(phi)"
    write(200,"(A)")"yf(r,phi) = r*sin(phi)"
    write(200,"(A)")"set style arrow 1 head filled size screen 0.01,15,45 fixed lc rgb 'black'"
    write(200,"(A)")"set palette model HSV defined ( 0 0 1 1, 1 1 1 1 )"
    write(200,"(A)")"do for [i=1:"//str(iter)//":1] {"
    write(200,"(A)")"file = sprintf('"//str(pfile)//"%012d.dat',i)"
    write(200,"(A)")"plot file u 1:2:3 with image, file u ($1):($2):(xf(0.5,$3)):(yf(0.5,$3)) with vectors as 1"
    write(200,"(A)")"}"
    close(200)
  end subroutine animate_lattice


  subroutine print_lattice(lattice,pfile)
    real(8),dimension(:,:) :: lattice
    character(len=*)       :: pfile
    integer                :: N,i,j
    !
    N=size(lattice,1)
    call assert_shape(lattice,[N,N])
    !
    open(200,file=str(pfile)//".dat")
    do i=1,N
       do j=1,N
          write(200,*)i,j,to_positive_angle(Lattice(i,j))
       enddo
       write(200,*)""
    enddo
    close(200)
    open(200,file="plot_"//str(pfile)//".gp")
    write(200,"(A)")"set terminal postscript eps enhanced color font 'Times-Roman,18'" 
    write(200,"(A)")"set output '|ps2pdf -dEPSCrop - "//str(pfile)//".pdf'"
    write(200,"(A)")"set size square"
    write(200,"(A)")"unset key"
    write(200,"(A)")"set xrange [0.5:"//str(N+0.5d0)//"]"
    write(200,"(A)")"set yrange [0.5:"//str(N+0.5d0)//"]"
    write(200,"(A)")"set cbrange [0:2*pi]"
    write(200,"(A)")"set xtics 5"
    write(200,"(A)")"set ytics 5"
    write(200,"(A)")"set cbtics ('0'0, '' pi/2, '{/Symbol-Italic p}' pi, '' 3*pi/2, '{/Symbol-Italic 2p}' 2*pi)"
    write(200,"(A)")"xf(r,phi) = r*cos(phi)"
    write(200,"(A)")"yf(r,phi) = r*sin(phi)"
    write(200,"(A)")"set style arrow 1 head filled size screen 0.01,15,45 fixed lc rgb 'black'"
    write(200,"(A)")"set palette model HSV defined ( 0 0 1 1, 1 1 1 1 )"    
    write(200,"(A)")"plot '"//str(pfile)//".dat"//"' u 1:2:3 with image,'"//str(pfile)//".dat"//"' u ($1):($2):(xf(0.5,$3)):(yf(0.5,$3)) with vectors as 1"
    close(200)
  end subroutine print_lattice


  function to_positive_angle(theta) result(angle)
    real(8) :: theta
    real(8) :: angle
    angle = mod(theta,pi2)
    if(angle<0d0)angle=angle+pi2
  end function to_positive_angle


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
