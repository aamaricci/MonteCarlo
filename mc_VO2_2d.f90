program vo2_2d
  USE SCIFOR
  USE MPI
  implicit none

  real(8)                          :: Jh
  integer                          :: Nx
  integer                          :: Nsweep
  integer                          :: Nwarm
  integer                          :: Nmeas
  integer                          :: Ntemp
  integer                          :: seed
  real(8)                          :: Drho,RhoLength
  real(8)                          :: Dtheta
  real(8)                          :: Temp,rnd
  !
  integer                          :: i,j,unit
  real(8),dimension(-22:22)        :: X1
  real(8),dimension(-30:29)        :: X2
  real(8),dimension(-22:22,-30:29) :: EnergyLocal
  type(finter2d_type)              :: vo2_ElocInter
  !
  integer                          :: MpiComm
  integer                          :: irank
  integer                          :: MpiRank
  integer                          :: MpiSize
  logical                          :: MpiMaster
  !
  call Init_MPI(MpiComm,.true.)
  MpiSize   = get_size_MPI(MpiComm)
  MpiRank   = get_rank_MPI(MpiComm)
  MpiMaster = get_master_MPI(MpiComm)
  !
  call parse_input_variable(Jh,"JH","inputXY2d.conf",default=1d0)
  call parse_input_variable(Nx,"Nx","inputXY2d.conf",default=10)
  call parse_input_variable(Nsweep,"Nsweep","inputXY2d.conf",default=100000)
  call parse_input_variable(Nwarm,"Nwarm","inputXY2d.conf",default=1000)
  call parse_input_variable(Nmeas,"Nmeas","inputXY2d.conf",default=100)
  call parse_input_variable(Temp,"Temp","inputXY2d.conf",default=1d0)
  call parse_input_variable(Drho,"Drho","inputXY2d.conf",default=0.1d0)
  call parse_input_variable(Dtheta,"Dtheta","inputXY2d.conf",default=0.1d0)
  call parse_input_variable(seed,"SEED","inputXY2d.conf",default=2342161)
  if(MpiMaster)call save_input("inputXY2d.conf")
  !
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
  !
  x1=0d0
  x2=0d0
  open(free_unit(unit),file="Data_Tot_VO2.conf")
  do i=0,22
     do j=0,29
        read(unit,*)x1(i),x2(j),EnergyLocal(i,j)
     enddo
  enddo
  close(unit)
  forall(i=1:22)x1(-i) = -x1(i)
  forall(j=1:30)x2(-j) = -x2(j-1)
  !mirror w/ to Y==X2 axis: x1-->-X1
  forall(i=1:22,j=0:29)EnergyLocal(-i,j) = EnergyLocal(i,j)
  !mirror w/ to X==X1 axis: x2-->-X2
  forall(i=0:22,j=1:30)EnergyLocal(i,-j) = EnergyLocal(i,j-1)
  !mirror w/ to Origin==(0,0) axis: x1-->-X1 && x2-->-X2
  forall(i=1:22,j=1:30)EnergyLocal(-i,-j) = EnergyLocal(i,j-1)

  call splot3d("VO2_x1x2_data.conf",X1,X2,EnergyLocal)

  call init_finter2d(vo2_ElocInter,X1,X2,EnergyLocal,3)

  call mersenne_init(seed)
  if(MpiMaster)open(100,file="vo2_2d.dat",position='append')
  call MC_vo2_2D(100)
  if(MpiMaster)close(100)

  call delete_finter2d(vo2_ElocInter)



  call Finalize_MPI()




contains



  function vo2_elocal(x,y) result(func)
    real(8) :: x,y
    real(8) :: func
    func = finter2d(vo2_ElocInter,x,y)
  end function vo2_elocal



  subroutine MC_vo2_2D(unit)
    integer                       :: unit
    real(8),dimension(Nx,Nx,2)    :: Lattice
    !
    real(8)                       :: rho,rnd_rho
    real(8)                       :: rho_flip
    !
    real(8)                       :: theta,rnd_theta
    real(8)                       :: theta_flip
    !
    real(8)                       :: v0(2),vf(2)
    real(8)                       :: E0
    real(8)                       :: Eflip
    real(8)                       :: Ediff
    real(8)                       :: P
    !
    integer                       :: Nacc
    integer                       :: Nave
    !
    real(8)                       :: Ene,Mag(2)
    real(8)                       :: CV,Chi
    real(8)                       :: Ang,Len
    !
    real(8)                       :: E_sum
    real(8)                       :: Esq_sum
    real(8)                       :: E_mean,E_mean_tmp
    real(8)                       :: Esq_mean,Esq_mean_tmp
    !
    real(8)                       :: aTheta,aRho
    real(8)                       :: A_sum,A_mean,A_mean_tmp
    real(8)                       :: L_sum,L_mean,L_mean_tmp
    !
    real(8)                       :: aMag
    real(8)                       :: M_sum
    real(8)                       :: Msq_sum
    real(8)                       :: M_mean,M_mean_tmp
    real(8)                       :: Msq_mean,Msq_mean_tmp
    !

    integer                       :: iter
    integer                       :: i,j,k,N,Nlat
    !
    !
    Nlat=Nx*Nx
    !
    call Init_Lattice(Lattice)
    !
    Ene = Lattice_Energy(Lattice)
    Mag = Lattice_Magnetization(Lattice)
    Ang = sum(Lattice(:,:,1))
    Len = sum(Lattice(:,:,2))
    ! !
    Nacc     = 0
    Nave     = 0
    E_sum    = 0d0
    M_sum    = 0d0
    A_sum    = 0d0
    L_sum    = 0d0
    Esq_sum  = 0d0
    Msq_sum  = 0d0    
    ! 
    if(MpiMaster)call start_timer()
    do iter=1,Nsweep
       !
       ! !Lattice Sweep
       ! do i=1,Nx
       !    do j=1,Nx
       i = int_mersenne(1,Nx)
       j = int_mersenne(1,Nx)
       !
       theta = Lattice(i,j,1)
       rho   = Lattice(i,j,2)
       !
       rnd   = 2d0*mersenne()-1d0 
       !
       if( mersenne() < 0.5d0)then !flip theta
          theta_flip = theta + Dtheta*pi*rnd
          rho_flip   = rho
       else
          theta_flip = theta
          rho_flip   = max(0d0,rho + Drho*rnd)
       endif
       !
       E0         = Lattice_Eloc(Lattice,i,j)
       Eflip      = Lattice_Eloc(Lattice,i,j,[theta_flip,rho_flip])
       Ediff      = Eflip - E0
       !
       P = exp(-Ediff/Temp)
       !
       if( min(1d0,P) > mersenne()) then
          Ene  = Ene + Ediff
          Mag  = Mag - rho*[cos(theta),sin(theta)] + rho_flip*[cos(theta_flip),sin(theta_flip)]
          Ang  = Ang - theta + theta_flip
          Len  = Len - rho   + rho_flip
          !
          Lattice(i,j,1) = theta_flip
          Lattice(i,j,2) = rho_flip
          !
          Nacc = Nacc + 1
       end if
       !
       if(iter>Nwarm.AND.mod(iter,Nmeas)==0)then
          Nave    = Nave + 1
          E_sum   = E_sum + Ene
          M_sum   = M_sum + sqrt(dot_product(Mag,Mag))
          A_sum   = A_sum + Ang
          L_sum   = L_sum + Len
          Esq_sum = Esq_sum + Ene*Ene
          Msq_Sum = Msq_Sum + dot_product(Mag,Mag)
       endif
       !
       !    enddo
       ! enddo
       if(MpiMaster)call eta(iter,Nsweep)
       !
       if(MpiMaster.AND.mod(iter,Nmeas)==0)&
            call animate_Lattice(Lattice,"Lattice_Animated",.false.)
    enddo
    if(MpiMaster)call stop_timer
    !
    E_mean_tmp = E_sum/Nave
    M_mean_tmp = M_sum/Nave
    A_mean_tmp = A_sum/Nave
    L_mean_tmp = L_sum/Nave
    Esq_mean_tmp = Esq_sum/Nave
    Msq_mean_tmp = Msq_sum/Nave
    !
    call AllReduce_MPI(MpiComm,E_mean_tmp,E_mean);E_mean=E_mean/MpiSize
    call AllReduce_MPI(MpiComm,M_mean_tmp,M_mean);M_mean=M_Mean/MpiSize
    call AllReduce_MPI(MpiComm,A_mean_tmp,A_mean);A_mean=A_Mean/MpiSize
    call AllReduce_MPI(MpiComm,L_mean_tmp,L_mean);L_mean=L_Mean/MpiSize
    call AllReduce_MPI(MpiComm,Esq_mean_tmp,Esq_mean);Esq_mean=Esq_mean/MpiSize
    call AllReduce_MPI(MpiComm,Msq_mean_tmp,Msq_mean);Msq_mean=Msq_mean/MpiSize
    !
    aMag = abs(M_mean)/Nlat
    Ene = E_mean/Nlat
    Ang = A_mean/Nlat
    Len = L_mean/Nlat
    Chi = (Msq_Mean - M_mean**2)/Temp/Nlat
    Cv  = (Esq_mean - E_mean**2)/Temp**2/Nlat
    !
    if(MpiMaster)then
       write(unit,*)temp,aMag,Ene,Ang,Len,Cv,Chi,dble(Nacc)/Nlat/Nsweep,Nave
       write(*,"(A5,F21.12)")"Temp=",temp
       write(*,"(A5,F21.12)")"Mat =",aMag
       write(*,"(A5,F21.12)")"Ene =",Ene
       write(*,"(A5,F21.12)")"<a> =",Ang
       write(*,"(A5,F21.12)")"<r> =",Len
       write(*,"(A5,F21.12)")"<X1>=",Len*Cos(Ang)
       write(*,"(A5,F21.12)")"<X2>=",Len*Sin(Ang)
       write(*,"(A5,F21.12)")"Cv  =",Cv
       write(*,"(A5,F21.12)")"Chi =",Chi
       !
       call print_Lattice(Lattice,"Lattice")
       call animate_Lattice(Lattice,"Lattice_Animated",.true.)
    endif
    !
  end subroutine MC_Vo2_2D






  subroutine Init_Lattice(lattice)
    real(8),dimension(Nx,Nx,2)    :: lattice
    integer                       :: N,i,j
    !
    do j=1,Nx
       do i=1,Nx
          lattice(i,j,1) = pi2*mersenne()
          lattice(i,j,2) = 1d0!mersenne()
       enddo
    enddo
  end subroutine Init_Lattice


  function Lattice_Neighbors(Lattice,i,j) result(neigh)
    real(8),dimension(Nx,Nx,2)    :: Lattice
    integer,dimension(4,2)        :: neigh
    integer                       :: i,j,k
    integer                       :: i_sx,i_dx
    integer                       :: j_up,j_dw
    !
    !PBC:
    i_dx = i+1 ;if(i_dx>Nx)i_dx=1
    j_up = j+1 ;if(j_up>Nx)j_up=1
    i_sx = i-1 ;if(i_sx<1)i_sx=Nx
    j_dw = j-1 ;if(j_dw<1)j_dw=Nx
    !
    neigh(1,:) = [i_dx,j]
    neigh(2,:) = [i,j_up]
    neigh(3,:) = [i_sx,j]
    neigh(4,:) = [i,j_dw]
  end function Lattice_Neighbors


  function Lattice_Eloc(Lattice,i,j,vec) result(E0)
    real(8),dimension(Nx,Nx,2)    :: Lattice
    integer                       :: i,j
    real(8),dimension(2),optional :: vec
    integer                       :: k
    real(8)                       :: E0
    real(8)                       :: Theta,Rho
    real(8)                       :: ThetaNN,RhoNN
    integer                       :: NstNbor(4,2)
    !
    NstNbor = Lattice_Neighbors(Lattice,i,j)
    !
    if(present(vec))then
       Theta = vec(1)
       Rho   = vec(2)
    else
       Theta = Lattice(i,j,1)
       Rho   = Lattice(i,j,2)
    endif
    !
    E0  =  vo2_elocal(rho*cos(theta),rho*sin(theta))
    do k=1,4
       ThetaNN = Lattice(NstNbor(k,1),NstNbor(k,2),1)
       RhoNN   = Lattice(NstNbor(k,1),NstNbor(k,2),2)
       E0 = E0 - Jh*RhoNN*Rho*cos(thetaNN - theta)
    enddo
  end function Lattice_Eloc


  function Lattice_Energy(Lattice) result(Ene)
    real(8),dimension(Nx,Nx,2)    :: Lattice
    real(8)                       :: Ene
    integer                       :: i,j
    !
    Ene=0d0
    do i=1,Nx
       do j=1,Nx
          Ene = Ene + Lattice_Eloc(Lattice,i,j)
       enddo
    enddo
  end function Lattice_Energy


  function Lattice_Magnetization(Lattice) result(Mag)
    real(8),dimension(Nx,Nx,2)    :: Lattice
    real(8)                       :: Mag(2),theta,rho
    integer                       :: i,j
    !
    Mag=0d0
    do i=1,Nx
       do j=1,Nx
          Theta = Lattice(i,j,1)
          Rho   = Lattice(i,j,2)
          Mag   = Mag + Rho*[cos(theta),sin(theta)]
       enddo
    enddo
  end function Lattice_Magnetization




  subroutine animate_lattice(lattice,pfile,last)
    real(8),dimension(Nx,Nx,2) :: lattice
    character(len=*)         :: pfile
    logical                  :: last
    integer                  :: i,j
    integer,save             :: iter=0
    !
    if(.not.last)then
       iter=iter+1           
       open(200,file=str(pfile)//str(iter,12)//".dat")
       do i=1,Nx
          do j=1,Nx
             write(200,*)i,j,to_positive_angle(Lattice(i,j,1))
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
    write(200,"(A)")"set xrange [0.5:"//str(Nx+0.5d0)//"]"
    write(200,"(A)")"set yrange [0.5:"//str(Nx+0.5d0)//"]"
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
    real(8),dimension(Nx,Nx,2) :: lattice
    character(len=*)           :: pfile
    integer                    :: i,j
    !
    open(200,file=str(pfile)//".dat")
    do i=1,Nx
       do j=1,Nx
          write(200,*)i,j,to_positive_angle(Lattice(i,j,1))
       enddo
       write(200,*)""
    enddo
    close(200)
    open(200,file="plot_"//str(pfile)//".gp")
    write(200,"(A)")"set size square"
    write(200,"(A)")"unset key"
    write(200,"(A)")"set xrange [0.5:"//str(Nx+0.5d0)//"]"
    write(200,"(A)")"set yrange [0.5:"//str(Nx+0.5d0)//"]"
    write(200,"(A)")"set cbrange [0:2*pi]"
    write(200,"(A)")"set xtics 5"
    write(200,"(A)")"set ytics 5"
    write(200,"(A)")"set cbtics pi/2"
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



end program vo2_2d


