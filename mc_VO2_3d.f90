program vo2_3d
  USE SCIFOR
  USE MPI
  implicit none

  real(8)                            :: Jh
  real(8)                            :: lambda
  integer                            :: Nx
  integer                            :: Nsweep
  integer                            :: Nwarm
  integer                            :: Nmeas
  integer                            :: seed
  integer                            :: MFit1
  integer                            :: Mfit2
  integer                            :: Rflip1
  integer                            :: Rflip2
  real(8)                            :: dx1
  real(8)                            :: dx2
  real(8)                            :: Temp
  !
  integer                            :: i,j,unit
  real(8),dimension(-22:22)          :: X1
  real(8),dimension(-30:29)          :: X2
  real(8),dimension(-22:22,-30:29)   :: EnergyLocal
  real(8)                            :: Emin,Emax
  real(8)                            :: Alat
  real(8)                            :: X1_min,X1_max
  real(8)                            :: X2_min,X2_max
  type(finter2d_type)                :: vo2_ElocInter
  !
  real(8),allocatable,dimension(:)   :: ArrayX1
  real(8),allocatable,dimension(:)   :: ArrayX2
  real(8),allocatable,dimension(:,:) :: VO2_Potential
  !
  !
  integer                            :: irank
  integer                            :: MpiComm
  integer                            :: MpiRank
  integer                            :: MpiSize
  logical                            :: MpiMaster
  real(8),allocatable,dimension(:)   :: MpiSeed

  !
  !
  call Init_MPI(MpiComm,.true.)
  MpiSize   = get_size_MPI(MpiComm)
  MpiRank   = get_rank_MPI(MpiComm)
  MpiMaster = get_master_MPI(MpiComm)
  !
  !
  call parse_input_variable(Jh,"JH","inputVO2.conf",default=1d0)
  call parse_input_variable(Lambda,"Lambda","inputVO2.conf",default=1d0)
  call parse_input_variable(Nx,"Nx","inputVO2.conf",default=20)
  call parse_input_variable(Nsweep,"Nsweep","inputVO2.conf",default=10000)
  call parse_input_variable(Nwarm,"Nwarm","inputVO2.conf",default=1000)
  call parse_input_variable(Nmeas,"Nmeas","inputVO2.conf",default=100)
  call parse_input_variable(Temp,"Temp","inputVO2.conf",default=4d0)
  call parse_input_variable(Mfit1,"MFit1","inputVO2.conf",default=1000)
  call parse_input_variable(Mfit2,"MFit2","inputVO2.conf",default=1000)
  call parse_input_variable(dx1,"dx1","inputVO2.conf",default=0.1d0)
  call parse_input_variable(dx2,"dx2","inputVO2.conf",default=0.1d0)
  call parse_input_variable(seed,"SEED","inputVO2.conf",default=2342161)
  if(MpiMaster)call save_input("inputVO2.conf")
  !
  allocate(MpiSeed(MpiSize))
  call random_number(MpiSeed)
  do irank=1,MpiSize
     call Barrier_MPI(MpiComm)
     if(irank==MpiRank+1)then
        seed = int(Seed*MpiSeed(irank))
        open(100,file="Seed_Rank"//str(irank,4)//".dat")
        write(100,*)seed
        write(*,*)"Irank",irank," seed=",seed
        close(100)
     endif
  end do
  !
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
  !
  Emin = minval(EnergyLocal)
  Emax = maxval(EnergyLocal)
  EnergyLocal = EnergyLocal - Emax
  EnergyLocal = EnergyLocal/(Emax-Emin)
  !
  call splot3d("VO2_x1x2_data.dat",X1,X2,EnergyLocal,nosurface=.true.)
  !
  !
  X1_min = minval(x1)
  X1_max = maxval(x1)
  !
  X2_min = minval(x2)
  X2_max = maxval(x2)
  !
  Alat = max(X1_max-X1_min,X2_max-X2_min)
  !
  call init_finter2d(vo2_ElocInter,X1,X2,EnergyLocal,3)
  !
  !Build the pre-evaluated local potential
  allocate(ArrayX1(Mfit1))
  allocate(ArrayX2(Mfit2))
  allocate(VO2_Potential(Mfit1,Mfit2))
  call build_VO2_Potential(ArrayX1,ArrayX2,VO2_Potential)

  call mersenne_init(seed)
  if(MpiMaster)open(free_unit(unit),file="vo2_2d.dat",position='append')
  call MC_vo2_3D(unit)
  if(MpiMaster)close(unit)
  !
  call delete_finter2d(vo2_ElocInter)
  !
  call Finalize_MPI()




contains




  subroutine build_VO2_potential(array1,array2,potential)
    real(8),dimension(Mfit1)       :: array1
    real(8),dimension(Mfit2)       :: array2
    real(8),dimension(Mfit1,Mfit2) :: potential
    integer                        :: i1,i2
    real(8)                        :: ddx1,ddx2
    !
    write(*,"(A)")"Building up the interpolated VO2-Phonon potential"
    array1 = linspace(X1_min,X1_max,Mfit1,mesh=ddx1)
    array2 = linspace(X2_min,X2_max,Mfit2,mesh=ddx2)
    Rflip1 = nint(dx1/ddx1)
    Rflip2 = nint(dx2/ddx2)
    write(*,"(A,I0,2x,I0)")"Rflips:",Rflip1,Rflip2
    do i1=1,Mfit1
       do i2=1,Mfit2
          potential(i1,i2) = lambda*vo2_elocal(array1(i1),array2(i2))
       enddo
    enddo
    write(*,"(A)")"Done"
    write(*,"(A)")""
  end subroutine build_VO2_potential

  function vo2_elocal(x,y) result(func)
    real(8) :: x,y
    real(8) :: func
    func = finter2d(vo2_ElocInter,x,y)
  end function vo2_elocal



  subroutine Init_Lattice(lattice)
    integer,dimension(Nx,Nx,Nx,2) :: lattice
    integer                       :: i,j,k
    !
    do k=1,Nx
       do j=1,Nx
          do i=1,Nx
             lattice(i,j,k,1) = int_mersenne(1,Mfit1)
             lattice(i,j,k,2) = int_mersenne(1,Mfit2)
             !
          enddo
       enddo
    enddo
  end subroutine Init_Lattice


  function Lattice_Neighbors(i,j,k) result(neigh)
    integer,dimension(6,3) :: neigh
    integer                :: i,j,k
    integer                :: i_sx,i_dx
    integer                :: j_up,j_dw
    integer                :: k_tp,k_bt
    !
    !PBC:
    i_dx = i+1 ;if(i_dx>Nx)i_dx=1
    j_up = j+1 ;if(j_up>Nx)j_up=1
    k_tp = k+i ;if(k_tp>Nx)k_tp=1
    i_sx = i-1 ;if(i_sx<1)i_sx=Nx
    j_dw = j-1 ;if(j_dw<1)j_dw=Nx
    k_bt = k-1 ;if(k_bt<1)k_bt=Nx
    !
    neigh(1,:) = [i_dx,j,k]
    neigh(2,:) = [i,j_up,k]
    neigh(3,:) = [i,j,k_tp]
    neigh(4,:) = [i_sx,j,k]
    neigh(5,:) = [i,j_dw,k]
    neigh(6,:) = [i,j,k_bt]
  end function Lattice_Neighbors


  function Lattice_Eloc(Lattice,i,j,k,flipd_indx) result(E0)
    integer,dimension(Nx,Nx,Nx,2) :: Lattice
    integer                       :: i,j,k
    integer,dimension(2),optional :: flipd_indx
    integer                       :: in
    real(8)                       :: E0,Wfield
    integer                       :: i1,i2
    real(8)                       :: x1,x2
    integer                       :: i1NN,i2NN
    real(8)                       :: x1NN,x2NN
    integer                       :: NstNbor(6,3)
    !
    NstNbor = Lattice_Neighbors(i,j,k)
    !
    if(present(flipd_indx))then
       i1 = flipd_indx(1)
       i2 = flipd_indx(2)
    else
       i1 = Lattice(i,j,k,1)
       i2 = Lattice(i,j,k,2)
    endif
    !
    x1 = arrayX1(i1)
    x2 = arrayX2(i2)
    !
    Wfield = 0d0
    do in=1,6
       i1NN = Lattice(NstNbor(in,1),NstNbor(in,2),NstNbor(in,3),1)
       i2NN = Lattice(NstNbor(in,1),NstNbor(in,2),NstNbor(in,3),2)
       x1NN = arrayX1(i1NN)
       x2NN = arrayX2(i2NN)
       !
       Wfield  = Wfield + Jh*dot_product([x1,x2]-[x1NN,x2NN],[x1,x2]-[x1NN,x2NN])
    enddo
    !
    E0     = Wfield  + ceiling(Alat**2)*VO2_potential(i1,i2)
    !
  end function Lattice_Eloc


  function Lattice_Energy(Lattice) result(Ene)
    integer,dimension(Nx,Nx,Nx,2) :: Lattice
    real(8)                       :: Ene
    integer                       :: i,j,k
    !
    Ene=0d0
    do i=1,Nx
       do j=1,Nx
          do k=1,Nx
             Ene = Ene + Lattice_Eloc(Lattice,i,j,k)
          enddo
       enddo
    enddo
  end function Lattice_Energy


  function Lattice_Magnetization(Lattice) result(Mag)
    integer,dimension(Nx,Nx,Nx,2)    :: Lattice
    real(8)                       :: x1,x2
    real(8)                       :: Mag(2)
    integer                       :: i,j,k
    !
    Mag=0d0
    do i=1,Nx
       do j=1,Nx
          do k=1,Nx
             x1 = arrayX1(Lattice(i,j,k,1))
             x2 = arrayX2(Lattice(i,j,k,2))
             Mag   = Mag + [abs(x1),abs(x2)]
          enddo
       enddo
    enddo
  end function Lattice_Magnetization




  !MAIN MC PROCEDURE:
  subroutine MC_vo2_3D(unit)
    integer                       :: unit
    integer,dimension(Nx,Nx,Nx,2) :: Lattice
    real(8)                       :: rnd(2)
    !
    integer                       :: i1,i2
    real(8)                       :: x1,x2
    integer                       :: i1_flip,i2_flip
    real(8)                       :: x1_flip,x2_flip
    real(8)                       :: E0
    real(8)                       :: Eflip
    real(8)                       :: Ediff
    real(8)                       :: P
    logical                       :: in_bool
    !    
    integer                       :: Nacc
    integer                       :: Nave
    integer                       :: Nbrk
    !
    real(8)                       :: Ene,Mag(2)
    real(8)                       :: CV,Chi
    !
    real(8)                       :: E_sum
    real(8)                       :: Esq_sum
    real(8)                       :: E_mean,E_mean_tmp
    real(8)                       :: Esq_mean,Esq_mean_tmp
    !
    real(8)                       :: aMag
    real(8)                       :: M_sum
    real(8)                       :: Msq_sum
    real(8)                       :: M_mean,M_mean_tmp
    real(8)                       :: Msq_mean,Msq_mean_tmp
    !
    integer                       :: iter,myI,myJ,myK
    integer                       :: i,j,k,Nlat
    !

    Nlat=Nx*Nx*Nx
    !
    call Init_Lattice(Lattice)
    !
    Ene = Lattice_Energy(Lattice)
    Mag = Lattice_Magnetization(Lattice)
    ! !
    Nacc     = 0
    Nave     = 0
    Nbrk     = 0
    E_sum    = 0d0
    M_sum    = 0d0
    Esq_sum  = 0d0
    Msq_sum  = 0d0
    myI = int_mersenne(1,Nx)
    myJ = int_mersenne(1,Nx)
    myK = int_mersenne(1,Nx)
    ! 
    if(MpiMaster)call start_timer()
    do iter=1,Nsweep
       !
       do i=1,Nx
          do j=1,Nx
             do k=1,Nx
                !
                i1 = Lattice(i,j,k,1)
                i2 = Lattice(i,j,k,2)
                x1 = arrayX1(i1)
                x2 = arrayX2(i2)
                !
                call mt_random(rnd)
                i1_flip = i1-Rflip1 + floor(rnd(1)*(2*Rflip1+1))
                i2_flip = i2-Rflip2 + floor(rnd(2)*(2*Rflip2+1))
                !
                in_bool = (i1_flip<Mfit1).AND.(1<i1_flip).AND.(i2_flip<Mfit2).AND.(1<i2_flip)
                if(.not.in_bool)then
                   Nbrk = Nbrk+1
                   cycle
                endif
                !
                x1_flip = arrayX1(i1_flip)
                x2_flip = arrayX2(i2_flip)
                !
                E0      = Lattice_Eloc(Lattice,i,j,k)
                Eflip   = Lattice_Eloc(Lattice,i,j,k,[i1_flip,i2_flip])
                Ediff   = Eflip - E0
                !
                P = exp(-Ediff/Temp)
                !
                if( min(1d0,P) > mersenne())then
                   Ene  = Ene + Ediff
                   Mag  = Mag - [abs(x1),abs(x2)] + [abs(x1_flip),abs(x2_flip)]
                   !
                   Lattice(i,j,k,1) = i1_flip
                   Lattice(i,j,k,2) = i2_flip
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
       enddo
       !
       if(MpiMaster)then
          call eta(iter,Nsweep)
          call Print_Point(Lattice(myI,myJ,myK,:),"mcVO2_PointGif",.false.)
       endif
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
       write(unit,*)temp,aMag,Ene,Cv,Chi,dble(Nacc)/Nlat/Nsweep,Nave,Nbrk
       write(*,"(A,F21.12)")"Temp=",temp
       write(*,"(A,F21.12)")"Lambda =",lambda
       write(*,"(A,I0)")"Alat^2 =",ceiling(Alat**2)
       write(*,"(A,I0)")"Nbrk=",Nbrk
       write(*,"(A,I0)")"Nave=",Nave
       write(*,"(A,I0)")"Nacc=",Nacc
       write(*,"(A,I0)")"Ntot=",(Nsweep-Nwarm)*Nlat
       write(*,"(A,2F21.12)")"Mag =",aMag,aMag/4d0
       write(*,"(A,F21.12)")"Ene =",Ene
       write(*,"(A,F21.12)")"Cv  =",Cv
       write(*,"(A,F21.12)")"Chi =",Chi
       !
       call print_Lattice(Lattice,"mcVO2_Lattice")
       call Print_Point(Lattice(myI,myJ,myK,:),"mcVO2_PointGif",.true.)
    endif
    !
  end subroutine MC_vo2_3D






  subroutine print_lattice(lattice,pfile)
    integer,dimension(Nx,Nx,Nx,2) :: lattice
    character(len=*)              :: pfile
    integer                       :: i,j,k
    real(8)                       :: rho,theta,x1,x2
    !
    open(200,file=str(pfile)//".dat")
    do i=1,Nx
       do j=1,Nx
          do k=1,Nx
             x1   = arrayX1(Lattice(i,j,k,1))
             x2   = arrayX2(Lattice(i,j,k,2))
             rho  = sqrt(x1**2 + x2**2)
             theta= get_theta(x1,x2)
             write(200,*)i,j,k,theta,rho
          enddo
          write(200,*)""
       enddo
       write(200,*)""
    enddo
    close(200)
    open(200,file="plot_"//str(pfile)//".gp")
    write(200,"(A)")"set terminal postscript eps enhanced color font 'Times-Roman'"
    write(200,"(A)")"set size square"
    write(200,"(A)")"unset key"
    !
    write(200,"(A)")"scale=0.6"
    write(200,"(A)")"xf(r,phi) = r*cos(phi)"
    write(200,"(A)")"yf(r,phi) = r*sin(phi)"
    write(200,"(A)")"set style arrow 1 head filled size screen 0.01,15,45 fixed lc rgb 'black'"
    !
    write(200,"(A)")"set output '|ps2pdf -dEPSCrop - "//str(pfile)//"AllPoints.pdf'"
    write(200,"(A)")"set pm3d map"
    write(200,"(A)")"plot 'VO2_x1x2_data.dat' u 1:2:3 with image, '"//str(pfile)//".dat' u (xf($5,$4)):(yf($5,$4)) w p pointtype 7 pointsize 0.5 lc rgb 'white'"
    !
    write(200,"(A)")"set xrange [0.5:"//str(Nx+0.5d0)//"]"
    write(200,"(A)")"set yrange [0.5:"//str(Nx+0.5d0)//"]"
    write(200,"(A)")"set cbrange [0:2*pi]"
    write(200,"(A)")"set xtics 5"
    write(200,"(A)")"set ytics 5"
    write(200,"(A)")"set cbtics ('0'0, '' pi/2, '{/Symbol-Italic p}' pi, '' 3*pi/2, '{/Symbol-Italic 2p}' 2*pi)"
    write(200,"(A)")"set palette model HSV defined ( 0 0 1 1, 1 1 1 1 )"
    write(200,"(A)")"set output '|ps2pdf -dEPSCrop - "//str(pfile)//"Vec.pdf'"
    write(200,"(A)")"k="//str(Nx/2)
    write(200,"(A)")"plot '"//str(pfile)//".dat"//"' every ::k::k u 1:2:4 with image, '"//str(pfile)//".dat' every ::k::k u 1:2:(xf($5*scale,$4)):(yf($5*scale,$4)) with vectors as 1"
    !
    close(200)
    !
  end subroutine print_lattice


  subroutine print_point(point,pfile,last)
    integer,dimension(2) :: point
    character(len=*)     :: pfile
    logical              :: last
    real(8)              :: rho,theta
    integer,save         :: iter=0
    !
    if(.not.last)then
       iter = iter+1
       rho  = sqrt(arrayX1(point(1))**2 + arrayX2(point(2))**2)
       theta= get_theta(arrayX1(point(1)),arrayX2(point(2)))
       open(200,file=str(pfile)//".dat",access='append')
       write(200,*)theta,rho,vo2_potential(point(1),point(2))
       close(200)
       return
    endif
    open(200,file="plot_"//str(pfile)//".gp")
    write(200,"(A)")"set term wxt"
    write(200,"(A)")"#set terminal gif size 450,450 nocrop animate delay 50 enhanced font 'Times-Roman'" 
    write(200,"(A)")"#set output '"//str(pfile)//".gif'"
    write(200,"(A)")"set size square"
    write(200,"(A)")"unset key"
    write(200,"(A)")"xf(r,phi) = r*cos(phi)"
    write(200,"(A)")"yf(r,phi) = r*sin(phi)"
    write(200,"(A)")""
    write(200,"(A)")""
    write(200,"(A)")"##Map plot"
    write(200,"(A)")"set pm3d map"
    write(200,"(A)")"do for [i=1:"//str(iter)//":10] {"
    write(200,"(A)")"set title 'i='.i"   
    write(200,"(A)")"plot 'VO2_x1x2_data.dat' u 1:2:3 with image,'"//str(pfile)//".dat' every ::i::i using (xf($2,$1)):(yf($2,$1)) w p pointtype 7 pointsize 1.5 lc rgb 'white','"//str(pfile)//".dat' every ::1::i using (xf($2,$1)):(yf($2,$1)) w l ls 1 lc rgb 'white'"
    write(200,"(A)")"}"
    write(200,"(A)")""
    write(200,"(A)")""
    write(200,"(A)")"##Surface plot"
    write(200,"(A)")"#set xrange [-2.2000:2.2000]"
    write(200,"(A)")"#set yrange [-2.8000:2.8000]"
    write(200,"(A)")"#set zrange [-0.2:0]"
    write(200,"(A)")"#do for [i=1:"//str(iter)//":1] {"
    write(200,"(A)")"#set title 'i='.i"   
    write(200,"(A)")"#splot 'VO2_x1x2_data.dat' u 1:2:3 with pm3d,'mcVO2_PointGif.dat' every ::i::i using (xf($2,$1)):(yf($2,$1)):3 w p pointtype 7 pointsize 1.5 lc rgb 'white', 'mcVO2_PointGif.dat' every ::1::i using (xf($2,$1)):(yf($2,$1)):3 w l ls 1 lc rgb 'white'"
    write(200,"(A)")"#}"
    write(200,"(A)")""
    write(200,"(A)")""
    close(200)
  end subroutine print_point






  function get_theta(x1,x2) result(theta)
    real(8) :: x1,x2
    real(8) :: theta,thetap
    thetap = atan(abs(x2/x1))
    if(x1>=0d0.AND.x2>=0d0)then
       theta = thetap           !Q1
    elseif(x1<0d0.AND.x2>0d0)then
       theta = pi - thetap      !Q2
    elseif(x1<0d0.AND.x2<0d0)then
       theta = pi + thetap      !Q3
    elseif(x1>0d0.AND.x2<0d0)then
       theta = 2*pi - thetap    !Q4
    end if
  end function get_theta


end program vo2_3d





