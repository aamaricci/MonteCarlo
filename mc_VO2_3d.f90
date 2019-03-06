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
  integer                            :: Npoints
  integer                            :: NprintLat
  real(8)                            :: dx1
  real(8)                            :: dx2
  real(8)                            :: Temp,p_global
  character(len=100)                 :: TempFile
  logical                            :: wPoint
  logical                            :: wPDF1d
  logical                            :: wPDF2d
  !
  real(8),dimension(:),allocatable   :: TempList
  integer                            :: i,j,unit,it,TempLen
  real(8),dimension(-22:22)          :: X1
  real(8),dimension(-30:29)          :: X2
  real(8),dimension(-22:22,-30:29)   :: EnergyLocal
  real(8)                            :: Emin,Emax
  real(8)                            :: Alat
  real(8)                            :: X1_min,X1_max
  real(8)                            :: X2_min,X2_max
  type(finter2d_type)                :: vo2_ElocInter
  logical                            :: iobool
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
  call parse_input_variable(p_global,"P_GLOBAL","inputVO2.conf",default=0.5d0)
  call parse_input_variable(Nsweep,"Nsweep","inputVO2.conf",default=5,comment="In units of 10**3")
  call parse_input_variable(Nwarm,"Nwarm","inputVO2.conf",default=1000)
  call parse_input_variable(Nmeas,"Nmeas","inputVO2.conf",default=100)
  call parse_input_variable(Npoints,"Npoints","inputVO2.conf",default=10)
  call parse_input_variable(NprintLat,"NprintLat","inputVO2.conf",default=10)
  call parse_input_variable(Temp,"Temp","inputVO2.conf",default=10d0)
  call parse_input_variable(TempFile,"TempFile","inputVO2.conf",default="list_temp.in")
  call parse_input_variable(Mfit1,"MFit1","inputVO2.conf",default=1000)
  call parse_input_variable(Mfit2,"MFit2","inputVO2.conf",default=1000)
  call parse_input_variable(dx1,"dx1","inputVO2.conf",default=0.1d0)
  call parse_input_variable(dx2,"dx2","inputVO2.conf",default=0.1d0)
  call parse_input_variable(seed,"SEED","inputVO2.conf",default=2342161)
  call parse_input_variable(wPoint,"WPOINT","inputVO2.conf",default=.false.)
  call parse_input_variable(wPDF1d,"wPDF1d","inputVO2.conf",default=.false.)
  call parse_input_variable(wPDF2d,"wPDF2d","inputVO2.conf",default=.false.)
  if(MpiMaster)call save_input("inputVO2.conf")
  if(MpiMaster)call print_input()
  call set_store_size(10)
  !
  if(Mfit1<=size(X1).OR.Mfit2<=size(X2))stop "Error: Mfit1 <= size(X1)=45 OR Mfit2 <= size(X2)=60 "
  !
  Nsweep=Nsweep*10**3
  !
  allocate(MpiSeed(MpiSize))
  call random_number(MpiSeed)
  do irank=0,MpiSize-1
     call Barrier_MPI(MpiComm)
     if(irank==MpiRank)then
        seed = int(Seed*MpiSeed(irank+1))
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
  !
  !Init the MT rng
  call mersenne_init(seed)
  !
  !Read Temps
  inquire(file=str(TempFile),exist=IObool)
  if(IObool)then
     if(MpiMaster)write(*,"(A)")'Reading Temperatures from file'//str(TempFile)
     TempLen = file_length(str(TempFile))
     open(free_unit(unit),file=str(TempFile))
     allocate(TempList(TempLen))
     do it=1,TempLen
        read(unit,*)TempList(it)
        if(MpiMaster)write(*,"(A)")"Temp="//str(TempList(it))
     enddo
     close(unit)
  else
     TempLen=1
     allocate(TempList(TempLen))
     TempList(1) = Temp
  endif
  if(MpiMaster)write(*,"(A)")""  

  !Start doing MC here:
  if(MpiMaster)open(free_unit(unit),file="vo2_3d.dat",position='append')
  do it=1,TempLen
     Temp = TempList(it)
     if(MpiMaster)write(*,"(A)")"Doing Temp="//str(Temp)
     call MC_vo2_3D(unit)
  enddo
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
    real(8)                        :: ddx1,ddx2,x,y
    !
    write(*,"(A)")"Building up the interpolated VO2-Phonon potential"
    array1 = linspace(X1_min,X1_max,Mfit1,mesh=ddx1)
    array2 = linspace(X2_min,X2_max,Mfit2,mesh=ddx2)
    Rflip1 = nint(dx1/ddx1)
    Rflip2 = nint(dx2/ddx2)
    write(*,"(A,I0,2x,I0)")"Rflips:",Rflip1,Rflip2
    do i1=1,Mfit1
       do i2=1,Mfit2
          x = array1(i1)
          y = array2(i2)
          potential(i1,i2) = ceiling(Alat**2)*lambda*finter2d(vo2_ElocInter,x,y)
       enddo
    enddo
    write(*,"(A)")"Done"
    write(*,"(A)")""
  end subroutine build_VO2_potential



  subroutine Init_Lattice(lattice)
    integer,dimension(Nx,Nx,Nx,2) :: lattice
    integer                       :: i,j,k,unit
    integer                       :: aI,bI
    logical                       :: bool
    !
    inquire(file="lattice_rank"//str(MpiRank,4)//".restart",exist=bool)
    if(bool)then
       write(*,*)"Reading Lattice from lattice_rank"//str(MpiRank,4)//".restart"
       open(free_unit(unit),file="lattice_rank"//str(MpiRank,4)//".restart")
       do i=1,Nx
          do j=1,Nx
             do k=1,Nx
                read(unit,*)aI,bI
                lattice(i,j,k,1) = aI
                lattice(i,j,k,2) = bI
             enddo
          enddo
       enddo
       close(unit)
    else
       do i=1,Nx
          do j=1,Nx
             do k=1,Nx
                lattice(i,j,k,1) = mt_uniform(1,Mfit1)
                lattice(i,j,k,2) = mt_uniform(1,Mfit2)
             enddo
          enddo
       enddo
    endif
  end subroutine Init_Lattice



  subroutine Save_Lattice(lattice)
    integer,dimension(Nx,Nx,Nx,2) :: lattice
    integer                       :: i,j,k,unit
    !
    open(free_unit(unit),file="lattice_rank"//str(MpiRank,4)//".restart")
    do i=1,Nx
       do j=1,Nx
          do k=1,Nx
             write(unit,*)lattice(i,j,k,1),lattice(i,j,k,2)
          enddo
       enddo
    enddo
    close(unit)
  end subroutine Save_Lattice



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
    k_tp = k+1 ;if(k_tp>Nx)k_tp=1
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
    integer                       :: ii,jj,kk
    real(8)                       :: x1,x2
    integer                       :: i1NN,i2NN
    real(8)                       :: x1NN,x2NN
    real(8)                       :: rho,rhoNN
    real(8)                       :: Ui(2),Uj(2)
    real(8)                       :: theta,thetaNN
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
    x1   = arrayX1(i1)
    x2   = arrayX2(i2)
    Ui   = [x1,x2]
    !
    Wfield = 0d0
    do in=1,6
       ii     = NstNbor(in,1)
       jj     = NstNbor(in,2)
       kk     = NstNbor(in,3)
       i1NN   = Lattice(ii,jj,kk,1)
       i2NN   = Lattice(ii,jj,kk,2)
       x1NN   = arrayX1(i1NN)
       x2NN   = arrayX2(i2NN)
       Uj     = [x1NN,x2NN]
       !
       Wfield  = Wfield + dot_product(Ui-Uj,Ui-Uj)
    enddo
    !
    E0     = Jh*Wfield + VO2_potential(i1,i2)
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
             Mag   = Mag + [x1,x2]
          enddo
       enddo
    enddo
  end function Lattice_Magnetization




  !MAIN MC PROCEDURE:
  subroutine MC_vo2_3D(unit)
    integer                        :: unit
    integer,dimension(Nx,Nx,Nx,2)  :: Lattice
    real(8)                        :: rnd(2)
    !
    integer                        :: i1,i2
    real(8)                        :: x1,x2
    integer                        :: i1_flip,i2_flip
    real(8)                        :: x1_flip,x2_flip
    real(8)                        :: E0
    real(8)                        :: Eflip
    real(8)                        :: Ediff
    real(8)                        :: P,ran
    logical                        :: in_bool
    !    
    integer                        :: Nacc
    integer                        :: Nave
    integer                        :: Nbrk
    !
    real(8)                        :: Ene,Mag(2)
    real(8)                        :: CV,Chi
    !
    real(8)                        :: E_sum
    real(8)                        :: Esq_sum
    real(8)                        :: E_mean
    real(8)                        :: Esq_mean
    !
    real(8)                        :: M_sum
    real(8)                        :: Msq_sum
    real(8)                        :: M_mean
    real(8)                        :: Msq_mean
    !
    real(8)                        :: X1_sum,X2_sum
    real(8)                        :: X1sq_sum,X2sq_sum
    real(8)                        :: X1_mean,X2_mean
    real(8)                        :: X1sq_mean,X2sq_mean
    real(8)                        :: sX1,sX2
    !
    integer                        :: iter
    integer,dimension(Npoints)     :: myI,myJ,myK
    integer                        :: i,ii,j,k,Nlat
    !
    real(8)                        :: data,Emin,Emax,ene_sigma
    integer,parameter              :: Npdf1=500,Npdf2=100
    real(8),dimension(Npdf1)       :: ene_pdf_tmp
    real(8),dimension(Npdf2,Npdf2) :: lat_pdf_tmp
    type(pdf_kernel)               :: ene_pdf
    type(pdf_kernel_2d)            :: lat_pdf
    real(8),dimension(2,2)         :: lat_sigma
    !
    !
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
    E_mean   = 0d0
    M_mean   = 0d0
    Esq_mean = 0d0
    Msq_mean = 0d0
    !
    X1_sum = 0d0
    X2_sum = 0d0
    !
    do i=1,Npoints
       myI(i) = mt_uniform(1,Nx)
       myJ(i) = mt_uniform(1,Nx)
       myK(i) = mt_uniform(1,Nx)
    enddo
    !
    Emin =  huge(1d0)
    Emax = -huge(1d0)
    !
    if(MpiMaster)call start_timer()
    MCsweep: do iter=1,Nsweep
       !
#ifdef _SITE
       i = mt_uniform(1,Nx)
       j = mt_uniform(1,Nx)
       k = mt_uniform(1,Nx)
#else
       iloop:do i=1,Nx
          jloop:do j=1,Nx
             kloop:do k=1,Nx
#endif
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
                if(mersenne()>p_global)then
                   ran = mersenne()
                   if(ran<=1d0/3d0)then
                      i1_flip = Mfit1-i1_flip+1
                   elseif(ran>2d0/3d0)then
                      i2_flip = Mfit2-i2_flip+1
                   else
                      i1_flip = Mfit1-i1_flip+1
                      i2_flip = Mfit2-i2_flip+1
                   endif
                endif
                !
                in_bool = (i1_flip<Mfit1).AND.(1<i1_flip).AND.(i2_flip<Mfit2).AND.(1<i2_flip)
                if(.not.in_bool)then
                   Nbrk = Nbrk+1
                   cycle kloop
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
                if( min(1d0,P) > mersenne() )then
                   Ene  = Ene + Ediff
                   Mag  = Mag - [x1,x2] + [x1_flip,x2_flip]
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
                   X1_sum  = X1_sum + Mag(1)
                   X2_sum  = X2_sum + Mag(2)
                   Esq_sum = Esq_sum + Ene*Ene
                   Msq_Sum = Msq_Sum + dot_product(Mag,Mag)
                   X1sq_sum= X1sq_sum + Mag(1)*Mag(1)
                   X2sq_sum= X2sq_sum + Mag(2)*Mag(2)
                   if(wPDF1d)write(999-MpiRank,*)Ene/Nlat
                   if(Ene/Nlat < Emin) Emin=Ene/Nlat
                   if(Ene/Nlat > Emax) Emax=Ene/Nlat
                   !
                endif
                !
#ifdef _SITE
                !
#else
             enddo kloop
          enddo jloop
       enddo iloop
#endif
       !
       if(MpiMaster)then
          call eta(iter,Nsweep)
       endif
       if(wPoint)then
          if(NprintLat>1.AND.mod(iter,Nsweep/NprintLat)==0)&
               call print_Lattice(Lattice,"mc_Lattice_Temp"//str(Temp)//"_"//str(MpiRank,4)//"_iter"//str(iter,12)//".dat")
          do ii=1,Npoints
             call Print_Point(Lattice(myI(ii),myJ(ii),myK(ii),:),&
                  "mc_PointGif_Temp"//str(Temp)//"_"//str(MpiRank,4)//"_"//str(ii,2)//".dat",.false.)
             !
             call Print_Angle(Lattice(myI(ii),myJ(ii),myK(ii),:),&
                  "mc_AngleDynamics_Temp"//str(Temp)//"_"//str(MpiRank,4)//"_"//str(ii,2)//".dat")
          enddo
       endif
    enddo MCsweep
    if(MpiMaster)call stop_timer
    !
    call Save_Lattice(Lattice)
    !
    call AllReduce_MPI(MpiComm,E_sum,E_mean);E_mean=E_mean/MpiSize/Nave/Nlat
    call AllReduce_MPI(MpiComm,M_sum,M_mean);M_mean=M_Mean/MpiSize/Nave/Nlat
    call AllReduce_MPI(MpiComm,X1_sum,X1_mean);X1_mean=X1_Mean/MpiSize/Nave/Nlat
    call AllReduce_MPI(MpiComm,X2_sum,X2_mean);X2_mean=X2_Mean/MpiSize/Nave/Nlat
    call AllReduce_MPI(MpiComm,Esq_sum,Esq_mean);Esq_mean=Esq_mean/MpiSize/Nave/Nlat/Nlat
    call AllReduce_MPI(MpiComm,Msq_sum,Msq_mean);Msq_mean=Msq_mean/MpiSize/Nave/Nlat/Nlat
    !
    Chi = (Msq_Mean - M_mean**2)/Temp
    Cv  = (Esq_mean - E_mean**2)/Temp/Temp
    sX1 = (X1sq_Mean - X1_mean**2)/Temp
    sX2 = (X2sq_Mean - X2_mean**2)/Temp
    !
    if(MpiMaster)then
       write(unit,*)temp,M_mean,X1_mean,X2_mean,E_mean,Cv,Chi,sX1,sX2,Nave,Nsweep,Nx
       write(*,"(A,I0)")     "Nx=",Nx
       write(*,"(A,I0)")     "Na=",Nave
       write(*,"(A,I0)")     "Ns=",Nsweep
       write(*,"(A,I0)")     "Nb=",Nbrk
       write(*,"(A,F21.12)") "T =",temp
       write(*,"(A,F21.12)") "M =",M_mean
       write(*,"(A,F21.12)") "X1=",X1_mean
       write(*,"(A,F21.12)") "X2=",X2_mean
       write(*,"(A,F21.12)") "E =",E_mean
       write(*,"(A,F21.12)") "C =",Cv
       write(*,"(A,F21.12)") "X =",Chi
    endif
    !
    call print_Lattice(Lattice,"mc_Lattice_Temp"//str(Temp))
    !
    if(wPoint)then
       do i=1,Npoints
          call Print_Point(Lattice(myI(i),myJ(i),myK(i),:),&
               "mc_PointGif_Temp"//str(Temp)//"_"//str(MpiRank,4)//"_"//str(i,2)//".dat",.true.)
       enddo
    endif
    !
    if(wPDF1d)then
       Emin = Emin-10d0
       Emax = Emax+10d0
       call Bcast_MPI(MpiComm,Emin)
       call Bcast_MPI(MpiComm,Emax)
       !
       call pdf_allocate(ene_pdf,Npdf1)
       call pdf_set_range(ene_pdf,Emin,Emax)
       call pdf_sigma(ene_pdf,sqrt(Esq_mean-E_mean**2),Nave,ene_sigma)
       call pdf_push_sigma(ene_pdf,ene_sigma)
       !
       rewind(999-MpiRank)
       if(MpiMaster)call start_timer()
       do i=1,Nave
          read(999-Mpirank,*)data
          call pdf_accumulate(ene_pdf,data)
       enddo
       call system("rm -fv fort."//str(999-MpiRank))
       if(MpiMaster)call stop_timer()
       !
       call pdf_normalize(ene_pdf)
       ene_pdf_tmp=0d0
       call AllReduce_MPI(MpiComm,ene_pdf%pdf,ene_pdf_tmp);ene_pdf_tmp=ene_pdf_tmp/MpiSize
       ene_pdf%pdf = ene_pdf_tmp
       !
       if(MpiMaster)then
          call pdf_print(ene_pdf,"mc_PDF_E_Temp"//str(Temp)//".dat")
          call pdf_save(ene_pdf,"mc_PDF_E_Temp"//str(Temp)//".save")
          !
          call pdf_print_moments(ene_pdf,"mc_Moments_PDF_E_Temp"//str(Temp)//".dat")
          write(*,"(A)")""
       endif
       !
       call pdf_deallocate(ene_pdf)
    endif
    !
    if(wPDF2d)then
       call pdf_allocate(lat_pdf,[Npdf2,Npdf2])
       call pdf_set_range(lat_pdf,[-3d0,-3d0],[3d0,3d0])
       ! lat_sigma = reshape([sqrt(Msq_mean-M_mean**2),0d0,0d0,sqrt(Msq_mean-M_mean**2)],[2,2])
       lat_sigma = reshape([0.02d0,0d0,0d0,0.02d0],[2,2])
       call pdf_push_sigma(lat_pdf,lat_sigma)
       !
       if(MpiMaster)call start_timer()
       do i=1,Nx
          do j=1,Nx
             do k=1,Nx
                mag(1)   = arrayX1(Lattice(i,j,k,1))
                mag(2)   = arrayX2(Lattice(i,j,k,2))
                call pdf_accumulate(lat_pdf,mag)
             enddo
          enddo
       enddo
       if(MpiMaster)call stop_timer()
       !
       call pdf_normalize(lat_pdf)
       lat_pdf_tmp=0d0
       call AllReduce_MPI(MpiComm,lat_pdf%pdf,lat_pdf_tmp);lat_pdf_tmp=lat_pdf_tmp/MpiSize
       lat_pdf%pdf = lat_pdf_tmp
       !
       if(MpiMaster)then
          call pdf_print(lat_pdf,"mc_PDF_X_Temp"//str(Temp)//".dat")
          call pdf_save(lat_pdf,"mc_PDF_X_Temp"//str(Temp)//".save")
          write(*,"(A)")""
       endif
       call pdf_deallocate(lat_pdf)
    endif
    !
  end subroutine MC_vo2_3D






  subroutine print_lattice(lattice,pfile)
    integer,dimension(Nx,Nx,Nx,2) :: lattice
    character(len=*)              :: pfile
    character(len=200)            :: pname
    integer                       :: i,j,k
    real(8)                       :: rho,theta,x1,x2
    !
    pname=str(pfile)//"_"//str(MpiRank,4)//".dat"
    open(200,file=str(pname))
    do k=1,Nx
       do i=1,Nx
          do j=1,Nx
             x1   = arrayX1(Lattice(i,j,k,1))
             x2   = arrayX2(Lattice(i,j,k,2))
             write(200,*)i,j,k,X1,X2,get_theta(x1,x2)
          enddo
       enddo
       write(200,*)
    enddo
    close(200)
    call file_gzip(str(pname))
    !
    if(MpiMaster)then
       open(200,file="plot_"//str(pfile)//".gp")
       write(200,"(A)")"reset"
       write(200,"(A)")"set terminal postscript eps enhanced color font 'Times-Roman'"
       write(200,"(A)")"reset"
       write(200,"(A)")""
       write(200,"(A)")"RANK='"//str(MpiRank,4)//"'"
       write(200,"(A)")"PFILE='"//str(pfile)//"_'.RANK.'.dat'"
       write(200,"(A)")"OUT_ALLPOINTS=PFILE.'_AllPoints.pdf'"
       write(200,"(A)")"OUT_VECPDF=PFILE.'_Vec.pdf'"
       write(200,"(A)")"OUT_VECGIF=PFILE.'_Vec.gif'"
       write(200,"(A)")"OUT_VEC3D=PFILE.'_3d.pdf'"
       write(200,"(A)")""
       write(200,"(A)")"system('gunzip -v '.PFILE.'.gz')"
       write(200,"(A)")""
       write(200,"(A)")"set size square"
       write(200,"(A)")"unset key"
       write(200,"(A)")"set style arrow 1 head filled size screen 0.005,15,45 fixed lc rgb 'black'"
       write(200,"(A)")""
       write(200,"(A)")"set output '|ps2pdf -dEPSCrop - '.OUT_ALLPOINTS"
       write(200,"(A)")"set pm3d map"
       write(200,"(A)")"plot 'VO2_x1x2_data.dat' u 1:2:3 with image, PFILE u 4:5 w p pointtype 7 pointsize 0.5 lc rgb 'white'"
       write(200,"(A)")"unset view"
       write(200,"(A)")"unset pm3d"
       write(200,"(A)")""
       write(200,"(A)")"set xrange [0.5:"//str(Nx+0.5d0)//"]"
       write(200,"(A)")"set yrange [0.5:"//str(Nx+0.5d0)//"]"
       write(200,"(A)")"set cbrange [0:2*pi]"
       write(200,"(A)")"set xtics 5"
       write(200,"(A)")"set ytics 5"
       write(200,"(A)")"set ztics 5"
       write(200,"(A)")"set cbtics ('0'0, '' pi/2, '{/Symbol-Italic p}' pi, '' 3*pi/2, '{/Symbol-Italic 2p}' 2*pi)"
       write(200,"(A)")"set palette model HSV defined ( 0 0 1 1, 1 1 1 1 )"
       write(200,"(A)")"set output '|ps2pdf -dEPSCrop - '.OUT_VECPDF"
       write(200,"(A)")"k="//str(mt_uniform(1,Nx)-1)
       write(200,"(A)")"plot PFILE every :::k::k u 1:2:6 with image, PFILE every :::k::k u 1:2:(0.6*$4):(0.6*$5) with vectors as 1"
       write(200,"(A)")""
       write(200,"(A)")"set output '|ps2pdf -dEPSCrop - '.OUT_VEC3D"
       write(200,"(A)")"set view 44,50"
       write(200,"(A)")"set xyplane at 0"
       write(200,"(A)")"splot PFILE u 1:2:3:6 with points palette pointsize 0.4 pointtype 7"
       write(200,"(A)")""
       write(200,"(A)")"set terminal gif size 450,450 nocrop animate delay 50 enhanced font 'Times-Roman'" 
       write(200,"(A)")"set output OUT_VECGIF"
       write(200,"(A)")"do for [k=0:"//str(Nx-1)//":1] {"
       write(200,"(A)")"set title 'k='.k"
       write(200,"(A)")"plot PFILE every :::k::k u 1:2:6 with image, PFILE every :::k::k u 1:2:(0.6*$4):(0.6*$5) with vectors as 1"
       write(200,"(A)")"}"
       write(200,"(A)")""
       write(200,"(A)")"system('gzip -v '.PFILE)"
       close(200)
    endif
    !
  end subroutine print_lattice


  subroutine print_point(point,pfile,last)
    integer,dimension(2) :: point
    character(len=*)     :: pfile
    logical              :: last
    real(8)              :: rho,theta,x,y
    integer,save         :: iter=0
    !
    if(.not.last)then
       iter = iter+1
       x = arrayX1(point(1))
       y = arrayX2(point(2))
       open(200,file=str(pfile),access='append')
       write(200,*)x,y,vo2_potential(point(1),point(2))/ceiling(Alat**2)/lambda
       close(200)
       return
    endif
    open(200,file="plot_"//str(pfile)//".gp")
    write(200,"(A)")"reset"
    write(200,"(A)")"set term wxt"
    write(200,"(A)")"#set terminal gif size 450,450 nocrop animate delay 50 enhanced font 'Times-Roman'" 
    write(200,"(A)")"#set output '"//str(pfile)//".gif'"
    write(200,"(A)")"reset"
    write(200,"(A)")"set size square"
    write(200,"(A)")"unset key"
    write(200,"(A)")"xf(r,phi) = r*cos(phi)"
    write(200,"(A)")"yf(r,phi) = r*sin(phi)"
    write(200,"(A)")""
    write(200,"(A)")""
    write(200,"(A)")"start=1"
    write(200,"(A)")"step=100"
    write(200,"(A)")""
    write(200,"(A)")"##Map plot"
    write(200,"(A)")"set pm3d map"
    write(200,"(A)")"do for [i=start:"//str(Nsweep)//":step] {"
    write(200,"(A)")"set title 'i='.i"   
    write(200,"(A)")"plot 'VO2_x1x2_data.dat' u 1:2:3 with image,'"//str(pfile)//"' every ::start::i using 1:2 w l ls 1  linewidth 0.3 lc rgb 'gray','"//str(pfile)//"' every ::i::i using 1:2 w p pointtype 7 pointsize 1.5 lc rgb 'green'"
    write(200,"(A)")"}"
    write(200,"(A)")""
    write(200,"(A)")""
    write(200,"(A)")"##Surface plot"
    write(200,"(A)")"#set xrange [-2.2000:2.2000]"
    write(200,"(A)")"#set yrange [-2.8000:2.8000]"
    write(200,"(A)")"##set zrange [-1:0]"
    write(200,"(A)")"#do for [i=start:"//str(Nsweep)//":step] {"
    write(200,"(A)")"#set title 'i='.i"   
    write(200,"(A)")"#splot 'VO2_x1x2_data.dat' u 1:2:3 with pm3d,'"//str(pfile)//"' every ::start::i using 1:2:3 w l ls 1  linewidth 0.3 lc rgb 'gray','"//str(pfile)//"' every ::i::i using 1:2:3 w p pointtype 7 pointsize 1.5 lc rgb 'green'"
    write(200,"(A)")"#}"
    write(200,"(A)")""
    write(200,"(A)")""
    close(200)
  end subroutine print_point



  subroutine print_angle(point,pfile)
    integer,dimension(2) :: point
    character(len=*)     :: pfile
    real(8)              :: rho,theta,x,y
    integer,save         :: iter=0
    !
    iter = iter+1
    x = arrayX1(point(1))
    y = arrayX2(point(2))
    rho  = sqrt(arrayX1(point(1))**2 + arrayX2(point(2))**2)
    theta= get_theta(arrayX1(point(1)),arrayX2(point(2)))
    open(200,file=pfile,access='append')
    write(200,*)iter,theta/pi*180
    close(200)
    return
  end subroutine print_angle






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


end program





