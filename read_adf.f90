program read_adf 
  !-----------------------------------------------------------
  !
  !     VERSION 1
  !     LAST UPDATE: 28.09.2019 M.S.
  !     updated to adf2018
  !     d.t. 29/06/2020:
  !     inserted implicit none and dynamic allocation of most
  !     arrays. Calculates nuclear dipole moment, Mulliken pop. analysis
  !     and generates various quantities needed by td-pdos
  !     tested through valgind for leaks: OK.
  !
  !     PURPOSE: EVALUATES dipOLE INTEGRALS BETWEEN
  !     EXCITED STATES
  !
  !     ONLY LENGTH GAUGE AVAILABLE
  !
  !     TAPES (INPUT)
  !           TAPE21   : GENERAL DATA 
  !           TAPE15   : ELECTRIC DIPOLE STO MATRIX ELEMENTS 
  !
  !     MUST BE LINKED WITH ADF LIBRARIES FOR KF (KEYED FILES) USAGE
  !
  !-----------------------------------------------------------
  use KF
  implicit none
  real*8, allocatable :: dipmatx(:,:) , dipmaty(:,:) , dipmatz(:,:)
  real*8, allocatable :: lmatx(:,:) , lmaty(:,:) , lmatz(:,:) 
  integer, allocatable :: nsymdav(:), ialpha(:)
  real*8, allocatable :: transmag(:,:)
  real*8, allocatable :: eigvks(:) , exciten(:), tddfteig(:,:,:), tddfteigl(:,:,:)
  real*8, allocatable :: frocin(:), frocfi(:)
  real*8, allocatable :: smat(:,:), smat_tri(:), mllkn_pop(:,:), ovrl_pop(:,:),buff(:)
  integer, parameter :: lbas =  13000, lnst = 3000, lnsym = 32
  integer, parameter :: larray = (lbas*(lbas+1))/2
  real*8, parameter :: evau = 27.21139628d0
  real*8, parameter :: twth = 2.d0/3.d0
  real*8, parameter :: two = 2.d0
  real*8, parameter :: one = 1.d0
  real*8, parameter :: zero = 0.d0
  integer, parameter :: lstidat = 3000
  integer, parameter :: lntyp =30, lnnuc = 2000
  integer :: ia, ib, ieff_sym, nnuc, ntyp, lenty, lennu, lenini, lenfin
  integer :: i_test, ityp, lalo, lahi, i, ii, iii, j, jj, jjj, k, kk, idx, kocc, kvirt, nocc, nunocc, ipair, mu, nu
  integer :: iinsy, ntotmo, nmo_in, indmoi, iinst, ipun
  integer :: indmoj, koccf, kvirtf, ifisy, naos, naosx, iu21, iu15, ifist, nmo
  integer :: ntoten, itoten, nener, isym, ios, lqtch, nsym
  integer :: ndimvx, iener, ispin
  integer :: nalloccd
  real*8 :: epsin, epsfi, dipx, dipy, dipz, exce, fvalue
  real*8 :: lx,ly,lz
  real*8 :: sig
  real*8, allocatable :: dip(:,:),eigin(:),eigfi(:),vectx(:),vecty(:),lm(:,:)
  real*8, allocatable :: vectz(:),xyznuc(:,:),e0(:),lvectx(:),lvecty(:),lvectz(:)
  !real*8 :: dip ( larray , 3 )
  !real*8 :: eigin ( lbas )
  !real*8 :: eigfi ( lbas )
  !real*8 :: vectx ( lbas ) , vecty ( lbas ) , vectz ( lbas )
  !real*8 :: xyznuc(3,lnnuc), dip_nuc(3)
  real*8 :: dip_nuc(3)
  real*8, allocatable :: qtch(:)
  integer, allocatable :: npartin(:), npartfi(:)
  integer, allocatable :: insy(:), inst(:), insp(:)
  !integer :: npartin ( lbas ) , npartfi ( lbas )
  !real*8 :: e0 ( lnst )
  !integer :: insy(lnst), inst(lnst), insp(lnst)
  integer :: jsymf(lnsym)
  integer :: nfitpt(lntyp), nqptr(lntyp), nbptr(lntyp)
  character*160 :: symrep(lnsym), secname, section, ieigstr
  !character*160, allocatable :: lab(:)
  character*11  :: eigspin
  character*5   :: eivspin, base_name
  character*4 :: excnr
  character*15 :: file_name
  logical :: locc,loccf
  logical, allocatable :: lrep2do(:)
  logical :: b3lyp,cdspectrum
  write(*,97) ' ************************************************ '
  write(*,99) ' *  DIP_STEC PROGRAM - OSCILLATOR STRENGTHS     * '
  write(*,99) ' * QUANTUM CHEMISTRY GROUP - TRIESTE UNIVERSITY * '
  write(*,99) ' *         LAST UPDATE: 28 SEPTEMBER 2019       * '
  write(*,99) ' ************************************************ '

97 format (//// , 20x , a )
99 format (     20x , a )
  allocate (dip(larray, 3), eigin(lbas), eigfi(lbas), vectx(lbas), &
           vecty (lbas) , vectz(lbas), xyznuc(3,lnnuc), e0(lnst))
  allocate (npartin(lbas), npartfi(lbas), insy(lnst), inst(lnst), insp(lnst))

  write(*,*) 'B3LYP or TDA calculation (.true. or .false.)'
  read(*,*) b3lyp
  write(*,*) 'ADF calculation with CD spectrum'
  read(*,*) cdspectrum

  if (cdspectrum) then
     allocate (lvectx(lbas),lvecty(lbas),lvectz(lbas))
     allocate(lm(larray,3))
  endif 

  !INPUT AND CHECK SECTION
  !START READ TAPES AND CHECK SECTION
  !open the kf-file fname and return the kf-unit number
  !the file should already exist
  call KFOPFL  (iu21, 'TAPE21')
  !Makes an existing variable the current variable
  call KFOPVR  (iu21, 'Symmetry%nsym')
  !Read the variable
  call KFREAD  (iu21, 'Symmetry%nsym', nsym)
  call KFOPVR  (iu21, 'Symmetry%symlab')
  !read character symrep(nsym)  
  call KFRDNS  (iu21, 'Symmetry%symlab', symrep , nsym , 1)
  call KFOPFL  (iu15, 'TAPE15')
  call KFREAD  (iu15, 'Basis%naos', naos)
  naosx=(naos*(naos+1))/2
  if (naosx.gt.larray) then
     write(*,*) 'naosx larray', naosx, larray
     stop ' CHANGE parameter lbas: naosx.GT.larray '
  end if
  !check dimensions
  lenty = KFLEN ( iu21 , 'Basis%nbaspt' )
  if ( lenty .gt. lntyp ) then
     stop ' lenty .GT. lntyp '
  end if
  lennu = KFLEN ( iu21 , 'Geometry%xyz' )
  if ( lennu .gt. ( lnnuc * 3 ) ) then
     stop ' lennu .GT. ( lnnuc * 3 ) '
  end if

  lqtch = KFLEN ( iu21 , 'Geometry%qtch' )
  !allocate (qtch(lqtch), lab(lqtch))
  allocate (qtch(lqtch))
  write(*,*) ' lqtch = ',lqtch

  call KFRDNR ( iu21 , 'Geometry%xyz' , xyznuc , lennu , 1 )
  call KFREAD ( iu21 ,         'nnuc' , nnuc   )
  call KFRDNI ( iu21,'nqptr', nqptr, lenty, 1 )
  call KFREAD ( iu21,'ntyp', ntyp )
  call KFRDNR ( iu21 ,         'qtch' , qtch   , lqtch , 1 )
  !call KFRDNS ( iu21 ,       'nuclab' , lab    , lqtch , 1 )
  call KFRDNI ( iu21 , 'Basis%nbptr', nbptr, lenty, 1 )
  !nuc. contribution to dipole moment
  do ityp = 1, ntyp
     write(*,*) ' ityp: ', ityp, ' qtch(ityp): ', qtch(ityp)  
     do j = nqptr(ityp), nqptr(ityp+1) - 1
        !write(*,*) ' j: ', j, ' X: ', xyznuc(1,j), &
        !   & ' Y: ', xyznuc(2,j), ' Z: ', xyznuc(3,j)
        dip_nuc(1) = dip_nuc(1) + qtch(ityp) * xyznuc(1,j)
        dip_nuc(2) = dip_nuc(2) + qtch(ityp) * xyznuc(2,j)
        dip_nuc(3) = dip_nuc(3) + qtch(ityp) * xyznuc(3,j)
      end do
  end do 

  write(*,*) ' dip_nuc: ', dip_nuc
  open(10,file='dip_nuc.dat')
  write(10,*) ' nuclear dipole contribution '
  write(10,*) dip_nuc
  close(10)
  !costruct overlap matrix over AOs
  write(*,*) ' naos, naosx = ' , naos, naosx
  allocate(smat_tri(naosx), smat(naos,naos))
  call KFRDNR (iu15 ,'Matrices%Smat', smat_tri, naosx, 1 )
  !write(*,*) ' smat_tri: ', smat_tri
  k = 0
  do j = 1, naos
     do i = 1, j
        k = k + 1
        smat(i,j) = smat_tri(k)
        smat(j,i) = smat_tri(k)
     enddo
  enddo
  
  deallocate (smat_tri)
  !read real dip(naosx,3)
  call KFRDNR  (iu15, 'Matrices%Dipmat_x', dip(1,1), naosx , 1)
  call KFRDNR  (iu15, 'Matrices%Dipmat_y', dip(1,2), naosx , 1)
  call KFRDNR  (iu15, 'Matrices%Dipmat_z', dip(1,3), naosx , 1)
  if (cdspectrum) then
  !read real lmat(naosx,3) for circular dichroism
     call KFRDNR  (iu15, 'Matrices%Lmat_x', lm(1,1), naosx , 1)
     call KFRDNR  (iu15, 'Matrices%Lmat_y', lm(1,2), naosx , 1)
     call KFRDNR  (iu15, 'Matrices%Lmat_z', lm(1,3), naosx , 1)
  endif

  !start loop over initial states
  eigspin = 'Eigen-Bas_A'
  eivspin = 'eps_A'
  ntotmo=0
  kocc =  0
  kvirt = 0
  do iinsy = 1 , nsym
     !open a section on the file iu and optionally create it
     call KFOPSC  (iu21, symrep (  iinsy ) )
     call KFREAD  (iu21, 'nmo_A', nmo_in)
     !close file
     call KFCLSC  (iu21)
     ntotmo=ntotmo+nmo_in
  enddo
  allocate(dipmatx(ntotmo,ntotmo))
  allocate(dipmaty(ntotmo,ntotmo))
  allocate(dipmatz(ntotmo,ntotmo))
  if (cdspectrum) then
     allocate(lmatx(ntotmo,ntotmo))
     allocate(lmaty(ntotmo,ntotmo))
     allocate(lmatz(ntotmo,ntotmo))
  endif
  allocate(eigvks(ntotmo))
  dipmatx = zero
  dipmaty = zero
  dipmatz = zero
  if (cdspectrum) then
     lmatx = zero
     lmaty = zero
     lmatz = zero
  endif
  open(50,file='mos_info.dat')
  write(50,*) ntotmo
  write(50,*) nsym
  write(*,*) 'ntotmo, nsym: ',  ntotmo, nsym
  !lookhere
  allocate(mllkn_pop(nnuc,ntotmo), ialpha(naos),&
          &buff(naos),ovrl_pop(naosx-naos,ntotmo))
  mllkn_pop = zero
  ialpha = 0
  idx = 0
  !run over atom type
  !naos -> total number of basis functions
  !nbos -> number of basis functions excluding repetitions of atoms of the same
  !        type
  !ntyp -> number of atom types
  !nbptr(ntyp) -> number of basis functions (nbos) per atom type [cumulative
  !               index]
  !nqptr(ntyp) -> number of atoms per atom type [cumulative index]
  !ialpha -> AO function idx belongs to atom j
  do ityp = 1, ntyp
     lalo = nbptr(ityp)
     lahi = nbptr(ityp+1) - 1
     !run over the atoms belonging to the same type
     do j = nqptr(ityp), nqptr(ityp+1) - 1
        do ia = lalo, lahi
           idx = idx + 1
           ialpha(idx)= j
         end do
     end do
  end do
  open(40,file='mulliken_pop.dat')
  write(40,*) nnuc, ntotmo
  !lookhere
  open(41,file='ovrl_pop.dat')
  write(41,*) naosx-naos, ntotmo
  open(42,file='AO_map.dat')
  write(42,*) naos
  write(42,*) '         nAO  ---  atom'
  !run over symmetries
  indmoi=0
  do iinsy = 1 , nsym
     call KFOPSC  (iu21, symrep (  iinsy ) ) 
     call KFREAD  (iu21, 'nmo_A', nmo_in)
     !length of the variable npart
     lenini = KFLEN ( iu21, 'npart')
     !read integer
     !npartin(lenini)
     call KFRDNI  (iu21, 'npart', npartin, lenini, 1 )
     allocate (frocin(nmo_in))
     !frocin(nmo_in)
     call KFRDNR  (iu21, 'froc_A', frocin, nmo_in, 1 )
     do i = 1, nmo_in
        if (frocin(i).gt.1.0d0) then
           nocc=i
        else
           exit
        endif
     enddo
     nunocc = nmo_in - nocc
     write(50,*) iinsy, nocc, nunocc
     call KFCLSC  (iu21) 
     !run over states
     do iinst = 1, nmo_in
        indmoi=indmoi+1
        if (frocin(iinst).gt.1.0d0) then
           locc = .true.
           kocc = kocc + 1
        else
           locc = .false.
           kvirt = kvirt + 1
        endif
        call KFOPSC  (iu21, symrep (  iinsy ) ) 
        ipun =  (iinst - 1 ) * lenini
        call KFOPVR  (iu21, eigspin)
        !Skips n elements of the current variable. The  variable
        !should exist and all skipped elements should have
        !been written previously
        call KFSKVR  (iu21, ipun)
        call KFRDNR  (iu21, '%', eigin, lenini , 1)
        call KFOPVR  (iu21, eivspin)
        call KFSKVR  (iu21, iinst-1 )
        call KFREAD  (iu21, '%', epsin)
        call KFCLSC  (iu21) 
        eigvks(indmoi)=epsin
        !calculate Mulliken population
        do i = 1, naos
           buff(i) = zero
           do k = 1, lenini
              kk = npartin(k)
              buff(i) = buff(i) + smat(kk,i)*eigin(k)
           end do
        end do
        do i = 1, lenini
           ii = npartin(i)
           iii = ialpha(ii)
           mllkn_pop(iii,indmoi)=mllkn_pop(iii,indmoi)+buff(ii)*eigin(i)
        end do
        !lookhere
        !here we can also define the OPDOS:
        !N_\mu,\nu(E) = 2*\sum_i C_\mu,iC_\nu,i * S_\mu,\nu
        !where indices \mu,\nu run over the basis functions
        !ovrl_pop(1:naosx,indmoi)=.....
        ipair=0
        do mu=1,naos
          do nu=mu+1,naos
            ipair=ipair+1
            ovrl_pop(ipair,indmoi)=2*eigin(mu)*eigin(nu)*smat(mu,nu)
          enddo
        enddo

10      format(//,' INITIAL SUBSPECIES: ',A15,' INITIAL ORBITAL: ',I5)
11      format(' INITIAL KS EIGENVALUE (eV) = ',F15.5)
12      format(' SPIN = ',I5,6X,' OCCUPATION = ', F15.8,/)
13      format(1X,A)
14      format(1X,A,3X,A,5X,A,7X,A,12X,A,12X,A,9X,A,/)
        !start loop over final symmetries
        indmoj=0
        koccf=0
        kvirtf=0
        do ifisy = 1, nsym
           call KFOPSC  (iu21, symrep ( ifisy ) )
           call KFREAD  (iu21, 'nmo_A', NMO)
           lenfin = KFLEN ( iu21, 'npart')
           call KFRDNI  (iu21, 'npart', npartfi, lenfin, 1 )
           allocate (frocfi(nmo))
           call KFRDNR  (iu21, 'froc_A', frocfi, nmo, 1 )
           call KFCLSC  (iu21)
           do j = 1 , lenfin
              vectx(j) = 0.0
              vecty(j) = 0.0
              vectz(j) = 0.0
              if (cdspectrum) then
                 lvectx(j) = 0.0
                 lvecty(j) = 0.0
                 lvectz(j) = 0.0
              endif
              jj = npartfi ( j )
              do i = 1 , lenini
                 ii = npartin ( i )
                 iii = MAX ( ii, jj )
                 jjj = MIN ( ii, jj )
                 sig = one
                 if (ii.lt.jj) then
                    sig = - one
                 endif
                 k = ( iii * ( iii - 1 ) ) / 2 + jjj
                 vectx(j) = vectx(j) +  eigin ( i ) * dip ( k , 1 )
                 vecty(j) = vecty(j) +  eigin ( i ) * dip ( k , 2 )
                 vectz(j) = vectz(j) +  eigin ( i ) * dip ( k , 3 )
                 if (cdspectrum) then
                    lvectx(j) = lvectx(j) +  eigin ( i ) * lm ( k , 1 ) * sig
                    lvecty(j) = lvecty(j) +  eigin ( i ) * lm ( k , 2 ) * sig
                    lvectz(j) = lvectz(j) +  eigin ( i ) * lm ( k , 3 ) * sig
                 endif
              end do
           end do

           !start loop over final states
           do ifist = 1, nmo
              indmoj=indmoj+1
              if (indmoj.gt.indmoi) cycle
              loccf=.false.
              if (frocfi(ifist).gt.1.0d0) then
                 koccf = koccf + 1
                 loccf=.true.
              else
                 kvirtf = kvirtf + 1
              endif
              ipun =  (ifist - 1 ) * lenfin
              call KFOPSC  (iu21, symrep ( ifisy ) )
              call KFOPVR  (iu21, eigspin)
              call KFSKVR  (iu21, ipun)
              call KFRDNR  (iu21, '%', eigfi, lenfin , 1)
              call KFOPVR  (iu21, eivspin)
              call KFSKVR  (iu21, ifist-1 )
              call KFREAD  (iu21, '%', epsfi)
              call KFCLSC  (iu21)
              !calculate dipole integrals < IN | D | FI >
              dipx = .0
              dipy = .0
              dipz = .0
              if (cdspectrum) then
                 lx = .0
                 ly = .0
                 lz = .0
              endif
              do j = 1, lenfin
                 dipx = dipx + eigfi(j)*vectx(j)
                 dipy = dipy + eigfi(j)*vecty(j)
                 dipz = dipz + eigfi(j)*vectz(j)
                 if (cdspectrum) then
                    lx = lx + eigfi(j)*lvectx(j)
                    ly = ly + eigfi(j)*lvecty(j)
                    lz = lz + eigfi(j)*lvectz(j)
                 endif 
              end do
              dipmatx(indmoi,indmoj) = dipx
              dipmatx(indmoj,indmoi) = dipx
              dipmaty(indmoi,indmoj) = dipy
              dipmaty(indmoj,indmoi) = dipy
              dipmatz(indmoi,indmoj) = dipz
              dipmatz(indmoj,indmoi) = dipz

              if (cdspectrum) then
                lmatx(indmoi,indmoj) = lx
                lmatx(indmoj,indmoi) = -lx
                lmaty(indmoi,indmoj) = ly
                lmaty(indmoj,indmoi) = -ly
                lmatz(indmoi,indmoj) = lz
                lmatz(indmoj,indmoi) = -lz
              endif

              exce = ( epsfi - epsin ) * evau
              fvalue = (dipx ** 2 + dipy ** 2 + dipz ** 2)
              fvalue = fvalue * ( epsfi - epsin ) * two * twth
15            format(2X,A12,X,I5,X,F12.6,4E14.6)
17            format(2x,f12.6,e14.6)
           end do !loop over final states
           deallocate ( frocfi )
        end do !loop over final symmetries
     end do !loop over initial states
     deallocate ( frocin )
  end do  !loop over initial symmetries
  do i = 1, ntotmo
     write(50,*) eigvks(i)
  end do
  close(50)  
  write(*,*) ' Mulliken population...'
  do i = 1, ntotmo
     write(40,*) (mllkn_pop(j,i), j=1, nnuc)
     !write(*,*) (mllkn_pop(j,i), j=1, nnuc)
  end do
  !lookhere
  write(*,*) ' Overlap population...'
  do i = 1, ntotmo
     write(41,*) (ovrl_pop(j,i), j=1, naosx-naos)
     !write(*,*) (ovrl_pop(j,i), j=1, naosx-naos)
  end do
  do i = 1, naos
       write(42,*) i, ialpha(i)
  end do
  close (40)
  close (41)
  close (42)
  !deallocate(mllkn_pop, buff, ialpha)
  nsym = 0
  call KFOPVR  (iu21, 'Symmetry%nsym excitations')
  call KFREAD  (iu21, 'Symmetry%nsym excitations', nsym)
  write(*,*) ' nsym excitations = ' , nsym
  call KFOPVR  (iu21, 'symlab excitations')
  call KFRDNS  (iu21, 'symlab excitations', symrep, nsym, 1)
  do i = 1, nsym
     write(*,*) ' isy symrep ', symrep(i)
  enddo
  call KFOPVR  (iu21, 'vecdimension excitations')
  call KFREAD  (iu21, 'vecdimension excitations', ndimvx)
  allocate(nsymdav(nsym))
  allocate(lrep2do(nsym))
  call KFOPVR  (iu21, 'nsymdav excitations')
  call KFRDNI  (iu21, 'nsymdav excitations', nsymdav, nsym, 1)
  call KFOPVR  (iu21, 'lrep2do excitations') 
  !read logical
  call kfrdnl (iu21, 'lrep2do excitations', lrep2do, nsym, 1) 
  write(*,*) '  nsymdav = ', nsymdav
  write(*,*) '  lrep2do  = ', lrep2do 
  section = 'SS '
  ntoten = 0
  do isym = 1 , nsym
     if (.not.lrep2do(isym)) cycle
     write (secname,'(A12,A2,1x,A)') 'Excitations ', trim(section), trim(symrep(isym))
     secname = trim(secname)
     call KFOPSC  (iu21, secname)
     call KFREAD  (iu21, 'nr of excenergies', nener )
     ntoten = ntoten + nener
  end do
  write(*,*) ' ntoten ' , ntoten
  allocate(tddfteig(kvirt, kocc, ntoten), tddfteigl(kvirt, kocc, ntoten))
  allocate(exciten(ntoten))
  open(80,file='eig.dat')
  open(81,file='eig_l.dat')
  open(60,file='tddft_info.dat')
  ieff_sym=0
  do isym =1,nsym
     if (.not.lrep2do(isym)) cycle
     ieff_sym=ieff_sym+1
  end do
  write(80,*) ieff_sym,kocc,kvirt,ntoten,ntotmo,cdspectrum
!  write(81,*) ieff_sym,kocc,kvirt,ntoten,ntotmo,cdspectrum
  write(60,*) ntoten
  close(60)
  itoten = 0
  base_name = 'o-cfc'
  do isym = 1 , nsym
     if (.not.lrep2do(isym)) cycle 
     !write(*,*) symrep(isym)
     !write(*,*) trim(symrep(isym) )
     write (secname,'(A12,A2,1x,A)') 'Excitations ',trim(section),trim(symrep(isym)) 
     !write (secname,'(A12,A2,1x,A)') 'Excitations ', trim(section)
     !write(*,*) secname
     secname = trim(secname)
     write(*,*) secname
     call KFOPSC  (iu21, secname)
     call KFREAD  (iu21, 'nr of excenergies', nener )
     write(*,*) ' symm nener', secname,nener
     call KFRDNR(iu21,'excenergies',exciten(itoten+1),nener,1)
     write(80,*) nener
!    write(81,*) nener
     do iener= 1, nener
        itoten = itoten + 1
        write(*,*) ' itoten: ', itoten
        ieigstr = ' '
        call csputi(ieigstr, iener)
        !v. 2018
        call KFRDNR(iu21,'eigenveceps '//trim(ieigstr),tddfteig(1,1,itoten),ndimvx,1)
        call KFRDNR(iu21,'eigenveceps_mag '//trim(ieigstr),tddfteigl(1,1,itoten),ndimvx,1)
        !write(*,*) 'eigenveceps_mag '// trim(ieigstr)
        !call KFRDNR(iu21,'eigenvector '//trim(ieigstr),tddfteig(1,1,itoten),ndimvx,1)
        !write(*,*) 'eigenvector '// trim(ieigstr)
        write(excnr,'(i4.4)') itoten
        write(file_name,*) adjustl(trim(base_name)), adjustl(trim(excnr)),".dat"
        write(*,*) ' file_name: ', file_name     
        open(40, file=file_name, form='formatted ', iostat=ios)
        write(40,*) ' #from ADF'
        write(40,*) ' #starting   #final  #ispin   Cij '
        ispin = 1 
        do i=1,kocc
           do j=1,kvirt
              if (.not.b3lyp) then
                 tddfteig(j,i,itoten)=tddfteig(j,i,itoten)/dsqrt(exciten(itoten))
                 tddfteigl(j,i,itoten)=tddfteigl(j,i,itoten)*dsqrt(exciten(itoten))
              endif
              !tddfteig(j,i,itoten)=tddfteig(j,i,itoten)*dsqrt((eigvks(nocc+j)-eigvks(i))/exciten(itoten)) 
              write(80,*) tddfteig(j,i,itoten)
              write(81,*) tddfteigl(j,i,itoten)
              !write(*,*) ' j, i, tddfteig: ', j, i, tddfteig(j,i,itoten)
              write(40, "(1x, i4, 2x, i4, 2x, i4, 2x, e21.13)") i, kocc+j, ispin, &
                   tddfteig(j,i,itoten)
           end do
        end do
        close(40)
        !debug
        !write(*,*) tddfteig
     enddo
  enddo
  write(*,*) ' kocc kvirt ' , kocc, kvirt, nocc
  write(*,*) ' excit ener ',exciten
  close(80)
  close(81)
  open(70,file='ene.dat')
  do i=1,ntoten
     write(70,*) exciten(i)
  enddo
  close(70)
  open(90,file='dipmat.dat')
  !add nuc. contribution
  !dipmatx = dipmatx
  !dipmaty = dipmaty
  !dipmatz = dipmatz
  do i=1,ntotmo
     do j=1,ntotmo
        write(90,*) dipmatx(j,i),dipmaty(j,i),dipmatz(j,i)
     enddo
  enddo      
  close(90)
  if (cdspectrum) then
     open(91,file='lmat.dat')
     do i=1,ntotmo
        do j=1,ntotmo
           write(91,*) lmatx(j,i),lmaty(j,i),lmatz(j,i)
        enddo
     enddo
     close(91)
  endif

  deallocate (dip,eigin,eigfi,vectx,vecty,vectz,xyznuc,e0)
  deallocate (npartin, npartfi, insy, inst, insp)
  deallocate(tddfteig, tddfteigl)
  deallocate(dipmatx)
  deallocate(dipmaty)
  deallocate(dipmatz)
  if (cdspectrum) then
     deallocate(lmatx)
     deallocate(lmaty)
     deallocate(lmatz)
     deallocate(lvectx,lvecty,lvectz)
     deallocate(lm)
  endif
  deallocate(eigvks)
  deallocate(nsymdav)
  deallocate(lrep2do)
  deallocate(qtch)
  deallocate(mllkn_pop, buff, ialpha)
  write(*,98) ' ************************************ '
  write(*,96) ' * REGULAR TERMIMATION OF READ_ADF  * '
  write(*,95) ' ************************************ '
98 format(//// , 30X , A     )
96 format(     30X , A     )
95 format(     30X , A ,//// )
  write(*,*) 'TERMINATION'
  stop
end program read_adf 
