program cippres_stieltjes
  use general
  use interpolation
  implicit none
  BEGIN_DOC
! CIPPRES stands for Configuration Interaction Plugin for Photoionized and Resonant Electronic States.
! ORMAS-CI calculations are performed to compute bound and approximate continuum states. Hamiltonian and Dipole coupling matrix elements can then be computed between these states. Finally, Stieltjes imaging technique is applied to recover the correct continuum coupling matrix elements.
! cippres_stieltjes computes the decay widths using Stieltjes imaging and the Fano couplings obtained from cippres_fano
  END_DOC

  integer :: npt
  integer, parameter :: QR_K = selected_real_kind (32)
  real (kind=QR_K), allocatable, dimension(:) :: e, g

  integer :: nmax = 30 ! according to Mueller-Plathe and Diercksen Stieltjes is inaccurate for n>=15
!real (kind=QR_K), dimension(0:nmax) :: sk
!real (kind=QR_K), dimension(nmax,nmax) ::  e1, g1, e2, g2
  real (kind=QR_K), allocatable, dimension(:) :: sk
  real (kind=QR_K), allocatable, dimension(:,:) ::  e1, g1
  real (kind=QR_K) :: shift1 = 0.1d0
  integer :: imax1, imax2
  integer :: inmax, ishift
  integer :: i, j
  integer :: exit_cycle
  character(len=60) :: fname

! Tsveta
  real (kind=QR_K) :: g_
  real (kind=QR_K) :: temp

! TODO Read couplings in EZFIO
! TODO Adapt Stieltjes code
! TODO Print results in EZFIO (some plots would be nice)

! GENERAL
! TODO Compute dipole matrix elements between two different CI runs
! TODO Include Stieltjes in qp

  if(ifcsf==3) then

    if(ifanosta==0) then
       print*, "Please set ifanosta"
       print*, "qp set cippres ifanosta X "
       stop
    endif

! reads and sorts energy and matrix elements

  npt = n_csf_cippres(ici2)
  allocate(e(npt),g(npt))
  e(:) = efano_cippres(1:n_csf_cippres(ici2),ifanosta)
  g(:) = cfano_cippres(1:n_csf_cippres(ici2),ifanosta)
 
  if(sum(g)==0d0) then
    write(*,*)"All matrix elements are zero, I stop"
    stop
  endif

  call sort2(npt,e,g)
  ishift = 2
  shift1 = ishift * 0.1d0
  write(fname, '(A8,I1,A4)')"gamma.sh",ishift,".dat"
  open(237,file=fname)

  do inmax = 15, nmax
     allocate(sk(0:inmax))
     allocate(e1(inmax,inmax), g1(inmax,inmax))       
     
! computes the nmax moments of the matrix elements
      open(unit=10,file='moments.txt')
      sk(:) = 0q0
      do i = 0, inmax
         do j = 1, npt
            sk(i) = sk(i) + e(j)**i*g(j)
         enddo
       write(10,*)i,abs(sk(i))
      enddo
      close(10)
    
! "images" the matrix elements for two different shift values, may be used to evaluate the highest accurate order
      shift1 = shift1 + abs(e(1))
      e(:) = e(:) + shift1
      call imaging(npt,e,g,inmax,imax1,e1,g1)
      e(:) = e(:) - shift1
       
      do i = 1, imax1
       do  j = 1, i-1
        write(100+i,'(2(f20.16,1X))')e1(j,i)-shift1,g1(j,i)
       enddo
      enddo
       
      call interp(e1(:,imax1),g1(:,imax1),imax1-1,shift1,g_)  
      g_ = 2.0d0*pi*g_
      shift1 = shift1 - abs(e(1))
      !write(*,'(1X,A9,f8.3)')"Shift 1: ",shift1-abs(e(1))
      !write(*,'(I2,A3,E23.15)')inmax," ",g_
     
      ! March 03, 2016 Tsveta
      write(237,'(I3,A3,E23.15)')inmax," ",g_
      
      deallocate(sk,e1,g1)
  end do
  close(237)
  deallocate(e,g)
    
  call ezfio_set_cippres_ifcsf(4)

  else 

    print*, "Please run cippres_fano first"

  endif

end program cippres_stieltjes


subroutine imaging(npt,e_point,g_point,nmax,maxord,eorder,gorder)
implicit none

 integer :: nmax

 integer, parameter :: QR_K = selected_real_kind (32) ! use quadruple precision
 double precision, parameter :: overmax=1d0 ! set arbitrary to 1 to determine the maximum order of pol.

 integer, intent(in) :: npt
 real (kind=QR_K), dimension(npt), intent(in) :: e_point, g_point

 real (kind=QR_K), dimension(0:nmax,npt) :: qpol
 real (kind=QR_K), dimension(nmax) :: acoef
 real (kind=QR_K), dimension(0:nmax) :: bcoef

 real (kind=QR_K), dimension(nmax) :: diag, offdiag
 real (kind=QR_K), dimension(nmax,nmax) :: abvec
 real (kind=QR_K) :: asum, bprod, qnorm, qoverlap
 integer :: iord, ierr, maxord, min, max

 real (kind=QR_K), dimension(nmax) :: enew, gnew
 real (kind=QR_K), dimension(nmax,nmax) :: eorder, gorder

 integer :: i, j

! initiate the recursive computation of the a,b coefficients and the orthogonal 
! polynomials according to (3.3.20-23) of Mueller-Plathe & Dierksen (1990)
       bcoef(0)=0.q0
       acoef(1)=0.q0
       do i=1,npt
          bcoef(0)=bcoef(0)+g_point(i)
          acoef(1)=acoef(1)+g_point(i)/e_point(i)
       end do
       acoef(1)=acoef(1)/bcoef(0)

       do i=1,npt
          qpol(0,i)=1.q0
          qpol(1,i)=1.q0/e_point(i)-acoef(1)
       end do

       bcoef(1)=0.q0
       acoef(2)=0.q0
       do i=1,npt
          bcoef(1)=bcoef(1)+qpol(1,i)*g_point(i)/e_point(i)
          acoef(2)=acoef(2)+qpol(1,i)*g_point(i)/(e_point(i)**2)
       end do
       bcoef(1)= bcoef(1)/bcoef(0)
       acoef(2)=acoef(2)/(bcoef(0)*bcoef(1))-acoef(1)

! calculate the higher-order coefficients and polynomials recursively
! up to the (NMAX-1)th order (total of NMAX polynomials)

       asum=acoef(1)
       do i=3,nmax

          asum=asum+acoef(i-1)

          do j=1,npt
             qpol(i-1,j)=(1.q0/e_point(j)-acoef(i-1))*qpol(i-2,j)-bcoef(i-2)*qpol(i-3,j)
          end do

          bprod=bcoef(0)
          do j=1,i-2
             bprod=bprod*bcoef(j)
          end do

          bcoef(i-1)=0.q0
          do j=1,npt
             bcoef(i-1)=bcoef(i-1)+qpol(i-1,j)*g_point(j)/(e_point(j)**(i-1))
          end do
          bcoef(i-1)=bcoef(i-1)/bprod

          bprod=bprod*bcoef(i-1)

          acoef(i)=0.q0
          do j=1,npt
             acoef(i)=acoef(i)+qpol(i-1,j)*g_point(j)/(e_point(j)**i)
          end do
          acoef(i)=acoef(i)/bprod-asum

       end do

! calculate the nmax-th order polynomial just for the orthogonality check 
       do j=1,npt
          qpol(nmax,j)=(1.q0/e_point(j)-acoef(nmax))*qpol(nmax-1,j)-bcoef(nmax-1)*qpol(nmax-2,j)
       end do

! check the orthogonality of the polynomials to define the maximal approximation order 
! if the orthogonality is preserved for all orders, MAXORD is set to NMAX
       maxord=nmax
       qnorm=bcoef(0)
       do i=1,nmax
          qnorm=0.q0
          qoverlap=0.q0
          do j=1,npt
             qnorm=qnorm+qpol(i,j)**2*g_point(j)
             qoverlap=qoverlap+qpol(i,j)*qpol(i-1,j)*g_point(j)
          end do
!nico change for qp2          if (qabs(qoverlap).lt.1.q-50) qoverlap=1.q-50
          if (abs(qoverlap).lt.1.q-50) qoverlap=1.q-50
          
!nico change for qp2          if (qnorm/qabs(qoverlap).le.overmax) then
          if (qnorm/abs(qoverlap).le.overmax) then
! MAXORD=I-1 is appropriate since the polynomial failing 
! the orthogonality check should not be used
             maxord=i-1
             go to 10
          end if
       end do
 10    continue

! look how many Stieltjes orders are available
       if (maxord.lt.5) then
          min=maxord
          max=maxord
          print*, '***WARNING*** Stieltjes:'
          print*, ' only very low-order approximation is available'
          print*, ' MAXORD=',maxord
       else
          min=5
          max=maxord
          print*, ' MAXORD=',maxord
       end if

! perform the gamma calculation using the successive approximations 
! n=5,...,nmax

   do iord=5,maxord

     write(*,*)"Performs Stieltjes at order",iord 

! fill the coefficients matrix
       do i=1,iord
          diag(i)=acoef(i)
!          write(*,*)"diag",i,diag(i)
       end do
       do i=2,iord
!nico change for qp2          offdiag(i)=-qsqrt(bcoef(i-1))
          offdiag(i)=-sqrt(bcoef(i-1))
!          write(*,*)"offdiag",i,offdiag(i)
       end do

! diagonalize the coefficients matrix
! initialize the arrays
       do i=1,nmax
          do j=1,nmax
             abvec(i,j)=0.q0
          end do
          abvec(i,i)=1.q0
       end do
       call tql2(nmax,iord,diag,offdiag,abvec,ierr)
       if (ierr.ne.0) then
          print*, '***WARNING*** Stieltjes:'
          print*, ' the eigenvalue no. ',ierr,' failed to converge'
       end if

! fill the Stieltjes energy and gamma arrays
! note that the eigenvalues are inverse energies and are given in ascending order 
       do i=1,iord
          enew(i)=1.q0/diag(iord+1-i)
          gnew(i)=bcoef(0)*abvec(1,iord+1-i)**2
       end do

       call eigsrtnico(enew,gnew,iord,nmax)

! calculate the gamma's by simple numerical differentiation at the middle 
! point of each [ENEW(I),ENEW(I+1)] interval
       do i=1,iord-1
          eorder(i,iord)=0.5d0*(enew(i)+enew(i+1))
          gorder(i,iord)=0.5d0*(gnew(i+1)+gnew(i))/(enew(i+1)-enew(i))
       end do

   enddo ! loop over iord  

end subroutine imaging

