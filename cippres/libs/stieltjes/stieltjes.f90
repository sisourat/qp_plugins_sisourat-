program stieltjes
use interpolation
implicit none

integer :: npt
integer, parameter :: QR_K = 16 !selected_real_kind (32)
real (kind=QR_K), allocatable, dimension(:) :: e, g

integer :: nmin=11, nmax = 26 ! according to Mueller-Plathe and Diercksen Stieltjes is inaccurate for n>=15
!real (kind=QR_K), dimension(0:nmax) :: sk
!real (kind=QR_K), dimension(nmax,nmax) ::  e1, g1, e2, g2
real (kind=QR_K), allocatable, dimension(:) :: sk, gord
real (kind=QR_K), allocatable, dimension(:,:) ::  e1, g1
real (kind=QR_K), allocatable, dimension(:) :: eallord,gallord
real (kind=QR_K) :: shift1 = 0.1d0
integer :: imax1, imax2
integer :: inmax, ishift
integer :: i, j, k, ichan
integer :: exit_cycle
character(len=60) :: fname

! Tsveta
real (kind=QR_K) :: g_
real (kind=QR_K) :: temp, pi, gav, stadev, gtot

integer :: nsta
real (kind=QR_K), allocatable, dimension(:,:) :: gpart
real (kind=QR_K), allocatable, dimension(:) :: st_gpart

pi = dacos(-1d0)

! reads and sorts energy and matrix elements

read(*,*)npt!, nsta!, nmax
nsta=0
allocate(e(npt),g(npt),gpart(npt,nsta))
allocate(eallord(nmax**2),gallord(nmax**2))
allocate(st_gpart(nsta))

do i = 1, npt
  read(*,*)e(i),g(i),(gpart(i,j),j=1,nsta)
enddo
 g(:) = abs(g(:))
 gpart(:,:) = abs(gpart(:,:))

!do i = 1, npt
!g(i) = sum(gpart(i,:))
!enddo
 
if(sum(g)==0d0) then
 write(*,*)"All matrix elements are zero, I stop"
 stop
endif

call sort3(npt,nsta,e,g,gpart)
open(unit=20,file='sorted.txt')
 do j = 1, npt
   write(20,'(2000(f23.15,1X))')e(j),g(j),(gpart(j,k),k=1,nsta)
 enddo
close(20)


ishift = 2
shift1 = ishift * 0.1d0
write(fname, '(A8,I1,A4)')"gamma.sh",ishift,".dat"

!! TOTAL RATES

open(237,file=fname)

 allocate(gord(0:nmax))
 allocate(sk(0:nmax))
 allocate(e1(nmax,nmax), g1(nmax,nmax))       
 gord = 0d0

 shift1 = shift1 + abs(e(1))
 e(:) = e(:) + shift1
 call imaging(npt,e,g,nmax,nmax,e1,g1)
 e(:) = e(:) - shift1

 do i = nmin, nmax
   do  j = nmin, i-1
      write(100+i,'(2(f20.16,1X))')e1(j,i)-shift1,g1(j,i)*2.0d0*pi*27211
   enddo
 enddo
!    gav = sum(gord(nmax-9:nmax))/10
       
 do imax1 = nmin, nmax
   write(*,*)"Interp at order", imax1
   call interp(e1(:,imax1),g1(:,imax1),imax1-1,shift1,g_)  
   g_ = 2.0d0*pi*g_
   gord(imax1) = g_
  write(237,'(I3,A3,F23.15,A)')imax1," ",g_*27211,' in meV'
 end do

 k = 0
 do i = nmin, nmax
   do  j = 1, i-1
      k = k + 1
      eallord(k) = e1(j,i)-shift1
      gallord(k) = g1(j,i)*2.0d0*pi
   enddo
 enddo
!    gav = sum(gord(nmax-9:nmax))/10
!    stadev = sqrt(sum( (gord(nmax-9:nmax)-gav)**2) )/10
 gav = sum(gord(nmin:nmax))/6
 stadev = sqrt(sum( (gord(nmin:nmax)-gav)**2) )/(nmax-nmin)

 write(*,'(2(f20.12,1X),a)')gav,stadev, 'in au'
 write(*,'(2(f20.12,1X),a)')gav*27210,stadev*27210, 'in meV'
 gtot = gav*27210

call sort2(k,eallord,gallord)

do i = 1, k
 write(222,'(2(f20.12,1X))')eallord(i),gallord(i)*27211
enddo
   call interp(eallord,gallord,k-1,0q0,g_)  
!   write(*,*)
   write(*,'(a)')"TOTAL RATE"
   write(*,'(1(f20.12,1X),a)')g_, 'in au'
   write(*,'(1(f20.12,1X),a)')g_*27210,'in meV'
!   write(*,*)

close(237)
deallocate(gord)
deallocate(sk,e1,g1)

stop

!! PARTIAL RATES

write(fname, '(A8,I1,A4)')"gamma_part.sh",ishift,".dat"

 do ichan = 1, nsta
 if(sum(gpart(:,ichan))==0d0) then
  gav = 0d0
  stadev = 0d0
  cycle
 endif

 write(*,*)
 write(*,*)"ISTATE=",ichan
 nmax = 28
! do k = 1, npt
! write(2111+ichan,'(2(f20.12,1X),a)')e(k),gpart(k,ichan)
! enddo

 allocate(gord(0:nmax))
 allocate(sk(0:nmax))
 allocate(e1(nmax,nmax), g1(nmax,nmax))
 gord = 0d0

 e(:) = e(:) + shift1
 call imaging(npt,e,gpart(:,ichan),nmax,nmax,e1,g1)
! write(*,*)ichan,nmax,sum(gpart(:,ichan))
 e(:) = e(:) - shift1

 do imax1 = 5, nmax
!   write(*,*)"Interp at order", imax1
   call interp(e1(:,imax1),g1(:,imax1),imax1-1,shift1,g_)
   g_ = 2.0d0*pi*g_
   gord(imax1) = g_
!  write(*,'(I3,A3,E23.15)')imax1," ",g_
 end do

 k = 0
 do i = 1, nmax
   do  j = 1, i-1
      write(100+i,'(2(f20.16,1X))')e1(j,i)-shift1,g1(j,i)*2.0d0*pi
      k = k + 1
      eallord(k) = e1(j,i)-shift1
      gallord(k) = g1(j,i)*2.0d0*pi
   enddo
 enddo

! gav = sum(gord(nmax-4:nmax-2))/3
! stadev = sqrt(sum( (gord(nmax-4:nmax-2)-gav)**2) )/3

 gav = sum(gord(20:25))/6
 stadev = sqrt(sum( (gord(20:25)-gav)**2) )/6

 st_gpart(ichan) = gav*27210

 write(*,'(2(f20.12,1X),a)')gav,stadev, 'in au'
 write(*,'(2(f20.12,1X),a)')gav*27210,stadev*27210, 'in meV'

call sort2(k,eallord,gallord)

do i = 1, k
 write(2220+ichan,'(2(f20.12,1X))')eallord(i),gallord(i)
enddo

 deallocate(gord)
 deallocate(sk,e1,g1)

 enddo ! ichan

! write(*,*)"absolute gamma in meV"
! write(*,'(2000(f20.12,1X))')gtot,sum(st_gpart),(st_gpart(i),i=1,nsta)
! write(*,*)"braching ratio in perc."
! write(*,'(2000(f20.12,1X))')gtot,100*sum(st_gpart)/gtot,(100*st_gpart(i)/gtot,i=1,nsta)

 deallocate(eallord,gallord)
 deallocate(st_gpart)
end program stieltjes
