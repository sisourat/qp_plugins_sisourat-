use bitmasks ! you need to include the bitmasks_module.f90 features
use general

BEGIN_PROVIDER [integer, ici1]
 implicit none
    call ezfio_get_cippres_ici1(ici1)
END_PROVIDER 

BEGIN_PROVIDER [integer, ici2]
 implicit none
    call ezfio_get_cippres_ici2(ici2)
END_PROVIDER 

BEGIN_PROVIDER [integer, n_ciruns_cippres]
 implicit none
    call ezfio_get_cippres_n_ciruns_cippres(n_ciruns_cippres)
END_PROVIDER 

BEGIN_PROVIDER [integer, n_csf_max]
 implicit none
    call ezfio_get_cippres_n_csf_max(n_csf_max)
END_PROVIDER 

BEGIN_PROVIDER [integer, n_det_max_csf]
 implicit none
    call ezfio_get_cippres_n_det_max_csf(n_det_max_csf)
END_PROVIDER 

 BEGIN_PROVIDER [integer(bit_kind), csf_basis, (N_int,2, n_det_max_csf, n_csf_max ,n_ciruns_cippres )]
&BEGIN_PROVIDER [double precision, coef_det_csf_basis, (n_det_max_csf, n_csf_max ,n_ciruns_cippres )]
&BEGIN_PROVIDER [integer, n_csf_cippres, (n_ciruns_cippres )]
&BEGIN_PROVIDER [integer, n_sta_cippres, (n_ciruns_cippres )]
&BEGIN_PROVIDER [double precision, prttol_cippres, (n_ciruns_cippres )]
&BEGIN_PROVIDER [integer, n_det_csf_cippres, (n_csf_max,n_ciruns_cippres )]
 implicit none
 logical :: exists

  PROVIDE ezfio_filename HF_bitmask mo_coef
  csf_basis = 0_bit_kind
  if (mpi_master) then
      call ezfio_has_determinants_N_int(exists)
      if (exists) then
        call ezfio_has_determinants_bit_kind(exists)
        if (exists) then
          call ezfio_has_cippres_n_det_max_csf(exists)
          if (exists) then
            call ezfio_has_cippres_n_csf_max(exists)
            if (exists) then
              call ezfio_has_cippres_n_ciruns_cippres(exists)
              if (exists) then
                call ezfio_has_cippres_csf_basis(exists)
              endif
            endif
          endif
        endif
      endif
  endif

  if(exists)then ! All data for the csf basis exist and therefore will be read
   call read_csf_basis(csf_basis)
   call ezfio_get_cippres_coef_det_csf_basis(coef_det_csf_basis)
   call ezfio_get_cippres_n_csf_cippres(n_csf_cippres)
   call ezfio_get_cippres_n_sta_cippres(n_sta_cippres)
   call ezfio_get_cippres_prttol_cippres(prttol_cippres)
   call ezfio_get_cippres_n_det_csf_cippres(n_det_csf_cippres)
  else 
   print*,'The basis for CSF is not set in the ezfio file !!'
   print*,'You should run cippres_gencsf !'
   stop
  endif


END_PROVIDER 


BEGIN_PROVIDER [double precision, H_matrix_cippres, (n_csf_max,n_csf_max,n_ciruns_cippres)]
 implicit none
 integer :: irun
 integer :: i,j,k,l
 double precision :: hij
 H_matrix_cippres = 0.d0
 do irun = 1, n_ciruns_cippres
  do i = 1, n_csf_cippres(irun) ! first loop on the csf of the space ispace 
   do j = 1, n_csf_cippres(irun)
    do k = 1, n_det_csf_cippres(i,irun) ! then on the determinants belonging to the ith CSF of space ispace
     do l = 1, n_det_csf_cippres(j,irun)
      call i_H_j( csf_basis(1,1,k,i,irun) , csf_basis(1,1,l,j,irun),N_int,hij)
      H_matrix_cippres(j,i,irun) += hij * coef_det_csf_basis(k,i,irun) * coef_det_csf_basis(l,j,irun) 
     enddo
    enddo
   enddo
  enddo
 enddo
END_PROVIDER 


 BEGIN_PROVIDER [double precision, eigvectors_cippres, (n_csf_max,n_csf_max,n_ciruns_cippres)]
&BEGIN_PROVIDER [double precision, eigvalues_cippres, (n_csf_max,n_ciruns_cippres)]
 use general
 implicit none
 integer :: irun, j, i
 double precision, allocatable :: eigval(:),eigvec(:,:),hmat(:,:)

  do irun = 1, n_ciruns_cippres
    allocate(eigval(n_csf_cippres(irun)),eigvec(n_csf_cippres(irun),n_csf_cippres(irun)),hmat(n_csf_cippres(irun),n_csf_cippres(irun)))
!  print*,'nico',irun,H_matrix_cippres(1:n_csf_cippres(irun),1:n_csf_cippres(irun),irun)
    hmat(:,:) = H_matrix_cippres(1:n_csf_cippres(irun),1:n_csf_cippres(irun),irun)

    call lapack_diagd(eigval,eigvec,hmat,n_csf_cippres(irun),n_csf_cippres(irun)) 
    eigvalues_cippres(:,irun) = eigval(:)

    do j = 1, n_csf_cippres(irun)
      print*,'eigval',irun,j,eigval(j)
    enddo

    eigvectors_cippres(:,:,irun) = eigvec(:,:)
    do j = 1, n_csf_cippres(irun)
     do i = 1, n_csf_cippres(irun)
       print*, j,i,eigvectors_cippres(j,i,irun)
     enddo
    enddo
    deallocate(eigval,eigvec,hmat)
  enddo
 
 call ezfio_set_cippres_eigvalues_cippres(eigvalues_cippres)
 call ezfio_set_cippres_eigvectors_cippres(eigvectors_cippres)
END_PROVIDER 


 BEGIN_PROVIDER [double precision, twoe_couplings_cippres, (n_csf_max,n_csf_max)]
&BEGIN_PROVIDER [double precision, e_couplings_cippres, (n_csf_max,n_csf_max)]
 use general
 implicit none
 integer :: i, j, k, l
 double precision, allocatable :: eigval1(:),eigvec1(:,:),eigval2(:),eigvec2(:,:),twoe_csf_mat(:,:), twoe_mat(:,:)

 double precision :: hij

    twoe_couplings_cippres(:,:) = 0d0
    e_couplings_cippres(:,:) = 0d0

 allocate(twoe_csf_mat(n_csf_cippres(ici2),n_csf_cippres(ici1)))
 twoe_csf_mat(:,:) = 0d0

  do i = 1, n_csf_cippres(ici1) ! first loop on the csf of the space ispace 
   do j = 1, n_csf_cippres(ici2)
    do k = 1, n_det_csf_cippres(i,ici1) ! then on the determinants belonging to the ith CSF of space ispace
     do l = 1, n_det_csf_cippres(j,ici2)
      call i_H_j( csf_basis(1,1,k,i,ici1) , csf_basis(1,1,l,j,ici2),N_int,hij)
      twoe_csf_mat(j,i) += hij * coef_det_csf_basis(k,i,ici1) * coef_det_csf_basis(l,j,ici2)
     enddo
    enddo
!    print*,j,i,twoe_csf_mat(j,i)
   enddo
  enddo


 allocate(eigval1(n_csf_cippres(ici1)),eigval2(n_csf_cippres(ici2)))
 eigval1(:) = eigvalues_cippres(:,ici1)
 eigval2(:) = eigvalues_cippres(:,ici2)

 allocate(eigvec1(n_csf_cippres(ici1),n_csf_cippres(ici1)))
 allocate(eigvec2(n_csf_cippres(ici2),n_csf_cippres(ici2)))

 eigvec1(:,:) = eigvectors_cippres(:,:,ici1)
 eigvec2(:,:) = eigvectors_cippres(:,:,ici2)

 allocate(twoe_mat(n_csf_cippres(ici2),n_csf_cippres(ici1)))
 twoe_mat(:,:) = 0d0

 do i = 1, n_csf_cippres(ici1) ! first loop on the first eigenvectors
   do j = 1, n_csf_cippres(ici2) ! then on the second eigenvectors
      do k = 1, n_csf_cippres(ici1) ! loop over the csfs of the ici1 run
        do l = 1, n_csf_cippres(ici2) ! then over the csfs of the ici2 run
          twoe_mat(j,i) += twoe_csf_mat(l,k) * eigvec1(k,i) * eigvec2(l,j)
     enddo
    enddo
    e_couplings_cippres(j,i) = eigval1(i)-eigval2(j)
!    print*,eigval1(i)-eigval2(j), twoe_mat(j,i)**2
   enddo
  enddo

 twoe_couplings_cippres(:,:) = twoe_mat(:,:)**2
! print*,twoe_couplings_cippres(:,:)
 call ezfio_set_cippres_cfano_cippres(twoe_couplings_cippres)
 call ezfio_set_cippres_efano_cippres(e_couplings_cippres)

 deallocate(twoe_csf_mat,eigval1,eigval2,eigvec1,eigvec2,twoe_mat) 

 END_PROVIDER


 BEGIN_PROVIDER [double precision, cfano_cippres, (n_csf_max,n_csf_max)]
    call ezfio_get_cippres_cfano_cippres(cfano_cippres)
 END_PROVIDER

 BEGIN_PROVIDER [double precision, efano_cippres, (n_csf_max,n_csf_max)]
    call ezfio_get_cippres_efano_cippres(efano_cippres)
 END_PROVIDER