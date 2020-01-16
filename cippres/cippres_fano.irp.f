

program cippres_fano
  use general
  implicit none
  BEGIN_DOC
! CIPPRES stands for Configuration Interaction Plugin for Photoionized and Resonant Electronic States.
! ORMAS-CI calculations are performed to compute bound and approximate continuum states. Hamiltonian and Dipole coupling matrix elements can then be computed between these states. Finally, Stieltjes imaging technique is applied to recover the correct continuum coupling matrix elements.
! cippres_fano computes the H matrice couplings between the CI eigenvectors of ici1 and ici2 runs
  END_DOC

  integer :: ici

! TODO Read eigenvalues/vectors in EZFIO (do not recompute them)
! TODO Read the info (ici1, ici2,...) from an input

! TODO Print couplings in EZFIO

! GENERAL
! MANU : how to get input filename from command line qp run??


! TODO Compute dipole matrix elements between two different CI runs
! TODO Include Stieltjes in qp

  if(ifcsf==2) then

    !read(10,*) ici
    call ezfio_set_cippres_ici1(4)
    call ezfio_set_cippres_ici2(3)

    print*,ici1,ici2,n_ciruns_cippres
    print*, n_csf_cippres(ici1)
    print*, n_csf_cippres(ici2)

    print*, twoe_couplings_cippres(:,:)
    call ezfio_set_cippres_ifcsf(3)

  else 

    print*, "Please run cippres_runci first"

  endif

end program cippres_fano
