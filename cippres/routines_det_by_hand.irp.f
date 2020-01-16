subroutine create_det(n_a,n_b,occ_a,occ_b,det)
 use bitmasks ! you need to include the bitmasks_module.f90 features
 integer, intent(in) :: n_a,n_b,occ_a(mo_num),occ_b(mo_num)
 integer(bit_kind), intent(out) :: det(N_int,2)
 integer :: iorb,i
 det = 0_bit_kind
 do i = 1, n_a
  iorb = occ_a(i)
  call set_bit_to_integer(iorb,det(1,1),N_int)
 enddo
 do i = 1, n_b
  iorb = occ_b(i)
  call set_bit_to_integer(iorb,det(1,2),N_int)
 enddo
end
