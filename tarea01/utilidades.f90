module utilidades

contains
  ! indice menor igual a 9 digitos
  subroutine make_name_of_file_1_id(filename, header, sz, n)
    implicit none
    character(len=90), intent(inout) :: filename
    character(len=sz), intent(in) :: header
    integer, intent(in) :: sz, n
    integer :: digits
    character(len=90) :: fmt

    digits = ceiling(log10(real(n)))
    write(fmt, '(A,I1,A)' ) '( ("'//trim(header)//'",I', digits, ' ) )'
    write(filename, fmt ) n
  end subroutine make_name_of_file_1_id

  ! dos indices menores o iguales a 9 digitos
  subroutine make_name_of_file_2_id(filename, header, sz, n, m)
    implicit none
    character(len=90), intent(inout) :: filename
    character(len=sz), intent(in) :: header
    integer, intent(in) :: sz, n, m
    integer :: digits_1, digits_2
    character(len=90) :: fmt, hdr

    write(hdr,*) header 
    digits_1 = ceiling(log10(real(n)))
    digits_2 = ceiling(log10(real(m)))
    write(fmt, '(A,I1,A,I1,A)' ) '( ("'//trim(hdr(2:))//'",I', digits_1, '"_",I', digits_2, ' ) )'
    write(filename, fmt ) n, m
  end subroutine make_name_of_file_2_id

end module utilidades
