module ErrorCalcs
  implicit none

contains
  function ArrayError(array1, array2) result (err)
    double precision, intent(in) :: array1(:), array2(:)
    double precision :: err(size(array1))
    call arraySizeCheck(array1, array2)
    err = abs(array1 - array2)
  end

  function absArrayError(array1, array2) result (err)
    double precision, intent(in) :: array1(:), array2(:)
    double precision :: err
    call arraySizeCheck(array1, array2)
    err = maxval(abs(array1 - array2))
  end function absarrayError

  function lTwoArrayError(array1, array2) result (l2err)
    double precision, intent(in) :: array1(:), array2(:)
    double precision :: l2err
    call arraySizeCheck(array1, array2)
    l2err = (sum((array1-array2)**2)/size(array1))**0.5
    !write(*,*) size(array1), l2err
  end function lTwoArrayError

  subroutine arraySizeCheck(array1, array2)
    double precision, intent(in) :: array1(:), array2(:)
    if (size(array1) .NE. size(array2)) then
        write(*,*) "Arrays are not equal size"
        return
    endif
  end subroutine arraySizeCheck



end module ErrorCalcs
