module BenchAxODE
  use VarPrecision
  use ODEs
  implicit none

  TYPE, extends(abstractODEClass)   ::    benchAxODEClass
    PRIVATE
    REAL(KIND=dp), POINTER   ::    a => NULL()
  CONTAINS
    PRIVATE
    PROCEDURE, PUBLIC, PASS   ::    SetODEParameters => BenchAxSetODEParameters
    PROCEDURE                 ::    BenchAxSetODEParameters
    PROCEDURE, PUBLIC, PASS   ::    CalcGrad => BenchAxCalcGrad
    PROCEDURE                 ::    BenchAxCalcGrad
    PROCEDURE, PUBLIC, PASS   ::    CalcTwoGrad => BenchAxCalcTwoGrad
    PROCEDURE                 ::    BenchAxCalcTwoGrad
    FINAL                     ::    DestroyBenchAxODEs
  END TYPE

contains

  !! -------------------------------------------------------------------------
  !! Set the input parameters for the Lorenz System
  !! -------------------------------------------------------------------------
  subroutine BenchAxSetODEParameters ( this, params_in )
    implicit none
    CLASS(benchAxODEClass), INTENT( OUT )                   ::   this
    REAL(kind=dp), INTENT( IN ) , DIMENSION( : ), TARGET    ::   params_in
    INTEGER                                                 ::   ioerr = 0

    !! check only three inputs to the parameters
    if ( size(params_in) /= 1 ) ioerr = 1

    !! call error message subroutine to ensure successful subroutine execution
    call ErrMsg(ioerr)

    !! assign the parameters to the respective values
    this%a => params_in(1)


  end subroutine

  !! -------------------------------------------------------------------------
  !! The gradient calculation of the set of ODEs that describe a Lorenz system
  !! -------------------------------------------------------------------------
  subroutine BenchAxCalcGrad ( this, curr_time, input_values, gradients_out )
    implicit none
    CLASS(benchAxODEClass), INTENT( IN )                          ::   this
    REAL(kind=dp)         , INTENT( IN )                          ::   curr_time
    REAL(kind=dp)         , INTENT( IN ) ,  DIMENSION( : )        ::   input_values
    REAL(kind=dp)         , INTENT( OUT ),  DIMENSION( : )        ::   gradients_out
    INTEGER                                                       ::   ioerr = 0

    if ( .NOT.ASSOCIATED( this%a ) ) ioerr = 2

    !! call error message subroutine to ensure successful subroutine execution
    call ErrMsg(ioerr)

    gradients_out(1) = this%a*input_values(1)
  end subroutine

  !! -------------------------------------------------------------------------
  !! The second derivative calculation of the set of ODEs that describe a Lorenz system
  !! -------------------------------------------------------------------------
  subroutine BenchAxCalcTwoGrad ( this, input_values, curr_time, second_deriv_out )
    implicit none
    CLASS(benchAxODEClass), INTENT( IN )                          ::   this
    REAL(kind=dp)         , INTENT( IN ) ,  DIMENSION( : )        ::   input_values
    REAL(kind=dp)         , INTENT( IN )                          ::   curr_time
    REAL(kind=dp), INTENT( OUT ), ALLOCATABLE, DIMENSION( : )     ::   second_deriv_out
    INTEGER                                                       ::   ioerr = 0

    allocate(second_deriv_out( size(input_values) ))

    if ( .NOT.ASSOCIATED( this%a ) ) ioerr = 2
    if ( .NOT.ALLOCATED( second_deriv_out ) ) ioerr = 3
    !! call error message subroutine to ensure successful subroutine execution
    call ErrMsg(ioerr)

    second_deriv_out(1) = this%a
  end subroutine

  !! -------------------------------------------------------------------------
  !! Nullify the input parameters to the Lorenz system
  !! -------------------------------------------------------------------------
  subroutine DestroyBenchAxODEs ( this )
    implicit none
    TYPE(benchAxODEClass)                   ::   this
    !! assign the parameters to the respective values
    ! if ( ASSOCIATED(this%a) ) DEALLOCATE(this%a)
  end subroutine

  !! -------------------------------------------------------------------------
  !! Check to see if any errors were encountered when executing the code
  !! -------------------------------------------------------------------------
  subroutine ErrMsg(err_type)
    INTEGER, INTENT( IN )   ::    err_type

    if ( err_type == 0 ) then
      return
    elseif (err_type == 1) then
      stop 'Error: Incorrect number of parameters input to BenchAx ODE Class'
    elseif (err_type == 2) then
      stop 'Error: Parameters have not been set for BenchAx ODE Class'
    elseif (err_type == 3) then
      stop 'Error: Gradient array input has not been allocated'
    endif
  end subroutine


end module
