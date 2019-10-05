module LorenzODES
  use VarPrecision
  use ODEs
  implicit none

  TYPE, extends(abstractODEClass)   ::    lorenzODEClass
    PRIVATE
    REAL(KIND=dp), POINTER   ::    sigma => NULL()
    REAL(KIND=dp), POINTER   ::    rho => NULL()
    REAL(KIND=dp), POINTER   ::    beta => NULL()
  CONTAINS
    PRIVATE
    PROCEDURE, PUBLIC, PASS   ::    SetODEParameters => LorenzSetODEParameters
    PROCEDURE                 ::    LorenzSetODEParameters
    PROCEDURE, PUBLIC, PASS   ::    CalcGrad => LorenzCalcGrad
    PROCEDURE                 ::    LorenzCalcGrad
    PROCEDURE, PUBLIC, PASS   ::    CalcTwoGrad => LorenzCalcTwoGrad
    PROCEDURE                 ::    LorenzCalcTwoGrad
    PROCEDURE, PUBLIC, PASS   ::    UpdateParameters => LorenzUpdateODEParameters
    PROCEDURE                 ::    LorenzUpdateODEParameters
    FINAL                     ::    DestroyLorenzODEs
  END TYPE

contains

  !! -------------------------------------------------------------------------
  !! Set the input parameters for the Lorenz System
  !! -------------------------------------------------------------------------
  subroutine LorenzSetODEParameters ( this, params_in )
    implicit none
    CLASS(lorenzODEClass), INTENT( OUT )                    ::   this
    REAL(kind=dp), INTENT( IN ) , DIMENSION( : )            ::   params_in
    INTEGER                                                 ::   ioerr = 0

    !! check only three inputs to the parameters
    if ( size(params_in) /= 3 ) ioerr = 1

    !! call error message subroutine to ensure successful subroutine execution
    call ErrMsg(ioerr)
    !! assign the parameters to the respective values
    this%sigma => params_in(1)
    this%rho => params_in(2)
    this%beta => params_in(3)
  end subroutine

  !! -------------------------------------------------------------------------
  !! Update the parameters in the Lorenz system if necessary
  !! -------------------------------------------------------------------------
  subroutine PoisUpdateODEParameters ( this, params_in )
    implicit none
    CLASS(lorenzODEClass), INTENT( OUT )                      ::   this
    REAL(kind=dp), INTENT( IN ) , DIMENSION( : ), TARGET    ::   params_in
    INTEGER                                                 ::   ioerr = 0

  end subroutine

  !! -------------------------------------------------------------------------
  !! The gradient calculation of the set of ODEs that describe a Lorenz system
  !! -------------------------------------------------------------------------
  subroutine LorenzCalcGrad ( this, curr_time, input_values, gradients_out )
    implicit none
    CLASS(lorenzODEClass), INTENT( IN )                          ::   this
    REAL(kind=dp)        , INTENT( IN )                          ::   curr_time
    REAL(kind=dp)        , INTENT( IN ) ,  DIMENSION( : )        ::   input_values
    REAL(kind=dp)        , INTENT( OUT ),  DIMENSION( : )        ::   gradients_out
    INTEGER                                                      ::   ioerr = 0

    if ( .NOT.ASSOCIATED( this%sigma ) .AND.                                  &
         .NOT.ASSOCIATED( this%rho )   .AND.                                  &
         .NOT.ASSOCIATED( this%beta ) ) ioerr = 2

    !! call error message subroutine to ensure successful subroutine execution
    call ErrMsg(ioerr)

    gradients_out(1) = this%sigma*(input_values(2)-input_values(1))
    gradients_out(2) = input_values(1)*(this%rho-input_values(3))-input_values(2)
    gradients_out(3) = input_values(1)*input_values(2)-this%beta*input_values(3)
  end subroutine

  !! -------------------------------------------------------------------------
  !! The second derivative calculation of the set of ODEs that describe a Lorenz system
  !! -------------------------------------------------------------------------
  subroutine LorenzCalcTwoGrad ( this, input_values, curr_time, second_deriv_out )
    implicit none
    CLASS(lorenzODEClass), INTENT( IN )                          ::   this
    REAL(kind=dp),         INTENT( IN ) ,  DIMENSION( : )        ::   input_values
    REAL(kind=dp),         INTENT( IN )                          ::   curr_time
    REAL(kind=dp), INTENT( OUT ), ALLOCATABLE, DIMENSION( : )    ::   second_deriv_out
    INTEGER                                                      ::   ioerr = 0

    allocate(second_deriv_out( size(input_values) ))

    if ( .NOT.ASSOCIATED( this%sigma ) .AND.                                  &
         .NOT.ASSOCIATED( this%rho )   .AND.                                  &
         .NOT.ASSOCIATED( this%beta ) ) ioerr = 2
    if ( .NOT.ALLOCATED( second_deriv_out ) ) ioerr = 3
    !! call error message subroutine to ensure successful subroutine execution
    call ErrMsg(ioerr)

    second_deriv_out(1) = - this%sigma
    second_deriv_out(2) = - 1
    second_deriv_out(3) = - this%beta
  end subroutine

  !! -------------------------------------------------------------------------
  !! Nullify the input parameters to the Lorenz system
  !! -------------------------------------------------------------------------
  subroutine DestroyLorenzODEs ( this )
    implicit none
    TYPE(lorenzODEClass)                    ::   this
    !! assign the parameters to the respective values
    if ( ASSOCIATED(this%sigma) ) DEALLOCATE(this%sigma)
    if ( ASSOCIATED(this%sigma) ) DEALLOCATE(this%rho)
    if ( ASSOCIATED(this%sigma) ) DEALLOCATE(this%beta)
  end subroutine

  !! -------------------------------------------------------------------------
  !! Check to see if any errors were encountered when executing the code
  !! -------------------------------------------------------------------------
  subroutine ErrMsg(err_type)
    INTEGER, INTENT( IN )   ::    err_type

    if ( err_type == 0 ) then
      return
    elseif (err_type == 1) then
      stop 'Error: Incorrect number of parameters input to Lorenz ODE Class'
    elseif (err_type == 2) then
      stop 'Error: Parameters have not been set for Lorenz ODE Class'
    elseif (err_type == 3) then
      stop 'Error: Gradient array input has not been allocated'
    endif
  end subroutine


end module
