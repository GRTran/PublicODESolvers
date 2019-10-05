module ODEs
  use VarPrecision
  implicit none


  TYPE, ABSTRACT    ::    abstractODEClass

  CONTAINS
    PROCEDURE(AbstractSetODEParameters), DEFERRED, PASS   ::    SetODEParameters
    PROCEDURE(AbstractCalcGrad),         DEFERRED, PASS   ::    CalcGrad
    PROCEDURE(AbstractCalcTwoGrad),      DEFERRED, PASS   ::    CalcTwoGrad
    PROCEDURE(AbstractUpdateParameters), DEFERRED, PASS   ::    UpdateParameters
  END TYPE

  !! this interface will execute the function associated with the type
  !! in which the dummy variable is dynamically allocated
  ABSTRACT INTERFACE

    !! -------------------------------------------------------------------------
    !! Set the input parameters for a system of ODEs
    !! -------------------------------------------------------------------------
    subroutine AbstractSetODEParameters ( this, params_in )
      import abstractODEClass
      implicit none
      INTEGER, parameter :: dp = selected_real_kind(8)
      CLASS(abstractODEClass),  INTENT( OUT )                    ::   this
      REAL(kind=dp),    INTENT( IN ) , DIMENSION( : ), TARGET    ::   params_in
    end subroutine

    !! -------------------------------------------------------------------------
    !! Update the parameters in the system if necessary
    !! -------------------------------------------------------------------------
    subroutine AbstractUpdateParameters ( this, params_in )
      import abstractODEClass
      implicit none
      INTEGER, parameter :: dp = selected_real_kind(8)
      CLASS(abstractODEClass), INTENT( OUT )                      ::   this
      REAL(kind=dp), INTENT( IN ) , DIMENSION( : )                ::   params_in
    end subroutine

    !! -------------------------------------------------------------------------
    !! The gradient calculation of a set of ODEs
    !! -------------------------------------------------------------------------
    subroutine AbstractCalcGrad ( this, curr_time,input_values, gradients_out )
      import abstractODEClass
      implicit none
      INTEGER, parameter :: dp = selected_real_kind(8)
      CLASS(abstractODEClass),  INTENT( IN )                     ::   this
      REAL(kind=dp),            INTENT( IN ) ,  DIMENSION( : )   ::   input_values
      REAL(kind=dp),            INTENT( IN )                     ::   curr_time
      REAL(kind=dp), INTENT( OUT ), DIMENSION( : )    ::   gradients_out
    end subroutine

    !! -------------------------------------------------------------------------
    !! The second derivative calculation for a set of ODEs
    !! -------------------------------------------------------------------------
    subroutine AbstractCalcTwoGrad ( this, input_values, curr_time, second_deriv_out )
      import abstractODEClass
      implicit none
      INTEGER, parameter :: dp = selected_real_kind(8)
      CLASS(abstractODEClass),  INTENT( IN )                     ::   this
      REAL(kind=dp),            INTENT( IN ) ,  DIMENSION( : )   ::   input_values
      REAL(kind=dp),            INTENT( IN )                     ::   curr_time
      REAL(kind=dp), INTENT( OUT ), ALLOCATABLE, DIMENSION( : )    ::   second_deriv_out
    end subroutine

  END INTERFACE

contains

end module
