module PoiseuilleFlowODE
!! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !!
!! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !!
!! Title:   Poiseuille tubular flow ode                     Date: 14/06/2019 !!
!! Author:  Greg Jones, Imperial College London                              !!
!! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !!
!! Description: This module uses a model of poiseuille flow to characterise  !!
!!              the height of a wetted powder region as time progresses      !!
!!              through a transient.                                         !!
!! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !!
!! Dependencies: ODE abstract class, general precision module.               !!
!! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !!
!! Revision: 1.0 (12/06/2019)                                                !!
!!           Original functionality implemented in module includes the       !!
!!           development subroutines to set poiseuille flow parameters.      !!
!!           The calculation of the gradient is also provided here.          !!
!! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !!
!! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !!
  use VarPrecision
  use ODEs
  implicit none

  TYPE, extends(abstractODEClass)   ::    poisODEClass
    PRIVATE
    REAL(KIND=dp)   ::    rhow
    REAL(KIND=dp)   ::    perm
    REAL(KIND=dp)   ::    mu
    REAL(KIND=dp)   ::    porosity
    REAL(KIND=dp)   ::    water_h
    REAL(KIND=dp)   ::    capillary_h
  CONTAINS
    PRIVATE
    PROCEDURE, PUBLIC, PASS   ::    SetODEParameters => PoisSetODEParameters
    PROCEDURE                 ::    PoisSetODEParameters
    PROCEDURE, PUBLIC, PASS   ::    CalcGrad => PoisCalcGrad
    PROCEDURE                 ::    PoisCalcGrad
    PROCEDURE, PUBLIC, PASS   ::    CalcTwoGrad => PoisCalcTwoGrad
    PROCEDURE                 ::    PoisCalcTwoGrad
    PROCEDURE, PUBLIC, PASS   ::    UpdateParameters => PoisUpdateODEParameters
    PROCEDURE                 ::    PoisUpdateODEParameters
    FINAL                     ::    DestroyPoisODEs
  END TYPE

contains

  !! -------------------------------------------------------------------------
  !! Set the input parameters for the Poiseuille flow ode. Params in order are
  !! as follows:
  !! (1) water density, (2) permeability
  !! (3) dynamic viscosity, (4) porosity, (5) starting height of water pond
  !! (6) capillary height
  !! -------------------------------------------------------------------------
  subroutine PoisSetODEParameters ( this, params_in )
    implicit none
    CLASS(poisODEClass), INTENT( OUT )                      ::   this
    REAL(kind=dp), INTENT( IN ) , DIMENSION( : ), TARGET    ::   params_in
    INTEGER                                                 ::   ioerr = 0
    INTEGER                                                 ::   i

    !! check only three inputs to the parameters
    if ( size(params_in) /= 6 ) ioerr = 1

    !! call error message subroutine to ensure successful subroutine execution
    call ErrMsg(ioerr)

    !! assign the parameters to the respective values
    this%rhow = params_in(1)
    this%perm = params_in(2)
    this%mu = params_in(3)
    this%porosity = params_in(4)
    this%water_h = params_in(5)
    this%capillary_h = params_in(6)
  end subroutine

  !! -------------------------------------------------------------------------
  !! Update the parameters in the Poiseuille system, the first argument is the
  !! height of the water above the wetted region.
  !! -------------------------------------------------------------------------
  subroutine PoisUpdateODEParameters ( this, params_in )
    implicit none
    CLASS(poisODEClass), INTENT( OUT )                      ::   this
    REAL(kind=dp), INTENT( IN ) , DIMENSION( : )            ::   params_in
    INTEGER                                                 ::   ioerr = 0

    this%water_h = params_in(1)
  end subroutine

  !! -------------------------------------------------------------------------
  !! The gradient calculation of the set of ODEs that describe the flow of
  !! water through the tubules of a porous medium. The gradient term represents
  !! the change in height of the wetted porous medium region
  !! -------------------------------------------------------------------------
  subroutine PoisCalcGrad ( this, curr_time, input_values, gradients_out )
    implicit none
    CLASS(poisODEClass), INTENT( IN )                            ::   this
    REAL(kind=dp)        , INTENT( IN )                          ::   curr_time
    REAL(kind=dp)        , INTENT( IN ) ,  DIMENSION( : )        ::   input_values
    REAL(kind=dp)        , INTENT( OUT ),  DIMENSION( : )        ::   gradients_out
    INTEGER                                                      ::   ioerr = 0
    INTEGER                                                      ::   n_inputs
    REAL(kind=dp)                                                ::   prec_conc


    !! calculate the gradient term for the changing wetted region height
    gradients_out(1) = ( ( PI * grav * this%rhow * this%perm ) / &
                       ( 8.0d0 * this%mu * this%porosity ) ) * &
                       ( ( this%water_h - input_values(1)*this%porosity &
                         + input_values(1) + this%capillary_h ) / input_values(1) )
  end subroutine

  !! -------------------------------------------------------------------------
  !! The second derivative calculation has not been implemented here
  !! -------------------------------------------------------------------------
  subroutine PoisCalcTwoGrad ( this, input_values, curr_time, second_deriv_out )
    implicit none
    CLASS(poisODEClass), INTENT( IN )                          ::   this
    REAL(kind=dp),         INTENT( IN ) ,  DIMENSION( : )        ::   input_values
    REAL(kind=dp),         INTENT( IN )                          ::   curr_time
    REAL(kind=dp), INTENT( OUT ), ALLOCATABLE, DIMENSION( : )    ::   second_deriv_out
    INTEGER                                                      ::   ioerr = 0

  end subroutine

  !! -------------------------------------------------------------------------
  !! Nullify the input parameters to the Lorenz system
  !! -------------------------------------------------------------------------
  subroutine DestroyPoisODEs ( this )
    implicit none
    TYPE(poisODEClass)                    ::   this
    !! assign the parameters to the respective values

  end subroutine

  !! -------------------------------------------------------------------------
  !! Check to see if any errors were encountered when executing the code
  !! -------------------------------------------------------------------------
  subroutine ErrMsg(err_type)
    INTEGER, INTENT( IN )   ::    err_type

    if ( err_type == 0 ) then
      return
    elseif (err_type == 1) then
      stop 'Error: Incorrect number of parameters input to Poiseuille ODE Class'
    elseif (err_type == 2) then
      stop 'Error: Parameters have not been set for Poiseuille ODE Class'
    elseif (err_type == 3) then
    endif
  end subroutine


end module
