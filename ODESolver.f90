module ODESolver
  use VarPrecision
  use ODEs
  implicit none

  TYPE, ABSTRACT    ::    odeSolverClass
    REAL(KIND=dp)                             ::    error_scale_out
    INTEGER                                   ::    order_accuracy = 1
    INTEGER                                   ::    neqn
  contains
    PROCEDURE( AbstractSolveTimestep  ), DEFERRED, PASS   ::    SolveTimestep
  end type

  ABSTRACT INTERFACE

    SUBROUTINE AbstractSolveTimestep( this, ode, input_values, curr_time, end_time, time_step, results_out )
      IMPORT odeSolverClass
      IMPORT abstractODEClass
      IMPORT dp
      CLASS(odeSolverClass), INTENT( IN )                       ::    this
      CLASS(abstractODEClass), INTENT( INOUT )                  ::    ode
      REAL(KIND=dp)          , INTENT( INOUT ), DIMENSION( : )  ::    input_values
      REAL(KIND=dp)          , INTENT( INOUT )                  ::    time_step
      REAL(kind=dp)          , INTENT( INOUT )                  ::    curr_time
      REAL(kind=dp)          , INTENT( IN )                     ::    end_time
      REAL(KIND=dp)          , DIMENSION( : )                   ::    results_out
    END SUBROUTINE

  END INTERFACE

  TYPE    odeContainerClass
    CLASS(odeSolverClass), pointer   ::    ode_solver
  END TYPE
contains


end module
