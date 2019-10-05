module ShampineGordonMethod
  use VarPrecision
  use ODESolver
  use ODEs
  use ShampGordon
  implicit none

  type, extends(odeSolverClass)   ::    shampGordSolverClass

  contains
    PROCEDURE, PUBLIC                         ::    SolveTimestep => ShampSolveTimestep

  end type

contains

  !! -------------------------------------------------------------------------
  !! This method implements the shampine gordon method to solve an ODE for given timestep
  !! using an externally sourced ode solver
  !! -------------------------------------------------------------------------
  subroutine ShampSolveTimestep( this, ode, input_values, curr_time, end_time, time_step, results_out )
    implicit none
    CLASS(shampGordSolverClass), INTENT( IN )                 ::    this
    CLASS(abstractODEClass), INTENT( INOUT )                  ::    ode
    REAL(KIND=dp)          , INTENT( INOUT ), DIMENSION( : )  ::    input_values
    REAL(KIND=dp)          , INTENT( INOUT )                  ::    time_step
    REAL(kind=dp)          , INTENT( INOUT )                  ::    curr_time
    REAL(kind=dp)          , INTENT( IN )                     ::    end_time
    REAL(KIND=dp)          , DIMENSION( : )                   ::    results_out
    REAL(KIND=dp)          , ALLOCATABLE, DIMENSION( : )      ::    gradients_out
    REAL(KIND=dp)                                             ::    in_time
    INTEGER                                                   ::    iflag
    INTEGER                                                   ::    iwork(5)
    REAL(KIND=dp)          , ALLOCATABLE, DIMENSION( : )      ::    work
    REAL(KIND=dp)                                             ::    abserr
    REAL(KIND=dp)                                             ::    relerr

    !! set the flag to -1 to catch any errors in the ode solver
    iflag =  1
    !! put in the error tolerances
    abserr = 1e-8
    relerr = 1e-8
    !! allocate memory for the solver
    ALLOCATE( work( 100+21*size(input_values) ) )
    !! allocate the gradients out array
    ALLOCATE(gradients_out(this%neqn))

    call shampODE( GradCalculatorInterface, ode, this%neqn, input_values, curr_time, end_time, relerr, abserr, iflag, work, iwork)
    results_out = input_values
    if (iflag /= 2) stop 'Error Shamp Gordon Solver not solved'


    DEALLOCATE(work)
    DEALLOCATE(gradients_out)
  end subroutine

  !! -------------------------------------------------------------------------
  !! This wrapper is used for the calcultion of the gradient
  !! -------------------------------------------------------------------------
  subroutine GradCalculatorInterface( ode, t, y, yp, neqn )
    implicit none
    INTEGER                         , INTENT( IN )            ::    neqn
    CLASS(abstractODEClass), INTENT( INOUT )                  ::    ode
    REAL(KIND=dp)                   , INTENT( IN )            ::    t
    REAL(KIND=dp), DIMENSION( neqn ), INTENT( INOUT )         ::    y
    REAL(KIND=dp) , DIMENSION( neqn )      ::    yp
    call ode%CalcGrad( t, y, yp )
  end subroutine

end module
