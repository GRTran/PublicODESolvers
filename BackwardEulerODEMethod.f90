module BackwardEulerODEMethod
  use VarPrecision
  use ODESolver
  use ODEs
  use NewtonRaphsonODEMod
  implicit none

  type, extends(odeSolverClass)   ::    bkwdEulerSolverClass

  contains
    PROCEDURE, PUBLIC                         ::    SolveTimestep => BackEulSolveTimestep
  end type

contains

  !! -------------------------------------------------------------------------
  !! This method implements and explicit forward euler method to solve an ODE for given timestep
  !! -------------------------------------------------------------------------
  subroutine BackEulSolveTimestep( this, ode, input_values, curr_time, end_time, time_step, results_out )
    CLASS(bkwdEulerSolverClass), INTENT( IN )                 ::    this
    CLASS(abstractODEClass), INTENT( INOUT )                  ::    ode
    TYPE(newtRaphODEClass)                                    ::    newt_raph
    REAL(KIND=dp)          , INTENT( INOUT ), DIMENSION( : )  ::    input_values
    REAL(KIND=dp)          , INTENT( INOUT )                  ::    time_step
    REAL(kind=dp)          , INTENT( INOUT )                  ::    curr_time
    REAL(kind=dp)          , INTENT( IN )                     ::    end_time
    REAL(KIND=dp)          , DIMENSION( : )                   ::    results_out
    REAL(KIND=dp)          , ALLOCATABLE, DIMENSION( : )      ::    gradients_out
    REAL(KIND=dp)                                             ::    in_time

    !! allocate the gradients out array
    ALLOCATE(gradients_out(this%neqn))

    !! set the in subroutine time step variable
    in_time = time_step

    do while ( abs( end_time-curr_time ) > 1e-10 )

      if ( curr_time + in_time > end_time ) then
        in_time = end_time - curr_time
      endif
      !! calculate the gradients of the ode functions
      results_out = newt_raph%NewtRaphODE ( input_values, ode, in_time, size(input_values), BackwardEulFunc )

      ! set the input values equal to the results for next recursive iteration
      input_values = results_out
      !! increment the start time for passing into next recursive iteration
      curr_time = curr_time + in_time
    end do
    DEALLOCATE(gradients_out)
  end subroutine

  function BackwardEulFunc( root, prev_val, ode, time, neqn ) result ( result_var )
    CLASS(abstractODEClass)         , INTENT( INOUT )         ::    ode
    INTEGER                         , INTENT( IN )            ::    neqn
    REAL(KIND=dp)                   , INTENT( IN )            ::    time
    REAL(KIND=dp), DIMENSION( neqn ), INTENT( IN )            ::    prev_val
    REAL(KIND=dp), DIMENSION( neqn ), INTENT( IN )            ::    root
    REAL(KIND=dp), DIMENSION( neqn )                          ::    result_var
    REAL(KIND=dp), DIMENSION( neqn )                          ::    grad

    call ode%CalcGrad( time, root, grad )
    result_var = root - prev_val - time*grad
  end function

end module
