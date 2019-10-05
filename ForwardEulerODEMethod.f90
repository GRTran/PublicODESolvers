module ForwardEulerODEMethod
  use VarPrecision
  use ODESolver
  use ODEs
  implicit none

  type, extends(odeSolverClass)   ::    fwdEulerSolverClass

  contains
    PROCEDURE, PUBLIC                         ::    SolveTimestep => FwdEulSolveTimestep
  end type

contains

  !! -------------------------------------------------------------------------
  !! This method implements and explicit forward euler method to solve an ODE for given timestep
  !! -------------------------------------------------------------------------
  subroutine FwdEulSolveTimestep( this, ode, input_values, curr_time, end_time, time_step, results_out )
    CLASS(fwdEulerSolverClass), INTENT( IN )                  ::    this
    CLASS(abstractODEClass), INTENT( INOUT )                  ::    ode
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
      !! calculate the gradients of the ode functions
      call ode%CalcGrad( curr_time, input_values, gradients_out )

      if ( curr_time + in_time > end_time ) then
        in_time = end_time - curr_time
      endif
      !! perform forward euler
      results_out = input_values + in_time*gradients_out
      !! set the input values equal to the results for next recursive iteration
      input_values = results_out
      !write(*,*) results_out
      !! increment the start time for passing into next recursive iteration
      curr_time = curr_time + in_time
    end do

    DEALLOCATE(gradients_out)

  end subroutine
end module
