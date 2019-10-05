module RK4ODEMethod
  use VarPrecision
  use ODESolver
  use ODEs
  implicit none

  type, extends(odeSolverClass)   ::    rk4SolverClass

  contains
    PROCEDURE, PUBLIC                         ::    SolveTimestep => RK4SolveTimestep
  end type

contains

  !! -------------------------------------------------------------------------
  !! This method implements and explicit forward euler method to solve an ODE for given timestep
  !! -------------------------------------------------------------------------
  subroutine RK4SolveTimestep( this, ode, input_values, curr_time, end_time, time_step, results_out )
    CLASS(rk4SolverClass)  , INTENT( IN )                     ::    this
    CLASS(abstractODEClass), INTENT( INOUT )                  ::    ode
    REAL(KIND=dp)          , INTENT( INOUT ), DIMENSION( : )  ::    input_values
    REAL(KIND=dp)          , INTENT( INOUT )                  ::    time_step
    REAL(kind=dp)          , INTENT( INOUT )                  ::    curr_time
    REAL(kind=dp)          , INTENT( IN )                     ::    end_time
    REAL(KIND=dp)                           , DIMENSION( : )  ::    results_out
    REAL(KIND=dp), ALLOCATABLE              , DIMENSION( : )  ::    gradients_out
    REAL(KIND=dp), ALLOCATABLE              , DIMENSION( :,: )::    k
    REAL(KIND=dp)                                             ::    in_time
    INTEGER                                                   ::    i

    !! allocate the gradients out array
    ALLOCATE(gradients_out(this%neqn))
    !! allocate the size of the slope array by the number of equations
    ALLOCATE( k( size(input_values) , 4) )

    !! set the in subroutine time step variable
    in_time = time_step

    do while ( abs( end_time-curr_time ) > 1e-10 )

      if ( curr_time + in_time > end_time ) then
        in_time = end_time - curr_time
      endif

      !! calculate slope at each position, weight them and multiply by the timestep
      slope_calc: do i = 1, 4
        if( i == 2 .OR. i == 3 ) then
          call ode%CalcGrad( curr_time+0.5d0*in_time, input_values+0.5d0*k(:,i-1), gradients_out )
          k(:,i) = 2.0d0 * in_time*gradients_out
        else
          call ode%CalcGrad( curr_time, input_values, gradients_out )
          k(:,i) = in_time*gradients_out
        endif
      end do slope_calc

      !! calculate the average slope and find the next term in the series
      results_out = input_values + SUM( k , 2) / 6.0d0
      !! set the input values equal to the results for next recursive iteration
      input_values = results_out
      !! increment the start time for passing into next recursive iteration
      curr_time = curr_time + in_time

    end do

    DEALLOCATE( k )
    DEALLOCATE(gradients_out)
  end subroutine
end module
