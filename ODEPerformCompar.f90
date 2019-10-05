program ODETester
  use LorenzODES
  use BenchAxODE
  use ForwardEulerODEMethod
  use BackwardEulerODEMethod
  use RK4ODEMethod
  use ShampineGordonMethod
  use dataOutput
  use ErrorCalcs

  implicit none

  TYPE(odeContainerClass), DIMENSION(4)          ::   ode_array
  CLASS(abstractODEClass), POINTER               ::   ode_set
  REAL(KIND=8), allocatable, dimension(:)        ::   input_values
  REAL(KIND=8), allocatable, dimension(:)        ::   gradients
  REAL(KIND=8), allocatable, dimension(:)        ::   outputs
  REAL(KIND=8)                                   ::   start_time
  REAL(KIND=8)                                   ::   end_time
  INTEGER                                        ::   num_steps
  REAL(KIND=8)                                   ::   data_out_t_step
  REAL(KIND=8)                                   ::   pref_t_step
  REAL(KIND=8)                                   ::   data_out_end
  REAL(KIND=8)                                   ::   timer
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:)      ::   solved_data
  INTEGER                                        ::   it_count, solver, i, j, end_runs, start_runs
  REAL(KIND=8), DIMENSION(1)                     ::   ode_params
  REAL(KIND=8), DIMENSION(4)                     ::   cpu_timer, l2
  REAL(KIND=8)                                   ::   cpu_start, cpu_end
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:)      ::   log_e, log_n, log_t

  ode_params = (/ 0.1d0 /)

  !! set the start runs and end runs remembering that they increase to power of 2
  start_runs = 1                !! 2**4 = 16
  end_runs = 2                 !! 2**16 = 65536

  !! allocate the ode set to the Lorenz system
  allocate(benchAxODEClass::ode_set)

  !! allocate each ode solver method to the container array
  allocate(fwdEulerSolverClass::ode_array(1)%ode_solver)
  allocate(bkwdEulerSolverClass::ode_array(2)%ode_solver)
  allocate(rk4SolverClass::ode_array(3)%ode_solver)
  allocate(shampGordSolverClass::ode_array(4)%ode_solver)

  call ode_set%SetODEParameters( ode_params )
  do solver = 1, 4
    ode_array(solver)%ode_solver%neqn = 1
  enddo

  !! allocate the arrays the appropriate size
  allocate(input_values(ode_array(1)%ode_solver%neqn))
  allocate(gradients(ode_array(1)%ode_solver%neqn))
  allocate(outputs(ode_array(1)%ode_solver%neqn))

  !! set ode parameters
  call ode_set%SetODEParameters( ode_params )

  allocate(log_e( end_runs-start_runs, 4 ))
  allocate(log_n( end_runs-start_runs, 4 ))
  allocate(log_t( end_runs-start_runs, 4 ))
  j = 0
  !! Loop through the steps to get an error plot
  step_loop: do i = start_runs, end_runs-1
    !! set the number of points to capture
    j = j + 1
    !num_steps = 50*i
    num_steps = 100

    !! allocate the full output data array, for all solution methods, and the analytical solution
    allocate(solved_data(num_steps+1,ode_array(1)%ode_solver%neqn+1+4))

    !! set ode parameters
    call ode_set%SetODEParameters( ode_params )

    !! Perform this operation for each ODE solution method FwdEul, BkwdEul, RK4, SG
    method_loop: do solver = 1, 4

      !! perform association
      ASSOCIATE ( smethod => ode_array(solver)%ode_solver )

        !! set the start and end time and the time step
        start_time = 0
        end_time = 100
        data_out_t_step = ( end_time - start_time ) / num_steps

        !! set the ode preferred time step
        pref_t_step = 2.0*(1.0d0/i)
        write(*,*) pref_t_step
        !pref_t_step = 0.0002

        !! initialise the start and current timer
        data_out_end = start_time + data_out_t_step
        timer = start_time

        !! specify the inital conditions
        input_values(1) = 1.0

        !! Start the ODE solver time count
        call cpu_time(cpu_start)

        !! perform a solution method
        solution_loop: do it_count = 1, num_steps+1


          !! set time
          solved_data(it_count,1) = timer
          !! carry out the forward Euler
          call smethod%solveTimestep(ode_set, input_values, timer, data_out_end, pref_t_step, solved_data(it_count, smethod%neqn+solver : smethod%neqn+solver ) )

          solved_data( it_count, smethod%neqn+1+4 ) = 2.71828**(ode_params(1)*timer)

          !! increment the end time, start time is already incremented
          data_out_end = data_out_end + data_out_t_step
        enddo solution_loop

        !! end the time count and calculate the length of time taken
        call cpu_time(cpu_end)
        cpu_timer(solver) = cpu_end - cpu_start
        l2(solver) = lTwoArrayError(solved_data( :,smethod%neqn+solver ), solved_data( :,smethod%neqn+1+4 ) )
        write(*,*)'---------------------------------------------------------------'
        write(*,*)'Time taken: '
        write(*,*) solver, cpu_timer(solver)
        write(*,*) 'L2-Norm Error: '
        write(*,*) solver, l2(solver)
        write(*,*)'---------------------------------------------------------------'

        !! record the log error and log time data for plot
        log_e( j, solver ) = log(l2(solver))
        log_n (j, solver) = log(real(num_steps))
        log_t (j, solver) = log(pref_t_step)

      END ASSOCIATE
    enddo method_loop


    !! write the data to a file
    call odeDataWriter(solved_data, 4, 145, 'output/bench/benchtest.dat')
    call EXECUTE_COMMAND_LINE('python PythonPlot.py')
    call odeGnuplotDataPlotter('output/bench/benchtest.dat', 146, 'output/bench/output_ploter.plot')

    !! deallocate the solved data array for next loop
    deallocate(solved_data)
  enddo step_loop

  call odeDataWriter(log_e, 4, 145, 'output/bench/err_dat.dat')
  call odeDataWriter(log_n, 4, 145, 'output/bench/step_dat.dat')
  call odeDataWriter(log_t, 4, 145, 'output/bench/time_dat.dat')
  if ( end_runs > 2 ) call EXECUTE_COMMAND_LINE('python PythonPlotErr.py')

  do solver = 1, 4
    deallocate( ode_array(solver)%ode_solver )
  enddo
  deallocate(input_values)
  deallocate(gradients)
  deallocate(outputs)
  deallocate(log_e)
  deallocate(log_n)
  deallocate(log_t)

end program
