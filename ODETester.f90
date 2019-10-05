program ODETester
  use LorenzODES
  use ForwardEulerODEMethod
  use BackwardEulerODEMethod
  use RK4ODEMethod
  use ShampineGordonMethod
  use dataOutput

  implicit none

  type(fwdEulerSolverClass)                      ::    fwd_eul
  type(bkwdEulerSolverClass)                     ::    bkwd_eul
  type(rk4SolverClass)                           ::    rk4
  type(shampGordSolverClass)                           ::    sg
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
  INTEGER                                        ::   it_count
  REAL(KIND=8), DIMENSION(3)                     ::   ode_params

  ode_params = (/ 10.0d0, 28.0d0, 8.0d0/3.0d0 /)

  !! allocate the arrays the appropriate size
  fwd_eul%neqn = 3
  bkwd_eul%neqn = 3
  rk4%neqn = 3
  sg%neqn = 3
  allocate(input_values(fwd_eul%neqn))
  allocate(gradients(fwd_eul%neqn))
  allocate(outputs(fwd_eul%neqn))

  !! allocate the ode set to the Lorenz system
  allocate(lorenzODEClass::ode_set)
  !! set ode parameters
  call ode_set%SetODEParameters( ode_params )

  !! set the start and end time and the number of points to capture and the time step
  start_time = 0
  end_time = 100
  num_steps = 4000
  data_out_t_step = ( end_time - start_time ) / num_steps

  !! set the ode preferred time step
  pref_t_step = 0.002

  !! initialise the start and current timer
  data_out_end = start_time + data_out_t_step
  timer = start_time


  !! allocate the full output data array
  allocate(solved_data(num_steps,fwd_eul%neqn+1))

  !! specify the inital conditions
  input_values(1) = 1.0
  input_values(2) = 1.0
  input_values(3) = 1.0

  !! perform forward Euler method
  do it_count = 1, num_steps

    !! set time
    solved_data(it_count,1) = timer
    !! carry out the forward Euler
    call bkwd_eul%solveTimestep(ode_set, input_values, timer, data_out_end, pref_t_step, solved_data(it_count,2:fwd_eul%neqn+1))

    !! increment the end time, start time is already incremented
    data_out_end = data_out_end + data_out_t_step
  enddo
  !! write the data to a file
  call odeDataWriter(solved_data, 4, 145, 'output/lorenz.dat')
  call odeGnuplotDataPlotter('output/lorenz.dat', 146, 'output_ploter.plot')
end program
