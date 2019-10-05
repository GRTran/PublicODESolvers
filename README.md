# PublicODESolvers
These are the suite of ODESolvers that have been developed in Fortran leveraging it's most recent OOP framework. 

The Lorenz system of ODEs are attached within this project directory as well as outputs from analysis to indicate the results of the chatotic system.

A simple validation was created using an ODE similar to dx/dt = ax which has a clearly defined analytic solution. The following ODE solving methods were implemented and compared for run-time and step size:

  - Forward Euler (explicit)
  - Backward Euler (implicit)
  - RK4 (explicit)
  - Shampine-Gordon (hybrid)
  
 The Shampine-Gordon solver performed optimally. This scheme is a hybrid scheme since it calculates the result for a series of future timesteps and determines the largest safe step size it can make whilst still producing a sensible solution.
