program NewtRaphTester
  use NewtonRaphsonMod
  implicit none

  TYPE(newtRaphClass)   ::    newt_raph
  REAL(KIND=dp), DIMENSION(1)     ::    x_guess, x_out
  INTEGER                         ::    neqn

  x_guess = 10000.0d0
  neqn = 1
  x_out = newt_raph%NewtRaph(x_guess, neqn, QuadFunc)
  write(*,'(a,1f10.5)') 'result: ',x_out
contains

  function QuadFunc( x , neqn ) result ( result_var )
    INTEGER                                                   ::    neqn
    REAL(KIND=dp), DIMENSION( neqn ), INTENT( IN )     ::    x
    REAL(KIND=dp), DIMENSION( neqn )    ::    result_var
    REAL(KIND=dp)   ::    mean, var, g

    mean = 82362231336.050201
    var = 3.1998599e20
    g = 2.902086e-6
    !result_var = x**2.0 + 5*x + 6    !! root is x = -2
    ! result_var = -4.0*x**2 + 12*x - 9 !! root is x = 1.5
    result_var = log( 1+ x*log( (1/(1+x)) * (var / mean + x ) )) - ( x / ((1+x)*mean))* (var/mean - 1)*log(1/g)
  end function

end program
