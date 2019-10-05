module NewtonRaphsonMod
  use VarPrecision
  implicit none

  TYPE newtRaphClass

  CONTAINS
    PROCEDURE, PUBLIC   ::    NewtRaph
  END TYPE

contains

  function NewtRaph ( this, input_guess, neqn, func ) result ( result_var )
    CLASS(newtRaphClass)                                      ::    this
    REAL(KIND=dp)          , DIMENSION( : ), INTENT( INOUT )  ::    input_guess
    INTEGER                                                   ::    neqn
    REAL(KIND=dp)          , DIMENSION( neqn )                ::    result_var
    REAL(KIND=dp)          , DIMENSION( neqn )                ::    g
    REAL(KIND=dp)          , DIMENSION( neqn )                ::    g_deriv
    REAL(KIND=dp)          , DIMENSION( neqn )                ::    delta
    REAL(KIND=dp)                                             ::    residual

    !! Create a function interface for a particular equation to find root
    INTERFACE
      function func( input_guess, neqn ) result ( func_result )
        INTEGER, parameter :: dp = selected_real_kind(8)
        INTEGER                                            ::    neqn
        REAL(KIND=dp), DIMENSION( neqn ), INTENT( IN )     ::    input_guess
        REAL(KIND=dp), DIMENSION( neqn )                   ::    func_result
      end function
    END INTERFACE

    residual = 1.0
    delta = 1e-16

    do while ( residual > 1e-8 )
      !! calculate the residual of the root finding function
      g = func( input_guess, neqn )
      !! calculate the gradient
      g_deriv = ( func( input_guess + delta, neqn ) - func( input_guess - delta, neqn ) ) / delta
      !! calculate the next set of values
      result_var = input_guess - g / g_deriv

      residual = MAXVAL( abs( ( result_var - input_guess ) / result_var ) )
      input_guess = result_var

      write(*,'(a,1f10.5,a,1f10.5)') 'NR Residual: ', residual, ' Value: ', result_var
    end do

  end function

end module
