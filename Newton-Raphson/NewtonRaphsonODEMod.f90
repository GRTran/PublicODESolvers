module NewtonRaphsonODEMod
  use NewtonRaphsonMod
  use ODEs
  implicit none

  TYPE, EXTENDS(newtRaphClass)    ::    newtRaphODEClass
  CONTAINS
    PROCEDURE, PUBLIC   ::    NewtRaphODE
  END TYPE

contains

  function NewtRaphODE ( this, prev_val, ode, time, neqn, func ) result ( result_var )
    CLASS(newtRaphODEClass)                                   ::    this
    CLASS(abstractODEClass)                , INTENT( INOUT )  ::    ode
    INTEGER                                , INTENT( IN )     ::    neqn
    REAL(KIND=dp)                          , INTENT( IN )     ::    time
    REAL(KIND=dp)          , DIMENSION( : ), INTENT( IN )     ::    prev_val
    REAL(KIND=dp)          , DIMENSION( neqn )                   ::    root
    REAL(KIND=dp)          , DIMENSION( neqn )                ::    result_var
    REAL(KIND=dp)          , DIMENSION( neqn )                ::    g
    REAL(KIND=dp)          , DIMENSION( neqn )                ::    g_deriv
    REAL(KIND=dp)          , DIMENSION( neqn )                ::    delta
    REAL(KIND=dp)                                             ::    residual

    !! Create a function interface for a particular equation to find root
    INTERFACE
      function func( root, prev_val, ode, time, neqn ) result ( func_result )
        import abstractODEClass
        INTEGER, parameter :: dp = selected_real_kind(8)
        CLASS(abstractODEClass)         , INTENT( INOUT )         ::    ode
        INTEGER                         , INTENT( IN )            ::    neqn
        REAL(KIND=dp)                   , INTENT( IN )            ::    time
        REAL(KIND=dp), DIMENSION( neqn ), INTENT( IN )            ::    prev_val
        REAL(KIND=dp), DIMENSION( neqn ), INTENT( IN )            ::    root
        REAL(KIND=dp), DIMENSION( neqn )                          ::    func_result
      end function
    END INTERFACE


    !! the previous value serves as the guess for the next
    root = prev_val

    residual = 1.0
    delta = 1e-08

    do while ( residual > 1e-8 )
      !! calculate the residual of the root finding function
      g = func( root, prev_val, ode, time, neqn )
      !! calculate the gradient
      g_deriv = ( func( root+delta, prev_val, ode, time, neqn  ) - func( root-delta, prev_val, ode, time, neqn  ) ) / (2.0*delta)
      !! calculate the next set of values
      result_var = root - g / g_deriv

      residual = MAXVAL( abs( ( result_var - root ) / result_var ) )
      !write(*,*) abs( ( result_var - root ) / result_var )
      root = result_var
    end do
    !write(*,'(a,1f10.5,a,1f10.5)') 'NR Residual: ', residual, ' Value: ', result_var
  end function
end module
