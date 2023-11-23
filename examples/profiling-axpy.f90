PROGRAM main
  USE mpi
  USE xns_ptimer
  IMPLICIT NONE

  INTEGER, PARAMETER :: N = 1e8
  REAL(dp), ALLOCATABLE :: x(:), y(:)
  INTEGER :: rank, size, ierr
  INTEGER :: i

  ! Initialize MPI
  CALL mpi_init(ierr)
  CALL chrono%start("setup")
  CALL chrono%start("allocation")
  ALLOCATE(x(N), y(N))
  CALL chrono%break("allocation")
  CALL chrono%start("initialization")
  DO i=1, N
    x(i) = REAL(i, dp);
  END DO
  CALL chrono%break("initialization")
  CALL chrono%break("setup")
  CALL chrono%start("compute")
  CALL axpy(N, REAL(1.2345, dp), x, y) 
  CALL axpy(N, REAL(0.1234, dp), y, x) 
  CALL axpy(N, REAL(0.0123, dp), x, y) 
  CALL axpy(N, REAL(0.0012, dp), y, x) 
  CALL chrono%break("compute")

  CALL chrono%report()

  DEALLOCATE(x, y)

  ! Finalize MPI
  CALL mpi_finalize(ierr)

  CONTAINS

  SUBROUTINE axpy(n, a, x, y)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n
    REAL(dp), INTENT(IN) :: a
    REAL(dp), INTENT(INOUT) :: x(n), y(n)
    INTEGER :: i
  
    CALL chrono%start("axpy")
    CALL chrono%add_flops(2*REAL(n, dp))
    CALL chrono%add_bytes(2*8*REAL(n, dp))
    ! manually count flops and memory accesses
    DO i = 1, n
      y(i) = a*x(i) + y(i)
    END DO
    CALL chrono%break("axpy")
  END SUBROUTINE axpy
END PROGRAM main
