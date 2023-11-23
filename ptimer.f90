! ptimer: a parallel timer for high-performance computing.
! Copyright (C) 2023 Xavier Alvarez Farre
! This file is part of ptimer <https://github.com/xafarre/ptimer>.
! 
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! 
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License
! along with ptimer. If not, see <http://www.gnu.org/licenses/>.

! Brief
!  This module is a parallel timer designed for distributed-memory performance analysis. It does not deal yet
!  with OpenMP multithreaded codes. Please, use it carefully to prevent data races or blocks.
!
! Definition
!  In scientific computing, the performance of the codes is a crucial characteristic. The aim of this object
!  is to provide developers with a complete but rather simple tool to evaluate the performance of the codes
!  at user-specified control points. Every control point is transformed into a channel within ptimer, whose
!  hierarchy within the code is also stored.
!
!  By design, the only dependencies of ptimer are MPI and OpenMP, two programming standards for shared- and
!  distributed-memory parallelism.
!
!  Two methods are provided in ptimer class to control channels in user-specified control points: start and
!  break.
!   · start(name): given a channel name, the method searches its channel list to find whether the chanel
!   exists, or not. Thus, inserting a control point for the first time will create a new channel, while
!   reaching the same control point multiple times will just increase the channel's calls counter.
!
!   · break(name): given a channel name, the method searches its channel list to find whether the channel
!   exists and is open, or not. If the channel exists and is open, the time, flops and memory traffic that
!   ocurred after starting the channel will be accumulated. Otherwise, the application will crash.
!
!  Inserting a control point within an existing channel will create a new channel and consider it a nested
!  child of the existing channel. Besides, if the name given to a nested control point is also used outside
!  its parent, it will create a new channel with a different hierarchy. For instance, following the code
!  below, the output will be given as shown in the right:
!
!    1. ptimer%start(one point)        |   PTIMER
!    2.   ! some calculations          |    one point........INFO
!    3. ptimer%start(another)          |     another.........INFO
!    4.   ! more calculations          |    another..........INFO
!    5. ptimer%break(another)          |
!    6. ptimer%break(one point)        |
!    7. ptimer%start(another)          |
!    8.   ! more stuff                 |
!    9. ptimer%break(another)          |
!
!  This allows for fine-grained evaluation of kernels and functions. For instance, inside SpMV kernel the
!  user may be interested in breaking up its total elapsed time in update and compute times.
!
!  Channels accumulate the timing information by default. To store additional information such as the number
!  of operations or the memory traffic, the information must be added manually inside the control points. To
!  do so, two functions are provided: add_flops and and_bytes, which increment the variables flops and bytes,
!  respectively. When calling start, the channel will store the current value of these variables to, later
!  in break, calculate their increment. For instance, the following code will account for 10 flops and 80
!  bytes in the control point:
!
!    1. ptimer%start(mypoint)
!    2.   ptimer%add_flops(10)
!    3.   ptimer%add_bytes(80)
!    4. ptimer%break(mypoint)
!
!  Users are encouraged to use the methods above to analyze their kernels by counting or estimating the
!  performance and throughput within a control point.
!
! Pending
!  Extend the parallel timer to allow for multithreaded performance analysis. Calculate the overhead of
!  storing one channel per thread.
!
!  Extend the member method ptimer%report to write in output files a detailed report process by process. A
!  collective CSV file would be great.
MODULE xns_ptimer
  USE mpi
  IMPLICIT NONE

  INTEGER, PUBLIC, PARAMETER :: dp = SELECTED_REAL_KIND(15, 307)
  INTEGER, PRIVATE, PARAMETER :: MAX_CHANNEL_NUM = 256 ! default max number of channels
  INTEGER, PRIVATE, PARAMETER :: MAX_CHANNEL_LEN = 32 ! default max length of channel name

  ! This type is designed to store all the required information of a channel (times, operation count, memory
  ! traffic) as an array within ptimer.
  TYPE chan_t
    REAL(dp) :: timer = 0.0 ! time variables
    REAL(dp) :: flops = 0.0 ! flop variables
    REAL(dp) :: bytes = 0.0 ! data variables
    INTEGER :: level = 0 ! channel's nested level
    INTEGER :: outer = -1 ! id of parent channel
    INTEGER :: calls = 0 ! number of calls to channel
    LOGICAL :: open = .false. ! status
    CHARACTER(LEN=MAX_CHANNEL_LEN) :: name = "" ! channel's name
  END TYPE chan_t

  ! Actual type for the ptimer object.
  TYPE, PUBLIC :: ptimer
    PRIVATE

    REAL(dp) :: flops = 0.0 ! floating-point operations counter
    REAL(dp) :: bytes = 0.0 ! bytes counter
    INTEGER :: nch = 0 ! number of channels created
    INTEGER :: nop = 0 ! current number of open channels
    INTEGER :: openid = 0 ! current open channel id
    TYPE(chan_t), DIMENSION(MAX_CHANNEL_NUM) :: channels ! array of channels

    CONTAINS

    PROCEDURE :: start
    PROCEDURE :: break
    PROCEDURE :: report
    PROCEDURE :: reset
    PROCEDURE :: add_flops
    PROCEDURE :: add_bytes
  END TYPE ptimer

  TYPE(ptimer), PUBLIC :: chrono

  CONTAINS

  ! The member method ptimer%start resumes an existing channel or creates a new one. To determine whether a
  ! channel exists, both the given name and the current open channel, openid, are evaluated. This allows
  ! to use the same name in different nested regions. If the channel name given is the same as the currently
  ! open channel, ptimer will crash.
  SUBROUTINE start(this, name)
    CLASS(ptimer), INTENT(INOUT) :: this
    CHARACTER(LEN=*), INTENT(IN) :: name

    INTEGER :: ch

    ! To prevent nesting channels with the same name.
    IF(this%openid .NE. 0) THEN
        IF(this%channels(this%openid)%name .EQ. name) THEN
          STOP "channel is already open"
        END IF
    END IF

    ch = search(this, name)

    ! If channel is not found, create new channel.
    IF(ch .EQ. 0) THEN
      ch = create(this, name)
    END IF

    this%nop = this%nop + 1
    this%openid = ch ! last opened channel

    ! Initial time is evaluated at the end to prevent overheads from the search routine.
    this%channels(ch)%open = .true.
    this%channels(ch)%flops = this%channels(ch)%flops - this%flops
    this%channels(ch)%bytes = this%channels(ch)%bytes - this%bytes
    this%channels(ch)%timer = this%channels(ch)%timer - gettime()
  END SUBROUTINE start

  ! The member method ptimer%break pauses an existing channel. To determine whether a chanel exists, both
  ! the given name and the current open channel, openid, are evaluated. In the case a channel is not found,
  ! or is found but was already closed, ptimer will crash.
  SUBROUTINE break(this, name)
    CLASS(ptimer), INTENT(INOUT) :: this
    CHARACTER(LEN=*), INTENT(IN) :: name

    IF(this%openid .EQ. 0) THEN
      STOP "no channel open"
    ELSE IF(this%channels(this%openid)%name .NE. name) THEN
      STOP "channel is not open"
    END IF

    ! Store the time, flop, and data variables.
    this%channels(this%openid)%timer = this%channels(this%openid)%timer + gettime()
    this%channels(this%openid)%bytes = this%channels(this%openid)%bytes + this%bytes
    this%channels(this%openid)%flops = this%channels(this%openid)%flops + this%flops
    this%channels(this%openid)%open = .false.
    this%channels(this%openid)%calls = this%channels(this%openid)%calls + 1

    ! The outer of the current channel becomes the current open channel.
    this%openid = this%channels(this%openid)%outer
    this%nop = this%nop - 1
  END SUBROUTINE break

  ! The member method ptimer%report prints a report of the timing on the terminal. The output is adjusted to
  ! 80 characters width.
  SUBROUTINE report(this)
    CLASS(ptimer), INTENT(INOUT) :: this

    REAL(dp), DIMENSION(:), ALLOCATABLE :: tmax, tmin, tavg, favg, bavg
    REAL(dp) :: global
    INTEGER :: nch, i, j
    INTEGER :: world, rank, ierr
    CHARACTER(LEN=MAX_CHANNEL_LEN) :: pname

    nch = this%nch
    global = 0.0

    ALLOCATE(tmax(nch), tmin(nch), tavg(nch), favg(nch), bavg(nch))

    DO i = 1, nch
      tmax(i) = this%channels(i)%timer
      tmin(i) = this%channels(i)%timer
      tavg(i) = this%channels(i)%timer
      favg(i) = this%channels(i)%flops
      bavg(i) = this%channels(i)%bytes
    END DO

    CALL MPI_ALLREDUCE(MPI_IN_PLACE, tmax, nch, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, tmin, nch, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, tavg, nch, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, favg, nch, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, bavg, nch, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, world, ierr)

    DO i = 1, nch
      tavg(i) = tavg(i) / REAL(world, dp)
      favg(i) = favg(i) / REAL(world, dp)
      bavg(i) = bavg(i) / REAL(world, dp)

      IF (this%channels(i)%level == 0) THEN
        global = global + tavg(i)
      END IF
    END DO

    CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

    IF (rank == 0) THEN
      WRITE (*, '(G0)') "--------------------------------------------------------------------------------"
      WRITE (*, '(G0)') " PARALLEL TIMER REPORT"
      WRITE (*, '(G0)') "--------------------------------------------------------------------------------"
      WRITE (*, '(G0)') " Channel                             N GFLOPS   GB/s   tavg   tmax   tmin      %"
      WRITE (*, '(G0)') "--------------------------------------------------------------------------------"

      DO i = 1, nch
        pname = this%channels(i)%name

        DO j = 1, this%channels(i)%level
          pname = " " // pname
        END DO

        WRITE(*, '(A32, 1X, I5, 1X, F6.2, 1X, F6.2, 1X, F6.2, 1X, F6.2, 1X, F6.2, 1X, F6.2)') &
            pname, &
            this%channels(i)%calls, &
            favg(i)/tavg(i)/1.0E9, &
            bavg(i)/tavg(i)/1.0E9, &
            tavg(i), &
            tmax(i), &
            tmin(i), &
            tavg(i)/MAX(global, 1.0e-30)*100
      END DO
    END IF

    DEALLOCATE(tmax, tmin, tavg, favg, bavg)
  END SUBROUTINE report

  ! To reset all channels.
  SUBROUTINE reset(this)
    CLASS(ptimer), INTENT(INOUT) :: this

    INTEGER :: i

    DO i=1, MAX_CHANNEL_NUM
      this%channels(i)%timer = 0.0
      this%channels(i)%flops = 0.0
      this%channels(i)%bytes = 0.0
      this%channels(i)%level = 0
      this%channels(i)%outer = -1
      this%channels(i)%calls = 0
      this%channels(i)%open = .false.
      this%channels(i)%name = ""
    END DO
    this%flops = 0.0
    this%bytes = 0.0
    this%nch = 0
    this%nop = 0
    this%openid = 0
  END SUBROUTINE reset

  ! To add flops to counter
  SUBROUTINE add_flops(this, flops)
    CLASS(ptimer), INTENT(INOUT) :: this
    REAL(dp), INTENT(IN) :: flops

    this%flops = this%flops + flops
  END SUBROUTINE add_flops

  ! To add bytes to counter
  SUBROUTINE add_bytes(this, bytes)
    CLASS(ptimer), INTENT(INOUT) :: this
    REAL(dp), INTENT(IN) :: bytes

    this%bytes = this%bytes + bytes
  END SUBROUTINE add_bytes

  ! The SYSTEM_CLOCK is a portable solution for time evaluations. Since we are only interested in elapsed
  ! time, the method below provides a floating point value with respect to epoch. We accumulate differences
  ! of such a value.
  REAL(dp) FUNCTION gettime()
    INTEGER :: clock_count, clock_rate

    CALL SYSTEM_CLOCK(COUNT=clock_count, COUNT_RATE=clock_rate)
    gettime = REAL(clock_count, dp)/REAL(clock_rate, dp)
  END FUNCTION gettime

  ! To create new channels.
  INTEGER FUNCTION create(this, name)
    CLASS(ptimer), INTENT(INOUT) :: this
    CHARACTER(LEN=*), INTENT(IN) :: name

    IF(LEN(name) .EQ. 0) THEN
      STOP "channel name is empty"
    ELSE IF(LEN(name) .GE. MAX_CHANNEL_LEN) THEN
      STOP "channel name is too long"
    ELSE IF(this%nch .GE. MAX_CHANNEL_NUM) THEN
      STOP "reached maximum number of channels"
    END IF

    ! Since Fortran indexes from 1, we increase the variable nch before accessing the channel. The opposite is
    ! implemented in C++.
    this%nch = this%nch + 1
    create = this%nch

    ! Store the already open channel into outer list; if no channel is open, outer is 0.
    this%channels(this%nch)%name = name
    this%channels(this%nch)%level = this%nop
    this%channels(this%nch)%outer = this%openid
  END FUNCTION create

  ! To search a given channel.
  INTEGER FUNCTION search(this, name)
    CLASS(ptimer), INTENT(INOUT) :: this
    CHARACTER(LEN=*), INTENT(IN) :: name

    INTEGER :: ch

    DO ch=this%openid + 1, this%nch
      IF(compare(this, name, ch)) THEN
        search = ch
        RETURN
      END IF
    END DO
    search = 0
  END FUNCTION search

  ! To compare a given channel.
  LOGICAL FUNCTION compare(this, name, ch)
    CLASS(ptimer), INTENT(INOUT) :: this
    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER, INTENT(IN)  :: ch

    compare = ((this%channels(ch)%name .EQ. name) .AND. (this%channels(ch)%outer .EQ. this%openid))
  END FUNCTION compare
END MODULE xns_ptimer
