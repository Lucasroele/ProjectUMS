**==readdat.spg  processed by SPAG 4.52O  at 18:54 on 27 Mar 1996
      SUBROUTINE READDAT(Equil, Prod, Nsamp, Ndispl, Dr, Iseed, nLambda, nGhosts,
     & nWidomCycle, runWidom, runTDI)
C     ---read input data and model parameters
c
c     ---input parameters: file: fort.15
c    ibeg  =  0  : initialize from a lattice
c             1  : read configuration from disk
c    Equil       : number of Monte Carlo cycles during equilibration
c    Prod        : number of Monte Carlo cycles during production
c    Nsamp       : number of Monte Carlo cycles between two sampling periods
c    Iseed       : seed random number generator
c    Dr          : maximum displacement
c    Ndispl      : number of attemps to displace a particle per MC cycle
c    NPART       : total numbero fo particles
c    TEMP        : temperature
c    rho         : density
c    nLambda     : number of distinct lambdas to simulate at
c    nWidomCycle : number of Monte Carlo cycles during Widom production
c    nGhosts     : number of particle insertion attempts
c    runWidom    : 0 dont run Widom
c                : 1 run Widom
c    runTDI      : 0 dont run TDI cycle
c                : 1 run TDI
c
c     ---input parameters: file: fort.25
c    TAILCO = .true. : tail corrections are applied
c             .false.: no-tail corrections
c    SHIFT  = .true. : potential is shifted
c             .false.: potential is NOT-shifted
c    eps    = epsilon Lennard-Jones potential
c    sig    = sigma Lennard-Jones potential
c    MASS   = mass of the particle
c    RC     = cut-off radius of the potential
c
c     ---input parameters: file: fort.11 (restart file
c                to continue a simulation from disk)
c    boxf   = Box length old configuration (if this one
c             does not correspond to the requested density, the positions
c             of the particles are rescaled!
c    NPART  = number of particles (over rules fort.15!!)
c    Dr     = optimized maximum displacement old configurations
c    X(1),Y(1),Z(1)            : position first particle 1
c        ...
c    X(NPART),Y(NPART),Z(NPART): position particle last particle
 
      IMPLICIT NONE
      INCLUDE 'parameter.inc'
      INCLUDE 'system.inc'
      INCLUDE 'potential.inc'
      INCLUDE 'conf.inc'
      INTEGER ibeg, Equil, Prod, i, Ndispl, Nsamp, Iseed, nLambda, nGhosts,
     &       nWidomCycle, runWidom, runTDI
      DOUBLE PRECISION eps, sig, CORU, CORP, vir, boxf, rhof, rho, Dr
     
c     ---read simulation data
      READ (15, *)
      READ (15, *) ibeg, Equil, Prod, Nsamp, Iseed
      READ (15, *)
      READ (15, *) Dr
      READ (15, *)
      READ (15, *) Ndispl, nLambda, nGhosts, nWidomCycle
      READ (15, *)
      READ (15, *) NPART, TEMP, rho
      READ (15, *)
      READ (15, *) runWidom, runTDI, sig
c     ---initialise and test random number generator
      CALL RANTEST(Iseed)
c     --- DOEN???
      Dr = Dr / sig
 
      IF (NPART.GT.NPMax) THEN
         WRITE (10, *) ' ERROR: number of particles too large'
         STOP
      END IF
      BOX = (NPART/rho)**(1.D0/3.D0)
      HBOX = BOX/2
c     ---read model parameters
      READ (25, *)
      READ (25, *) TAILCO, SHIFT
      READ (25, *)
      READ (25, *) eps, MASS, RC
c     ---read/generate configuration
      IF (ibeg.EQ.0) THEN
c        ---generate configuration form lattice
         CALL LATTICE

c     ---skip###########
      ELSE
         WRITE (10, *) ' read conf from disk '
         READ (11, *) boxf
         READ (11, *) NPART
         READ (11, *) Dr
         WRITE (10,*) 'boxf ', boxf, ' NPART ', NPART, ' Dr ', Dr
         rhof = NPART/boxf**3
         IF (ABS(boxf-BOX).GT.1D-6) THEN
            WRITE (10, 99007) rho, rhof
         END IF
         DO i = 1, NPART
            READ (11, *) X(i), Y(i), Z(i)
            X(i) = X(i)*BOX/boxf
            Y(i) = Y(i)*BOX/boxf
            Z(i) = Z(i)*BOX/boxf
         END DO
         REWIND (11)
      END IF
c     ---skip###########

c     ---write input data
      WRITE (10, 99001) Equil, Prod, Nsamp
      WRITE (10, 99002) Ndispl, Dr
      WRITE (10, 99003) NPART, TEMP, rho, BOX
      WRITE (10, 99004) eps, sig, MASS
      WRITE (10, 99008) nLambda
c     ---calculate parameters:
      BETA = 1/TEMP
      PI = 4*ATAN(1.D0)
c     ---calculate cut-off radius potential
      RC = MIN(RC, HBOX)
      RC2 = RC*RC
      EPS4 = 4*eps
      EPS48 = 48*eps
      SIG2 = sig*sig
      IF (SHIFT) THEN
c     ---calculate energy of the shift
         ECUT = 0
         CALL ENER(ECUT, vir, RC2, 1.0d0, 1.0d0)
         WRITE (10, 99005) RC, ECUT
      END IF
      IF (TAILCO) THEN
         WRITE (10, 99006) RC, CORU(RC, rho, 1.0d0, 1.0d0), CORP(RC, rho)
      END IF
      RETURN
99001 FORMAT ('  Number of equilibration cycles             :', i10, /, 
     &        '  Number of production cycles                :', i10, /, 
     &        '  Sample frequency                           :',i10, /)
99002 FORMAT ('  Number of att. to displ. a part. per cycle :', i10, /, 
     &        '  Maximum displacement                       :', f10.3, 
     &        //)
99003 FORMAT ('  Number of particles                        :', i10, /, 
     &        '  Temperature                                :', f10.3, 
     &        /, '  Density                                    :', 
     &        f10.3, /, '  Box length                                 :'
     &        , f10.3, /)
99004 FORMAT ('  Model parameters: ', /, '     epsilon: ', f5.3, /, 
     &        '     sigma  : ', f5.3, /, '     mass   : ', f5.3)
99005 FORMAT (' Simulations with TRUNCATED AND SHIFTED potential: ', /, 
     &        ' Potential truncated at :', f10.3, /, 
     &        ' Energy shif            :', f10.3, //)
99006 FORMAT (' Simulations with tail correction: ', /, 
     &        ' Potential truncated at  Rc =:', f10.3, /, 
     &        ' Tail corrections:   energy = ', f10.3, ' pressure ', 
     &        f10.3, /, /)
99007 FORMAT (' Requested density: ', f5.2, 
     &        ' different from density on disk: ', f5.2, /, 
     &        ' Rescaling of coordinates!',
     &         //)
99008 FORMAT (' number of lambda values: ', i10,
     &        //)
      END
