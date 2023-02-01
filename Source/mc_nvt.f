**==mc_nvt.spg  processed by SPAG 4.52O  at 18:54 on 27 Mar 1996
      PROGRAM MC_NVT
c________________________________________________________________________
c
c   Understanding Molecular Simulations: From Algorithms to Applications
c
c                 Daan Frenkel  and  Berend Smit
c
c  We make no warranties, express or implied, that the programs contained
c  in this work are free of error, or that they will meet your requirements
c  for any particular application. They should not be relied on for solving
c  problems whose incorrect solution could results in injury, damage, or
c  loss of property. The authors and publishers disclaim all liability for
c  direct or consequential damages resulting from your use of the programs
c
c__________________________________________________________________________
c
c
c   Case Study 1: Equation of state of the Lennard-Jones fluid
c
c__________________________________________________________________________
 
      IMPLICIT NONE
      INCLUDE 'system.inc'
      INCLUDE 'parameter.inc'
      INCLUDE 'conf.inc'
      INCLUDE 'potential.inc'
      INTEGER iseed, equil, prod, nsamp, ii, icycl, ndispl, attempt, 
     &        nacc, ncycl, nmoves, imove, nLambda, I, J, SampleCount,
     &        nGhosts, nWidomCycle, runWidom, runTDI
      DOUBLE PRECISION en, ent, vir, virt, dr, Lambda , Press, PressSum, 
     &        ChemicalPotentialSum, RANF, widomPress,
     &        Xi, Yi, Zi, EnSum, EnSquaredSum, rho, CORU,
     &        EnDummy, VirDummy, EnStart, chempotId, chempotEx, chempotL
      WRITE (10, *) '**************** MC_NVT ***************'
c     ---initialize system
      CALL READDAT(equil, prod, nsamp, ndispl, dr, iseed, nLambda, nGhosts, nWidomCycle, runWidom,
     & runTDI)
      nmoves = ndispl
c     --- initializing system Widom     
      PressSum = 0.0d0
      SampleCount = 0
      EnSum = 0.0d0
      EnSquaredSum = 0.0d0
      ChemicalPotentialSum = 0.0d0
      Lambda = 1.0d0
      EnStart = 0.0d0
c     ---total energy of the system
      CALL TOTERG(en, vir, Lambda)
      WRITE (10, 99001) en, vir

c     ---start MC-cycle
      WRITE (10, *) en, vir

c     --- Widom
      IF (runWidom.ne.0) THEN
      DO ii = 1, 2

c        --- ii=1 equilibration
c        --- ii=2 production
         IF (ii.EQ.1) THEN
            ncycl = equil
            IF (ncycl.NE.0) WRITE (10, *) ' Start equilibration '
         ELSE
            IF (ncycl.NE.0) WRITE (10, *) ' Start production '
            ncycl = nWidomCycle
         END IF

         attempt = 0
         nacc = 0

c        Mathmatically describe ADJUST
c        Mathmatically describe ADJUST
c        Mathmatically describe ADJUST
c        Mathmatically describe ADJUST
c        Mathmatically describe ADJUST

c        ---intialize the subroutine that adjusts the maximum displacement
         CALL ADJUST(attempt, nacc, dr)

         DO icycl = 1, ncycl
            DO imove = 1, nmoves
c              ---attempt to displace a particle
               CALL MCMOVE(en, vir, attempt, nacc, dr, iseed, Lambda)
            END DO

c           --- sample the system every nsamp times
c           --- assumes decorellation after nsamp steps
            IF (ii.EQ.2) THEN
               IF (MOD(icycl,nsamp).EQ.0) Then
                  CALL SAMPLE(icycl, en, vir, press, Lambda)
                  SampleCount = SampleCount + 1
                  PressSum = PressSum + press
                  EnSquaredSum = EnSquaredSum + en * en
                  EnSum = EnSum + en
c                 --- Calculation of chemical potential with nGhosts trial chains
                  DO J = 1, nGhosts
                     Xi = BOX * RANF()
                     Yi = BOX * RANF()
                     Zi = BOX * RANF()
c                    --- Taking lambda as 1
                     CALL ENERI(Xi, Yi, Zi, 0, 1, EnDummy, VirDummy, Lambda, Lambda)
                     IF (TAILCO) THEN
c                       --- verify tail correction
                        rho = NPART/(BOX**3)
                        EnDummy = EnDummy + 2*CORU(RC, rho, Lambda, Lambda)
                     END IF
                     ChemicalPotentialSum = ChemicalPotentialSum + Exp(-BETA * EnDummy)
                  END DO
               END IF
            END IF

c           --- Runs every 1/5th of ncycl (prod and equil)
            IF (MOD(icycl,ncycl/5).EQ.0) THEN
               WRITE (10, *) '======>> Done ', icycl, ' out of ', ncycl
c              ---write intermediate configuration to file
c              CALL STORE(8, dr)
c              ---adjust maximum displacements
               CALL ADJUST(attempt, nacc, dr)
            END IF
         END DO
         IF (ncycl.NE.0) THEN
            IF (attempt.NE.0) WRITE (10, 99003) attempt, nacc,
     &                               100.*FLOAT(nacc)/FLOAT(attempt)
c           ---test total energy
            CALL TOTERG(ent, virt, Lambda)

            IF (ABS(ent-en).GT.1.D-6) THEN
                WRITE (10, *)
     &         ' ######### PROBLEMS ENERGY ################ '
            END IF
            IF (ABS(virt-vir).GT.1.D-6) THEN
               WRITE (10, *)
     &        ' ######### PROBLEMS VIRIAL ################ '
            END IF
            WRITE (10, 99002) ent, en, ent - en, virt, vir, virt - vir


            IF(ii .Eq. 2) THEN
               widomPress = (PressSum/DBLE(SampleCount))/DBLE(nGhosts)
               chempotEx = -Log((ChemicalPotentialSum/DBLE(SampleCount))
     &         / DBLE(nGhosts))/BETA
               chempotId = -log(BOX*BOX*BOX/DBLE(npart+1))/BETA
               Write(10,99004) widomPress, nGhosts*SampleCount, chempotEx
               Write(10,99004) widomPress, nGhosts*SampleCount, chempotEx
               Write(40,*) rho, chempotEx
               Write(44, 99005) chempotEx, chempotId
            END IF
         END IF
      END DO

      WRITE (10, *)
     &        ' ################################################ '
      WRITE (10, *)
     &        ' ################################################ '
      WRITE (10, *)
     &        ' ################ WIDOM FINISHED ################ '
      WRITE (10, *)
     &        ' ################################################ '
      WRITE (10, *)
     &        ' ################################################ '

      END IF
      

      IF (runTDI.NE.0) THEN
c     --- TDI
      DO I = 0, 2 * nlambda
         IF (MOD(I,nlambda/10).eq.0) write(6,*)'Running lambda step', I+1, 'out of', 2*nlambda/10
         IF (I .EQ. nlambda + 1) THEN
            chempotL = (en - EnStart) / DBLE(NPART)
            write (44,99006) Lambda, en - Enstart, chempotL
         END IF
         IF (I .LE. nlambda) THEN
            Lambda =  1.0d0 - (DBLE(I) / DBLE(nLambda))
         ELSE
            Lambda = (DBLE(I) / DBLE(nLambda)) - 1.0d0
         END IF
         WRITE (10,*) 'lambda = ', I
         CALL TOTERG(en, vir, Lambda)
            DO ii = 1, 2
c           --- ii=1 equilibration
c           --- ii=2 production
            IF (ii.EQ.1) THEN
               ncycl = equil
               IF (ncycl.NE.0) WRITE (10, *) ' Start equilibration '
            ELSE
               IF (ncycl.NE.0) WRITE (10, *) ' Start production '
               ncycl = prod
            END IF

            attempt = 0
            nacc = 0

c           ---intialize the subroutine that adjusts the maximum displacement
            CALL ADJUST(attempt, nacc, dr)

            DO icycl = 1, ncycl

               DO imove = 1, nmoves
c                 ---attempt to displace a particle
                  CALL MCMOVE(en, vir, attempt, nacc, dr, iseed, Lambda)
               END DO

c              --- sample the system every nsamp times
c              --- assumes decorellation after nsamp steps
c               IF (ii.EQ.2) THEN
c                  IF (MOD(icycl,nsamp).EQ.0) Then
c                     CALL SAMPLE(icycl, en, vir, press, lambda)
c                  END IF
c               END IF

c              --- Runs every 1/5th of ncycl (prod and equil)
               IF (MOD(icycl,ncycl/5).EQ.0) THEN
                  WRITE (10, *) '======>> Done ', icycl, ' out of ', ncycl
c                 ---write intermediate configuration to file
c                  CALL STORE(8, dr)
c                 ---adjust maximum displacements
                  CALL ADJUST(attempt, nacc, dr)
               END IF
            END DO


            IF (I .eq. 0) THEN
               IF (ii .eq. 1) THEN
                  EnStart = en
               END IF
            END IF

            IF (ncycl.NE.0) THEN
               IF (attempt.NE.0) WRITE (10, 99003) attempt, nacc,
     &                               100.*FLOAT(nacc)/FLOAT(attempt)

c           ---test total energy
               CALL TOTERG(ent, virt, Lambda)
               IF (II .EQ. 2) WRITE(44,*) lambda, ent 
               IF (ABS(ent-en).GT.1.D-6) THEN
                  WRITE (10, *)
     &                    ' ######### PROBLEMS ENERGY ################ '
               END IF
               IF (ABS(virt-vir).GT.1.D-6) THEN
                  WRITE (10, *)
     &                    ' ######### PROBLEMS VIRIAL ################ '
               END IF
               WRITE (10, 99002) ent, en, ent - en, virt, vir, virt - vir
            END IF
         END DO
      END DO
      END IF
      chempotL = (EnStart - en) / DBLE(NPART) + chempotL
      write (44, 99006) Lambda, EnStart + DBLE(NPART) * chempotL - en, chempotL
      CALL STORE(21, dr)

      WRITE (10, *)
     &        ' ################################################ '
      WRITE (10, *)
     &        ' ################################################ '
      WRITE (10, *)
     &        ' ################ LAMBDA FINISHED ############### '
      WRITE (10, *)
     &        ' ################################################ '
      WRITE (10, *)
     &        ' ################################################ '
      
      STOP
 
99001 FORMAT (' Total energy initial configuration: ', f12.5, /, 
     &        ' Total virial initial configuration: ', f12.5)
99002 FORMAT (' Total energy end of simulation    : ', f12.5, /, 
     &        '       running energy              : ', f12.5, /, 
     &        '       difference                  :  ', e12.5, /, 
     &        ' Total virial end of simulation    : ', f12.5, /, 
     &        '       running virial              : ', f12.5, /, 
     &        '       difference                  :  ', e12.5)
99003 FORMAT (' Number of att. to displ. a part.  : ', i10, /, 
     &        ' success: ', i10, '(= ', f5.2, '%)')
99004 FORMAT (' Results Widom test particle method : '/,
     &        ' Average pressure          :', f12.5, /,
     &        ' Number of samples         :', i12, /,
     &        ' Excess chemical potential :', f12.5)
99005 FORMAT ('Finished a WIDOM cycle', /,
     &        'Ex Chemical potential      :', f12.5, /,
     &        'Id Chemical potential      :', f12.5)
99006 FORMAT ('Finished a LAMBDA cycle ending at lambda =', f8.6, /,
     &        'Free energy                :', f12.5, /,
     &        'Chemical potential         :', f12.5)
      END
