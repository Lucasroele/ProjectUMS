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
      INTEGER iseed, equil, prod, nsamp, ii, icycl, ndispl, attempt, 
     &        nacc, ncycl, nmoves, imove, nLambda, I, J, SampleCount,
     &        nWidom, nWidomCycle
      DOUBLE PRECISION en, ent, vir, virt, dr, Lambda , Press, PressSum, 
     &        ChemicalPotentialSum, RANF,
     &        Xi, Yi, Zi, EnSum, EnSquaredSum,
     &        EnDummy, VirDummy
      WRITE (6, *) '**************** MC_NVT ***************'
c     ---initialize system
      CALL READDAT(equil, prod, nsamp, ndispl, dr, iseed, nLambda, nWidom, nWidomCycle)
      nmoves = ndispl
c     --- initializing system Widom     
      PressSum = 0.0d0
      SampleCount = 0
      EnSum = 0.0d0
      EnSquaredSum = 0.0d0
      ChemicalPotentialSum = 0.0d0

c     ---total energy of the system
      Lambda = 1

      CALL TOTERG(en, vir, Lambda)
      WRITE (6, 99001) en, vir
c     ---start MC-cycle
      WRITE (6, *) en, vir

c     --- Widom
      DO ii = 1, 2

c        --- ii=1 equilibration
c        --- ii=2 production
         IF (ii.EQ.1) THEN
            ncycl = equil
            IF (ncycl.NE.0) WRITE (6, *) ' Start equilibration '
         ELSE
            IF (ncycl.NE.0) WRITE (6, *) ' Start production '
            ncycl = nWidomCycle
         END IF

         attempt = 0
         nacc = 0

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
                  CALL SAMPLE(icycl, en, vir, press, lambda)
                  SampleCount = SampleCount + 1
                  PressSum = PressSum + Press
                  EnSquaredSum = EnSquaredSum + en * en
                  EnSum = EnSum + en
c                 --- Calculation of chemical potential with 10 trial chains
                  DO J = 1, nWidom
                     Xi = BOX * RANF()
                     Yi = BOX * RANF()
                     Zi = BOX * RANF()
c                    --- Taking lambda as 1
c                    --- ENERI needs checking
                     CALL ENERI(Xi, Yi, Zi, 0, NPART+1, EnDummy, VirDummy, Lambda)
                     ChemicalPotentialSum = ChemicalPotentialSum + Exp(-BETA * EnDummy)
                  END DO
               END IF
            END IF

            IF (MOD(icycl,ncycl/5).EQ.0) THEN
               WRITE (6, *) '======>> Done ', icycl, ' out of ', ncycl
c              ---write intermediate configuration to file
               CALL STORE(8, dr)
c              ---adjust maximum displacements
               CALL ADJUST(attempt, nacc, dr)
            END IF
            WRITE (66, *)
         END DO

         IF (ncycl.NE.0) THEN
            IF (attempt.NE.0) WRITE (6, 99003) attempt, nacc,
     &                               100.*FLOAT(nacc)/FLOAT(attempt)
c           ---test total energy
            CALL TOTERG(ent, virt, Lambda)

            IF (ABS(ent-en).GT.1.D-6) THEN
                WRITE (6, *)
     &         ' ######### PROBLEMS ENERGY ################ '
            END IF
            IF (ABS(virt-vir).GT.1.D-6) THEN
               WRITE (6, *)
     &        ' ######### PROBLEMS VIRIAL ################ '
            END IF
            WRITE (6, 99002) ent, en, ent - en, virt, vir, virt - vir

c           --- Print Chemical Potential and Pressure
            IF(ii .Eq. 2) THEN
               Write(10,*) 'Heat Capacity                      :',0.0
               Write(10,*) 'Average Pressure                  : ',
     &         PressSum/DBLE(SampleCount)
               Write(10,*) 'Chemical Potential                : ',
     &         -Log(((ChemicalPotentialSum / DBLE(SampleCount)) / DBLE(nWidom))*
     &         (BOX*BOX*BOX/DBLE(npart)))/BETA
            END IF
         END IF
      END DO
      WRITE (6, *)
     &        ' ################################################ '
      WRITE (6, *)
     &        ' ################################################ '
      WRITE (6, *)
     &        ' ################ WIDOM FINISHED ################ '
      WRITE (6, *)
     &        ' ################################################ '
      WRITE (6, *)
     &        ' ################################################ '

c     --- TDI
      DO I = nLambda, 0, -1
         Lambda =  DBLE(I) / DBLE(nLambda)
         write (6,*) I
         DO ii = 1, 2

c           --- ii=1 equilibration
c           --- ii=2 production
            IF (ii.EQ.1) THEN
               ncycl = equil
               IF (ncycl.NE.0) WRITE (6, *) ' Start equilibration '
            ELSE
               IF (ncycl.NE.0) WRITE (6, *) ' Start production '
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
               IF (ii.EQ.2) THEN
                  IF (MOD(icycl,nsamp).EQ.0) Then
                     CALL SAMPLE(icycl, en, vir, press, lambda)
                  END IF
               END IF

               IF (MOD(icycl,ncycl/5).EQ.0) THEN
                  WRITE (6, *) '======>> Done ', icycl, ' out of ', ncycl
c                 ---write intermediate configuration to file
                  CALL STORE(8, dr)
c                 ---adjust maximum displacements
                  CALL ADJUST(attempt, nacc, dr)
               END IF
               WRITE (66, *)
            END DO


            IF (ncycl.NE.0) THEN
               IF (attempt.NE.0) WRITE (6, 99003) attempt, nacc,
     &                               100.*FLOAT(nacc)/FLOAT(attempt)
c           ---test total energy
               CALL TOTERG(ent, virt, Lambda)
               IF (ABS(ent-en).GT.1.D-6) THEN
                  WRITE (6, *)
     &                    ' ######### PROBLEMS ENERGY ################ '
               END IF
               IF (ABS(virt-vir).GT.1.D-6) THEN
                  WRITE (6, *)
     &                    ' ######### PROBLEMS VIRIAL ################ '
               END IF
               WRITE (6, 99002) ent, en, ent - en, virt, vir, virt - vir
            END IF
         END DO
      END DO
      CALL STORE(21, dr)
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
      END
      
