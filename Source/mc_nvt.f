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
      INTEGER iseed, equil, prod, nsamp, ii, icycl, ndispl, attempt, 
     &        nacc, ncycl, nmoves, imove, nLambda, I, J
      DOUBLE PRECISION en, ent, vir, virt, dr, Lambda , Press, PressSum, 
     &        PressCount, ChemicalPotentialSum, ChemicalPotentialCount, 
     &        Xi, Yi, Zi, EnSum, EnSquaredSum, EnCount, RandomNumber,
     &        EnDummy, VirDummy
 
      WRITE (6, *) '**************** MC_NVT ***************'
c     ---initialize system
      CALL READDAT(equil, prod, nsamp, ndispl, dr, iseed)
      nmoves = ndispl
c     --- initializing system Widom     
      PressSum = 0.0d0
      PressCount = 0.0d0
      EnSum = 0.0
      EnSquaredSum = 0.0
      ChemicalPotentialSum = 0.0
      ChemicalPotentialCount = 0.0d0

c     ---total energy of the system
      Lambda = 1
      nLambda = 25000
      CALL TOTERG(en, vir, Lambda)
c      CALL TOTERG(en, vir)
      WRITE (6, 99001) en, vir
c     ---start MC-cycle
      WRITE (6, *) en, vir
      DO ii = 1, 2
c        --- ii=1 equilibration
c        --- ii=2 production
         IF (ii.EQ.1) THEN
            ncycl = equil
            IF (ncycl.NE.0) WRITE (6, *) ' Start equilibration '
         ELSE
            IF (ncycl.NE.0) WRITE (6, *) ' Start production '
            ncycl = prod
         END IF
         attempt = 0
         nacc = 0
c        ---intialize the subroutine that adjust the maximum displacement
         CALL ADJUST(attempt, nacc, dr)


         DO I = nLambda, 0
            If (I==0) then
               Lambda = 0
            Else
               Lambda = DBLE(nLambda) / DBLE(I)
               WRITE (66, *) "Lambda = ", Lambda   
            Endif
             DO icycl = 1, ncycl
                DO imove = 1, nmoves
c                  ---attempt to displace a particle
                   CALL MCMOVE(en, vir, attempt, nacc, dr, iseed, Lambda)
                END DO
                IF (ii.EQ.2) THEN
c                  ---sample averages
                   IF (MOD(icycl,nsamp).EQ.0) Then
                     CALL SAMPLE(icycl, en, vir, press, lambda)
                     PressSum = PressSum + Press
                     PressCount = PressCount + 1.0d0

                     EnSquaredSum = EnSquaredSum + en * en
                     EnSum = EnSum + en 
                     EnCount = EnCount + 1.0d0

c              --- Calculation of chemical potential with 10 trial chains                     
                     DO J = 1,10
                        Xi = Box * RandomNumber()
                        Yi = Box * RandomNumber()
                        Zi = Box * RandomNumber()
c                    --- Taking lambda as 1?
                        CALL ENERI(Xi, Yi, Zi, npart +1, 1, EnDummy, VirDummy,
     &                  1)
                        ChemicalPotentialCount = ChemicalPotentialCount + 1.0d0
                        ChemicalPotentialSum = ChemicalPotentialSum + Dble(Exp(-Beta * EnDummy))
                     END DO
                  END IF
                END IF

                IF (MOD(icycl,ncycl/5).EQ.0) THEN
                   WRITE (6, *) '======>> Done ', icycl, ' out of ', ncycl
c                  ---write intermediate configuration to file
                   CALL STORE(8, dr)
c                  ---adjust maximum displacements
                   CALL ADJUST(attempt, nacc, dr)
                END IF
             END DO
             WRITE (66, *)
         END DO



         IF (ncycl.NE.0) THEN
            IF (attempt.NE.0) WRITE (6, 99003) attempt, nacc, 
     &                               100.*FLOAT(nacc)/FLOAT(attempt)
c           ---test total energy
            CALL TOTERG(ent, virt, Lambda)
c            CALL TOTERG(ent, virt)

            IF (ABS(ent-en).GT.1.D-6) THEN
               WRITE (6, *) 
     &                    ' ######### PROBLEMS ENERGY ################ '
            END IF
            IF (ABS(virt-vir).GT.1.D-6) THEN
               WRITE (6, *) 
     &                    ' ######### PROBLEMS VIRIAL ################ '
            END IF
            WRITE (6, 99002) ent, en, ent - en, virt, vir, virt - vir

c           --- Print Chemical Potential and Pressure
            IF(ii .Eq. 2) THEN
               Write(10,*) 'Heat Capacity                      :',0.0
               Write(10,*) 'Average Pressure                  : ',
     &           PressSum/PressCount
               Write(10,*) 'Chemical Potential                : ',
     &           -Log((ChemicalPotentialSum/ChemicalPotentialCount)*
     &           (Box*Box*Box/Dble(npart)))/Beta
            END IF
         END IF
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
