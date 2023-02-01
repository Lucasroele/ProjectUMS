**==rantest.spg  processed by SPAG 4.52O  at 18:54 on 27 Mar 1996
 
      SUBROUTINE RANTEST(Iseed)
c
c     test and initialize the random number generator
c
      IMPLICIT NONE
      INTEGER Iseed, i
      DOUBLE PRECISION RANF
 
      CALL RANSET(Iseed)
      WRITE(10, *) ' ******** test random numbers ***********'
      DO i = 1, 5
         WRITE(10,*) ' i,ranf() ', i, RANF(Iseed)
      END DO
      RETURN
      END
