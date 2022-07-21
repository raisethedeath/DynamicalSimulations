C FILE: HamiltonianF.f90

        subroutine Hamiltonian(JMAX, DIMEN, ENERGY, ENDIM)
        implicit none
C
C       Populates the potential energy terms for the 
C       Hamiltonian matrix
C
        
        INTEGER :: JMAX, DIMEN, ENDIM

        INTEGER :: J, K, M
        INTEGER :: JJ, KK, MM
        INTEGER :: BRAMIN, BRAMAX
        INTEGER :: KETMIN, KETMAX

        INTEGER :: COL, ROW
        INTEGER :: ANGLE

        DOUBLE PRECISION :: BRANORM, KETNORM

        DOUBLE PRECISION :: PREHAM(DIMEN, DIMEN, 12)
        DOUBLE PRECISION :: TMPHAM(DIMEN, DIMEN)
        DOUBLE PRECISION :: HAMIL(DIMEN, DIMEN)
        DOUBLE PRECISION :: ENERGY(:,:)

        DOUBLE PRECISION :: ALPHA, BETA, GAMA       
        DOUBLE PRECISION :: WA, WB, WG, WEIGHT

        DOUBLE PRECISION :: BRALILD, KETLILD

        COMPLEX :: BRABIGD, KETBIGD
        COMPLEX :: BRA, KET
        COMPLEX :: H

        LOGICAL :: TIMECHECK
        DOUBLE PRECISION :: TIMEVAL1, TIMEVAL2, TIMEVAL
        DOUBLE PRECISION :: CORNER
        INTEGER :: HC

        DOUBLE PRECISION :: WIGNERD

        DOUBLE PRECISION :: PI
        PARAMETER (PI = 4 * ATAN(1.d0))

        WRITE (*,*)
        WRITE (*,*) "       Constructing pre-Hamiltonian arrays"
        WRITE (*,*)

        COL = 1

        DO J = 0, JMAX
          BRANORM = SQRT((2 * J + 1) / (8 * PI**2))
          DO K = -J, J
            DO M = -J, J
              BRAMIN = MAX(0, K-M)
              BRAMAX = MIN(J-M, J+K)

              PREHAM(COL, :, 1) = J
              PREHAM(COL, :, 2) = K
              PREHAM(COL, :, 3) = M
              PREHAM(COL, :, 4) = BRANORM
              PREHAM(COL, :, 5) = BRAMIN
              PREHAM(COL, :, 6) = BRAMAX

              PREHAM(:, COL, 7 ) = J
              PREHAM(:, COL, 8 ) = K
              PREHAM(:, COL, 9 ) = M
              PREHAM(:, COL, 10) = BRANORM
              PREHAM(:, COL, 11) = BRAMIN
              PREHAM(:, COL, 12) = BRAMAX


              COL = COL + 1

            ENDDO
          ENDDO
        ENDDO


        OPEN (unit=12, file='tmphamil.txt', status='OLD')
        READ (12,*) TMPHAM
        CLOSE (unit=12)


        HC = 1
        TIMECHECK = .FALSE.

        DO COL = 1, DIMEN
          J = INT(PREHAM(COL, 1, 1))
          K = INT(PREHAM(COL, 1, 2))
          M = INT(PREHAM(COL, 1, 3))
          BRANORM = PREHAM(COL, 1, 4)
          BRAMIN  = INT(PREHAM(COL, 1, 5))
          BRAMAX  = INT(PREHAM(COL, 1, 6))

          DO ROW = COL, DIMEN   
            JJ = INT(PREHAM(1, ROW, 7))
            KK = INT(PREHAM(1, ROW, 8))
            MM = INT(PREHAM(1, ROW, 9))
            KETNORM = PREHAM(1, ROW, 10)
            KETMIN  = INT(PREHAM(1, ROW, 11))
            KETMAX  = INT(PREHAM(1, ROW, 12))

            IF (ABS(TMPHAM(COL, ROW)) .GT. 0.d0) THEN
              HAMIL(COL, ROW) = TMPHAM(COL, ROW)
            ELSE

              CALL CPU_TIME(TIMEVAL1)

              H = 0.d0
              
              DO ANGLE = 1, ENDIM
                ALPHA = ENERGY(angle, 1)
                BETA  = ENERGY(angle, 2)
                GAMA  = ENERGY(angle, 3)

                WA = ENERGY(angle, 4)
                WB = ENERGY(angle, 5)
                WG = ENERGY(angle, 6)

                WEIGHT = WA * WB * WG

                BRALILD = WIGNERD(J, K, M, BRAMIN, BRAMAX, BETA)
                KETLILD = WIGNERD(JJ, KK, MM, KETMIN, KETMAX, BETA)

                BRABIGD = EXP(-(M*ALPHA + K*GAMA) * COMPLEX(0,1))
                BRABIGD = BRABIGD * BRALILD

                KETBIGD = EXP(-(MM*ALPHA + KK*GAMA) * COMPLEX(0,1))
                KETBIGD = KETBIGD * KETLILD

                BRA = CONJG(BRABIGD) * BRANORM
                KET = KETBIGD * KETNORM

                H = H + WEIGHT * ENERGY(ANGLE,7) * BRA * KET

              ENDDO

              CALL CPU_TIME(TIMEVAL2)

              IF (TIMECHECK .EQV. .FALSE.) THEN
                TIMEVAL = (TIMEVAL2-TIMEVAL1)
                TIMEVAL = TIMEVAL * DIMEN*DIMEN

                CORNER = (DIMEN*DIMEN) / 2. + DIMEN/2.
                TIMEVAL = TIMEVAL * (CORNER / (DIMEN*DIMEN))

                WRITE (*,*) "       Estimated time to populate"  , 
     1                   " Hamiltonian matrix - ", 
     2                    TIMEVAL, " seconds"
                WRITE (*,*)

                TIMECHECK = .TRUE.
                  
              ENDIF

              HAMIL(COL, ROW) = REAL(H)
              HAMIL(ROW, COL) = REAL(H)

            ENDIF

          ENDDO

        WRITE (*,*) "  ", (COL*DIMEN), " out of ", DIMEN*DIMEN, 
     1         "entries in the Rotational Hamiltonian Matrix populated"

        ENDDO

        OPEN (unit=13, file='hamil.tmp', status='replace')
        WRITE (13,*) HAMIL
        CLOSE (unit=13)


        END subroutine Hamiltonian


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

        DOUBLE PRECISION function WIGNERD(J, K, M, NMIN, NMAX, BETA)
C 
C       Uses the rotational quantum numbers to 
C       construct the Wigner-D rotational 
C       matrix elements following the work of Tajima
C       https://doi.org/10.1103/PhysRevC.91.014320
C
        DOUBLE PRECISION :: FACT
        INTEGER :: J, K, M
        INTEGER :: NMIN, NMAX
        DOUBLE PRECISION :: LILD, LILW, BIGW
        DOUBLE PRECISION :: BETA

        LILD = 0.d0

        DO N=NMIN, NMAX
            LILW = SQRT(FACT(J+M) * FACT(J-M) * FACT(J+K) * FACT(J-K))
     1                 / (FACT(J-M-N) * FACT(J+K-N) * FACT(N+M-K)
     2                 * FACT(N))

            BIGW = LILW * ((COS(BETA/2.d0))**(2.d0*J+K-M-2*N))
     1                  * ((-SIN(BETA/2.d0))**(M-K+2*N))

            LILD = LILD + ((-1)**N) * BIGW

        ENDDO

        WIGNERD = lilD

        END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

        DOUBLE PRECISION function fact(n)
C      
C       Calculates the factorial of n
C
        INTEGER :: n, p
        
        p = 1

        DO I=1,n
            p = p*I
        ENDDO

        fact = p

        END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C END FILE Hamiltonian_Fort.F
