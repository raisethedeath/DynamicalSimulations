C FILE: InteractionEnergyF.f90

        subroutine InteractionEnergy(ENERGY, H2POS, H2POLY, H2ABSV, 
     1             H2ABSW, H2SMEAR, H2SIG, HPOS, HPOLY, HABSV, 
     2             HABSW, HSMEAR, HSIG, THRESH)
        implicit none
C
C       Populates the potential energy terms for the 
C       Hamiltonian matrix
C

        interface
          function Rotate(A,B,G) RESULT(MAT)
            DOUBLE PRECISION :: A, B, G
            DOUBLE PRECISION :: MAT(3,3)
          end function
          function IntEn(POS, FIT, ABSV, ABSW,
     1                   QUAD, SIG, THRESH) RESULT(E)
            DOUBLE PRECISION :: POS(3), FIT(:,:), ABSV(:), ABSW(:)
            DOUBLE PRECISION :: SIG, THRESH, E
            INTEGER :: QUAD
          end function
        end interface

        DOUBLE PRECISION :: ENERGY(:,:)

        DOUBLE PRECISION :: H2POS(:,:), H2POLY(:,:)
        DOUBLE PRECISION :: H2ABSV(:), H2ABSW(:)
        INTEGER :: H2SMEAR
        DOUBLE PRECISION :: H2SIG

        DOUBLE PRECISION :: HPOS(:,:), HPOLY(:,:)
        DOUBLE PRECISION :: HABSV(:), HABSW(:)
        INTEGER :: HSMEAR
        DOUBLE PRECISION :: HSIG

        DOUBLE PRECISION :: THRESH

        DOUBLE PRECISION :: ALPHA, BETA, GAMA
        INTEGER :: NUMH2(2), NUMH(2), NUMANG(2), ANGLE, SITE, TOTMOL
        INTEGER :: ANGC, ANGCC

        DOUBLE PRECISION :: E, ROTMAT(3,3), NEWPOS(3)
        
        DOUBLE PRECISION :: START, FINISH, TIMEVAL

        NUMH2  = SHAPE(H2POS)
        NUMH   = SHAPE(HPOS)
        NUMANG = SHAPE(ENERGY)

        ANGC = NUMANG(1) / 10
        ANGCC = 1

        IF (NUMH(2) .NE. 1) THEN
          TOTMOL = NUMH2(1) + NUMH(1)
        ELSE
          TOTMOL = NUMH2(1)
        ENDIF

        WRITE(*,*)

        DO ANGLE=1,NUMANG(1)
          ALPHA = ENERGY(ANGLE,1)
          BETA  = ENERGY(ANGLE,2)
          GAMA  = ENERGY(ANGLE,3)

          E = 0.d0

          ROTMAT = Rotate(ALPHA, BETA, GAMA)

          IF (ANGLE .EQ. 1) THEN
            CALL CPU_TIME(START)
          ENDIF

          DO SITE=1,NUMH2(1)
            NEWPOS = MATMUL(H2POS(SITE,:), ROTMAT)

            E = E + IntEn(NEWPOS, H2POLY, H2ABSV, H2ABSW,
     1                    H2SMEAR, H2SIG, THRESH)

            WRITE(*,*) E

            CALL EXIT(0)

            IF (ANGLE .EQ. 1) THEN
              IF (SITE .EQ. 1) THEN
                CALL CPU_TIME(FINISH)
                TIMEVAL = FINISH - START
                TIMEVAL = TIMEVAL * NUMANG(1) * TOTMOL

                WRITE(*,*) "       Estimated time for calculating"  ,
     1                     " all interaction energies - ",
     2                      TIMEVAL, " seconds"
                WRITE (*,*)

              ENDIF
            ENDIF

          ENDDO
        
          IF (NUMH(2) .NE. 1) THEN
            DO SITE=1,NUMH(1)
              NEWPOS = MATMUL(HPOS(SITE,:), ROTMAT)

              E = E + IntEn(NEWPOS, HPOLY, HABSV, HABSW,
     1                      HSMEAR, HSIG, THRESH)
            ENDDO
          ENDIF

          ENERGY(ANGLE,7) = E

          if (ANGLE .GT. ANGC*ANGCC) THEN
            WRITE(*,*) "  ", (ANGLE*TOTMOL), " out of ",
     1                       (TOTMOL*NUMANG(1)),
     2         "interaction energies calculated"
            ANGCC = ANGCC + 1
          ENDIF

        ENDDO

        OPEN (unit=13, file='E.en', status='replace')
        WRITE (13,*) ENERGY(:,7)
        CLOSE (unit=13)

        

        END subroutine InteractionEnergy


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

        DOUBLE PRECISION function IntEn(POS, FIT, ABSV, ABSW,
     1          QUAD, SIG, THRESH) RESULT(E)

C
C       Calculate interaction energy
C

        interface 
          function Lebedev(COEF,T,P,THRESH) RESULT(E)
            DOUBLE PRECISION :: COEF(81), T, P, THRESH, E
          end function
          function LebCoef(POLY,R) RESULT(C)
            DOUBLE PRECISION :: POLY(:), R, C
          end function
        end interface

        DOUBLE PRECISION :: POS(3), FIT(:,:), ABSV(:), ABSW(:)
        DOUBLE PRECISION :: SIG, THRESH
        INTEGER :: QUAD

        DOUBLE PRECISION :: WX, WY, WZ
        INTEGER :: JX, JY, JZ
        DOUBLE PRECISION :: NEWPOS(3), NEWPOSr(3), R, T, P
        DOUBLE PRECISION :: COEF(81)
        INTEGER :: C

        DOUBLE PRECISION :: PI
        PARAMETER (PI = 4 * ATAN(1.d0))

        E = 0.d0


        DO JX=1,QUAD
          NEWPOS(1) = SIG * ABSV(JX) + POS(1)
          WX = ABSW(JX)
          DO JY=1,QUAD
            NEWPOS(2) = SIG * ABSV(JY) + POS(2)
            WY = ABSW(JY)
            DO JZ=1,QUAD
              NEWPOS(3) = SIG * ABSV(JZ) + POS(3)
              WZ = ABSW(JZ)

              R = NORM2(NEWPOS)

              DO C=1,81
                COEF(C) = LebCoef(FIT(C,:), R)
              ENDDO

              T = ACOS(NEWPOSr(3)/R)
              P = ATAN(NEWPOSr(2), NEWPOSr(1))

              IF (P .LT. 0.d0) THEN
                P = P + 2.d0*PI
              ENDIF

              E = E + Lebedev(COEF,T,P,THRESH) * WX*WY*WZ

              WRITE(*,*) E / (WX*WY*WZ)
              CALL EXIT(0)

            ENDDO
          ENDDO
        ENDDO

        E = E / SUM(ABSW)**3

        End function

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

        DOUBLE PRECISION function LebCoef(POLY,R) RESULT(C)

C
C       Calculate the interpolated lebedev coefficients
C

        DOUBLE PRECISION :: POLY(:), R
        INTEGER :: N(1), X

        N = SHAPE(POLY)

        C = 0.d0

        DO X=0,N(1)-1
          C = C + POLY(N(1)-X) * R**X
        ENDDO

        END function


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

        DOUBLE PRECISION function Lebedev(COEF,T,P,THRESH) RESULT(E)

C
C       Calculate interaction energy using Lebedev quadrature
C

        interface 
          function SphereHarm(J,M,T,P) RESULT(VAL)
            DOUBLE PRECISION :: T, P, Val
            INTEGER :: J, M
          end function
        end interface

        DOUBLE PRECISION :: COEF(81), T, P, THRESH
        DOUBLE PRECISION :: SH
        INTEGER :: C, J, M

        DOUBLE PRECISION :: PI
        PARAMETER (PI = 4 * ATAN(1.d0))

        C = 1
        E = 0.d0

        DO J=0,8
          DO M=-J,J
            IF (ABS(COEF(C)) .GT. THRESH) THEN
              SH = SphereHarm(J,M,T,P)
              E = E + SH * COEF(C)
              WRITE(*,*) J, M, C, SH, COEF(C), E
            ENDIF
            C = C + 1
          ENDDO
        ENDDO

        END function

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

        DOUBLE PRECISION function SphereHarm(J,M,T,P) RESULT(VAL)

C      
C       Compute Spherical Harmonics
C

        interface
          function Fac(N) RESULT(P)
            INTEGER :: N
            DOUBLE PRECISION :: P
          end function
        end interface

        DOUBLE PRECISION :: T, P
        INTEGER :: J, M, K, N, KK
        PARAMETER (K=1, KK=3)

        DOUBLE PRECISION :: PREFAC
        COMPLEX :: EX, SPHERE
        DOUBLE PRECISION :: LEG(0:J,0:J), LEGN1(0:J), LEGN2(0:J)
        DOUBLE PRECISION :: DLEG(0:J), ALEG, TEMPVAL
        INTEGER :: D, DD

        DOUBLE PRECISION :: PI
        PARAMETER (PI = 4 * ATAN(1.d0))




        ! Coefficients of Legendre Polynomials

        LEG = 0

        IF (J .EQ. 0) THEN
         LEG(0,0) = 1
        ELSE IF (J .EQ. 1) THEN
         LEG(0,0) = 1
         LEG(1,1) = 1
        ELSE
           LEG(0,0) = 1
           LEG(1,1) = 1
          DO N=2,J
            LEGN2 = LEG(N-2,:) * (N-1)
            
            LEGN1 = 0
            LEGN1(1:J) = LEG(N-1,0:J-1) * (2 * N - 1)

            LEG(N,:) = (LEGN1 - LEGN2) / N
          ENDDO
        ENDIF

        ! Mth Derivative of Legendre Polynomials

        DLEG = LEG(J,:)

        DO D=1,ABS(M)
          DLEG(0:J-D) = DLEG(1:J)
          DLEG(J-D+1) = 0

          DO DD=0,J-D
            DLEG(DD) = DLEG(DD) * (DD+1)
          ENDDO
        ENDDO

        ! Value of Associate Legendre Polynomial at COS(T)

        ALEG = (-1)**ABS(M) * (1 - COS(T)**2)**(ABS(M)/2.)

        TEMPVAL = 0.d0

        DO D=0,J
          TEMPVAL = TEMPVAL + DLEG(D) * COS(T)**D
        ENDDO

        ALEG = ALEG * TEMPVAL

        ! Change for M < 0

        IF (M .LT. 0) THEN
          ALEG = ALEG * (-1)**M * Fac(J-ABS(M)) / Fac(J+ABS(M))
        ENDIF

        ! Prefactor for spherical harmonics

        PREFAC = SQRT(((2.d0*J+1.d0) * Fac(J-M)) /
     1                ((4.d0*PI) * Fac(J+M)))

        ! Complex eponential for spherical harmonic

        EX = EXP(COMPLEX(0,1) * M * P)

        ! Value of spherical harmonic

        SPHERE = PREFAC * EX * ALEG

        VAL = REAL(SPHERE)

        END function

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

        DOUBLE PRECISION function Fac(N) RESULT(P)
C      
C       Calculates the factorial of n
C
        INTEGER :: N, I

        P = 1

        DO I=1,N
            P = P*I
        ENDDO

        END function

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

        DOUBLE PRECISION function Rotate(A, B, G) RESULT(MAT)

C
C       Constructs a rotation matrix using Euler angles
C

        DOUBLE PRECISION :: A, B, G
        DIMENSION :: MAT(3,3)

        DOUBLE PRECISION :: COSA, COSB, COSG
        DOUBLE PRECISION :: SINA, SINB, SING

        COSA = COS(A)
        COSB = COS(B)
        COSG = COS(G)

        SINA = SIN(A)
        SINB = SIN(B)
        SING = SIN(G)

        
        MAT(1,1) =  COSA * COSB * COSG - SINA * SING
        MAT(1,2) = -COSA * COSB * SING - SINA * COSG
        MAT(1,3) =  COSA * SINB

        MAT(2,1) =  SINA * COSB * COSG + COSA * SING
        MAT(2,2) = -SINA * COSB * SING + COSA * COSG
        MAT(2,3) =  SINA * SINB

        MAT(3,1) = -SINB * COSG
        MAT(3,2) =  SINB * SING
        MAT(3,3) =  COSB

        END function

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C END FILE InteractionEnergyF.f90
