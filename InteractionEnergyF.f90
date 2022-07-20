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
     1                   QUAD, SIG, THRESH) RESULT(V)
            DOUBLE PRECISION :: POS(3), FIT(:,:), ABSV(:), ABSW(:)
            DOUBLE PRECISION :: SIG, THRESH, V
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
        INTEGER :: NUMH2(2), NUMH(2), NUMANG(2), ANGLE, SITE

        DOUBLE PRECISION :: E, ROTMAT(3,3), NEWPOS(3)

        NUMH2  = SHAPE(H2POS)
        NUMH   = SHAPE(HPOS)
        NUMANG = SHAPE(ENERGY)

        DO ANGLE=1,NUMANG(1)
          ALPHA = ENERGY(ANGLE,1)
          BETA  = ENERGY(ANGLE,2)
          GAMA  = ENERGY(ANGLE,3)

          E = 0.d0

          ROTMAT = Rotate(ALPHA, BETA, GAMA)

          DO SITE=1,NUMH2(1)
            NEWPOS = MATMUL(H2POS(SITE,:), ROTMAT)

            E = E + IntEn(NEWPOS, H2POLY, H2ABSV, H2ABSW,
     1                   H2SMEAR, H2SIG, THRESH)

            WRITE(*,*) E

            CALL EXIT(0)


          ENDDO



          



          CALL EXIT(0)

        

        ENDDO



        

        END subroutine InteractionEnergy


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

        DOUBLE PRECISION function IntEn(POS, FIT, ABSV, ABSW,
     1          QUAD, SIG, THRESH) RESULT(V)

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

        DOUBLE PRECISION :: WX, WY, WZ, W
        INTEGER :: JX, JY, JZ
        DOUBLE PRECISION :: NEWPOS(3), R, T, P
        DOUBLE PRECISION :: COEF(81)
        INTEGER :: C

        DOUBLE PRECISION :: PI
        PARAMETER (PI = 4 * ATAN(1.d0))

        V = 0.d0


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

              NEWPOS = NEWPOS / R

              DO C=1,81
                COEF(C) = LebCoef(FIT(C,:), R)
              ENDDO

              T = ACOS(NEWPOS(3)/R)
              P = ATAN(NEWPOS(2), NEWPOS(1))

              IF (P .LT. 0.d0) THEN
                P = P + 2.d0*PI
              ENDIF

              E = E + Lebedev(COEF,T,P,THRESH) * WX*WY*WZ

              WRITE (*,*) E
              

              CALL EXIT(0)

            ENDDO
          ENDDO
        ENDDO


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
          function Fac(N) RESULT(P)
            INTEGER :: N
            DOUBLE PRECISION :: P
          end function
          function Legendre(K,X) RESULT(VAL)
            DOUBLE PRECISION :: X
            INTEGER :: K
            DOUBLE PRECISION :: VAL
          end function
        end interface

        DOUBLE PRECISION :: COEF(81), T, P, THRESH
        INTEGER :: C, J, M
        DOUBLE PRECISION :: EX, PREFAC, LEG

        DOUBLE PRECISION :: PI
        PARAMETER (PI = 4 * ATAN(1.d0))

        C = 1

        DO J=1,8
          DO M=-J,J
            IF (ABS(COEF(C)) .GT. THRESH) THEN
              PREFAC = SQRT(((2.d0*J+1.d0)/(4.d0*PI)) * 
     1                      (Fac(J-M)/Fac(J+M)))
              EX = EXP(COMPLEX(0,1) * M * P)
              LEG = Legendre(J,COS(T))

              WRITE (*,*) PREFAC, EX, LEG

              CALL EXIT(0)

            ENDIF
            C = C + 1

          ENDDO
        ENDDO

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

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

        DOUBLE PRECISION function Legendre(K,X) RESULT(VAL)

C
C       Calculates Legendre polynomials using recursion formula
C

        DOUBLE PRECISION :: X
        INTEGER :: K, N, M
        PARAMETER (M=2)
        DOUBLE PRECISION :: LEG(0:M)

        K = M

        if (K .EQ. 0) THEN
          LEG(0) = 1.d0
        
        ELSE IF (K .EQ. 1) THEN
          LEG(0) = 0.d0
          LEG(1) = 1.d0

        ELSE
          DO 

        ENDIF


        !IF (K .GE. 2) THEN
        !  DO N=1,K
        !    LEG(N+1) = 2.d0 * N * LEG(N) - N * LEG(N-1)
        !    LEG(N+1) = LEG(N+1) * (N+1.d0)
        !  ENDDO
        !ENDIF

        WRITE (*,*) LEG

        CALL EXIT(0)

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
        MAT(1,3) =  COSA * SING

        MAT(2,1) =  SINA * COSB * COSG + COSA * SING
        MAT(2,2) = -SINA * COSB * SING + COSA * COSG
        MAT(2,3) =  SINA * SINB

        MAT(3,1) = -SINB * COSG
        MAT(3,2) =  SINB * SING
        MAT(3,3) =  COSB

        END function

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C END FILE InteractionEnergyF.f90
