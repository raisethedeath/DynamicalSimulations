C FILE: Dynamical.F

        subroutine WignerD(J, K, M, beta)

C 
C       Uses the rotational quantum numbers to 
C       construct the Wigner-D rotational 
C       matrix elements following the work of Tajima
C       https://doi.org/10.1103/PhysRevC.91.014320
C
        DOUBLE PRECISION fact
        INTEGER J, K, M
        INTEGER nmin, nmax
        DOUBLE PRECISION lilD, lilW, bigW
        DOUBLE PRECISION beta

        nmin = MAX(0, K-M)
        nmax = MIN(J-M, J+K)

        lilD = 0

        DO N=nmin,nmax
            lilW = SQRT(fact(J+M) * fact(J-M) * fact(J+K) * fact(J-K))
     1                 / (fact(J-M-N) * fact(J+K-N) * fact(N+M-K)
     2                 * fact(N))

            bigW = lilW * ((COS(beta/2.d0))**(2.d0*J+K-M-2*N))
     1                  * ((-SIN(beta/2.d0))**(M-K+2*N))

            lilD = lilD + ((-1)**N) * BigW

            write (*,*) N, lilW, bigW, lilD
        ENDDO

        END





        DOUBLE PRECISION function fact(n)

C      
C       Calculates the factorial of n
C

        INTEGER n, p
        
        p = 1

        DO I=1,n
            p = p*I
        ENDDO

        fact = p

        END

C END FILE Dynamical.F
