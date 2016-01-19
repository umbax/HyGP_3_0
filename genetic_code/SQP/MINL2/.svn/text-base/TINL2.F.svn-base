C
C      TEST OF MINL2     21.11.1991
C
      CALL OPTI(0)
      CALL OPTI(1)
      END
C
      SUBROUTINE OPTI(METHOD)
      INTEGER METHOD
C
      PARAMETER (N=2, M=3, IW=M*(2*N+4)+(N/2+1)*(3*N+33)+93)
      DOUBLE PRECISION X(N),W(IW),DX,EPS
      INTEGER INDEX(4,2),ICONTR
      LOGICAL OPTIM
      EXTERNAL FDF
C     SET PAPAMETERS
      EPS=1D-10
      MAXFUN=25
C     SET INITIAL GUESS
      X(1)=1.D0
      X(2)=1.D0
      ICONTR=METHOD
      OPTIM=(ICONTR.GT.0)
      IF (.NOT.OPTIM) DX=0.001D0
C
      CALL MINL2(FDF,N,M,X,DX,EPS,MAXFUN,W,IW,ICONTR)
C
      IF (ICONTR.LT.0) GOTO 100
      IF (.NOT.OPTIM ) GOTO 200
                       GOTO 300
C
C     PARAMETER OUTSIDE RANGE
  100 CONTINUE
        WRITE(6,10) -ICONTR
   10   FORMAT(' INPUT ERROR. PARAMETER NUMBER ',I1,
     -         ' IS OUTSIDE ITS RANGE.')
        GOTO 400
C
C     RESULTS FROM GRADIENT TEST
  200 CONTINUE
        DO 20 K=2,4
          INDEX(K,1)=W(2*K+1)
          INDEX(K,2)=W(2*K+2)
   20   CONTINUE
        WRITE(6,21) (W(K),INDEX(K,1),INDEX(K,2),K=2,4)
   21   FORMAT(' TEST OF GRADIENTS '//
     -         ' MAXIMUM FORWARD  DIFFERENCE: ',D8.2,
     -         ' AT FUNCTION NO ',I1,' AND VARIABLE NO ',I1/
     -         ' MAXIMUM BACKWARD DIFFERENCE: ',D8.2,
     -         ' AT FUNCTION NO ',I1,' AND VARIABLE NO ',I1/
     -         ' MAXIMUM CENTRAL  DIFFERENCE: ',D8.2,
     -         ' AT FUNCTION NO ',I1,' AND VARIABLE NO ',I1)
        WRITE(6,22) W(1)
   22   FORMAT(//' MAXIMUM ELEMENT IN DF: ',D8.2)
        GOTO 400
C
C     RESULTS FROM OPTIMIZATION
  300 CONTINUE
        IF (ICONTR.EQ.2) WRITE(6,30)
   30   FORMAT(' UPPER LIMIT FOR FUNCTION EVALUATIONS EXCEEDED.'/)
        WRITE(6,31) (X(I),I=1,N),MAXFUN,(W(J),J=1,M)
   31   FORMAT(23X,' SOLUTION: ',D18.10/D52.10//
     -         ' NUMBER OF CALLS OF FDF: ',I4//
     -         ' FUNCTION VALUES AT THE SOLUTION: ',D18.10,2(/D52.10))
C
  400 CONTINUE
      RETURN
C
      END
      SUBROUTINE FDF(N,M,X,DF,F)
      DOUBLE PRECISION X(N),DF(M,N),F(M)
      F(1)=1.5D0   - X(1)*(1D0-X(2))
      F(2)=2.25D0  - X(1)*(1D0-X(2)**2)
      F(3)=2.625D0 - X(1)*(1D0-X(2)**3)
      DF(1,1)=X(2)-1D0
      DF(1,2)=X(1)
      DF(2,1)=X(2)**2-1D0
      DF(2,2)=2D0*X(1)*X(2)
      DF(3,1)=X(2)**3-1D0
      DF(3,2)=3D0*X(1)*X(2)**2
      RETURN
      END
C
