C
C      TEST OF MI0L2   17.5.1991
C
c      CALL OPTI(0)
c      CALL OPTI(1)
c      STOP
c      END 
C
c      SUBROUTINE opti (method)

      SUBROUTINE opti (method, x,n,m, iw)

      INTEGER METHOD
	   INTEGER IW
		DOUBLE PRECISION X(*),W(IW),DX,EPS
C      DOUBLE PRECISION X(*),W(13560),DX,EPS

c      DOUBLE PRECISION X(3),W(3560),DX,EPS
      LOGICAL GRAD
      EXTERNAL FDF
      INTEGER N,M,MAXFUN,I,J,ICONTR
C	    TO PRINT THE SET OF X INHERITED FROM Population::tuning_individual(...)
C      WRITE(*,05) (X(I),I=1,N)
C   05  FORMAT(23X,' TI0L2 : PARAMETERS: ',D18.10/D52.10//)	
C     SET PARAMETERS
c      N=3
c      M=3
      EPS=1D-7
      DX=0.10D0
      MAXFUN=100
C	   IW=13560
C     INITIAL GUESS
c      X(1)=1D0
c      X(2)=1D0
c      X(3)=1D0
      ICONTR=METHOD
      GRAD=(ICONTR.GT.0)
      IF (.NOT.GRAD) DX=0.001
C      WRITE(*,*)
C     IF (.NOT.GRAD) WRITE(*,*) 'OPTIMIZE : NON-GRADIENT METHOD'
C      IF (GRAD)
C     1WRITE(*,*) 'OPTIMIZE : GRADIENT METHOD WITH APPROX. GRADIENT'
C     WRITE(*,*)
C

      CALL MI0L2(FDF,N,M,X,DX,EPS,MAXFUN,W,IW,ICONTR)
C
      IF (ICONTR.GE.0) GOTO 300
C
C     PARAMETER OUTSIDE RANGE
        WRITE(*,10) -ICONTR
   10   FORMAT(' INPUT ERROR. PARAMETER NUMBER ',I1,
     -         ' IS OUTSIDE ITS RANGE.')
        GOTO 400
C
C     RESULTS FROM OPTIMIZATION
  300 CONTINUE
C        IF (ICONTR.EQ.1) WRITE(*,40)
C   40   FORMAT(' SUM OF SQUARES FAILS TO DECREASE'/)
C        IF (ICONTR.EQ.2) WRITE(*,30)
C   30   FORMAT(' UPPER LIMIT FOR FUNCTION EVALUATIONS EXCEEDED.'/)
C        WRITE(*,31) (X(I),I=1,N)
C   31   FORMAT(23X,' SOLUTION: ',D18.10/D52.10//)
C        IF (GRAD) WRITE(*,35) MAXFUN
C   35   FORMAT(' NUMBER OF ITERATIONS: ',I4//)
C        IF (.NOT.GRAD) WRITE(*,33) MAXFUN
C  33   FORMAT(' NUMBER OF CALLS OF FDF: ',I4//)
C        WRITE(*,34) (W(J),J=1,M)
C   34   FORMAT(' FUNCTION VALUES AT THE SOLUTION: ',D18.10,2(/D52.10))
C       IF (.NOT.GRAD) THEN
C          WRITE(*,*) 'APPROX. OF PARTIAL DERIVATIVES'
C         DO 50 I=1,M
C            WRITE(*,'(D18.10)') (W(M+N*(I-1)+J),J=1,N)
C   50     CONTINUE
C       END IF
C
  400 CONTINUE
		method=ICONTR
      RETURN
C
      END
      SUBROUTINE FDF(N,M,X,F)
      DOUBLE PRECISION X(N),F(M)
		integer i
	


c      F(1)=1.5D0   - X(1)*(1D0-X(2))
c      F(2)=2.25D0  - X(1)*(1D0-X(2)**2)
c      F(3)=2.625D0 - X(1)*(1D0-X(2)**3)

c       Z1=0D0
c       DO I=1,M    
c        F(I)=  2.0+8.0*(Z1**2)-(X(1)+X(2)**2*x(3)*z1**3)
c        Z1=Z1+0.1
c       ENDDO	

	i=fdf_c(f,x,m, n) !call for c-function


	
       RETURN
       END
C
