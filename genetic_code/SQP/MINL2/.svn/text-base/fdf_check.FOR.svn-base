c TO DEFINE IN THE MAIN FUNCTION:(variables used)
c
C X(N), N = 3 (three parameters to be tuned)
C Initial guess of X = 1, 1, 1
C
C M = 20 (20 points where the functions are evaluated). Those points are:
C 		0.0 0.1 0.2 0.3 0.4 .......... z(i)=z(i-1)+.1 ...... 1.9
C
C IW = 356  ( IW = m*(n+3)+93+n*(n*3+33)/2+n*m+m )   is it right?
C
c the function to minimize is the sum on the Z1 points of:
C
C	2+8*Z1^2 - [ x1 + x2^2 * x3 * z1^3]  
C      
      SUBROUTINE FDF(N,M,X,DF,F)
      DOUBLE PRECISION X(N),DF(M,N),F(M)
c     
      DOUBLE PRECISION Z1
      INTEGER I
c
      Z1=0D0  	
c
      DO 10 I=1, M
c       FUNCTION       
	F(I)=2D0+8DO*(Z1**2) - (X(1)+(X(2)**2)*(Z1**2)*X(3)*Z1)
c       DERIVATIVE
	DF(I,1)=-1D0
        DF(I,2)=-2D0*X(2)*X(3)*(Z1**3)
 	DF(I,3)=-(X(2)**2)*(Z1**3)
c 	UPDATE Z1   
        Z1=Z1+1D-1
10    CONTINUE
      RETURN
      END
C
