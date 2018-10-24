       SUBROUTINE VUMAT(
! READ ONLY - DO NOT MODIFY
     . NBLOCK, NDIR, NSHR, NSTATEV, NFIELDV, NPROPS, LANNEAL, STEPTIME,
     . TOTALTIME, DT, CMNAME, COORDMP, CHARLENGTH, PROPS, DENSITY,
     . STRAININC, RELSPININC, TEMPOLD, STRETCHOLD, DEFGRADOLD, FIELDOLD,
     . STRESSOLD, STATEOLD, ENERINTERNOLD, ENERINELASOLD, TEMPNEW,
     . STRETCHNEW, DEFGRADNEW, FIELDNEW,
! WRITE ONLY - DO NOT READ
     . STRESSNEW, STATENEW, ENERINTERNNEW, ENERINELASNEW)
!-----------------------------------------------------------------------
!     ABAQUS implicit variable declaration included in VABA_PARAM.INC
!     states the following:
!     a to h are real variables
!     o to z are real variables
!     i to n are integer variables
!-----------------------------------------------------------------------
      INCLUDE 'VABA_PARAM.INC'
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!     ABAQUS variables 
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      DIMENSION PROPS(NPROPS), DENSITY(NBLOCK), COORDMP(NBLOCK,*),
     . CHARLENGTH(*), STRAININC(NBLOCK,NDIR+NSHR), RELSPININC(*),
     . TEMPOLD(*), STRETCHOLD(*), DEFGRADOLD(*), FIELDOLD(*),
     . STRESSOLD(NBLOCK,NDIR+NSHR), STATEOLD(NBLOCK,NSTATEV),
     . ENERINTERNOLD(NBLOCK),  ENERINELASOLD(NBLOCK), TEMPNEW(*),
     . STRETCHNEW(*), DEFGRADNEW(*), FIELDNEW(*),
     . STRESSNEW(NBLOCK,NDIR+NSHR), STATENEW(NBLOCK,NSTATEV),
     . ENERINTERNNEW(NBLOCK), ENERINELASNEW(NBLOCK)
C
      CHARACTER*80 CMNAME
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!     Internal UMAT variables 
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      DIMENSION STRESS(NDIR+NSHR) ! Stress tensor inside UMAT
      DIMENSION DFDS(NDIR+NSHR)   ! Derivative of the yield function
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!     Internal UMAT variables could be declared as follow
!     but it is not required due to the VABA_PARAM.INC file
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!   Material parameters
!      REAL YOUNG       ! Young's modulus 
!      REAL POISS       ! Poisson's ratio
!      REAL C11,C12,C44 ! Elasticity matrix C components
!      REAL SIGMA0      ! Initial yield stress
!      REAL ET          ! Tangent modulus
!   Internal variables
!      REAL P           ! Equivalent plastic strain
!   Plasticity variables
!      REAL DLAMBDA     ! Plastic multiplier
!      REAL DDLAMBDA    ! Increment in plastic multiplier
!   Yield function variables
!      REAL F           ! Yield function
!      REAL PHI         ! Equivalent stress
!      REAL SIGMAY      ! Yield stress
!   Computational variables
!      REAL RESNOR      ! Convergence criterion
!      REAL TOL         ! Tolerance for the RMAP algorithm
!      INTEGER MXITER   ! Maximum number of iteration for the RMAP
!      INTEGER ITER     ! Number of iteration for the RMAP
!-----------------------------------------------------------------------
!     Read parameters from ABAQUS material card
!-----------------------------------------------------------------------
      YOUNG  = props(1)
      POISS  = props(2)
      SIGMA0 = props(3)
      ET     = props(4)
      TOL    = props(5)
      MXITER = props(6)
!-----------------------------------------------------------------------
!     Compute elasticity matrix
!-----------------------------------------------------------------------
      C11    = YOUNG*(1.0-POISS)/((1.0+POISS)*(1.0-2.0*POISS))
      C12    = POISS*C11/(1.0-POISS)
      C44    = 0.5*YOUNG/(1.0+POISS)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!     Loop over integration points
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      DO i=1,NBLOCK
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!        If time = 0 then pure elastic computation
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
         IF(TOTALTIME.eq.0.0)THEN
            STRESSNEW(i,1) = STRESSOLD(i,1)+C11*STRAININC(i,1)
     .                       +C12*STRAININC(i,2)+C12*STRAININC(i,3)
            STRESSNEW(i,2) = STRESSOLD(i,2)+C12*STRAININC(i,1)
     .                       +C11*STRAININC(i,2)+C12*STRAININC(i,3)
            STRESSNEW(i,3) = STRESSOLD(i,3)+C12*STRAININC(i,1)
     .                       +C12*STRAININC(i,2)+C11*STRAININC(i,3)
            STRESSNEW(i,4) = STRESSOLD(i,4)+C44*STRAININC(i,4)*2.0
            STRESSNEW(i,5) = STRESSOLD(i,5)+C44*STRAININC(i,5)*2.0
            STRESSNEW(i,6) = STRESSOLD(i,6)+C44*STRAININC(i,6)*2.0
         ELSE
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!        Elastic-predictor-corrector scheme
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!-----------------------------------------------------------------------
!           Elastic prediction
!-----------------------------------------------------------------------
            STRESS(1) = STRESSOLD(i,1)+C11*STRAININC(i,1)
     .                                +C12*STRAININC(i,2)
     .                                +C12*STRAININC(i,3)
            STRESS(2) = STRESSOLD(i,2)+C12*STRAININC(i,1)
     .                                +C11*STRAININC(i,2)
     .                                +C12*STRAININC(i,3)
            STRESS(3) = STRESSOLD(i,3)+C12*STRAININC(i,1)
     .                                +C12*STRAININC(i,2)
     .                                +C11*STRAININC(i,3)
            STRESS(4) = STRESSOLD(i,4)+C44*STRAININC(i,4)*2.0
            STRESS(5) = STRESSOLD(i,5)+C44*STRAININC(i,5)*2.0
            STRESS(6) = STRESSOLD(i,6)+C44*STRAININC(i,6)*2.0   
!-----------------------------------------------------------------------
!           Equivalent stress
!-----------------------------------------------------------------------
            PHI       = sqrt(STRESS(1)*STRESS(1)
     .                      +STRESS(2)*STRESS(2)
     .                      +STRESS(3)*STRESS(3)
     .                      -STRESS(1)*STRESS(2)
     .                      -STRESS(2)*STRESS(3)
     .                      -STRESS(3)*STRESS(1)
     .                 +3.0*(STRESS(4)*STRESS(4)
     .                      +STRESS(5)*STRESS(5)
     .                      +STRESS(6)*STRESS(6)))
!-----------------------------------------------------------------------
!           Equivalent plastic strain from previous time step
!-----------------------------------------------------------------------
            P        = STATEOLD(i,1)
!-----------------------------------------------------------------------
!           Update yield stress
!-----------------------------------------------------------------------
            SIGMAY   = SIGMA0+ET*P
!-----------------------------------------------------------------------
!           Compute yield function
!-----------------------------------------------------------------------
            F        = PHI-SIGMAY
!-----------------------------------------------------------------------
!           Initialize the plastic multiplier
!-----------------------------------------------------------------------
            DLAMBDA  = 0.0
!-----------------------------------------------------------------------
!           Compute the derivative of the yield function
!-----------------------------------------------------------------------
            IF(PHI.eq.0)THEN
               DENOM = 1.0
            ELSE
               DENOM = PHI
            ENDIF        
c
            DFDS(1) = (STRESS(1)-0.5*(STRESS(2)+STRESS(3)))/DENOM
            DFDS(2) = (STRESS(2)-0.5*(STRESS(3)+STRESS(1)))/DENOM
            DFDS(3) = (STRESS(3)-0.5*(STRESS(1)+STRESS(2)))/DENOM
            DFDS(4) = 3.0*STRESS(4)/DENOM
            DFDS(5) = 3.0*STRESS(5)/DENOM
            DFDS(6) = 3.0*STRESS(6)/DENOM
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!           Check for plasticity
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!-----------------------------------------------------------------------
!           Plastic flow
!-----------------------------------------------------------------------
            IF(F.GT.0.0)THEN ! Plastic flow
               DO ITER=1,MXITER
!-----------------------------------------------------------------------
!                 Compute increment in plastic multiplier
!-----------------------------------------------------------------------
                  DDLAMBDA = -F/(-3.0*C44-ET)
!-----------------------------------------------------------------------
!                 Update plastic multiplier
!-----------------------------------------------------------------------
                  DLAMBDA  = DLAMBDA+DDLAMBDA
!-----------------------------------------------------------------------
!                 Update equivalent plastic strain
!-----------------------------------------------------------------------
                  P      = P+DDLAMBDA
!-----------------------------------------------------------------------
!                 Update Yield stress
!-----------------------------------------------------------------------
                  SIGMAY = SIGMA0+ET*P
!-----------------------------------------------------------------------
!                 Update Yield function
!-----------------------------------------------------------------------
                  F      = PHI-3.0*C44*DLAMBDA-SIGMAY
!-----------------------------------------------------------------------
!                 Compute convergence criterion
!-----------------------------------------------------------------------
                  RESNOR = ABS(F/SIGMAY)
!-----------------------------------------------------------------------
!                 Check for convergence
!-----------------------------------------------------------------------
                  IF(RESNOR.LE.TOL)THEN ! RMAP has converged
!-----------------------------------------------------------------------
!                    Update the stress tensor
!-----------------------------------------------------------------------
                     STRESSNEW(i,1)= STRESS(1)-DLAMBDA*(C11*DFDS(1)
     .                                                 +C12*DFDS(2)
     .                                                 +C12*DFDS(3))
                     STRESSNEW(i,2)= STRESS(2)-DLAMBDA*(C12*DFDS(1)
     .                                                 +C11*DFDS(2)
     .                                                 +C12*DFDS(3))
                     STRESSNEW(i,3)= STRESS(3)-DLAMBDA*(C12*DFDS(1)
     .                                                 +C12*DFDS(2)
     .                                                 +C11*DFDS(3))
                     STRESSNEW(i,4)= STRESS(4)-DLAMBDA*C44*DFDS(4)
                     STRESSNEW(i,5)= STRESS(5)-DLAMBDA*C44*DFDS(5)
                     STRESSNEW(i,6)= STRESS(6)-DLAMBDA*C44*DFDS(6) 
!-----------------------------------------------------------------------
!                    Update the history variables
!-----------------------------------------------------------------------
                     STATENEW(i,1) = P
                     STATENEW(i,2) = PHI-3.0*C44*DLAMBDA
                     STATENEW(i,3) = F
                     STATENEW(i,4) = SIGMAY
                     STATENEW(i,5) = DGAMA
                     STATENEW(i,6) = ITER 
                     GOTO 90
                  ELSE ! RMAP has not converged yet
                     IF(ITER.eq.MXITER)THEN
                        write(*,*) 'RMAP has not converged'
                        write(*,*) 'Integration point',i
                        write(*,*) 'Convergence',RESNOR
                        write(*,*) 'dlambda',DLAMBDA,P
                        STOP
                     ENDIF
                  ENDIF          
               ENDDO
!-----------------------------------------------------------------------
!           Elastic point
!-----------------------------------------------------------------------
            ELSE
!-----------------------------------------------------------------------
!              Update the stress tensor
!-----------------------------------------------------------------------
               STRESSNEW(i,1) = STRESS(1)
               STRESSNEW(i,2) = STRESS(2)
               STRESSNEW(i,3) = STRESS(3)
               STRESSNEW(i,4) = STRESS(4)
               STRESSNEW(i,5) = STRESS(5)
               STRESSNEW(i,6) = STRESS(6)
!-----------------------------------------------------------------------
!              Update the history variables
!-----------------------------------------------------------------------
               STATENEW(i,1) = P
               STATENEW(i,2) = PHI
               STATENEW(i,3) = F
               STATENEW(i,4) = SIGMAY
               STATENEW(i,5) = 0.0
               STATENEW(i,6) = 0
            ENDIF
!-----------------------------------------------------------------------
!        End of loop over integration points
!-----------------------------------------------------------------------
         ENDIF           
  90     CONTINUE  
      ENDDO
!-----------------------------------------------------------------------
!     END SUBROUTINE
!-----------------------------------------------------------------------
      RETURN
      END