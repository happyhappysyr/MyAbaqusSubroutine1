      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     1 PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
     2 KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,NPREDF,
     3 LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,PERIOD)
C
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*),
     1 SVARS(*),ENERGY(8),COORDS(MCRD,NNODE),U(NDOFEL),
     2 DU(MLVARX,*),V(NDOFEL),A(NDOFEL),TIME(2),PARAMS(*),
     3 JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*),
     4 PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)
      REAL BMATRIX(3,6)
      REAL BTRANS(6,3)
      REAL DMATRIX(3,3)


C      Checking for correct element types
      IF(NNODE.EQ.3) THEN
          WRITE(*,*) "ELEMENT HAS CORRECT NUMBER OF NODES (",NNODE,")"
        ELSE
          WRITE(*,*) "ELEMENT HAS WRONG NUMBER OF NODES (",NNODE,")"
        END IF
C      Assiging the Material Coordinates from COORD
      E = PROPS(1)
      PR = PROPS(2)
      THICKNESS = PROPS(3)
C      AREA = PROPS(4)
      WRITE(*,*) COORDS(1,1)
      WRITE(*,*) COORDS(2,1)
      WRITE(*,*) COORDS(3,1)
      WRITE(*,*) COORDS(1,2)
      WRITE(*,*) COORDS(2,2)
      WRITE(*,*) COORDS(3,2)
C      Assiging the Material Coordinates from COORD
      X_I = COORDS(1,1)
      X_J = COORDS(2,1)
      X_M = COORDS(3,1)
      Y_I = COORDS(1,2)
      Y_J = COORDS(2,2)
      Y_M = COORDS(3,2)

      BETA_I = Y_J - Y_M
      BETA_J = Y_M - Y_I
      BETA_M = Y_I - Y_J
C
      GAMMA_I = X_M - X_J
      GAMMA_J = X_I - X_M
      GAMMA_M = X_J - X_I
C
      AREA = ( X_I*(Y_J - Y_M) + X_J*(Y_M - Y_I) + X_M*(Y_I - Y_J))/2

C     Assembling B Matrix
      BMATRIX(1,1) = BETA_I/(2*AREA)
      BMATRIX(1,2) = 0
      BMATRIX(1,3) = BETA_J/(2*AREA)
      BMATRIX(1,4) = 0
      BMATRIX(1,5) = BETA_M/(2*AREA)
      BMATRIX(1,6) = 0
      BMATRIX(2,1) = 0
      BMATRIX(2,2) = GAMMA_I/(2*AREA)
      BMATRIX(2,3) = 0
      BMATRIX(2,4) = GAMMA_J/(2*AREA)
      BMATRIX(2,5) = 0
      BMATRIX(2,6) = GAMMA_M/(2*AREA)
      BMATRIX(3,1) = GAMMA_I/(2*AREA)
      BMATRIX(3,2) = BETA_I/(2*AREA)
      BMATRIX(3,3) = GAMMA_J/(2*AREA)
      BMATRIX(3,4) = BETA_J/(2*AREA)
      BMATRIX(3,5) = GAMMA_M/(2*AREA)
      BMATRIX(3,6) = BETA_M/(2*AREA)

C       Assembling [D] Matrix
      DMATRIX(1,1) = (E/(1-PR**2))*(1)
      DMATRIX(2,1) = (E/(1-PR**2))*(PR)
      DMATRIX(3,1) = 0
      DMATRIX(1,2) = (E/(1-PR**2))*(PR)
      DMATRIX(2,2) = (E/(1-PR**2))*(1)
      DMATRIX(3,2) = 0
      DMATRIX(1,3) = 0
      DMATRIX(2,3) = 0
      DMATRIX(3,3) = (E/(1-PR**2))*((1-2*PR)/2)

C      MATRIX MULTIPLICATION FOR STIFFNESS MATRIX
      BTRANS = TRANSPOSE(BMATRIX)
c      BTRANS = TRANSPOSE(BMATRIX)*(THICKNESS*AREA)
      AMATRX = MATMUL(MATMUL(BTRANS,DMATRIX),BMATRIX)

      RETURN
      END
