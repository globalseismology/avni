      SUBROUTINE PDAZ(EPLA,EPLO,AZIM,DELTA,XLAT,XLON)
C
C     INPUT PARAMETERS:
C        EPLA : EPICENTRAL LATITUDE
C        EPLO : EPICENTRAL LONGITUDE
C        AZIM : AZIMUTH OF RECEIVER AS SEEN FROM EPICENTER
C        DELTA: ANGULAR DISTANCE
C     OUTPUT PARAMETERS:
C        XLAT : RECEIVER LATITUDE
C        XLON : RECEIVER LONGITUDE
C
C     AZIMUTH IS MEASURED COUNTERCLOCKWISE FROM THE NORTH
C     ALL ANGLES ARE IN DEGREES
C
C       [HOFFMAN TRADITIONAL ROUTINE]
C
      DATA RADIAN/57.29578/
      PH21=(180.-AZIM)/RADIAN
      TH1=(90.-EPLA)/RADIAN
      PH1=EPLO/RADIAN
      DEL=DELTA/RADIAN
      STH1=SIN(TH1)
      CTH1=COS(TH1)
      SPH1=SIN(PH1)
      CPH1=COS(PH1)
      SDEL=SIN(DEL)
      CDEL=COS(DEL)
      CPH21=COS(PH21)
      SPH21=SIN(PH21)
      CTH2=-SDEL*CPH21*STH1+CDEL*CTH1
      STH2=SQRT(1.-CTH2*CTH2)
      CPH2=(SDEL*(CPH21*CTH1*CPH1-SPH21*SPH1)+CDEL*STH1*CPH1)
      SPH2=(SDEL*(CPH21*CTH1*SPH1+SPH21*CPH1)+CDEL*STH1*SPH1)
      TH2=ATAN2(STH2,CTH2)
      PH2=ATAN2(SPH2,CPH2)
      XLAT=90.-TH2*RADIAN
      XLON=PH2*RADIAN
      RETURN
      END
