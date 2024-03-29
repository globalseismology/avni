      SUBROUTINE RAYSEQ (NP,NS,NN,ISS,IRR,ICNSx,ICNR,ICN,ISEQ,
     #idirseg,NSEG)
      INTEGER NP(NN),NS(NN),ICN(NN),ISEQ(NN),NPS(20,2),idirseg(nn)
      common /plevel/iprtlv
C SUBROUTINE TO FORM ARRAY DESCRIBING RAY SEQUENCE THROUGH LAYERS .
C OUTPUT IS ISEQ(NSEG) WHICH IS A SEQUENCE OF LAYER NUMBERS , POSITIVE
C FOR P SEGMENT , NEGATIVE FOR S SEGMENT . NOTE THE REFLECTION/
C TRANSMISSION COEFFICIENTS ARE NOT USED , SO THE RAY FORMED IS ONLY A
C KINEMATIC ANALOGUE . IT IS DIFFICULT TO SEQUENCE A GENERAL DYNAMIC
C ANALOGUE SO THIS IS NOT ATTEMPTED. THE RAY MINIMIZES THE CONVERSIONS .
C
C FIND TOTAL NUMBER OF SEGMENTS AND TRANSFER TO INTERNAL ARRAYS .
      icns=icnsx
      if(iss.gt.4) icns=icns+1
      NSEQ=0
      DO 1 I=1,NN
      NSEQ=NSEQ+NP(I)+NS(I)
      NPS(I,1)=NP(I)
1     NPS(I,2)=NS(I)
      NSEG=NSEQ
C FIND RECEIVER LAYER .
      DO 20 IL=2,NN
      IF (ICNR.LT.ICN(IL)) GO TO 21
20    CONTINUE
21    IPR=MOD(IRR,4)
      IL=IL-1
      IF (IRR.LE.4.OR.IRR.GE.9) GO TO 22
C DOWNGOING RAY AT RECEIVER .
      IF (ICNR.EQ.ICN(IL)) GO TO 25
C SUBTRACT RECEIVER SEGMENT .
24    NPS(IL,IPR)=NPS(IL,IPR)-1
      ISEQ(NSEQ)=ISIGN(IL,3-IPR-IPR)
      idirseg(nseq)=-1
      IF (NSEQ.EQ.1) GO TO 101
      NSEQ=NSEQ-1
      GO TO 25
C UPGOING RAY AT RECEIVER .
22    IF (ICNR.EQ.ICN(IL)) GO TO 24
25    N0=1
C FIND SOURCE LAYER .
      DO 2 IL=2,NN
      IF (ICNS.LT.ICN(IL)) GO TO 3
2     CONTINUE
3     IPS=MOD(ISS,4)
      IL=IL-1
      IF (ISS.LE.4) GO TO 4
      IDIR=-1
      IF (ICNS.EQ.ICN(IL)) GO TO 5
7     N0=2
      ISEQ(1)=ISIGN(IL,3-IPS-IPS)
      idirseg(1)=idir
      NPS(IL,IPS)=NPS(IL,IPS)-1
      GO TO 5
4     IDIR=1
      IF (ICNS.EQ.ICN(IL)) GO TO 7
C
5     IF (N0.GT.NSEQ) GO TO 101
      DO 8 N=N0,NSEQ
      K=0
      IL=IL+IDIR
      IF (IDIR.EQ.1) GO TO 10
C UPWARD TRAVELLING RAY .
15    IF (IL.EQ.0) GO TO 11
      K=K+1
      IF (K.GT.2) GO TO 100
      IF (NPS(IL,1)+NPS(IL,2).NE.1) GO TO 16
      IF (IL.LT.ILMAX) GO TO 11
      ILMAX=IL-1
16    IF (NPS(IL,IPS).NE.0) GO TO 13
      IPS=3-IPS
      IF (NPS(IL,IPS).NE.0) GO TO 13
      IPS=3-IPS
C CHANGE DIRECTION TO DOWNWARD TRAVELLING .
11    IDIR=1
      IPS=3-IPR
      IL=IL+1
C DOWNWARD TRAVELLING RAY .
10    K=K+1
      IF (K.GT.2) GO TO 100
      IF (NPS(IL,IPS).NE.0) GO TO 13
      IPS=3-IPS
      IF (NPS(IL,IPS).NE.0) GO TO 13
C CHANGE DIRECTION TO UPWARD TRAVELLING .
      IL=IL-1
      ILMAX=IL
      IDIR=-1
      IPS=3-IPS
      GO TO 15
13    NPS(IL,IPS)=NPS(IL,IPS)-1
      if(iprtlv.gt.2) then
        write(6,"('il,ips',3i5)") il,ips,isign(il,3-ips-ips)
      endif
      idirseg(n)=idir
8     ISEQ(N)=ISIGN(IL,3-IPS-IPS)
C DECODING COMPLETED - CHECK ALL ZEROS .
101   DO 103 I=1,NN
      IF (NPS(I,1).NE.0.OR.NPS(I,2).NE.0) GO TO 100
103   CONTINUE
      RETURN
100   WRITE (6,600)
      STOP
600   FORMAT (' ERROR IN RAY SEQUENCING')
      END
