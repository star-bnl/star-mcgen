        INTEGER FUNCTION ISAHEP(IDISA)
C Returns the HEP ID of stable particles starting with ISEJET ID
      PARAMETER (NID=103)
      INTEGER IDISA, IOFF, ISIGN, ITRAN(NID)
      DATA ITRAN/22,12,11,14,13,16,15,111,113,221,223,0,333,441,443
     1,551,553,211,213,311,313,-431,-433,541,543,3122,4232,321,323
     2,-411,-413,531,533,0,4132,4122,5232,-421,-423,511,513,3*0,5132
     3,5122,0,521,523,2224,2112,2114,3322,4424,0,1114,3312,3314,2212
     4,2214,3212,0,4322,4324,3112,0,4312,4314,3222,0,4212,4214,5322
     5,5324,4112,4114,5312,5314,4222,4224,5212,5214,3334,0,5112,5114
     6,2*0,5222,5224,0,4332,4334,8*0,5332,5334/
        IOFF = 0
        ISIGN = 1
        IAISA = IABS(IDISA)
        IF(IAISA.NE.IDISA) ISIGN = -1
        IF(IAISA.LT.19) THEN
          IOFF = -9
          GO TO 1000
        ENDIF
        IF(IAISA.EQ.20) THEN
          ISAHEP=130
          IF(ISIGN.LT.0) ISAHEP=310
          RETURN
        ENDIF
        IF(IAISA.LT.199) THEN
          IOFF = -102
          GO TO 1000
        ENDIF
        IF(IAISA.LT.299) THEN
          IOFF = -210
          GO TO 1000
        ENDIF
        IF(IAISA.LT.399) THEN
          IOFF = -318
          GO TO 1000
        ENDIF
        IF(IAISA.LT.499) THEN
          IOFF = -426
          GO TO 1000
        ENDIF
        IF(IAISA.LT.599) THEN
          IOFF = -534
          GO TO 1000
        ENDIF
        IF(IAISA.LT.1199) THEN
          IOFF = -1061
          GO TO 1000
        ENDIF
        IF(IAISA.LT.1299) THEN
          IOFF = -1169
          GO TO 1000
        ENDIF
        IF(IAISA.LT.1999) THEN
          IOFF = -1277
          GO TO 1000
        ENDIF
        IF(IAISA.LT.2199) THEN
          IOFF = -2104
          GO TO 1000
        ENDIF
        IF(IAISA.LT.2299) THEN
          IOFF = -2165
          GO TO 1000
        ENDIF
        IF(IAISA.LT.2999) THEN
          IOFF = -2273
          GO TO 1000
        ENDIF
        IF(IAISA.LT.3199) THEN
          IOFF = -3113
          GO TO 1000
        ENDIF
        IF(IAISA.LT.3299) THEN
          IOFF = -3205
          GO TO 1000
        ENDIF
        IF(IAISA.LT.3999) THEN
          IOFF = -3248
          GO TO 1000
        ENDIF
1000    ISAHEP = ITRAN(IAISA+IOFF)*ISIGN
        RETURN
        END
