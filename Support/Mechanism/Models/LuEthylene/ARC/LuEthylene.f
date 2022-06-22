C----------------------------------------------------------------------C
C                                                                  
C   A 22-Species Reduced Mechanism for C2H4-ai

C   Z. Luo, C.S. Yoo, E.S. Richardson, J.H. Chen, C.K. Law, and T.F. Lu,
C   "Chemical explosive mode analysis for a turbulent lifted ethylene
C    jet flame in highly-heated coflow,"
C   Combustion and Flame,
C   doi:10.1016/j.combustflame.2011.05.023, 2011
C                                                                 
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKWYP  (P, T, Y, ICKWRK, RCKWRK, WDOT)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION Y(*),WDOT(*),ICKWRK(*),RCKWRK(*)
      DIMENSION RF(206),RB(206),RKLOW(21),C(22)
      DIMENSION XQ(10)
C
      CALL YTCP(P, T, Y, C)
      CALL RATT(T, RF, RB, RKLOW)
      CALL RATX(T, C, RF, RB, RKLOW)
      CALL QSSA(RF, RB, XQ)
      CALL RDOT(RF, RB, WDOT)
C
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE YTCP (P, T, Y, C)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION Y(*), C(*)
      DATA SMALL/1D-50/
C
      C(1) = Y(1)/2.01593995D0
      C(2) = Y(2)/1.00796998D0
      C(3) = Y(3)/1.59994001D1
      C(4) = Y(4)/3.19988003D1
      C(5) = Y(5)/1.70073701D1
      C(6) = Y(6)/1.80153401D1
      C(7) = Y(7)/3.30067703D1
      C(8) = Y(8)/3.40147402D1
      C(9) = Y(9)/1.50350603D1
      C(10) = Y(10)/1.60430303D1
      C(11) = Y(11)/2.80105505D1
      C(12) = Y(12)/4.40099506D1
      C(13) = Y(13)/3.00264904D1
      C(14) = Y(14)/2.60382407D1
      C(15) = Y(15)/2.80541806D1
      C(16) = Y(16)/3.00701206D1
      C(17) = Y(17)/4.10296708D1
      C(18) = Y(18)/4.20376408D1
      C(19) = Y(19)/4.40535808D1
      C(20) = Y(20)/4.1073301D1
      C(21) = Y(21)/4.20812709D1
      C(22) = Y(22)/2.80133991D1
C
      SUM = 0D0
      DO K = 1, 22
         SUM = SUM + C(K)
      ENDDO
      SUM = P/(SUM*T*8.314510D7)
C
      DO K = 1, 22
         C(K) = C(K) * SUM
      ENDDO
C
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE RATT (T, RF, RB, RKLOW)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      PARAMETER (RU=8.31451D7, SMALL=1.D-200, PATM=1.01325D6)
      DIMENSION RF(*), RB(*), RKLOW(*)
      DIMENSION SMH(31), EG(31)
C
      ALOGT = LOG(T)
      TI = 1D0/T
      TI2 = TI*TI
C
      CALL RDSMH (T, SMH)
      EG(1) = EXP(SMH(1))
      EG(2) = EXP(SMH(2))
      EG(3) = EXP(SMH(3))
      EG(4) = EXP(SMH(4))
      EG(5) = EXP(SMH(5))
      EG(6) = EXP(SMH(6))
      EG(7) = EXP(SMH(7))
      EG(8) = EXP(SMH(8))
      EG(9) = EXP(SMH(9))
      EG(10) = EXP(SMH(10))
      EG(11) = EXP(SMH(11))
      EG(12) = EXP(SMH(12))
      EG(13) = EXP(SMH(13))
      EG(14) = EXP(SMH(14))
      EG(15) = EXP(SMH(15))
      EG(16) = EXP(SMH(16))
      EG(17) = EXP(SMH(17))
      EG(18) = EXP(SMH(18))
      EG(19) = EXP(SMH(19))
      EG(20) = EXP(SMH(20))
      EG(21) = EXP(SMH(21))
      EG(22) = EXP(SMH(22))
      EG(23) = EXP(SMH(23))
      EG(24) = EXP(SMH(24))
      EG(25) = EXP(SMH(25))
      EG(26) = EXP(SMH(26))
      EG(27) = EXP(SMH(27))
      EG(28) = EXP(SMH(28))
      EG(29) = EXP(SMH(29))
      EG(30) = EXP(SMH(30))
      EG(31) = EXP(SMH(31))
      PFAC1 = PATM / (RU*T)
      PFAC2 = PFAC1*PFAC1
      PFAC3 = PFAC2*PFAC1
C
C
      RF(1) = EXP(3.20498617D1 -7.25286183D3*TI)
      EQK = EG(3)*EG(5)/EG(2)/EG(4)
      RB(1) = RF(1) / MAX(EQK, SMALL)
      RF(2) = EXP(1.08197783D1 +2.67D0*ALOGT -3.16523284D3*TI)
      EQK = EG(2)*EG(5)/EG(1)/EG(3)
      RB(2) = RF(2) / MAX(EQK, SMALL)
      RF(3) = EXP(1.9190789D1 +1.51D0*ALOGT -1.72603317D3*TI)
      EQK = EG(2)*EG(6)/EG(1)/EG(5)
      RB(3) = RF(3) / MAX(EQK, SMALL)
      RF(4) = EXP(1.0482906D1 +2.4D0*ALOGT +1.06178717D3*TI)
      EQK = EG(3)*EG(6)/EG(5)/EG(5)
      RB(4) = RF(4) / MAX(EQK, SMALL)
      RF(5) = EXP(4.14465317D1 -1D0*ALOGT)
      EQK = EG(1)/EG(2)/EG(2)/PFAC1
      RB(5) = RF(5) / MAX(EQK, SMALL)
      RF(6) = EXP(3.90385861D1 -6D-1*ALOGT)
      EQK = EG(1)/EG(2)/EG(2)/PFAC1
      RB(6) = RF(6) / MAX(EQK, SMALL)
      RF(7) = EXP(4.55408762D1 -1.25D0*ALOGT)
      EQK = EG(1)/EG(2)/EG(2)/PFAC1
      RB(7) = RF(7) / MAX(EQK, SMALL)
      RF(8) = EXP(4.775645D1 -2D0*ALOGT)
      EQK = EG(1)/EG(2)/EG(2)/PFAC1
      RB(8) = RF(8) / MAX(EQK, SMALL)
      RF(9) = 4D1*RF(8)
      EQK = EG(6)/EG(2)/EG(5)/PFAC1
      RB(9) = RF(9) / MAX(EQK, SMALL)
      RF(10) = 5D-1*RF(5)
      EQK = EG(5)/EG(2)/EG(3)/PFAC1
      RB(10) = RF(10) / MAX(EQK, SMALL)
      RF(11) = 1.2D-1*RF(5)
      EQK = EG(4)/EG(3)/EG(3)/PFAC1
      RB(11) = RF(11) / MAX(EQK, SMALL)
      RF(12) = EXP(4.24761511D1 -8.6D-1*ALOGT)
      EQK = EG(7)/EG(2)/EG(4)/PFAC1
      RB(12) = RF(12) / MAX(EQK, SMALL)
      RF(13) = EXP(4.71503141D1 -1.72D0*ALOGT)
      EQK = EG(7)/EG(2)/EG(4)/PFAC1
      RB(13) = RF(13) / MAX(EQK, SMALL)
      RF(14) = EXP(4.42511034D1 -7.6D-1*ALOGT)
      EQK = EG(7)/EG(2)/EG(4)/PFAC1
      RB(14) = RF(14) / MAX(EQK, SMALL)
      RF(15) = EXP(4.47046282D1 -1.24D0*ALOGT)
      EQK = EG(7)/EG(2)/EG(4)/PFAC1
      RB(15) = RF(15) / MAX(EQK, SMALL)
      RF(16) = EXP(3.19350862D1 -3.7D-1*ALOGT)
      EQK = EG(8)/EG(5)/EG(5)/PFAC1
      RB(16) = RF(16) / MAX(EQK, SMALL)
      RF(17) = EXP(2.90097872D1 -3.37658384D2*TI)
      EQK = EG(3)*EG(6)/EG(2)/EG(7)
      RB(17) = RF(17) / MAX(EQK, SMALL)
      RF(18) = EXP(3.04404238D1 -4.12637667D2*TI)
      EQK = EG(1)*EG(4)/EG(2)/EG(7)
      RB(18) = RF(18) / MAX(EQK, SMALL)
      RF(19) = EXP(3.18908801D1 -1.50965D2*TI)
      EQK = EG(5)*EG(5)/EG(2)/EG(7)
      RB(19) = RF(19) / MAX(EQK, SMALL)
      RF(20) = 2D13
      EQK = EG(4)*EG(5)/EG(3)/EG(7)
      RB(20) = RF(20) / MAX(EQK, SMALL)
      RF(21) = EXP(3.14683206D1 +2.51608334D2*TI)
      EQK = EG(4)*EG(6)/EG(5)/EG(7)
      RB(21) = RF(21) / MAX(EQK, SMALL)
      RF(22) = EXP(2.55908003D1 +8.20243168D2*TI)
      EQK = EG(4)*EG(8)/EG(7)/EG(7)
      RB(22) = RF(22) / MAX(EQK, SMALL)
      RF(23) = EXP(3.36712758D1 -6.03860001D3*TI)
      EQK = EG(4)*EG(8)/EG(7)/EG(7)
      RB(23) = RF(23) / MAX(EQK, SMALL)
      RF(24) = EXP(1.6308716D1 +2D0*ALOGT -2.61672667D3*TI)
      EQK = EG(1)*EG(7)/EG(2)/EG(8)
      RB(24) = RF(24) / MAX(EQK, SMALL)
      RF(25) = EXP(2.99336062D1 -1.81158D3*TI)
      EQK = EG(5)*EG(6)/EG(2)/EG(8)
      RB(25) = RF(25) / MAX(EQK, SMALL)
      RF(26) = EXP(1.60803938D1 +2D0*ALOGT -2.01286667D3*TI)
      EQK = EG(5)*EG(7)/EG(3)/EG(8)
      RB(26) = RF(26) / MAX(EQK, SMALL)
      RF(27) = EXP(2.81906369D1 -1.61029334D2*TI)
      EQK = EG(6)*EG(7)/EG(5)/EG(8)
      RB(27) = RF(27) / MAX(EQK, SMALL)
      RF(28) = EXP(3.39940492D1 -4.81075134D3*TI)
      EQK = EG(6)*EG(7)/EG(5)/EG(8)
      RB(28) = RF(28) / MAX(EQK, SMALL)
      RF(29) = EXP(3.40312786D1 -1.50965D3*TI)
      EQK = EG(15)/EG(3)/EG(14)/PFAC1
      RB(29) = RF(29) / MAX(EQK, SMALL)
      RF(30) = EXP(1.76783433D1 +1.228D0*ALOGT -3.52251667D1*TI)
      EQK = EG(2)*EG(15)/EG(5)/EG(14)
      RB(30) = RF(30) / MAX(EQK, SMALL)
      RF(31) = EXP(1.75767107D1 +1.5D0*ALOGT -4.00560467D4*TI)
      EQK = EG(17)/EG(1)/EG(14)/PFAC1
      RB(31) = RF(31) / MAX(EQK, SMALL)
      RF(32) = EXP(2.85473118D1 -2.40537567D4*TI)
      EQK = EG(3)*EG(15)/EG(4)/EG(14)
      RB(32) = RF(32) / MAX(EQK, SMALL)
      RF(33) = EXP(3.26416564D1 -1.18759134D4*TI)
      EQK = EG(5)*EG(15)/EG(7)/EG(14)
      RB(33) = RF(33) / MAX(EQK, SMALL)
      RF(34) = 5.7D13
      EQK = EG(2)*EG(14)/EG(3)/EG(9)
      RB(34) = RF(34) / MAX(EQK, SMALL)
      RF(35) = 3D13
      EQK = EG(2)*EG(16)/EG(5)/EG(9)
      RB(35) = RF(35) / MAX(EQK, SMALL)
      RF(36) = EXP(1.85223344D1 +1.79D0*ALOGT -8.40371835D2*TI)
      EQK = EG(2)*EG(10)/EG(1)/EG(9)
      RB(36) = RF(36) / MAX(EQK, SMALL)
      RF(37) = EXP(2.93732401D1 +3.79928584D2*TI)
      EQK = EG(2)*EG(17)/EG(6)/EG(9)
      RB(37) = RF(37) / MAX(EQK, SMALL)
      RF(38) = 3.3D13
      EQK = EG(3)*EG(16)/EG(4)/EG(9)
      RB(38) = RF(38) / MAX(EQK, SMALL)
      RF(39) = 5D13
      EQK = EG(25)/EG(9)/EG(14)/PFAC1
      RB(39) = RF(39) / MAX(EQK, SMALL)
      RF(40) = EXP(2.88547965D1 -3.47219501D2*TI)
      EQK = EG(14)*EG(16)/EG(9)/EG(15)
      RB(40) = RF(40) / MAX(EQK, SMALL)
      RF(41) = EXP(2.77171988D1 +4.8D-1*ALOGT +1.30836334D2*TI)
      EQK = EG(17)/EG(2)/EG(16)/PFAC1
      RB(41) = RF(41) / MAX(EQK, SMALL)
      RF(42) = 7.34D13
      EQK = EG(1)*EG(14)/EG(2)/EG(16)
      RB(42) = RF(42) / MAX(EQK, SMALL)
      RF(43) = 3D13
      EQK = EG(5)*EG(14)/EG(3)/EG(16)
      RB(43) = RF(43) / MAX(EQK, SMALL)
      RF(44) = 3D13
      EQK = EG(2)*EG(15)/EG(3)/EG(16)
      RB(44) = RF(44) / MAX(EQK, SMALL)
      RF(45) = 5D13
      EQK = EG(6)*EG(14)/EG(5)/EG(16)
      RB(45) = RF(45) / MAX(EQK, SMALL)
      RF(46) = EXP(3.9769885D1 -1D0*ALOGT -8.55468335D3*TI)
      EQK = EG(2)*EG(14)/EG(16)*PFAC1
      RB(46) = RF(46) / MAX(EQK, SMALL)
      RF(47) = EXP(2.96591694D1 -2.01286667D2*TI)
      EQK = EG(7)*EG(14)/EG(4)/EG(16)
      RB(47) = RF(47) / MAX(EQK, SMALL)
      RF(48) = EXP(3.77576522D1 -8D-1*ALOGT)
      EQK = EG(12)/EG(2)/EG(10)/PFAC1
      RB(48) = RF(48) / MAX(EQK, SMALL)
      RF(49) = EXP(1.31223634D1 +2D0*ALOGT -3.63825651D3*TI)
      EQK = EG(2)*EG(12)/EG(1)/EG(10)
      RB(49) = RF(49) / MAX(EQK, SMALL)
      RF(50) = 8D13
      EQK = EG(2)*EG(16)/EG(3)/EG(10)
      RB(50) = RF(50) / MAX(EQK, SMALL)
      RF(51) = EXP(2.99880944D1 -7.54825001D2*TI)
      EQK = EG(5)*EG(16)/EG(4)/EG(10)
      RB(51) = RF(51) / MAX(EQK, SMALL)
      RF(52) = 2.5D-1*RF(51)
      EQK = EG(2)*EG(2)*EG(15)/EG(4)/EG(10)*PFAC1
      RB(52) = RF(52) / MAX(EQK, SMALL)
      RF(53) = 2D13
      EQK = EG(2)*EG(17)/EG(5)/EG(10)
      RB(53) = RF(53) / MAX(EQK, SMALL)
      RF(54) = EXP(1.62403133D1 +2D0*ALOGT -1.50965D3*TI)
      EQK = EG(6)*EG(9)/EG(5)/EG(10)
      RB(54) = RF(54) / MAX(EQK, SMALL)
      RF(55) = 2D13
      EQK = EG(5)*EG(17)/EG(7)/EG(10)
      RB(55) = RF(55) / MAX(EQK, SMALL)
      RF(56) = EXP(2.74203001D1 +5D-1*ALOGT -2.26950717D3*TI)
      EQK = EG(26)/EG(10)/EG(14)/PFAC1
      RB(56) = RF(56) / MAX(EQK, SMALL)
      RF(57) = 4D13
      EQK = EG(2)*EG(19)/EG(9)/EG(10)
      RB(57) = RF(57) / MAX(EQK, SMALL)
      RF(58) = 3.2D13
      EQK = EG(1)*EG(19)/EG(10)/EG(10)
      RB(58) = RF(58) / MAX(EQK, SMALL)
      RF(59) = EXP(3.03390713D1 -3.01930001D2*TI)
      EQK = EG(10)/EG(11)
      RB(59) = RF(59) / MAX(EQK, SMALL)
      RF(60) = 3D13
      EQK = EG(1)*EG(9)/EG(2)/EG(11)
      RB(60) = RF(60) / MAX(EQK, SMALL)
      RF(61) = 1.5D13
      EQK = EG(1)*EG(14)/EG(3)/EG(11)
      RB(61) = RF(61) / MAX(EQK, SMALL)
      RF(62) = 1.5D13
      EQK = EG(2)*EG(16)/EG(3)/EG(11)
      RB(62) = RF(62) / MAX(EQK, SMALL)
      RF(63) = 3D13
      EQK = EG(2)*EG(17)/EG(5)/EG(11)
      RB(63) = RF(63) / MAX(EQK, SMALL)
      RF(64) = 7D13
      EQK = EG(2)*EG(12)/EG(1)/EG(11)
      RB(64) = RF(64) / MAX(EQK, SMALL)
      RF(65) = 2.8D13
      EQK = EG(2)*EG(5)*EG(14)/EG(4)/EG(11)*PFAC1
      RB(65) = RF(65) / MAX(EQK, SMALL)
      RF(66) = 1.2D13
      EQK = EG(6)*EG(14)/EG(4)/EG(11)
      RB(66) = RF(66) / MAX(EQK, SMALL)
      RF(67) = 3D13
      EQK = EG(10)/EG(11)
      RB(67) = RF(67) / MAX(EQK, SMALL)
      RF(68) = 9D12
      EQK = EG(10)/EG(11)
      RB(68) = RF(68) / MAX(EQK, SMALL)
      RF(69) = 7D12
      EQK = EG(10)/EG(11)
      RB(69) = RF(69) / MAX(EQK, SMALL)
      RF(70) = 1.4D13
      EQK = EG(14)*EG(17)/EG(11)/EG(15)
      RB(70) = RF(70) / MAX(EQK, SMALL)
      RF(71) = EXP(2.7014835D1 +4.54D-1*ALOGT -1.30836334D3*TI)
      EQK = EG(18)/EG(2)/EG(17)/PFAC1
      RB(71) = RF(71) / MAX(EQK, SMALL)
      RF(72) = EXP(2.38587601D1 +1.05D0*ALOGT -1.64803459D3*TI)
      EQK = EG(1)*EG(16)/EG(2)/EG(17)
      RB(72) = RF(72) / MAX(EQK, SMALL)
      RF(73) = EXP(3.12945828D1 -1.781387D3*TI)
      EQK = EG(5)*EG(16)/EG(3)/EG(17)
      RB(73) = RF(73) / MAX(EQK, SMALL)
      RF(74) = EXP(2.19558261D1 +1.18D0*ALOGT +2.2493785D2*TI)
      EQK = EG(6)*EG(16)/EG(5)/EG(17)
      RB(74) = RF(74) / MAX(EQK, SMALL)
      RF(75) = EXP(3.22361913D1 -2.01286667D4*TI)
      EQK = EG(7)*EG(16)/EG(4)/EG(17)
      RB(75) = RF(75) / MAX(EQK, SMALL)
      RF(76) = EXP(2.76310211D1 -4.02573334D3*TI)
      EQK = EG(8)*EG(16)/EG(7)/EG(17)
      RB(76) = RF(76) / MAX(EQK, SMALL)
      RF(77) = EXP(3.21806786D1 +2.59156584D2*TI)
      EQK = EG(2)*EG(26)/EG(9)/EG(17)
      RB(77) = RF(77) / MAX(EQK, SMALL)
      RF(78) = EXP(3.70803784D1 -6.3D-1*ALOGT -1.92731984D2*TI)
      EQK = EG(13)/EG(2)/EG(12)/PFAC1
      RB(78) = RF(78) / MAX(EQK, SMALL)
      RF(79) = 8.43D13
      EQK = EG(2)*EG(17)/EG(3)/EG(12)
      RB(79) = RF(79) / MAX(EQK, SMALL)
      RF(80) = EXP(1.78408622D1 +1.6D0*ALOGT -2.72743434D3*TI)
      EQK = EG(6)*EG(10)/EG(5)/EG(12)
      RB(80) = RF(80) / MAX(EQK, SMALL)
      RF(81) = 2.501D13
      EQK = EG(6)*EG(11)/EG(5)/EG(12)
      RB(81) = RF(81) / MAX(EQK, SMALL)
      RF(82) = EXP(3.10595094D1 -1.449264D4*TI)
      EQK = EG(3)*EG(18)/EG(4)/EG(12)
      RB(82) = RF(82) / MAX(EQK, SMALL)
      RF(83) = EXP(2.43067848D1 -4.49875701D3*TI)
      EQK = EG(5)*EG(17)/EG(4)/EG(12)
      RB(83) = RF(83) / MAX(EQK, SMALL)
      RF(84) = 1D12
      EQK = EG(4)*EG(13)/EG(7)/EG(12)
      RB(84) = RF(84) / MAX(EQK, SMALL)
      RF(85) = 1.34D13
      EQK = EG(5)*EG(18)/EG(7)/EG(12)
      RB(85) = RF(85) / MAX(EQK, SMALL)
      RF(86) = EXP(1.01064284D1 +2.47D0*ALOGT -2.60666234D3*TI)
      EQK = EG(7)*EG(13)/EG(8)/EG(12)
      RB(86) = RF(86) / MAX(EQK, SMALL)
      RF(87) = 3D13
      EQK = EG(2)*EG(21)/EG(9)/EG(12)
      RB(87) = RF(87) / MAX(EQK, SMALL)
      RF(88) = 8.48D12
      EQK = EG(13)*EG(14)/EG(12)/EG(16)
      RB(88) = RF(88) / MAX(EQK, SMALL)
      RF(89) = 1.8D13
      EQK = EG(28)/EG(12)/EG(16)/PFAC1
      RB(89) = RF(89) / MAX(EQK, SMALL)
      RF(90) = EXP(8.10772006D0 +2.81D0*ALOGT -2.94884967D3*TI)
      EQK = EG(13)*EG(16)/EG(12)/EG(17)
      RB(90) = RF(90) / MAX(EQK, SMALL)
      RF(91) = 4D13
      EQK = EG(2)*EG(22)/EG(10)/EG(12)
      RB(91) = RF(91) / MAX(EQK, SMALL)
      RF(92) = EXP(3.01159278D1 +2.86833501D2*TI)
      EQK = EG(2)*EG(22)/EG(11)/EG(12)
      RB(92) = RF(92) / MAX(EQK, SMALL)
      RF(93) = EXP(3.75927776D1 -9.7D-1*ALOGT -3.11994334D2*TI)
      EQK = EG(24)/EG(12)/EG(12)/PFAC1
      RB(93) = RF(93) / MAX(EQK, SMALL)
      RF(94) = EXP(2.9238457D1 +1D-1*ALOGT -5.33409668D3*TI)
      EQK = EG(2)*EG(23)/EG(12)/EG(12)
      RB(94) = RF(94) / MAX(EQK, SMALL)
      RF(95) = 5D13
      EQK = EG(14)*EG(22)/EG(12)/EG(25)
      RB(95) = RF(95) / MAX(EQK, SMALL)
      RF(96) = 2D13
      EQK = EG(1)*EG(17)/EG(2)/EG(18)
      RB(96) = RF(96) / MAX(EQK, SMALL)
      RF(97) = 3.2D13
      EQK = EG(5)*EG(12)/EG(2)/EG(18)
      RB(97) = RF(97) / MAX(EQK, SMALL)
      RF(98) = 1.6D13
      EQK = EG(6)*EG(11)/EG(2)/EG(18)
      RB(98) = RF(98) / MAX(EQK, SMALL)
      RF(99) = 1D13
      EQK = EG(5)*EG(17)/EG(3)/EG(18)
      RB(99) = RF(99) / MAX(EQK, SMALL)
      RF(100) = 5D12
      EQK = EG(6)*EG(17)/EG(5)/EG(18)
      RB(100) = RF(100) / MAX(EQK, SMALL)
      RF(101) = EXP(-2.84796532D1 +7.6D0*ALOGT +1.77635484D3*TI)
      EQK = EG(7)*EG(17)/EG(4)/EG(18)
      RB(101) = RF(101) / MAX(EQK, SMALL)
      RF(102) = EXP(2.03077504D1 +1.62D0*ALOGT -5.45486868D3*TI)
      EQK = EG(1)*EG(12)/EG(2)/EG(13)
      RB(102) = RF(102) / MAX(EQK, SMALL)
      RF(103) = EXP(2.07430685D1 +1.5D0*ALOGT -4.32766334D3*TI)
      EQK = EG(5)*EG(12)/EG(3)/EG(13)
      RB(103) = RF(103) / MAX(EQK, SMALL)
      RF(104) = EXP(1.84206807D1 +1.6D0*ALOGT -1.570036D3*TI)
      EQK = EG(6)*EG(12)/EG(5)/EG(13)
      RB(104) = RF(104) / MAX(EQK, SMALL)
      RF(105) = 6D13
      EQK = EG(2)*EG(22)/EG(9)/EG(13)
      RB(105) = RF(105) / MAX(EQK, SMALL)
      RF(106) = EXP(1.47156719D1 +2D0*ALOGT -4.16160184D3*TI)
      EQK = EG(12)*EG(12)/EG(10)/EG(13)
      RB(106) = RF(106) / MAX(EQK, SMALL)
      RF(107) = 1.33333333D0*RF(92)
      EQK = EG(12)*EG(12)/EG(11)/EG(13)
      RB(107) = RF(107) / MAX(EQK, SMALL)
      RF(108) = 1D14
      EQK = EG(11)*EG(14)/EG(2)/EG(25)
      RB(108) = RF(108) / MAX(EQK, SMALL)
      RF(109) = 1D14
      EQK = EG(2)*EG(14)*EG(14)/EG(3)/EG(25)*PFAC1
      RB(109) = RF(109) / MAX(EQK, SMALL)
      RF(110) = EXP(2.81010247D1 -4.29747034D2*TI)
      EQK = EG(5)*EG(14)*EG(14)/EG(4)/EG(25)*PFAC1
      RB(110) = RF(110) / MAX(EQK, SMALL)
      RF(111) = 5D13
      EQK = EG(14)*EG(19)/EG(9)/EG(25)
      RB(111) = RF(111) / MAX(EQK, SMALL)
      RF(112) = 3D13
      EQK = EG(14)*EG(21)/EG(10)/EG(25)
      RB(112) = RF(112) / MAX(EQK, SMALL)
      RF(113) = 1D13
      EQK = EG(14)*EG(14)*EG(19)/EG(25)/EG(25)*PFAC1
      RB(113) = RF(113) / MAX(EQK, SMALL)
      RF(114) = EXP(3.43156328D1 -5.2D-1*ALOGT -2.55382459D4*TI)
      EQK = EG(20)/EG(19)
      RB(114) = RF(114) / MAX(EQK, SMALL)
      RF(115) = EXP(1.97713479D1 +1.62D0*ALOGT -1.86432818D4*TI)
      EQK = EG(2)*EG(19)/EG(21)*PFAC1
      RB(115) = RF(115) / MAX(EQK, SMALL)
      RF(116) = EXP(1.66079019D1 +2D0*ALOGT -9.56111669D2*TI)
      EQK = EG(2)*EG(25)/EG(3)/EG(19)
      RB(116) = RF(116) / MAX(EQK, SMALL)
      RF(117) = 2.5D-1*RF(116)
      EQK = EG(10)*EG(14)/EG(3)/EG(19)
      RB(117) = RF(117) / MAX(EQK, SMALL)
      RF(118) = EXP(-8.4310155D0 +4.5D0*ALOGT +5.03216668D2*TI)
      EQK = EG(2)*EG(26)/EG(5)/EG(19)
      RB(118) = RF(118) / MAX(EQK, SMALL)
      RF(119) = EXP(-7.6354939D0 +4D0*ALOGT +1.00643334D3*TI)
      EQK = EG(12)*EG(14)/EG(5)/EG(19)
      RB(119) = RF(119) / MAX(EQK, SMALL)
      RF(120) = EXP(1.61180957D1 +2D0*ALOGT -3.01930001D3*TI)
      EQK = EG(14)*EG(21)/EG(16)/EG(19)
      RB(120) = RF(120) / MAX(EQK, SMALL)
      RF(121) = EXP(1.27430637D2 -1.182D1*ALOGT -1.79799315D4*TI)
      EQK = EG(29)/EG(12)/EG(19)/PFAC1
      RB(121) = RF(121) / MAX(EQK, SMALL)
      RF(122) = 1D14
      EQK = EG(19)/EG(20)
      RB(122) = RF(122) / MAX(EQK, SMALL)
      RF(123) = 1D14
      EQK = EG(10)*EG(14)/EG(3)/EG(20)
      RB(123) = RF(123) / MAX(EQK, SMALL)
      RF(124) = 2D13
      EQK = EG(2)*EG(26)/EG(5)/EG(20)
      RB(124) = RF(124) / MAX(EQK, SMALL)
      RF(125) = 1D13
      EQK = EG(10)*EG(15)/EG(4)/EG(20)
      RB(125) = RF(125) / MAX(EQK, SMALL)
      RF(126) = EXP(3.34301138D1 -6D-2*ALOGT -4.27734167D3*TI)
      EQK = EG(27)/EG(2)/EG(26)/PFAC1
      RB(126) = RF(126) / MAX(EQK, SMALL)
      RF(127) = 5D1*RF(76)
      EQK = EG(1)*EG(25)/EG(2)/EG(26)
      RB(127) = RF(127) / MAX(EQK, SMALL)
      RF(128) = EXP(2.11287309D1 +1.43D0*ALOGT -1.35365284D3*TI)
      EQK = EG(12)*EG(14)/EG(2)/EG(26)
      RB(128) = RF(128) / MAX(EQK, SMALL)
      RF(129) = 1D1*RF(76)
      EQK = EG(5)*EG(25)/EG(3)/EG(26)
      RB(129) = RF(129) / MAX(EQK, SMALL)
      RF(130) = EXP(2.81906369D1 -6.79342501D2*TI)
      EQK = EG(10)*EG(15)/EG(3)/EG(26)
      RB(130) = RF(130) / MAX(EQK, SMALL)
      RF(131) = EXP(2.96459241D1 -1.00643334D3*TI)
      EQK = EG(6)*EG(25)/EG(5)/EG(26)
      RB(131) = RF(131) / MAX(EQK, SMALL)
      RF(132) = EXP(2.94360258D1 +2.7D-1*ALOGT -1.40900667D2*TI)
      EQK = EG(22)/EG(2)/EG(21)/PFAC1
      RB(132) = RF(132) / MAX(EQK, SMALL)
      RF(133) = 3D13
      EQK = EG(1)*EG(19)/EG(2)/EG(21)
      RB(133) = RF(133) / MAX(EQK, SMALL)
      RF(134) = 6D13
      EQK = EG(1)*EG(20)/EG(2)/EG(21)
      RB(134) = RF(134) / MAX(EQK, SMALL)
      RF(135) = 4.8D13
      EQK = EG(2)*EG(26)/EG(3)/EG(21)
      RB(135) = RF(135) / MAX(EQK, SMALL)
      RF(136) = 4.8D13
      EQK = EG(12)*EG(14)/EG(3)/EG(21)
      RB(136) = RF(136) / MAX(EQK, SMALL)
      RF(137) = 3.011D13
      EQK = EG(6)*EG(19)/EG(5)/EG(21)
      RB(137) = RF(137) / MAX(EQK, SMALL)
      RF(138) = EXP(1.41081802D1 +1.61D0*ALOGT +1.9293327D2*TI)
      EQK = EG(7)*EG(19)/EG(4)/EG(21)
      RB(138) = RF(138) / MAX(EQK, SMALL)
      RF(139) = EXP(2.64270483D1 +2.9D-1*ALOGT -5.53538334D0*TI)
      EQK = EG(3)*EG(27)/EG(4)/EG(21)
      RB(139) = RF(139) / MAX(EQK, SMALL)
      RF(140) = EXP(3.83674178D1 -1.39D0*ALOGT -5.08248834D2*TI)
      EQK = EG(16)*EG(17)/EG(4)/EG(21)
      RB(140) = RF(140) / MAX(EQK, SMALL)
      RF(141) = 1D13
      EQK = EG(5)*EG(27)/EG(7)/EG(21)
      RB(141) = RF(141) / MAX(EQK, SMALL)
      RF(142) = EXP(2.32164713D1 +2.99917134D2*TI)
      EQK = EG(7)*EG(22)/EG(8)/EG(21)
      RB(142) = RF(142) / MAX(EQK, SMALL)
      RF(143) = 9.033D13
      EQK = EG(14)*EG(22)/EG(16)/EG(21)
      RB(143) = RF(143) / MAX(EQK, SMALL)
      RF(144) = 3.92D11
      EQK = EG(13)*EG(19)/EG(12)/EG(21)
      RB(144) = RF(144) / MAX(EQK, SMALL)
      RF(145) = 2.5D13
      EQK = EG(30)/EG(12)/EG(21)/PFAC1
      RB(145) = RF(145) / MAX(EQK, SMALL)
      RF(146) = EXP(5.56675073D1 -2.83D0*ALOGT -9.36888792D3*TI)
      EQK = EG(2)*EG(29)/EG(12)/EG(21)
      RB(146) = RF(146) / MAX(EQK, SMALL)
      RF(147) = EXP(9.64601125D1 -9.147D0*ALOGT -2.36008617D4*TI)
      EQK = EG(12)*EG(14)/EG(27)*PFAC1
      RB(147) = RF(147) / MAX(EQK, SMALL)
      RF(148) = 1D14
      EQK = EG(28)/EG(2)/EG(27)/PFAC1
      RB(148) = RF(148) / MAX(EQK, SMALL)
      RF(149) = 9D13
      EQK = EG(12)*EG(16)/EG(2)/EG(27)
      RB(149) = RF(149) / MAX(EQK, SMALL)
      RF(150) = EXP(3.06267534D1 -2.01286667D3*TI)
      EQK = EG(1)*EG(26)/EG(2)/EG(27)
      RB(150) = RF(150) / MAX(EQK, SMALL)
      RF(151) = RF(150)
      EQK = EG(5)*EG(26)/EG(3)/EG(27)
      RB(151) = RF(151) / MAX(EQK, SMALL)
      RF(152) = 1.33333333D0*RF(131)
      EQK = EG(6)*EG(26)/EG(5)/EG(27)
      RB(152) = RF(152) / MAX(EQK, SMALL)
      RF(153) = 1.4D11
      EQK = EG(7)*EG(26)/EG(4)/EG(27)
      RB(153) = RF(153) / MAX(EQK, SMALL)
      RF(154) = 1.8D10
      EQK = EG(5)*EG(14)*EG(17)/EG(4)/EG(27)*PFAC1
      RB(154) = RF(154) / MAX(EQK, SMALL)
      RF(155) = EXP(2.97104627D1 +4.4D-1*ALOGT -4.46705436D4*TI)
      EQK = EG(1)*EG(20)/EG(22)*PFAC1
      RB(155) = RF(155) / MAX(EQK, SMALL)
      RF(156) = EXP(2.77079822D1 +4.54D-1*ALOGT -9.15854335D2*TI)
      EQK = EG(23)/EG(2)/EG(22)/PFAC1
      RB(156) = RF(156) / MAX(EQK, SMALL)
      RF(157) = EXP(1.77414365D1 +1.93D0*ALOGT -6.51665585D3*TI)
      EQK = EG(1)*EG(21)/EG(2)/EG(22)
      RB(157) = RF(157) / MAX(EQK, SMALL)
      RF(158) = EXP(1.65302053D1 +1.91D0*ALOGT -1.88203034D3*TI)
      EQK = EG(5)*EG(21)/EG(3)/EG(22)
      RB(158) = RF(158) / MAX(EQK, SMALL)
      RF(159) = EXP(1.67704208D1 +1.83D0*ALOGT -1.10707667D2*TI)
      EQK = EG(12)*EG(16)/EG(3)/EG(22)
      RB(159) = RF(159) / MAX(EQK, SMALL)
      RF(160) = 2D-2*RF(159)
      EQK = EG(10)*EG(17)/EG(3)/EG(22)
      RB(160) = RF(160) / MAX(EQK, SMALL)
      RF(161) = EXP(1.50964444D1 +2D0*ALOGT -1.25804167D3*TI)
      EQK = EG(6)*EG(21)/EG(5)/EG(22)
      RB(161) = RF(161) / MAX(EQK, SMALL)
      RF(162) = EXP(3.13734413D1 -3.05955734D4*TI)
      EQK = EG(7)*EG(21)/EG(4)/EG(22)
      RB(162) = RF(162) / MAX(EQK, SMALL)
      RF(163) = EXP(2.83241683D1 -7.04503335D3*TI)
      EQK = EG(5)*EG(28)/EG(7)/EG(22)
      RB(163) = RF(163) / MAX(EQK, SMALL)
      RF(164) = EXP(1.61180957D1 +2D0*ALOGT -4.02573334D3*TI)
      EQK = EG(14)*EG(23)/EG(16)/EG(22)
      RB(164) = RF(164) / MAX(EQK, SMALL)
      RF(165) = EXP(3.06267534D1 -3.01930001D3*TI)
      EQK = EG(2)*EG(29)/EG(10)/EG(22)
      RB(165) = RF(165) / MAX(EQK, SMALL)
      RF(166) = 5D13
      EQK = EG(13)*EG(20)/EG(11)/EG(22)
      RB(166) = RF(166) / MAX(EQK, SMALL)
      RF(167) = 5D13
      EQK = EG(2)*EG(29)/EG(11)/EG(22)
      RB(167) = RF(167) / MAX(EQK, SMALL)
      RF(168) = EXP(1.23327053D1 +2D0*ALOGT -4.62959334D3*TI)
      EQK = EG(13)*EG(21)/EG(12)/EG(22)
      RB(168) = RF(168) / MAX(EQK, SMALL)
      RF(169) = EXP(2.65223585D1 -3.87476834D3*TI)
      EQK = EG(31)/EG(12)/EG(22)/PFAC1
      RB(169) = RF(169) / MAX(EQK, SMALL)
      RF(170) = EXP(4.07945264D1 -9.9D-1*ALOGT -7.95082335D2*TI)
      EQK = EG(24)/EG(2)/EG(23)/PFAC1
      RB(170) = RF(170) / MAX(EQK, SMALL)
      RF(171) = 2D12
      EQK = EG(1)*EG(22)/EG(2)/EG(23)
      RB(171) = RF(171) / MAX(EQK, SMALL)
      RF(172) = 1.604D13
      EQK = EG(12)*EG(17)/EG(3)/EG(23)
      RB(172) = RF(172) / MAX(EQK, SMALL)
      RF(173) = 8.02D13
      EQK = EG(2)*EG(28)/EG(3)/EG(23)
      RB(173) = RF(173) / MAX(EQK, SMALL)
      RF(174) = 2D10
      EQK = EG(7)*EG(22)/EG(4)/EG(23)
      RB(174) = RF(174) / MAX(EQK, SMALL)
      RF(175) = 3D11
      EQK = EG(4)*EG(24)/EG(7)/EG(23)
      RB(175) = RF(175) / MAX(EQK, SMALL)
      RF(176) = 3D11
      EQK = EG(8)*EG(22)/EG(7)/EG(23)
      RB(176) = RF(176) / MAX(EQK, SMALL)
      RF(177) = 2.4D13
      EQK = EG(5)*EG(12)*EG(17)/EG(7)/EG(23)*PFAC1
      RB(177) = RF(177) / MAX(EQK, SMALL)
      RF(178) = EXP(2.28865889D1 -4.90133034D2*TI)
      EQK = EG(7)*EG(24)/EG(8)/EG(23)
      RB(178) = RF(178) / MAX(EQK, SMALL)
      RF(179) = 1.2D14
      EQK = EG(14)*EG(24)/EG(16)/EG(23)
      RB(179) = RF(179) / MAX(EQK, SMALL)
      RF(180) = EXP(1.85604427D1 +1.9D0*ALOGT -3.78922151D3*TI)
      EQK = EG(1)*EG(23)/EG(2)/EG(24)
      RB(180) = RF(180) / MAX(EQK, SMALL)
      RF(181) = EXP(1.83130955D1 +1.92D0*ALOGT -2.86330284D3*TI)
      EQK = EG(5)*EG(23)/EG(3)/EG(24)
      RB(181) = RF(181) / MAX(EQK, SMALL)
      RF(182) = EXP(1.50796373D1 +2.12D0*ALOGT -4.37798501D2*TI)
      EQK = EG(6)*EG(23)/EG(5)/EG(24)
      RB(182) = RF(182) / MAX(EQK, SMALL)
      RF(183) = EXP(3.13199006D1 +2.76769167D2*TI)
      EQK = EG(12)*EG(23)/EG(11)/EG(24)
      RB(183) = RF(183) / MAX(EQK, SMALL)
      RF(184) = EXP(1.56303353D1 +1.74D0*ALOGT -5.25861418D3*TI)
      EQK = EG(13)*EG(23)/EG(12)/EG(24)
      RB(184) = RF(184) / MAX(EQK, SMALL)
      RF(185) = 2D14
      EQK = EG(30)/EG(2)/EG(29)/PFAC1
      RB(185) = RF(185) / MAX(EQK, SMALL)
      RF(186) = 2.66666667D0*RF(131)
      EQK = EG(13)*EG(20)/EG(2)/EG(29)
      RB(186) = RF(186) / MAX(EQK, SMALL)
      RF(187) = 2.66D12
      EQK = EG(4)*EG(30)/EG(7)/EG(29)
      RB(187) = RF(187) / MAX(EQK, SMALL)
      RF(188) = 6.6D12
      EQK = EG(5)*EG(17)*EG(21)/EG(7)/EG(29)*PFAC1
      RB(188) = RF(188) / MAX(EQK, SMALL)
      RF(189) = 6D13
      EQK = EG(14)*EG(30)/EG(16)/EG(29)
      RB(189) = RF(189) / MAX(EQK, SMALL)
      RF(190) = EXP(3.02187852D1 -1.64083859D3*TI)
      EQK = EG(31)/EG(2)/EG(30)/PFAC1
      RB(190) = RF(190) / MAX(EQK, SMALL)
      RF(191) = EXP(5.11268757D1 -2.39D0*ALOGT -5.62596234D3*TI)
      EQK = EG(12)*EG(22)/EG(2)/EG(30)
      RB(191) = RF(191) / MAX(EQK, SMALL)
      RF(192) = EXP(1.20435537D1 +2.5D0*ALOGT -1.2530095D3*TI)
      EQK = EG(1)*EG(29)/EG(2)/EG(30)
      RB(192) = RF(192) / MAX(EQK, SMALL)
      RF(193) = EXP(1.86030023D1 +1.65D0*ALOGT -1.6455185D2*TI)
      EQK = EG(2)*EG(12)*EG(26)/EG(3)/EG(30)*PFAC1
      RB(193) = RF(193) / MAX(EQK, SMALL)
      RF(194) = EXP(1.73708586D1 +1.65D0*ALOGT +4.89126601D2*TI)
      EQK = EG(16)*EG(23)/EG(3)/EG(30)
      RB(194) = RF(194) / MAX(EQK, SMALL)
      RF(195) = EXP(2.59162227D1 +7D-1*ALOGT -2.95891401D3*TI)
      EQK = EG(5)*EG(29)/EG(3)/EG(30)
      RB(195) = RF(195) / MAX(EQK, SMALL)
      RF(196) = EXP(1.49469127D1 +2D0*ALOGT +1.49958567D2*TI)
      EQK = EG(6)*EG(29)/EG(5)/EG(30)
      RB(196) = RF(196) / MAX(EQK, SMALL)
      RF(197) = EXP(9.16951838D0 +2.6D0*ALOGT -6.99974385D3*TI)
      EQK = EG(8)*EG(29)/EG(7)/EG(30)
      RB(197) = RF(197) / MAX(EQK, SMALL)
      RF(198) = EXP(7.8845736D-1 +3.5D0*ALOGT -2.85575459D3*TI)
      EQK = EG(13)*EG(29)/EG(12)/EG(30)
      RB(198) = RF(198) / MAX(EQK, SMALL)
      RF(199) = EXP(5.65703751D1 -2.92D0*ALOGT -6.29272443D3*TI)
      EQK = EG(12)*EG(23)/EG(2)/EG(31)
      RB(199) = RF(199) / MAX(EQK, SMALL)
      RF(200) = 1.8D12
      EQK = EG(1)*EG(30)/EG(2)/EG(31)
      RB(200) = RF(200) / MAX(EQK, SMALL)
      RF(201) = 9.6D13
      EQK = EG(17)*EG(23)/EG(3)/EG(31)
      RB(201) = RF(201) / MAX(EQK, SMALL)
      RF(202) = 2.4D13
      EQK = EG(6)*EG(30)/EG(5)/EG(31)
      RB(202) = RF(202) / MAX(EQK, SMALL)
      RF(203) = 9D10
      EQK = EG(7)*EG(30)/EG(4)/EG(31)
      RB(203) = RF(203) / MAX(EQK, SMALL)
      RF(204) = 2.4D13
      EQK = EG(5)*EG(17)*EG(23)/EG(7)/EG(31)*PFAC1
      RB(204) = RF(204) / MAX(EQK, SMALL)
      RF(205) = 1.1D13
      EQK = EG(13)*EG(30)/EG(12)/EG(31)
      RB(205) = RF(205) / MAX(EQK, SMALL)
      RF(206) = EXP(7.50436995D1 -5.22D0*ALOGT -9.93701954D3*TI)
      EQK = EG(12)*EG(29)/EG(21)/EG(23)
      RB(206) = RF(206) / MAX(EQK, SMALL)
C
      RKLOW(1) = EXP(4.22794408D1 -9D-1*ALOGT +8.55468335D2*TI)
      RKLOW(2) = EXP(6.37931383D1 -3.42D0*ALOGT -4.24463259D4*TI)
      RKLOW(3) = EXP(6.54619238D1 -3.74D0*ALOGT -9.74227469D2*TI)
      RKLOW(4) = EXP(5.55621468D1 -2.57D0*ALOGT -7.17083751D2*TI)
      RKLOW(5) = EXP(6.33329483D1 -3.14D0*ALOGT -6.18956501D2*TI)
      RKLOW(6) = EXP(7.69748493D1 -5.11D0*ALOGT -3.57032226D3*TI)
      RKLOW(7) = EXP(6.98660102D1 -4.8D0*ALOGT -2.79788467D3*TI)
      RKLOW(8) = EXP(7.68923562D1 -4.76D0*ALOGT -1.22784867D3*TI)
      RKLOW(9) = EXP(1.11312542D2 -9.588D0*ALOGT -2.566405D3*TI)
      RKLOW(10) = EXP(1.15700234D2 -9.67D0*ALOGT -3.13000767D3*TI)
      RKLOW(11) = EXP(3.54348644D1 -6.4D-1*ALOGT -2.50098684D4*TI)
      RKLOW(12) = EXP(6.3111756D1 -3.4D0*ALOGT -1.80145126D4*TI)
      RKLOW(13) = EXP(9.57409899D1 -7.64D0*ALOGT -5.98827834D3*TI)
      RKLOW(14) = EXP(6.9414025D1 -3.86D0*ALOGT -1.67067934D3*TI)
      RKLOW(15) = EXP(1.35001549D2 -1.194D1*ALOGT -4.9163262D3*TI)
      RKLOW(16) = EXP(9.14494773D1 -7.297D0*ALOGT -2.36511834D3*TI)
      RKLOW(17) = EXP(1.17075165D2 -9.31D0*ALOGT -5.02512164D4*TI)
      RKLOW(18) = EXP(9.68908955D1 -7.62D0*ALOGT -3.50742017D3*TI)
      RKLOW(19) = EXP(9.50941235D1 -7.08D0*ALOGT -3.36400342D3*TI)
      RKLOW(20) = EXP(1.38440285D2 -1.2D1*ALOGT -3.00309643D3*TI)
      RKLOW(21) = EXP(8.93324137D1 -6.66D0*ALOGT -3.52251667D3*TI)
C
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE RDSMH  (T, SMH)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION SMH(*), TN(5)
C
      TLOG = LOG(T)
      TI = 1D0/T
C
      TN(1) = TLOG - 1D0
      TN(2) = T
      TN(3) = TN(2)*T
      TN(4) = TN(3)*T
      TN(5) = TN(4)*T
C H2
      IF (T .GT. 1000) THEN
      SMH(1) = -3.20502331D0 +9.50158922D2*TI 
     *         +3.3372792D0*TN(1) -2.47012366D-5*TN(2) 
     *         +8.32427963D-8*TN(3) -1.49638662D-11*TN(4) 
     *         +1.00127688D-15*TN(5) 
      ELSE
      SMH(1) = 6.83010238D-1 +9.17935173D2*TI 
     *         +2.34433112D0*TN(1) +3.99026037D-3*TN(2) 
     *         -3.2463585D-6*TN(3) +1.67976745D-9*TN(4) 
     *         -3.68805881D-13*TN(5) 
      ENDIF
C H
      IF (T .GT. 1000) THEN
      SMH(2) = -4.46682914D-1 -2.54736599D4*TI 
     *         +2.50000001D0*TN(1) -1.15421486D-11*TN(2) 
     *         +2.69269913D-15*TN(3) -3.94596029D-19*TN(4) 
     *         +2.49098679D-23*TN(5) 
      ELSE
      SMH(2) = -4.46682853D-1 -2.54736599D4*TI 
     *         +2.5D0*TN(1) +3.5266641D-13*TN(2) 
     *         -3.32653273D-16*TN(3) +1.91734693D-19*TN(4) 
     *         -4.63866166D-23*TN(5) 
      ENDIF
C O
      IF (T .GT. 1000) THEN
      SMH(3) = 4.78433864D0 -2.92175791D4*TI 
     *         +2.56942078D0*TN(1) -4.29870568D-5*TN(2) 
     *         +6.99140982D-9*TN(3) -8.34814992D-13*TN(4) 
     *         +6.14168455D-17*TN(5) 
      ELSE
      SMH(3) = 2.05193346D0 -2.91222592D4*TI 
     *         +3.1682671D0*TN(1) -1.63965942D-3*TN(2) 
     *         +1.10717733D-6*TN(3) -5.10672187D-10*TN(4) 
     *         +1.05632986D-13*TN(5) 
      ENDIF
C O2
      IF (T .GT. 1000) THEN
      SMH(4) = 5.45323129D0 +1.08845772D3*TI 
     *         +3.28253784D0*TN(1) +7.4154377D-4*TN(2) 
     *         -1.26327778D-7*TN(3) +1.74558796D-11*TN(4) 
     *         -1.08358897D-15*TN(5) 
      ELSE
      SMH(4) = 3.65767573D0 +1.06394356D3*TI 
     *         +3.78245636D0*TN(1) -1.49836708D-3*TN(2) 
     *         +1.641217D-6*TN(3) -8.06774591D-10*TN(4) 
     *         +1.62186419D-13*TN(5) 
      ENDIF
C OH
      IF (T .GT. 1000) THEN
      SMH(5) = 4.4766961D0 -3.858657D3*TI 
     *         +3.09288767D0*TN(1) +2.74214858D-4*TN(2) 
     *         +2.10842047D-8*TN(3) -7.3288463D-12*TN(4) 
     *         +5.8706188D-16*TN(5) 
      ELSE
      SMH(5) = -1.03925458D-1 -3.61508056D3*TI 
     *         +3.99201543D0*TN(1) -1.20065876D-3*TN(2) 
     *         +7.69656402D-7*TN(3) -3.23427778D-10*TN(4) 
     *         +6.8205735D-14*TN(5) 
      ENDIF
C H2O
      IF (T .GT. 1000) THEN
      SMH(6) = 4.9667701D0 +3.00042971D4*TI 
     *         +3.03399249D0*TN(1) +1.08845902D-3*TN(2) 
     *         -2.73454197D-8*TN(3) -8.08683225D-12*TN(4) 
     *         +8.4100496D-16*TN(5) 
      ELSE
      SMH(6) = -8.49032208D-1 +3.02937267D4*TI 
     *         +4.19864056D0*TN(1) -1.01821705D-3*TN(2) 
     *         +1.08673369D-6*TN(3) -4.57330885D-10*TN(4) 
     *         +8.85989085D-14*TN(5) 
      ENDIF
C HO2
      IF (T .GT. 1000) THEN
      SMH(7) = 3.78510215D0 -1.11856713D2*TI 
     *         +4.0172109D0*TN(1) +1.11991007D-3*TN(2) 
     *         -1.05609692D-7*TN(3) +9.52053083D-12*TN(4) 
     *         -5.39542675D-16*TN(5) 
      ELSE
      SMH(7) = 3.71666245D0 -2.9480804D2*TI 
     *         +4.30179801D0*TN(1) -2.37456026D-3*TN(2) 
     *         +3.52638152D-6*TN(3) -2.02303245D-9*TN(4) 
     *         +4.64612562D-13*TN(5) 
      ENDIF
C H2O2
      IF (T .GT. 1000) THEN
      SMH(8) = 2.91615662D0 +1.78617877D4*TI 
     *         +4.16500285D0*TN(1) +2.45415847D-3*TN(2) 
     *         -3.16898708D-7*TN(3) +3.09321655D-11*TN(4) 
     *         -1.43954153D-15*TN(5) 
      ELSE
      SMH(8) = 3.43505074D0 +1.77025821D4*TI 
     *         +4.27611269D0*TN(1) -2.71411208D-4*TN(2) 
     *         +2.78892835D-6*TN(3) -1.79809011D-9*TN(4) 
     *         +4.31227181D-13*TN(5) 
      ENDIF
C CH
      IF (T .GT. 1000) THEN
      SMH(9) = 5.48497999D0 -7.10124364D4*TI 
     *         +2.87846473D0*TN(1) +4.85456841D-4*TN(2) 
     *         +2.40742758D-8*TN(3) -1.08906541D-11*TN(4) 
     *         +8.80396915D-16*TN(5) 
      ELSE
      SMH(9) = 2.08401108D0 -7.07972934D4*TI 
     *         +3.48981665D0*TN(1) +1.61917771D-4*TN(2) 
     *         -2.81498442D-7*TN(3) +2.63514439D-10*TN(4) 
     *         -7.03045335D-14*TN(5) 
      ENDIF
C CH2
      IF (T .GT. 1000) THEN
      SMH(10) = 6.17119324D0 -4.6263604D4*TI 
     *         +2.87410113D0*TN(1) +1.82819646D-3*TN(2) 
     *         -2.34824328D-7*TN(3) +2.16816291D-11*TN(4) 
     *         -9.38637835D-16*TN(5) 
      ELSE
      SMH(10) = 1.56253185D0 -4.60040401D4*TI 
     *         +3.76267867D0*TN(1) +4.84436072D-4*TN(2) 
     *         +4.65816402D-7*TN(3) -3.20909294D-10*TN(4) 
     *         +8.43708595D-14*TN(5) 
      ENDIF
C CH2*
      IF (T .GT. 1000) THEN
      SMH(11) = 8.62650169D0 -5.09259997D4*TI 
     *         +2.29203842D0*TN(1) +2.32794319D-3*TN(2) 
     *         -3.35319912D-7*TN(3) +3.48255D-11*TN(4) 
     *         -1.69858183D-15*TN(5) 
      ELSE
      SMH(11) = -7.69118967D-1 -5.04968163D4*TI 
     *         +4.19860411D0*TN(1) -1.1833071D-3*TN(2) 
     *         +1.37216037D-6*TN(3) -5.57346651D-10*TN(4) 
     *         +9.71573685D-14*TN(5) 
      ENDIF
C CH3
      IF (T .GT. 1000) THEN
      SMH(12) = 8.48007179D0 -1.67755843D4*TI 
     *         +2.28571772D0*TN(1) +3.61995019D-3*TN(2) 
     *         -4.97857247D-7*TN(3) +4.9640387D-11*TN(4) 
     *         -2.33577197D-15*TN(5) 
      ELSE
      SMH(12) = 1.60456433D0 -1.64449988D4*TI 
     *         +3.6735904D0*TN(1) +1.00547587D-3*TN(2) 
     *         +9.55036427D-7*TN(3) -5.72597854D-10*TN(4) 
     *         +1.27192867D-13*TN(5) 
      ENDIF
C CH4
      IF (T .GT. 1000) THEN
      SMH(13) = 1.8437318D1 +9.46834459D3*TI 
     *         +7.4851495D-2*TN(1) +6.69547335D-3*TN(2) 
     *         -9.55476348D-7*TN(3) +1.01910446D-10*TN(4) 
     *         -5.0907615D-15*TN(5) 
      ELSE
      SMH(13) = -4.64130376D0 +1.02466476D4*TI 
     *         +5.14987613D0*TN(1) -6.8354894D-3*TN(2) 
     *         +8.19667665D-6*TN(3) -4.03952522D-9*TN(4) 
     *         +8.3346978D-13*TN(5) 
      ENDIF
C CO
      IF (T .GT. 1000) THEN
      SMH(14) = 7.81868772D0 +1.41518724D4*TI 
     *         +2.71518561D0*TN(1) +1.03126372D-3*TN(2) 
     *         -1.66470962D-7*TN(3) +1.9171084D-11*TN(4) 
     *         -1.01823858D-15*TN(5) 
      ELSE
      SMH(14) = 3.50840928D0 +1.4344086D4*TI 
     *         +3.57953347D0*TN(1) -3.0517684D-4*TN(2) 
     *         +1.69469055D-7*TN(3) +7.55838237D-11*TN(4) 
     *         -4.52212249D-14*TN(5) 
      ENDIF
C CO2
      IF (T .GT. 1000) THEN
      SMH(15) = 2.27163806D0 +4.8759166D4*TI 
     *         +3.85746029D0*TN(1) +2.20718513D-3*TN(2) 
     *         -3.69135673D-7*TN(3) +4.36241823D-11*TN(4) 
     *         -2.36042082D-15*TN(5) 
      ELSE
      SMH(15) = 9.90105222D0 +4.83719697D4*TI 
     *         +2.35677352D0*TN(1) +4.49229838D-3*TN(2) 
     *         -1.18726045D-6*TN(3) +2.04932518D-10*TN(4) 
     *         -7.1849774D-15*TN(5) 
      ENDIF
C HCO
      IF (T .GT. 1000) THEN
      SMH(16) = 9.79834492D0 -4.01191815D3*TI 
     *         +2.77217438D0*TN(1) +2.47847763D-3*TN(2) 
     *         -4.14076022D-7*TN(3) +4.90968148D-11*TN(4) 
     *         -2.66754356D-15*TN(5) 
      ELSE
      SMH(16) = 3.39437243D0 -3.83956496D3*TI 
     *         +4.22118584D0*TN(1) -1.62196266D-3*TN(2) 
     *         +2.29665743D-6*TN(3) -1.10953411D-9*TN(4) 
     *         +2.16884432D-13*TN(5) 
      ENDIF
C CH2O
      IF (T .GT. 1000) THEN
      SMH(17) = 1.3656323D1 +1.39958323D4*TI 
     *         +1.76069008D0*TN(1) +4.60000041D-3*TN(2) 
     *         -7.37098022D-7*TN(3) +8.38676767D-11*TN(4) 
     *         -4.4192782D-15*TN(5) 
      ELSE
      SMH(17) = 6.028129D-1 +1.43089567D4*TI 
     *         +4.79372315D0*TN(1) -4.95416684D-3*TN(2) 
     *         +6.22033347D-6*TN(3) -3.16071051D-9*TN(4) 
     *         +6.5886326D-13*TN(5) 
      ENDIF
C CH3O
      IF (T .GT. 1000) THEN
      SMH(18) = 2.929575D0 -1.2783252D2*TI 
     *         +3.770799D0*TN(1) +3.9357485D-3*TN(2) 
     *         -4.42730667D-7*TN(3) +3.28702583D-11*TN(4) 
     *         -1.056308D-15*TN(5) 
      ELSE
      SMH(18) = 1.3152177D1 -9.786011D2*TI 
     *         +2.106204D0*TN(1) +3.6082975D-3*TN(2) 
     *         +8.89745333D-7*TN(3) -6.14803D-10*TN(4) 
     *         +1.037805D-13*TN(5) 
      ENDIF
C C2H2
      IF (T .GT. 1000) THEN
      SMH(19) = -1.23028121D0 -2.59359992D4*TI 
     *         +4.14756964D0*TN(1) +2.98083332D-3*TN(2) 
     *         -3.9549142D-7*TN(3) +3.89510143D-11*TN(4) 
     *         -1.80617607D-15*TN(5) 
      ELSE
      SMH(19) = 1.39397051D1 -2.64289807D4*TI 
     *         +8.08681094D-1*TN(1) +1.16807815D-2*TN(2) 
     *         -5.91953025D-6*TN(3) +2.33460364D-9*TN(4) 
     *         -4.25036487D-13*TN(5) 
      ENDIF
C H2CC
      IF (T .GT. 1000) THEN
      SMH(20) = 6.4023701D-1 -4.8316688D4*TI 
     *         +4.278034D0*TN(1) +2.3781402D-3*TN(2) 
     *         -2.71683483D-7*TN(3) +2.1219005D-11*TN(4) 
     *         -7.4431895D-16*TN(5) 
      ELSE
      SMH(20) = 5.920391D0 -4.8621794D4*TI 
     *         +3.2815483D0*TN(1) +3.48823955D-3*TN(2) 
     *         -3.975874D-7*TN(3) -1.00870267D-10*TN(4) 
     *         +4.90947725D-14*TN(5) 
      ENDIF
C C2H3
      IF (T .GT. 1000) THEN
      SMH(21) = 7.78732378D0 -3.46128739D4*TI 
     *         +3.016724D0*TN(1) +5.1651146D-3*TN(2) 
     *         -7.80137248D-7*TN(3) +8.480274D-11*TN(4) 
     *         -4.3130352D-15*TN(5) 
      ELSE
      SMH(21) = 8.51054025D0 -3.48598468D4*TI 
     *         +3.21246645D0*TN(1) +7.5739581D-4*TN(2) 
     *         +4.32015687D-6*TN(3) -2.98048206D-9*TN(4) 
     *         +7.35754365D-13*TN(5) 
      ENDIF
C C2H4
      IF (T .GT. 1000) THEN
      SMH(22) = 1.03053693D1 -4.93988614D3*TI 
     *         +2.03611116D0*TN(1) +7.32270755D-3*TN(2) 
     *         -1.11846319D-6*TN(3) +1.22685769D-10*TN(4) 
     *         -6.28530305D-15*TN(5) 
      ELSE
      SMH(22) = 4.09733096D0 -5.08977593D3*TI 
     *         +3.95920148D0*TN(1) -3.78526124D-3*TN(2) 
     *         +9.51650487D-6*TN(3) -5.76323961D-9*TN(4) 
     *         +1.34942187D-12*TN(5) 
      ENDIF
C C2H5
      IF (T .GT. 1000) THEN
      SMH(23) = 1.34624343D1 -1.285752D4*TI 
     *         +1.95465642D0*TN(1) +8.6986361D-3*TN(2) 
     *         -1.33034445D-6*TN(3) +1.46014741D-10*TN(4) 
     *         -7.4820788D-15*TN(5) 
      ELSE
      SMH(23) = 4.70720924D0 -1.28416265D4*TI 
     *         +4.30646568D0*TN(1) -2.09329446D-3*TN(2) 
     *         +8.28571345D-6*TN(3) -4.99272172D-9*TN(4) 
     *         +1.15254502D-12*TN(5) 
      ENDIF
C C2H6
      IF (T .GT. 1000) THEN
      SMH(24) = 1.51156107D1 +1.14263932D4*TI 
     *         +1.0718815D0*TN(1) +1.08426339D-2*TN(2) 
     *         -1.67093445D-6*TN(3) +1.84510001D-10*TN(4) 
     *         -9.5001445D-15*TN(5) 
      ELSE
      SMH(24) = 2.66682316D0 +1.15222055D4*TI 
     *         +4.29142492D0*TN(1) -2.75077135D-3*TN(2) 
     *         +9.99063813D-6*TN(3) -5.90388571D-9*TN(4) 
     *         +1.34342886D-12*TN(5) 
      ENDIF
C HCCO
      IF (T .GT. 1000) THEN
      SMH(25) = -3.9302595D0 -1.9327215D4*TI 
     *         +5.6282058D0*TN(1) +2.04267005D-3*TN(2) 
     *         -2.65575783D-7*TN(3) +2.38550433D-11*TN(4) 
     *         -9.703916D-16*TN(5) 
      ELSE
      SMH(25) = 1.2490417D1 -2.0059449D4*TI 
     *         +2.2517214D0*TN(1) +8.8275105D-3*TN(2) 
     *         -3.95485017D-6*TN(3) +1.43964658D-9*TN(4) 
     *         -2.53324055D-13*TN(5) 
      ENDIF
C CH2CO
      IF (T .GT. 1000) THEN
      SMH(26) = 6.32247205D-1 +7.55105311D3*TI 
     *         +4.51129732D0*TN(1) +4.50179872D-3*TN(2) 
     *         -6.94899392D-7*TN(3) +7.69454902D-11*TN(4) 
     *         -3.974191D-15*TN(5) 
      ELSE
      SMH(26) = 1.2215648D1 +7.04291804D3*TI 
     *         +2.1358363D0*TN(1) +9.05943605D-3*TN(2) 
     *         -2.89912457D-6*TN(3) +7.7866464D-10*TN(4) 
     *         -1.00728807D-13*TN(5) 
      ENDIF
C CH2CHO
      IF (T .GT. 1000) THEN
      SMH(27) = -5.0320879D0 -4.9032178D2*TI 
     *         +5.9756699D0*TN(1) +4.0652957D-3*TN(2) 
     *         -4.5727075D-7*TN(3) +3.39192008D-11*TN(4) 
     *         -1.08800855D-15*TN(5) 
      ELSE
      SMH(27) = 9.5714535D0 -1.5214766D3*TI 
     *         +3.4090624D0*TN(1) +5.369287D-3*TN(2) 
     *         +3.1524875D-7*TN(3) +5.96548592D-10*TN(4) 
     *         +1.43369255D-13*TN(5) 
      ENDIF
C CH3CHO
      IF (T .GT. 1000) THEN
      SMH(28) = -3.4807917D0 +2.2593122D4*TI 
     *         +5.4041108D0*TN(1) +5.8615295D-3*TN(2) 
     *         -7.04385617D-7*TN(3) +5.69770425D-11*TN(4) 
     *         -2.04924315D-15*TN(5) 
      ELSE
      SMH(28) = 4.1030159D0 +2.1572878D4*TI 
     *         +4.7294595D0*TN(1) -1.5966429D-3*TN(2) 
     *         +7.92248683D-6*TN(3) -4.78821758D-9*TN(4) 
     *         +1.0965556D-12*TN(5) 
      ENDIF
C aC3H5
      IF (T .GT. 1000) THEN
      SMH(29) = -1.124305D1 -1.7482449D4*TI 
     *         +6.5007877D0*TN(1) +7.1623655D-3*TN(2) 
     *         -9.46360533D-7*TN(3) +9.23400083D-11*TN(4) 
     *         -4.51819435D-15*TN(5) 
      ELSE
      SMH(29) = 1.7173214D1 -1.9245629D4*TI 
     *         +1.3631835D0*TN(1) +9.9069105D-3*TN(2) 
     *         +2.08284333D-6*TN(3) -2.77962958D-9*TN(4) 
     *         +7.9232855D-13*TN(5) 
      ENDIF
C C3H6
      IF (T .GT. 1000) THEN
      SMH(30) = -1.331335D1 +9.235703D2*TI 
     *         +6.732257D0*TN(1) +7.45417D-3*TN(2) 
     *         -8.24983167D-7*TN(3) +6.01001833D-11*TN(4) 
     *         -1.883102D-15*TN(5) 
      ELSE
      SMH(30) = 1.614534D1 -1.074826D3*TI 
     *         +1.493307D0*TN(1) +1.046259D-2*TN(2) 
     *         +7.47799D-7*TN(3) -1.39076D-9*TN(4) 
     *         +3.579073D-13*TN(5) 
      ENDIF
C nC3H7
      IF (T .GT. 1000) THEN
      SMH(31) = -1.5515297D1 -7.9762236D3*TI 
     *         +7.7097479D0*TN(1) +8.0157425D-3*TN(2) 
     *         -8.78670633D-7*TN(3) +6.32402933D-11*TN(4) 
     *         -1.94313595D-15*TN(5) 
      ELSE
      SMH(31) = 2.1136034D1 -1.0312346D4*TI 
     *         +1.0491173D0*TN(1) +1.30044865D-2*TN(2) 
     *         +3.92375267D-7*TN(3) -1.63292767D-9*TN(4) 
     *         +4.68601035D-13*TN(5) 
      ENDIF
C
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE RATX (T, C, RF, RB, RKLOW)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      PARAMETER (SMALL = 1D-200)
      DIMENSION C(*), RF(*), RB(*), RKLOW(*)
C
      ALOGT = LOG(T)
      CTOT = 0.0
      DO K = 1, 22
         CTOT = CTOT + C(K)
      ENDDO
C
      RF(1) = RF(1)*C(2)*C(4)
      RB(1) = RB(1)*C(3)*C(5)
      RF(2) = RF(2)*C(1)*C(3)
      RB(2) = RB(2)*C(2)*C(5)
      RF(3) = RF(3)*C(1)*C(5)
      RB(3) = RB(3)*C(2)*C(6)
      RF(4) = RF(4)*C(5)*C(5)
      RB(4) = RB(4)*C(3)*C(6)
      CTB = CTOT-C(1)-C(6)+C(10)-C(12)
     * +2D0*C(16)+2D0*C(14)+2D0*C(15)
      RF(5) = RF(5)*CTB*C(2)*C(2)
      RB(5) = RB(5)*CTB*C(1)
      RF(6) = RF(6)*C(1)*C(2)*C(2)
      RB(6) = RB(6)*C(1)*C(1)
      RF(7) = RF(7)*C(2)*C(2)*C(6)
      RB(7) = RB(7)*C(1)*C(6)
      RF(8) = RF(8)*C(2)*C(2)*C(12)
      RB(8) = RB(8)*C(1)*C(12)
      CTB = CTOT-2.7D-1*C(1)+2.65D0*C(6)+C(10)+2D0*C(16)
     * +2D0*C(14)+2D0*C(15)
      RF(9) = RF(9)*CTB*C(2)*C(5)
      RB(9) = RB(9)*CTB*C(6)
      CTB = CTOT+C(1)+5D0*C(6)+C(10)+5D-1*C(11)
     * +C(12)+2D0*C(16)+2D0*C(14)+2D0*C(15)
      RF(10) = RF(10)*CTB*C(2)*C(3)
      RB(10) = RB(10)*CTB*C(5)
      CTB = CTOT+1.4D0*C(1)+1.44D1*C(6)+C(10)+7.5D-1*C(11)
     * +2.6D0*C(12)+2D0*C(16)+2D0*C(14)+2D0*C(15)
      RF(11) = RF(11)*CTB*C(3)*C(3)
      RB(11) = RB(11)*CTB*C(4)
      CTB = CTOT-C(4)-C(6)-2.5D-1*C(11)+5D-1*C(12)
     * +5D-1*C(16)-C(22)+2D0*C(14)+2D0*C(15)
      RF(12) = RF(12)*CTB*C(2)*C(4)
      RB(12) = RB(12)*CTB*C(7)
      RF(13) = RF(13)*C(2)*C(4)*C(4)
      RB(13) = RB(13)*C(4)*C(7)
      RF(14) = RF(14)*C(2)*C(4)*C(6)
      RB(14) = RB(14)*C(6)*C(7)
      RF(15) = RF(15)*C(2)*C(4)*C(22)
      RB(15) = RB(15)*C(7)*C(22)
      CTB = CTOT+C(1)+5D0*C(6)+C(10)+5D-1*C(11)
     * +C(12)+2D0*C(16)+2D0*C(14)+2D0*C(15)
      PR = RKLOW(1) * CTB / RF(16)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 2.654D-1*EXP(-T/9.4D1) + 7.346D-1*EXP(-T/1.756D3)
     *     + EXP(-5.182D3/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(16) = RF(16) * PCOR
      RB(16) = RB(16) * PCOR
      RF(16) = RF(16)*C(5)*C(5)
      RB(16) = RB(16)*C(8)
      RF(17) = RF(17)*C(2)*C(7)
      RB(17) = RB(17)*C(3)*C(6)
      RF(18) = RF(18)*C(2)*C(7)
      RB(18) = RB(18)*C(1)*C(4)
      RF(19) = RF(19)*C(2)*C(7)
      RB(19) = RB(19)*C(5)*C(5)
      RF(20) = RF(20)*C(3)*C(7)
      RB(20) = RB(20)*C(4)*C(5)
      RF(21) = RF(21)*C(5)*C(7)
      RB(21) = RB(21)*C(4)*C(6)
      RF(22) = RF(22)*C(7)*C(7)
      RB(22) = RB(22)*C(4)*C(8)
      RF(23) = RF(23)*C(7)*C(7)
      RB(23) = RB(23)*C(4)*C(8)
      RF(24) = RF(24)*C(2)*C(8)
      RB(24) = RB(24)*C(1)*C(7)
      RF(25) = RF(25)*C(2)*C(8)
      RB(25) = RB(25)*C(5)*C(6)
      RF(26) = RF(26)*C(3)*C(8)
      RB(26) = RB(26)*C(5)*C(7)
      RF(27) = RF(27)*C(5)*C(8)
      RB(27) = RB(27)*C(6)*C(7)
      RF(28) = RF(28)*C(5)*C(8)
      RB(28) = RB(28)*C(6)*C(7)
      CTB = CTOT+C(1)+5D0*C(4)+5D0*C(6)+C(10)
     * +5D-1*C(11)+2.5D0*C(12)+2D0*C(16)+2D0*C(14)
     * +2D0*C(15)
      RF(29) = RF(29)*CTB*C(3)*C(11)
      RB(29) = RB(29)*CTB*C(12)
      RF(30) = RF(30)*C(5)*C(11)
      RB(30) = RB(30)*C(2)*C(12)
      CTB = CTOT+C(1)+5D0*C(6)+C(10)+5D-1*C(11)
     * +C(12)+2D0*C(16)+2D0*C(14)+2D0*C(15)
      PR = RKLOW(2) * CTB / RF(31)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 6.8D-2*EXP(-T/1.97D2) + 9.32D-1*EXP(-T/1.54D3)
     *     + EXP(-1.03D4/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(31) = RF(31) * PCOR
      RB(31) = RB(31) * PCOR
      RF(31) = RF(31)*C(1)*C(11)
      RB(31) = RB(31)*C(13)
      RF(32) = RF(32)*C(4)*C(11)
      RB(32) = RB(32)*C(3)*C(12)
      RF(33) = RF(33)*C(7)*C(11)
      RB(33) = RB(33)*C(5)*C(12)
      RF(34) = RF(34)*C(3)
      RB(34) = RB(34)*C(2)*C(11)
      RF(35) = RF(35)*C(5)
      RB(35) = RB(35)*C(2)
      RF(36) = RF(36)*C(1)
      RB(36) = RB(36)*C(2)
      RF(37) = RF(37)*C(6)
      RB(37) = RB(37)*C(2)*C(13)
      RF(38) = RF(38)*C(4)
      RB(38) = RB(38)*C(3)
      CTB = CTOT+C(1)+5D0*C(6)+C(10)+5D-1*C(11)
     * +C(12)+2D0*C(16)+2D0*C(14)+2D0*C(15)
      PR = RKLOW(3) * CTB / RF(39)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 4.243D-1*EXP(-T/2.37D2) + 5.757D-1*EXP(-T/1.652D3)
     *     + EXP(-5.069D3/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(39) = RF(39) * PCOR
      RB(39) = RB(39) * PCOR
      RF(39) = RF(39)*C(11)
      RB(39) = RB(39)*C(17)
      RF(40) = RF(40)*C(12)
      RB(40) = RB(40)*C(11)
      CTB = CTOT+C(1)+5D0*C(6)+C(10)+5D-1*C(11)
     * +C(12)+2D0*C(16)+2D0*C(14)+2D0*C(15)
      PR = RKLOW(4) * CTB / RF(41)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 2.176D-1*EXP(-T/2.71D2) + 7.824D-1*EXP(-T/2.755D3)
     *     + EXP(-6.57D3/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(41) = RF(41) * PCOR
      RB(41) = RB(41) * PCOR
      RF(41) = RF(41)*C(2)
      RB(41) = RB(41)*C(13)
      RF(42) = RF(42)*C(2)
      RB(42) = RB(42)*C(1)*C(11)
      RF(43) = RF(43)*C(3)
      RB(43) = RB(43)*C(5)*C(11)
      RF(44) = RF(44)*C(3)
      RB(44) = RB(44)*C(2)*C(12)
      RF(45) = RF(45)*C(5)
      RB(45) = RB(45)*C(6)*C(11)
      CTB = CTOT+C(1)+5D0*C(6)+C(10)+5D-1*C(11)
     * +C(12)+2D0*C(16)+2D0*C(14)+2D0*C(15)
      RF(46) = RF(46)*CTB
      RB(46) = RB(46)*CTB*C(2)*C(11)
      RF(47) = RF(47)*C(4)
      RB(47) = RB(47)*C(7)*C(11)
      CTB = CTOT+C(1)+5D0*C(6)+C(10)+5D-1*C(11)
     * +C(12)+2D0*C(16)+2D0*C(14)+2D0*C(15)
      PR = RKLOW(5) * CTB / RF(48)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 3.2D-1*EXP(-T/7.8D1) + 6.8D-1*EXP(-T/1.995D3)
     *     + EXP(-5.59D3/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(48) = RF(48) * PCOR
      RB(48) = RB(48) * PCOR
      RF(48) = RF(48)*C(2)
      RB(48) = RB(48)*C(9)
      RF(49) = RF(49)*C(1)
      RB(49) = RB(49)*C(2)*C(9)
      RF(50) = RF(50)*C(3)
      RB(50) = RB(50)*C(2)
      RF(51) = RF(51)*C(4)
      RB(51) = RB(51)*C(5)
      RF(52) = RF(52)*C(4)
      RB(52) = RB(52)*C(2)*C(2)*C(12)
      RF(53) = RF(53)*C(5)
      RB(53) = RB(53)*C(2)*C(13)
      RF(54) = RF(54)*C(5)
      RB(54) = RB(54)*C(6)
      RF(55) = RF(55)*C(7)
      RB(55) = RB(55)*C(5)*C(13)
      CTB = CTOT+C(1)+5D0*C(6)+C(10)+5D-1*C(11)
     * +C(12)+2D0*C(16)+2D0*C(14)+2D0*C(15)
      PR = RKLOW(6) * CTB / RF(56)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 4.093D-1*EXP(-T/2.75D2) + 5.907D-1*EXP(-T/1.226D3)
     *     + EXP(-5.185D3/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(56) = RF(56) * PCOR
      RB(56) = RB(56) * PCOR
      RF(56) = RF(56)*C(11)
      RB(56) = RB(56)*C(18)
      RB(57) = RB(57)*C(2)*C(14)
      RB(58) = RB(58)*C(1)*C(14)
      RF(59) = RF(59)*C(22)
      RB(59) = RB(59)*C(22)
      RF(60) = RF(60)*C(2)
      RB(60) = RB(60)*C(1)
      RF(61) = RF(61)*C(3)
      RB(61) = RB(61)*C(1)*C(11)
      RF(62) = RF(62)*C(3)
      RB(62) = RB(62)*C(2)
      RF(63) = RF(63)*C(5)
      RB(63) = RB(63)*C(2)*C(13)
      RF(64) = RF(64)*C(1)
      RB(64) = RB(64)*C(2)*C(9)
      RF(65) = RF(65)*C(4)
      RB(65) = RB(65)*C(2)*C(5)*C(11)
      RF(66) = RF(66)*C(4)
      RB(66) = RB(66)*C(6)*C(11)
      RF(67) = RF(67)*C(6)
      RB(67) = RB(67)*C(6)
      RF(68) = RF(68)*C(11)
      RB(68) = RB(68)*C(11)
      RF(69) = RF(69)*C(12)
      RB(69) = RB(69)*C(12)
      RF(70) = RF(70)*C(12)
      RB(70) = RB(70)*C(11)*C(13)
      CTB = CTOT+C(1)+5D0*C(6)+C(10)+5D-1*C(11)
     * +C(12)+2D0*C(16)+2D0*C(14)+2D0*C(15)
      PR = RKLOW(7) * CTB / RF(71)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 2.42D-1*EXP(-T/9.4D1) + 7.58D-1*EXP(-T/1.555D3)
     *     + EXP(-4.2D3/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(71) = RF(71) * PCOR
      RB(71) = RB(71) * PCOR
      RF(71) = RF(71)*C(2)*C(13)
      RF(72) = RF(72)*C(2)*C(13)
      RB(72) = RB(72)*C(1)
      RF(73) = RF(73)*C(3)*C(13)
      RB(73) = RB(73)*C(5)
      RF(74) = RF(74)*C(5)*C(13)
      RB(74) = RB(74)*C(6)
      RF(75) = RF(75)*C(4)*C(13)
      RB(75) = RB(75)*C(7)
      RF(76) = RF(76)*C(7)*C(13)
      RB(76) = RB(76)*C(8)
      RF(77) = RF(77)*C(13)
      RB(77) = RB(77)*C(2)*C(18)
      CTB = CTOT+C(1)+5D0*C(6)+C(10)+5D-1*C(11)
     * +C(12)+2D0*C(16)+2D0*C(14)+2D0*C(15)
      PR = RKLOW(8) * CTB / RF(78)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 2.17D-1*EXP(-T/7.4D1) + 7.83D-1*EXP(-T/2.941D3)
     *     + EXP(-6.964D3/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(78) = RF(78) * PCOR
      RB(78) = RB(78) * PCOR
      RF(78) = RF(78)*C(2)*C(9)
      RB(78) = RB(78)*C(10)
      RF(79) = RF(79)*C(3)*C(9)
      RB(79) = RB(79)*C(2)*C(13)
      RF(80) = RF(80)*C(5)*C(9)
      RB(80) = RB(80)*C(6)
      RF(81) = RF(81)*C(5)*C(9)
      RB(81) = RB(81)*C(6)
      RF(82) = RF(82)*C(4)*C(9)
      RB(82) = RB(82)*C(3)
      RF(83) = RF(83)*C(4)*C(9)
      RB(83) = RB(83)*C(5)*C(13)
      RF(84) = RF(84)*C(7)*C(9)
      RB(84) = RB(84)*C(4)*C(10)
      RF(85) = RF(85)*C(7)*C(9)
      RB(85) = RB(85)*C(5)
      RF(86) = RF(86)*C(8)*C(9)
      RB(86) = RB(86)*C(7)*C(10)
      RF(87) = RF(87)*C(9)
      RB(87) = RB(87)*C(2)
      RF(88) = RF(88)*C(9)
      RB(88) = RB(88)*C(10)*C(11)
      CTB = CTOT+C(1)+5D0*C(6)+C(10)+5D-1*C(11)
     * +C(12)+2D0*C(16)+2D0*C(14)+2D0*C(15)
      PR = RKLOW(9) * CTB / RF(89)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 3.827D-1*EXP(-T/1.3076D1) + 6.173D-1*EXP(-T/2.078D3)
     *     + EXP(-5.093D3/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(89) = RF(89) * PCOR
      RB(89) = RB(89) * PCOR
      RF(89) = RF(89)*C(9)
      RB(89) = RB(89)*C(19)
      RF(90) = RF(90)*C(9)*C(13)
      RB(90) = RB(90)*C(10)
      RF(91) = RF(91)*C(9)
      RB(91) = RB(91)*C(2)*C(15)
      RF(92) = RF(92)*C(9)
      RB(92) = RB(92)*C(2)*C(15)
      CTB = CTOT+C(1)+5D0*C(6)+C(10)+5D-1*C(11)
     * +C(12)+2D0*C(16)+2D0*C(14)+2D0*C(15)
      PR = RKLOW(10) * CTB / RF(93)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 4.675D-1*EXP(-T/1.51D2) + 5.325D-1*EXP(-T/1.038D3)
     *     + EXP(-4.97D3/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(93) = RF(93) * PCOR
      RB(93) = RB(93) * PCOR
      RF(93) = RF(93)*C(9)*C(9)
      RB(93) = RB(93)*C(16)
      RF(94) = RF(94)*C(9)*C(9)
      RB(94) = RB(94)*C(2)
      RF(95) = RF(95)*C(9)*C(17)
      RB(95) = RB(95)*C(11)*C(15)
      RF(96) = RF(96)*C(2)
      RB(96) = RB(96)*C(1)*C(13)
      RF(97) = RF(97)*C(2)
      RB(97) = RB(97)*C(5)*C(9)
      RF(98) = RF(98)*C(2)
      RB(98) = RB(98)*C(6)
      RF(99) = RF(99)*C(3)
      RB(99) = RB(99)*C(5)*C(13)
      RF(100) = RF(100)*C(5)
      RB(100) = RB(100)*C(6)*C(13)
      RF(101) = RF(101)*C(4)
      RB(101) = RB(101)*C(7)*C(13)
      RF(102) = RF(102)*C(2)*C(10)
      RB(102) = RB(102)*C(1)*C(9)
      RF(103) = RF(103)*C(3)*C(10)
      RB(103) = RB(103)*C(5)*C(9)
      RF(104) = RF(104)*C(5)*C(10)
      RB(104) = RB(104)*C(6)*C(9)
      RF(105) = RF(105)*C(10)
      RB(105) = RB(105)*C(2)*C(15)
      RF(106) = RF(106)*C(10)
      RB(106) = RB(106)*C(9)*C(9)
      RF(107) = RF(107)*C(10)
      RB(107) = RB(107)*C(9)*C(9)
      RF(108) = RF(108)*C(2)*C(17)
      RB(108) = RB(108)*C(11)
      RF(109) = RF(109)*C(3)*C(17)
      RB(109) = RB(109)*C(2)*C(11)*C(11)
      RF(110) = RF(110)*C(4)*C(17)
      RB(110) = RB(110)*C(5)*C(11)*C(11)
      RF(111) = RF(111)*C(17)
      RB(111) = RB(111)*C(11)*C(14)
      RF(112) = RF(112)*C(17)
      RB(112) = RB(112)*C(11)
      RF(113) = RF(113)*C(17)*C(17)
      RB(113) = RB(113)*C(11)*C(11)*C(14)
      CTB = CTOT+C(1)+5D0*C(6)+C(10)+5D-1*C(11)
     * +C(12)+2D0*C(16)+1.5D0*C(14)+1.5D0*C(15)
      PR = RKLOW(11) * CTB / RF(114)
      PCOR = PR / (1.0 + PR)
      RF(114) = RF(114) * PCOR
      RB(114) = RB(114) * PCOR
      RF(114) = RF(114)*C(14)
      CTB = CTOT+C(1)+5D0*C(6)+C(10)+5D-1*C(11)
     * +C(12)+2D0*C(16)+2D0*C(14)+2D0*C(15)
      PR = RKLOW(12) * CTB / RF(115)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = -9.816D-1*EXP(-T/5.3837D3) + 1.9816D0*EXP(-T/4.2932D0)
     *     + EXP(7.95D-2/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(115) = RF(115) * PCOR
      RB(115) = RB(115) * PCOR
      RB(115) = RB(115)*C(2)*C(14)
      RF(116) = RF(116)*C(3)*C(14)
      RB(116) = RB(116)*C(2)*C(17)
      RF(117) = RF(117)*C(3)*C(14)
      RB(117) = RB(117)*C(11)
      RF(118) = RF(118)*C(5)*C(14)
      RB(118) = RB(118)*C(2)*C(18)
      RF(119) = RF(119)*C(5)*C(14)
      RB(119) = RB(119)*C(9)*C(11)
      RF(120) = RF(120)*C(14)
      RB(120) = RB(120)*C(11)
      CTB = CTOT
      RF(121) = RF(121)*CTB*C(9)*C(14)
      RB(121) = RB(121)*CTB*C(20)
      RF(122) = RF(122)*C(2)
      RB(122) = RB(122)*C(2)*C(14)
      RF(123) = RF(123)*C(3)
      RB(123) = RB(123)*C(11)
      RF(124) = RF(124)*C(5)
      RB(124) = RB(124)*C(2)*C(18)
      RF(125) = RF(125)*C(4)
      RB(125) = RB(125)*C(12)
      CTB = CTOT+C(1)+5D0*C(6)+C(10)+5D-1*C(11)
     * +C(12)+2D0*C(16)+2D0*C(14)+2D0*C(15)
      PR = RKLOW(13) * CTB / RF(126)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 6.63D-1*EXP(-T/1.707D3) + 3.37D-1*EXP(-T/3.2D3)
     *     + EXP(-4.131D3/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(126) = RF(126) * PCOR
      RB(126) = RB(126) * PCOR
      RF(126) = RF(126)*C(2)*C(18)
      RF(127) = RF(127)*C(2)*C(18)
      RB(127) = RB(127)*C(1)*C(17)
      RF(128) = RF(128)*C(2)*C(18)
      RB(128) = RB(128)*C(9)*C(11)
      RF(129) = RF(129)*C(3)*C(18)
      RB(129) = RB(129)*C(5)*C(17)
      RF(130) = RF(130)*C(3)*C(18)
      RB(130) = RB(130)*C(12)
      RF(131) = RF(131)*C(5)*C(18)
      RB(131) = RB(131)*C(6)*C(17)
      CTB = CTOT+C(1)+5D0*C(6)+C(10)+5D-1*C(11)
     * +C(12)+2D0*C(16)+2D0*C(14)+2D0*C(15)
      PR = RKLOW(14) * CTB / RF(132)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 2.18D-1*EXP(-T/2.075D2) + 7.82D-1*EXP(-T/2.663D3)
     *     + EXP(-6.095D3/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(132) = RF(132) * PCOR
      RB(132) = RB(132) * PCOR
      RF(132) = RF(132)*C(2)
      RB(132) = RB(132)*C(15)
      RF(133) = RF(133)*C(2)
      RB(133) = RB(133)*C(1)*C(14)
      RF(134) = RF(134)*C(2)
      RB(134) = RB(134)*C(1)
      RF(135) = RF(135)*C(3)
      RB(135) = RB(135)*C(2)*C(18)
      RF(136) = RF(136)*C(3)
      RB(136) = RB(136)*C(9)*C(11)
      RF(137) = RF(137)*C(5)
      RB(137) = RB(137)*C(6)*C(14)
      RF(138) = RF(138)*C(4)
      RB(138) = RB(138)*C(7)*C(14)
      RF(139) = RF(139)*C(4)
      RB(139) = RB(139)*C(3)
      RF(140) = RF(140)*C(4)
      RB(140) = RB(140)*C(13)
      RF(141) = RF(141)*C(7)
      RB(141) = RB(141)*C(5)
      RF(142) = RF(142)*C(8)
      RB(142) = RB(142)*C(7)*C(15)
      RB(143) = RB(143)*C(11)*C(15)
      RF(144) = RF(144)*C(9)
      RB(144) = RB(144)*C(10)*C(14)
      CTB = CTOT+C(1)+5D0*C(6)+C(10)+5D-1*C(11)
     * +C(12)+2D0*C(16)+2D0*C(14)+2D0*C(15)
      PR = RKLOW(15) * CTB / RF(145)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 8.25D-1*EXP(-T/1.3406D3) + 1.75D-1*EXP(-T/6D4)
     *     + EXP(-1.01398D4/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(145) = RF(145) * PCOR
      RB(145) = RB(145) * PCOR
      RF(145) = RF(145)*C(9)
      RB(145) = RB(145)*C(21)
      RF(146) = RF(146)*C(9)
      RB(146) = RB(146)*C(2)*C(20)
      RB(147) = RB(147)*C(9)*C(11)
      CTB = CTOT+C(1)+5D0*C(6)+C(10)+5D-1*C(11)
     * +C(12)+2D0*C(16)+2D0*C(14)+2D0*C(15)
      PR = RKLOW(16) * CTB / RF(148)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 4.5D-1*EXP(-T/8.9D3) + 5.5D-1*EXP(-T/4.35D3)
     *     + EXP(-7.244D3/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(148) = RF(148) * PCOR
      RB(148) = RB(148) * PCOR
      RF(148) = RF(148)*C(2)
      RB(148) = RB(148)*C(19)
      RF(149) = RF(149)*C(2)
      RB(149) = RB(149)*C(9)
      RF(150) = RF(150)*C(2)
      RB(150) = RB(150)*C(1)*C(18)
      RF(151) = RF(151)*C(3)
      RB(151) = RB(151)*C(5)*C(18)
      RF(152) = RF(152)*C(5)
      RB(152) = RB(152)*C(6)*C(18)
      RF(153) = RF(153)*C(4)
      RB(153) = RB(153)*C(7)*C(18)
      RF(154) = RF(154)*C(4)
      RB(154) = RB(154)*C(5)*C(11)*C(13)
      CTB = CTOT+C(1)+5D0*C(6)+C(10)+5D-1*C(11)
     * +C(12)+2D0*C(16)+2D0*C(14)+2D0*C(15)
      PR = RKLOW(17) * CTB / RF(155)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 2.655D-1*EXP(-T/1.8D2) + 7.345D-1*EXP(-T/1.035D3)
     *     + EXP(-5.417D3/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(155) = RF(155) * PCOR
      RB(155) = RB(155) * PCOR
      RF(155) = RF(155)*C(15)
      RB(155) = RB(155)*C(1)
      CTB = CTOT+C(1)+5D0*C(6)+C(10)+5D-1*C(11)
     * +C(12)+2D0*C(16)+2D0*C(14)+2D0*C(15)
      PR = RKLOW(18) * CTB / RF(156)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 2.47D-2*EXP(-T/2.1D2) + 9.753D-1*EXP(-T/9.84D2)
     *     + EXP(-4.374D3/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(156) = RF(156) * PCOR
      RB(156) = RB(156) * PCOR
      RF(156) = RF(156)*C(2)*C(15)
      RF(157) = RF(157)*C(2)*C(15)
      RB(157) = RB(157)*C(1)
      RF(158) = RF(158)*C(3)*C(15)
      RB(158) = RB(158)*C(5)
      RF(159) = RF(159)*C(3)*C(15)
      RB(159) = RB(159)*C(9)
      RF(160) = RF(160)*C(3)*C(15)
      RB(160) = RB(160)*C(13)
      RF(161) = RF(161)*C(5)*C(15)
      RB(161) = RB(161)*C(6)
      RF(162) = RF(162)*C(4)*C(15)
      RB(162) = RB(162)*C(7)
      RF(163) = RF(163)*C(7)*C(15)
      RB(163) = RB(163)*C(5)*C(19)
      RF(164) = RF(164)*C(15)
      RB(164) = RB(164)*C(11)
      RF(165) = RF(165)*C(15)
      RB(165) = RB(165)*C(2)*C(20)
      RF(166) = RF(166)*C(15)
      RB(166) = RB(166)*C(10)
      RF(167) = RF(167)*C(15)
      RB(167) = RB(167)*C(2)*C(20)
      RF(168) = RF(168)*C(9)*C(15)
      RB(168) = RB(168)*C(10)
      RF(169) = RF(169)*C(9)*C(15)
      CTB = CTOT+C(1)+5D0*C(6)+C(10)+5D-1*C(11)
     * +C(12)+2D0*C(16)+2D0*C(14)+2D0*C(15)
      PR = RKLOW(19) * CTB / RF(170)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 1.578D-1*EXP(-T/1.25D2) + 8.422D-1*EXP(-T/2.219D3)
     *     + EXP(-6.882D3/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(170) = RF(170) * PCOR
      RB(170) = RB(170) * PCOR
      RF(170) = RF(170)*C(2)
      RB(170) = RB(170)*C(16)
      RF(171) = RF(171)*C(2)
      RB(171) = RB(171)*C(1)*C(15)
      RF(172) = RF(172)*C(3)
      RB(172) = RB(172)*C(9)*C(13)
      RF(173) = RF(173)*C(3)
      RB(173) = RB(173)*C(2)*C(19)
      RF(174) = RF(174)*C(4)
      RB(174) = RB(174)*C(7)*C(15)
      RF(175) = RF(175)*C(7)
      RB(175) = RB(175)*C(4)*C(16)
      RF(176) = RF(176)*C(7)
      RB(176) = RB(176)*C(8)*C(15)
      RF(177) = RF(177)*C(7)
      RB(177) = RB(177)*C(5)*C(9)*C(13)
      RF(178) = RF(178)*C(8)
      RB(178) = RB(178)*C(7)*C(16)
      RB(179) = RB(179)*C(11)*C(16)
      RF(180) = RF(180)*C(2)*C(16)
      RB(180) = RB(180)*C(1)
      RF(181) = RF(181)*C(3)*C(16)
      RB(181) = RB(181)*C(5)
      RF(182) = RF(182)*C(5)*C(16)
      RB(182) = RB(182)*C(6)
      RF(183) = RF(183)*C(16)
      RB(183) = RB(183)*C(9)
      RF(184) = RF(184)*C(9)*C(16)
      RB(184) = RB(184)*C(10)
      CTB = CTOT+C(1)+5D0*C(6)+C(10)+5D-1*C(11)
     * +C(12)+2D0*C(16)+2D0*C(14)+2D0*C(15)
      PR = RKLOW(20) * CTB / RF(185)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 9.8D-1*EXP(-T/1.0966D3) + 2D-2*EXP(-T/1.0966D3)
     *     + EXP(-6.8595D3/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(185) = RF(185) * PCOR
      RB(185) = RB(185) * PCOR
      RF(185) = RF(185)*C(2)*C(20)
      RB(185) = RB(185)*C(21)
      RF(186) = RF(186)*C(2)*C(20)
      RB(186) = RB(186)*C(10)
      RF(187) = RF(187)*C(7)*C(20)
      RB(187) = RB(187)*C(4)*C(21)
      RF(188) = RF(188)*C(7)*C(20)
      RB(188) = RB(188)*C(5)*C(13)
      RF(189) = RF(189)*C(20)
      RB(189) = RB(189)*C(11)*C(21)
      CTB = CTOT+C(1)+5D0*C(6)+C(10)+5D-1*C(11)
     * +C(12)+2D0*C(16)
      PR = RKLOW(21) * CTB / RF(190)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 0D0*EXP(-T/1D3) + 1D0*EXP(-T/1.31D3)
     *     + EXP(-4.8097D4/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(190) = RF(190) * PCOR
      RB(190) = RB(190) * PCOR
      RF(190) = RF(190)*C(2)*C(21)
      RF(191) = RF(191)*C(2)*C(21)
      RB(191) = RB(191)*C(9)*C(15)
      RF(192) = RF(192)*C(2)*C(21)
      RB(192) = RB(192)*C(1)*C(20)
      RF(193) = RF(193)*C(3)*C(21)
      RB(193) = RB(193)*C(2)*C(9)*C(18)
      RF(194) = RF(194)*C(3)*C(21)
      RF(195) = RF(195)*C(3)*C(21)
      RB(195) = RB(195)*C(5)*C(20)
      RF(196) = RF(196)*C(5)*C(21)
      RB(196) = RB(196)*C(6)*C(20)
      RF(197) = RF(197)*C(7)*C(21)
      RB(197) = RB(197)*C(8)*C(20)
      RF(198) = RF(198)*C(9)*C(21)
      RB(198) = RB(198)*C(10)*C(20)
      RF(199) = RF(199)*C(2)
      RB(199) = RB(199)*C(9)
      RF(200) = RF(200)*C(2)
      RB(200) = RB(200)*C(1)*C(21)
      RF(201) = RF(201)*C(3)
      RB(201) = RB(201)*C(13)
      RF(202) = RF(202)*C(5)
      RB(202) = RB(202)*C(6)*C(21)
      RF(203) = RF(203)*C(4)
      RB(203) = RB(203)*C(7)*C(21)
      RF(204) = RF(204)*C(7)
      RB(204) = RB(204)*C(5)*C(13)
      RF(205) = RF(205)*C(9)
      RB(205) = RB(205)*C(10)*C(21)
      RB(206) = RB(206)*C(9)*C(20)
C
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE RDOT(RF, RB, WDOT)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION RF(*), RB(*), WDOT(*)
C
      DO K = 1, 22
         WDOT(K) = 0D0
      ENDDO
C
      ROP = RF(1)-RB(1)
      WDOT(2) = WDOT(2) -ROP
      WDOT(3) = WDOT(3) +ROP
      WDOT(4) = WDOT(4) -ROP
      WDOT(5) = WDOT(5) +ROP
      ROP = RF(2)-RB(2)
      WDOT(1) = WDOT(1) -ROP
      WDOT(2) = WDOT(2) +ROP
      WDOT(3) = WDOT(3) -ROP
      WDOT(5) = WDOT(5) +ROP
      ROP = RF(3)-RB(3)
      WDOT(1) = WDOT(1) -ROP
      WDOT(2) = WDOT(2) +ROP
      WDOT(5) = WDOT(5) -ROP
      WDOT(6) = WDOT(6) +ROP
      ROP = RF(4)-RB(4)
      WDOT(3) = WDOT(3) +ROP
      WDOT(5) = WDOT(5) -ROP -ROP
      WDOT(6) = WDOT(6) +ROP
      ROP = RF(5)-RB(5)
      WDOT(1) = WDOT(1) +ROP
      WDOT(2) = WDOT(2) -ROP -ROP
      ROP = RF(6)-RB(6)
      WDOT(1) = WDOT(1) +ROP
      WDOT(2) = WDOT(2) -ROP -ROP
      ROP = RF(7)-RB(7)
      WDOT(1) = WDOT(1) +ROP
      WDOT(2) = WDOT(2) -ROP -ROP
      ROP = RF(8)-RB(8)
      WDOT(1) = WDOT(1) +ROP
      WDOT(2) = WDOT(2) -ROP -ROP
      ROP = RF(9)-RB(9)
      WDOT(2) = WDOT(2) -ROP
      WDOT(5) = WDOT(5) -ROP
      WDOT(6) = WDOT(6) +ROP
      ROP = RF(10)-RB(10)
      WDOT(2) = WDOT(2) -ROP
      WDOT(3) = WDOT(3) -ROP
      WDOT(5) = WDOT(5) +ROP
      ROP = RF(11)-RB(11)
      WDOT(3) = WDOT(3) -ROP -ROP
      WDOT(4) = WDOT(4) +ROP
      ROP = RF(12)-RB(12)
      WDOT(2) = WDOT(2) -ROP
      WDOT(4) = WDOT(4) -ROP
      WDOT(7) = WDOT(7) +ROP
      ROP = RF(13)-RB(13)
      WDOT(2) = WDOT(2) -ROP
      WDOT(4) = WDOT(4) -ROP
      WDOT(7) = WDOT(7) +ROP
      ROP = RF(14)-RB(14)
      WDOT(2) = WDOT(2) -ROP
      WDOT(4) = WDOT(4) -ROP
      WDOT(7) = WDOT(7) +ROP
      ROP = RF(15)-RB(15)
      WDOT(2) = WDOT(2) -ROP
      WDOT(4) = WDOT(4) -ROP
      WDOT(7) = WDOT(7) +ROP
      ROP = RF(16)-RB(16)
      WDOT(5) = WDOT(5) -ROP -ROP
      WDOT(8) = WDOT(8) +ROP
      ROP = RF(17)-RB(17)
      WDOT(2) = WDOT(2) -ROP
      WDOT(3) = WDOT(3) +ROP
      WDOT(6) = WDOT(6) +ROP
      WDOT(7) = WDOT(7) -ROP
      ROP = RF(18)-RB(18)
      WDOT(1) = WDOT(1) +ROP
      WDOT(2) = WDOT(2) -ROP
      WDOT(4) = WDOT(4) +ROP
      WDOT(7) = WDOT(7) -ROP
      ROP = RF(19)-RB(19)
      WDOT(2) = WDOT(2) -ROP
      WDOT(5) = WDOT(5) +ROP +ROP
      WDOT(7) = WDOT(7) -ROP
      ROP = RF(20)-RB(20)
      WDOT(3) = WDOT(3) -ROP
      WDOT(4) = WDOT(4) +ROP
      WDOT(5) = WDOT(5) +ROP
      WDOT(7) = WDOT(7) -ROP
      ROP = RF(21)-RB(21)
      WDOT(4) = WDOT(4) +ROP
      WDOT(5) = WDOT(5) -ROP
      WDOT(6) = WDOT(6) +ROP
      WDOT(7) = WDOT(7) -ROP
      ROP = RF(22)-RB(22)
      WDOT(4) = WDOT(4) +ROP
      WDOT(7) = WDOT(7) -ROP -ROP
      WDOT(8) = WDOT(8) +ROP
      ROP = RF(23)-RB(23)
      WDOT(4) = WDOT(4) +ROP
      WDOT(7) = WDOT(7) -ROP -ROP
      WDOT(8) = WDOT(8) +ROP
      ROP = RF(24)-RB(24)
      WDOT(1) = WDOT(1) +ROP
      WDOT(2) = WDOT(2) -ROP
      WDOT(7) = WDOT(7) +ROP
      WDOT(8) = WDOT(8) -ROP
      ROP = RF(25)-RB(25)
      WDOT(2) = WDOT(2) -ROP
      WDOT(5) = WDOT(5) +ROP
      WDOT(6) = WDOT(6) +ROP
      WDOT(8) = WDOT(8) -ROP
      ROP = RF(26)-RB(26)
      WDOT(3) = WDOT(3) -ROP
      WDOT(5) = WDOT(5) +ROP
      WDOT(7) = WDOT(7) +ROP
      WDOT(8) = WDOT(8) -ROP
      ROP = RF(27)-RB(27)
      WDOT(5) = WDOT(5) -ROP
      WDOT(6) = WDOT(6) +ROP
      WDOT(7) = WDOT(7) +ROP
      WDOT(8) = WDOT(8) -ROP
      ROP = RF(28)-RB(28)
      WDOT(5) = WDOT(5) -ROP
      WDOT(6) = WDOT(6) +ROP
      WDOT(7) = WDOT(7) +ROP
      WDOT(8) = WDOT(8) -ROP
      ROP = RF(29)-RB(29)
      WDOT(3) = WDOT(3) -ROP
      WDOT(11) = WDOT(11) -ROP
      WDOT(12) = WDOT(12) +ROP
      ROP = RF(30)-RB(30)
      WDOT(2) = WDOT(2) +ROP
      WDOT(5) = WDOT(5) -ROP
      WDOT(11) = WDOT(11) -ROP
      WDOT(12) = WDOT(12) +ROP
      ROP = RF(31)-RB(31)
      WDOT(1) = WDOT(1) -ROP
      WDOT(11) = WDOT(11) -ROP
      WDOT(13) = WDOT(13) +ROP
      ROP = RF(32)-RB(32)
      WDOT(3) = WDOT(3) +ROP
      WDOT(4) = WDOT(4) -ROP
      WDOT(11) = WDOT(11) -ROP
      WDOT(12) = WDOT(12) +ROP
      ROP = RF(33)-RB(33)
      WDOT(5) = WDOT(5) +ROP
      WDOT(7) = WDOT(7) -ROP
      WDOT(11) = WDOT(11) -ROP
      WDOT(12) = WDOT(12) +ROP
      ROP = RF(34)-RB(34)
      WDOT(2) = WDOT(2) +ROP
      WDOT(3) = WDOT(3) -ROP
      WDOT(11) = WDOT(11) +ROP
      ROP = RF(35)-RB(35)
      WDOT(2) = WDOT(2) +ROP
      WDOT(5) = WDOT(5) -ROP
      ROP = RF(36)-RB(36)
      WDOT(1) = WDOT(1) -ROP
      WDOT(2) = WDOT(2) +ROP
      ROP = RF(37)-RB(37)
      WDOT(2) = WDOT(2) +ROP
      WDOT(6) = WDOT(6) -ROP
      WDOT(13) = WDOT(13) +ROP
      ROP = RF(38)-RB(38)
      WDOT(3) = WDOT(3) +ROP
      WDOT(4) = WDOT(4) -ROP
      ROP = RF(39)-RB(39)
      WDOT(11) = WDOT(11) -ROP
      WDOT(17) = WDOT(17) +ROP
      ROP = RF(40)-RB(40)
      WDOT(11) = WDOT(11) +ROP
      WDOT(12) = WDOT(12) -ROP
      ROP = RF(41)-RB(41)
      WDOT(2) = WDOT(2) -ROP
      WDOT(13) = WDOT(13) +ROP
      ROP = RF(42)-RB(42)
      WDOT(1) = WDOT(1) +ROP
      WDOT(2) = WDOT(2) -ROP
      WDOT(11) = WDOT(11) +ROP
      ROP = RF(43)-RB(43)
      WDOT(3) = WDOT(3) -ROP
      WDOT(5) = WDOT(5) +ROP
      WDOT(11) = WDOT(11) +ROP
      ROP = RF(44)-RB(44)
      WDOT(2) = WDOT(2) +ROP
      WDOT(3) = WDOT(3) -ROP
      WDOT(12) = WDOT(12) +ROP
      ROP = RF(45)-RB(45)
      WDOT(5) = WDOT(5) -ROP
      WDOT(6) = WDOT(6) +ROP
      WDOT(11) = WDOT(11) +ROP
      ROP = RF(46)-RB(46)
      WDOT(2) = WDOT(2) +ROP
      WDOT(11) = WDOT(11) +ROP
      ROP = RF(47)-RB(47)
      WDOT(4) = WDOT(4) -ROP
      WDOT(7) = WDOT(7) +ROP
      WDOT(11) = WDOT(11) +ROP
      ROP = RF(48)-RB(48)
      WDOT(2) = WDOT(2) -ROP
      WDOT(9) = WDOT(9) +ROP
      ROP = RF(49)-RB(49)
      WDOT(1) = WDOT(1) -ROP
      WDOT(2) = WDOT(2) +ROP
      WDOT(9) = WDOT(9) +ROP
      ROP = RF(50)-RB(50)
      WDOT(2) = WDOT(2) +ROP
      WDOT(3) = WDOT(3) -ROP
      ROP = RF(51)-RB(51)
      WDOT(4) = WDOT(4) -ROP
      WDOT(5) = WDOT(5) +ROP
      ROP = RF(52)-RB(52)
      WDOT(2) = WDOT(2) +ROP +ROP
      WDOT(4) = WDOT(4) -ROP
      WDOT(12) = WDOT(12) +ROP
      ROP = RF(53)-RB(53)
      WDOT(2) = WDOT(2) +ROP
      WDOT(5) = WDOT(5) -ROP
      WDOT(13) = WDOT(13) +ROP
      ROP = RF(54)-RB(54)
      WDOT(5) = WDOT(5) -ROP
      WDOT(6) = WDOT(6) +ROP
      ROP = RF(55)-RB(55)
      WDOT(5) = WDOT(5) +ROP
      WDOT(7) = WDOT(7) -ROP
      WDOT(13) = WDOT(13) +ROP
      ROP = RF(56)-RB(56)
      WDOT(11) = WDOT(11) -ROP
      WDOT(18) = WDOT(18) +ROP
      ROP = RF(57)-RB(57)
      WDOT(2) = WDOT(2) +ROP
      WDOT(14) = WDOT(14) +ROP
      ROP = RF(58)-RB(58)
      WDOT(1) = WDOT(1) +ROP
      WDOT(14) = WDOT(14) +ROP
      ROP = RF(59)-RB(59)
      ROP = RF(60)-RB(60)
      WDOT(1) = WDOT(1) +ROP
      WDOT(2) = WDOT(2) -ROP
      ROP = RF(61)-RB(61)
      WDOT(1) = WDOT(1) +ROP
      WDOT(3) = WDOT(3) -ROP
      WDOT(11) = WDOT(11) +ROP
      ROP = RF(62)-RB(62)
      WDOT(2) = WDOT(2) +ROP
      WDOT(3) = WDOT(3) -ROP
      ROP = RF(63)-RB(63)
      WDOT(2) = WDOT(2) +ROP
      WDOT(5) = WDOT(5) -ROP
      WDOT(13) = WDOT(13) +ROP
      ROP = RF(64)-RB(64)
      WDOT(1) = WDOT(1) -ROP
      WDOT(2) = WDOT(2) +ROP
      WDOT(9) = WDOT(9) +ROP
      ROP = RF(65)-RB(65)
      WDOT(2) = WDOT(2) +ROP
      WDOT(4) = WDOT(4) -ROP
      WDOT(5) = WDOT(5) +ROP
      WDOT(11) = WDOT(11) +ROP
      ROP = RF(66)-RB(66)
      WDOT(4) = WDOT(4) -ROP
      WDOT(6) = WDOT(6) +ROP
      WDOT(11) = WDOT(11) +ROP
      ROP = RF(67)-RB(67)
      ROP = RF(68)-RB(68)
      ROP = RF(69)-RB(69)
      ROP = RF(70)-RB(70)
      WDOT(11) = WDOT(11) +ROP
      WDOT(12) = WDOT(12) -ROP
      WDOT(13) = WDOT(13) +ROP
      ROP = RF(71)-RB(71)
      WDOT(2) = WDOT(2) -ROP
      WDOT(13) = WDOT(13) -ROP
      ROP = RF(72)-RB(72)
      WDOT(1) = WDOT(1) +ROP
      WDOT(2) = WDOT(2) -ROP
      WDOT(13) = WDOT(13) -ROP
      ROP = RF(73)-RB(73)
      WDOT(3) = WDOT(3) -ROP
      WDOT(5) = WDOT(5) +ROP
      WDOT(13) = WDOT(13) -ROP
      ROP = RF(74)-RB(74)
      WDOT(5) = WDOT(5) -ROP
      WDOT(6) = WDOT(6) +ROP
      WDOT(13) = WDOT(13) -ROP
      ROP = RF(75)-RB(75)
      WDOT(4) = WDOT(4) -ROP
      WDOT(7) = WDOT(7) +ROP
      WDOT(13) = WDOT(13) -ROP
      ROP = RF(76)-RB(76)
      WDOT(7) = WDOT(7) -ROP
      WDOT(8) = WDOT(8) +ROP
      WDOT(13) = WDOT(13) -ROP
      ROP = RF(77)-RB(77)
      WDOT(2) = WDOT(2) +ROP
      WDOT(13) = WDOT(13) -ROP
      WDOT(18) = WDOT(18) +ROP
      ROP = RF(78)-RB(78)
      WDOT(2) = WDOT(2) -ROP
      WDOT(9) = WDOT(9) -ROP
      WDOT(10) = WDOT(10) +ROP
      ROP = RF(79)-RB(79)
      WDOT(2) = WDOT(2) +ROP
      WDOT(3) = WDOT(3) -ROP
      WDOT(9) = WDOT(9) -ROP
      WDOT(13) = WDOT(13) +ROP
      ROP = RF(80)-RB(80)
      WDOT(5) = WDOT(5) -ROP
      WDOT(6) = WDOT(6) +ROP
      WDOT(9) = WDOT(9) -ROP
      ROP = RF(81)-RB(81)
      WDOT(5) = WDOT(5) -ROP
      WDOT(6) = WDOT(6) +ROP
      WDOT(9) = WDOT(9) -ROP
      ROP = RF(82)-RB(82)
      WDOT(3) = WDOT(3) +ROP
      WDOT(4) = WDOT(4) -ROP
      WDOT(9) = WDOT(9) -ROP
      ROP = RF(83)-RB(83)
      WDOT(4) = WDOT(4) -ROP
      WDOT(5) = WDOT(5) +ROP
      WDOT(9) = WDOT(9) -ROP
      WDOT(13) = WDOT(13) +ROP
      ROP = RF(84)-RB(84)
      WDOT(4) = WDOT(4) +ROP
      WDOT(7) = WDOT(7) -ROP
      WDOT(9) = WDOT(9) -ROP
      WDOT(10) = WDOT(10) +ROP
      ROP = RF(85)-RB(85)
      WDOT(5) = WDOT(5) +ROP
      WDOT(7) = WDOT(7) -ROP
      WDOT(9) = WDOT(9) -ROP
      ROP = RF(86)-RB(86)
      WDOT(7) = WDOT(7) +ROP
      WDOT(8) = WDOT(8) -ROP
      WDOT(9) = WDOT(9) -ROP
      WDOT(10) = WDOT(10) +ROP
      ROP = RF(87)-RB(87)
      WDOT(2) = WDOT(2) +ROP
      WDOT(9) = WDOT(9) -ROP
      ROP = RF(88)-RB(88)
      WDOT(9) = WDOT(9) -ROP
      WDOT(10) = WDOT(10) +ROP
      WDOT(11) = WDOT(11) +ROP
      ROP = RF(89)-RB(89)
      WDOT(9) = WDOT(9) -ROP
      WDOT(19) = WDOT(19) +ROP
      ROP = RF(90)-RB(90)
      WDOT(9) = WDOT(9) -ROP
      WDOT(10) = WDOT(10) +ROP
      WDOT(13) = WDOT(13) -ROP
      ROP = RF(91)-RB(91)
      WDOT(2) = WDOT(2) +ROP
      WDOT(9) = WDOT(9) -ROP
      WDOT(15) = WDOT(15) +ROP
      ROP = RF(92)-RB(92)
      WDOT(2) = WDOT(2) +ROP
      WDOT(9) = WDOT(9) -ROP
      WDOT(15) = WDOT(15) +ROP
      ROP = RF(93)-RB(93)
      WDOT(9) = WDOT(9) -ROP -ROP
      WDOT(16) = WDOT(16) +ROP
      ROP = RF(94)-RB(94)
      WDOT(2) = WDOT(2) +ROP
      WDOT(9) = WDOT(9) -ROP -ROP
      ROP = RF(95)-RB(95)
      WDOT(9) = WDOT(9) -ROP
      WDOT(11) = WDOT(11) +ROP
      WDOT(15) = WDOT(15) +ROP
      WDOT(17) = WDOT(17) -ROP
      ROP = RF(96)-RB(96)
      WDOT(1) = WDOT(1) +ROP
      WDOT(2) = WDOT(2) -ROP
      WDOT(13) = WDOT(13) +ROP
      ROP = RF(97)-RB(97)
      WDOT(2) = WDOT(2) -ROP
      WDOT(5) = WDOT(5) +ROP
      WDOT(9) = WDOT(9) +ROP
      ROP = RF(98)-RB(98)
      WDOT(2) = WDOT(2) -ROP
      WDOT(6) = WDOT(6) +ROP
      ROP = RF(99)-RB(99)
      WDOT(3) = WDOT(3) -ROP
      WDOT(5) = WDOT(5) +ROP
      WDOT(13) = WDOT(13) +ROP
      ROP = RF(100)-RB(100)
      WDOT(5) = WDOT(5) -ROP
      WDOT(6) = WDOT(6) +ROP
      WDOT(13) = WDOT(13) +ROP
      ROP = RF(101)-RB(101)
      WDOT(4) = WDOT(4) -ROP
      WDOT(7) = WDOT(7) +ROP
      WDOT(13) = WDOT(13) +ROP
      ROP = RF(102)-RB(102)
      WDOT(1) = WDOT(1) +ROP
      WDOT(2) = WDOT(2) -ROP
      WDOT(9) = WDOT(9) +ROP
      WDOT(10) = WDOT(10) -ROP
      ROP = RF(103)-RB(103)
      WDOT(3) = WDOT(3) -ROP
      WDOT(5) = WDOT(5) +ROP
      WDOT(9) = WDOT(9) +ROP
      WDOT(10) = WDOT(10) -ROP
      ROP = RF(104)-RB(104)
      WDOT(5) = WDOT(5) -ROP
      WDOT(6) = WDOT(6) +ROP
      WDOT(9) = WDOT(9) +ROP
      WDOT(10) = WDOT(10) -ROP
      ROP = RF(105)-RB(105)
      WDOT(2) = WDOT(2) +ROP
      WDOT(10) = WDOT(10) -ROP
      WDOT(15) = WDOT(15) +ROP
      ROP = RF(106)-RB(106)
      WDOT(9) = WDOT(9) +ROP +ROP
      WDOT(10) = WDOT(10) -ROP
      ROP = RF(107)-RB(107)
      WDOT(9) = WDOT(9) +ROP +ROP
      WDOT(10) = WDOT(10) -ROP
      ROP = RF(108)-RB(108)
      WDOT(2) = WDOT(2) -ROP
      WDOT(11) = WDOT(11) +ROP
      WDOT(17) = WDOT(17) -ROP
      ROP = RF(109)-RB(109)
      WDOT(2) = WDOT(2) +ROP
      WDOT(3) = WDOT(3) -ROP
      WDOT(11) = WDOT(11) +ROP +ROP
      WDOT(17) = WDOT(17) -ROP
      ROP = RF(110)-RB(110)
      WDOT(4) = WDOT(4) -ROP
      WDOT(5) = WDOT(5) +ROP
      WDOT(11) = WDOT(11) +ROP +ROP
      WDOT(17) = WDOT(17) -ROP
      ROP = RF(111)-RB(111)
      WDOT(11) = WDOT(11) +ROP
      WDOT(14) = WDOT(14) +ROP
      WDOT(17) = WDOT(17) -ROP
      ROP = RF(112)-RB(112)
      WDOT(11) = WDOT(11) +ROP
      WDOT(17) = WDOT(17) -ROP
      ROP = RF(113)-RB(113)
      WDOT(11) = WDOT(11) +ROP +ROP
      WDOT(14) = WDOT(14) +ROP
      WDOT(17) = WDOT(17) -ROP -ROP
      ROP = RF(114)-RB(114)
      WDOT(14) = WDOT(14) -ROP
      ROP = RF(115)-RB(115)
      WDOT(2) = WDOT(2) +ROP
      WDOT(14) = WDOT(14) +ROP
      ROP = RF(116)-RB(116)
      WDOT(2) = WDOT(2) +ROP
      WDOT(3) = WDOT(3) -ROP
      WDOT(14) = WDOT(14) -ROP
      WDOT(17) = WDOT(17) +ROP
      ROP = RF(117)-RB(117)
      WDOT(3) = WDOT(3) -ROP
      WDOT(11) = WDOT(11) +ROP
      WDOT(14) = WDOT(14) -ROP
      ROP = RF(118)-RB(118)
      WDOT(2) = WDOT(2) +ROP
      WDOT(5) = WDOT(5) -ROP
      WDOT(14) = WDOT(14) -ROP
      WDOT(18) = WDOT(18) +ROP
      ROP = RF(119)-RB(119)
      WDOT(5) = WDOT(5) -ROP
      WDOT(9) = WDOT(9) +ROP
      WDOT(11) = WDOT(11) +ROP
      WDOT(14) = WDOT(14) -ROP
      ROP = RF(120)-RB(120)
      WDOT(11) = WDOT(11) +ROP
      WDOT(14) = WDOT(14) -ROP
      ROP = RF(121)-RB(121)
      WDOT(9) = WDOT(9) -ROP
      WDOT(14) = WDOT(14) -ROP
      WDOT(20) = WDOT(20) +ROP
      ROP = RF(122)-RB(122)
      WDOT(14) = WDOT(14) +ROP
      ROP = RF(123)-RB(123)
      WDOT(3) = WDOT(3) -ROP
      WDOT(11) = WDOT(11) +ROP
      ROP = RF(124)-RB(124)
      WDOT(2) = WDOT(2) +ROP
      WDOT(5) = WDOT(5) -ROP
      WDOT(18) = WDOT(18) +ROP
      ROP = RF(125)-RB(125)
      WDOT(4) = WDOT(4) -ROP
      WDOT(12) = WDOT(12) +ROP
      ROP = RF(126)-RB(126)
      WDOT(2) = WDOT(2) -ROP
      WDOT(18) = WDOT(18) -ROP
      ROP = RF(127)-RB(127)
      WDOT(1) = WDOT(1) +ROP
      WDOT(2) = WDOT(2) -ROP
      WDOT(17) = WDOT(17) +ROP
      WDOT(18) = WDOT(18) -ROP
      ROP = RF(128)-RB(128)
      WDOT(2) = WDOT(2) -ROP
      WDOT(9) = WDOT(9) +ROP
      WDOT(11) = WDOT(11) +ROP
      WDOT(18) = WDOT(18) -ROP
      ROP = RF(129)-RB(129)
      WDOT(3) = WDOT(3) -ROP
      WDOT(5) = WDOT(5) +ROP
      WDOT(17) = WDOT(17) +ROP
      WDOT(18) = WDOT(18) -ROP
      ROP = RF(130)-RB(130)
      WDOT(3) = WDOT(3) -ROP
      WDOT(12) = WDOT(12) +ROP
      WDOT(18) = WDOT(18) -ROP
      ROP = RF(131)-RB(131)
      WDOT(5) = WDOT(5) -ROP
      WDOT(6) = WDOT(6) +ROP
      WDOT(17) = WDOT(17) +ROP
      WDOT(18) = WDOT(18) -ROP
      ROP = RF(132)-RB(132)
      WDOT(2) = WDOT(2) -ROP
      WDOT(15) = WDOT(15) +ROP
      ROP = RF(133)-RB(133)
      WDOT(1) = WDOT(1) +ROP
      WDOT(2) = WDOT(2) -ROP
      WDOT(14) = WDOT(14) +ROP
      ROP = RF(134)-RB(134)
      WDOT(1) = WDOT(1) +ROP
      WDOT(2) = WDOT(2) -ROP
      ROP = RF(135)-RB(135)
      WDOT(2) = WDOT(2) +ROP
      WDOT(3) = WDOT(3) -ROP
      WDOT(18) = WDOT(18) +ROP
      ROP = RF(136)-RB(136)
      WDOT(3) = WDOT(3) -ROP
      WDOT(9) = WDOT(9) +ROP
      WDOT(11) = WDOT(11) +ROP
      ROP = RF(137)-RB(137)
      WDOT(5) = WDOT(5) -ROP
      WDOT(6) = WDOT(6) +ROP
      WDOT(14) = WDOT(14) +ROP
      ROP = RF(138)-RB(138)
      WDOT(4) = WDOT(4) -ROP
      WDOT(7) = WDOT(7) +ROP
      WDOT(14) = WDOT(14) +ROP
      ROP = RF(139)-RB(139)
      WDOT(3) = WDOT(3) +ROP
      WDOT(4) = WDOT(4) -ROP
      ROP = RF(140)-RB(140)
      WDOT(4) = WDOT(4) -ROP
      WDOT(13) = WDOT(13) +ROP
      ROP = RF(141)-RB(141)
      WDOT(5) = WDOT(5) +ROP
      WDOT(7) = WDOT(7) -ROP
      ROP = RF(142)-RB(142)
      WDOT(7) = WDOT(7) +ROP
      WDOT(8) = WDOT(8) -ROP
      WDOT(15) = WDOT(15) +ROP
      ROP = RF(143)-RB(143)
      WDOT(11) = WDOT(11) +ROP
      WDOT(15) = WDOT(15) +ROP
      ROP = RF(144)-RB(144)
      WDOT(9) = WDOT(9) -ROP
      WDOT(10) = WDOT(10) +ROP
      WDOT(14) = WDOT(14) +ROP
      ROP = RF(145)-RB(145)
      WDOT(9) = WDOT(9) -ROP
      WDOT(21) = WDOT(21) +ROP
      ROP = RF(146)-RB(146)
      WDOT(2) = WDOT(2) +ROP
      WDOT(9) = WDOT(9) -ROP
      WDOT(20) = WDOT(20) +ROP
      ROP = RF(147)-RB(147)
      WDOT(9) = WDOT(9) +ROP
      WDOT(11) = WDOT(11) +ROP
      ROP = RF(148)-RB(148)
      WDOT(2) = WDOT(2) -ROP
      WDOT(19) = WDOT(19) +ROP
      ROP = RF(149)-RB(149)
      WDOT(2) = WDOT(2) -ROP
      WDOT(9) = WDOT(9) +ROP
      ROP = RF(150)-RB(150)
      WDOT(1) = WDOT(1) +ROP
      WDOT(2) = WDOT(2) -ROP
      WDOT(18) = WDOT(18) +ROP
      ROP = RF(151)-RB(151)
      WDOT(3) = WDOT(3) -ROP
      WDOT(5) = WDOT(5) +ROP
      WDOT(18) = WDOT(18) +ROP
      ROP = RF(152)-RB(152)
      WDOT(5) = WDOT(5) -ROP
      WDOT(6) = WDOT(6) +ROP
      WDOT(18) = WDOT(18) +ROP
      ROP = RF(153)-RB(153)
      WDOT(4) = WDOT(4) -ROP
      WDOT(7) = WDOT(7) +ROP
      WDOT(18) = WDOT(18) +ROP
      ROP = RF(154)-RB(154)
      WDOT(4) = WDOT(4) -ROP
      WDOT(5) = WDOT(5) +ROP
      WDOT(11) = WDOT(11) +ROP
      WDOT(13) = WDOT(13) +ROP
      ROP = RF(155)-RB(155)
      WDOT(1) = WDOT(1) +ROP
      WDOT(15) = WDOT(15) -ROP
      ROP = RF(156)-RB(156)
      WDOT(2) = WDOT(2) -ROP
      WDOT(15) = WDOT(15) -ROP
      ROP = RF(157)-RB(157)
      WDOT(1) = WDOT(1) +ROP
      WDOT(2) = WDOT(2) -ROP
      WDOT(15) = WDOT(15) -ROP
      ROP = RF(158)-RB(158)
      WDOT(3) = WDOT(3) -ROP
      WDOT(5) = WDOT(5) +ROP
      WDOT(15) = WDOT(15) -ROP
      ROP = RF(159)-RB(159)
      WDOT(3) = WDOT(3) -ROP
      WDOT(9) = WDOT(9) +ROP
      WDOT(15) = WDOT(15) -ROP
      ROP = RF(160)-RB(160)
      WDOT(3) = WDOT(3) -ROP
      WDOT(13) = WDOT(13) +ROP
      WDOT(15) = WDOT(15) -ROP
      ROP = RF(161)-RB(161)
      WDOT(5) = WDOT(5) -ROP
      WDOT(6) = WDOT(6) +ROP
      WDOT(15) = WDOT(15) -ROP
      ROP = RF(162)-RB(162)
      WDOT(4) = WDOT(4) -ROP
      WDOT(7) = WDOT(7) +ROP
      WDOT(15) = WDOT(15) -ROP
      ROP = RF(163)-RB(163)
      WDOT(5) = WDOT(5) +ROP
      WDOT(7) = WDOT(7) -ROP
      WDOT(15) = WDOT(15) -ROP
      WDOT(19) = WDOT(19) +ROP
      ROP = RF(164)-RB(164)
      WDOT(11) = WDOT(11) +ROP
      WDOT(15) = WDOT(15) -ROP
      ROP = RF(165)-RB(165)
      WDOT(2) = WDOT(2) +ROP
      WDOT(15) = WDOT(15) -ROP
      WDOT(20) = WDOT(20) +ROP
      ROP = RF(166)-RB(166)
      WDOT(10) = WDOT(10) +ROP
      WDOT(15) = WDOT(15) -ROP
      ROP = RF(167)-RB(167)
      WDOT(2) = WDOT(2) +ROP
      WDOT(15) = WDOT(15) -ROP
      WDOT(20) = WDOT(20) +ROP
      ROP = RF(168)-RB(168)
      WDOT(9) = WDOT(9) -ROP
      WDOT(10) = WDOT(10) +ROP
      WDOT(15) = WDOT(15) -ROP
      ROP = RF(169)-RB(169)
      WDOT(9) = WDOT(9) -ROP
      WDOT(15) = WDOT(15) -ROP
      ROP = RF(170)-RB(170)
      WDOT(2) = WDOT(2) -ROP
      WDOT(16) = WDOT(16) +ROP
      ROP = RF(171)-RB(171)
      WDOT(1) = WDOT(1) +ROP
      WDOT(2) = WDOT(2) -ROP
      WDOT(15) = WDOT(15) +ROP
      ROP = RF(172)-RB(172)
      WDOT(3) = WDOT(3) -ROP
      WDOT(9) = WDOT(9) +ROP
      WDOT(13) = WDOT(13) +ROP
      ROP = RF(173)-RB(173)
      WDOT(2) = WDOT(2) +ROP
      WDOT(3) = WDOT(3) -ROP
      WDOT(19) = WDOT(19) +ROP
      ROP = RF(174)-RB(174)
      WDOT(4) = WDOT(4) -ROP
      WDOT(7) = WDOT(7) +ROP
      WDOT(15) = WDOT(15) +ROP
      ROP = RF(175)-RB(175)
      WDOT(4) = WDOT(4) +ROP
      WDOT(7) = WDOT(7) -ROP
      WDOT(16) = WDOT(16) +ROP
      ROP = RF(176)-RB(176)
      WDOT(7) = WDOT(7) -ROP
      WDOT(8) = WDOT(8) +ROP
      WDOT(15) = WDOT(15) +ROP
      ROP = RF(177)-RB(177)
      WDOT(5) = WDOT(5) +ROP
      WDOT(7) = WDOT(7) -ROP
      WDOT(9) = WDOT(9) +ROP
      WDOT(13) = WDOT(13) +ROP
      ROP = RF(178)-RB(178)
      WDOT(7) = WDOT(7) +ROP
      WDOT(8) = WDOT(8) -ROP
      WDOT(16) = WDOT(16) +ROP
      ROP = RF(179)-RB(179)
      WDOT(11) = WDOT(11) +ROP
      WDOT(16) = WDOT(16) +ROP
      ROP = RF(180)-RB(180)
      WDOT(1) = WDOT(1) +ROP
      WDOT(2) = WDOT(2) -ROP
      WDOT(16) = WDOT(16) -ROP
      ROP = RF(181)-RB(181)
      WDOT(3) = WDOT(3) -ROP
      WDOT(5) = WDOT(5) +ROP
      WDOT(16) = WDOT(16) -ROP
      ROP = RF(182)-RB(182)
      WDOT(5) = WDOT(5) -ROP
      WDOT(6) = WDOT(6) +ROP
      WDOT(16) = WDOT(16) -ROP
      ROP = RF(183)-RB(183)
      WDOT(9) = WDOT(9) +ROP
      WDOT(16) = WDOT(16) -ROP
      ROP = RF(184)-RB(184)
      WDOT(9) = WDOT(9) -ROP
      WDOT(10) = WDOT(10) +ROP
      WDOT(16) = WDOT(16) -ROP
      ROP = RF(185)-RB(185)
      WDOT(2) = WDOT(2) -ROP
      WDOT(20) = WDOT(20) -ROP
      WDOT(21) = WDOT(21) +ROP
      ROP = RF(186)-RB(186)
      WDOT(2) = WDOT(2) -ROP
      WDOT(10) = WDOT(10) +ROP
      WDOT(20) = WDOT(20) -ROP
      ROP = RF(187)-RB(187)
      WDOT(4) = WDOT(4) +ROP
      WDOT(7) = WDOT(7) -ROP
      WDOT(20) = WDOT(20) -ROP
      WDOT(21) = WDOT(21) +ROP
      ROP = RF(188)-RB(188)
      WDOT(5) = WDOT(5) +ROP
      WDOT(7) = WDOT(7) -ROP
      WDOT(13) = WDOT(13) +ROP
      WDOT(20) = WDOT(20) -ROP
      ROP = RF(189)-RB(189)
      WDOT(11) = WDOT(11) +ROP
      WDOT(20) = WDOT(20) -ROP
      WDOT(21) = WDOT(21) +ROP
      ROP = RF(190)-RB(190)
      WDOT(2) = WDOT(2) -ROP
      WDOT(21) = WDOT(21) -ROP
      ROP = RF(191)-RB(191)
      WDOT(2) = WDOT(2) -ROP
      WDOT(9) = WDOT(9) +ROP
      WDOT(15) = WDOT(15) +ROP
      WDOT(21) = WDOT(21) -ROP
      ROP = RF(192)-RB(192)
      WDOT(1) = WDOT(1) +ROP
      WDOT(2) = WDOT(2) -ROP
      WDOT(20) = WDOT(20) +ROP
      WDOT(21) = WDOT(21) -ROP
      ROP = RF(193)-RB(193)
      WDOT(2) = WDOT(2) +ROP
      WDOT(3) = WDOT(3) -ROP
      WDOT(9) = WDOT(9) +ROP
      WDOT(18) = WDOT(18) +ROP
      WDOT(21) = WDOT(21) -ROP
      ROP = RF(194)-RB(194)
      WDOT(3) = WDOT(3) -ROP
      WDOT(21) = WDOT(21) -ROP
      ROP = RF(195)-RB(195)
      WDOT(3) = WDOT(3) -ROP
      WDOT(5) = WDOT(5) +ROP
      WDOT(20) = WDOT(20) +ROP
      WDOT(21) = WDOT(21) -ROP
      ROP = RF(196)-RB(196)
      WDOT(5) = WDOT(5) -ROP
      WDOT(6) = WDOT(6) +ROP
      WDOT(20) = WDOT(20) +ROP
      WDOT(21) = WDOT(21) -ROP
      ROP = RF(197)-RB(197)
      WDOT(7) = WDOT(7) -ROP
      WDOT(8) = WDOT(8) +ROP
      WDOT(20) = WDOT(20) +ROP
      WDOT(21) = WDOT(21) -ROP
      ROP = RF(198)-RB(198)
      WDOT(9) = WDOT(9) -ROP
      WDOT(10) = WDOT(10) +ROP
      WDOT(20) = WDOT(20) +ROP
      WDOT(21) = WDOT(21) -ROP
      ROP = RF(199)-RB(199)
      WDOT(2) = WDOT(2) -ROP
      WDOT(9) = WDOT(9) +ROP
      ROP = RF(200)-RB(200)
      WDOT(1) = WDOT(1) +ROP
      WDOT(2) = WDOT(2) -ROP
      WDOT(21) = WDOT(21) +ROP
      ROP = RF(201)-RB(201)
      WDOT(3) = WDOT(3) -ROP
      WDOT(13) = WDOT(13) +ROP
      ROP = RF(202)-RB(202)
      WDOT(5) = WDOT(5) -ROP
      WDOT(6) = WDOT(6) +ROP
      WDOT(21) = WDOT(21) +ROP
      ROP = RF(203)-RB(203)
      WDOT(4) = WDOT(4) -ROP
      WDOT(7) = WDOT(7) +ROP
      WDOT(21) = WDOT(21) +ROP
      ROP = RF(204)-RB(204)
      WDOT(5) = WDOT(5) +ROP
      WDOT(7) = WDOT(7) -ROP
      WDOT(13) = WDOT(13) +ROP
      ROP = RF(205)-RB(205)
      WDOT(9) = WDOT(9) -ROP
      WDOT(10) = WDOT(10) +ROP
      WDOT(21) = WDOT(21) +ROP
      ROP = RF(206)-RB(206)
      WDOT(9) = WDOT(9) +ROP
      WDOT(20) = WDOT(20) +ROP
C
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE QSSA(RF, RB, XQ)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION RF(*), RB(*), XQ(*)
C
      RF(57) = 0.D0
      RF(58) = 0.D0
      RF(143) = 0.D0
      RF(179) = 0.D0
      RB(194) = 0.D0
      RF(206) = 0.D0
C
C     CH
      DEN = +RF( 34) +RF( 35) +RF( 36) +RF( 37) +RF( 38) 
     *  +RF( 39) +RF( 40) +RF( 77) +RF( 87) +RF(105) +RF(111) 
     *  +RB( 54) +RB( 60) 
      A1_0 = ( +RB( 34) +RB( 37) +RB( 39) +RB( 57) +RB( 77) 
     *  +RB(105) +RB(111) )/DEN
      A1_2 = ( +RB( 36) +RF( 54) )/DEN
      A1_3 = ( +RF( 60) )/DEN
      A1_4 = ( +RB( 35) +RB( 38) +RB( 40) )/DEN
      A1_7 = ( +RB( 87) )/DEN
C     CH2
      DEN = +RF( 48) +RF( 49) +RF( 50) +RF( 51) +RF( 52) 
     *  +RF( 53) +RF( 54) +RF( 55) +RF( 56) +RF( 91) +RF(106) 
     *  +RF(112) +RF(165) +RB( 36) +RB( 59) +RB( 67) +RB( 68) 
     *  +RB( 69) +RB( 80) +RB(117) +RB(123) +RB(125) +RB(130) 
     *  +RB(160) 
      A2_0 = ( +RB( 48) +RB( 49) +RB( 52) +RB( 53) +RB( 55) 
     *  +RB( 56) +RB( 57) +RB( 58) +RB( 58) +RF( 80) +RB( 91) 
     *  +RB(106) +RF(117) +RF(130) +RF(160) +RB(165) )/DEN
      A2_1 = ( +RF( 36) +RB( 54) )/DEN
      A2_3 = ( +RF( 59) +RF( 67) +RF( 68) +RF( 69) )/DEN
      A2_4 = ( +RB( 50) +RB( 51) )/DEN
      A2_6 = ( +RF(123) +RF(125) )/DEN
      A2_7 = ( +RB(112) )/DEN
C     CH2*
      DEN = +RF( 59) +RF( 60) +RF( 61) +RF( 62) +RF( 63) 
     *  +RF( 64) +RF( 65) +RF( 66) +RF( 67) +RF( 68) +RF( 69) 
     *  +RF( 70) +RF( 92) +RF(107) +RF(166) +RF(167) +RF(183) 
     *  +RB( 81) +RB( 98) +RB(108) 
      A3_0 = ( +RB( 61) +RB( 63) +RB( 64) +RB( 65) +RB( 66) 
     *  +RB( 70) +RF( 81) +RB( 92) +RB(107) +RF(108) +RB(167) )/DEN
      A3_1 = ( +RB( 60) )/DEN
      A3_2 = ( +RB( 59) +RB( 67) +RB( 68) +RB( 69) )/DEN
      A3_4 = ( +RB( 62) )/DEN
      A3_5 = ( +RF( 98) )/DEN
      A3_6 = ( +RB(166) )/DEN
      A3_8 = ( +RB(183) )/DEN
C     HCO
      DEN = +RF( 41) +RF( 42) +RF( 43) +RF( 44) +RF( 45) 
     *  +RF( 46) +RF( 47) +RF( 88) +RF( 89) +RF(120) +RF(164) 
     *  +RF(189) +RB( 35) +RB( 38) +RB( 40) +RB( 50) +RB( 51) 
     *  +RB( 62) +RB( 72) +RB( 73) +RB( 74) +RB( 75) +RB( 76) 
     *  +RB( 90) +RB(140) +RB(149) +RB(159) 
      A4_0 = ( +RB( 41) +RB( 42) +RB( 43) +RB( 44) +RB( 45) 
     *  +RB( 46) +RB( 47) +RF( 72) +RF( 73) +RF( 74) +RF( 75) 
     *  +RF( 76) +RB( 88) +RB( 89) +RF( 90) +RB(143) +RF(159) 
     *  +RB(179) +RB(189) +RF(194) )/DEN
      A4_1 = ( +RF( 35) +RF( 38) +RF( 40) )/DEN
      A4_2 = ( +RF( 50) +RF( 51) )/DEN
      A4_3 = ( +RF( 62) )/DEN
      A4_7 = ( +RB(120) +RF(140) )/DEN
      A4_8 = ( +RB(164) )/DEN
      A4_9 = ( +RF(149) )/DEN
C     CH3O
      DEN = +RF( 96) +RF( 97) +RF( 98) +RF( 99) +RF(100) 
     *  +RF(101) +RB( 71) +RB( 82) +RB( 85) 
      A5_0 = ( +RF( 71) +RF( 82) +RF( 85) +RB( 96) +RB( 97) 
     *  +RB( 99) +RB(100) +RB(101) )/DEN
      A5_3 = ( +RB( 98) )/DEN
C     H2CC
      DEN = +RF(122) +RF(123) +RF(124) +RF(125) +RB(114) 
     *  +RB(134) +RB(155) +RB(166) +RB(186) 
      A6_0 = ( +RF(114) +RB(122) +RB(124) +RF(155) +RF(186) )/DEN
      A6_2 = ( +RB(123) +RB(125) )/DEN
      A6_3 = ( +RF(166) )/DEN
      A6_7 = ( +RF(134) )/DEN
C     C2H3
      DEN = +RF(115) +RF(132) +RF(133) +RF(134) +RF(135) 
     *  +RF(136) +RF(137) +RF(138) +RF(139) +RF(140) +RF(141) 
     *  +RF(142) +RF(144) +RF(145) +RF(146) +RB( 87) +RB(112) 
     *  +RB(120) +RB(157) +RB(158) +RB(161) +RB(162) +RB(168) 
     *  +RB(188) 
      A7_0 = ( +RB(115) +RB(132) +RB(133) +RB(135) +RB(136) 
     *  +RB(137) +RB(138) +RB(142) +RB(143) +RB(144) +RB(145) 
     *  +RB(146) +RF(157) +RF(158) +RF(161) +RF(162) +RF(168) 
     *  +RF(188) +RB(206) )/DEN
      A7_1 = ( +RF( 87) )/DEN
      A7_2 = ( +RF(112) )/DEN
      A7_4 = ( +RF(120) +RB(140) )/DEN
      A7_6 = ( +RB(134) )/DEN
      A7_9 = ( +RB(139) +RB(141) )/DEN
C     C2H5
      DEN = +RF(170) +RF(171) +RF(172) +RF(173) +RF(174) 
     *  +RF(175) +RF(176) +RF(177) +RF(178) +RB( 94) +RB(156) 
     *  +RB(164) +RB(180) +RB(181) +RB(182) +RB(183) +RB(184) 
     *  +RB(199) +RB(201) +RB(204) 
      A8_0 = ( +RF( 94) +RF(156) +RB(170) +RB(171) +RB(172) 
     *  +RB(173) +RB(174) +RB(175) +RB(176) +RB(177) +RB(178) 
     *  +RB(179) +RF(180) +RF(181) +RF(182) +RF(184) +RF(194) 
     *  +RB(206) )/DEN
      A8_3 = ( +RF(183) )/DEN
      A8_4 = ( +RF(164) )/DEN
      A8_10 = ( +RF(199) +RF(201) +RF(204) )/DEN
C     CH2CHO
      DEN = +RF(147) +RF(148) +RF(149) +RF(150) +RF(151) 
     *  +RF(152) +RF(153) +RF(154) +RB(126) +RB(139) +RB(141) 
      A9_0 = ( +RF(126) +RB(147) +RB(148) +RB(150) +RB(151) 
     *  +RB(152) +RB(153) +RB(154) )/DEN
      A9_4 = ( +RB(149) )/DEN
      A9_7 = ( +RF(139) +RF(141) )/DEN
C     nC3H7
      DEN = +RF(199) +RF(200) +RF(201) +RF(202) +RF(203) 
     *  +RF(204) +RF(205) +RB(169) +RB(190) 
      A10_0 = ( +RF(169) +RF(190) +RB(200) +RB(202) +RB(203) 
     *  +RB(205) )/DEN
      A10_8 = ( +RB(199) +RB(201) +RB(204) )/DEN
C
      A3_0 = A3_0 + A3_5*A5_0
      DEN = 1 -A3_5*A5_3
      A3_0 = A3_0/DEN
      A3_4 = A3_4/DEN
      A3_2 = A3_2/DEN
      A3_1 = A3_1/DEN
      A3_6 = A3_6/DEN
      A3_8 = A3_8/DEN
      A8_0 = A8_0 + A8_10*A10_0
      DEN = 1 -A8_10*A10_8
      A8_0 = A8_0/DEN
      A8_3 = A8_3/DEN
      A8_4 = A8_4/DEN
      A4_0 = A4_0 + A4_9*A9_0
      A4_7 = A4_7 + A4_9*A9_7
      DEN = 1 -A4_9*A9_4
      A4_0 = A4_0/DEN
      A4_3 = A4_3/DEN
      A4_7 = A4_7/DEN
      A4_2 = A4_2/DEN
      A4_1 = A4_1/DEN
      A4_8 = A4_8/DEN
      A7_0 = A7_0 + A7_9*A9_0
      A7_4 = A7_4 + A7_9*A9_4
      DEN = 1 -A7_9*A9_7
      A7_0 = A7_0/DEN
      A7_4 = A7_4/DEN
      A7_2 = A7_2/DEN
      A7_1 = A7_1/DEN
      A7_6 = A7_6/DEN
      A3_0 = A3_0 + A3_8*A8_0
      A3_4 = A3_4 + A3_8*A8_4
      DEN = 1 -A3_8*A8_3
      A3_0 = A3_0/DEN
      A3_4 = A3_4/DEN
      A3_2 = A3_2/DEN
      A3_1 = A3_1/DEN
      A3_6 = A3_6/DEN
      A4_0 = A4_0 + A4_8*A8_0
      A4_3 = A4_3 + A4_8*A8_3
      DEN = 1 -A4_8*A8_4
      A4_0 = A4_0/DEN
      A4_3 = A4_3/DEN
      A4_7 = A4_7/DEN
      A4_2 = A4_2/DEN
      A4_1 = A4_1/DEN
      A3_0 = A3_0 + A3_6*A6_0
      A3_7 = A3_6*A6_7
      A3_2 = A3_2 + A3_6*A6_2
      DEN = 1 -A3_6*A6_3
      A3_0 = A3_0/DEN
      A3_4 = A3_4/DEN
      A3_7 = A3_7/DEN
      A3_2 = A3_2/DEN
      A3_1 = A3_1/DEN
      A7_0 = A7_0 + A7_6*A6_0
      A7_3 = A7_6*A6_3
      A7_2 = A7_2 + A7_6*A6_2
      DEN = 1 -A7_6*A6_7
      A7_0 = A7_0/DEN
      A7_3 = A7_3/DEN
      A7_4 = A7_4/DEN
      A7_2 = A7_2/DEN
      A7_1 = A7_1/DEN
      A2_0 = A2_0 + A2_6*A6_0
      A2_3 = A2_3 + A2_6*A6_3
      A2_7 = A2_7 + A2_6*A6_7
      DEN = 1 -A2_6*A6_2
      A2_0 = A2_0/DEN
      A2_3 = A2_3/DEN
      A2_4 = A2_4/DEN
      A2_7 = A2_7/DEN
      A2_1 = A2_1/DEN
      A3_0 = A3_0 + A3_1*A1_0
      A3_4 = A3_4 + A3_1*A1_4
      A3_7 = A3_7 + A3_1*A1_7
      A3_2 = A3_2 + A3_1*A1_2
      DEN = 1 -A3_1*A1_3
      A3_0 = A3_0/DEN
      A3_4 = A3_4/DEN
      A3_7 = A3_7/DEN
      A3_2 = A3_2/DEN
      A4_0 = A4_0 + A4_1*A1_0
      A4_3 = A4_3 + A4_1*A1_3
      A4_7 = A4_7 + A4_1*A1_7
      A4_2 = A4_2 + A4_1*A1_2
      DEN = 1 -A4_1*A1_4
      A4_0 = A4_0/DEN
      A4_3 = A4_3/DEN
      A4_7 = A4_7/DEN
      A4_2 = A4_2/DEN
      A7_0 = A7_0 + A7_1*A1_0
      A7_3 = A7_3 + A7_1*A1_3
      A7_4 = A7_4 + A7_1*A1_4
      A7_2 = A7_2 + A7_1*A1_2
      DEN = 1 -A7_1*A1_7
      A7_0 = A7_0/DEN
      A7_3 = A7_3/DEN
      A7_4 = A7_4/DEN
      A7_2 = A7_2/DEN
      A2_0 = A2_0 + A2_1*A1_0
      A2_3 = A2_3 + A2_1*A1_3
      A2_4 = A2_4 + A2_1*A1_4
      A2_7 = A2_7 + A2_1*A1_7
      DEN = 1 -A2_1*A1_2
      A2_0 = A2_0/DEN
      A2_3 = A2_3/DEN
      A2_4 = A2_4/DEN
      A2_7 = A2_7/DEN
      A3_0 = A3_0 + A3_2*A2_0
      A3_4 = A3_4 + A3_2*A2_4
      A3_7 = A3_7 + A3_2*A2_7
      DEN = 1 -A3_2*A2_3
      A3_0 = A3_0/DEN
      A3_4 = A3_4/DEN
      A3_7 = A3_7/DEN
      A4_0 = A4_0 + A4_2*A2_0
      A4_3 = A4_3 + A4_2*A2_3
      A4_7 = A4_7 + A4_2*A2_7
      DEN = 1 -A4_2*A2_4
      A4_0 = A4_0/DEN
      A4_3 = A4_3/DEN
      A4_7 = A4_7/DEN
      A7_0 = A7_0 + A7_2*A2_0
      A7_3 = A7_3 + A7_2*A2_3
      A7_4 = A7_4 + A7_2*A2_4
      DEN = 1 -A7_2*A2_7
      A7_0 = A7_0/DEN
      A7_3 = A7_3/DEN
      A7_4 = A7_4/DEN
      A3_0 = A3_0 + A3_7*A7_0
      A3_4 = A3_4 + A3_7*A7_4
      DEN = 1 -A3_7*A7_3
      A3_0 = A3_0/DEN
      A3_4 = A3_4/DEN
      A4_0 = A4_0 + A4_7*A7_0
      A4_3 = A4_3 + A4_7*A7_3
      DEN = 1 -A4_7*A7_4
      A4_0 = A4_0/DEN
      A4_3 = A4_3/DEN
      A3_0 = A3_0 + A3_4*A4_0
      DEN = 1 -A3_4*A4_3
      A3_0 = A3_0/DEN
      XQ(3) = A3_0
      XQ(4) = A4_0 +A4_3*XQ(3)
      XQ(7) = A7_0 +A7_3*XQ(3) +A7_4*XQ(4)
      XQ(2) = A2_0 +A2_3*XQ(3) +A2_4*XQ(4) +A2_7*XQ(7)
      XQ(1) = A1_0 +A1_3*XQ(3) +A1_4*XQ(4) +A1_7*XQ(7) +A1_2*XQ(2)
      XQ(6) = A6_0 +A6_3*XQ(3) +A6_7*XQ(7) +A6_2*XQ(2)
      XQ(8) = A8_0 +A8_3*XQ(3) +A8_4*XQ(4)
      XQ(9) = A9_0 +A9_4*XQ(4) +A9_7*XQ(7)
      XQ(10) = A10_0 +A10_8*XQ(8)
      XQ(5) = A5_0 +A5_3*XQ(3)
C
      RF( 34) = RF( 34)*XQ( 1)
      RF( 35) = RF( 35)*XQ( 1)
      RB( 35) = RB( 35)*XQ( 4)
      RF( 36) = RF( 36)*XQ( 1)
      RB( 36) = RB( 36)*XQ( 2)
      RF( 37) = RF( 37)*XQ( 1)
      RF( 38) = RF( 38)*XQ( 1)
      RB( 38) = RB( 38)*XQ( 4)
      RF( 39) = RF( 39)*XQ( 1)
      RF( 40) = RF( 40)*XQ( 1)
      RB( 40) = RB( 40)*XQ( 4)
      RF( 41) = RF( 41)*XQ( 4)
      RF( 42) = RF( 42)*XQ( 4)
      RF( 43) = RF( 43)*XQ( 4)
      RF( 44) = RF( 44)*XQ( 4)
      RF( 45) = RF( 45)*XQ( 4)
      RF( 46) = RF( 46)*XQ( 4)
      RF( 47) = RF( 47)*XQ( 4)
      RF( 48) = RF( 48)*XQ( 2)
      RF( 49) = RF( 49)*XQ( 2)
      RF( 50) = RF( 50)*XQ( 2)
      RB( 50) = RB( 50)*XQ( 4)
      RF( 51) = RF( 51)*XQ( 2)
      RB( 51) = RB( 51)*XQ( 4)
      RF( 52) = RF( 52)*XQ( 2)
      RF( 53) = RF( 53)*XQ( 2)
      RF( 54) = RF( 54)*XQ( 2)
      RB( 54) = RB( 54)*XQ( 1)
      RF( 55) = RF( 55)*XQ( 2)
      RF( 56) = RF( 56)*XQ( 2)
      RF( 59) = RF( 59)*XQ( 3)
      RB( 59) = RB( 59)*XQ( 2)
      RF( 60) = RF( 60)*XQ( 3)
      RB( 60) = RB( 60)*XQ( 1)
      RF( 61) = RF( 61)*XQ( 3)
      RF( 62) = RF( 62)*XQ( 3)
      RB( 62) = RB( 62)*XQ( 4)
      RF( 63) = RF( 63)*XQ( 3)
      RF( 64) = RF( 64)*XQ( 3)
      RF( 65) = RF( 65)*XQ( 3)
      RF( 66) = RF( 66)*XQ( 3)
      RF( 67) = RF( 67)*XQ( 3)
      RB( 67) = RB( 67)*XQ( 2)
      RF( 68) = RF( 68)*XQ( 3)
      RB( 68) = RB( 68)*XQ( 2)
      RF( 69) = RF( 69)*XQ( 3)
      RB( 69) = RB( 69)*XQ( 2)
      RF( 70) = RF( 70)*XQ( 3)
      RB( 71) = RB( 71)*XQ( 5)
      RB( 72) = RB( 72)*XQ( 4)
      RB( 73) = RB( 73)*XQ( 4)
      RB( 74) = RB( 74)*XQ( 4)
      RB( 75) = RB( 75)*XQ( 4)
      RB( 76) = RB( 76)*XQ( 4)
      RF( 77) = RF( 77)*XQ( 1)
      RB( 80) = RB( 80)*XQ( 2)
      RB( 81) = RB( 81)*XQ( 3)
      RB( 82) = RB( 82)*XQ( 5)
      RB( 85) = RB( 85)*XQ( 5)
      RF( 87) = RF( 87)*XQ( 1)
      RB( 87) = RB( 87)*XQ( 7)
      RF( 88) = RF( 88)*XQ( 4)
      RF( 89) = RF( 89)*XQ( 4)
      RB( 90) = RB( 90)*XQ( 4)
      RF( 91) = RF( 91)*XQ( 2)
      RF( 92) = RF( 92)*XQ( 3)
      RB( 94) = RB( 94)*XQ( 8)
      RF( 96) = RF( 96)*XQ( 5)
      RF( 97) = RF( 97)*XQ( 5)
      RF( 98) = RF( 98)*XQ( 5)
      RB( 98) = RB( 98)*XQ( 3)
      RF( 99) = RF( 99)*XQ( 5)
      RF(100) = RF(100)*XQ( 5)
      RF(101) = RF(101)*XQ( 5)
      RF(105) = RF(105)*XQ( 1)
      RF(106) = RF(106)*XQ( 2)
      RF(107) = RF(107)*XQ( 3)
      RB(108) = RB(108)*XQ( 3)
      RF(111) = RF(111)*XQ( 1)
      RF(112) = RF(112)*XQ( 2)
      RB(112) = RB(112)*XQ( 7)
      RB(114) = RB(114)*XQ( 6)
      RF(115) = RF(115)*XQ( 7)
      RB(117) = RB(117)*XQ( 2)
      RF(120) = RF(120)*XQ( 4)
      RB(120) = RB(120)*XQ( 7)
      RF(122) = RF(122)*XQ( 6)
      RF(123) = RF(123)*XQ( 6)
      RB(123) = RB(123)*XQ( 2)
      RF(124) = RF(124)*XQ( 6)
      RF(125) = RF(125)*XQ( 6)
      RB(125) = RB(125)*XQ( 2)
      RB(126) = RB(126)*XQ( 9)
      RB(130) = RB(130)*XQ( 2)
      RF(132) = RF(132)*XQ( 7)
      RF(133) = RF(133)*XQ( 7)
      RF(134) = RF(134)*XQ( 7)
      RB(134) = RB(134)*XQ( 6)
      RF(135) = RF(135)*XQ( 7)
      RF(136) = RF(136)*XQ( 7)
      RF(137) = RF(137)*XQ( 7)
      RF(138) = RF(138)*XQ( 7)
      RF(139) = RF(139)*XQ( 7)
      RB(139) = RB(139)*XQ( 9)
      RF(140) = RF(140)*XQ( 7)
      RB(140) = RB(140)*XQ( 4)
      RF(141) = RF(141)*XQ( 7)
      RB(141) = RB(141)*XQ( 9)
      RF(142) = RF(142)*XQ( 7)
      RF(144) = RF(144)*XQ( 7)
      RF(145) = RF(145)*XQ( 7)
      RF(146) = RF(146)*XQ( 7)
      RF(147) = RF(147)*XQ( 9)
      RF(148) = RF(148)*XQ( 9)
      RF(149) = RF(149)*XQ( 9)
      RB(149) = RB(149)*XQ( 4)
      RF(150) = RF(150)*XQ( 9)
      RF(151) = RF(151)*XQ( 9)
      RF(152) = RF(152)*XQ( 9)
      RF(153) = RF(153)*XQ( 9)
      RF(154) = RF(154)*XQ( 9)
      RB(155) = RB(155)*XQ( 6)
      RB(156) = RB(156)*XQ( 8)
      RB(157) = RB(157)*XQ( 7)
      RB(158) = RB(158)*XQ( 7)
      RB(159) = RB(159)*XQ( 4)
      RB(160) = RB(160)*XQ( 2)
      RB(161) = RB(161)*XQ( 7)
      RB(162) = RB(162)*XQ( 7)
      RF(164) = RF(164)*XQ( 4)
      RB(164) = RB(164)*XQ( 8)
      RF(165) = RF(165)*XQ( 2)
      RF(166) = RF(166)*XQ( 3)
      RB(166) = RB(166)*XQ( 6)
      RF(167) = RF(167)*XQ( 3)
      RB(168) = RB(168)*XQ( 7)
      RB(169) = RB(169)*XQ(10)
      RF(170) = RF(170)*XQ( 8)
      RF(171) = RF(171)*XQ( 8)
      RF(172) = RF(172)*XQ( 8)
      RF(173) = RF(173)*XQ( 8)
      RF(174) = RF(174)*XQ( 8)
      RF(175) = RF(175)*XQ( 8)
      RF(176) = RF(176)*XQ( 8)
      RF(177) = RF(177)*XQ( 8)
      RF(178) = RF(178)*XQ( 8)
      RB(180) = RB(180)*XQ( 8)
      RB(181) = RB(181)*XQ( 8)
      RB(182) = RB(182)*XQ( 8)
      RF(183) = RF(183)*XQ( 3)
      RB(183) = RB(183)*XQ( 8)
      RB(184) = RB(184)*XQ( 8)
      RB(186) = RB(186)*XQ( 6)
      RB(188) = RB(188)*XQ( 7)
      RF(189) = RF(189)*XQ( 4)
      RB(190) = RB(190)*XQ(10)
      RF(199) = RF(199)*XQ(10)
      RB(199) = RB(199)*XQ( 8)
      RF(200) = RF(200)*XQ(10)
      RF(201) = RF(201)*XQ(10)
      RB(201) = RB(201)*XQ( 8)
      RF(202) = RF(202)*XQ(10)
      RF(203) = RF(203)*XQ(10)
      RF(204) = RF(204)*XQ(10)
      RB(204) = RB(204)*XQ( 8)
      RF(205) = RF(205)*XQ(10)
C
      END
