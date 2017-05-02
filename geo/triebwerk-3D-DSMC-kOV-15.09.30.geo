//Variablen //Nicht verändern! Unten existiert eine Variable unabhängig von diesen Werten
D1 =  	2.4/1000; 	// Durchmesser Zündkerze [m]
D2 = 	1.8/1000; 		// Durchmesser Verbindung Eintritt-Defusor
D3 = 	1/1000; 		// Durchmesser Kathode-Spitze [m] (Aus Origianl)
D4 = 	16/1000; 		//Durchmesser Inlet [m]
D5 = 	12/1000; 		// Durchmesser Inlet nach Senkung [m]
D6 = 	80/1000;	// Durchmesser Austrittsvolumen [m]

L1 = 	4.84529946162/1000; 	//Eintrittslänge vor Beginn der Krümmung [m]
L2 = 	5.93211598/1000; 		//Länge in X-Richtung der Krümmung [m]
L3 = 2.63/1000; 				// Länge der Kathode vom Beginn der Krümmung (ohne Spitze) [m]	
L4 = 	60/1000;				// Länge Austrittsvolumen [m]
L5 = 	2/1000;					// Länge Düsenhals [m]
L6 = 	20/1000;				//Länge Defusor
L8 = 	1.15470053838/1000;			//Tiefe der Senkung in X-Richtung
L9 =	8/1000;					//Zusatzlänge Eintritt
L10 = 1/1000; // Kathodenspitze Länge
Alpha = 20; 					// Öffnungswinkel [°]

dsmcFac = 0.5;					//Anpassung der Eintrittslänge 0<dsmcFac<1; 1--> Vollständiges Inlet

//Variablen End



KathodenSpitzenWinkel=ArcTan(L10/((D1-D3)/2))*360/(2*Pi);
md=0.001;


X1=0;
Y1=0;
Z1=0;

AX1=X1;
AY1=Y1;
AZ1=Z1;

AY2=AY1+D4/2;
AZ2=AZ1+D4/2;

//A außen
Point(1)={AX1,AY1,AZ1,md};
Point(2)={0,0,AZ2,md};
Point(3)={0,AY2,0,md};
Point(4)={0,0,-AZ2,md};
Point(5)={0,-AY2,0,md};

BX1=AX1+L9;
BY1=AY1;
BZ1=AZ1;

BX2=BX1;
BY2=BY1+D4/2;
BZ2=BZ1+D4/2;

//B außen
Point(6)={BX1,BY1,BZ1,md};
Point(7)={BX2,0,BZ2,md};
Point(8)={BX2,BY2,0,md};
Point(9)={BX2,0,-BZ2,md};
Point(10)={BX2,-BY2,0,md};


CX1=BX1+L8;
CY1=BY1;
CZ1=BZ1;

//Inlet_Mod DSMC 15.09.10
CXDSMC=CX1+L1-dsmcFac*L1;

CX2=CX1;
CY2=CY1+D5/2;
CZ2=CZ1+D5/2;
//DSMC-Mod	30.09.15
Point(11)={CXDSMC,CY1,CZ1,md};
Point(12)={CXDSMC,0,CZ2,md};
Point(13)={CXDSMC,CY2,0,md};
Point(14)={CXDSMC,0,-CZ2,md};
Point(15)={CXDSMC,-CY2,0,md};

/*
Point(11)={CX1,CY1,CZ1,md};
Point(12)={CX2,0,CZ2,md};
Point(13)={CX2,CY2,0,md};
Point(14)={CX2,0,-CZ2,md};
Point(15)={CX2,-CY2,0,md};
*/
DX1=CX1+L1;
DY1=CY1;
DZ1=CZ1;

DX2=DX1;
DY2=DY1+D5/2;
DZ2=DZ1+D5/2;

//D außen
Point(16)={DX1,DY1,DZ1,md};
Point(17)={DX2,0,DZ2,md};
Point(18)={DX2,DY2,0,md};
Point(19)={DX2,0,-DZ2,md};
Point(20)={DX2,-DY2,0,md};


EX1=DX1+L2;
EY1=DY1;
EZ1=DZ1;

EX2=EX1;
EY2=EY1+D2/2;
EZ2=EZ1+D2/2;

//E außen
Point(21)={EX1,EY1,EZ1,md};
Point(22)={EX2,0,EZ2,md};
Point(23)={EX2,EY2,0,md};
Point(24)={EX2,0,-EZ2,md};
Point(25)={EX2,-EY2,0,md};


FX1=EX1+L5;
FY1=EY1;
FZ1=EZ1;

FX2=FX1;
FY2=FY1+D2/2;
FZ2=FZ1+D2/2;

//F außen
Point(26)={FX1,FY1,FZ1,md};
Point(27)={FX2,0,FZ2,md};
Point(28)={FX2,FY2,0,md};
Point(29)={FX2,0,-FZ2,md};
Point(30)={FX2,-FY2,0,md};


TanAlpha=Tan(Alpha*Pi/180);

GX1=FX1+L6;
GY1=FY1;
GZ1=FZ1;

GX2=GX1;

GY2=FY2+TanAlpha*L6;
GZ2=FZ2+TanAlpha*L6;

//G außen
Point(31)={GX1,GY1,GZ1,md};
Point(32)={GX2,0,GZ2,md};
Point(33)={GX2,GY2,0,md};
Point(34)={GX2,0,-GZ2,md};
Point(35)={GX2,-GY2,0,md};

HX1=GX1+L4;
HY1=GY1;
HZ1=GZ1;

HX2=HX1;
HY2=HY1+D6/2;
HZ2=HZ1+D6/2;

//Auf Höhe H
Point(36)={GX1,0,HZ2,md};
Point(37)={GX1,HY2,0,md};
Point(38)={GX1,0,-HZ2,md};
Point(39)={GX1,-HY2,0,md};

//H außen
Point(40)={HX1,HY1,HZ1,md};
Point(41)={HX2,0,HZ2,md};
Point(42)={HX2,HY2,0,md};
Point(43)={HX2,0,-HZ2,md};
Point(44)={HX2,-HY2,0,md};


//Linien außen

//Ring A-B
/*				DSMC - Mod 30.09.15
Line(1) = {2,7};
Line(2) = {3,8};
Line(3) = {4,9};
Line(4) = {5,10};
*/
//Ring B-C
/*				DSMC - Mod 30.09.15
Line(5) = {7,12};
Line(6) = {8,13};
Line(7) = {9,14};
Line(8) = {10,15};
*/
//Ring C-D
Line(9) = {12,17};
Line(10) = {13,18};
Line(11) = {14,19};
Line(12) = {15,20};
//Ring E-F
Line(13) = {22,27};
Line(14) = {23,28};
Line(15) = {24,29};
Line(16) = {25,30};
//Ring F-G
Line(17) = {27,32};
Line(18) = {28,33};
Line(19) = {29,34};
Line(20) = {30,35};
//G auf Höhe H
Line(21) = {32,36};
Line(22) = {33,37};
Line(23) = {34,38};
Line(24) = {35,39};
// => Ring G-H
Line(25) = {36,41};
Line(26) = {37,42};
Line(27) = {38,43};
Line(28) = {39,44};


//Circles außen

//Ebene A
/*				DSMC - Mod 30.09.15
Circle(29)={2,1,3};
Circle(30)={3,1,4};
Circle(31)={4,1,5};
Circle(32)={5,1,2};
*/
//Ebene B
/*				DSMC - Mod 30.09.15
Circle(33)={7,6,8};
Circle(34)={8,6,9};
Circle(35)={9,6,10};
Circle(36)={10,6,7};
*/
//Ebene C
Circle(37)={12,11,13};
Circle(38)={13,11,14};
Circle(39)={14,11,15};
Circle(40)={15,11,12};

//Ebene D
Circle(41)={17,16,18};
Circle(42)={18,16,19};
Circle(43)={19,16,20};
Circle(44)={20,16,17};

//Krümmung D-E
//Verschoben nach unten
/*
Circle(45)={17,16,22};
Circle(46)={18,16,23};
Circle(47)={19,16,24};
Circle(48)={20,16,25};
*/

//Ebene E
Circle(49)={22,21,23};
Circle(50)={23,21,24};
Circle(51)={24,21,25};
Circle(52)={25,21,22};

//Ebene F
Circle(53)={27,26,28};
Circle(54)={28,26,29};
Circle(55)={29,26,30};
Circle(56)={30,26,27};

//Ebene G
Circle(57)={32,31,33};
Circle(58)={33,31,34};
Circle(59)={34,31,35};
Circle(60)={35,31,32};

//Ebene G auf Höhe H
Circle(61)={36,31,37};
Circle(62)={37,31,38};
Circle(63)={38,31,39};
Circle(64)={39,31,36};

//Ebene H
Circle(65)={41,40,42};
Circle(66)={42,40,43};
Circle(67)={43,40,44};
Circle(68)={44,40,41};

////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
//Kathode
KAX=AX1;
KBX=BX1;
KCX=CX1;
KDX=DX1+L3;
KEX=KDX+L10;

KY=Y1+D1/2;
KZ=Z1+D1/2;

KSY=Y1+D3/2;
KSZ=Z1+D3/2;
// A
Point(45)={KAX,0,KZ,md};
Point(46)={KAX,KY,0,md};
Point(47)={KAX,0,-KZ,md};
Point(48)={KAX,-KY,0,md};
// B
Point(49)={KBX,0,KZ,md};
Point(50)={KBX,KY,0,md};
Point(51)={KBX,0,-KZ,md};
Point(52)={KBX,-KY,0,md};
// C
/*
Point(53)={KCX,0,KZ,md};
Point(54)={KCX,KY,0,md};
Point(55)={KCX,0,-KZ,md};
Point(56)={KCX,-KY,0,md};
*/
//DSMC-Mod 30.09.15
Point(53)={CXDSMC,0,KZ,md};
Point(54)={CXDSMC,KY,0,md};
Point(55)={CXDSMC,0,-KZ,md};
Point(56)={CXDSMC,-KY,0,md};

//D in der Ebene
Point(83)={KDX-L3,0,KZ,md};
Point(84)={KDX-L3,KY,0,md};
Point(85)={KDX-L3,0,-KZ,md};
Point(86)={KDX-L3,-KY,0,md};

// D bzw Beginn Spitze
Point(57)={KDX,0,0,md};
Point(58)={KDX,0,KZ,md};
Point(59)={KDX,KY,0,md};
Point(60)={KDX,0,-KZ,md};
Point(61)={KDX,-KY,0,md};
// E bzw Spitze
Point(62)={KEX,0,0,md};
Point(63)={KEX,0,KSZ,md};
Point(64)={KEX,KSY,0,md};
Point(65)={KEX,0,-KSZ,md};
Point(66)={KEX,-KSY,0,md};

//Linien Kathode
//=>A-B
/*				DSMC - Mod 30.09.15
Line(69)={45,49};
Line(70)={46,50};
Line(71)={47,51};
Line(72)={48,52};
*/
//=>B-C
/*
Line(73)={49,53};
Line(74)={50,54};
Line(75)={51,55};
Line(76)={52,56};
*/
//=>C-D
Line(77)={53,83};
Line(78)={54,84};
Line(79)={55,85};
Line(80)={56,86};

//=> D-D(Spitze)
Line(169)={83,58};
Line(170)={84,59};
Line(171)={85,60};
Line(172)={86,61};

/* Alt C-D(Spitze)
Line(77)={53,58};
Line(78)={54,59};
Line(79)={55,60};
Line(80)={56,61};
*/

//Circles Kathode
// A
/*				DSMC - Mod 30.09.15
Circle(81)={45,1,46};
Circle(82)={46,1,47};
Circle(83)={47,1,48};
Circle(84)={48,1,45};
*/
// B
/*				DSMC - Mod 30.09.15
Circle(85)={49,6,50};
Circle(86)={50,6,51};
Circle(87)={51,6,52};
Circle(88)={52,6,49};
*/
// C
Circle(89)={53,11,54};
Circle(90)={54,11,55};
Circle(91)={55,11,56};
Circle(92)={56,11,53};
//D-Ebene

Circle(149)={83,16,84};
Circle(150)={84,16,85};
Circle(151)={85,16,86};
Circle(152)={86,16,83};

// D-Spitze
Circle(93)={58,57,59};
Circle(94)={59,57,60};
Circle(95)={60,57,61};
Circle(96)={61,57,58};
// E
Circle(97)={63,62,64};
Circle(98)={64,62,65};
Circle(99)={65,62,66};
Circle(100)={66,62,63};

//Spitze
Line(101)={58,63};
Line(102)={59,64};
Line(103)={60,65};
Line(104)={61,66};

////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
//Mitte/Inneres


//Ebene A
AX3=AX1;
AY3=KY+0.25*(D4-D1);
AZ3=KZ+0.25*(D4-D1);
Point(67)={AX3,0,AZ3,md};
Point(68)={AX3,AY3,0,md};
Point(69)={AX3,0,-AZ3,md};
Point(70)={AX3,-AY3,0,md};

//Ebene B
BX3=BX1;
BY3=KY+0.25*(D4-D1);
BZ3=KZ+0.25*(D4-D1);
Point(71)={BX3,0,BZ3,md};
Point(72)={BX3,BY3,0,md};
Point(73)={BX3,0,-BZ3,md};
Point(74)={BX3,-BY3,0,md};

//Ebene C
CX3=CX1;
CY3=KY+0.25*(D5-D1);
CZ3=KZ+0.25*(D5-D1);

/*	DSMC Mod 30.09.15
Point(75)={CX3,0,CZ3,md};
Point(76)={CX3,CY3,0,md};
Point(77)={CX3,0,-CZ3,md};
Point(78)={CX3,-CY3,0,md};
*/
Point(75)={CXDSMC,0,CZ3,md};
Point(76)={CXDSMC,CY3,0,md};
Point(77)={CXDSMC,0,-CZ3,md};
Point(78)={CXDSMC,-CY3,0,md};

//Ebene D
DX3=DX1;
DY3=KY+0.25*(D5-D1);
DZ3=KZ+0.25*(D5-D1);
Point(79)={DX3,0,DZ3,md};
Point(80)={DX3,DY3,0,md};
Point(81)={DX3,0,-DZ3,md};
Point(82)={DX3,-DY3,0,md};

//Linien

//1.Verbindung in X-Richtung
//A-B
/*				DSMC - Mod 30.09.15
Line(105)={67,71};
Line(106)={68,72};
Line(107)={69,73};
Line(108)={70,74};
*/
//B-C
/*				DSMC - Mod 30.09.15
Line(109)={71,75};
Line(110)={72,76};
Line(111)={73,77};
Line(112)={74,78};
*/
//C-D
Line(113)={75,79};
Line(114)={76,80};
Line(115)={77,81};
Line(116)={78,82};
//2.Verbindung in YZ-Richtung
//Ebene A
/*				DSMC - Mod 30.09.15
Line(117)={2,67};
Line(118)={3,68};
Line(119)={4,69};
Line(120)={5,70};

Line(121)={67,45};
Line(122)={68,46};
Line(123)={69,47};
Line(124)={70,48};
*/

//Ebene B
/*				DSMC - Mod 30.09.15
Line(125)={7,71};
Line(126)={8,72};
Line(127)={9,73};
Line(128)={10,74};

Line(129)={71,49};
Line(130)={72,50};
Line(131)={73,51};
Line(132)={74,52};
*/
//Ebene C
Line(133)={12,75};
Line(134)={13,76};
Line(135)={14,77};
Line(136)={15,78};

Line(137)={75,53};
Line(138)={76,54};
Line(139)={77,55};
Line(140)={78,56};

//Ebene D
Line(141)={17,79};
Line(142)={18,80};
Line(143)={19,81};
Line(144)={20,82};

Line(145)={79,83};
Line(146)={80,84};
Line(147)={81,85};
Line(148)={82,86};

//Circles
//A mittig
/*				DSMC - Mod 30.09.15
Circle(153)={67,1,68};
Circle(154)={68,1,69};
Circle(155)={69,1,70};
Circle(156)={70,1,67};
*/

//B mittig
/*				DSMC - Mod 30.09.15
Circle(157)={71,6,72};
Circle(158)={72,6,73};
Circle(159)={73,6,74};
Circle(160)={74,6,71};
*/

//C mittig
Circle(161)={75,11,76};
Circle(162)={76,11,77};
Circle(163)={77,11,78};
Circle(164)={78,11,75};

//D mittig
Circle(165)={79,16,80};
Circle(166)={80,16,81};
Circle(167)={81,16,82};
Circle(168)={82,16,79};



//Rechtecke A-B
/*				DSMC - Mod 30.09.15
Line Loop(1) = {106, 130, -70, -122};
Ruled Surface(1) = {1};
Line Loop(2) = {2, 126, -106, -118};
Ruled Surface(2) = {2};
Line Loop(3) = {72, -132, -108, 124};
Ruled Surface(3) = {3};
Line Loop(4) = {120, 108, -128, -4};
Ruled Surface(4) = {4};
Line Loop(5) = {69, -129, -105, 121};
Ruled Surface(5) = {5};
Line Loop(6) = {105, -125, -1, 117};
Ruled Surface(6) = {6};
Line Loop(7) = {107, 131, -71, -123};
Ruled Surface(7) = {7};
Line Loop(8) = {3, 127, -107, -119};
Ruled Surface(8) = {8};
*/

//Krümmung um Kathode A-B außen
/*				DSMC - Mod 30.09.15
Line Loop(9) = {2, 34, -3, -30};
Ruled Surface(9) = {9};
Line Loop(11) = {31, 4, -35, -3};
Ruled Surface(11) = {11};
Line Loop(13) = {1, -36, -4, 32};
Ruled Surface(13) = {13};
Line Loop(15) = {2, -33, -1, 29};
Ruled Surface(15) = {15};
*/

//Krümmung um Kathode A-B mittig
/*				DSMC - Mod 30.09.15
Line Loop(14) = {105, -160, -108, 156};
Ruled Surface(14) = {14};
Line Loop(16) = {106, -157, -105, 153};
Ruled Surface(16) = {16};
Line Loop(12) = {155, 108, -159, -107};
Ruled Surface(12) = {12};
Line Loop(10) = {106, 158, -107, -154};
Ruled Surface(10) = {10};
*/


//Ebene A Kreisringsegmente
/*				DSMC - Mod 30.09.15
Line Loop(17) = {153, 122, -81, -121};
Ruled Surface(17) = {17};
Line Loop(18) = {154, 123, -82, -122};
Ruled Surface(18) = {18};
Line Loop(19) = {123, 83, -124, -155};
Ruled Surface(19) = {19};
Line Loop(20) = {84, -121, -156, 124};
Ruled Surface(20) = {20};
Line Loop(21) = {29, 118, -153, -117};
Ruled Surface(21) = {21};
Line Loop(22) = {30, 119, -154, -118};
Ruled Surface(22) = {22};
Line Loop(23) = {155, -120, -31, 119};
Ruled Surface(23) = {23};
Line Loop(24) = {156, -117, -32, 120};
Ruled Surface(24) = {24};
*/

//Ebene B Kreisringsegmente
/*				DSMC - Mod 30.09.15
Line Loop(25) = {36, 125, -160, -128};
Ruled Surface(25) = {25};
Line Loop(26) = {125, 157, -126, -33};
Ruled Surface(26) = {26};
Line Loop(27) = {34, 127, -158, -126};
Ruled Surface(27) = {27};
Line Loop(28) = {35, 128, -159, -127};
Ruled Surface(28) = {28};
Line Loop(29) = {160, 129, -88, -132};
Ruled Surface(29) = {29};
Line Loop(30) = {132, -87, -131, 159};
Ruled Surface(30) = {30};
Line Loop(31) = {131, -86, -130, 158};
Ruled Surface(31) = {31};
Line Loop(32) = {157, 130, -85, -129};
Ruled Surface(32) = {32};
*/				
//Krümmung um Kathode A-B

/*				DSMC - Mod 30.09.15
Line Loop(33) = {82, 71, -86, -70};
Ruled Surface(33) = {33};
Line Loop(34) = {81, 70, -85, -69};
Ruled Surface(34) = {34};
Line Loop(35) = {84, 69, -88, -72};
Ruled Surface(35) = {35};
Line Loop(36) = {72, -87, -71, 83};
Ruled Surface(36) = {36};
*/

// Flanken B-C innen
/*				DSMC - Mod 30.09.15
Line Loop(37)={-138,-110,130,74};
Ruled Surface(37) = {37};
Line Loop(38) = {76, -140, -112, 132};
Ruled Surface(38) = {38};
Line Loop(39) = {109, 137, -73, -129};
Ruled Surface(39) = {39};
Line Loop(40) = {111, 139, -75, -131};
Ruled Surface(40) = {40};
*/
//Flanken B-C außen
/*				DSMC - Mod 30.09.15
Line Loop(41) = {6, 134, -110, -126};
Ruled Surface(41) = {41};
Line Loop(42) = {128, 112, -136, -8};
Ruled Surface(42) = {42};
Line Loop(43) = {109, -133, -5, 125};
Ruled Surface(43) = {43};
Line Loop(44) = {7, 135, -111, -127};
Ruled Surface(44) = {44};

//Krümmung um Kathode B-C mittig
Line Loop(45) = {110, -161, -109, 157};
Ruled Surface(45) = {45};
Line Loop(46) = {158, 111, -162, -110};
Ruled Surface(46) = {46};
Line Loop(47) = {160, 109, -164, -112};
Ruled Surface(47) = {47};
Line Loop(48) = {159, 112, -163, -111};
Ruled Surface(48) = {48};
*/
//Krümmung um Kathode B-C
/*				DSMC - Mod 30.09.15
Line Loop(49) = {92, -73, -88, 76};
Ruled Surface(49) = {49};
Line Loop(50) = {87, 76, -91, -75};
Ruled Surface(50) = {50};
Line Loop(51) = {86, 75, -90, -74};
Ruled Surface(51) = {51};
Line Loop(52) = {74, -89, -73, 85};
Ruled Surface(52) = {52};


//Krümmung um Kathode B-C außen
Line Loop(53) = {34, 7, -38, -6};
Ruled Surface(53) = {53};
Line Loop(54) = {7, 39, -8, -35};
Ruled Surface(54) = {54};
Line Loop(55) = {8, 40, -5, -36};
Ruled Surface(55) = {55};
Line Loop(56) = {5, 37, -6, -33};
Ruled Surface(56) = {56};
*/
//Ebene C Kreisringsegmente
Line Loop(57) = {91, -140, -163, 139};
Ruled Surface(57) = {57};
Line Loop(58) = {139, -90, -138, 162};
Ruled Surface(58) = {58};
Line Loop(59) = {89, -138, -161, 137};
Ruled Surface(59) = {59};
Line Loop(60) = {137, -92, -140, 164};
Ruled Surface(60) = {60};
Line Loop(61) = {39, 136, -163, -135};
Ruled Surface(61) = {61};
Line Loop(62) = {38, 135, -162, -134};
Ruled Surface(62) = {62};
Line Loop(63) = {134, -161, -133, 37};
Ruled Surface(63) = {63};
Line Loop(64) = {133, -164, -136, 40};
Ruled Surface(64) = {64};

//Krümmung um Kathode C-D
Line Loop(65) = {79, 151, -80, -91};
Ruled Surface(65) = {65};
Line Loop(66) = {80, 152, -77, -92};
Ruled Surface(66) = {66};
Line Loop(67) = {78, 150, -79, -90};
Ruled Surface(67) = {67};
Line Loop(68) = {77, 149, -78, -89};
Ruled Surface(68) = {68};

//Krümmung um Kathode c-D mittig
Line Loop(69) = {113, 165, -114, -161};
Ruled Surface(69) = {69};
Line Loop(70) = {114, 166, -115, -162};
Ruled Surface(70) = {70};
Line Loop(71) = {115, 167, -116, -163};
Ruled Surface(71) = {71};
Line Loop(72) = {113, -168, -116, 164};
Ruled Surface(72) = {72};

//Krümmung um Kathode C-D außen
Line Loop(73) = {38, 11, -42, -10};
Ruled Surface(73) = {73};
Line Loop(74) = {10, -41, -9, 37};
Ruled Surface(74) = {74};
Line Loop(75) = {9, -44, -12, 40};
Ruled Surface(75) = {75};
Line Loop(76) = {12, -43, -11, 39};
Ruled Surface(76) = {76};


//Rechtecke C-D
Line Loop(78) = {9, 141, -113, -133};
Ruled Surface(78) = {78};
Line Loop(79) = {145, -77, -137, 113};
Ruled Surface(79) = {79};
Line Loop(80) = {147, -79, -139, 115};
Ruled Surface(80) = {80};
Line Loop(81) = {11, 143, -115, -135};
Ruled Surface(81) = {81};
Line Loop(82) = {148, -80, -140, 116};
Ruled Surface(82) = {82};
Line Loop(83) = {116, -144, -12, 136};
Ruled Surface(83) = {83};
Line Loop(84) = {10, 142, -114, -134};
Ruled Surface(84) = {84};
Line Loop(85) = {114, 146, -78, -138};
Ruled Surface(85) = {85};

//Ebene D Kreisringsegmente
Line Loop(86) = {149, -146, -165, 145};
Ruled Surface(86) = {86};
Line Loop(87) = {146, 150, -147, -166};
Ruled Surface(87) = {87};
Line Loop(88) = {147, 151, -148, -167};
Ruled Surface(88) = {88};
Line Loop(89) = {145, -152, -148, 168};
Ruled Surface(89) = {89};
Line Loop(90) = {165, -142, -41, 141};
Ruled Surface(90) = {90};
Line Loop(91) = {168, -141, -44, 144};
Ruled Surface(91) = {91};
Line Loop(92) = {144, -167, -143, 43};
Ruled Surface(92) = {92};
Line Loop(93) = {143, -166, -142, 42};
Ruled Surface(93) = {93};


//KrümmungsKrams
//erster Bereich
//HKR=(((D5)*Cos((Pi*(K1X-DX1))/(2*L2))-D1)/4)+D1/2; //halber Spalt... von Kathode aus zur Krümmung
K1Laenge=(3.667/1000)/2.5;//manuell gefunden ....
K1Winkel=90-((90-KathodenSpitzenWinkel)/2);
K1Y=Tan(K1Winkel*Pi/180)*Cos(K1Winkel*Pi/180)*K1Laenge+KY;// Formel für Y, wenn die Länge über K1Laenge definiert sein soll
K1X=KDX+Cos(K1Winkel*Pi/180)*K1Laenge;
K1Z=K1Y;
Point(87)={K1X,K1Y,0,md};
Point(88)={K1X,0,K1Z,md};
Point(89)={K1X,-K1Y,0,md};
Point(90)={K1X,0,-K1Z,md};

Hilfswinkel=(ArcTan(K1Y/(K1X-DX1))*180)/Pi;
K1AX=(D5/2)*Cos((Hilfswinkel*Pi)/180);
K1AY=(D5/2)*Sin((Hilfswinkel*Pi)/180);
K1AZ=K1AY;
Point(91)={K1AX+DX1,K1AY,0,md};
Point(92)={K1AX+DX1,0,K1AZ,md};
Point(93)={K1AX+DX1,-K1AY,0,md};
Point(94)={K1AX+DX1,0,-K1AZ,md};

Point(95)={K1X,0,0,md};
Point(96)={K1AX+DX1,0,0,md};
A=Cos(2*K1Winkel*(Pi/180));
B=2*KDX*Cos(K1Winkel*Pi/180)-2*DX1*Cos(K1Winkel*Pi/180)-2*KY*Sin(K1Winkel*Pi/180);
C=KDX^2-2*KDX*DX1+DX1^2-(D5^2/4)-KY^2;
Bdurch2A=B/(2*A);
CdurchA=C/A;
NULLA1=-Bdurch2A*1000+Sqrt((Bdurch2A*1000)^2-CdurchA*1000);
NULLA1REAL=NULLA1/1000;
laenge=Sqrt((KDX+NULLA1REAL*Cos(K1Winkel*Pi/180))+(KY+NULLA1REAL*Sin(K1Winkel*Pi/180))^2);
l=Sqrt((KDX+(-(B/2*A)+Sqrt((B/2*A)^2-(C/A)))*Cos(K1Winkel*Pi/180))^2+(KY+(-(B/2*A)+Sqrt((B/2*A)^2-(C/A)))*Sin(K1Winkel*Pi/180))^2);
Printf("A = %f,B = %f,c = %f,Bdurch2A = %f,CdurchA = %f, NULLA1REAL = %f,laenge = %F",A,B,C,Bdurch2A,CdurchA,NULLA1REAL,laenge);
//Linien nach Außen "Krakenarme"
Line(173)={59,87};
Line(174)={87,91};
Line(175)={58,88};
Line(176)={88,92};
Line(177)={61,89};
Line(178)={89,93};
Line(179)={60,90};
Line(180)={90,94};

//Verbindung D-K1
Line(181)={80,87};
Line(182)={79,88};
Line(183)={82,89};
Line(184)={81,90};

//Kreise
Circle(185)={87,95,88};
Circle(186)={88,95,89};
Circle(187)={89,95,90};
Circle(188)={90,95,87};
//Außen
Circle(189)={91,96,92};
Circle(190)={92,96,93};
Circle(191)={93,96,94};
Circle(192)={94,96,91};

//zweiter Bereich
K2Laenge=(2.8725/1000)/2.5;//manuell gefunden ....
K2Winkel=KathodenSpitzenWinkel;
K2Y=Tan(K2Winkel*Pi/180)*Cos(K2Winkel*Pi/180)*K2Laenge+KSY;// Formel für Y, wenn die Länge über K1Laenge definiert sein soll
K2X=KEX+Cos(K2Winkel*Pi/180)*K2Laenge;
K2Z=K2Y;

Point(97)={K2X,0,0,md};
Point(98)={K2X,K2Y,0,md};
Point(99)={K2X,0,K2Z,md};
Point(100)={K2X,-K2Y,0,md};
Point(101)={K2X,0,-K2Z,md};

Hilfswinkel2=(ArcTan(K2Y/(K2X-DX1))*180)/Pi;

K2AX=(D5/2)*Cos((Hilfswinkel2*Pi)/180);
K2AY=(D5/2)*Sin((Hilfswinkel2*Pi)/180);
K2AZ=K2AY;
Point(102)={K2AX+DX1,0,0,md};

Point(103)={K2AX+DX1,K2AY,0,md};
Point(104)={K2AX+DX1,0,K2AZ,md};
Point(105)={K2AX+DX1,-K2AY,0,md};
Point(106)={K2AX+DX1,0,-K2AZ,md};


//Linien nach Außen "Krakenarme"
Line(193)={64,98};
Line(194)={98,103};
Line(195)={65,101};
Line(196)={101,106};
Line(197)={66,100};
Line(198)={100,105};
Line(199)={63,99};
Line(200)={99,104};

//Verbindung K1-K2

Line(201)={98,87};
Line(202)={101,90};
Line(203)={100,89};
Line(204)={99,88};


//Kreise
Circle(205)={98,97,99};
Circle(206)={99,97,100};
Circle(207)={100,97,101};
Circle(208)={101,97,98};
//Außen
Circle(209)={103,102,106};
Circle(210)={106,102,105};
Circle(211)={105,102,104};
Circle(212)={104,102,103};

//KathodenSpitze Viereck
KSY2=KSY/2;
KSZ2=KSZ/2;
Point(107)={KEX,KSY2,0,md};
Point(108)={KEX,0,KSZ2,md};
Point(109)={KEX,-KSY2,0,md};
Point(110)={KEX,0,-KSZ2,md};
//Viereck
Line(213)={107,108};
Line(214)={108,109};
Line(215)={109,110};
Line(216)={110,107};
//Linien nach außen
Line(217)={107,64};
Line(218)={108,63};
Line(219)={109,66};
Line(220)={110,65};

//DüsenHals Zusätze
DHZLenX=(K2AX+DX1)-K2X;
DHZLenY=K2AY-K2Y;
DHZLen=Sqrt(DHZLenX^2 + DHZLenY^2);
DHZWinkel=ArcTan(DHZLenX/DHZLenY);

EX3=EX2-(0.45*DHZLen*Sin(DHZWinkel));

EY3=EY2-(0.45*DHZLen*Cos(DHZWinkel));
EZ3=EY3;
Point(111)={EX3,EY3,0,md};
Point(112)={EX3,0,EZ3,md};
Point(113)={EX3,-EY3,0,md};
Point(114)={EX3,0,-EZ3,md};
Point(115)={EX3,0,0,md};
//Verbindung zum Hals
Line(221)={23,111};
Line(222)={22,112};
Line(223)={25,113};
Line(224)={24,114};
//Kreise
Circle(225)={111,115,112};
Circle(226)={112,115,113};
Circle(227)={113,115,114};
Circle(228)={114,115,111};

//Verbindung von Zusatz nach hinten -X
Line(229)={111,98};
Line(230)={112,99};
Line(231)={113,100};
Line(232)={114,101};

//Restliches Inneres
EX4=EX3;
EY4=EY3/4;
EZ4=EY4;
//Viereck
Point(116)={EX4,EY4,0,md};
Point(117)={EX4,0,EZ4,md};
Point(118)={EX4,-EY4,0,md};
Point(119)={EX4,0,-EZ4,md};

EX5=EX4;
EY5=EY3/2.5;
EZ5=EY5;
//Kreis
Point(120)={EX5,EY5,0,md};
Point(121)={EX5,0,EZ5,md};
Point(122)={EX5,-EY5,0,md};
Point(123)={EX5,0,-EZ5,md};



//Verbinung nach hinten
//>Viereck
Line(233)={116,107};
Line(234)={117,108};
Line(235)={118,109};
Line(236)={119,110};
//>Kreis
Line(237)={120,64};
Line(238)={121,63};
Line(239)={122,66};
Line(240)={123,65};
//Linien im Kreis zu dem innern Kreis
Line(241)={120,111};
Line(242)={121,112};
Line(243)={122,113};
Line(244)={123,114};

//Linien im inneren Kreis zu dem Viereck
Line(245)={120,116};
Line(246)={121,117};
Line(247)={122,118};
Line(248)={123,119};
//Bilde Viereck
Line(249)={116,117};
Line(250)={117,118};
Line(251)={118,119};
Line(252)={119,116};
//Kreis
Circle(253)={120,115,121};
Circle(254)={121,115,122};
Circle(255)={122,115,123};
Circle(256)={123,115,120};

//Restliches Inneres in Ebene F
RIFX=FX2;
//Viereck
Point(124)={RIFX,EY4,0,md};
Point(125)={RIFX,0,EZ4,md};
Point(126)={RIFX,-EY4,0,md};
Point(127)={RIFX,0,-EZ4,md};

//Innerer Kreis
Point(128)={RIFX,EY5,0,md};
Point(129)={RIFX,0,EZ5,md};
Point(130)={RIFX,-EY5,0,md};
Point(131)={RIFX,0,-EZ5,md};

//Zusatz->Kreis
Point(132)={RIFX,EY3,0,md};
Point(133)={RIFX,0,EZ3,md};
Point(134)={RIFX,-EY3,0,md};
Point(135)={RIFX,0,-EZ3,md};
Printf("EZ3=%f, EZ5=%f\n",EZ3,EZ5);
//Bilde Viereck
Line(257)={124,125};
Line(258)={125,126};
Line(259)={126,127};
Line(260)={127,124};
//Innerer Kreis
Circle(261)={128,26,129};
Circle(262)={129,26,130};
Circle(263)={130,26,131};
Circle(264)={131,26,128};
//Zusatz Kreis
Circle(265)={132,26,133};
Circle(266)={133,26,134};
Circle(267)={134,26,135};
Circle(268)={135,26,132};

//Linien nach hinten -X
//>Viereck
Line(269)={124,116};
Line(270)={125,117};
Line(271)={126,118};
Line(272)={127,119};
//>Innerer Kreis
Line(273)={128,120};
Line(274)={129,121};
Line(275)={130,122};
Line(276)={131,123};
//>Zusatz Kreis
Line(277)={132,111};
Line(278)={133,112};
Line(279)={134,113};
Line(280)={135,114};


/////////////////// Ebene G
RIGX=GX1;
GY3=GY2*(EY3/FY2);
GY4=GY2*(EY4/FY2);
GY5=GY2*(EY5/FY2);
GZ3=GY3;
GZ4=GY4;
GZ5=GY5;


//Viereck
Point(136)={RIGX,GY4,0,md};
Point(137)={RIGX,0,GZ4,md};
Point(138)={RIGX,-GY4,0,md};
Point(139)={RIGX,0,-GZ4,md};

//Innerer Kreis
Point(140)={RIGX,GY5,0,md};
Point(141)={RIGX,0,GZ5,md};
Point(142)={RIGX,-GY5,0,md};
Point(143)={RIGX,0,-GZ5,md};
//Zusatz Kreis
Point(144)={RIGX,GY3,0,md};
Point(145)={RIGX,0,GZ3,md};
Point(146)={RIGX,-GY3,0,md};
Point(147)={RIGX,0,-GZ3,md};


//Linien nach hinten -X
//>Viereck
Line(281)={124,136};
Line(282)={125,137};
Line(283)={126,138};
Line(284)={127,139};
//>Innerer Kreis
Line(285)={128,140};
Line(286)={129,141};
Line(287)={130,142};
Line(288)={131,143};
//>Zusatz Kreis
Line(289)={132,144};
Line(290)={133,145};
Line(291)={134,146};
Line(292)={135,147};

//Baue Viereck
Line(293)={136,137};
Line(294)={137,138};
Line(295)={138,139};
Line(296)={139,136};
//Kreise
//innen
Circle(297)={140,31,141};
Circle(298)={141,31,142};
Circle(299)={142,31,143};
Circle(300)={143,31,140};
//zusatz
Circle(301)={144,31,145};
Circle(302)={145,31,146};
Circle(303)={146,31,147};
Circle(304)={147,31,144};

//Linien nach außen angefangen von innen
Line(305)={136,140};
Line(306)={137,141};
Line(307)={138,142};
Line(308)={139,143};

Line(309)={140,144};
Line(310)={141,145};
Line(311)={142,146};
Line(312)={143,147};

Line(313)={144,33};
Line(314)={145,32};
Line(315)={146,35};
Line(316)={147,34};
////////////////////////// Ebene H
RIHX=HX1;
HY3=GY3;
HY4=GY4;
HY5=GY5;
HZ3=HY3;
HZ4=HY4;
HZ5=HY5;


//Viereck
Point(148)={RIHX,HY4,0,md};
Point(149)={RIHX,0,HZ4,md};
Point(150)={RIHX,-HY4,0,md};
Point(151)={RIHX,0,-HZ4,md};

//Innerer Kreis
Point(152)={RIHX,HY5,0,md};
Point(153)={RIHX,0,HZ5,md};
Point(154)={RIHX,-HY5,0,md};
Point(155)={RIHX,0,-HZ5,md};
//Zusatz Kreis
Point(156)={RIHX,HY3,0,md};
Point(157)={RIHX,0,HZ3,md};
Point(158)={RIHX,-HY3,0,md};
Point(159)={RIHX,0,-HZ3,md};


//Linien nach hinten -X
//>Viereck
Line(317)={136,148};
Line(318)={137,149};
Line(319)={138,150};
Line(320)={139,151};
//>Innerer Kreis
Line(321)={140,152};
Line(322)={141,153};
Line(323)={142,154};
Line(324)={143,155};
//>Zusatz Kreis
Line(325)={144,156};
Line(326)={145,157};
Line(327)={146,158};
Line(328)={147,159};

//Baue Viereck
Line(329)={148,149};
Line(330)={149,150};
Line(331)={150,151};
Line(332)={151,148};
//Kreise
//innen
Circle(333)={152,40,153};
Circle(334)={153,40,154};
Circle(335)={154,40,155};
Circle(336)={155,40,152};
//zusatz
Circle(337)={156,40,157};
Circle(338)={157,40,158};
Circle(339)={158,40,159};
Circle(340)={159,40,156};

//Linien nach außen angefangen von innen
Line(341)={148,152};
Line(342)={149,153};
Line(343)={150,154};
Line(344)={151,155};

Line(345)={152,156};
Line(346)={153,157};
Line(347)={154,158};
Line(348)={155,159};





//Letzer Kreis aus Ebene G
Point(160)={RIHX,GY2,0,md};
Point(161)={RIHX,0,GZ2,md};
Point(162)={RIHX,-GY2,0,md};
Point(163)={RIHX,0,-GZ2,md};

Circle(353)={160,40,161};
Circle(354)={161,40,162};
Circle(355)={162,40,163};
Circle(356)={163,40,160};
//Linien nach Hinten -X
Line(357)={33,160};
Line(358)={32,161};
Line(359)={35,162};
Line(360)={34,163};

//Krümmung D-E
Circle(45)={17,16,92};
Circle(46)={18,16,91};
Circle(47)={19,16,94};
Circle(48)={20,16,93};

Circle(361)={92,16,104};
Circle(362)={91,16,103};
Circle(363)={94,16,106};
Circle(364)={93,16,105};

Circle(365)={104,16,22};
Circle(366)={103,16,23};
Circle(367)={106,16,24};
Circle(368)={105,16,25};

///////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////

//Flächen D-E
//> Flügel innen um Kathode
Line Loop(94) = {181, -173, -170, -146};
Ruled Surface(94) = {94};
Line Loop(95) = {172, 177, -183, 148};
Ruled Surface(95) = {95};
Line Loop(96) = {182, -175, -169, -145};
Ruled Surface(96) = {96};
Line Loop(97) = {184, -179, -171, -147};
Ruled Surface(97) = {97};
//>2.Ebene Flügel um Kathode
Line Loop(98) = {182, 176, -45, 141};
Ruled Surface(98) = {98};
Line Loop(99) = {47, -180, -184, -143};
Ruled Surface(99) = {99};
Line Loop(100) = {46, -174, -181, -142};
Ruled Surface(100) = {100};
Line Loop(101) = {183, 178, -48, 144};
Ruled Surface(101) = {101};
//>Kathode
Line Loop(102) = {170, -93, -169, 149};
Ruled Surface(102) = {102};
Line Loop(103) = {169, -96, -172, 152};
Ruled Surface(103) = {103};
Line Loop(104) = {95, -172, -151, 171};
Ruled Surface(104) = {104};
Line Loop(105) = {171, -94, -170, 150};
Ruled Surface(105) = {105};
//>2.Ebene Zylinder um Kathode
Line Loop(106) = {181, 185, -182, 165};
Ruled Surface(106) = {106};
Line Loop(107) = {181, -188, -184, -166};
Ruled Surface(107) = {107};
Line Loop(108) = {184, -187, -183, -167};
Ruled Surface(108) = {108};
Line Loop(109) = {183, -186, -182, -168};
Ruled Surface(109) = {109};

//>3.Ebene Zylinder um Kathode (Außenseite)
Line Loop(110) = {46, 189, -45, 41};
Ruled Surface(110) = {110};
Line Loop(111) = {45, 190, -48, 44};
Ruled Surface(111) = {111};
Line Loop(112) = {48, 191, -47, 43};
Ruled Surface(112) = {112};
Line Loop(113) = {46, -192, -47, -42};
Ruled Surface(113) = {113};

//Verbindungsflächen in +X Richtung
Line Loop(114) = {93, 173, 185, -175};
Ruled Surface(114) = {114};
Line Loop(115) = {173, -188, -179, -94};
Ruled Surface(115) = {115};
Line Loop(116) = {179, -187, -177, -95};
Ruled Surface(116) = {116};
Line Loop(117) = {177, -186, -175, -96};
Ruled Surface(117) = {117};
//Verbindungsflächen in +X Richtung 2.Ebene
Line Loop(118) = {176, -189, -174, 185};
Ruled Surface(118) = {118};
Line Loop(119) = {176, 190, -178, -186};
Ruled Surface(119) = {119};
Line Loop(120) = {174, -192, -180, 188};
Ruled Surface(120) = {120};
Line Loop(121) = {180, -191, -178, 187};
Ruled Surface(121) = {121};

//>Flügel im KathodenSpitze
Line Loop(122) = {201, -173, 102, 193};
Ruled Surface(122) = {122};
Line Loop(123) = {104, 197, 203, -177};
Ruled Surface(123) = {123};
Line Loop(124) = {204, -175, 101, 199};
Ruled Surface(124) = {124};
Line Loop(125) = {202, -179, 103, 195};
Ruled Surface(125) = {125};
//>Flügel im KathodenSpitze 2.Ebene
Line Loop(126) = {363, -196, 202, 180};
Ruled Surface(126) = {126};
Line Loop(127) = {200, -361, -176, -204};
Ruled Surface(127) = {127};
Line Loop(128) = {362, -194, 201, 174};
Ruled Surface(128) = {128};
Line Loop(129) = {203, 178, 364, -198};
Ruled Surface(129) = {129};

////KathodenSpitze Fläche
Line Loop(130) = {102, -97, -101, 93};
Ruled Surface(130) = {130};
Line Loop(131) = {101, -100, -104, 96};
Ruled Surface(131) = {131};
Line Loop(132) = {104, -99, -103, 95};
Ruled Surface(132) = {132};
Line Loop(133) = {103, -98, -102, 94};
Ruled Surface(133) = {133};

//Verbindungsflächen in +X Richtung 2.Ebene offen in -X
Line Loop(134) = {201, -188, -202, 208};
Ruled Surface(134) = {134};
Line Loop(135) = {202, -187, -203, 207};
Ruled Surface(135) = {135};
Line Loop(136) = {203, -186, -204, 206};
Ruled Surface(136) = {136};
Line Loop(137) = {204, -185, -201, 205};
Ruled Surface(137) = {137};

//Flächen an der Spitze Kathoden
Line Loop(138) = {97, -217, 213, 218};
Ruled Surface(138) = {138};
Line Loop(139) = {218, -100, -219, -214};
Ruled Surface(139) = {139};
Line Loop(140) = {219, -99, -220, -215};
Ruled Surface(140) = {140};
Line Loop(141) = {220, -98, -217, -216};
Ruled Surface(141) = {141};
Line Loop(142) = {216, 213, 214, 215};
Ruled Surface(142) = {142};

//Verbindungsflächen in +X Richtung offen in +X
Line Loop(143) = {205, -199, 97, 193};
Ruled Surface(143) = {143};
Line Loop(144) = {199, 206, -197, 100};
Ruled Surface(144) = {144};
Line Loop(145) = {207, -195, 99, 197};
Ruled Surface(145) = {145};
Line Loop(146) = {195, 208, -193, 98};
Ruled Surface(146) = {146};
//Verbindungsflächen in +X Richtung offen in +X 2.Ebene
Line Loop(147) = {194, -212, -200, -205};
Ruled Surface(147) = {147};
Line Loop(148) = {200, -211, -198, -206};
Ruled Surface(148) = {148};
Line Loop(149) = {198, -210, -196, -207};
Ruled Surface(149) = {149};
Line Loop(150) = {196, -209, -194, -208};
Ruled Surface(150) = {150};

//Außenhülle 2. Teil
Line Loop(151) = {189, 361, 212, -362};
Ruled Surface(151) = {151};
Line Loop(152) = {361, -211, -364, -190};
Ruled Surface(152) = {152};
Line Loop(153) = {362, 209, -363, 192};
Ruled Surface(153) = {153};
Line Loop(154) = {363, 210, -364, 191};
Ruled Surface(154) = {154};

//Flügel um Verengung hin zum Düsenhals
Line Loop(155) = {238, 199, -230, -242};
Ruled Surface(155) = {155};
Line Loop(156) = {193, -229, -241, 237};
Ruled Surface(156) = {156};
Line Loop(157) = {231, -197, -239, 243};
Ruled Surface(157) = {157};
Line Loop(158) = {232, -195, -240, 244};
Ruled Surface(158) = {158};
//Außenhülle 3.Teil
Line Loop(159) = {209, 367, -50, -366};
Ruled Surface(159) = {159};
Line Loop(160) = {211, 365, -52, -368};
Ruled Surface(160) = {160};
Line Loop(161) = {368, -51, -367, 210};
Ruled Surface(161) = {161};
Line Loop(162) = {212, 366, -49, -365};
Ruled Surface(162) = {162};

//Verbindungsfläche in +X offen in +X zum Düsenhals
Line Loop(163) = {224, -227, -223, -51};
Ruled Surface(163) = {163};
Line Loop(164) = {52, 222, 226, -223};
Ruled Surface(164) = {164};
Line Loop(165) = {225, -222, 49, 221};
Ruled Surface(165) = {165};
Line Loop(166) = {224, 228, -221, 50};
Ruled Surface(166) = {166};

//Kasten von Kathodenspize aus
Line Loop(167) = {236, 216, -233, -252};
Ruled Surface(167) = {167};
Line Loop(168) = {214, -235, -250, 234};
Ruled Surface(168) = {168};
Line Loop(169) = {233, 213, -234, -249};
Ruled Surface(169) = {169};
Line Loop(170) = {236, -215, -235, 251};
Ruled Surface(170) = {170};
//Flügel auf dem Kasten
Line Loop(171) = {238, -218, -234, -246};
Ruled Surface(171) = {171};
Line Loop(172) = {219, -239, 247, 235};
Ruled Surface(172) = {172};
Line Loop(173) = {217, -237, 245, 233};
Ruled Surface(173) = {173};
Line Loop(174) = {220, -240, 248, 236};
Ruled Surface(174) = {174};

//Zylinder um dem Kasten
Line Loop(175) = {237, -97, -238, -253};
Ruled Surface(175) = {175};
Line Loop(176) = {100, -238, 254, 239};
Ruled Surface(176) = {176};
Line Loop(177) = {98, -240, 256, 237};
Ruled Surface(177) = {177};
Line Loop(178) = {99, -239, 255, 240};
Ruled Surface(178) = {178};

//Verbindungsfläche offen in -x
Line Loop(179) = {207, -232, -227, 231};
Ruled Surface(179) = {179};
Line Loop(180) = {231, -206, -230, 226};
Ruled Surface(180) = {180};
Line Loop(181) = {230, -205, -229, 225};
Ruled Surface(181) = {181};
Line Loop(182) = {229, -208, -232, 228};
Ruled Surface(182) = {182};


//Flügel um den Zusatz
Line Loop(183) = {232, 196, 367, 224};
Ruled Surface(183) = {183};
Line Loop(184) = {222, 230, 200, 365};
Ruled Surface(184) = {184};
Line Loop(185) = {194, 366, 221, 229};
Ruled Surface(185) = {185};
Line Loop(186) = {231, 198, 368, 223};
Ruled Surface(186) = {186};

//Kreisfläche zum Düsenhals
Line Loop(187) = {251, 252, 249, 250};
Ruled Surface(187) = {187};
Line Loop(188) = {256, 245, -252, -248};
Ruled Surface(188) = {188};
Line Loop(189) = {245, 249, -246, -253};
Ruled Surface(189) = {189};
Line Loop(190) = {246, 250, -247, -254};
Ruled Surface(190) = {190};
Line Loop(191) = {251, -248, -255, 247};
Ruled Surface(191) = {191};

Line Loop(192) = {225, -242, -253, 241};
Ruled Surface(192) = {192};
Line Loop(193) = {242, 226, -243, -254};
Ruled Surface(193) = {193};
Line Loop(194) = {243, 227, -244, -255};
Ruled Surface(194) = {194};
Line Loop(195) = {244, 228, -241, -256};
Ruled Surface(195) = {195};

//Verbindungslinien Ebene F
Line(369) = {124,128};
Line(370) = {125,129};
Line(371) = {126,130};
Line(372) = {127,131};

Line(373) = {128,132};
Line(374) = {129,133};
Line(375) = {130,134};
Line(376) = {131,135};

Line(377) = {132,28};
Line(378) = {133,27};
Line(379) = {134,30};
Line(380) = {135,29};

//Flächen um Kasten in Düsenhals
Line Loop(196) = {257, 270, -249, -269};
Ruled Surface(196) = {196};
Line Loop(197) = {269, -252, -272, 260};
Ruled Surface(197) = {197};
Line Loop(198) = {270, 250, -271, -258};
Ruled Surface(198) = {198};
Line Loop(199) = {251, -272, -259, 271};
Ruled Surface(199) = {199};

//Flügel am Kasten
Line Loop(200) = {273, 245, -269, 369};
Ruled Surface(200) = {200};
Line Loop(201) = {246, -270, 370, 274};
Ruled Surface(201) = {201};
Line Loop(202) = {271, -247, -275, -371};
Ruled Surface(202) = {202};
Line Loop(203) = {372, 276, 248, -272};
Ruled Surface(203) = {203};
//Zylinder um Kasten
Line Loop(204) = {256, -273, -264, 276};
Ruled Surface(204) = {204};
Line Loop(205) = {261, 274, -253, -273};
Ruled Surface(205) = {205};
Line Loop(206) = {274, 254, -275, -262};
Ruled Surface(206) = {206};
Line Loop(207) = {275, 255, -276, -263};
Ruled Surface(207) = {207};
//Flügel um Zylinder
Line Loop(208) = {241, -277, -373, 273};
Ruled Surface(208) = {208};
Line Loop(209) = {376, 280, -244, -276};
Ruled Surface(209) = {209};
Line Loop(210) = {375, 279, -243, -275};
Ruled Surface(210) = {210};
Line Loop(211) = {242, -278, -374, 274};
Ruled Surface(211) = {211};
//2. Zylinder

Line Loop(212) = {277, 225, -278, -265};
Ruled Surface(212) = {212};
Line Loop(213) = {277, -228, -280, 268};
Ruled Surface(213) = {213};
Line Loop(214) = {227, -280, -267, 279};
Ruled Surface(214) = {214};
Line Loop(215) = {266, 279, -226, -278};
Ruled Surface(215) = {215};
//Flügel mit Schräge um 2. Zylinder
Line Loop(216) = {221, -277, 377, -14};
Ruled Surface(216) = {216};
Line Loop(217) = {379, -16, 223, -279};
Ruled Surface(217) = {217};
Line Loop(218) = {222, -278, 378, -13};
Ruled Surface(218) = {218};
Line Loop(219) = {-224, 15, -380, 280};
Ruled Surface(219) = {219};
//Düsenhals Hülle
Line Loop(220) = {13, 53, -14, -49};
Ruled Surface(220) = {220};
Line Loop(221) = {13, -56, -16, 52};
Ruled Surface(221) = {221};
Line Loop(222) = {51, 16, -55, -15};
Ruled Surface(222) = {222};
Line Loop(223) = {15, -54, -14, 50};
Ruled Surface(223) = {223};

//Ebene F Kreismitte
Line Loop(224) = {260, 257, 258, 259};
Ruled Surface(224) = {224};
Line Loop(225) = {257, 370, -261, -369};
Ruled Surface(225) = {225};
Line Loop(226) = {262, -371, -258, 370};
Ruled Surface(226) = {226};
Line Loop(227) = {371, 263, -372, -259};
Ruled Surface(227) = {227};
Line Loop(228) = {372, 264, -369, -260};
Ruled Surface(228) = {228};

//Ebene F Kreisrand
Line Loop(229) = {265, -374, -261, 373};
Ruled Surface(229) = {229};
Line Loop(230) = {373, -268, -376, 264};
Ruled Surface(230) = {230};
Line Loop(231) = {376, -267, -375, 263};
Ruled Surface(231) = {231};
Line Loop(232) = {375, -266, -374, 262};
Ruled Surface(232) = {232};
Line Loop(233) = {266, 379, 56, -378};
Ruled Surface(233) = {233};
Line Loop(234) = {267, 380, 55, -379};
Ruled Surface(234) = {234};
Line Loop(235) = {380, -54, -377, -268};
Ruled Surface(235) = {235};
Line Loop(236) = {377, -53, -378, -265};
Ruled Surface(236) = {236};

//Kasten F-G
Line Loop(237) = {257, 282, -293, -281};
Ruled Surface(237) = {237};
Line Loop(238) = {294, -283, -258, 282};
Ruled Surface(238) = {238};
Line Loop(239) = {-283, 259, 284, -295};
Ruled Surface(239) = {239};
Line Loop(240) = {296, -281, -260, 284};
Ruled Surface(240) = {240};
//Flügel um Kasten
Line Loop(241) = {-285, -369, 281, 305};
Ruled Surface(241) = {241};
Line Loop(242) = {308, -288, -372, 284};
Ruled Surface(242) = {242};
Line Loop(243) = {370, 286, -306, -282};
Ruled Surface(243) = {243};
Line Loop(244) = {307, -287, -371, 283};
Ruled Surface(244) = {244};
//Zylinder um Kasten
Line Loop(245) = {298, -287, -262, 286};
Ruled Surface(245) = {245};
Line Loop(246) = {261, 286, -297, -285};
Ruled Surface(246) = {246};
Line Loop(247) = {300, -285, -264, 288};
Ruled Surface(247) = {247};
Line Loop(248) = {263, 288, -299, -287};
Ruled Surface(248) = {248};
//Flügel um Zylinder
Line Loop(249) = {373, 289, -309, -285};
Ruled Surface(249) = {249};
Line Loop(250) = {374, 290, -310, -286};
Ruled Surface(250) = {250};
Line Loop(251) = {311, -291, -375, 287};
Ruled Surface(251) = {251};
Line Loop(252) = {312, -292, -376, 288};
Ruled Surface(252) = {252};
//Flügel auf 2. Zylinder
Line Loop(253) = {315, -20, -379, 291};
Ruled Surface(253) = {253};
Line Loop(254) = {18, -313, -289, 377};
Ruled Surface(254) = {254};
Line Loop(255) = {314, -17, -378, 290};
Ruled Surface(255) = {255};
Line Loop(256) = {380, 19, -316, -292};
Ruled Surface(256) = {256};
//3. Zylinder
Line Loop(257) = {304, -289, -268, 292};
Ruled Surface(257) = {257};
Line Loop(258) = {265, 290, -301, -289};
Ruled Surface(258) = {258};
Line Loop(259) = {-290, 266, 291, -302};
Ruled Surface(259) = {259};
Line Loop(260) = {-291, 267, 292, -303};
Ruled Surface(260) = {260};
//Hülle
Line Loop(261) = {53, 18, -57, -17};
Ruled Surface(261) = {261};
Line Loop(262) = {18, 58, -19, -54};
Ruled Surface(262) = {262};
Line Loop(263) = {59, -20, -55, 19};
Ruled Surface(263) = {263};
Line Loop(264) = {56, 17, -60, -20};
Ruled Surface(264) = {264};
//Ebene G Kreisfläche
Line Loop(265) = {296, 293, 294, 295};
Ruled Surface(265) = {265};
Line Loop(266) = {300, -305, -296, 308};
Ruled Surface(266) = {266};
Line Loop(267) = {293, 306, -297, -305};
Ruled Surface(267) = {267};
Line Loop(268) = {306, 298, -307, -294};
Ruled Surface(268) = {268};
Line Loop(269) = {307, 299, -308, -295};
Ruled Surface(269) = {269};
Line Loop(270) = {309, -304, -312, 300};
Ruled Surface(270) = {270};
Line Loop(271) = {309, 301, -310, -297};
Ruled Surface(271) = {271};
Line Loop(272) = {302, -311, -298, 310};
Ruled Surface(272) = {272};
Line Loop(273) = {311, 303, -312, -299};
Ruled Surface(273) = {273};
Line Loop(274) = {57, -313, 301, 314};
Ruled Surface(274) = {274};
Line Loop(275) = {314, -60, -315, -302};
Ruled Surface(275) = {275};
Line Loop(278) = {315, -59, -316, -303};
Ruled Surface(278) = {278};
Line Loop(279) = {316, -58, -313, -304};
Ruled Surface(279) = {279};

//Kasten H-H
Line Loop(280) = {296, 317, -332, -320};
Ruled Surface(280) = {280};
Line Loop(281) = {329, -318, -293, 317};
Ruled Surface(281) = {281};
Line Loop(282) = {294, 319, -330, -318};
Ruled Surface(282) = {282};
Line Loop(283) = {331, -320, -295, 319};
Ruled Surface(283) = {283};
//Flügel auf Kasten
Line Loop(284) = {305, 321, -341, -317};
Ruled Surface(284) = {284};
Line Loop(285) = {344, -324, -308, 320};
Ruled Surface(285) = {285};
Line Loop(286) = {322, -342, -318, 306};
Ruled Surface(286) = {286};
Line Loop(287) = {307, 323, -343, -319};
Ruled Surface(287) = {287};
//Zylinder um Kasten
Line Loop(288) = {298, 323, -334, -322};
Ruled Surface(288) = {288};
Line Loop(289) = {300, 321, -336, -324};
Ruled Surface(289) = {289};
Line Loop(290) = {333, -322, -297, 321};
Ruled Surface(290) = {290};
Line Loop(291) = {299, 324, -335, -323};
Ruled Surface(291) = {291};
//Flügel auf Zylinder
Line Loop(292) = {309, 325, -345, -321};
Ruled Surface(292) = {292};
Line Loop(293) = {348, -328, -312, 324};
Ruled Surface(293) = {293};
Line Loop(294) = {346, -326, -310, 322};
Ruled Surface(294) = {294};
Line Loop(295) = {311, 327, -347, -323};
Ruled Surface(295) = {295};
//2. Zylinder
Line Loop(296) = {337, -326, -301, 325};
Ruled Surface(296) = {296};
Line Loop(297) = {304, 325, -340, -328};
Ruled Surface(297) = {297};
Line Loop(298) = {339, -328, -303, 327};
Ruled Surface(298) = {298};
Line Loop(299) = {302, 327, -338, -326};
Ruled Surface(299) = {299};

Line(381) = {160,42};
Line(382) = {161,41};
Line(383) = {162,44};
Line(384) = {163,43};


Line(349)={156,160};
Line(350)={157,161};
Line(351)={158,162};
Line(352)={159,163};

//2. Flügel
Line Loop(300) = {349, -357, -313, 325};
Ruled Surface(300) = {300};
Line Loop(301) = {316, 360, -352, -328};
Ruled Surface(301) = {301};
Line Loop(302) = {351, -359, -315, 327};
Ruled Surface(302) = {302};
Line Loop(303) = {314, 358, -350, -326};
Ruled Surface(303) = {303};
//3.Zylinder
Line Loop(304) = {-58,357,-356,-360};
Ruled Surface(304) = {304};
Line Loop(305) = {57, 357, 353, -358};
Ruled Surface(305)={305};
Line Loop(306) = {-358, -60, 359, -354};
Ruled Surface(306)={306};
Line Loop(307) = {59, 359, 355, -360};
Ruled Surface(307)={307};
//3.Flügel
Line Loop(308) = {23, 27, -384, -360};
Ruled Surface(308) = {308};
Line Loop(309) = {381, -26, -22, 357};
Ruled Surface(309) = {309};
Line Loop(310) = {21, 25, -382, -358};
Ruled Surface(310) = {310};
Line Loop(311) = {383, -28, -24, 359};
Ruled Surface(311) = {311};
//Ebene G große Aussenfläche
Line Loop(312) = {22, 62, -23, -58};
Ruled Surface(312) = {312};
Line Loop(313) = {22, -61, -21, 57};
Ruled Surface(313) = {313};
Line Loop(314) = {21, -64, -24, 60};
Ruled Surface(314) = {314};
Line Loop(315) = {63, -24, -59, 23};
Ruled Surface(315) = {315};
//Ebene H alles
Line Loop(316) = {329, 330, 331, 332};
Ruled Surface(316) = {316};
Line Loop(317) = {341, 333, -342, -329};
Ruled Surface(317) = {317};
Line Loop(318) = {330, 343, -334, -342};
Ruled Surface(318) = {318};
Line Loop(319) = {332, 341, -336, -344};
Ruled Surface(319) = {319};
Line Loop(320) = {344, -335, -343, 331};
Ruled Surface(320) = {320};
Line Loop(321) = {340, -345, -336, 348};
Ruled Surface(321) = {321};
Line Loop(322) = {348, -339, -347, 335};
Ruled Surface(322) = {322};
Line Loop(323) = {347, -338, -346, 334};
Ruled Surface(323) = {323};
Line Loop(324) = {346, -337, -345, 333};
Ruled Surface(324) = {324};
Line Loop(325) = {353, -350, -337, 349};
Ruled Surface(325) = {325};
Line Loop(326) = {350, 354, -351, -338};
Ruled Surface(326) = {326};
Line Loop(327) = {351, 355, -352, -339};
Ruled Surface(327) = {327};
Line Loop(328) = {352, 356, -349, -340};
Ruled Surface(328) = {328};
Line Loop(329) = {381, 66, -384, 356};
Ruled Surface(329) = {329};
Line Loop(330) = {65, -381, 353, 382};
Ruled Surface(330) = {330};
Line Loop(331) = {354, 383, 68, -382};
Ruled Surface(331) = {331};
Line Loop(332) = {383, -67, -384, -355};
Ruled Surface(332) = {332};
//Hülle
Line Loop(333) = {26, -65, -25, 61};
Ruled Surface(333) = {333};
Line Loop(334) = {64, 25, -68, -28};
Ruled Surface(334) = {334};
Line Loop(335) = {28, -67, -27, 63};
Ruled Surface(335) = {335};
Line Loop(336) = {27, -66, -26, 62};
Ruled Surface(336) = {336};
/*
//Volumen A-B		DSMC-Mod -30.09.15

Surface Loop(500) = {19, 36, 7, 30, 3, 12};
Volume(500) = {500};
Surface Loop(501) = {20, 35, 29, 14, 5, 3};
Volume(501) = {501};
Surface Loop(502) = {18, 33, 31, 1, 10, 7};
Volume(502) = {502};
Surface Loop(503) = {17, 34, 32, 1, 5, 16};
Volume(503) = {503};
Surface Loop(504) = {22, 9, 27, 10, 2, 8};
Volume(504) = {504};
Surface Loop(505) = {8, 23, 11, 28, 12, 4};
Volume(505) = {505};
Surface Loop(506) = {4, 24, 13, 25, 6, 14};
Volume(506) = {506};
Surface Loop(507) = {21, 15, 26, 6, 16, 2};
Volume(507) = {507};
//Volumen B-C

Surface Loop(508) = {32, 52, 59, 45, 37, 39};
Volume(508) = {508};
Surface Loop(509) = {29, 49, 60, 38, 39, 47};
Volume(509) = {509};
Surface Loop(510) = {50, 30, 57, 40, 38, 48};
Volume(510) = {510};
Surface Loop(511) = {51, 31, 58, 40, 46, 37};
Volume(511) = {511};
Surface Loop(512) = {56, 63, 26, 45, 43, 41};
Volume(512) = {512};
Surface Loop(513) = {41, 53, 27, 62, 46, 44};
Volume(513) = {513};
Surface Loop(514) = {44, 54, 61, 28, 48, 42};
Volume(514) = {514};
Surface Loop(515) = {42, 25, 55, 64, 43, 47};
Volume(515) = {515};
*/
//Volumen C-D
Surface Loop(516) = {65, 88, 57, 71, 80, 82};
Volume(516) = {516};
Surface Loop(517) = {80, 67, 87, 58, 85, 70};
Volume(517) = {517};
Surface Loop(518) = {85, 68, 86, 59, 79, 69};
Volume(518) = {518};
Surface Loop(519) = {66, 89, 60, 82, 72, 79};
Volume(519) = {519};
Surface Loop(520) = {76, 92, 61, 81, 83, 71};
Volume(520) = {520};
Surface Loop(521) = {73, 62, 93, 70, 81, 84};
Volume(521) = {521};
Surface Loop(522) = {75, 91, 64, 72, 78, 83};
Volume(522) = {522};
Surface Loop(523) = {74, 90, 63, 84, 78, 69};
Volume(523) = {523};

//
Surface Loop(524) = {102, 86, 106, 94, 96, 114};
Volume(524) = {524};
Surface Loop(525) = {115, 105, 87, 94, 97, 107};
Volume(525) = {525};
Surface Loop(526) = {116, 104, 88, 108, 95, 97};
Volume(526) = {526};
Surface Loop(527) = {109, 96, 117, 89, 103, 95};
Volume(527) = {527};

//
Surface Loop(528) = {110, 90, 106, 118, 100, 98};
Volume(528) = {528};
Surface Loop(529) = {99, 101, 92, 112, 108, 121};
Volume(529) = {529};
Surface Loop(530) = {111, 91, 98, 109, 101, 119};
Volume(530) = {530};
Surface Loop(531) = {93, 113, 107, 100, 99, 120};
Volume(531) = {531};
//
Surface Loop(532) = {122, 124, 114, 137, 143, 130};
Volume(532) = {532};
Surface Loop(533) = {125, 133, 115, 134, 146, 122};
Volume(533) = {533};
Surface Loop(534) = {123, 131, 124, 117, 144, 136};
Volume(534) = {534};
Surface Loop(535) = {123, 132, 116, 125, 135, 145};
Volume(535) = {535};
//
Surface Loop(536) = {120, 153, 128, 126, 134, 150};
Volume(536) = {536};
Surface Loop(537) = {121, 154, 126, 135, 129, 149};
Volume(537) = {537};
Surface Loop(538) = {119, 152, 129, 136, 148, 127};
Volume(538) = {538};
Surface Loop(539) = {118, 151, 128, 137, 147, 127};
Volume(539) = {539};

//Volumen von Kathoden-Spitze aus
Surface Loop(540) = {142, 167, 169, 170, 168, 187};
Volume(540) = {540};
Surface Loop(541) = {141, 177, 174, 167, 173, 188};
Volume(541) = {541};
Surface Loop(542) = {173, 175, 169, 171, 189, 138};
Volume(542) = {542};
Surface Loop(543) = {172, 176, 171, 168, 139, 190};
Volume(543) = {543};
Surface Loop(544) = {178, 174, 170, 172, 140, 191};
Volume(544) = {544};
//1.Ebene um Volumen von Spitze aus
Surface Loop(545) = {175, 156, 181, 155, 192, 143};
Volume(545) = {545};
Surface Loop(546) = {157, 178, 145, 158, 179, 194};
Volume(546) = {546};
Surface Loop(547) = {146, 182, 156, 177, 195, 158};
Volume(547) = {547};
Surface Loop(548) = {144, 176, 155, 193, 180, 157};
Volume(548) = {548};
//2. Ebene um Volumen von Spitze aus
Surface Loop(549) = {185, 147, 162, 165, 184, 181};
Volume(549) = {549};
Surface Loop(550) = {150, 159, 166, 185, 183, 182};
Volume(550) = {550};
Surface Loop(551) = {148, 160, 164, 186, 180, 184};
Volume(551) = {551};
Surface Loop(552) = {183, 149, 161, 163, 186, 179};
Volume(552) = {552};
//Düsenhals Kasten und 1. Zylinder
Surface Loop(553) = {197, 196, 198, 199, 187, 224};
Volume(553) = {553};
Surface Loop(554) = {197, 204, 203, 200, 228, 188};
Volume(554) = {554};
Surface Loop(555) = {191, 202, 207, 203, 199, 227};
Volume(555) = {555};
Surface Loop(556) = {226, 206, 198, 201, 202, 190};
Volume(556) = {556};
Surface Loop(557) = {189, 196, 205, 200, 201, 225};
Volume(557) = {557};
//2.Zylinder
Surface Loop(558) = {204, 213, 195, 208, 209, 230};
Volume(558) = {558};
Surface Loop(559) = {208, 211, 205, 212, 192, 229};
Volume(559) = {559};
Surface Loop(560) = {206, 215, 193, 210, 211, 232};
Volume(560) = {560};
Surface Loop(561) = {214, 194, 207, 209, 210, 231};
Volume(561) = {561};
//3.Zylinder außen
Surface Loop(562) = {165, 220, 236, 216, 218, 212};
Volume(562) = {562};
Surface Loop(563) = {218, 215, 164, 221, 233, 217};
Volume(563) = {563};
Surface Loop(564) = {234, 222, 163, 214, 217, 219};
Volume(564) = {564};
Surface Loop(565) = {219, 235, 223, 166, 213, 216};
Volume(565) = {565};
//F-G Kasten+1.Zylinder
Surface Loop(566) = {224, 265, 237, 240, 238, 239};
Volume(566) = {566};
Surface Loop(567) = {239, 248, 242, 244, 227, 269};
Volume(567) = {567};
Surface Loop(568) = {245, 238, 244, 243, 226, 268};
Volume(568) = {568};
Surface Loop(569) = {246, 237, 241, 243, 267, 225};
Volume(569) = {569};
Surface Loop(570) = {240, 241, 247, 242, 266, 228};
Volume(570) = {570};
//2.Zylinder
Surface Loop(571) = {249, 257, 247, 252, 270, 230};
Volume(571) = {571};
Surface Loop(572) = {249, 258, 250, 246, 271, 229};
Volume(572) = {572};
Surface Loop(573) = {245, 250, 251, 259, 272, 232};
Volume(573) = {573};
Surface Loop(574) = {231, 248, 251, 260, 252, 273};
Volume(574) = {574};
//3.Zylinder
Surface Loop(575) = {258, 261, 236, 274, 254, 255};
Volume(575) = {575};
Surface Loop(576) = {262, 279, 235, 254, 257, 256};
Volume(576) = {576};
Surface Loop(577) = {256, 260, 263, 278, 234, 253};
Volume(577) = {577};
Surface Loop(578) = {253, 264, 233, 275, 255, 259};
Volume(578) = {578};
//Austrittsvolumen
Surface Loop(579) = {265, 280, 283, 282, 281, 316};
Volume(579) = {579};
Surface Loop(580) = {317, 284, 286, 290, 281, 267};
Volume(580) = {580};
Surface Loop(581) = {266, 284, 289, 285, 280, 319};
Volume(581) = {581};
Surface Loop(582) = {286, 282, 287, 288, 318, 268};
Volume(582) = {582};
Surface Loop(583) = {283, 291, 285, 287, 320, 269};
Volume(583) = {583};
Surface Loop(584) = {290, 292, 296, 294, 324, 271};
Volume(584) = {584};
Surface Loop(585) = {270, 289, 297, 292, 293, 321};
Volume(585) = {585};
Surface Loop(586) = {293, 298, 291, 295, 273, 322};
Volume(586) = {586};
Surface Loop(587) = {323, 294, 295, 299, 288, 272};
Volume(587) = {587};
Surface Loop(588) = {305, 296, 303, 300, 325, 274};
Volume(588) = {588};
Surface Loop(589) = {300, 279, 301, 297, 304, 328};
Volume(589) = {589};
Surface Loop(590) = {327, 298, 307, 301, 302, 278};
Volume(590) = {590};
Surface Loop(591) = {326, 302, 303, 299, 306, 275};
Volume(591) = {591};
Surface Loop(592) = {305, 309, 310, 333, 330, 313};
Volume(592) = {592};
Surface Loop(593) = {336, 329, 312, 309, 308, 304};
Volume(593) = {593};
Surface Loop(594) = {308, 335, 332, 315, 307, 311};
Volume(594) = {594};
Surface Loop(595) = {311, 306, 310, 334, 314, 331};
Volume(595) = {595};

/*
//Bereich A-B in YZ-Richtung Gerade
Transfinite Line {121,122,123,124,117,118,119,120} = 11 Using Progression 1;
//Bereich A-B in YZ-Richtung Krümmung
Transfinite Line {153,154,155,156,84,81,82,83,29,30,31,32} = 15 Using Progression 1;
//Bereich A-B in X-Richtung
Transfinite Line {2,106,70,3,107,71,72,108,4,1,105,69} = 21 Using Progression 1;

//Bereich B-C in YZ-Richtung Gerade
Transfinite Line {129,125,130,126,131,127,128,132} = 11 Using Progression 1;
//Bereich B-C in YZ-Richtung Krümmung
Transfinite Line {36,160,33,157,158,134,159,35,88,85,86,87} = 15 Using Progression 1;
//Bereich B-C in X-Richtung
Transfinite Line {8,112,76,73,109,5,74,110,6,7,111,75} = 21 Using Progression 1;
*/

//Bereich C-D in YZ-Richtung Gerade
Transfinite Line {136,140,137,133,138,134,135,139} = 11 Using Progression 1;
//Bereich C-D in YZ-Richtung Krümmung
Transfinite Line {89,92,91,90,161,162,163,164,40,37,38,39} = 15 Using Progression 1;
//Bereich C-D in X-Richtung
Transfinite Line {77,113,9,79,115,11,78,114,10,80,116,12} = 11 Using Progression 1;


//Bereich D-K1 in YZ-Richtung Gerade
Transfinite Line {146,142,147,143,145,141,148,144} = 11 Using Progression 1;
//Bereich D-K1 in YZ-Richtung Krümmung
Transfinite Line {168,152,44,150,42,166,149,165,41,151,167,43} = 15 Using Progression 1;
//Bereich D-K1 in X-Richtung
Transfinite Line {172,183,48,170,181,46,171,184,47,182,169,45} = 21 Using Progression 1;

//Bereich K1-K2 in YZ-Richtung Gerade
Transfinite Line {179,180,175,176,173,174,177,178} = 11 Using Progression 1;
//Bereich K1-K2 in YZ-Richtung Krümmung
Transfinite Line {188,185,186,187,96,93,95,94,189,192,191,190} = 15 Using Progression 1;
//Bereich K1-K2 in X-Richtung
Transfinite Line {362,201,102,104,203,364,101,204,361,363,202,103} = 21 Using Progression 1;


//Bereich K2-E in YZ-Richtung Gerade
Transfinite Line {193,194,197,198,199,200,196,195,217,219,218,220} = 11 Using Progression 1;
//Bereich K2-E in YZ-Richtung Krümmung
Transfinite Line {205,208,206,207,212,209,211,210,99,100,98,97,213,216,214,215} = 15 Using Progression 1;
//Bereich K2-E in X-Richtung
Transfinite Line {366,229,237,233,235,239,231,368,365,230,238,234,236,240,232,367} = 21 Using Progression 1;


//Bereich E-F in YZ-Richtung Gerade
Transfinite Line {245,241,221,247,243,223,246,242,222,248,244,224} = 11 Using Progression 1;
//Bereich E-F in YZ-Richtung Krümmung
Transfinite Line {249,253,225,49,252,256,228,50,250,254,226,52,251,255,227,51} = 15 Using Progression 1;
//Bereich E-F in X-Richtung
Transfinite Line {13,14,15,16,269,270,271,272,273,274,275,276,277,278,279,280} = 21 Using Progression 1;

//Bereich F-G in YZ-Richtung Gerade
Transfinite Line {370,378,372,380,369,377,371,379} = 11 Using Progression 1;
Transfinite Line {374,376,375,373} = 11 Using Progression 1;
//Bereich F-G in YZ-Richtung Krümmung
Transfinite Line {258,262,266,56,257,261,265,53,54,268,264,260,259,263,267,55} = 15 Using Progression 1;
//Bereich F-G in X-Richtung
Transfinite Line {282,281,284,283,286,285,288,287,290,289,292,291,20,18,19,17} = 50 Using Progression 1.05;

//Bereich G-H in YZ-Richtung Gerade
Transfinite Line {306,314,305,313,308,316,315,307} = 11 Using Progression 1;
Transfinite Line {309,310,311,312} = 11 Using Progression 1;
//Bereich G-H in YZ-Richtung Austrittsvolumen-Gerade
Transfinite Line {21,23,22,24} = 21 Using Progression 1.15;	//Mod
//Bereich G-H in YZ-Richtung Krümmung
Transfinite Line {294,298,302,60,293,297,301,57,296,300,304,58,295,299,303,59} = 15 Using Progression 1;
//Bereich G-H in YZ-Richtung Austrittsvolumen-Krümmung
Transfinite Line {61,62,63,64} = 15 Using Progression 1;
//Bereich G-H in X-Richtung
Transfinite Line {317,321,325,26,318,322,326,25,319,323,327,28,320,324,328,27,357,358,359,360} = 50 Using Progression 1.015;//Mod


//Ebene H in YZ-Richtung Gerade
Transfinite Line {341,345,349,342,346,350,343,347,351,344,348,352} = 11 Using Progression 1;
//Ebene H in YZ-Richtung Austrittsvolumen-Gerade
Transfinite Line {381,384,382,383} = 21 Using Progression 1.15;	//Mod
//Ebene H in YZ-Richtung Krümmung
Transfinite Line {329,333,337,353,356,340,336,332,331,335,339,355,330,334,338,354} = 15 Using Progression 1;
//Ebene H in YZ-Richtung Austrittsvolumen-Krümmung
Transfinite Line {65,66,67,68} = 15 Using Progression 1;

For i In {57:76}
Transfinite Surface {i} = {};
EndFor

For i In {78:275}
Transfinite Surface {i} = {};
EndFor

For i In {278:336}
Transfinite Surface {i} = {};
EndFor



For j In {57:76}
Recombine Surface {j};
EndFor

For j In {78:275}
Recombine Surface {j};
EndFor

For j In {278:336}
Recombine Surface {j};
EndFor


For k In {516:595}
Transfinite Volume{k}={};
EndFor

Physical Surface("inlet") = {57,58,59,60,61,62,63,64};//done
//Physical Surface("fixedWalls") = {142,141,140,139,138,102,103,104,105,65,66,67,68,49,50,51,52,33,34,35,36,9,11,13,15,53,54,55,56,73,74,75,76,159,160,161,162,220,221,222,223,261,262,263,264,312,313,314,315};//done
Physical Surface("fixedWalls") = {142,141,140,139,138,102,103,104,105,65,66,67,68,73,74,75,76,159,160,161,162,220,221,222,223,261,262,263,264,312,313,314,315};//done
Physical Surface("Kathode") = {130,131,132,133};//done 
Physical Surface("Anode") = {110,111,112,113,151,152,153,154};//done 
Physical Surface("seiten") = {333,334,335,336}; //done
Physical Surface("outlet") = {329,330,331,332,325,326,327,328,321,322,323,324,317,318,319,320,316}; //done


Physical Volume("internalField") = {516,517,518,519,520,521,522,523,524,525,526,527,528,529,530,531,532,533,534,535,536,537,538,539,540,541,542,543,544,545,546,547,548,549,550,551,552,553,554,555,556,557,558,559,560,561,562,563,564,565,566,567,568,569,570,571,572,573,574,575,576,577,578,579,580,581,582,583,584,585,586,587,588,589,590,591,592,593,594,595};//done


