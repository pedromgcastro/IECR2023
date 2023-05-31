$Title Global optimization of QCPs arising from the design of water using networks

$Ontext
          Code for of the global optimization algorithm featured in the article submitted to Ind. Eng. Chem. Res. in June 2023.
          "Global Optimization of QCPs using MIP relaxations with a base-2 logarithmic partitioning scheme"

          WHATM=    7         BARON
                    8         GUROBI
                    9         Mixed-radix MDT algorithm with spatial B&B
         
          WHATSBB=  1         LP-based bound contracting in spatial B&B
                    2         MIP-based bound contracting in spatial B&B
$Offtext

$eolcom #

SETS
I         Operation units /O1*O30/
K         Fresh Water Sources /W1*W7/
C         Process Contaminants          /C1*C6/
J         Numbers   /0*29/
KD        Significant digits   /KD1*KD21/
N         Domain partitions   /N1*N1000/ #N1048576
IT        Iterations        /IT0*IT100/
SP        Solutions in pool /file1*file100/
ND        Nodes for spatial B&B         /ND1*ND20000/
ACTOP(I)  Active operations
ACTWS(K)  Active water sources
ACTCN(C)  Active contaminants
SPOOL(SP) Actual solutions in pool
ACTJ(J,KD)          Active numbers J in digit KD
ACTK(KD,I,C)        Active significant digits KD for outlet concentration of treatment I contaminant C
FML(I)    Operations with fixed mass load
FFO(I)    Operations with fixed flowrate
BNDFTS(I) Flowrate from unit I to treatment system being bound contracted
BRCHFTS(I)          Flowrate from unit I to treatment system being branched
BNDFOO(I,I)         Flowrate from unit I to I being bound contracted
BRCHFOO(I,I)        Flowrate from unit I to I being branched
BNDCOUT(I,C)        Concentration from unit I contaminant C being bound contracted
CONNECT(I,I)        Connected units
OUTLETC(I,C)        Outlet concentrations to consider
FIXCOUT(I,C)        Range of outlet concentration from unit I contaminant C has been reduced to zero by OBBT
FIXFOO(I,I)         Range of flow from unit I to unit IL has been reduced to zero by OBBT
FIXFTS(I) Range of flow from unit I to to treatment system has been reduced to zero by OBBT
ACTICJK(I,C,J,KD)   For operation I concentration C value J in digit KD is active
CURRENT(ND)         Current node
SRCNODE(ND)         Nodes search in spatial branch and bound
NEWNODE(ND)         New node
WAITING(ND)         Waiting node list
NODEBRN(ND,ND)      Relation between parent and child nodes
NEXTFOO(ND,I,I)     Next FOO variable for branching
NEXTFTS(ND,I)       Next FTS variable for branching
;

ALIAS(I,IL,ILL,ILLL);ALIAS(KD,KDL,KDLL);ALIAS(ND,NDL);

SCALAR    WHATM     What Model                    /9/;
SCALAR    WHATSBB   What OBBT in SBB              /2/;
SCALAR    TCPU      Total computational time of algorithm (s)         /3600/;
SCALAR    CPUOBBT   Approximate maximum computational time in OBBT (s)  /300/;
SCALAR    CPURELX   Time limit in relaxation for trigerring SBB (s)   /100/;
SCALAR    TRGGAP    Optimality gap (%)     /0.0001/;
SCALAR    NDIGITS   Number of digits used in numeric representation
SCALAR    NWATS     Number of water sources
SCALAR    PARTVAR   Number of partitioned variables
SCALAR    DISGVAR   Number of disaggregated variables
SCALAR    BILINTR   Number of bilinear variables
SCALAR    LBOUND    Lower Bound from relaxation problem
SCALAR    UBOUND    Upper bound from original problem
SCALAR    OPTGAP    Optimality gap (%)
SCALAR    CURIT     Current iteration /1/;
SCALAR    OBBTCPU   Computational time of OBBT (s)
SCALAR    MIPCPU    Computational time in Global McCormick or NMDT lower bounding
SCALAR    TSTART    Time at the start of a certain procedure /0/;
SCALAR    NBOUNDP   Number of bounding problems
SCALAR    ESTOBBT   Estimate of OBBT computational time (s)
SCALAR    ESTRELX   Estimate of relaxation computational time (s)
SCALAR    SBBOBBT   Number of digits used for performing OBBT in spatial branch and bound /0/;
SCALAR    SBBRELX   Number of digits used for performing relaxation in spatial branch and bound
SCALAR    CARDSP    Number of solutions in pool
SCALAR    CLOOP     Continue loop /1/;
SCALAR    DONE      Termination of spatial B&B    /0/;
SCALAR    BISECT    Bisection of variable domain
SCALAR    FIRST     Controller for loop
SCALAR    FIRST2    Controller for loop
SCALAR    STATUS    Another controller;

PARAMETER
BASIS(KD) Basis used in digit KD        /KD1*KD1 2/
CINMAX(I,C)         Maximum Inlet Concentration in Operation I of Contaminant C (ppm)
COUTMAX(I,C)        Maximum Outlet Concentration in Operation I of Contaminant C (ppm)
PFIN(I)   Inlet flowrate to fixed flowrate operation I (in th-1)
PFOUT(I)  Outlet flowrate to fixed flowrate operation I (in th-1)
FLIM(I)   Limiting Flowrate Water in Operation I (in th-1)
FMIN(I)   Minimum flowrate in operation I (th-1)
FMAX(I)   Maximum flowrate in operation I (th-1)
COFW(K,C) Concentration of Contaminant C on FreshWater K
MLFC(I,C) Mass to remove for operation I contaminant C
CINMIN(I,C)         Lower bound on inlet concentration operation I contaminant C
COUTMIN(I,C)        Lower bound on outlet concentration operation I contaminant C
FOOMIN(I,IL)        Minimum flowrate from fixed load operation I to operation IL
FOOMAX(I,IL)        Maximum flowrate from fixed load operation I to operation IL
FOOMINT(I,IL,C,J,KD)          Minimum flowrate from fixed load operation I to operation IL when constrained to values of contaminant C number J digit KD
FOOMAXT(I,IL,C,J,KD)          Maximum flowrate from fixed load operation I to operation IL when constrained to values of contaminant C number J digit KD
FTSMIN(I) Minimum flowrate from fixed load unit I to treatment system
FTSMAX(I) Maximum flowrate from fixed load unit I to treatment system
FTSMINT(I,C,J,KD)   Minimum flowrate from fixed load unit I to treatment system when constrained to values of contaminant C number J digit KD
FTSMAXT(I,C,J,KD)   Maximum flowrate from fixed load unit I to treatment system when constrained to values of contaminant C number J digit KD
COUTJL(I,C,J)       Upper bound of outlet concentration from unit I contaminant C value J in first digit
COUTJU(I,C,J)       Lower bound of outlet concentration from unit I contaminant C value J in first digit
RPART(I,C)          Real partitions for the outlet concentration from unit I contaminant C
PRICE(K)  Price of water source K
CARDK(I,C)          Required number of positions for outlet concentration of unit I contaminant C
KVAL(KD,I,C)        Position of significant digit KD in numerical basis representation for concentration of unit I contaminant C
JVAL(J,KD)          Value to use for number J in digit KD
LBFTS(I,ND)         Lower bound for FTS variable from unit I in node ND
UBFTS(I,ND)         Upper bound for FTS variable from unit I in node ND
LBFOO(I,IL,ND)      Lower bound for FOO variable from unit I to IL in node ND
UBFOO(I,IL,ND)      Upper bound for FOO variable from unit I to IL in node ND
LBCOUT(I,C,ND)      Lower bound for COUT variable form unit I contaminant C node ND
UBCOUT(I,C,ND)      Upper bound for COUT variable form unit I contaminant C node ND
COUTRNG(I,C,IT)     Range of COUT for operation I contaminant C in iteration IT of OBBT
COUTRR(I,C,IT)      Range reduction of COUT variables for operation I contaminant C in iteration IT of OBBT (%)
FOORNG(I,IL,IT)     Range of FOO from unit I to IL in iteration IT of OBBT
FOORR(I,IL,IT)      Range reduction of FOO variables from unit I to IL in iteration IT of OBBT (%)
FTSRNG(I,IT)        Range of FTS from unit I in iteration IT of OBBT
FTSRR(I,IT)  Range reduction of FTS variables from unit I in iteration IT of OBBT (%)
NFIXV(IT) Number of fixed variables after OBBT in iteration IT
AVERR(IT) Average range reduction in iteration IT
LBSP(SP)  Lower bound from solution SP in pool
UBSP(SP)  Upper bound from solution SP in pool
ERRFOO(IL,I,C)      Error of bilinear term associated to FOO from unit IL to I contaminant C
ERRFTS(I,C)         Error of bilinear term associated to FTS from unit I contaminant C
BOUND(ND) Lower bound of node ND
LOGBB(ND,*)         Logging information for branch and bound tree
LOGIT(ND,*)         Logging information for iterations
LOGND(ND,I,I,*)     Logging information for nodes FOO variables
LOGND2(ND,I,*)      Logging information for nodes FTS variables
LOGND3(ND,I,C,*)    Logging information for nodes COUT variables
;

*Change location of Data files
$include  C:\Users\Castro\Documents\GAMS files\Aguas\Teles\Exemplo9.txt

ACTOP(I)=yes$(SMAX(C,COUTMAX(I,C)) GT 0);
ACTWS(K)=yes$(ord(K) LE NWATS);
ACTCN(C)=yes$(SMAX(I$(ACTOP(I)),COUTMAX(I,C)) GT 0);
FML(I)=yes$(FLIM(I) GT 0);
FFO(I)=yes$(PFIN(I) GT 0 or PFOUT(I) GT 0);
MLFC(I,C) = FLIM(I)*(COUTMAX(I,C)-CINMAX(I,C));
CINMIN(I,C)$(FML(I)) = SMIN(K$(ACTWS(K)), COFW(K,C));
COUTMIN(I,C)$(FML(I)) = SMIN(K$(ACTWS(K)), COFW(K,C))+MLFC(I,C)/FLIM(I);
FOOMIN(FML,ACTOP)=0;
FOOMAX(FML,IL)$(ACTOP(IL))=MIN(FLIM(FML), FLIM(IL)$(FML(IL))+PFIN(IL)$(FFO(IL)));
FTSMIN(FML)=0;
FTSMAX(FML)=FLIM(FML);
FMIN(FML)=SUM(ACTOP,FOOMIN(FML,ACTOP))+FTSMIN(FML);
FMAX(FML)=MIN(FLIM(FML),SUM(ACTOP,FOOMAX(FML,ACTOP))+FTSMAX(FML));
CONNECT(FML,ACTOP)=yes;OUTLETC(FML,ACTCN)=yes;
BNDFOO(CONNECT)=no;BNDFTS(FML)=no;BNDCOUT(OUTLETC)=no;
CURRENT('ND1')=yes;NEWNODE('ND1')=yes;WAITING('ND1')=yes;
LBCOUT(OUTLETC,CURRENT)=COUTMIN(OUTLETC);UBCOUT(OUTLETC,CURRENT)=COUTMAX(OUTLETC);
LBFOO(CONNECT,CURRENT)=FOOMIN(CONNECT);UBFOO(CONNECT,CURRENT)=FOOMAX(CONNECT);
LBFTS(FML,CURRENT)=FTSMIN(FML);UBFTS(FML,CURRENT)=FTSMAX(FML);
COUTRNG(OUTLETC,'IT0')=SUM(CURRENT,UBCOUT(OUTLETC,CURRENT)-LBCOUT(OUTLETC,CURRENT));
FOORNG(CONNECT,'IT0')=SUM(CURRENT,UBFOO(CONNECT,CURRENT)-LBFOO(CONNECT,CURRENT));
FTSRNG(FML,'IT0')=SUM(CURRENT,UBFTS(FML,CURRENT)-LBFTS(FML,CURRENT));
PARTVAR=card(OUTLETC);DISGVAR=card(CONNECT)+card(FML);BILINTR=card(CONNECT)*card(ACTCN)+card(OUTLETC);
FIXCOUT(OUTLETC)=no;FIXFOO(CONNECT)=no;FIXFTS(FML)=no;
NDIGITS=SUM(KD$(BASIS(KD) GT 0),1);
CARDK(I,C)$(FML(I) and ACTCN(C))=NDIGITS;
ACTK(KD,I,C)$(FML(I) and ACTCN(C))=yes$(ord(KD) LE CARDK(I,C));
KVAL(KD,I,C)$(ACTK(KD,I,C))=-BASIS(KD);
ACTJ(J,KD)$(ord(KD) LE NDIGITS)=yes$(ord(J) LE BASIS(KD));
JVAL(J,KD)$(ACTJ(J,KD))=ord(J)-1;

Display ACTOP,ACTWS,ACTCN,FML,FFO,MLFC,FOOMAX,COUTMIN,COUTMAX;
Display PARTVAR,DISGVAR,BILINTR,NDIGITS,CARDK,ACTK,KVAL,JVAL;

BINARY VARIABLES
Y(I,C,J,KD)         Binary variable for operation I contaminant C with digit J is associated to position KD

POSITIVE VARIABLES
FFW(K,I)  Flowrate From Fresh Water K to Operation I
FIO(I)    Flowrate Inlet Operation I
FTS(I)    Flowrate From Operation I to Treatment Systems
FOO(I,IL) Flowrate from Operation I to operation IL
CIN(I,C)  Inlet Concentration of Contaminant C to Operation I
COUT(I,C) Outlet Concentration of Contaminant C to Operation I
UMIn(I,C) Mass into operation unit I for contaminant C
UMOut(I,C)          Mass from operation unit I for contaminant C
ZFOO(I,IL,C)        Variable associated to bilinear term involving FOO I IL and COUT I C
ZFTS(I,C)           Variable associated to bilinear term involving FTS I and COUT I C
FOODS(I,IL,C,J,KD)  Disaggregated Flowrate from operation unit I to unit IL associated to contaminant C digit J in position KD
FTSDS(I,C,J,KD)     Disaggregated Flowrate from operation unit I to treatment system associated to contaminant C digit J in position KD
LCOUT(I,C)          Discretized variable related to concentration from unit I contaminant C
LFOOC(IL,I,C)       Discretized variable related to mass from unit IL into I contaminant C
LFTSC(I,C)          Discretized variable related to mass from unit I contaminant C into treatment system
LFOOCR(I,IL,C)      Residual variable from unit I to IL contaminant C
LFTSCR(I,C)         Residual variable from unit I contaminant C into treatment system
LCOUTR(I,C)         Residual variable related to concentration from unit I contaminant C

VARIABLES
OBJ       Variable of the objective function;

EQUATIONS
OBJ_NLP   Objective function for NLP (minimization of total flow through treatment units)
NLP_1(I)  Flow balance over mixer at unit I
NLP_2(I,C)          Mass balance over mixer in unit I contaminant C for fixed load operations
NLP_3(I,C)          Maximum inlet concentration constraint for fixed flowrate I for contaminant C
NLP_4(I,C)          Mass balance in inlet and outlet operation unit I for contaminant C
NLP_5(I)  Flow balance over splitter at operation unit I

LP_1(I,C) Maximum inlet concentration constraint for unit I for contaminant C
LP_2(I,C) Mass balance over mixer in operation unit I for contaminant C
LP_3(I,C) Mass balance over splitter in operation unit I for contaminant C
LP_4(I,C) Mass balance in inlet and outlet operation unit I for contaminant C

GLMC_1(I,IL,C)      Mc Cormick underestimators 1 for contaminant mass from I into IL contaminant C
GLMC_2(I,IL,C)      Mc Cormick underestimators 2 for contaminant mass from I into IL contaminant C
GLMC_3(I,IL,C)      Mc Cormick overestimators 1 for contaminant mass from I into IL contaminant C
GLMC_4(I,IL,C)      Mc Cormick overestimators 2 for contaminant mass from I into IL contaminant C
GLMC_5(I,C)         Mc Cormick underestimators 1 for contaminant mass from I into treatment system contaminant C
GLMC_6(I,C)         Mc Cormick underestimators 2 for contaminant mass from I into treatment system contaminant C
GLMC_7(I,C)         Mc Cormick overestimators 1 for contaminant mass from I into treatment system contaminant C
GLMC_8(I,C)         Mc Cormick overestimators 2 for contaminant mass from I into treatment system contaminant C

NMDT1(I,C)          Relation between outlet concentration unit I contaminant C and discretized variable
NMDT2(I,IL,C)       Definition of bilinear term involving mass from unit I to IL contaminant C
NMDT3(I,C)          Definition of bilinear term involving mass from unit I contaminant C into treatment system
MDT_1(I,IL,C)       Definition of bilinear term involving mass from unit I to IL contaminant C
MDT_2(I,C)          Definition of bilinear term involving mass from unit I to treatment system contaminant C
MDT_3(I,C)          Definition of outlet concentration from unit I to contaminant C as a multiparametric sum
MDT_4(I,IL,C,KD)    Outlet flow between I and IL must be equal to disaggregated flow over C for every KD
MDT_5(I,C,KD)       Flow from unit I to treatment system must be equal to disaggregated flow over C for every KD
MDT_6(I,C,KD,J)     Relation between disaggregated flowrate variables and binary variables (upper bound)
MDT_7(I,C,KD,J)     Relation between disaggregated flowrate variables and binary variables (lower bound)
MDT_8(I,C,KD)       For a given unit I contaminant C and significant digit KD a single value is selected
MDT_9(IL,I,C)       McCormick relaxation of bilinear term involving mass from unit I to IL contaminant C and its residual
MDT_10(IL,I,C)      McCormick relaxation of bilinear term involving mass from unit I to IL contaminant C and its residual
MDT_11(IL,I,C)      McCormick relaxation of bilinear term involving mass from unit I to IL contaminant C and its residual
MDT_12(IL,I,C)      McCormick relaxation of bilinear term involving mass from unit I to IL contaminant C and its residual
MDT_13(I,C)         McCormick relaxation of bilinear term involving mass from unit I to contaminant C and its residual
MDT_14(I,C)         McCormick relaxation of bilinear term involving mass from unit I to contaminant C and its residual
MDT_15(I,C)         McCormick relaxation of bilinear term involving mass from unit I to contaminant C and its residual
MDT_16(I,C)         McCormick relaxation of bilinear term involving mass from unit I to contaminant C and its residual

OBJ_BTN   Objective function for optimality based bound tightening
BTNEQ1    Objective function of original NLP must be lower than current best solution;

OBJ_NLP.. OBJ=e=Sum((K,I)$(ACTOP(I) and ACTWS(K)), FFW(K,I)*PRICE(K));
NLP_1(I)$(ACTOP(I))..         FIO(i)$(FML(i))+PFIN(i)$(FFO(i))=e=Sum(k$(ACTWS(k)), FFW(k,i)) + Sum(IL$(ACTOP(IL)), FOO(IL,I))  ;
NLP_2(I,C)$(FML(I) and ACTCN(C))..      FIO(i)*CIN(i,c)=e=Sum(k$(ACTWS(k)), FFW(k,i)*COFW(k,c)) + Sum(IL$(ACTOP(IL)), FOO(IL,I)*(COUT(IL,c)$(FML(IL))+COUTMAX(IL,c)$(FFO(IL))))  ;
NLP_3(I,C)$(FFO(I) and ACTCN(C))..      PFIN(i)*CINMAX(i,c)=g=Sum(k$(ACTWS(k)), FFW(k, i)*COFW(k,c)) + Sum(IL$(ACTOP(IL)), FOO(IL,I)*(COUT(IL,c)$(FML(IL))+COUTMAX(IL,c)$(FFO(IL))))  ;
NLP_4(I,C)$(FML(I) and ACTCN(C))..      MLFC(i,c)=e=FIO(i)*(COUT(i,c) - CIN(i,c))  ;
NLP_5(I)$(ACTOP(I))..         FIO(i)$(FML(i))+PFOUT(i)$(FFO(i))=e=FTS(i) + Sum(IL$(ACTOP(IL)), FOO(I,IL));

LP_1(I,C)$(ACTOP(I) and ACTCN(C))..     (PFIN(I)$(FFO(I))+FIO(I)$(FML(I)))*CINMAX(I,C)=g=Sum(k$(ACTWS(k)),FFW(k,I)*COFW(k,C)) + Sum(IL$(ACTOP(IL)),ZFOO(IL,I,C)$(FML(IL))+FOO(IL,I)*COUTMAX(IL,C)$(FFO(IL)));
LP_2(I,C)$(FML(I) and ACTCN(C))..       UMIn(I,C)=e=Sum(k$(ACTWS(k)),FFW(k,I)*COFW(k,C)) + Sum(IL$(ACTOP(IL)),ZFOO(IL,I,C)$(FML(IL))+FOO(IL,I)*COUTMAX(IL,C)$(FFO(IL)));
LP_3(I,C)$(FML(I) and ACTCN(C))..       UMOut(I,C)=e=Sum(IL$(ACTOP(IL)),ZFOO(I,IL,C))+ZFTS(I,C);
LP_4(I,C)$(FML(I) and ACTCN(C))..       MLFC(I,C)=e=UMOut(I,C) - UMIn(I,C)   ;

GLMC_1(IL,I,C)$(FML(IL) and ACTOP(I) and ACTCN(C))..        ZFOO(IL,I,C)=g=FOO(IL,I)*COUTMIN(IL,C)+FOOMIN(IL,I)*COUT(IL,C)-FOOMIN(IL,I)*COUTMIN(IL,C);
GLMC_2(IL,I,C)$(FML(IL) and ACTOP(I) and ACTCN(C))..        ZFOO(IL,I,C)=g=FOO(IL,I)*COUTMAX(IL,C)+FOOMAX(IL,I)*COUT(IL,C)-FOOMAX(IL,I)*COUTMAX(IL,C);
GLMC_3(IL,I,C)$(FML(IL) and ACTOP(I) and ACTCN(C))..        ZFOO(IL,I,C)=l=FOO(IL,I)*COUTMIN(IL,C)+FOOMAX(IL,I)*COUT(IL,C)-FOOMAX(IL,I)*COUTMIN(IL,C);
GLMC_4(IL,I,C)$(FML(IL) and ACTOP(I) and ACTCN(C))..        ZFOO(IL,I,C)=l=FOO(IL,I)*COUTMAX(IL,C)+FOOMIN(IL,I)*COUT(IL,C)-FOOMIN(IL,I)*COUTMAX(IL,C);
GLMC_5(I,C)$(FML(I) and ACTCN(C))..     ZFTS(I,C)=g=FTS(I)*COUTMIN(I,C)+FTSMIN(I)*COUT(I,C)-FTSMIN(I)*COUTMIN(I,C);
GLMC_6(I,C)$(FML(I) and ACTCN(C))..     ZFTS(I,C)=g=FTS(I)*COUTMAX(I,C)+FTSMAX(I)*COUT(I,C)-FTSMAX(I)*COUTMAX(I,C);
GLMC_7(I,C)$(FML(I) and ACTCN(C))..     ZFTS(I,C)=l=FTS(I)*COUTMIN(I,C)+FTSMAX(I)*COUT(I,C)-FTSMAX(I)*COUTMIN(I,C);
GLMC_8(I,C)$(FML(I) and ACTCN(C))..     ZFTS(I,C)=l=FTS(I)*COUTMAX(I,C)+FTSMIN(I)*COUT(I,C)-FTSMIN(I)*COUTMAX(I,C);

NMDT1(I,C)$(FML(I) and ACTCN(C))..      COUT(I,C)=e=COUTMIN(I,C)+LCOUT(I,C)*(COUTMAX(I,C)-COUTMIN(I,C));
NMDT2(IL,I,C)$(FML(IL) and ACTOP(I) and ACTCN(C))..         ZFOO(IL,I,C)=e=FOO(IL,I)*COUTMIN(IL,C)+LFOOC(IL,I,C)*(COUTMAX(IL,C)-COUTMIN(IL,C));
NMDT3(I,C)$(FML(I) and ACTCN(C))..      ZFTS(I,C)=e=FTS(I)*COUTMIN(I,C)+LFTSC(I,C)*(COUTMAX(I,C)-COUTMIN(I,C));
MDT_1(IL,I,C)$(FML(IL) and ACTOP(I) and ACTCN(C))..         LFOOC(IL,I,C)=e=SUM((J,KD)$(ACTJ(J,KD) and ACTK(KD,IL,C)),FOODS(IL,I,C,J,KD)*JVAL(J,KD)*PROD(KDL$(ord(KDL) LE ord(KD)),1/BASIS(KDL)))+LFOOCR(IL,I,C);
MDT_2(I,C)$(FML(I) and ACTCN(C))..      LFTSC(I,C)=e=SUM((J,KD)$(ACTJ(J,KD) and ACTK(KD,I,C)),FTSDS(I,C,J,KD)*JVAL(J,KD)*PROD(KDL$(ord(KDL) LE ord(KD)),1/BASIS(KDL)))+LFTSCR(I,C);
MDT_3(I,C)$(FML(I) and ACTCN(C))..      LCOUT(I,C)=e=SUM((J,KD)$(ACTJ(J,KD) and ACTK(KD,I,C)),JVAL(J,KD)*Y(I,C,J,KD)*PROD(KDL$(ord(KDL) LE ord(KD)),1/BASIS(KDL)))+LCOUTR(I,C);
MDT_4(I,IL,C,KD)$(FML(I) and ACTOP(IL) and ACTCN(C) and ACTK(KD,I,C))..         FOO(I,IL)=e=SUM(J$(ACTJ(J,KD)),FOODS(I,IL,C,J,KD));
MDT_5(I,C,KD)$(ACTK(KD,I,C))..          FTS(I)=e=SUM(J$(ACTJ(J,KD)),FTSDS(I,C,J,KD));
MDT_6(I,C,KD,J)$(ACTJ(J,KD) and ACTK(KD,I,C))..   SUM(IL$(ACTOP(IL)),FOODS(I,IL,C,J,KD))+FTSDS(I,C,J,KD)=l=FMAX(I)*Y(I,C,J,KD);
MDT_7(I,C,KD,J)$(ACTJ(J,KD) and ACTK(KD,I,C))..   SUM(IL$(ACTOP(IL)),FOODS(I,IL,C,J,KD))+FTSDS(I,C,J,KD)=g=FMIN(I)*Y(I,C,J,KD);
MDT_8(I,C,KD)$(ACTK(KD,I,C))..          SUM(J$(ACTJ(J,KD)),Y(I,C,J,KD))=e=1;
MDT_9(IL,I,C)$(FML(IL) and ACTOP(I) and ACTCN(C))..         LFOOCR(IL,I,C)=l=FOOMAX(IL,I)*LCOUTR(IL,C);
MDT_10(IL,I,C)$(FML(IL) and ACTOP(I) and ACTCN(C))..        LFOOCR(IL,I,C)=g=FOOMIN(IL,I)*LCOUTR(IL,C);
MDT_11(IL,I,C)$(FML(IL) and ACTOP(I) and ACTCN(C))..        LFOOCR(IL,I,C)=g=(FOO(IL,I)-FOOMAX(IL,I))*MIN(1,PROD(KD$(ACTK(KD,IL,C)),1/BASIS(KD)))+FOOMAX(IL,I)*LCOUTR(IL,C);
MDT_12(IL,I,C)$(FML(IL) and ACTOP(I) and ACTCN(C))..        LFOOCR(IL,I,C)=l=(FOO(IL,I)-FOOMIN(IL,I))*MIN(1,PROD(KD$(ACTK(KD,IL,C)),1/BASIS(KD)))+FOOMIN(IL,I)*LCOUTR(IL,C);
MDT_13(I,C)$(FML(I) and ACTCN(C))..     LFTSCR(I,C)=l=FTSMAX(I)*LCOUTR(I,C);
MDT_14(I,C)$(FML(I) and ACTCN(C))..     LFTSCR(I,C)=g=FTSMIN(I)*LCOUTR(I,C);
MDT_15(I,C)$(FML(I) and ACTCN(C))..     LFTSCR(I,C)=g=(FTS(I)-FTSMAX(I))*MIN(1,PROD(KD$(ACTK(KD,I,C)),1/BASIS(KD)))+FTSMAX(I)*LCOUTR(I,C);
MDT_16(I,C)$(FML(I) and ACTCN(C))..     LFTSCR(I,C)=l=(FTS(I)-FTSMIN(I))*MIN(1,PROD(KD$(ACTK(KD,I,C)),1/BASIS(KD)))+FTSMIN(I)*LCOUTR(I,C);

OBJ_BTN.. OBJ=e=SUM(I$(BNDFTS(I)),FTS(I))+SUM((I,IL)$(BNDFOO(I,IL)),FOO(I,IL))+SUM((I,C)$(BNDCOUT(I,C)),COUT(I,C));
BTNEQ1..  SUM((K,I)$(ACTOP(I) and ACTWS(K)), FFW(K,I)*PRICE(K))=l=UBOUND;


CIN.up(I,C)$(ACTOP(I) and ACTCN(C)) = CINMAX(I,C);
CIN.lo(I,C)$(ACTOP(I) and ACTCN(C)) = CINMIN(I,C);
COUT.up(I,C)$(ACTOP(I) and ACTCN(C)) = COUTMAX(I,C);
COUT.lo(I,C)$(ACTOP(I) and ACTCN(C)) = COUTMIN(I,C);
FOO.up(I,IL)$(FML(I) and ACTOP(IL)) = FOOMAX(I,IL);
FTS.up(i)$(FML(I)) = FTSMAX(I);
FIO.up(I)$(FML(I)) = FLIM(I);
FFW.up(K,I)$(ACTWS(K) and FML(I)) = 10000;
LCOUT.up(I,C)$(FML(I) and ACTCN(C))=1;
LCOUTR.up(I,C)$(FML(I) and ACTCN(C))=MIN(1,PROD(KD$(ACTK(KD,I,C)),1/BASIS(KD)));

OPTION optcr=1E-6;
OPTION limrow=0;
OPTION limcol=0;
OPTION Solprint=Off;
OPTION iterlim=2000000000;
OPTION decimals=6;
OPTION threads=0;
OPTION NLP=CONOPT;
OPTION QCP=GloMIQO;
OPTION MIP=CPLEX;
OPTION reslim=3600;

MODEL WUNNLP using /OBJ_NLP,NLP_1,NLP_2,NLP_3,NLP_4,NLP_5/;
MODEL WUNLP using /OBJ_NLP,NLP_1,NLP_5,LP_1,LP_2,LP_3,LP_4,GLMC_1,GLMC_2,GLMC_3,GLMC_4,GLMC_5,GLMC_6,GLMC_7,GLMC_8/;
MODEL MDT /OBJ_NLP,NLP_1,NLP_5,LP_1,LP_2,LP_3,LP_4,NMDT1,NMDT2,NMDT3,MDT_1,MDT_2,MDT_3,MDT_4,MDT_5,MDT_6,MDT_7,MDT_8,MDT_9,MDT_10,MDT_11,MDT_12,MDT_13,MDT_14,MDT_15,MDT_16/;
MODEL OBBT_LP /OBJ_BTN,BTNEQ1,NLP_1,NLP_5,LP_1,LP_2,LP_3,LP_4,GLMC_1,GLMC_2,GLMC_3,GLMC_4,GLMC_5,GLMC_6,GLMC_7,GLMC_8/;
MODEL OBBT_MDT /OBJ_BTN,BTNEQ1,NLP_1,NLP_5,LP_1,LP_2,LP_3,LP_4,NMDT1,NMDT2,NMDT3,MDT_1,MDT_2,MDT_3,MDT_4,MDT_5,MDT_6,MDT_7,MDT_8,MDT_9,MDT_10,MDT_11,MDT_12,MDT_13,MDT_14,MDT_15,MDT_16/;

*Change location of files
file OptFile /C:\Users\Castro\Documents\GAMS files\GitHub\WUN\Cplex.op9/; #Directory must match location of WUN.gms file
file Results   /C:\Users\Castro\Documents\My Data Sources\GAMS Output\OBBT.txt/;
file Search   /C:\Users\Castro\Documents\My Data Sources\GAMS Output\GlobalSearchWUN.txt/;
file SPoolRes /C:\Users\Castro\Documents\My Data Sources\GAMS Output\SPoolRes.txt/;
file fsoln;
Results.pw=700;Search.pw=700;SPoolRes.pw=700;
Results.pc=6;Search.pc=6;SPoolRes.pc=6; #Separa com tabs os valores das variaveis


if(WHATM EQ 7,
OPTION NLP=BARON;
SOLVE WUNNLP using NLP minimizing OBJ;
LBOUND=WUNNLP.objest;
UBOUND=WUNNLP.objval;
OPTGAP=ABS((UBOUND-LBOUND)/LBOUND)*100;
Display OBJ.l,FFW.l,FIO.l,FOO.l,FTS.l,CIN.l,COUT.l,UBOUND,LBOUND,OPTGAP;
);
if(WHATM EQ 8,
OPTION QCP=GUROBI;
WUNNLP.optfile=2; #option file has single line: nonconvex 2
SOLVE WUNNLP using QCP minimizing OBJ;
LBOUND=WUNNLP.objest;
UBOUND=WUNNLP.objval;
OPTGAP=ABS((UBOUND-LBOUND)/LBOUND)*100;
Display OBJ.l,FFW.l,FIO.l,FOO.l,FTS.l,CIN.l,COUT.l,UBOUND,LBOUND,OPTGAP;
);

if(WHATM EQ 9,
MDT.optcr=1E-7;
OBBT_MDT.reslim=60;
OBBT_MDT.optcr=1E-5;
put Search;
put 'CPU','LB','UB','Gap (%)','Method','Digits'/;

#Quick computation of lower and upper bound
SOLVE WUNLP using LP minimizing OBJ;
MIPCPU=MAX(WUNLP.resusd,0.001);
CIN.l(I,C)$(ACTOP(I) and FML(I) and ACTCN(C) and FIO.l(I) NE 0) = (Sum(k$(ACTWS(k)), FFW.l(k,i)*COFW(k,c)) + Sum(IL$(ACTOP(IL)), FOO.l(IL,i)*(COUT.l(IL,c)$(FML(IL))+COUTMAX(IL,c)$(FFO(IL))))) / FIO.l(i)  ;
LBOUND=WUNLP.objval;
SOLVE WUNNLP using NLP minimizing OBJ;
    if(WUNNLP.modelstat EQ 2,
    UBOUND=WUNNLP.objval;
    else
    UBOUND=+INF;
    );
OPTGAP=ABS((UBOUND-LBOUND)/LBOUND)*100;
put TimeElapsed,LBOUND:9:5,UBOUND:9:5,OPTGAP:9:5,'LP relaxation' /;
putclose Search;Search.ap=1; #Appends text to file instead of overwriting
NBOUNDP=2*(PARTVAR+DISGVAR);
ESTOBBT=MIPCPU*NBOUNDP;TSTART=TimeElapsed;
CLOOP$(OPTGAP LE TRGGAP)=0;
#Global McCormick OBBT
    loop(IT$(CLOOP and ESTOBBT LT CPUOBBT and ord(IT) EQ CURIT+1),
        loop(OUTLETC$(not FIXCOUT(OUTLETC)),
        BNDCOUT(OUTLETC)=yes;
        SOLVE OBBT_LP using LP minimizing OBJ;
        LBCOUT(OUTLETC,CURRENT)=MIN(OBJ.l,COUTMAX(OUTLETC)); #Avoids numerical errors
        COUTMIN(OUTLETC)=SUM(CURRENT,LBCOUT(OUTLETC,CURRENT));COUT.lo(OUTLETC)=COUTMIN(OUTLETC);
        SOLVE OBBT_LP using LP maximizing OBJ;
        UBCOUT(OUTLETC,CURRENT)=MAX(OBJ.l,COUTMIN(OUTLETC)); #Avoids numerical errors
        COUTMAX(OUTLETC)=SUM(CURRENT,UBCOUT(OUTLETC,CURRENT));COUT.up(OUTLETC)=COUTMAX(OUTLETC);
        BNDCOUT(OUTLETC)=no;
        );
        loop(CONNECT$(not FIXFOO(CONNECT)),
        BNDFOO(CONNECT)=yes;
        SOLVE OBBT_LP using LP minimizing OBJ;
        LBFOO(CONNECT,CURRENT)=MIN(OBJ.l,FOOMAX(CONNECT)); #Avoids numerical errors
        FOOMIN(CONNECT)=SUM(CURRENT,LBFOO(CONNECT,CURRENT));FOO.lo(CONNECT)=FOOMIN(CONNECT);
        FMIN(FML)=SUM(ACTOP,FOOMIN(FML,ACTOP))+FTSMIN(FML);
        SOLVE OBBT_LP using LP maximizing OBJ;
        UBFOO(CONNECT,CURRENT)=MAX(OBJ.l,FOOMIN(CONNECT)); #Avoids numerical errors
        FOOMAX(CONNECT)=SUM(CURRENT,UBFOO(CONNECT,CURRENT));FOO.up(CONNECT)=FOOMAX(CONNECT);
        FMAX(FML)=MIN(FLIM(FML),SUM(ACTOP,FOOMAX(FML,ACTOP))+FTSMAX(FML));
        BNDFOO(CONNECT)=no;
        );
        loop(FML$(not FIXFTS(FML)),
        BNDFTS(FML)=yes;
        SOLVE OBBT_LP using LP minimizing OBJ;
        LBFTS(FML,CURRENT)=MIN(OBJ.l,FTSMAX(FML)); #Avoids numerical errors
        FTSMIN(FML)=SUM(CURRENT,LBFTS(FML,CURRENT));FTS.lo(FML)=FTSMIN(FML);
        FMIN(FML)=SUM(ACTOP,FOOMIN(FML,ACTOP))+FTSMIN(FML);
        SOLVE OBBT_LP using LP maximizing OBJ;
        UBFTS(FML,CURRENT)=MAX(OBJ.l,FTSMIN(FML)); #Avoids numerical errors
        FTSMAX(FML)=SUM(CURRENT,UBFTS(FML,CURRENT));FTS.up(FML)=FTSMAX(FML);
        FMAX(FML)=MIN(FLIM(FML),SUM(ACTOP,FOOMAX(FML,ACTOP))+FTSMAX(FML));
        BNDFTS(FML)=no;
        );
    #Report results of OBBT
    OBBTCPU=TimeElapsed-TSTART;CPUOBBT=MAX(0,CPUOBBT-OBBTCPU);
    COUTRNG(OUTLETC,IT)=SUM(CURRENT,UBCOUT(OUTLETC,CURRENT)-LBCOUT(OUTLETC,CURRENT));
    FOORNG(CONNECT,IT)=SUM(CURRENT,UBFOO(CONNECT,CURRENT)-LBFOO(CONNECT,CURRENT));
    FTSRNG(FML,IT)=SUM(CURRENT,UBFTS(FML,CURRENT)-LBFTS(FML,CURRENT));
    COUTRR(OUTLETC,IT)=100$(COUTRNG(OUTLETC,'IT0') EQ 0)+((COUTRNG(OUTLETC,'IT0')-COUTRNG(OUTLETC,IT))/COUTRNG(OUTLETC,'IT0')*100)$(COUTRNG(OUTLETC,'IT0') GT 0);
    FOORR(CONNECT,IT)=(FOORNG(CONNECT,'IT0')-FOORNG(CONNECT,IT))/FOORNG(CONNECT,'IT0')*100;
    FTSRR(FML,IT)=(FTSRNG(FML,'IT0')-FTSRNG(FML,IT))/FTSRNG(FML,'IT0')*100;
    AVERR(IT)=(SUM(OUTLETC,(COUTRR(OUTLETC,IT))$(not FIXCOUT(OUTLETC))+100$(FIXCOUT(OUTLETC)))+SUM(CONNECT,(FOORR(CONNECT,IT))$(not FIXFOO(CONNECT))+100$(FIXFOO(CONNECT)))+SUM(FML,(FTSRR(FML,IT))$(not FIXFTS(FML))+100$(FIXFTS(FML))))/(PARTVAR+DISGVAR);
    FIXCOUT(OUTLETC)=yes$(COUTRR(OUTLETC,IT) GT 99.9999);
    FIXFOO(CONNECT)=yes$(FOORR(CONNECT,IT) GT 99.9999);
    FIXFTS(FML)=yes$(FTSRR(FML,IT) GT 99.9999);
    NFIXV(IT)=SUM(OUTLETC$(FIXCOUT(OUTLETC)),1)+SUM(CONNECT$(FIXFOO(CONNECT)),1)+SUM(FML$(FIXFTS(FML)),1);
    NBOUNDP=2*(PARTVAR+DISGVAR-card(FIXCOUT)-card(FIXFOO)-card(FIXFTS));
    #Recompute Bounds
    SOLVE WUNLP using LP minimizing OBJ;
    MIPCPU=MAX(WUNLP.resusd,0.001);
    CIN.l(I,C)$(ACTOP(I) and FML(I) and ACTCN(C) and FIO.l(I) NE 0) = (Sum(k$(ACTWS(k)), FFW.l(k,i)*COFW(k,c)) + Sum(IL$(ACTOP(IL)), FOO.l(IL,i)*(COUT.l(IL,c)$(FML(IL))+COUTMAX(IL,c)$(FFO(IL))))) / FIO.l(i)  ;
    LBOUND=WUNLP.objval;
    SOLVE WUNNLP using NLP minimizing OBJ;
        if(WUNNLP.modelstat EQ 2,
            if(WUNNLP.objval LT UBOUND*0.9999,
            UBOUND=WUNNLP.objval;
            CURIT=CURIT+1; #One more round of LP-based OBBT
            );
        );
    OPTGAP=ABS((UBOUND-LBOUND)/LBOUND)*100;
    put TimeElapsed,LBOUND:9:5,UBOUND:9:5,OPTGAP:9:5,'LP relaxation after LP-based OBBT' /;
    putclose Search;
    );
CLOOP$(OPTGAP LE TRGGAP or TimeElapsed GT TCPU)=0;
BASIS(KD)=0;KVAL(KD,I,C)=0;ACTJ(J,KD)=no;JVAL(J,KD)=0;
#Relaxation with mixed-radix MDT
    loop(KDLL$(CLOOP),    
    BASIS(KDLL)=2;
    NDIGITS=SUM(KD$(BASIS(KD) GT 0),1);
    CARDK(I,C)$(FML(I) and ACTCN(C))=NDIGITS;
    ACTK(KD,I,C)$(FML(I) and ACTCN(C))=yes$(ord(KD) LE CARDK(I,C));
    KVAL(KD,I,C)$(ACTK(KD,I,C))=-BASIS(KD);
    ACTJ(J,KD)$(ord(KD) LE NDIGITS)=yes$(ord(J) LE BASIS(KD));
    JVAL(J,KD)$(ACTJ(J,KD))=ord(J)-1;
    LCOUTR.up(I,C)$(FML(I) and ACTCN(C))=MIN(1,PROD(KD$(ACTK(KD,I,C)),1/BASIS(KD)));
    put OptFile;
    put 'tilim ', (TCPU-TimeElapsed):<9:0 /;
    put 'Solnpool solnpool.gdx' /;
    put 'SolnPoolPop 1' /;
    put 'SolnPoolIntensity 1' /;
    put 'SolnPoolCapacity 10' /;
    put 'SolnPoolGap 0.01' /;
    put 'SolnPoolReplace 2' /;
    putclose OptFile;
    MDT.optfile=9;
    SOLVE MDT using MIP minimizing OBJ;
    ESTRELX=sqr(MDT.resusd)/MIPCPU;
    MIPCPU=MDT.resusd;
    ESTOBBT=MIPCPU*NBOUNDP;
        if(MDT.modelstat NE 10 and MDT.modelstat NE 14,
        LBOUND=MAX(LBOUND,MDT.objest);
        OPTGAP=ABS((UBOUND-LBOUND)/LBOUND)*100;
        );
        if(OPTGAP GT TRGGAP,
        #Solution pool to generate multiple starting points for NLP
        put SPoolRes;
        put 'Solutions from CPLEX pool' /;
        put 'File','LB','UB' /;
        execute_load 'solnpool.gdx', SPOOL=Index; CARDSP=card(SPOOL);
            loop(SPOOL(SP),
            put_utility fsoln 'gdxin' / SPOOL.te(SP):0:0;
            execute_loadpoint;
            LBSP(SP)=OBJ.l;
            put SPoolRes SP.te(SP),LBSP(SP):9:5;
            CIN.l(I,C)$(ACTOP(I) and FML(I) and ACTCN(C) and FIO.l(I) NE 0) = (Sum(k$(ACTWS(k)), FFW.l(k,i)*COFW(k,c)) + Sum(IL$(ACTOP(IL)), FOO.l(IL,i)*(COUT.l(IL,c)$(FML(IL))+COUTMAX(IL,c)$(FFO(IL))))) / FIO.l(i)  ;
            SOLVE WUNNLP using NLP minimizing OBJ;
                if(WUNNLP.modelstat EQ 2,
                UBSP(SP)=OBJ.l; put UBSP(SP):9:5 /;
                UBOUND=MIN(WUNNLP.objval,UBOUND);
                else
                put 'Infeasible' /;
                );
            );
        putclose SPoolRes;
        OPTGAP=ABS((UBOUND-LBOUND)/LBOUND)*100;
        );        
    put Search;
    put TimeElapsed,LBOUND:9:5,UBOUND:9:5,OPTGAP:9:5,'MIP Relaxation',NDIGITS:9:0 /;
    putclose Search;
    CLOOP$(OPTGAP LE TRGGAP or TimeElapsed GT TCPU or MIPCPU GT CPURELX)=0;
        if(CLOOP and ESTOBBT LT CPUOBBT and ESTRELX*NBOUNDP GT CPUOBBT/2,
        #NMDT OBBT with current settings
        TSTART=TimeElapsed;
            loop(OUTLETC$(not FIXCOUT(OUTLETC)),
            BNDCOUT(OUTLETC)=yes;
            SOLVE OBBT_MDT using MIP minimizing OBJ;
                if(OBBT_MDT.modelstat NE 10 and OBBT_MDT.modelstat NE 14 and OBBT_MDT.solvestat NE 4,
                LBCOUT(OUTLETC,CURRENT)=MIN(OBBT_MDT.objest,COUTMAX(OUTLETC)); #Avoids numerical errors
                else
                SOLVE OBBT_LP using LP minimizing OBJ;
                    if(OBBT_LP.modelstat NE 4,LBCOUT(OUTLETC,CURRENT)=MIN(OBJ.l,COUTMAX(OUTLETC)););
                );
            COUTMIN(OUTLETC)=SUM(CURRENT,LBCOUT(OUTLETC,CURRENT));COUT.lo(OUTLETC)=COUTMIN(OUTLETC);
            SOLVE OBBT_MDT using MIP maximizing OBJ;
                if(OBBT_MDT.modelstat NE 10 and OBBT_MDT.modelstat NE 14 and OBBT_MDT.solvestat NE 4,
                UBCOUT(OUTLETC,CURRENT)=MAX(OBBT_MDT.objest,COUTMIN(OUTLETC)); #Avoids numerical errors
                else
                SOLVE OBBT_LP using LP maximizing OBJ;
                    if(OBBT_LP.modelstat NE 4,UBCOUT(OUTLETC,CURRENT)=MAX(OBJ.l,COUTMIN(OUTLETC)););
                );
            COUTMAX(OUTLETC)=SUM(CURRENT,UBCOUT(OUTLETC,CURRENT));COUT.up(OUTLETC)=COUTMAX(OUTLETC);
            BNDCOUT(OUTLETC)=no;
            );
            loop(CONNECT$(not FIXFOO(CONNECT)),
            BNDFOO(CONNECT)=yes;
            SOLVE OBBT_MDT using MIP minimizing OBJ;
                if(OBBT_MDT.modelstat NE 10 and OBBT_MDT.modelstat NE 14 and OBBT_MDT.solvestat NE 4,
                LBFOO(CONNECT,CURRENT)=MIN(OBBT_MDT.objest,FOOMAX(CONNECT)); #Avoids numerical errors
                else
                SOLVE OBBT_LP using LP minimizing OBJ;
                    if(OBBT_LP.modelstat NE 4,LBFOO(CONNECT,CURRENT)=MIN(OBJ.l,FOOMAX(CONNECT)););
                );
            FOOMIN(CONNECT)=SUM(CURRENT,LBFOO(CONNECT,CURRENT));FOO.lo(CONNECT)=FOOMIN(CONNECT);
            FMIN(FML)=SUM(ACTOP,FOOMIN(FML,ACTOP))+FTSMIN(FML);
            SOLVE OBBT_MDT using MIP maximizing OBJ;
                if(OBBT_MDT.modelstat NE 10 and OBBT_MDT.modelstat NE 14 and OBBT_MDT.solvestat NE 4,
                UBFOO(CONNECT,CURRENT)=MAX(OBBT_MDT.objest,FOOMIN(CONNECT)); #Avoids numerical errors
                else
                SOLVE OBBT_LP using LP maximizing OBJ;
                    if(OBBT_LP.modelstat NE 4,UBFOO(CONNECT,CURRENT)=MAX(OBJ.l,FOOMIN(CONNECT)););
                );
            FOOMAX(CONNECT)=SUM(CURRENT,UBFOO(CONNECT,CURRENT));FOO.up(CONNECT)=FOOMAX(CONNECT);
            FMAX(FML)=MIN(FLIM(FML),SUM(ACTOP,FOOMAX(FML,ACTOP))+FTSMAX(FML));
            BNDFOO(CONNECT)=no;
            );
            loop(FML$(not FIXFTS(FML)),
            BNDFTS(FML)=yes;
            SOLVE OBBT_MDT using MIP minimizing OBJ;
                if(OBBT_MDT.modelstat NE 10 and OBBT_MDT.modelstat NE 14 and OBBT_MDT.solvestat NE 4,
                LBFTS(FML,CURRENT)=MIN(OBBT_MDT.objest,FTSMAX(FML)); #Avoids numerical errors
                else
                SOLVE OBBT_LP using LP minimizing OBJ;
                    if(OBBT_LP.modelstat NE 4,LBFTS(FML,CURRENT)=MIN(OBJ.l,FTSMAX(FML)););
                );
            FTSMIN(FML)=SUM(CURRENT,LBFTS(FML,CURRENT));FTS.lo(FML)=FTSMIN(FML);
            FMIN(FML)=SUM(ACTOP,FOOMIN(FML,ACTOP))+FTSMIN(FML);
            SOLVE OBBT_MDT using MIP maximizing OBJ;
                if(OBBT_MDT.modelstat NE 10 and OBBT_MDT.modelstat NE 14 and OBBT_MDT.solvestat NE 4,
                UBFTS(FML,CURRENT)=MAX(OBBT_MDT.objest,FTSMIN(FML)); #Avoids numerical errors
                else
                SOLVE OBBT_LP using LP maximizing OBJ;
                    if(OBBT_LP.modelstat NE 4,UBFTS(FML,CURRENT)=MAX(OBJ.l,FTSMIN(FML)););
                );
            FTSMAX(FML)=SUM(CURRENT,UBFTS(FML,CURRENT));FTS.up(FML)=FTSMAX(FML);
            FMAX(FML)=MIN(FLIM(FML),SUM(ACTOP,FOOMAX(FML,ACTOP))+FTSMAX(FML));
            BNDFTS(FML)=no;
            );
        #Report results of OBBT
        OBBTCPU=TimeElapsed-TSTART;CPUOBBT=MAX(0,CPUOBBT-OBBTCPU);
        CURIT=CURIT+1;
        COUTRNG(OUTLETC,IT)$(ord(IT) EQ CURIT+1)=SUM(CURRENT,UBCOUT(OUTLETC,CURRENT)-LBCOUT(OUTLETC,CURRENT));
        FOORNG(CONNECT,IT)$(ord(IT) EQ CURIT+1)=SUM(CURRENT,UBFOO(CONNECT,CURRENT)-LBFOO(CONNECT,CURRENT));
        FTSRNG(FML,IT)$(ord(IT) EQ CURIT+1)=SUM(CURRENT,UBFTS(FML,CURRENT)-LBFTS(FML,CURRENT));
        COUTRR(OUTLETC,IT)$(ord(IT) EQ CURIT+1)=100$(COUTRNG(OUTLETC,'IT0') EQ 0)+((COUTRNG(OUTLETC,'IT0')-COUTRNG(OUTLETC,IT))/COUTRNG(OUTLETC,'IT0')*100)$(COUTRNG(OUTLETC,'IT0') GT 0);
        FOORR(CONNECT,IT)$(ord(IT) EQ CURIT+1)=(FOORNG(CONNECT,'IT0')-FOORNG(CONNECT,IT))/FOORNG(CONNECT,'IT0')*100;
        FTSRR(FML,IT)$(ord(IT) EQ CURIT+1)=(FTSRNG(FML,'IT0')-FTSRNG(FML,IT))/FTSRNG(FML,'IT0')*100;
        FIXCOUT(OUTLETC)=yes$(SUM(IT$(ord(IT) EQ CURIT+1),COUTRR(OUTLETC,IT)) GT 99.9999);
        AVERR(IT)$(ord(IT) EQ CURIT+1)=(SUM(OUTLETC,(COUTRR(OUTLETC,IT))$(not FIXCOUT(OUTLETC))+100$(FIXCOUT(OUTLETC)))+SUM(CONNECT,(FOORR(CONNECT,IT))$(not FIXFOO(CONNECT))+100$(FIXFOO(CONNECT)))+SUM(FML,(FTSRR(FML,IT))$(not FIXFTS(FML))+100$(FIXFTS(FML))))/(PARTVAR+DISGVAR);
        FIXFOO(CONNECT)=yes$(SUM(IT$(ord(IT) EQ CURIT+1),FOORR(CONNECT,IT)) GT 99.9999);
        FIXFTS(FML)=yes$(SUM(IT$(ord(IT) EQ CURIT+1),FTSRR(FML,IT)) GT 99.9999);
        NFIXV(IT)$(ord(IT) EQ CURIT+1)=SUM(OUTLETC$(FIXCOUT(OUTLETC)),1)+SUM(CONNECT$(FIXFOO(CONNECT)),1)+SUM(FML$(FIXFTS(FML)),1);
        NBOUNDP=2*(PARTVAR+DISGVAR-card(FIXCOUT)-card(FIXFOO)-card(FIXFTS));
        #Check if with results from McCormick relaxation we can already stop
        SOLVE WUNLP using LP minimizing OBJ;
        LBOUND=MAX(LBOUND,WUNLP.objval);
        OPTGAP=ABS((UBOUND-LBOUND)/LBOUND)*100;
        CLOOP$(OPTGAP LE TRGGAP)=0;
        put TimeElapsed,LBOUND:9:5,UBOUND:9:5,OPTGAP:9:5,'LP relaxation after MIP-based OBBT',NDIGITS:9:0 /;
        putclose Search;
        SBBOBBT=ord(KDLL);
        );
    SBBRELX=ord(KDLL);
    );
WHATSBB$(SBBOBBT EQ 0)=1;
#Spatial Branch and Bound
    if(TimeElapsed LT TCPU and OPTGAP GT TRGGAP,
    put 'Spatial B&B' /;
    put 'CPU','LB','UB','Gap (%)','Expl.','Open'/;
    #Compute the error that will be used to select branching variavel (assumes MIP relaxation)
    ERRFOO(IL,I,C)$(FML(IL) and ACTOP(I) and ACTCN(C))=ABS(LFOOC.l(IL,I,C)-FOO.l(IL,I)*LCOUT.l(IL,C));
    ERRFTS(FML,C)$(ACTCN(C))=ABS(LFTSC.l(FML,C)-FTS.l(FML)*LCOUT.l(FML,C));
    display ERRFOO,ERRFTS;
    BOUND('ND1')=LBOUND;
    SRCNODE('ND1')=yes;
    LOGBB('ND1','Time (s)')=TimeElapsed;LOGBB('ND1','Explored')=card(SRCNODE);LOGBB('ND1','Remaining')=card(WAITING);
    LOGBB('ND1','UBOUND')=UBOUND;LOGBB('ND1','LBOUND')=LBOUND;LOGBB('ND1','OPTGAP (%)')=OPTGAP;
    LOGND(CURRENT,CONNECT,'XLO')=LBFOO(CONNECT,CURRENT);LOGND(CURRENT,CONNECT,'XUP')=UBFOO(CONNECT,CURRENT);
    LOGND2(CURRENT,FML,'XLO')=LBFTS(FML,CURRENT);LOGND2(CURRENT,FML,'XUP')=UBFTS(FML,CURRENT);
    LOGND3(CURRENT,OUTLETC,'XLO')=LBCOUT(OUTLETC,CURRENT);LOGND3(CURRENT,OUTLETC,'XUP')=UBCOUT(OUTLETC,CURRENT);
    put LOGBB('ND1','Time (s)'):9:2, put LOGBB('ND1','LBOUND'):9:5, put LOGBB('ND1','UBOUND'):9:5, put LOGBB('ND1','OPTGAP (%)'):9:5, put LOGBB('ND1','Explored'):9:0 ; put LOGBB('ND1','Remaining'):9:0 /;
    putclose Search;
    #Choose branching variavel
    FIRST2=1;
        if(SMAX((CONNECT,C)$(ACTCN(C)),ERRFOO(CONNECT,C)) GE SMAX((FML,C)$(ACTCN(C)),ERRFTS(FML,C)),
            loop(CONNECT$(FIRST2 and SMAX(C$(ACTCN(C)),ERRFOO(CONNECT,C)) EQ SMAX((I,IL,C)$(FML(I) and ACTOP(IL) and ACTCN(C)),ERRFOO(I,IL,C))),
            NEXTFOO(CURRENT,CONNECT)=yes;
            FIRST2=0;
            );
        else
            loop(FML$(FIRST2 and SMAX(C$(ACTCN(C)),ERRFTS(FML,C)) EQ SMAX((IL,C)$(FML(IL) and ACTCN(C)),ERRFTS(IL,C))),
            NEXTFTS(CURRENT,FML)=yes;
            FIRST2=0;
            );
        );
        loop(CURRENT,BRCHFOO(CONNECT)=NEXTFOO(CURRENT,CONNECT);BRCHFTS(FML)=NEXTFTS(CURRENT,FML););
    LOGND(CURRENT,BRCHFOO,'LBOUND')=LBOUND;LOGND(CURRENT,BRCHFOO,'UBOUND')=UBOUND;LOGND(CURRENT,BRCHFOO,'Next From')=SUM((I,IL)$(NEXTFOO(CURRENT,I,IL)),ord(I));LOGND(CURRENT,BRCHFOO,'Next To')=SUM((I,IL)$(NEXTFOO(CURRENT,I,IL)),ord(IL));
    LOGND2(CURRENT,BRCHFTS,'LBOUND')=LBOUND;LOGND2(CURRENT,BRCHFTS,'UBOUND')=UBOUND;LOGND2(CURRENT,BRCHFTS,'Next From')=SUM(I$(NEXTFTS(CURRENT,I)),ord(I));LOGND2(CURRENT,BRCHFTS,'Next To')=SUM(I$(NEXTFTS(CURRENT,I)),ord(I));
    #Perform branching
        loop(ND$(not DONE), # and ord(ND) LE 1
        #Node selection
        CURRENT(NDL)=no;
        CURRENT(WAITING(NDL))$(BOUND(NDL) EQ LBOUND)=yes;
        FIRST=1;
            loop(CURRENT$FIRST, #Select only one node for branching
            WAITING(CURRENT)=no;
                loop(NDL$(CURRENT(NDL)),BRCHFOO(I,IL)=NEXTFOO(NDL,I,IL);BRCHFTS(I)=NEXTFTS(NDL,I);); #Assign next variavel to branch
            LOGIT(ND,'Node_Sel.')=SUM(NDL$(CURRENT(NDL)),ord(NDL));
            LOGIT(ND,'From_Unit')=SUM((I,IL)$(BRCHFOO(I,IL)),ord(I))+SUM(I$(BRCHFTS(I)),ord(I));
            LOGIT(ND,'To_Unit')=SUM((I,IL)$(BRCHFOO(I,IL)),ord(IL))+SUM(I$(BRCHFTS(I)),ord(I));
            #Branch
                for(BISECT=1 to 2,
                NEWNODE(NDL)=NEWNODE(NDL-1);
                NODEBRN(CURRENT,NEWNODE)=yes;
                SRCNODE(NEWNODE)=yes;
                WAITING(NEWNODE)=yes;
                LBFOO(CONNECT,NEWNODE)=LBFOO(CONNECT,CURRENT);UBFOO(CONNECT,NEWNODE)=UBFOO(CONNECT,CURRENT);
                LBFTS(FML,NEWNODE)=LBFTS(FML,CURRENT);UBFTS(FML,NEWNODE)=UBFTS(FML,CURRENT);
                LBCOUT(OUTLETC,NEWNODE)=LBCOUT(OUTLETC,CURRENT);UBCOUT(OUTLETC,NEWNODE)=UBCOUT(OUTLETC,CURRENT);
                    if(BISECT EQ 1, #Left branch
                    LBFOO(BRCHFOO,NEWNODE)=LBFOO(BRCHFOO,CURRENT);
                    UBFOO(BRCHFOO,NEWNODE)=LBFOO(BRCHFOO,CURRENT)+(UBFOO(BRCHFOO,CURRENT)-LBFOO(BRCHFOO,CURRENT))/2;
                    LBFTS(BRCHFTS,NEWNODE)=LBFTS(BRCHFTS,CURRENT);
                    UBFTS(BRCHFTS,NEWNODE)=LBFTS(BRCHFTS,CURRENT)+(UBFTS(BRCHFTS,CURRENT)-LBFTS(BRCHFTS,CURRENT))/2;
                    else
                    LBFOO(BRCHFOO,NEWNODE)=LBFOO(BRCHFOO,CURRENT)+(UBFOO(BRCHFOO,CURRENT)-LBFOO(BRCHFOO,CURRENT))/2;
                    UBFOO(BRCHFOO,NEWNODE)=UBFOO(BRCHFOO,CURRENT);
                    LBFTS(BRCHFTS,NEWNODE)=LBFTS(BRCHFTS,CURRENT)+(UBFTS(BRCHFTS,CURRENT)-LBFTS(BRCHFTS,CURRENT))/2;
                    UBFTS(BRCHFTS,NEWNODE)=UBFTS(BRCHFTS,CURRENT);
                    );
                FOOMIN(CONNECT)=SUM(NDL$(NEWNODE(NDL)),LBFOO(CONNECT,NDL));FOOMAX(CONNECT)=SUM(NDL$(NEWNODE(NDL)),UBFOO(CONNECT,NDL));
                FTSMIN(FML)=SUM(NDL$(NEWNODE(NDL)),LBFTS(FML,NDL));FTSMAX(FML)=SUM(NDL$(NEWNODE(NDL)),UBFTS(FML,NDL));
                FMIN(FML)=SUM(ACTOP,FOOMIN(FML,ACTOP))+FTSMIN(FML);FMAX(FML)=MIN(FLIM(FML),SUM(ACTOP,FOOMAX(FML,ACTOP))+FTSMAX(FML));
                COUTMIN(OUTLETC)=SUM(NDL$(NEWNODE(NDL)),LBCOUT(OUTLETC,NDL));COUTMAX(OUTLETC)=SUM(NDL$(NEWNODE(NDL)),UBCOUT(OUTLETC,NDL));
                FOO.lo(CONNECT)=FOOMIN(CONNECT);FOO.up(CONNECT)=FOOMAX(CONNECT);
                FTS.lo(FML)=FTSMIN(FML);FTS.up(FML)=FTSMAX(FML);
                COUT.lo(OUTLETC)=COUTMIN(OUTLETC);COUT.up(OUTLETC)=COUTMAX(OUTLETC);
                LOGND(NEWNODE,CONNECT,'XLO')=FOO.lo(CONNECT);LOGND(NEWNODE,CONNECT,'XUP')=FOO.up(CONNECT);
                LOGND2(NEWNODE,FML,'XLO')=FTS.lo(FML);LOGND2(NEWNODE,FML,'XUP')=FTS.up(FML);
                LOGND3(NEWNODE,OUTLETC,'XLO')=COUT.lo(OUTLETC);LOGND3(NEWNODE,OUTLETC,'XUP')=COUT.up(OUTLETC);
                #OBBT
                STATUS=0;
                    loop(OUTLETC$(not STATUS and not FIXCOUT(OUTLETC)),
                    BNDCOUT(OUTLETC)=yes;
                        if(WHATSBB EQ 1,
                        #LP-based OBBT                        
                        SOLVE OBBT_LP using LP minimizing OBJ;
                            if(OBBT_LP.modelstat EQ 4, #Branch can be fathomed
                            STATUS=1;
                            WAITING(NEWNODE)=no;
                            else
                            LBCOUT(OUTLETC,NEWNODE)=MIN(OBJ.l,COUTMAX(OUTLETC)); #Avoids numerical errors
                            COUTMIN(OUTLETC)=SUM(NEWNODE,LBCOUT(OUTLETC,NEWNODE));
                            COUT.lo(OUTLETC)=COUTMIN(OUTLETC);
                            SOLVE OBBT_LP using LP maximizing OBJ;
                            UBCOUT(OUTLETC,NEWNODE)=MAX(OBJ.l,COUTMIN(OUTLETC)); #Avoids numerical errors
                            COUTMAX(OUTLETC)=SUM(NEWNODE,UBCOUT(OUTLETC,NEWNODE));
                            COUT.up(OUTLETC)=COUTMAX(OUTLETC);
                            LOGND3(NEWNODE,OUTLETC,'XLO')=COUT.lo(OUTLETC);LOGND3(NEWNODE,OUTLETC,'XUP')=COUT.up(OUTLETC);
                            );
                        else
                        #MIP-based OBBT
                        KVAL(KD,I,C)=0;ACTJ(J,KD)=no;JVAL(J,KD)=0;
                        BASIS(KD)=2$(ord(KD) LE SBBOBBT);
                        NDIGITS=SUM(KD$(BASIS(KD) GT 0),1);
                        CARDK(I,C)$(FML(I) and ACTCN(C))=NDIGITS;
                        ACTK(KD,I,C)$(FML(I) and ACTCN(C))=yes$(ord(KD) LE CARDK(I,C));
                        KVAL(KD,I,C)$(ACTK(KD,I,C))=-BASIS(KD);
                        ACTJ(J,KD)$(ord(KD) LE NDIGITS)=yes$(ord(J) LE BASIS(KD));
                        JVAL(J,KD)$(ACTJ(J,KD))=ord(J)-1;
                        LCOUTR.up(I,C)$(FML(I) and ACTCN(C))=MIN(1,PROD(KD$(ACTK(KD,I,C)),1/BASIS(KD)));
                        SOLVE OBBT_MDT using MIP minimizing OBJ;
                            if(OBBT_MDT.modelstat NE 10 and OBBT_MDT.modelstat NE 14 and OBBT_MDT.solvestat NE 4,
                            LBCOUT(OUTLETC,NEWNODE)=MIN(OBBT_MDT.objest,COUTMAX(OUTLETC)); #Avoids numerical errors
                            COUTMIN(OUTLETC)=SUM(NEWNODE,LBCOUT(OUTLETC,NEWNODE));COUT.lo(OUTLETC)=COUTMIN(OUTLETC);
                            else
                            SOLVE OBBT_LP using LP minimizing OBJ;
                                if(OBBT_LP.modelstat EQ 4, #Branch can be fathomed)
                                STATUS=1;
                                WAITING(NEWNODE)=no;
                                else
                                LBCOUT(OUTLETC,NEWNODE)=MIN(OBJ.l,COUTMAX(OUTLETC)); #Avoids numerical errors
                                COUTMIN(OUTLETC)=SUM(NEWNODE,LBCOUT(OUTLETC,NEWNODE));COUT.lo(OUTLETC)=COUTMIN(OUTLETC);
                                );
                            );
                            if(not STATUS,
                            SOLVE OBBT_MDT using MIP maximizing OBJ;
                                if(OBBT_MDT.modelstat NE 10 and OBBT_MDT.modelstat NE 14 and OBBT_MDT.solvestat NE 4,
                                UBCOUT(OUTLETC,NEWNODE)=MAX(OBBT_MDT.objest,COUTMIN(OUTLETC)); #Avoids numerical errors
                                COUTMAX(OUTLETC)=SUM(NEWNODE,UBCOUT(OUTLETC,NEWNODE));COUT.up(OUTLETC)=COUTMAX(OUTLETC);
                                LOGND3(NEWNODE,OUTLETC,'XLO')=COUT.lo(OUTLETC);LOGND3(NEWNODE,OUTLETC,'XUP')=COUT.up(OUTLETC);
                                else
                                SOLVE OBBT_LP using LP maximizing OBJ;
                                    if(OBBT_LP.modelstat EQ 4, #Branch can be fathomed)
                                    STATUS=1;
                                    WAITING(NEWNODE)=no;
                                    else
                                    UBCOUT(OUTLETC,NEWNODE)=MAX(OBJ.l,COUTMIN(OUTLETC)); #Avoids numerical errors
                                    COUTMAX(OUTLETC)=SUM(NEWNODE,UBCOUT(OUTLETC,NEWNODE));COUT.up(OUTLETC)=COUTMAX(OUTLETC);
                                    LOGND3(NEWNODE,OUTLETC,'XLO')=COUT.lo(OUTLETC);LOGND3(NEWNODE,OUTLETC,'XUP')=COUT.up(OUTLETC);
                                    );
                                );
                            );
                        );
                    BNDCOUT(OUTLETC)=no;
                    );
                    loop((ILL,ILLL)$(CONNECT(ILL,ILLL) and not BRCHFOO(ILL,ILLL) and not STATUS and not FIXFOO(ILL,ILLL)),
                    BNDFOO(ILL,ILLL)=yes;
                        if(WHATSBB EQ 1,
                        #LP-based OBBT                        
                        SOLVE OBBT_LP using LP minimizing OBJ;
                            if(OBBT_LP.modelstat EQ 4, #Branch can be fathomed
                            STATUS=1;
                            WAITING(NEWNODE)=no;
                            else
                            LBFOO(ILL,ILLL,NEWNODE)=MIN(OBJ.l,FOOMAX(ILL,ILLL)); #Avoids numerical errors
                            FOOMIN(ILL,ILLL)=SUM(NEWNODE,LBFOO(ILL,ILLL,NEWNODE));FOO.lo(ILL,ILLL)=FOOMIN(ILL,ILLL);
                            SOLVE OBBT_LP using LP maximizing OBJ;
                            UBFOO(ILL,ILLL,NEWNODE)=MAX(OBJ.l,FOOMIN(ILL,ILLL)); #Avoids numerical errors
                            FOOMAX(ILL,ILLL)=SUM(NEWNODE,UBFOO(ILL,ILLL,NEWNODE));FOO.up(ILL,ILLL)=FOOMAX(ILL,ILLL);
                            LOGND(NEWNODE,ILL,ILLL,'XLO')=FOO.lo(ILL,ILLL);LOGND(NEWNODE,ILL,ILLL,'XUP')=FOO.up(ILL,ILLL);
                            );
                        else
                        #MIP-based OBBT
                        SOLVE OBBT_MDT using MIP minimizing OBJ;
                            if(OBBT_MDT.modelstat NE 10 and OBBT_MDT.modelstat NE 14 and OBBT_MDT.solvestat NE 4,
                            LBFOO(ILL,ILLL,NEWNODE)=MIN(OBBT_MDT.objest,FOOMAX(ILL,ILLL)); #Avoids numerical errors
                            FOOMIN(ILL,ILLL)=SUM(NEWNODE,LBFOO(ILL,ILLL,NEWNODE));FOO.lo(ILL,ILLL)=FOOMIN(ILL,ILLL);
                            else
                            SOLVE OBBT_LP using LP minimizing OBJ;
                                if(OBBT_LP.modelstat EQ 4, #Branch can be fathomed)
                                STATUS=1;
                                WAITING(NEWNODE)=no;
                                else
                                LBFOO(ILL,ILLL,NEWNODE)=MIN(OBJ.l,FOOMAX(ILL,ILLL)); #Avoids numerical errors
                                FOOMIN(ILL,ILLL)=SUM(NEWNODE,LBFOO(ILL,ILLL,NEWNODE));FOO.lo(ILL,ILLL)=FOOMIN(ILL,ILLL);
                                );
                            );
                            if(not STATUS,
                            SOLVE OBBT_MDT using MIP maximizing OBJ;
                                if(OBBT_MDT.modelstat NE 10 and OBBT_MDT.modelstat NE 14 and OBBT_MDT.solvestat NE 4,
                                UBFOO(ILL,ILLL,NEWNODE)=MAX(OBBT_MDT.objest,FOOMIN(ILL,ILLL)); #Avoids numerical errors
                                FOOMAX(ILL,ILLL)=SUM(NEWNODE,UBFOO(ILL,ILLL,NEWNODE));FOO.up(ILL,ILLL)=FOOMAX(ILL,ILLL);
                                LOGND(NEWNODE,ILL,ILLL,'XLO')=FOO.lo(ILL,ILLL);LOGND(NEWNODE,ILL,ILLL,'XUP')=FOO.up(ILL,ILLL);                                
                                else
                                SOLVE OBBT_LP using LP maximizing OBJ;
                                    if(OBBT_LP.modelstat EQ 4, #Branch can be fathomed)
                                    STATUS=1;
                                    WAITING(NEWNODE)=no;
                                    else
                                    UBFOO(ILL,ILLL,NEWNODE)=MAX(OBJ.l,FOOMIN(ILL,ILLL)); #Avoids numerical errors
                                    FOOMAX(ILL,ILLL)=SUM(NEWNODE,UBFOO(ILL,ILLL,NEWNODE));FOO.up(ILL,ILLL)=FOOMAX(ILL,ILLL);
                                    LOGND(NEWNODE,ILL,ILLL,'XLO')=FOO.lo(ILL,ILLL);LOGND(NEWNODE,ILL,ILLL,'XUP')=FOO.up(ILL,ILLL);
                                    );
                                );
                            );                       
                        );
                    BNDFOO(ILL,ILLL)=no;
                    );
                    loop(ILL$(FML(ILL) and not BRCHFTS(ILL) and not STATUS and not FIXFTS(ILL)),
                    BNDFTS(ILL)=yes;
                        if(WHATSBB EQ 1,
                        #LP-based OBBT                        
                        SOLVE OBBT_LP using LP minimizing OBJ;
                            if(OBBT_LP.modelstat EQ 4, #Branch can be fathomed
                            STATUS=1;
                            WAITING(NEWNODE)=no;
                            else
                            LBFTS(ILL,NEWNODE)=MIN(OBJ.l,FTSMAX(ILL)); #Avoids numerical errors
                            FTSMIN(ILL)=SUM(NEWNODE,LBFTS(ILL,NEWNODE));FTS.lo(ILL)=FTSMIN(ILL);
                            SOLVE OBBT_LP using LP maximizing OBJ;
                            UBFTS(ILL,NEWNODE)=MAX(OBJ.l,FTSMIN(ILL)); #Avoids numerical errors
                            FTSMAX(ILL)=SUM(NEWNODE,UBFTS(ILL,NEWNODE));FTS.up(ILL)=FTSMAX(ILL);
                            LOGND2(NEWNODE,ILL,'XLO')=FTS.lo(ILL);LOGND2(NEWNODE,ILL,'XUP')=FTS.up(ILL);
                            );
                        else
                        #MIP-based OBBT
                        SOLVE OBBT_MDT using MIP minimizing OBJ;
                            if(OBBT_MDT.modelstat NE 10 and OBBT_MDT.modelstat NE 14 and OBBT_MDT.solvestat NE 4,
                            LBFTS(ILL,NEWNODE)=MIN(OBBT_MDT.objest,FTSMAX(ILL)); #Avoids numerical errors
                            FTSMIN(ILL)=SUM(NEWNODE,LBFTS(ILL,NEWNODE));FTS.lo(ILL)=FTSMIN(ILL);
                            else
                            SOLVE OBBT_LP using LP minimizing OBJ;
                                if(OBBT_LP.modelstat EQ 4, #Branch can be fathomed)
                                STATUS=1;
                                WAITING(NEWNODE)=no;
                                else
                                LBFTS(ILL,NEWNODE)=MIN(OBJ.l,FTSMAX(ILL)); #Avoids numerical errors
                                FTSMIN(ILL)=SUM(NEWNODE,LBFTS(ILL,NEWNODE));FTS.lo(ILL)=FTSMIN(ILL);
                                );
                            );
                            if(not STATUS,
                            SOLVE OBBT_MDT using MIP maximizing OBJ;
                                if(OBBT_MDT.modelstat NE 10 and OBBT_MDT.modelstat NE 14 and OBBT_MDT.solvestat NE 4,
                                UBFTS(ILL,NEWNODE)=MAX(OBBT_MDT.objest,FTSMIN(ILL)); #Avoids numerical errors
                                FTSMAX(ILL)=SUM(NEWNODE,UBFTS(ILL,NEWNODE));FTS.up(ILL)=FTSMAX(ILL);
                                LOGND2(NEWNODE,ILL,'XLO')=FTS.lo(ILL);LOGND2(NEWNODE,ILL,'XUP')=FTS.up(ILL);                                
                                else
                                SOLVE OBBT_LP using LP maximizing OBJ;
                                    if(OBBT_LP.modelstat EQ 4, #Branch can be fathomed)
                                    STATUS=1;
                                    WAITING(NEWNODE)=no;
                                    else
                                    UBFTS(ILL,NEWNODE)=MAX(OBJ.l,FTSMIN(ILL)); #Avoids numerical errors
                                    FTSMAX(ILL)=SUM(NEWNODE,UBFTS(ILL,NEWNODE));FTS.up(ILL)=FTSMAX(ILL);
                                    LOGND2(NEWNODE,ILL,'XLO')=FTS.lo(ILL);LOGND2(NEWNODE,ILL,'XUP')=FTS.up(ILL);
                                    );
                                );
                            );
                        );
                    BNDFTS(ILL)=no;
                    );
                #MDT Relaxation
                    if(not STATUS,
                    KVAL(KD,I,C)=0;ACTJ(J,KD)=no;JVAL(J,KD)=0;
                    BASIS(KD)=2$(ord(KD) LE SBBRELX);
                    NDIGITS=SUM(KD$(BASIS(KD) GT 0),1);
                    CARDK(I,C)$(FML(I) and ACTCN(C))=NDIGITS;
                    ACTK(KD,I,C)$(FML(I) and ACTCN(C))=yes$(ord(KD) LE CARDK(I,C));
                    KVAL(KD,I,C)$(ACTK(KD,I,C))=-BASIS(KD);
                    ACTJ(J,KD)$(ord(KD) LE NDIGITS)=yes$(ord(J) LE BASIS(KD));
                    JVAL(J,KD)$(ACTJ(J,KD))=ord(J)-1;
                    LCOUTR.up(I,C)$(FML(I) and ACTCN(C))=MIN(1,PROD(KD$(ACTK(KD,I,C)),1/BASIS(KD)));
                    #No need for solution pool at this stage
                    put OptFile;
                    put 'tilim ', (MAX(0,TCPU-TimeElapsed)):<9:0 /;
                    putclose OptFile;
                    MDT.cutoff=UBOUND; #*(1-TRGGAP*0.01); #Avoids issues with MILP solver
                    SOLVE MDT using MIP minimizing OBJ;
                        if(MDT.modelstat EQ 10,BOUND(NEWNODE)=+INF;
                        else BOUND(NEWNODE)=MAX(MDT.objest,BOUND(CURRENT));
                            if(MDT.modelstat NE 14,
                            ERRFOO(IL,I,C)$(FML(IL) and ACTOP(I) and ACTCN(C))=ABS(LFOOC.l(IL,I,C)-FOO.l(IL,I)*LCOUT.l(IL,C));
                            ERRFTS(FML,C)$(ACTCN(C))=ABS(LFTSC.l(FML,C)-FTS.l(FML)*LCOUT.l(FML,C));
                            else  #Use LP relaxation to compute the error that will be used to select branching variable
                            SOLVE WUNLP using LP minimizing OBJ;
                            ERRFOO(IL,I,C)$(FML(IL) and ACTOP(I) and ACTCN(C))=ABS(ZFOO.l(IL,I,C)-FOO.l(IL,I)*COUT.l(IL,C));
                            ERRFTS(FML,C)$(ACTCN(C))=ABS(ZFTS.l(FML,C)-FTS.l(FML)*COUT.l(FML,C));
                            );
                        );
                    LOGND(NEWNODE,BRCHFOO,'LBOUND')=BOUND(NEWNODE);LOGND2(NEWNODE,BRCHFTS,'LBOUND')=BOUND(NEWNODE);
                        if(SUM(NEWNODE,BOUND(NEWNODE)) GT UBOUND/(1+TRGGAP*0.01),
                        WAITING(NEWNODE)=no; #Fathom node
                        else
                        #Choose branching variavel
                        FIRST2=1;
                            if(SMAX((CONNECT,C)$(ACTCN(C)),ERRFOO(CONNECT,C)) GE SMAX((FML,C)$(ACTCN(C)),ERRFTS(FML,C)),
                                loop(CONNECT$(FIRST2 and SMAX(C$(ACTCN(C)),ERRFOO(CONNECT,C)) EQ SMAX((I,IL,C)$(FML(I) and ACTOP(IL) and ACTCN(C)),ERRFOO(I,IL,C))),
                                NEXTFOO(NEWNODE,CONNECT)=yes;
                                FIRST2=0;
                                );
                            else
                                loop(FML$(FIRST2 and SMAX(C$(ACTCN(C)),ERRFTS(FML,C)) EQ SMAX((IL,C)$(FML(IL) and ACTCN(C)),ERRFTS(IL,C))),
                                NEXTFTS(NEWNODE,FML)=yes;
                                FIRST2=0;
                                );
                            );
                        LOGND(NEWNODE,BRCHFOO,'Next From')=SUM((I,IL)$(NEXTFOO(NEWNODE,I,IL)),ord(I))+SUM(I$(NEXTFTS(NEWNODE,I)),ord(I));
                        LOGND(NEWNODE,BRCHFOO,'Next To')=SUM((I,IL)$(NEXTFOO(NEWNODE,I,IL)),ord(IL))+SUM(I$(NEXTFTS(NEWNODE,I)),ord(I));
                        LOGND2(NEWNODE,BRCHFTS,'Next From')=SUM((I,IL)$(NEXTFOO(NEWNODE,I,IL)),ord(I))+SUM(I$(NEXTFTS(NEWNODE,I)),ord(I));
                        LOGND2(NEWNODE,BRCHFTS,'Next To')=SUM((I,IL)$(NEXTFOO(NEWNODE,I,IL)),ord(IL))+SUM(I$(NEXTFTS(NEWNODE,I)),ord(I));
                        FOO.lo(CONNECT)=LBFOO(CONNECT,'ND1');FOO.up(CONNECT)=UBFOO(CONNECT,'ND1');
                        FTS.lo(FML)=LBFTS(FML,'ND1');FTS.up(FML)=UBFTS(FML,'ND1');
                        CIN.l(I,C)$(ACTOP(I)  and FML(I) and ACTCN(C) and FIO.l(I) NE 0) = (Sum(k$(ACTWS(k)), FFW.l(k,i)*COFW(k,c)) + Sum(IL$(ACTOP(IL)), FOO.l(IL,i)*(COUT.l(IL,c)$(FML(IL))+COUTMAX(IL,c)$(FFO(IL))))) / FIO.l(i)  ;
                        SOLVE WUNNLP using NLP minimizing OBJ;
                            if(WUNNLP.modelstat EQ 2,
                            LOGND(NEWNODE,BRCHFOO,'UBOUND')=WUNNLP.objval;LOGND2(NEWNODE,BRCHFTS,'UBOUND')=WUNNLP.objval;
                                if(WUNNLP.objval LT UBOUND,
                                UBOUND=WUNNLP.objval;
                                    loop(NDL$(WAITING(NDL)), #Eliminate nodes with lower bound greater than improved upper bound
                                        if(BOUND(NDL) GT UBOUND/(1+TRGGAP*0.01),WAITING(NDL)=no;);
                                    );
                                );
                            else
                            LOGND(NEWNODE,BRCHFOO,'UBOUND')=UBOUND;LOGND2(NEWNODE,BRCHFTS,'UBOUND')=UBOUND;
                            );
                        );
                    else
                    LOGND(NEWNODE,BRCHFOO,'LBOUND')=+INF;
                    LOGND2(NEWNODE,BRCHFTS,'LBOUND')=+INF;
                    );
                );
            FIRST=0;
            );
        LBOUND=MIN(SMIN(WAITING(NDL),BOUND(NDL)),UBOUND*(1-TRGGAP*0.01));
        OPTGAP=ABS((UBOUND-LBOUND)/LBOUND)*100;
        LOGBB(ND+1,'Time (s)')=TimeElapsed;LOGBB(ND+1,'Explored')=card(SRCNODE);LOGBB(ND+1,'Remaining')=card(WAITING);
        LOGBB(ND+1,'UBOUND')=UBOUND;LOGBB(ND+1,'LBOUND')=LBOUND;LOGBB(ND+1,'OPTGAP (%)')=OPTGAP;
        LOGIT(ND,'LBOUND')=LBOUND;
        LOGIT(ND,'UBOUND')=UBOUND;
        LOGIT(ND,'OPTGAP (%)')=OPTGAP;
        LOGIT(ND,'WAITING')=card(WAITING);
        LOGIT(ND,'CPUs')=TimeElapsed;
        put Search;
        put LOGBB(ND+1,'Time (s)'):9:2, put LOGBB(ND+1,'LBOUND'):9:5, put LOGBB(ND+1,'UBOUND'):9:5, put LOGBB(ND+1,'OPTGAP (%)'):9:5, put LOGBB(ND+1,'Explored'):9:0, put LOGBB(ND+1,'Remaining'):9:0 /;
        putclose Search;
        DONE$(card(WAITING) EQ 0 or OPTGAP LE TRGGAP or TimeElapsed GE TCPU)=1;
        );
    display LOGBB,LOGND,LOGND2,LOGND3,LOGIT,NODEBRN,WAITING;
    );
put Results;
put 'OBBT'/;
put 'Var.', 'RangeR (%)'/;
          loop(OUTLETC,
          put OUTLETC.te(OUTLETC);
                    loop(IT$(ord(IT) GT 1 and ord(IT) LE CURIT+1),
                    put COUTRR(OUTLETC,IT):8:4;
                    );
          put/;
          );
          loop(CONNECT,
          put CONNECT.te(CONNECT);
                    loop(IT$(ord(IT) GT 1 and ord(IT) LE CURIT+1),
                    put FOORR(CONNECT,IT):8:4;
                    );
          put/;
          );
          loop(FML,
          put FML.te(FML);
                    loop(IT$(ord(IT) GT 1 and ord(IT) LE CURIT+1),
                    put FTSRR(FML,IT):8:4;
                    );
          put/;
          );
put 'Fixed V',loop(IT$(ord(IT) GT 1 and ord(IT) LE CURIT+1), put NFIXV(IT);); put/;
put 'ARR',loop(IT$(ord(IT) GT 1 and ord(IT) LE CURIT+1), put AVERR(IT);); put/;
display NBOUNDP,UBOUND,LBOUND,OPTGAP;
);