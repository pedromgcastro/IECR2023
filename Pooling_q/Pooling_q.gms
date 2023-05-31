$Title Global optimization of QCPs arising from the pooling problem (q-formulation)

$Ontext
          Code for of the global optimization algorithm featured in the article submitted to Ind. Eng. Chem. Res. in June 2023.
          "Global Optimization of QCPs using MIP relaxations with a base-2 logarithmic partitioning scheme"
          NLP model taken from Ben-Tal et al. Math Programming 63 (1994), 193-212.
          This is the q-formulation. Relaxation values match PQ-relaxation in Alfaki and Haugland JOGO (2013) 56: 897-916.

          WHATM=    7         BARON
                    8         GUROBI
                    9         Mixed-radix MDT algorithm with spatial B&B

          WHATSBB=  1         LP-based bound contracting in spatial B&B
                    2         MIP-based bound contracting in spatial B&B
$Offtext

$eolcom #

SETS
S         Sources   /S1*S60/
L         Pools     /L1*L30/
P         Products  /P1*P40/
K         Qualities /K1*K20/
J         Numbers   /0*29/
KD        Digits  /KD1*KD21/
IT        Iterations        /IT0*IT100/
SP        Solutions in pool /file1*file100/
ND        Nodes for spatial B&B         /ND1*ND20000/
ACTS(S)   Active sources
ACTL(L)   Active pools
ACTP(P)   Active products
ACTK(K)   Active qualities
SPOOL(SP) Actual solutions in pool
ACTJ(J,KD)          Active numbers J in digit KD
ACTSL(S,L)          Active connections from source S to pool L
ACTSP(S,P)          Active connections from source S to product P
ACTLP(L,P)          Active connections from pool L to product P
ACTKD(KD,S,L)       Active significant digits KD for fraction of source S in pool L
ACTSLJK(S,L,J,KD)   Active connections from source S to pool L number J digit KD
BNDQ(S,L) Bound contract is done for fraction from source S in pool L
BNDF(L,P) Bound contracting is done for flow from pool L to product P
BRCHF(L,P)          Flowrate from pool L product P being branched
FRAC(S,L) Fractions to consider from source S in pool L
FLOW(L,P) Flows to consider from pool L to product P
FIXQ(S,L) Range of fraction from source S in pool L has been reduced to zero by OBBT
FIXF(L,P) Range of flow from pool L to product P has been reduced to zero by OBBT
CURRENT(ND)         Current node
SRCNODE(ND)         Nodes search in spatial branch and bound
NEWNODE(ND)         New node
WAITING(ND)         Waiting node list
NODEBRN(ND,ND)      Relation between parent and child nodes
NEXTF(ND,L,P)       Next flow variable for branching;

ALIAS(IT,ITL);ALIAS(KD,KDL,KDLL);ALIAS(L,LL);ALIAS(P,PL);ALIAS(ND,NDL);

SCALAR    WHATM     What Model                    /9/;
SCALAR    WHATSBB   What OBBT in SBB              /2/;
SCALAR    TCPU      Total computational time of algorithm (s)         /3600/;
SCALAR    CPUOBBT   Remaining computational time in OBBT (s)  /300/;
SCALAR    CPURELX   Time limit in relaxation for trigerring SBB (s)   /100/;
SCALAR    TRGGAP    Optimality gap (%)     /0.0001/;
SCALAR    NDIGITS   Number of digits used in numeric representation
SCALAR    PARTVAR   Number of partitioned variables
SCALAR    DISGVAR   Number of disaggregated variables
SCALAR    BILINTR   Number of bilinear variables
SCALAR    LBOUND    Lower Bound from original problem
SCALAR    UBOUND    Upper bound from relaxation problem
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
SCALAR    STATUS    Another controller

PARAMETERS
BASIS(KD) Basis used in digit KD        /KD1*KD1 2/
SUPMAX(S) Maximum supply of source S
CAPMAX(L) Maximum size of pool L
DEMMIN(P) Minimum demand for product P
DEMMAX(P) Maximum demand for product P
PRICE(S)  Unit price of source S
VALUE(P)  Unit price of product P
CS(S,K)   Level for source S of quality K
CMIN(P,K) Lower bound for product P quality K
CMAX(P,K) Upper bound for product P quality K
LINKSL(S,L)   Connectivity source S pool L
LINKSP(S,P)   Connectivity source S product P
LINKLP(L,P)   Connectivity pool L product P
QMIN(S,L) Minimum fraction of source S in pool L
QMAX(S,L) Maximum fraction of source S in pool L
QMINJ(S,L,J)        Minimum fraction of source S in pool L number J
QMAXJ(S,L,J)        Maximum fraction of source S in pool L number J
FLPMIN(L,P)         Minimum flow from pool L to product P
FLPMAX(L,P)         Maximum flow from pool L to product P
FLPMINT(S,L,P,J,KD) Minimum flow from source S through pool L to product P when fraction is constrained to number J in digit KD
FLPMAXT(S,L,P,J,KD) Maximum flow from source S through pool L to product P when fraction is constrained to number J in digit KD
RPART(S,L)          Real partitions for the fraction from source S to pool L
CARDK(S,L)          Required number of positions for fraction of source S in pool L
KVAL(KD,S,L)        Position of significant digit KD in numerical basis representation for fraction of source S in pool L
JVAL(J,KD)          Value to use for number J in digit KD
LBF(L,P,ND)         Lower bound for FLP variable from pool L to product P in node ND
UBF(L,P,ND)         Upper bound for FLP variable from pool L to product P in node ND
LBQ(S,L,ND)         Lower bound for Q variable from source S to pool L in node ND
UBQ(S,L,ND)         Lower bound for Q variable from source S to pool L in node ND
QRNG(S,L,IT)        Range of fraction variables for source S in pool L in iteration IT of OBBT
QRR(S,L,IT)         Range reduction of fraction variables for source S in pool L in iteration IT (%)
FLPRNG(L,P,IT)      Range of flow variables from pool L to product P in iteration IT of OBBT
FLPRR(L,P,IT)       Range reduction of flow variables from pool L to pool K in iteration IT (%)
NFIXV(IT) Number of fixed variables after OBBT in iteration IT
AVERR(IT) Average range reduction in iteration IT
LBSP(SP)  Lower bound from solution SP in pool
UBSP(SP)  Upper bound from solution SP in pool
ERRF(S,L,P)         Error of bilinear term associated to flow from source S through pool L to product P
BOUND(ND) Lower bound of node ND
LOGBB(ND,*)         Logging information for branch and bound tree
LOGIT(ND,*)         Logging information for iterations
LOGND(ND,L,P,*)     Logging information for nodes FLP variables
LOGND2(ND,S,L,*)    Logging information for nodes Q variables

*Change location of Data files
$include C:\Users\Castro\Documents\GAMS files\Aulas\Data Files\PoolingA2.txt

ACTS(S)=yes$(SUPMAX(S) GT 0);
ACTL(L)=yes$(CAPMAX(L) GT 0);
ACTP(P)=yes$(DEMMAX(P) GT 0);
ACTK(K)=yes$(SMAX(P$(ACTP(P)),CMAX(P,K)) GT 0);
ACTSL(S,L)$(ACTS(S) and ACTL(L))=yes$(LINKSL(S,L) EQ 1);
ACTSP(S,P)$(ACTS(S) and ACTP(P))=yes$(LINKSP(S,P) EQ 1);
ACTLP(L,P)$(ACTL(L) and ACTP(P))=yes$(LINKLP(L,P) EQ 1);
FLPMIN(L,P)$(ACTLP(L,P))=0;
FLPMAX(L,P)$(ACTLP(L,P))=DEMMAX(P);
QMIN(S,L)$(ACTSL(S,L))=0;
QMAX(S,L)$(ACTSL(S,L))=1;
FLOW(ACTLP)=yes;FRAC(ACTSL)=yes;
FIXQ(FRAC)=no;FIXF(FLOW)=no;
QRNG(FRAC,'IT0')=QMAX(FRAC)-QMIN(FRAC);
FLPRNG(FLOW,'IT0')=FLPMAX(FLOW)-FLPMIN(FLOW);
CURRENT('ND1')=yes;NEWNODE('ND1')=yes;WAITING('ND1')=yes;
LBQ(FRAC,CURRENT)=QMIN(FRAC);UBQ(FRAC,CURRENT)=QMAX(FRAC);
LBF(FLOW,CURRENT)=FLPMIN(FLOW);UBF(FLOW,CURRENT)=FLPMAX(FLOW);
PARTVAR=card(FRAC);DISGVAR=card(FLOW);BILINTR=SUM((S,L,P)$(ACTSL(S,L) and ACTLP(L,P)),1);

Display ACTSL,ACTSP,ACTLP,ACTK,FLPMIN,FLPMAX,QMIN,QMAX;
Display PARTVAR,DISGVAR,BILINTR;

VARIABLES
PROFIT,VBOUND;

BINARY VARIABLES
Y(S,L,J,KD)         Identifies that fraction of source S in pool L has digit J in position KD

POSITIVE VARIABLES
FLP(L,P)  Flow from pool L to product P
FSP(S,P)  Flow from source S going directly to product P
Q(S,L)    Fraction of source S in pool L
ZFLP(S,L,P)         Variable associated to bilinear term involving variables FLP and Q
QR(S,L)   Residual variable related to fraction of source S in pool L
FLPQR(S,L,P)        Residual variable related to flow from source S through pool L to product P
FLPDS(S,L,P,J,KD)   Disaggregated flowrate from source S through pool L to product P associated to digit J in position KD
LQ(S,L)   Discretized variable related to fraction of source S in pool L
LQR(S,L)  Residual variable related to fraction of source S in pool L
LFLP(S,L,P)         Discretized variable related to flow from source S through pool L to product P
LFLPQR(S,L,P)       Residual variable related to flow from source S through pool L to product P

EQUATIONS
OBJ       Objective function (profit maximization)
EQ1(S)    Balance on availability of source S
EQ3(L)    Upper bound on size of pool L
EQ4(P)    Supply of product P must exceed its minimum demand
EQ5(P)    Supply of product P must not exceed its demand
EQ6(L)    Sum of source fractions in pool L must be equal to 1
EQ7(P,K)  The level for product P of quality K cannot exceed maximum value
EQ8(P,K)  The level for product P of quality K must exceed minimum value

OBJLP     Linear objective function (profit maximization)
LP1(S)    Linear balance on availability of source S
LP7(P,K)  Linear constraint for ensuring the level for product P of quality K cannot exceed maximum value
LP8(P,K)  Linear constraint for ensuring the level for product P of quality K must exceed minimum value
LP9(S,L)  Redundant constraint to improve quality of relaxation for relating bilinear variable Z with flow from source S to pool L
LP10(L,P) Redundant constraint to improve quality of relaxation for relating bilinear variable Z with flow from pool L to product P
GLMC_1(S,L,P)       Mc Cormick underestimators 1 for flow from source S through pool L to product P
GLMC_2(S,L,P)       Mc Cormick underestimators 2 for flow from source S through pool L to product P
GLMC_3(S,L,P)       Mc Cormick overestimators 1 for flow from source S through pool L to product P
GLMC_4(S,L,P)       Mc Cormick overestimators 2 for flow from source S through pool L to product P

NMDT1(S,L)      Relation between fraction of source S in pool L and its discretized variable
NMDT2(S,L,P)    Definition of bilinear term involving fraction of source S in pool L and flow from L to P
NMDT3(S,L,P)    Relation between bilinear term involving flow from source S through pool L to product P and its residual variable
NMDT4(S,L)      Relation between discretized variable of fraction of source S in pool L and its residual variable
NMDT5(S,L,P,KD) Relation between disaggregated flow from source S through pool L to product P and position KD with original FLP variable
NMDT6(S,L,KD,J) Disaggregated flow variable from source S to pool L can only take positive values for position KD value J if the binary has been selected
NMDT7(S,L,KD,J) Disaggregated flow variable from source S to pool L can only take positive values for position KD value J if the binary has been selected
NMDT8(S,L,KD)   To compute fraction of source S in pool L exactly one value must be selected for position KD
NMDT9(S,L,P)    McCormick relaxation of bilinear term involving flow from source S through pool L to product P and its residual
NMDT10(S,L,P)   McCormick relaxation of bilinear term involving flow from source S through pool L to product P and its residual
NMDT11(S,L,P)   McCormick relaxation of bilinear term involving flow from source S through pool L to product P and its residual
NMDT12(S,L,P)   McCormick relaxation of bilinear term involving flow from source S through pool L to product P and its residual

OBJOBBT   Objective function for optimality-based bound tightening
OBBTEQ1   Profit must be greater than current lower bound
;

OBJ..     PROFIT=e=SUM((L,P)$(ACTLP(L,P)),VALUE(P)*FLP(L,P))+SUM((S,P)$(ACTSP(S,P)),(VALUE(P)-PRICE(S))*FSP(S,P))-SUM((S,L,P)$(ACTSL(S,L) and ACTLP(L,P)),PRICE(S)*Q(S,L)*FLP(L,P));
EQ1(S)$(ACTS(S))..  SUM((L,P)$(ACTSL(S,L) and ACTLP(L,P)),Q(S,L)*FLP(L,P))+SUM(P$(ACTSP(S,P)),FSP(S,P))=l=SUPMAX(S);
EQ3(L)$(ACTL(L))..  SUM(P$(ACTLP(L,P)),FLP(L,P))=l=CAPMAX(L);
EQ4(P)$(ACTP(P) and DEMMIN(P) GT 0)..   SUM(L$(ACTLP(L,P)),FLP(L,P))+SUM(S$(ACTSP(S,P)),FSP(S,P))=g=DEMMIN(P);
EQ5(P)$(ACTP(P))..  SUM(L$(ACTLP(L,P)),FLP(L,P))+SUM(S$(ACTSP(S,P)),FSP(S,P))=l=DEMMAX(P);
EQ6(L)$(ACTL(L))..  SUM(S$(ACTSL(S,L)),Q(S,L))=e=1;
EQ7(P,K)$(ACTP(P) and ACTK(K))..        SUM((S,L)$(ACTSL(S,L) and ACTLP(L,P)),CS(S,K)*Q(S,L)*FLP(L,P))+SUM(S$(ACTSP(S,P)),CS(S,K)*FSP(S,P))=l=(SUM(L$(ACTLP(L,P)),FLP(L,P))+SUM(S$(ACTSP(S,P)),FSP(S,P)))*CMAX(P,K);
EQ8(P,K)$(ACTP(P) and ACTK(K) and CMIN(P,K) GT 0)..         SUM((S,L)$(ACTSL(S,L) and ACTLP(L,P)),CS(S,K)*Q(S,L)*FLP(L,P))+SUM(S$(ACTSP(S,P)),CS(S,K)*FSP(S,P))=g=(SUM(L$(ACTLP(L,P)),FLP(L,P))+SUM(S$(ACTSP(S,P)),FSP(S,P)))*CMIN(P,K);

OBJLP..   PROFIT=e=SUM((L,P)$(ACTLP(L,P)),VALUE(P)*FLP(L,P))+SUM((S,P)$(ACTSP(S,P)),(VALUE(P)-PRICE(S))*FSP(S,P))-SUM((S,L,P)$(ACTSL(S,L) and ACTLP(L,P)),PRICE(S)*ZFLP(S,L,P));
LP1(S)$(ACTS(S))..        SUM((L,P)$(ACTSL(S,L) and ACTLP(L,P)),ZFLP(S,L,P))+SUM(P$(ACTSP(S,P)),FSP(S,P))=l=SUPMAX(S);
LP7(P,K)$(ACTP(P) and ACTK(K))..        SUM((S,L)$(ACTSL(S,L) and ACTLP(L,P)),CS(S,K)*ZFLP(S,L,P))+SUM(S$(ACTSP(S,P)),CS(S,K)*FSP(S,P))=l=(SUM(L$(ACTLP(L,P)),FLP(L,P))+SUM(S$(ACTSP(S,P)),FSP(S,P)))*CMAX(P,K);
LP8(P,K)$(ACTP(P) and ACTK(K) and CMIN(P,K) GT 0)..         SUM((S,L)$(ACTSL(S,L) and ACTLP(L,P)),CS(S,K)*ZFLP(S,L,P))+SUM(S$(ACTSP(S,P)),CS(S,K)*FSP(S,P))=g=(SUM(L$(ACTLP(L,P)),FLP(L,P))+SUM(S$(ACTSP(S,P)),FSP(S,P)))*CMIN(P,K);
LP9(S,L)$(ACTSL(S,L))..       SUM(P$(ACTLP(L,P)),ZFLP(S,L,P))=l=Q(S,L)*CAPMAX(L);
LP10(L,P)$(ACTLP(L,P))..      SUM(S$(ACTSL(S,L)),ZFLP(S,L,P))=e=FLP(L,P);
GLMC_1(S,L,P)$(ACTSL(S,L) and ACTLP(L,P))..       ZFLP(S,L,P)=g=FLP(L,P)*QMIN(S,L)+FLPMIN(L,P)*Q(S,L)-FLPMIN(L,P)*QMIN(S,L);
GLMC_2(S,L,P)$(ACTSL(S,L) and ACTLP(L,P))..       ZFLP(S,L,P)=g=FLP(L,P)*QMAX(S,L)+FLPMAX(L,P)*Q(S,L)-FLPMAX(L,P)*QMAX(S,L);
GLMC_3(S,L,P)$(ACTSL(S,L) and ACTLP(L,P))..       ZFLP(S,L,P)=l=FLP(L,P)*QMIN(S,L)+FLPMAX(L,P)*Q(S,L)-FLPMAX(L,P)*QMIN(S,L);
GLMC_4(S,L,P)$(ACTSL(S,L) and ACTLP(L,P))..       ZFLP(S,L,P)=l=FLP(L,P)*QMAX(S,L)+FLPMIN(L,P)*Q(S,L)-FLPMIN(L,P)*QMAX(S,L);

NMDT1(S,L)$(ACTSL(S,L))..     Q(S,L)=e=QMIN(S,L)+LQ(S,L)*(QMAX(S,L)-QMIN(S,L));
NMDT2(S,L,P)$(ACTSL(S,L) and ACTLP(L,P))..        ZFLP(S,L,P)=e=FLP(L,P)*QMIN(S,L)+LFLP(S,L,P)*(QMAX(S,L)-QMIN(S,L));
NMDT3(S,L,P)$(ACTSL(S,L) and ACTLP(L,P))..        LFLP(S,L,P)=e=SUM((J,KD)$(ACTJ(J,KD) and ACTKD(KD,S,L)),FLPDS(S,L,P,J,KD)*JVAL(J,KD)*PROD(KDL$(ord(KDL) LE ord(KD)),1/BASIS(KDL)))+LFLPQR(S,L,P);
NMDT4(S,L)$(ACTSL(S,L))..     LQ(S,L)=e=SUM((J,KD)$(ACTJ(J,KD) and ACTKD(KD,S,L)),JVAL(J,KD)*Y(S,L,J,KD)*PROD(KDL$(ord(KDL) LE ord(KD)),1/BASIS(KDL)))+LQR(S,L);
NMDT5(S,L,P,KD)$(ACTLP(L,P) and ACTKD(KD,S,L))..  FLP(L,P)=e=SUM(J$(ACTJ(J,KD)),FLPDS(S,L,P,J,KD));
NMDT6(S,L,KD,J)$(ACTJ(J,KD) and ACTKD(KD,S,L))..  SUM(P$(ACTLP(L,P)),FLPDS(S,L,P,J,KD))=l=SUM(P$(ACTLP(L,P)),FLPMAX(L,P))*Y(S,L,J,KD);
NMDT7(S,L,KD,J)$(ACTJ(J,KD) and ACTKD(KD,S,L))..  SUM(P$(ACTLP(L,P)),FLPDS(S,L,P,J,KD))=g=SUM(P$(ACTLP(L,P)),FLPMIN(L,P))*Y(S,L,J,KD);
NMDT8(S,L,KD)$(ACTKD(KD,S,L))..         SUM(J$(ACTJ(J,KD)),Y(S,L,J,KD))=e=1;
NMDT9(S,L,P)$(ACTSL(S,L) and ACTLP(L,P))..        LFLPQR(S,L,P)=l=FLPMAX(L,P)*LQR(S,L);
NMDT10(S,L,P)$(ACTSL(S,L) and ACTLP(L,P))..       LFLPQR(S,L,P)=g=FLPMIN(L,P)*LQR(S,L);
NMDT11(S,L,P)$(ACTSL(S,L) and ACTLP(L,P))..       LFLPQR(S,L,P)=g=(FLP(L,P)-FLPMAX(L,P))*MIN(1,PROD(KD$(ACTKD(KD,S,L)),1/BASIS(KD)))+FLPMAX(L,P)*LQR(S,L);
NMDT12(S,L,P)$(ACTSL(S,L) and ACTLP(L,P))..       LFLPQR(S,L,P)=l=(FLP(L,P)-FLPMIN(L,P))*MIN(1,PROD(KD$(ACTKD(KD,S,L)),1/BASIS(KD)))+FLPMIN(L,P)*LQR(S,L);

OBJOBBT.. VBOUND=e=SUM((L,P)$(BNDF(L,P)),FLP(L,P))+SUM((S,L)$(BNDQ(S,L)),Q(S,L));
OBBTEQ1.. SUM((L,P)$(ACTLP(L,P)),VALUE(P)*FLP(L,P))+SUM((S,P)$(ACTSP(S,P)),(VALUE(P)-PRICE(S))*FSP(S,P))-SUM((S,L,P)$(ACTSL(S,L) and ACTLP(L,P)),PRICE(S)*ZFLP(S,L,P))=g=LBOUND;


Q.lo(S,L)$(ACTSL(S,L))=QMIN(S,L);
Q.up(S,L)$(ACTSL(S,L))=QMAX(S,L);
FLP.lo(L,P)$(ACTLP(L,P))=FLPMIN(L,P);
FLP.up(L,P)$(ACTLP(L,P))=FLPMAX(L,P);
LQ.up(S,L)$(ACTSL(S,L))=1;


OPTION optcr=1E-6;
OPTION limrow=0;
OPTION limcol=0;
OPTION Solprint=Off;
OPTION decimals=6;
OPTION threads=0;
OPTION reslim=3600;
OPTION QCP=GloMIQO;
OPTION NLP=CONOPT;

MODEL Poolq using /OBJ,EQ1,EQ3,EQ4,EQ5,EQ6,EQ7,EQ8/;
MODEL PoolqLP using /OBJLP,LP1,EQ3,EQ4,EQ5,EQ6,LP7,LP8,LP9,LP10,GLMC_1,GLMC_2,GLMC_3,GLMC_4/; #
MODEL MDT using /OBJLP,LP1,EQ3,EQ4,EQ5,EQ6,LP7,LP8,LP9,LP10,NMDT1,NMDT2,NMDT3,NMDT4,NMDT5,NMDT6,NMDT7,NMDT8,NMDT9,NMDT10,NMDT11,NMDT12/;
MODEL OBBT_LP using /OBJOBBT,OBBTEQ1,LP1,EQ3,EQ4,EQ5,EQ6,LP7,LP8,LP9,LP10,GLMC_1,GLMC_2,GLMC_3,GLMC_4/;
MODEL OBBT_MDT using /OBJOBBT,OBBTEQ1,LP1,EQ3,EQ4,EQ5,EQ6,LP7,LP8,LP9,LP10,NMDT1,NMDT2,NMDT3,NMDT4,NMDT5,NMDT6,NMDT7,NMDT8,NMDT9,NMDT10,NMDT11,NMDT12/;

file OptFile /C:\Users\Castro\Documents\GAMS files\GitHub\Pooling_q\Cplex.op9/; #Directory must match location of Pooling_q.gms file
file Results   /C:\Users\Castro\Documents\My Data Sources\GAMS Output\OBBTPool.txt/;
file Search   /C:\Users\Castro\Documents\My Data Sources\GAMS Output\GlobalSearchPool.txt/;
file SPoolRes /C:\Users\Castro\Documents\My Data Sources\GAMS Output\SPoolRes.txt/;
file fsoln;
Results.pw=700;Search.pw=700;SPoolRes.pw=700;
Results.pc=6;Search.pc=6;SPoolRes.pc=6; #Separa com tabs os valores das variaveis

if(WHATM EQ 7,
OPTION NLP=BARON;
SOLVE Poolq using NLP maximizing PROFIT;
LBOUND=Poolq.objval;
UBOUND=Poolq.objest;
OPTGAP=ABS((UBOUND-LBOUND)/UBOUND)*100;
Display FLP.l,FSP.l,Q.l,UBOUND,LBOUND,OPTGAP;
);
if(WHATM EQ 8,
OPTION QCP=GUROBI;
Poolq.optfile=2;
SOLVE Poolq using QCP maximizing PROFIT;
LBOUND=Poolq.objval;
UBOUND=Poolq.objest;
OPTGAP=ABS((UBOUND-LBOUND)/UBOUND)*100;
Display FLP.l,FSP.l,Q.l,UBOUND,LBOUND,OPTGAP;
);

BNDQ(FRAC)=no;BNDF(FLOW)=no;

if(WHATM EQ 9,
MDT.optcr=1E-7;
OBBT_MDT.reslim=60;
OBBT_MDT.optcr=1E-5;
put Search;
put 'CPU','LB','UB','Gap (%)','Method','Intervals'/;

#Quick computation of lower and upper bound
SOLVE PoolqLP using LP maximizing PROFIT;
MIPCPU=MAX(PoolqLP.resusd,0.001);
UBOUND=PoolqLP.objval;
SOLVE Poolq using NLP maximizing PROFIT;
          if(Poolq.modelstat EQ 2,
          LBOUND=Poolq.objval;
          else
          LBOUND=-INF;
          );
OPTGAP=ABS((UBOUND-LBOUND)/UBOUND)*100;
put TimeElapsed,LBOUND:9:3,UBOUND:9:3,OPTGAP:9:5,'Global McCormick' /;
putclose Search;Search.ap=1; #Appends text to file instead of overwriting
NBOUNDP=2*(PARTVAR+DISGVAR);
ESTOBBT=MIPCPU*NBOUNDP;TSTART=TimeElapsed;
CLOOP$(OPTGAP LE TRGGAP)=0;
#Global McCormick OBBT
          loop(IT$(CLOOP and ESTOBBT LT CPUOBBT and ord(IT) EQ CURIT+1),
                    loop(FRAC$(not FIXQ(FRAC)),
                    BNDQ(FRAC)=yes;
                    SOLVE OBBT_LP using LP minimizing VBOUND;
                    LBQ(FRAC,CURRENT)=MIN(VBOUND.l,QMAX(FRAC)); #Avoids numerical errors
                    QMIN(FRAC)=SUM(CURRENT,LBQ(FRAC,CURRENT));Q.lo(FRAC)=QMIN(FRAC);
                    SOLVE OBBT_LP using LP maximizing VBOUND;
                    UBQ(FRAC,CURRENT)=MAX(VBOUND.l,QMIN(FRAC)); #Avoids numerical errors
                    QMAX(FRAC)=SUM(CURRENT,UBQ(FRAC,CURRENT));Q.up(FRAC)=QMAX(FRAC);
                    BNDQ(FRAC)=no;
                    );
                    loop(FLOW$(not FIXF(FLOW)),
                    BNDF(FLOW)=yes;
                    SOLVE OBBT_LP using LP minimizing VBOUND;
                    LBF(FLOW,CURRENT)=MIN(VBOUND.l,FLPMAX(FLOW)); #Avoids numerical errors
                    FLPMIN(FLOW)=SUM(CURRENT,LBF(FLOW,CURRENT));FLP.lo(FLOW)=FLPMIN(FLOW);
                    SOLVE OBBT_LP using LP maximizing VBOUND;
                    UBF(FLOW,CURRENT)=MAX(VBOUND.l,FLPMIN(FLOW)); #Avoids numerical errors
                    FLPMAX(FLOW)=SUM(CURRENT,UBF(FLOW,CURRENT));FLP.up(FLOW)=FLPMAX(FLOW);
                    BNDF(FLOW)=no;
                    );
          #Report results of OBBT
          OBBTCPU=TimeElapsed-TSTART;CPUOBBT=MAX(0,CPUOBBT-OBBTCPU);
          QRNG(FRAC,IT)=QMAX(FRAC)-QMIN(FRAC);FLPRNG(FLOW,IT)=FLPMAX(FLOW)-FLPMIN(FLOW);
          QRR(FRAC,IT)=(QRNG(FRAC,'IT0')-QRNG(FRAC,IT))/QRNG(FRAC,'IT0')*100;
          FLPRR(FLOW,IT)=(FLPRNG(FLOW,'IT0')-FLPRNG(FLOW,IT))/FLPRNG(FLOW,'IT0')*100;
          AVERR(IT)=(SUM(FRAC,(QRR(FRAC,IT))$(not FIXQ(FRAC))+100$(FIXQ(FRAC)))+SUM(FLOW,(FLPRR(FLOW,IT))$(not FIXF(FLOW))+100$(FIXF(FLOW))))/(PARTVAR+DISGVAR);
          FIXQ(FRAC)=yes$(QRR(FRAC,IT) GT 99.9999);
          FIXF(FLOW)=yes$(FLPRR(FLOW,IT) GT 99.9999);
          NFIXV(IT)=SUM(FRAC$(FIXQ(FRAC)),1)+SUM(FLOW$(FIXF(FLOW)),1);
          NBOUNDP=2*(PARTVAR+DISGVAR-card(FIXQ)-card(FIXF));
          #Recompute Bounds
          SOLVE PoolqLP using LP maximizing PROFIT;
          MIPCPU=MAX(PoolqLP.resusd,0.001);
          UBOUND=PoolqLP.objval;
          SOLVE Poolq using NLP maximizing PROFIT;
                    if(Poolq.modelstat EQ 2,
                              if(Poolq.objval GT LBOUND*1.0001,
                              LBOUND=Poolq.objval;
                              CURIT=CURIT+1;
                              );
                    );
          OPTGAP=ABS((UBOUND-LBOUND)/UBOUND)*100;
          put TimeElapsed,LBOUND:9:3,UBOUND:9:3,OPTGAP:9:5,'Global McCormick after LP-based OBBT' /;
          putclose Search;
          );
CLOOP$(OPTGAP LE TRGGAP or TimeElapsed GT TCPU)=0;
BASIS(KD)=0;KVAL(KD,S,L)=0;ACTJ(J,KD)=no;JVAL(J,KD)=0;
#Relaxation with mixed-radix MDT
          loop(KDLL$(CLOOP),
          BASIS(KDLL)=2;
          NDIGITS=SUM(KD$(BASIS(KD) GT 0),1);
          CARDK(S,L)$(ACTSL(S,L))=NDIGITS;
          ACTKD(KD,S,L)$(ACTSL(S,L))=yes$(ord(KD) LE CARDK(S,L));
          KVAL(KD,S,L)$(ACTKD(KD,S,L))=-BASIS(KD);
          ACTJ(J,KD)$(ord(KD) LE NDIGITS)=yes$(ord(J) LE BASIS(KD));
          JVAL(J,KD)$(ACTJ(J,KD))=ord(J)-1;
          LQR.up(S,L)$(ACTSL(S,L))=MIN(1,PROD(KD$(ACTKD(KD,S,L)),1/BASIS(KD)));
          put OptFile;
          put 'tilim ', (TCPU-TimeElapsed):<9:0 /;
          put 'Solnpool solnpool.gdx' /;
          put 'SolnPoolPop 1' /;
          put 'SolnPoolIntensity 1' /;
          put 'SolnPoolCapacity 10' /;
          put 'SolnPoolGap 0.01' /;
          put 'SolnpoolReplace 2' /;
          putclose OptFile;
          MDT.optfile=9;
          SOLVE MDT using MIP maximizing PROFIT;
          ESTRELX=sqr(MDT.resusd)/MIPCPU;
          MIPCPU=MDT.resusd;
          ESTOBBT=MIPCPU*NBOUNDP;
                    if(MDT.modelstat NE 10 and MDT.modelstat NE 14,
                    UBOUND=MIN(UBOUND,MDT.objest);
                    OPTGAP=ABS((UBOUND-LBOUND)/UBOUND)*100;
                    );
                    if(OPTGAP GT TRGGAP,
                    #Solution pool to generate multiple starting points for NLP
                    put SPoolRes;
                    put 'Solutions from CPLEX pool' /;
                    put 'File','UB','LB' /;
                    execute_load 'solnpool.gdx', SPOOL=Index; CARDSP=card(SPOOL);
                              loop(SPOOL(SP),
                              put_utility fsoln 'gdxin' / SPOOL.te(SP):0:0;
                              execute_loadpoint;
                              UBSP(SP)=PROFIT.l;
                              put SPoolRes SP.te(SP),UBSP(SP):9:5;
                              SOLVE Poolq using NLP maximizing PROFIT;
                                        if(Poolq.modelstat EQ 2,
                                        LBSP(SP)=PROFIT.l; put LBSP(SP):9:5 /;
                                        LBOUND=MAX(LBOUND,Poolq.objval);
                                        else
                                        put 'Infeasible' /;
                                        );
                              );
                    putclose SPoolRes;
                    OPTGAP=ABS((UBOUND-LBOUND)/UBOUND)*100;
                    );
          put Search;
          put TimeElapsed,LBOUND:9:3,UBOUND:9:3,OPTGAP:9:5,'NMDT Relaxation',(2**NDIGITS):9:0 /;
          putclose Search;
          CLOOP$(OPTGAP LE TRGGAP or TimeElapsed GT TCPU or MIPCPU GT CPURELX)=0;
                    if(CLOOP and ESTOBBT LT CPUOBBT and ESTRELX*NBOUNDP GT CPUOBBT/2,
                    #NMDT OBBT with current settings
                    TSTART=TimeElapsed;
                              loop(FRAC$(not FIXQ(FRAC)),
                              BNDQ(FRAC)=yes;
                              SOLVE OBBT_MDT using MIP minimizing VBOUND;
                                        if(OBBT_MDT.modelstat NE 10 and OBBT_MDT.modelstat NE 14 and OBBT_MDT.solvestat NE 4,
                                        LBQ(FRAC,CURRENT)=MIN(OBBT_MDT.objest,QMAX(FRAC)); #Avoids numerical errors
                                        else
                                        SOLVE OBBT_LP using LP minimizing VBOUND;
                                                  if(OBBT_LP.modelstat NE 4,
                                                  LBQ(FRAC,CURRENT)=MIN(VBOUND.l,QMAX(FRAC));
                                                  );
                                        );
                              QMIN(FRAC)=SUM(CURRENT,LBQ(FRAC,CURRENT));Q.lo(FRAC)=QMIN(FRAC);
                              SOLVE OBBT_MDT using MIP maximizing VBOUND;
                                        if(OBBT_MDT.modelstat NE 10 and OBBT_MDT.modelstat NE 14 and OBBT_MDT.solvestat NE 4,
                                        UBQ(FRAC,CURRENT)=MAX(OBBT_MDT.objest,QMIN(FRAC)); #Avoids numerical errors
                                        else
                                        SOLVE OBBT_LP using LP maximizing VBOUND;
                                                  if(OBBT_LP.modelstat NE 4,
                                                  UBQ(FRAC,CURRENT)=MAX(VBOUND.l,QMIN(FRAC));
                                                  );
                                        );
                              QMAX(FRAC)=SUM(CURRENT,UBQ(FRAC,CURRENT));Q.up(FRAC)=QMAX(FRAC);
                              BNDQ(FRAC)=no;
                              );
                              loop(FLOW$(not FIXF(FLOW)),
                              BNDF(FLOW)=yes;
                              SOLVE OBBT_MDT using MIP minimizing VBOUND;
                                        if(OBBT_MDT.modelstat NE 10 and OBBT_MDT.modelstat NE 14 and OBBT_MDT.solvestat NE 4,
                                        LBF(FLOW,CURRENT)=MIN(OBBT_MDT.objest,FLPMAX(FLOW)); #Avoids numerical errors
                                        else
                                        SOLVE OBBT_LP using LP minimizing VBOUND;
                                                  if(OBBT_LP.modelstat NE 4,
                                                  LBF(FLOW,CURRENT)=MIN(VBOUND.l,FLPMAX(FLOW));
                                                  );
                                        );
                              FLPMIN(FLOW)=SUM(CURRENT,LBF(FLOW,CURRENT));FLP.lo(FLOW)=FLPMIN(FLOW);
                              SOLVE OBBT_MDT using MIP maximizing VBOUND;
                                        if(OBBT_MDT.modelstat NE 10 and OBBT_MDT.modelstat NE 14 and OBBT_MDT.solvestat NE 4,
                                        UBF(FLOW,CURRENT)=MAX(OBBT_MDT.objest,FLPMIN(FLOW)); #Avoids numerical errors
                                        else
                                                  if(OBBT_LP.modelstat NE 4,
                                                  SOLVE OBBT_LP using LP maximizing VBOUND;
                                                  UBF(FLOW,CURRENT)=MAX(VBOUND.l,FLPMIN(FLOW));
                                                  );
                                        );
                              FLPMAX(FLOW)=SUM(CURRENT,UBF(FLOW,CURRENT));FLP.up(FLOW)=FLPMAX(FLOW);
                              BNDF(FLOW)=no;
                              );
                    #Report results of OBBT
                    OBBTCPU=TimeElapsed-TSTART;CPUOBBT=MAX(0,CPUOBBT-OBBTCPU);
                    CURIT=CURIT+1;
                    QRNG(FRAC,IT)$(ord(IT) EQ CURIT+1)=QMAX(FRAC)-QMIN(FRAC);FLPRNG(FLOW,IT)$(ord(IT) EQ CURIT+1)=FLPMAX(FLOW)-FLPMIN(FLOW);
                    QRR(FRAC,IT)$(ord(IT) EQ CURIT+1)=(QRNG(FRAC,'IT0')-QRNG(FRAC,IT))/QRNG(FRAC,'IT0')*100;
                    FLPRR(FLOW,IT)$(ord(IT) EQ CURIT+1)=(FLPRNG(FLOW,'IT0')-FLPRNG(FLOW,IT))/FLPRNG(FLOW,'IT0')*100;
                    AVERR(IT)$(ord(IT) EQ CURIT+1)=(SUM(FRAC,(QRR(FRAC,IT))$(not FIXQ(FRAC))+100$(FIXQ(FRAC)))+SUM(FLOW,(FLPRR(FLOW,IT))$(not FIXF(FLOW))+100$(FIXF(FLOW))))/(PARTVAR+DISGVAR);
                    FIXQ(FRAC)=yes$(SUM(IT$(ord(IT) EQ CURIT+1),QRR(FRAC,IT)) GT 99.9999);
                    FIXF(FLOW)=yes$(SUM(IT$(ord(IT) EQ CURIT+1),FLPRR(FLOW,IT)) GT 99.9999);
                    NFIXV(IT)$(ord(IT) EQ CURIT+1)=SUM(FRAC$(FIXQ(FRAC)),1)+SUM(FLOW$(FIXF(FLOW)),1);
                    NBOUNDP=2*(PARTVAR+DISGVAR-card(FIXQ)-card(FIXF));
                    #Check if with results from McCormick relaxation we can already stop
                    SOLVE PoolqLP using LP maximizing PROFIT;
                    UBOUND=MIN(UBOUND,PoolqLP.objval);
                    OPTGAP=ABS((UBOUND-LBOUND)/UBOUND)*100;
                    CLOOP$(OPTGAP LE TRGGAP)=0;
                    put TimeElapsed,LBOUND:9:3,UBOUND:9:3,OPTGAP:9:5,'NMDT OBBT',(2**NDIGITS):9:0 /;
                    putclose Search;
                    SBBOBBT=ord(KDLL);
                    );
          SBBRELX=ord(KDLL);
          );
WHATSBB$(SBBOBBT EQ 0)=1;
#Spatial Branch and Bound
          if(TimeElapsed LT TCPU and OPTGAP GT TRGGAP,
          put 'Spatial B&B' /;
          put 'CPU','LB','UB','Gap (%)','Expl.','Left'/;
          #Compute the error that will be used to select branching variavel
          ERRF(S,L,P)$(ACTSL(S,L) and ACTLP(L,P))=ABS(LFLP.l(S,L,P)-FLP.l(L,P)*LQ.l(S,L));
          BOUND('ND1')=UBOUND;
          SRCNODE('ND1')=yes;
          LOGBB('ND1','Time (s)')=TimeElapsed;LOGBB('ND1','Explored')=card(SRCNODE);LOGBB('ND1','Remaining')=card(WAITING);
          LOGBB('ND1','UBOUND')=UBOUND;LOGBB('ND1','LBOUND')=LBOUND;LOGBB('ND1','OPTGAP (%)')=OPTGAP;
          LOGND(CURRENT,FLOW,'XLO')=LBF(FLOW,CURRENT);LOGND(CURRENT,FLOW,'XUP')=UBF(FLOW,CURRENT);
          LOGND2(CURRENT,FRAC,'XLO')=LBQ(FRAC,CURRENT);LOGND2(CURRENT,FRAC,'XUP')=UBQ(FRAC,CURRENT);
          put LOGBB('ND1','Time (s)'):9:2, put LOGBB('ND1','LBOUND'):9:3, put LOGBB('ND1','UBOUND'):9:3, put LOGBB('ND1','OPTGAP (%)'):9:5, put LOGBB('ND1','Explored'):9:0 ; put LOGBB('ND1','Remaining'):9:0 /;
          putclose Search;
          #Choose branching variavel
                    loop((L,P)$(ACTLP(L,P) and SMAX(S$(ACTSL(S,L)),ERRF(S,L,P)) EQ SMAX((S,LL,PL)$(ACTSL(S,LL) and ACTLP(LL,PL)),ERRF(S,LL,PL))),
                    NEXTF(CURRENT,L,P)=yes;
                    );
          loop(CURRENT,BRCHF(FLOW)=NEXTF(CURRENT,FLOW););
          LOGND(CURRENT,BRCHF,'LBOUND')=LBOUND;LOGND(CURRENT,BRCHF,'UBOUND')=UBOUND;LOGND(CURRENT,BRCHF,'Next Pool')=SUM((L,P)$(NEXTF(CURRENT,L,P)),ord(L));LOGND(CURRENT,BRCHF,'Next Prod')=SUM((L,P)$(NEXTF(CURRENT,L,P)),ord(P));
          #Perform branching
                    loop(ND$(not DONE), # and ord(ND) LE 1
                    #Node selection
                    CURRENT(NDL)=no;
                    CURRENT(WAITING(NDL))$(BOUND(NDL) EQ UBOUND)=yes;
                    FIRST=1;
                              loop(CURRENT$FIRST, #Select only one node for branching
                              WAITING(CURRENT)=no;
                              loop(NDL$(CURRENT(NDL)),BRCHF(L,P)=NEXTF(NDL,L,P);); #Assign next variavel to branch
                              LOGIT(ND,'Node_Sel.')=SUM(NDL$(CURRENT(NDL)),ord(NDL));
                              LOGIT(ND,'From_Pool')=SUM((L,P)$(BRCHF(L,P)),ord(L));
                              LOGIT(ND,'To_Prod')=SUM((L,P)$(BRCHF(L,P)),ord(P));
                              #Branch
                                        for(BISECT=1 to 2,
                                        NEWNODE(NDL)=NEWNODE(NDL-1);
                                        NODEBRN(CURRENT,NEWNODE)=yes;
                                        SRCNODE(NEWNODE)=yes;
                                        WAITING(NEWNODE)=yes;
                                        LBF(FLOW,NEWNODE)=LBF(FLOW,CURRENT);UBF(FLOW,NEWNODE)=UBF(FLOW,CURRENT);
                                        LBQ(FRAC,NEWNODE)=LBQ(FRAC,CURRENT);UBQ(FRAC,NEWNODE)=UBQ(FRAC,CURRENT);
                                                  if(BISECT EQ 1, #Left branch
                                                  LBF(BRCHF,NEWNODE)=LBF(BRCHF,CURRENT);
                                                  UBF(BRCHF,NEWNODE)=LBF(BRCHF,CURRENT)+(UBF(BRCHF,CURRENT)-LBF(BRCHF,CURRENT))/2;
                                                  else
                                                  LBF(BRCHF,NEWNODE)=LBF(BRCHF,CURRENT)+(UBF(BRCHF,CURRENT)-LBF(BRCHF,CURRENT))/2;
                                                  UBF(BRCHF,NEWNODE)=UBF(BRCHF,CURRENT);
                                                  );
                                        FLPMIN(FLOW)=SUM(NDL$(NEWNODE(NDL)),LBF(FLOW,NDL));FLPMAX(FLOW)=SUM(NDL$(NEWNODE(NDL)),UBF(FLOW,NDL));
                                        QMIN(FRAC)=SUM(NDL$(NEWNODE(NDL)),LBQ(FRAC,NDL));QMAX(FRAC)=SUM(NDL$(NEWNODE(NDL)),UBQ(FRAC,NDL));
                                        FLP.lo(FLOW)=FLPMIN(FLOW);FLP.up(FLOW)=FLPMAX(FLOW);
                                        Q.lo(FRAC)=QMIN(FRAC);Q.up(FRAC)=QMAX(FRAC);
                                        LOGND(NEWNODE,FLOW,'XLO')=FLP.lo(FLOW);LOGND(NEWNODE,FLOW,'XUP')=FLP.up(FLOW);
                                        LOGND2(NEWNODE,FRAC,'XLO')=Q.lo(FRAC);LOGND2(NEWNODE,FRAC,'XUP')=Q.up(FRAC);
                                        #OBBT
                                        STATUS=0;
                                                  loop(FRAC$(not STATUS and not FIXQ(FRAC)),
                                                  BNDQ(FRAC)=yes;
                                                            if(WHATSBB EQ 1,
                                                            #LP-based OBBT
                                                            SOLVE OBBT_LP using LP minimizing VBOUND;
                                                                      if(OBBT_LP.modelstat EQ 4, #Branch can be fathomed
                                                                      STATUS=1;
                                                                      WAITING(NEWNODE)=no;
                                                                      else
                                                                      LBQ(FRAC,NEWNODE)=MIN(VBOUND.l,QMAX(FRAC)); #Avoids numerical errors
                                                                      QMIN(FRAC)=SUM(NEWNODE,LBQ(FRAC,NEWNODE));Q.lo(FRAC)=QMIN(FRAC);
                                                                      SOLVE OBBT_LP using LP maximizing VBOUND;
                                                                      UBQ(FRAC,NEWNODE)=MAX(VBOUND.l,QMIN(FRAC)); #Avoids numerical errors
                                                                      QMAX(FRAC)=SUM(NEWNODE,UBQ(FRAC,NEWNODE));Q.up(FRAC)=QMAX(FRAC);
                                                                      LOGND2(NEWNODE,FRAC,'XLO')=Q.lo(FRAC);LOGND2(NEWNODE,FRAC,'XUP')=Q.up(FRAC);
                                                                      );
                                                            else
                                                            #MIP-based OBBT
                                                            KVAL(KD,S,L)=0;ACTJ(J,KD)=no;JVAL(J,KD)=0;
                                                            BASIS(KD)=2$(ord(KD) LE SBBOBBT);
                                                            NDIGITS=SUM(KD$(BASIS(KD) GT 0),1);
                                                            CARDK(S,L)$(ACTSL(S,L))=NDIGITS;
                                                            ACTKD(KD,S,L)$(ACTSL(S,L))=yes$(ord(KD) LE CARDK(S,L));
                                                            KVAL(KD,S,L)$(ACTKD(KD,S,L))=-BASIS(KD);
                                                            ACTJ(J,KD)$(ord(KD) LE NDIGITS)=yes$(ord(J) LE BASIS(KD));
                                                            JVAL(J,KD)$(ACTJ(J,KD))=ord(J)-1;
                                                            LQR.up(S,L)$(ACTSL(S,L))=MIN(1,PROD(KD$(ACTKD(KD,S,L)),1/BASIS(KD)));
                                                            SOLVE OBBT_MDT using MIP minimizing VBOUND;
                                                                      if(OBBT_MDT.modelstat NE 10 and OBBT_MDT.modelstat NE 14 and OBBT_MDT.solvestat NE 4,
                                                                      LBQ(FRAC,NEWNODE)=MIN(OBBT_MDT.objest,QMAX(FRAC)); #Avoids numerical errors
                                                                      QMIN(FRAC)=SUM(NEWNODE,LBQ(FRAC,NEWNODE));Q.lo(FRAC)=QMIN(FRAC);
                                                                      else
                                                                      SOLVE OBBT_LP using LP minimizing VBOUND;
                                                                                if(OBBT_LP.modelstat EQ 4, #Branch can be fathomed
                                                                                STATUS=1;
                                                                                WAITING(NEWNODE)=no;
                                                                                else
                                                                                LBQ(FRAC,NEWNODE)=MIN(VBOUND.l,QMAX(FRAC)); #Avoids numerical errors
                                                                                QMIN(FRAC)=SUM(NEWNODE,LBQ(FRAC,NEWNODE));Q.lo(FRAC)=QMIN(FRAC);
                                                                                );
                                                                      );
                                                                      if(not STATUS,
                                                                      SOLVE OBBT_MDT using MIP maximizing VBOUND;
                                                                                if(OBBT_MDT.modelstat NE 10 and OBBT_MDT.modelstat NE 14 and OBBT_MDT.solvestat NE 4,
                                                                                UBQ(FRAC,NEWNODE)=MAX(OBBT_MDT.objest,QMIN(FRAC)); #Avoids numerical errors
                                                                                QMAX(FRAC)=SUM(NEWNODE,UBQ(FRAC,NEWNODE));Q.up(FRAC)=QMAX(FRAC);
                                                                                LOGND2(NEWNODE,FRAC,'XLO')=Q.lo(FRAC);LOGND2(NEWNODE,FRAC,'XUP')=Q.up(FRAC);
                                                                                else
                                                                                SOLVE OBBT_LP using LP maximizing VBOUND;
                                                                                          if(OBBT_LP.modelstat EQ 4, #Branch can be fathomed
                                                                                          STATUS=1;
                                                                                          WAITING(NEWNODE)=no;
                                                                                          else
                                                                                          UBQ(FRAC,NEWNODE)=MAX(VBOUND.l,QMIN(FRAC)); #Avoids numerical errors
                                                                                          QMAX(FRAC)=SUM(NEWNODE,UBQ(FRAC,NEWNODE));Q.up(FRAC)=QMAX(FRAC);
                                                                                          LOGND2(NEWNODE,FRAC,'XLO')=Q.lo(FRAC);LOGND2(NEWNODE,FRAC,'XUP')=Q.up(FRAC);
                                                                                          );
                                                                                );
                                                                      );
                                                            );
                                                  BNDQ(FRAC)=no;
                                                  );
                                                  loop((LL,PL)$(FLOW(LL,PL) and not BRCHF(LL,PL) and not STATUS and not FIXF(LL,PL)),
                                                  BNDF(LL,PL)=yes;
                                                            if(WHATSBB EQ 1,
                                                            #LP-based OBBT
                                                            SOLVE OBBT_LP using LP minimizing VBOUND;
                                                                      if(OBBT_LP.modelstat EQ 4, #Branch can be fathomed
                                                                      STATUS=1;
                                                                      WAITING(NEWNODE)=no;
                                                                      else
                                                                      LBF(LL,PL,NEWNODE)=MIN(VBOUND.l,FLPMAX(LL,PL)); #Avoids numerical errors
                                                                      FLPMIN(LL,PL)=SUM(NEWNODE,LBF(LL,PL,NEWNODE));FLP.lo(LL,PL)=FLPMIN(LL,PL);
                                                                      SOLVE OBBT_LP using LP maximizing VBOUND;
                                                                      UBF(LL,PL,NEWNODE)=MAX(VBOUND.l,FLPMIN(LL,PL)); #Avoids numerical errors
                                                                      FLPMAX(LL,PL)=SUM(NEWNODE,UBF(LL,PL,NEWNODE));FLP.up(LL,PL)=FLPMAX(LL,PL);
                                                                      LOGND(NEWNODE,LL,PL,'XLO')=FLP.lo(LL,PL);LOGND(NEWNODE,LL,PL,'XUP')=FLP.up(LL,PL);
                                                                      );
                                                            else
                                                            #MIP-based OBBT
                                                            SOLVE OBBT_MDT using MIP minimizing VBOUND;
                                                                      if(OBBT_MDT.modelstat NE 10 and OBBT_MDT.modelstat NE 14 and OBBT_MDT.solvestat NE 4,
                                                                      LBF(LL,PL,NEWNODE)=MIN(OBBT_MDT.objest,FLPMAX(LL,PL));
                                                                      FLPMIN(LL,PL)=SUM(NEWNODE,LBF(LL,PL,NEWNODE));FLP.lo(LL,PL)=FLPMIN(LL,PL);
                                                                      else
                                                                      SOLVE OBBT_LP using LP minimizing VBOUND;
                                                                                if(OBBT_LP.modelstat EQ 4, #Branch can be fathomed
                                                                                STATUS=1;
                                                                                WAITING(NEWNODE)=no;
                                                                                else
                                                                                LBF(LL,PL,NEWNODE)=MIN(VBOUND.l,FLPMAX(LL,PL)); #Avoids numerical errors
                                                                                FLPMIN(LL,PL)=SUM(NEWNODE,LBF(LL,PL,NEWNODE));FLP.lo(LL,PL)=FLPMIN(LL,PL);
                                                                                );
                                                                      );
                                                                      if(not STATUS,
                                                                      SOLVE OBBT_MDT using MIP maximizing VBOUND;
                                                                                if(OBBT_MDT.modelstat NE 10 and OBBT_MDT.modelstat NE 14 and OBBT_MDT.solvestat NE 4,
                                                                                UBF(LL,PL,NEWNODE)=MAX(OBBT_MDT.objest,FLPMIN(LL,PL));
                                                                                FLPMAX(LL,PL)=SUM(NEWNODE,UBF(LL,PL,NEWNODE));FLP.up(LL,PL)=FLPMAX(LL,PL);
                                                                                LOGND(NEWNODE,LL,PL,'XLO')=FLP.lo(LL,PL);LOGND(NEWNODE,LL,PL,'XUP')=FLP.up(LL,PL);
                                                                                else
                                                                                SOLVE OBBT_LP using LP maximizing VBOUND;
                                                                                          if(OBBT_LP.modelstat EQ 4, #Branch can be fathomed
                                                                                          STATUS=1;
                                                                                          WAITING(NEWNODE)=no;
                                                                                          else
                                                                                          UBF(LL,PL,NEWNODE)=MAX(VBOUND.l,FLPMIN(LL,PL)); #Avoids numerical errors
                                                                                          FLPMAX(LL,PL)=SUM(NEWNODE,UBF(LL,PL,NEWNODE));FLP.up(LL,PL)=FLPMAX(LL,PL);
                                                                                          LOGND(NEWNODE,LL,PL,'XLO')=FLP.lo(LL,PL);LOGND(NEWNODE,LL,PL,'XUP')=FLP.up(LL,PL);
                                                                                          );
                                                                                );
                                                                      );
                                                            );
                                                  BNDF(LL,PL)=no;
                                                  );
                                        #MDT Relaxation
                                                  if(not STATUS,
                                                  KVAL(KD,S,L)=0;ACTJ(J,KD)=no;JVAL(J,KD)=0;
                                                  BASIS(KD)=2$(ord(KD) LE SBBRELX);
                                                  NDIGITS=SUM(KD$(BASIS(KD) GT 0),1);
                                                  CARDK(S,L)$(ACTSL(S,L))=NDIGITS;
                                                  ACTKD(KD,S,L)$(ACTSL(S,L))=yes$(ord(KD) LE CARDK(S,L));
                                                  KVAL(KD,S,L)$(ACTKD(KD,S,L))=-BASIS(KD);
                                                  ACTJ(J,KD)$(ord(KD) LE NDIGITS)=yes$(ord(J) LE BASIS(KD));
                                                  JVAL(J,KD)$(ACTJ(J,KD))=ord(J)-1;
                                                  LQR.up(S,L)$(ACTSL(S,L))=MIN(1,PROD(KD$(ACTKD(KD,S,L)),1/BASIS(KD)));
                                                  #No need for solution pool at this stage
                                                  put OptFile;
                                                  put 'tilim ', (MAX(0,TCPU-TimeElapsed)):<9:0 /;
                                                  putclose OptFile;
                                                  MDT.cutoff=LBOUND;
                                                  SOLVE MDT using MIP maximizing PROFIT;
                                                            if(MDT.modelstat EQ 10,BOUND(NEWNODE)=-INF;
                                                            else BOUND(NEWNODE)=MIN(MDT.objest,BOUND(CURRENT));
                                                                      if(MDT.modelstat NE 14,
                                                                      ERRF(S,L,P)$(ACTSL(S,L) and ACTLP(L,P))=ABS(LFLP.l(S,L,P)-FLP.l(L,P)*LQ.l(S,L));
                                                                      else  #Use LP relaxation to compute the error that will be used to select branching variable
                                                                      SOLVE PoolqLP using LP maximizing PROFIT;
                                                                      ERRF(S,L,P)$(ACTSL(S,L) and ACTLP(L,P))=ABS(ZFLP.l(S,L,P)-FLP.l(L,P)*Q.l(S,L));
                                                                      );
                                                            );
                                                  LOGND(NEWNODE,BRCHF,'UBOUND')=BOUND(NEWNODE);
                                                            if(SUM(NEWNODE,BOUND(NEWNODE)) LT LBOUND/(1-TRGGAP*0.01),
                                                            WAITING(NEWNODE)=no; #Fathom node
                                                            else
                                                            #Choose branching variavel
                                                                      loop((L,P)$(ACTLP(L,P) and SMAX(S$(ACTSL(S,L)),ERRF(S,L,P)) EQ SMAX((S,LL,PL)$(ACTSL(S,LL) and ACTLP(LL,PL)),ERRF(S,LL,PL))),
                                                                      NEXTF(NEWNODE,L,P)=yes;
                                                                      );
                                                            LOGND(NEWNODE,BRCHF,'Next Pool')=SUM((L,P)$(NEXTF(NEWNODE,L,P)),ord(L));
                                                            LOGND(NEWNODE,BRCHF,'Next Prod')=SUM((L,P)$(NEXTF(NEWNODE,L,P)),ord(P));
                                                            FLP.lo(FLOW)=LBF(FLOW,'ND1');FLP.up(FLOW)=UBF(FLOW,'ND1');
                                                            SOLVE Poolq using NLP maximizing PROFIT;
                                                                      if(Poolq.modelstat EQ 2,
                                                                      LOGND(NEWNODE,BRCHF,'LBOUND')=Poolq.objval;
                                                                                if(Poolq.objval GT LBOUND,
                                                                                LBOUND=Poolq.objval;
                                                                                          loop(NDL$(WAITING(NDL)), #Eliminate nodes with upper bound lower than improved lower bound
                                                                                                    if(BOUND(NDL) LT LBOUND/(1-TRGGAP*0.01),WAITING(NDL)=no;);
                                                                                          );
                                                                                );
                                                                      else
                                                                      LOGND(NEWNODE,BRCHF,'LBOUND')=LBOUND;
                                                                      );
                                                            );
                                                  else
                                                  LOGND(NEWNODE,BRCHF,'UBOUND')=-INF;
                                                  );
                                        );
                              FIRST=0;
                              );
                    UBOUND=MAX(SMAX(WAITING(NDL),BOUND(NDL)),LBOUND/(1-TRGGAP*0.01));
                    OPTGAP=ABS((UBOUND-LBOUND)/UBOUND)*100;
                    LOGBB(ND+1,'Time (s)')=TimeElapsed;LOGBB(ND+1,'Explored')=card(SRCNODE);LOGBB(ND+1,'Remaining')=card(WAITING);
                    LOGBB(ND+1,'UBOUND')=UBOUND;LOGBB(ND+1,'LBOUND')=LBOUND;LOGBB(ND+1,'OPTGAP (%)')=OPTGAP;
                    LOGIT(ND,'LBOUND')=LBOUND;
                    LOGIT(ND,'UBOUND')=UBOUND;
                    LOGIT(ND,'OPTGAP (%)')=OPTGAP;
                    LOGIT(ND,'WAITING')=card(WAITING);
                    LOGIT(ND,'CPUs')=TimeElapsed;
                    put Search;
                    put LOGBB(ND+1,'Time (s)'):9:2, put LOGBB(ND+1,'LBOUND'):9:3, put LOGBB(ND+1,'UBOUND'):9:3, put LOGBB(ND+1,'OPTGAP (%)'):9:5, put LOGBB(ND+1,'Explored'):9:0, put LOGBB(ND+1,'Remaining'):9:0 /;
                    putclose Search;
                    DONE$(card(WAITING) EQ 0 or OPTGAP LE TRGGAP or TimeElapsed GE TCPU)=1;
                    );
          display LOGBB,LOGND,LOGND2,LOGIT,NODEBRN,WAITING;
          );
put Results;
put 'OBBT';
put 'Var.', 'RangeR (%)'/;
          loop(FRAC,
          put FRAC.te(FRAC);
                    loop(IT$(ord(IT) GT 1 and ord(IT) LE CURIT+1),
                    put QRR(FRAC,IT):8:4;
                    );
          put/;
          );
          loop(FLOW,
          put FLOW.te(FLOW);
                    loop(IT$(ord(IT) GT 1 and ord(IT) LE CURIT+1),
                    put FLPRR(FLOW,IT):8:4;
                    );
          put/;
          );
put 'Fixed V',loop(IT$(ord(IT) GT 1 and ord(IT) LE CURIT+1), put NFIXV(IT);); put/;
put 'ARR',loop(IT$(ord(IT) GT 1 and ord(IT) LE CURIT+1), put AVERR(IT);); put/;
display NBOUNDP,UBOUND,LBOUND,OPTGAP;
);