$Title Global optimization of QCPs arising from the pooling problem (tp-formulation)

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
S         Sources   /S1*S60/
L         Pools     /L1*L30/
P         Products  /P1*P40/
K         Qualities /K1*K20/
J         Numbers   /0*29/
KD        Digits   /KD1*KD21/
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
ACTKD(KD,L,P)       Active significant digits KD for fraction from pool L to product P
BNDX(L,P) Bound contract is done for fraction from pool L to product P
BNDF(S,L) Bound contracting is done for flow from source S to pool L
BRCHF(S,L)          Flowrate from source S to pool L being branched
FLOW(S,L) Flows to consider from source S to pool L
FRAC(L,P) Fractions to consider from pool L to product P
FIXF(S,L) Range of flow from source S to pool L has been reduced to zero by OBBT
FIXX(L,P) Range of fraction from pool L to product P has been reduced to zero by OBBT
CURRENT(ND)         Current node
SRCNODE(ND)         Nodes search in spatial branch and bound
NEWNODE(ND)         New node
WAITING(ND)         Waiting node list
NODEBRN(ND,ND)      Relation between parent and child nodes
NEXTF(ND,S,L)       Next flow variable for branching;

ALIAS(KD,KDL,KDLL);;ALIAS(S,SL);ALIAS(L,LL);ALIAS(ND,NDL);

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
SCALAR    STATUS    Another controller;

PARAMETERS
BASIS(KD) Basis used in digit KD        /KD1*KD3 2/
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
XMIN(L,P) Minimum fraction from pool L to product P
XMAX(L,P) Maximum fraction from pool L to product P
FSLMIN(S,L)         Minimum flow from source S to pool L
FSLMAX(S,L)         Maximum flow from source S to pool L
CARDK(L,P)          Required number of positions for fraction from pool L to product P
KVAL(KD,L,P)        Position of significant digit KD in numerical basis representation for fraction from pool L to product P
JVAL(J,KD)          Value to use for number J in digit KD
LBF(S,L,ND)         Lower bound for FSL variable from source S to pool L in node ND
UBF(S,L,ND)         Upper bound for FSL variable from source S to pool L in node ND
LBX(L,P,ND)         Lower bound for X variable from pool L to product P in node ND
UBX(L,P,ND)         Lower bound for X variable from pool L to product P in node ND
XRNG(L,P,IT)        Range of fraction variables for pool L to product P in iteration IT of OBBT
XRR(L,P,IT)         Range reduction of fraction variables for pool L to product P in iteration IT (%)
FSLRNG(S,L,IT)      Range of flow variables from source S to pool L in iteration IT of OBBT
FSLRR(S,L,IT)       Range reduction of flow variables from source S to pool L in iteration IT (%)
NFIXV(IT) Number of fixed variables after OBBT in iteration IT
AVERR(IT) Average range reduction in iteration IT
LBSP(SP)  Lower bound from solution SP in pool
UBSP(SP)  Upper bound from solution SP in pool
ERRF(S,L,P)         Error of bilinear term associated to flow from source S through pool L to product P
BOUND(ND) Lower bound of node ND
LOGBB(ND,*)         Logging information for branch and bound tree
LOGIT(ND,*)         Logging information for iterations
LOGND(ND,S,L,*)     Logging information for nodes FSL variables
LOGND2(ND,L,P,*)    Logging information for nodes X variables;

*Change location of Data files
$include C:\Users\Castro\Documents\GAMS files\Aulas\Data Files\PoolingA2.txt

ACTS(S)=yes$(SUPMAX(S) GT 0);
ACTL(L)=yes$(CAPMAX(L) GT 0);
ACTP(P)=yes$(DEMMAX(P) GT 0);
ACTK(K)=yes$(SMAX(P$(ACTP(P)),CMAX(P,K)) GT 0);
ACTSL(S,L)$(ACTS(S) and ACTL(L))=yes$(LINKSL(S,L) EQ 1);
ACTSP(S,P)$(ACTS(S) and ACTP(P))=yes$(LINKSP(S,P) EQ 1);
ACTLP(L,P)$(ACTL(L) and ACTP(P))=yes$(LINKLP(L,P) EQ 1);
FSLMIN(S,L)$(ACTSL(S,L))=0;
FSLMAX(S,L)$(ACTSL(S,L))=MIN(SUPMAX(S),CAPMAX(L));
XMIN(L,P)$(ACTLP(L,P))=0;
XMAX(L,P)$(ACTLP(L,P))=1;
FLOW(ACTSL)=yes;FRAC(ACTLP)=yes;
FIXX(FRAC)=no;FIXF(FLOW)=no;
XRNG(FRAC,'IT0')=XMAX(FRAC)-XMIN(FRAC);
FSLRNG(FLOW,'IT0')=FSLMAX(FLOW)-FSLMIN(FLOW);
CURRENT('ND1')=yes;NEWNODE('ND1')=yes;WAITING('ND1')=yes;
LBX(FRAC,CURRENT)=XMIN(FRAC);UBX(FRAC,CURRENT)=XMAX(FRAC);
LBF(FLOW,CURRENT)=FSLMIN(FLOW);UBF(FLOW,CURRENT)=FSLMAX(FLOW);
PARTVAR=card(FRAC);DISGVAR=card(FLOW);BILINTR=SUM((S,L,P)$(ACTSL(S,L) and ACTLP(L,P)),1);

display ACTSL,ACTSP,ACTLP,ACTK,FSLMIN,FSLMAX,XMIN,XMAX;
Display PARTVAR,DISGVAR,BILINTR;

VARIABLES
PROFIT,VBOUND;

BINARY VARIABLES
Y(L,P,J,KD)         Identifies that fraction from pool L to product P has digit J in position KD

POSITIVE VARIABLES
FSL(S,L)  Flow from source S to pool L
FSP(S,P)  Flow from source S going directly to product P
FSLP(S,L,P)         Flow from source S through pool L to product P
X(L,P)    Split fraction from pool L to product P
XR(L,P)   Residual variable related to split fraction from pool L to product P
FSLXR(S,L,P)        Residual variable related to flow from source S through pool L to product P
FSLDS(S,L,P,J,KD)   Disaggregated flowrate from source S through pool L to product P associated to digit J in position KD
LX(L,P)   Discretized variable related to split fraction from pool L to product P
LXR(L,P)  Residual variable related to split fraction from pool L to product P
LFSLP(S,L,P)        Discretized variable related to flow from source S through pool L to product P
LFSLXR(S,L,P)       Residual variable related to flow from source S through pool L to product P

EQUATIONS
OBJ       Objective function (profit maximization)
EQ1(S)    Balance on availability of source S
EQ2(L)    Upper bound on size of pool L
EQ3(L)    Sum of split fractions out of pool L must be equal to 1
EQ4(P)    Supply of product P must exceed its minimum demand
EQ5(P)    Supply of product P must not exceed its demand
EQ6(P,K)  The level for product P of quality K cannot exceed maximum value
EQ7(P,K)  The level for product P of quality K must exceed minimum value

OBJLP     Linear objective function (profit maximization)
LP4(P)    Linear constraint for ensuring supply of product P must exceed its minimum demand
LP5(P)    Linear constraint for ensuring supply of product P must not exceed its demand
LP6(P,K)  Linear constraint for ensuring the level for product P of quality K cannot exceed maximum value
LP7(P,K)  Linear constraint for ensuring the level for product P of quality K must exceed minimum value
LP8(L,P)  Redundant constraint to improve quality of relaxation
LP9(S,L)  Redundant constraint to improve quality of relaxation
GLMC_1(S,L,P)       Mc Cormick underestimators 1 for flow from source S through pool L to product P
GLMC_2(S,L,P)       Mc Cormick underestimators 2 for flow from source S through pool L to product P
GLMC_3(S,L,P)       Mc Cormick overestimators 1 for flow from source S through pool L to product P
GLMC_4(S,L,P)       Mc Cormick overestimators 2 for flow from source S through pool L to product P

NMDT1(L,P)          Relation between split fraction from pool L to product P and its discretized variable
NMDT2(S,L,P)        Definition of bilinear term involving flow from S to L and split fraction into P
NMDT3(S,L,P)        Relation between bilinear term involving flow from source S through pool L to product P and its residual variable
NMDT4(L,P)          Relation between discretized variable of fraction from pool L to product P and its residual variable
NMDT5(S,L,P,KD)     Relation between disaggregated flow from source S through pool L to product P and position KD with original FLP variable
NMDT6(L,P,KD,J)     Disaggregated flow variable from pool L to product P can only take positive values for position KD value J if the binary has been selected
NMDT7(L,P,KD,J)     Disaggregated flow variable from pool L to product P can only take positive values for position KD value J if the binary has been selected
NMDT8(L,P,KD)       To compute fraction from pool L to product P exactly one value must be selected for position KD
NMDT9(S,L,P)        McCormick relaxation of bilinear term involving flow from source S through pool L to product P and its residual
NMDT10(S,L,P)       McCormick relaxation of bilinear term involving flow from source S through pool L to product P and its residual
NMDT11(S,L,P)       McCormick relaxation of bilinear term involving flow from source S through pool L to product P and its residual
NMDT12(S,L,P)       McCormick relaxation of bilinear term involving flow from source S through pool L to product P and its residual

OBJOBBT   Objective function for optimality-based bound tightening
OBBTEQ1   Profit must be greater than current lower bound
;

OBJ..     PROFIT=e=SUM((S,L,P)$(ACTSL(S,L) and ACTLP(L,P)),(VALUE(P)-PRICE(S))*X(L,P)*FSL(S,L))+SUM((S,P)$(ACTSP(S,P)),(VALUE(P)-PRICE(S))*FSP(S,P));
EQ1(S)$(ACTS(S))..  SUM(L$(ACTSL(S,L)),FSL(S,L))+SUM(P$(ACTSP(S,P)),FSP(S,P))=l=SUPMAX(S);
EQ2(L)$(ACTL(L))..  SUM(S$(ACTSL(S,L)),FSL(S,L))=l=CAPMAX(L);
EQ3(L)$(ACTL(L))..  SUM(P$(ACTLP(L,P)),X(L,P))=e=1;
EQ4(P)$(ACTP(P))..  SUM((S,L)$(ACTSL(S,L) and ACTLP(L,P)),X(L,P)*FSL(S,L))+SUM(S$(ACTSP(S,P)),FSP(S,P))=g=DEMMIN(P);
EQ5(P)$(ACTP(P))..  SUM((S,L)$(ACTSL(S,L) and ACTLP(L,P)),X(L,P)*FSL(S,L))+SUM(S$(ACTSP(S,P)),FSP(S,P))=l=DEMMAX(P);
EQ6(P,K)$(ACTP(P) and ACTK(K))..        SUM((S,L)$(ACTSL(S,L) and ACTLP(L,P)),CS(S,K)*X(L,P)*FSL(S,L))+SUM(S$(ACTSP(S,P)),CS(S,K)*FSP(S,P))=l=(SUM((S,L)$(ACTSL(S,L) and ACTLP(L,P)),X(L,P)*FSL(S,L))+SUM(S$(ACTSP(S,P)),FSP(S,P)))*CMAX(P,K);
EQ7(P,K)$(ACTP(P) and ACTK(K))..        SUM((S,L)$(ACTSL(S,L) and ACTLP(L,P)),CS(S,K)*X(L,P)*FSL(S,L))+SUM(S$(ACTSP(S,P)),CS(S,K)*FSP(S,P))=g=(SUM((S,L)$(ACTSL(S,L) and ACTLP(L,P)),X(L,P)*FSL(S,L))+SUM(S$(ACTSP(S,P)),FSP(S,P)))*CMIN(P,K);

OBJLP..   PROFIT=e=SUM((S,L,P)$(ACTSL(S,L) and ACTLP(L,P)),(VALUE(P)-PRICE(S))*FSLP(S,L,P))+SUM((S,P)$(ACTSP(S,P)),(VALUE(P)-PRICE(S))*FSP(S,P));
LP4(P)$(ACTP(P))..  SUM((S,L)$(ACTSL(S,L) and ACTLP(L,P)),FSLP(S,L,P))+SUM(S$(ACTSP(S,P)),FSP(S,P))=g=DEMMIN(P);
LP5(P)$(ACTP(P))..  SUM((S,L)$(ACTSL(S,L) and ACTLP(L,P)),FSLP(S,L,P))+SUM(S$(ACTSP(S,P)),FSP(S,P))=l=DEMMAX(P);
LP6(P,K)$(ACTP(P) and ACTK(K))..        SUM((S,L)$(ACTSL(S,L) and ACTLP(L,P)),CS(S,K)*FSLP(S,L,P))+SUM(S$(ACTSP(S,P)),CS(S,K)*FSP(S,P))=l=(SUM((S,L)$(ACTSL(S,L) and ACTLP(L,P)),FSLP(S,L,P))+SUM(S$(ACTSP(S,P)),FSP(S,P)))*CMAX(P,K);
LP7(P,K)$(ACTP(P) and ACTK(K))..        SUM((S,L)$(ACTSL(S,L) and ACTLP(L,P)),CS(S,K)*FSLP(S,L,P))+SUM(S$(ACTSP(S,P)),CS(S,K)*FSP(S,P))=g=(SUM((S,L)$(ACTSL(S,L) and ACTLP(L,P)),FSLP(S,L,P))+SUM(S$(ACTSP(S,P)),FSP(S,P)))*CMIN(P,K);
LP8(L,P)$(ACTLP(L,P))..       SUM(S$(ACTSL(S,L)),FSLP(S,L,P))=l=CAPMAX(L)*X(L,P);
LP9(S,L)$(ACTSL(S,L))..       FSL(S,L)=e=SUM(P$(ACTLP(L,P)),FSLP(S,L,P));
GLMC_1(S,L,P)$(ACTSL(S,L) and ACTLP(L,P))..       FSLP(S,L,P)=g=FSL(S,L)*XMIN(L,P)+FSLMIN(S,L)*X(L,P)-FSLMIN(S,L)*XMIN(L,P);
GLMC_2(S,L,P)$(ACTSL(S,L) and ACTLP(L,P))..       FSLP(S,L,P)=g=FSL(S,L)*XMAX(L,P)+FSLMAX(S,L)*X(L,P)-FSLMAX(S,L)*XMAX(L,P);
GLMC_3(S,L,P)$(ACTSL(S,L) and ACTLP(L,P))..       FSLP(S,L,P)=l=FSL(S,L)*XMIN(L,P)+FSLMAX(S,L)*X(L,P)-FSLMAX(S,L)*XMIN(L,P);
GLMC_4(S,L,P)$(ACTSL(S,L) and ACTLP(L,P))..       FSLP(S,L,P)=l=FSL(S,L)*XMAX(L,P)+FSLMIN(S,L)*X(L,P)-FSLMIN(S,L)*XMAX(L,P);

NMDT1(L,P)$(ACTLP(L,P))..     X(L,P)=e=XMIN(L,P)+LX(L,P)*(XMAX(L,P)-XMIN(L,P));
NMDT2(S,L,P)$(ACTSL(S,L) and ACTLP(L,P))..        FSLP(S,L,P)=e=FSL(S,L)*XMIN(L,P)+LFSLP(S,L,P)*(XMAX(L,P)-XMIN(L,P));
NMDT3(S,L,P)$(ACTSL(S,L) and ACTLP(L,P))..        LFSLP(S,L,P)=e=SUM((J,KD)$(ACTJ(J,KD) and ACTKD(KD,L,P)),FSLDS(S,L,P,J,KD)*JVAL(J,KD)*PROD(KDL$(ord(KDL) LE ord(KD)),1/BASIS(KDL)))+LFSLXR(S,L,P);
NMDT4(L,P)$(ACTLP(L,P))..     LX(L,P)=e=SUM((J,KD)$(ACTJ(J,KD) and ACTKD(KD,L,P)),JVAL(J,KD)*Y(L,P,J,KD)*PROD(KDL$(ord(KDL) LE ord(KD)),1/BASIS(KDL)))+LXR(L,P);
NMDT5(S,L,P,KD)$(ACTSL(S,L) and ACTKD(KD,L,P))..  FSL(S,L)=e=SUM(J$(ACTJ(J,KD)),FSLDS(S,L,P,J,KD));
NMDT6(L,P,KD,J)$(ACTJ(J,KD) and ACTKD(KD,L,P))..  SUM(S$(ACTSL(S,L)),FSLDS(S,L,P,J,KD))=l=SUM(S$(ACTSL(S,L)),FSLMAX(S,L))*Y(L,P,J,KD);
NMDT7(L,P,KD,J)$(ACTJ(J,KD) and ACTKD(KD,L,P))..  SUM(S$(ACTSL(S,L)),FSLDS(S,L,P,J,KD))=g=SUM(S$(ACTSL(S,L)),FSLMIN(S,L))*Y(L,P,J,KD);
NMDT8(L,P,KD)$(ACTKD(KD,L,P))..         SUM(J$(ACTJ(J,KD)),Y(L,P,J,KD))=e=1;
NMDT9(S,L,P)$(ACTSL(S,L) and ACTLP(L,P))..        LFSLXR(S,L,P)=l=FSLMAX(S,L)*LXR(L,P);
NMDT10(S,L,P)$(ACTSL(S,L) and ACTLP(L,P))..       LFSLXR(S,L,P)=g=FSLMIN(S,L)*LXR(L,P);
NMDT11(S,L,P)$(ACTSL(S,L) and ACTLP(L,P))..       LFSLXR(S,L,P)=g=(FSL(S,L)-FSLMAX(S,L))*MIN(1,PROD(KD$(ACTKD(KD,L,P)),1/BASIS(KD)))+FSLMAX(S,L)*LXR(L,P);
NMDT12(S,L,P)$(ACTSL(S,L) and ACTLP(L,P))..       LFSLXR(S,L,P)=l=(FSL(S,L)-FSLMIN(S,L))*MIN(1,PROD(KD$(ACTKD(KD,L,P)),1/BASIS(KD)))+FSLMIN(S,L)*LXR(L,P);

OBJOBBT.. VBOUND=e=SUM((S,L)$(BNDF(S,L)),FSL(S,L))+SUM((L,P)$(BNDX(L,P)),X(L,P));
OBBTEQ1.. SUM((S,L,P)$(ACTSL(S,L) and ACTLP(L,P)),(VALUE(P)-PRICE(S))*FSLP(S,L,P))+SUM((S,P)$(ACTSP(S,P)),(VALUE(P)-PRICE(S))*FSP(S,P))=g=LBOUND;


X.lo(L,P)$(ACTLP(L,P))=XMIN(L,P);
X.up(L,P)$(ACTLP(L,P))=XMAX(L,P);
FSL.lo(S,L)$(ACTSL(S,L))=FSLMIN(S,L);
FSL.up(S,L)$(ACTSL(S,L))=FSLMAX(S,L);


OPTION optcr=1E-6;
OPTION limrow=0;
OPTION limcol=0;
OPTION Solprint=Off;
OPTION decimals=6;
OPTION threads=0;
OPTION reslim=3600;
OPTION QCP=GloMIQO;
OPTION NLP=CONOPT;

MODEL Pooltp using /OBJ,EQ1,EQ2,EQ3,EQ4,EQ5,EQ6,EQ7/;#
MODEL PooltpLP using /OBJLP,EQ1,EQ2,EQ3,LP4,LP5,LP6,LP7,LP8,LP9,GLMC_1,GLMC_2,GLMC_3,GLMC_4/; #
MODEL MDT using /OBJLP,EQ1,EQ2,EQ3,LP4,LP5,LP6,LP7,LP8,LP9,NMDT1,NMDT2,NMDT3,NMDT4,NMDT5,NMDT6,NMDT7,NMDT8,NMDT9,NMDT10,NMDT11,NMDT12/;
MODEL OBBT_LP using /OBJOBBT,OBBTEQ1,EQ1,EQ2,EQ3,LP4,LP5,LP6,LP7,LP8,LP9,GLMC_1,GLMC_2,GLMC_3,GLMC_4/;
MODEL OBBT_MDT using /OBJOBBT,OBBTEQ1,EQ1,EQ2,EQ3,LP4,LP5,LP6,LP7,LP8,LP9,NMDT1,NMDT2,NMDT3,NMDT4,NMDT5,NMDT6,NMDT7,NMDT8,NMDT9,NMDT10,NMDT11,NMDT12/;

file OptFile /C:\Users\Castro\Documents\GAMS files\GitHub\Pooling_tp\Cplex.op9/; #Directory must match location of Pooling_tp.gms file
file Results   /C:\Users\Castro\Documents\My Data Sources\GAMS Output\OBBTPool.txt/;
file Search   /C:\Users\Castro\Documents\My Data Sources\GAMS Output\GlobalSearchPool.txt/;
file SPoolRes /C:\Users\Castro\Documents\My Data Sources\GAMS Output\SPoolRes.txt/;
file fsoln;
Results.pw=700;Search.pw=700;SPoolRes.pw=700;
Results.pc=6;Search.pc=6;SPoolRes.pc=6; #Separa com tabs os valores das variaveis

if(WHATM EQ 7,
OPTION NLP=BARON;
SOLVE Pooltp using NLP maximizing PROFIT;
LBOUND=Pooltp.objval;
UBOUND=Pooltp.objest;
OPTGAP=ABS((UBOUND-LBOUND)/UBOUND)*100;
Display FSL.l,FSP.l,X.l,UBOUND,LBOUND,OPTGAP;
);
if(WHATM EQ 8,
OPTION QCP=GUROBI;
Pooltp.optfile=2; #option file has single line: nonconvex 2
SOLVE Pooltp using QCP maximizing PROFIT;
LBOUND=Pooltp.objval;
UBOUND=Pooltp.objest;
OPTGAP=ABS((UBOUND-LBOUND)/UBOUND)*100;
Display FSL.l,FSP.l,X.l,UBOUND,LBOUND,OPTGAP;
);

BNDF(ACTSL)=no;BNDX(ACTLP)=no;

if(WHATM EQ 9,
MDT.optcr=1E-7;
OBBT_MDT.reslim=60;
OBBT_MDT.optcr=1E-5;
put Search;
put 'CPU','LB','UB','Gap (%)','Method','Intervals'/;

#Quick computation of lower and upper bound
SOLVE PooltpLP using LP maximizing PROFIT;
MIPCPU=MAX(PooltpLP.resusd,0.001);
UBOUND=PooltpLP.objval;
SOLVE Pooltp using NLP maximizing PROFIT;
          if(Pooltp.modelstat EQ 2,
          LBOUND=Pooltp.objval;
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
                    loop(FRAC$(not FIXX(FRAC)),
                    BNDX(FRAC)=yes;
                    SOLVE OBBT_LP using LP minimizing VBOUND;
                    LBX(FRAC,CURRENT)=MIN(VBOUND.l,XMAX(FRAC)); #Avoids numerical errors
                    XMIN(FRAC)=SUM(CURRENT,LBX(FRAC,CURRENT));X.lo(FRAC)=XMIN(FRAC);
                    SOLVE OBBT_LP using LP maximizing VBOUND;
                    UBX(FRAC,CURRENT)=MAX(VBOUND.l,XMIN(FRAC)); #Avoids numerical errors
                    XMAX(FRAC)=SUM(CURRENT,UBX(FRAC,CURRENT));X.up(FRAC)=XMAX(FRAC);
                    BNDX(FRAC)=no;
                    );
                    loop(FLOW$(not FIXF(FLOW)),
                    BNDF(FLOW)=yes;
                    SOLVE OBBT_LP using LP minimizing VBOUND;
                    LBF(FLOW,CURRENT)=MIN(VBOUND.l,FSLMAX(FLOW)); #Avoids numerical errors
                    FSLMIN(FLOW)=SUM(CURRENT,LBF(FLOW,CURRENT));FSL.lo(FLOW)=FSLMIN(FLOW);
                    SOLVE OBBT_LP using LP maximizing VBOUND;
                    UBF(FLOW,CURRENT)=MAX(VBOUND.l,FSLMIN(FLOW)); #Avoids numerical errors
                    FSLMAX(FLOW)=SUM(CURRENT,UBF(FLOW,CURRENT));FSL.up(FLOW)=FSLMAX(FLOW);
                    BNDF(FLOW)=no;
                    );
          #Report results of OBBT
          OBBTCPU=TimeElapsed-TSTART;CPUOBBT=MAX(0,CPUOBBT-OBBTCPU);
          XRNG(FRAC,IT)=XMAX(FRAC)-XMIN(FRAC);FSLRNG(FLOW,IT)=FSLMAX(FLOW)-FSLMIN(FLOW);
          XRR(FRAC,IT)=(XRNG(FRAC,'IT0')-XRNG(FRAC,IT))/XRNG(FRAC,'IT0')*100;
          FSLRR(FLOW,IT)=(FSLRNG(FLOW,'IT0')-FSLRNG(FLOW,IT))/FSLRNG(FLOW,'IT0')*100;
          AVERR(IT)=(SUM(FRAC,(XRR(FRAC,IT))$(not FIXX(FRAC))+100$(FIXX(FRAC)))+SUM(FLOW,(FSLRR(FLOW,IT))$(not FIXF(FLOW))+100$(FIXF(FLOW))))/(PARTVAR+DISGVAR);
          FIXX(FRAC)=yes$(XRR(FRAC,IT) GT 99.9999);
          FIXF(FLOW)=yes$(FSLRR(FLOW,IT) GT 99.9999);
          NFIXV(IT)=SUM(FRAC$(FIXX(FRAC)),1)+SUM(FLOW$(FIXF(FLOW)),1);
          NBOUNDP=2*(PARTVAR+DISGVAR-card(FIXX)-card(FIXF));
          #Recompute Bounds
          SOLVE PooltpLP using LP maximizing PROFIT;
          MIPCPU=PooltpLP.resusd;
          UBOUND=PooltpLP.objval;
          SOLVE Pooltp using NLP maximizing PROFIT;
                    if(Pooltp.modelstat EQ 2,
                              if(Pooltp.objval GT LBOUND*1.0001,
                              LBOUND=Pooltp.objval;
                              CURIT=CURIT+1;
                              );
                    );
          OPTGAP=ABS((UBOUND-LBOUND)/UBOUND)*100;
          put TimeElapsed,LBOUND:9:3,UBOUND:9:3,OPTGAP:9:5,'Global McCormick after LP-based OBBT' /;
          putclose Search;
          );
CLOOP$(OPTGAP LE TRGGAP or TimeElapsed GT TCPU)=0;
BASIS(KD)=0;KVAL(KD,L,P)=0;ACTJ(J,KD)=no;JVAL(J,KD)=0;
#Relaxation with mixed-radix MDT
          loop(KDLL$(CLOOP),
          BASIS(KDLL)=2;
          NDIGITS=SUM(KD$(BASIS(KD) GT 0),1);
          CARDK(L,P)$(ACTLP(L,P))=NDIGITS;
          ACTKD(KD,L,P)$(ACTLP(L,P))=yes$(ord(KD) LE CARDK(L,P));
          KVAL(KD,L,P)$(ACTKD(KD,L,P))=-BASIS(KD);
          ACTJ(J,KD)$(ord(KD) LE NDIGITS)=yes$(ord(J) LE BASIS(KD));
          JVAL(J,KD)$(ACTJ(J,KD))=ord(J)-1;
          LXR.up(L,P)$(ACTLP(L,P))=MIN(1,PROD(KD$(ACTKD(KD,L,P)),1/BASIS(KD)));
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
                              SOLVE Pooltp using NLP maximizing PROFIT;
                                        if(Pooltp.modelstat EQ 2,
                                        LBSP(SP)=PROFIT.l; put LBSP(SP):9:5 /;
                                        LBOUND=MAX(LBOUND,Pooltp.objval);
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
                              loop(FRAC$(not FIXX(FRAC)),
                              BNDX(FRAC)=yes;
                              SOLVE OBBT_MDT using MIP minimizing VBOUND;
                                        if(OBBT_MDT.modelstat NE 10 and OBBT_MDT.modelstat NE 14 and OBBT_MDT.solvestat NE 4,
                                        LBX(FRAC,CURRENT)=MIN(OBBT_MDT.objest,XMAX(FRAC)); #Avoids numerical errors
                                        else
                                        SOLVE OBBT_LP using LP minimizing VBOUND;
                                                  if(OBBT_LP.modelstat NE 4,
                                                  LBX(FRAC,CURRENT)=MIN(VBOUND.l,XMAX(FRAC));
                                                  );
                                        );
                              XMIN(FRAC)=SUM(CURRENT,LBX(FRAC,CURRENT));X.lo(FRAC)=XMIN(FRAC);
                              SOLVE OBBT_MDT using MIP maximizing VBOUND;
                                        if(OBBT_MDT.modelstat NE 10 and OBBT_MDT.modelstat NE 14 and OBBT_MDT.solvestat NE 4,
                                        UBX(FRAC,CURRENT)=MAX(OBBT_MDT.objest,XMIN(FRAC)); #Avoids numerical errors
                                        else
                                        SOLVE OBBT_LP using LP maximizing VBOUND;
                                                  if(OBBT_LP.modelstat NE 4,
                                                  UBX(FRAC,CURRENT)=MAX(VBOUND.l,XMIN(FRAC));
                                                  );
                                        );
                              XMAX(FRAC)=SUM(CURRENT,UBX(FRAC,CURRENT));X.up(FRAC)=XMAX(FRAC);
                              BNDX(FRAC)=no;
                              );
                              loop(FLOW$(not FIXF(FLOW)),
                              BNDF(FLOW)=yes;
                              SOLVE OBBT_MDT using MIP minimizing VBOUND;
                                        if(OBBT_MDT.modelstat NE 10 and OBBT_MDT.modelstat NE 14 and OBBT_MDT.solvestat NE 4,
                                        LBF(FLOW,CURRENT)=MIN(OBBT_MDT.objest,FSLMAX(FLOW)); #Avoids numerical errors
                                        else
                                        SOLVE OBBT_LP using LP minimizing VBOUND;
                                                  if(OBBT_LP.modelstat NE 4,
                                                  LBF(FLOW,CURRENT)=MIN(VBOUND.l,FSLMAX(FLOW));
                                                  );
                                        );
                              FSLMIN(FLOW)=SUM(CURRENT,LBF(FLOW,CURRENT));FSL.lo(FLOW)=FSLMIN(FLOW);
                              SOLVE OBBT_MDT using MIP maximizing VBOUND;
                                        if(OBBT_MDT.modelstat NE 10 and OBBT_MDT.modelstat NE 14 and OBBT_MDT.solvestat NE 4,
                                        UBF(FLOW,CURRENT)=MAX(OBBT_MDT.objest,FSLMIN(FLOW)); #Avoids numerical errors
                                        else
                                        SOLVE OBBT_LP using LP maximizing VBOUND;
                                                  if(OBBT_LP.modelstat NE 4,
                                                  UBF(FLOW,CURRENT)=MAX(VBOUND.l,FSLMIN(FLOW));
                                                  );
                                        );
                              FSLMAX(FLOW)=SUM(CURRENT,UBF(FLOW,CURRENT));FSL.up(FLOW)=FSLMAX(FLOW);
                              BNDF(FLOW)=no;
                              );
                    #Report results of OBBT
                    OBBTCPU=TimeElapsed-TSTART;CPUOBBT=MAX(0,CPUOBBT-OBBTCPU);
                    CURIT=CURIT+1;
                    XRNG(FRAC,IT)$(ord(IT) EQ CURIT+1)=XMAX(FRAC)-XMIN(FRAC);FSLRNG(FLOW,IT)$(ord(IT) EQ CURIT+1)=FSLMAX(FLOW)-FSLMIN(FLOW);
                    XRR(FRAC,IT)$(ord(IT) EQ CURIT+1)=(XRNG(FRAC,'IT0')-XRNG(FRAC,IT))/XRNG(FRAC,'IT0')*100;
                    FSLRR(FLOW,IT)$(ord(IT) EQ CURIT+1)=(FSLRNG(FLOW,'IT0')-FSLRNG(FLOW,IT))/FSLRNG(FLOW,'IT0')*100;
                    AVERR(IT)$(ord(IT) EQ CURIT+1)=(SUM(FRAC,(XRR(FRAC,IT))$(not FIXX(FRAC))+100$(FIXX(FRAC)))+SUM(FLOW,(FSLRR(FLOW,IT))$(not FIXF(FLOW))+100$(FIXF(FLOW))))/(PARTVAR+DISGVAR);
                    FIXX(FRAC)=yes$(SUM(IT$(ord(IT) EQ CURIT+1),XRR(FRAC,IT)) GT 99.9999);
                    FIXF(FLOW)=yes$(SUM(IT$(ord(IT) EQ CURIT+1),FSLRR(FLOW,IT)) GT 99.9999);
                    NFIXV(IT)$(ord(IT) EQ CURIT+1)=SUM(FRAC$(FIXX(FRAC)),1)+SUM(FLOW$(FIXF(FLOW)),1);
                    NBOUNDP=2*(PARTVAR+DISGVAR-card(FIXX)-card(FIXF));
                    #Check if with results from McCormick relaxation we can already stop
                    SOLVE PooltpLP using LP maximizing PROFIT;
                    UBOUND=MIN(UBOUND,PooltpLP.objval);
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
          ERRF(S,L,P)$(ACTSL(S,L) and ACTLP(L,P))=ABS(LFSLP.l(S,L,P)-FSL.l(S,L)*LX.l(L,P));
          BOUND('ND1')=UBOUND;
          SRCNODE('ND1')=yes;
          LOGBB('ND1','Time (s)')=TimeElapsed;LOGBB('ND1','Explored')=card(SRCNODE);LOGBB('ND1','Remaining')=card(WAITING);
          LOGBB('ND1','UBOUND')=UBOUND;LOGBB('ND1','LBOUND')=LBOUND;LOGBB('ND1','OPTGAP (%)')=OPTGAP;
          LOGND(CURRENT,FLOW,'XLO')=LBF(FLOW,CURRENT);LOGND(CURRENT,FLOW,'XUP')=UBF(FLOW,CURRENT);
          LOGND2(CURRENT,FRAC,'XLO')=LBX(FRAC,CURRENT);LOGND2(CURRENT,FRAC,'XUP')=UBX(FRAC,CURRENT);
          put LOGBB('ND1','Time (s)'):9:2, put LOGBB('ND1','LBOUND'):9:3, put LOGBB('ND1','UBOUND'):9:3, put LOGBB('ND1','OPTGAP (%)'):9:5, put LOGBB('ND1','Explored'):9:0 ; put LOGBB('ND1','Remaining'):9:0 /;
          putclose Search;
          #Choose branching variavel
                    loop((S,L)$(ACTSL(S,L) and SMAX(P$(ACTLP(L,P)),ERRF(S,L,P)) EQ SMAX((SL,LL,P)$(ACTSL(SL,LL) and ACTLP(LL,P)),ERRF(SL,LL,P))),
                    NEXTF(CURRENT,S,L)=yes;
                    );
          loop(CURRENT,BRCHF(FLOW)=NEXTF(CURRENT,FLOW););
          LOGND(CURRENT,BRCHF,'LBOUND')=LBOUND;LOGND(CURRENT,BRCHF,'UBOUND')=UBOUND;LOGND(CURRENT,BRCHF,'Next Src.')=SUM((S,L)$(NEXTF(CURRENT,S,L)),ord(S));LOGND(CURRENT,BRCHF,'Next Pool')=SUM((S,L)$(NEXTF(CURRENT,S,L)),ord(L));
          #Perform branching
                    loop(ND$(not DONE), # and ord(ND) LE 1
                    #Node selection
                    CURRENT(NDL)=no;
                    CURRENT(WAITING(NDL))$(BOUND(NDL) EQ UBOUND)=yes;
                    FIRST=1;
                              loop(CURRENT$FIRST, #Select only one node for branching
                              WAITING(CURRENT)=no;
                              loop(NDL$(CURRENT(NDL)),BRCHF(S,L)=NEXTF(NDL,S,L);); #Assign next variavel to branch
                              LOGIT(ND,'Node_Sel.')=SUM(NDL$(CURRENT(NDL)),ord(NDL));
                              LOGIT(ND,'From_Src.')=SUM((S,L)$(BRCHF(S,L)),ord(S));
                              LOGIT(ND,'To_Pool')=SUM((S,L)$(BRCHF(S,L)),ord(L));
                              #Branch
                                        for(BISECT=1 to 2,
                                        NEWNODE(NDL)=NEWNODE(NDL-1);
                                        NODEBRN(CURRENT,NEWNODE)=yes;
                                        SRCNODE(NEWNODE)=yes;
                                        WAITING(NEWNODE)=yes;
                                        LBF(FLOW,NEWNODE)=LBF(FLOW,CURRENT);UBF(FLOW,NEWNODE)=UBF(FLOW,CURRENT);
                                        LBX(FRAC,NEWNODE)=LBX(FRAC,CURRENT);UBX(FRAC,NEWNODE)=UBX(FRAC,CURRENT);
                                                  if(BISECT EQ 1, #Left branch
                                                  LBF(BRCHF,NEWNODE)=LBF(BRCHF,CURRENT);
                                                  UBF(BRCHF,NEWNODE)=LBF(BRCHF,CURRENT)+(UBF(BRCHF,CURRENT)-LBF(BRCHF,CURRENT))/2;
                                                  else
                                                  LBF(BRCHF,NEWNODE)=LBF(BRCHF,CURRENT)+(UBF(BRCHF,CURRENT)-LBF(BRCHF,CURRENT))/2;
                                                  UBF(BRCHF,NEWNODE)=UBF(BRCHF,CURRENT);
                                                  );
                                        FSLMIN(FLOW)=SUM(NDL$(NEWNODE(NDL)),LBF(FLOW,NDL));FSLMAX(FLOW)=SUM(NDL$(NEWNODE(NDL)),UBF(FLOW,NDL));
                                        XMIN(FRAC)=SUM(NDL$(NEWNODE(NDL)),LBX(FRAC,NDL));XMAX(FRAC)=SUM(NDL$(NEWNODE(NDL)),UBX(FRAC,NDL));
                                        FSL.lo(FLOW)=FSLMIN(FLOW);FSL.up(FLOW)=FSLMAX(FLOW);
                                        X.lo(FRAC)=XMIN(FRAC);X.up(FRAC)=XMAX(FRAC);
                                        LOGND(NEWNODE,FLOW,'XLO')=FSL.lo(FLOW);LOGND(NEWNODE,FLOW,'XUP')=FSL.up(FLOW);
                                        LOGND2(NEWNODE,FRAC,'XLO')=X.lo(FRAC);LOGND2(NEWNODE,FRAC,'XUP')=X.up(FRAC);
                                        #Bound contract
                                        STATUS=0;
                                                  loop(FRAC$(not STATUS and not FIXX(FRAC)),
                                                  BNDX(FRAC)=yes;
                                                            if(WHATSBB EQ 1,
                                                            #LP-based OBBT
                                                            SOLVE OBBT_LP using LP minimizing VBOUND;
                                                                      if(OBBT_LP.modelstat EQ 4, #Branch can be fathomed
                                                                      STATUS=1;
                                                                      WAITING(NEWNODE)=no;
                                                                      else
                                                                      LBX(FRAC,NEWNODE)=MIN(VBOUND.l,XMAX(FRAC)); #Avoids numerical errors
                                                                      XMIN(FRAC)=SUM(NEWNODE,LBX(FRAC,NEWNODE));X.lo(FRAC)=XMIN(FRAC);
                                                                      SOLVE OBBT_LP using LP maximizing VBOUND;
                                                                      UBX(FRAC,NEWNODE)=MAX(VBOUND.l,XMIN(FRAC)); #Avoids numerical errors
                                                                      XMAX(FRAC)=SUM(NEWNODE,UBX(FRAC,NEWNODE));X.up(FRAC)=XMAX(FRAC);
                                                                      LOGND2(NEWNODE,FRAC,'XLO')=X.lo(FRAC);LOGND2(NEWNODE,FRAC,'XUP')=X.up(FRAC);
                                                                      );
                                                            else
                                                            #MIP-based OBBT
                                                            KVAL(KD,L,P)=0;ACTJ(J,KD)=no;JVAL(J,KD)=0;
                                                            BASIS(KD)=2$(ord(KD) LE SBBOBBT);
                                                            NDIGITS=SUM(KD$(BASIS(KD) GT 0),1);
                                                            CARDK(L,P)$(ACTLP(L,P))=NDIGITS;
                                                            ACTKD(KD,L,P)$(ACTLP(L,P))=yes$(ord(KD) LE CARDK(L,P));
                                                            KVAL(KD,L,P)$(ACTKD(KD,L,P))=-BASIS(KD);
                                                            ACTJ(J,KD)$(ord(KD) LE NDIGITS)=yes$(ord(J) LE BASIS(KD));
                                                            JVAL(J,KD)$(ACTJ(J,KD))=ord(J)-1;
                                                            LXR.up(L,P)$(ACTLP(L,P))=MIN(1,PROD(KD$(ACTKD(KD,L,P)),1/BASIS(KD)));
                                                            SOLVE OBBT_MDT using MIP minimizing VBOUND;
                                                                      if(OBBT_MDT.modelstat NE 10 and OBBT_MDT.modelstat NE 14 and OBBT_MDT.solvestat NE 4,
                                                                      LBX(FRAC,NEWNODE)=MIN(OBBT_MDT.objest,XMAX(FRAC)); #Avoids numerical errors
                                                                      XMIN(FRAC)=SUM(NEWNODE,LBX(FRAC,NEWNODE));X.lo(FRAC)=XMIN(FRAC);
                                                                      else
                                                                      SOLVE OBBT_LP using LP minimizing VBOUND;
                                                                                if(OBBT_LP.modelstat EQ 4, #Branch can be fathomed
                                                                                STATUS=1;
                                                                                WAITING(NEWNODE)=no;
                                                                                else
                                                                                LBX(FRAC,NEWNODE)=MIN(VBOUND.l,XMAX(FRAC)); #Avoids numerical errors
                                                                                XMIN(FRAC)=SUM(NEWNODE,LBX(FRAC,NEWNODE));X.lo(FRAC)=XMIN(FRAC);
                                                                                );
                                                                      );
                                                                      if(not STATUS,
                                                                      SOLVE OBBT_MDT using MIP maximizing VBOUND;
                                                                                if(OBBT_MDT.modelstat NE 10 and OBBT_MDT.modelstat NE 14 and OBBT_MDT.solvestat NE 4,
                                                                                UBX(FRAC,NEWNODE)=MAX(OBBT_MDT.objest,XMIN(FRAC)); #Avoids numerical errors
                                                                                XMAX(FRAC)=SUM(NEWNODE,UBX(FRAC,NEWNODE));X.up(FRAC)=XMAX(FRAC);
                                                                                LOGND2(NEWNODE,FRAC,'XLO')=X.lo(FRAC);LOGND2(NEWNODE,FRAC,'XUP')=X.up(FRAC);
                                                                                else
                                                                                SOLVE OBBT_LP using LP maximizing VBOUND;
                                                                                          if(OBBT_LP.modelstat EQ 4, #Branch can be fathomed
                                                                                          STATUS=1;
                                                                                          WAITING(NEWNODE)=no;
                                                                                          else
                                                                                          UBX(FRAC,NEWNODE)=MAX(VBOUND.l,XMIN(FRAC)); #Avoids numerical errors
                                                                                          XMAX(FRAC)=SUM(NEWNODE,UBX(FRAC,NEWNODE));X.up(FRAC)=XMAX(FRAC);
                                                                                          LOGND2(NEWNODE,FRAC,'XLO')=X.lo(FRAC);LOGND2(NEWNODE,FRAC,'XUP')=X.up(FRAC);
                                                                                          );
                                                                                );
                                                                      );
                                                            );
                                                  BNDX(FRAC)=no;
                                                  );
                                                  loop((SL,LL)$(FLOW(SL,LL) and not BRCHF(SL,LL) and not STATUS and not FIXF(SL,LL)),
                                                  BNDF(SL,LL)=yes;
                                                            if(WHATSBB EQ 1,
                                                            #LP-based OBBT
                                                            SOLVE OBBT_LP using LP minimizing VBOUND;
                                                                      if(OBBT_LP.modelstat EQ 4, #Branch can be fathomed
                                                                      STATUS=1;
                                                                      WAITING(NEWNODE)=no;
                                                                      else
                                                                      LBF(SL,LL,NEWNODE)=MIN(VBOUND.l,FSLMAX(SL,LL)); #Avoids numerical errors
                                                                      FSLMIN(SL,LL)=SUM(NEWNODE,LBF(SL,LL,NEWNODE));FSL.lo(SL,LL)=FSLMIN(SL,LL);
                                                                      SOLVE OBBT_LP using LP maximizing VBOUND;
                                                                      UBF(SL,LL,NEWNODE)=MAX(VBOUND.l,FSLMIN(SL,LL)); #Avoids numerical errors
                                                                      FSLMAX(SL,LL)=SUM(NEWNODE,UBF(SL,LL,NEWNODE));FSL.up(SL,LL)=FSLMAX(SL,LL);
                                                                      LOGND(NEWNODE,SL,LL,'XLO')=FSL.lo(SL,LL);LOGND(NEWNODE,SL,LL,'XUP')=FSL.up(SL,LL);
                                                                      );
                                                            else
                                                            #MIP-based OBBT
                                                            SOLVE OBBT_MDT using MIP minimizing VBOUND;
                                                                      if(OBBT_MDT.modelstat NE 10 and OBBT_MDT.modelstat NE 14 and OBBT_MDT.solvestat NE 4,
                                                                      LBF(SL,LL,NEWNODE)=MIN(OBBT_MDT.objest,FSLMAX(SL,LL));
                                                                      FSLMIN(SL,LL)=SUM(NEWNODE,LBF(SL,LL,NEWNODE));FSL.lo(SL,LL)=FSLMIN(SL,LL);
                                                                      else
                                                                      SOLVE OBBT_LP using LP minimizing VBOUND;
                                                                                if(OBBT_LP.modelstat EQ 4, #Branch can be fathomed
                                                                                STATUS=1;
                                                                                WAITING(NEWNODE)=no;
                                                                                else
                                                                                LBF(SL,LL,NEWNODE)=MIN(VBOUND.l,FSLMAX(SL,LL)); #Avoids numerical errors
                                                                                FSLMIN(SL,LL)=SUM(NEWNODE,LBF(SL,LL,NEWNODE));FSL.lo(SL,LL)=FSLMIN(SL,LL);
                                                                                );
                                                                      );
                                                                      if(not STATUS,
                                                                      SOLVE OBBT_MDT using MIP maximizing VBOUND;
                                                                                if(OBBT_MDT.modelstat NE 10 and OBBT_MDT.modelstat NE 14 and OBBT_MDT.solvestat NE 4,
                                                                                UBF(SL,LL,NEWNODE)=MAX(OBBT_MDT.objest,FSLMIN(SL,LL));
                                                                                FSLMAX(SL,LL)=SUM(NEWNODE,UBF(SL,LL,NEWNODE));FSL.up(SL,LL)=FSLMAX(SL,LL);
                                                                                LOGND(NEWNODE,SL,LL,'XLO')=FSL.lo(SL,LL);LOGND(NEWNODE,SL,LL,'XUP')=FSL.up(SL,LL);
                                                                                else
                                                                                SOLVE OBBT_LP using LP maximizing VBOUND;
                                                                                          if(OBBT_LP.modelstat EQ 4, #Branch can be fathomed
                                                                                          STATUS=1;
                                                                                          WAITING(NEWNODE)=no;
                                                                                          else
                                                                                          UBF(SL,LL,NEWNODE)=MAX(VBOUND.l,FSLMIN(SL,LL)); #Avoids numerical errors
                                                                                          FSLMAX(SL,LL)=SUM(NEWNODE,UBF(SL,LL,NEWNODE));FSL.up(SL,LL)=FSLMAX(SL,LL);
                                                                                          LOGND(NEWNODE,SL,LL,'XLO')=FSL.lo(SL,LL);LOGND(NEWNODE,SL,LL,'XUP')=FSL.up(SL,LL);
                                                                                          );
                                                                                );
                                                                      );
                                                            );
                                                  BNDF(SL,LL)=no;
                                                  );
                                        #MDT Relaxation
                                                  if(not STATUS,
                                                  KVAL(KD,L,P)=0;ACTJ(J,KD)=no;JVAL(J,KD)=0;
                                                  BASIS(KD)=2$(ord(KD) LE SBBRELX);
                                                  NDIGITS=SUM(KD$(BASIS(KD) GT 0),1);
                                                  CARDK(L,P)$(ACTLP(L,P))=NDIGITS;
                                                  ACTKD(KD,L,P)$(ACTLP(L,P))=yes$(ord(KD) LE CARDK(L,P));
                                                  KVAL(KD,L,P)$(ACTKD(KD,L,P))=-BASIS(KD);
                                                  ACTJ(J,KD)$(ord(KD) LE NDIGITS)=yes$(ord(J) LE BASIS(KD));
                                                  JVAL(J,KD)$(ACTJ(J,KD))=ord(J)-1;
                                                  LXR.up(L,P)$(ACTLP(L,P))=MIN(1,PROD(KD$(ACTKD(KD,L,P)),1/BASIS(KD)));
                                                  #No need for solution pool at this stage
                                                  put OptFile;
                                                  put 'tilim ', (MAX(0,TCPU-TimeElapsed)):<9:0 /;
                                                  putclose OptFile;
                                                  MDT.cutoff=LBOUND;
                                                  SOLVE MDT using MIP maximizing PROFIT;
                                                            if(MDT.modelstat EQ 10,BOUND(NEWNODE)=-INF;
                                                            else BOUND(NEWNODE)=MIN(MDT.objest,BOUND(CURRENT));
                                                                      if(MDT.modelstat NE 14,
                                                                      ERRF(S,L,P)$(ACTSL(S,L) and ACTLP(L,P))=ABS(LFSLP.l(S,L,P)-FSL.l(S,L)*LX.l(L,P));
                                                                      else  #Use LP relaxation to compute the error that will be used to select branching variable
                                                                      SOLVE PooltpLP using LP maximizing PROFIT;
                                                                      ERRF(S,L,P)$(ACTSL(S,L) and ACTLP(L,P))=ABS(FSLP.l(S,L,P)-FSL.l(S,L)*X.l(L,P));
                                                                      );
                                                            );
                                                  LOGND(NEWNODE,BRCHF,'UBOUND')=BOUND(NEWNODE);
                                                            if(SUM(NEWNODE,BOUND(NEWNODE)) LT LBOUND/(1-TRGGAP*0.01),
                                                            WAITING(NEWNODE)=no; #Fathom node
                                                            else
                                                            #Choose branching variavel
                                                                      loop((S,L)$(ACTSL(S,L) and SMAX(P$(ACTLP(L,P)),ERRF(S,L,P)) EQ SMAX((SL,LL,P)$(ACTSL(SL,LL) and ACTLP(LL,P)),ERRF(SL,LL,P))),
                                                                      NEXTF(NEWNODE,S,L)=yes;
                                                                      );
                                                            LOGND(NEWNODE,BRCHF,'Next Src.')=SUM((S,L)$(NEXTF(NEWNODE,S,L)),ord(S));
                                                            LOGND(NEWNODE,BRCHF,'Next Pool')=SUM((S,L)$(NEXTF(NEWNODE,S,L)),ord(L));
                                                            FSL.lo(FLOW)=LBF(FLOW,'ND1');FSL.up(FLOW)=UBF(FLOW,'ND1');
                                                            SOLVE Pooltp using NLP maximizing PROFIT;
                                                                      if(Pooltp.modelstat EQ 2,
                                                                      LOGND(NEWNODE,BRCHF,'LBOUND')=Pooltp.objval;
                                                                                if(Pooltp.objval GT LBOUND,
                                                                                LBOUND=Pooltp.objval;
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
put 'OBBT'/;
put 'Var.', 'RangeR (%)'/;
          loop(FRAC,
          put FRAC.te(FRAC);
                    loop(IT$(ord(IT) GT 1 and ord(IT) LE CURIT+1),
                    put XRR(FRAC,IT):8:4;
                    );
          put/;
          );
          loop(FLOW,
          put FLOW.te(FLOW);
                    loop(IT$(ord(IT) GT 1 and ord(IT) LE CURIT+1),
                    put FSLRR(FLOW,IT):8:4;
                    );
          put/;
          );
put 'Fixed V',loop(IT$(ord(IT) GT 1 and ord(IT) LE CURIT+1), put NFIXV(IT);); put/;
put 'ARR',loop(IT$(ord(IT) GT 1 and ord(IT) LE CURIT+1), put AVERR(IT);); put/;
display NBOUNDP,UBOUND,LBOUND,OPTGAP;
);
