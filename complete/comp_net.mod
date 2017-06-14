/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Model under complete network structure  %
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

@#define shockconfig = "SC2"
@#define liqinj      = "on"


/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Declaration of variables: Total = 139
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

var d_A d_B d_C d_D c_A c_B c_C c_D w_A w_B w_C w_D n_A n_B n_C n_D            // Household variables (16)

k_A k_B k_C k_D y_A y_B y_C y_D x_A x_B x_C x_D alp_A alp_B alp_C alp_D
fpr_A fpr_B fpr_C fpr_D svk_A svk_B svk_C svk_D                                // Firm variables      (24)

del_AB del_AC del_AD l_AB l_AC l_AD b_AB b_AC b_AD
del_BA del_BC del_BD l_BA l_BC l_BD b_BA b_BC b_BD 
del_CA del_CB del_CD l_CA l_CB l_CD b_CA b_CB b_CD
del_DA del_DB del_DC l_DA l_DB l_DC b_DA b_DB b_DC 

F_A F_B F_C F_D  bpr_A bpr_B bpr_C bpr_D
svp_A svp_B svp_C svp_D svf_A svf_B svf_C svf_D                                // Bank variables      (52)

rd_A rd_B rd_C rd_D rl_A rl_B rl_C rl_D rd_av rl_av rb_av       
rb_AB rb_AC rb_AD
rb_BA rb_BC rb_BD
rb_CA rb_CB rb_CD
rb_DA rb_DB rb_DC

IBspread_A IBspread_D IBspread_B IBspread_C IBspread_AB IBspread_BD ratespread        

T_A T_B T_C T_D tot_bpr tot_fpr tot_y tot_x  tot_l tot_b                             // Other              (9)
 
@#include "liqinjvardefs_comp.mod" 
@#include "shockvardefs.mod"

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Declaration of stochastic shocks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

@#include "liqinjvarexo.mod"

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Declaration of parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

parameters chi beta dpb dpf rwf rwb rws nu etaf etab tau  mu k ykr num_links num_banks
rhob_Gamma rho_Gamma rho_tfp  alp_ss del_ss d_ss n_ss S_ss rb_ss rl_ss rd_ss         // Calibration: Parameters and steady states

omegab db dF xi l_ss b_ss bpr_ss F_ss Gamma_ss svf_ss svp_ss omegaf df
tfp_ss  k_ss x_ss y_ss w_ss svk_ss fpr_ss T_ss c_ss                          // Implied steady state parameters 

@#include "liqinjpardefs.mod"
@#include "shockpardefs.mod"
                                                                                      
/*----------------------------------------------------------------------------------------
Calibrations
----------------------------------------------------------------------------------------*/

chi    = 0.01;
dpb    = 0.80;
dpf    = 0.80;
rwf    = 0.8;
rwb    = 0.05;
rws    = 1.2;
nu     = 0.5;
etaf   = 100;
etab   = 200;
tau    = 0.03;
mu     = 1/3;
k      = 0.08;
ykr    = 0.1;

rhob_Gamma = 0.03;

rho_Gamma  = 0.5; 
rho_tfp    = 0.95;

alp_ss = 0.95;
del_ss = 0.99;

n_ss = 0.20;

//rb_ss    = 0.007;
//rl_ss    = 0.016;
//rd_ss    = 0.0035;

//rb_ss = 0.01;
//rl_ss = 0.008;
//rd_ss = 0.004;

//rb_ss = 0.01;
//rl_ss = 0.005;
//rd_ss = 0.005;

rb_ss = 0.012;
rl_ss = 0.001;
rd_ss = 0.005;

beta     = 1/(1+rd_ss);

k_ss     = n_ss/(ykr)^(1/(1-mu));
x_ss     = tau*k_ss*(1+rl_ss);
S_ss   = 1*x_ss;

tfp_ss   = 1;
Gamma_ss = rhob_Gamma*S_ss; 

@#include "liqinjpareqs.mod"
@#include "shockpareqs.mod"

num_links = 12;
num_banks = 4;

/*----------------------------------------------------------------------------------------
Parameters specific to the firms
----------------------------------------------------------------------------------------*/

y_ss     = k_ss*ykr;
w_ss     = (1-mu)*(k_ss/n_ss)^mu;

svk_ss   = mu*(k_ss/n_ss)^(mu-1)/(1-beta*(1-tau));
omegaf   = (svk_ss/(1+rl_ss)-beta*alp_ss)/(beta^2*(1-alp_ss)^2*x_ss);
df       = x_ss-beta*omegaf*(1-alp_ss)*x_ss^2;

fpr_ss   = tfp_ss*y_ss-w_ss*n_ss-alp_ss*x_ss-omegaf/2*((1-alp_ss)*x_ss)^2; 
        
/*----------------------------------------------------------------------------------------
Parameters specific to the banks
----------------------------------------------------------------------------------------*/

l_ss   = 0.5*x_ss;
b_ss   = l_ss;
d_ss   = 2*x_ss;
F_ss   = k*(rwb*l_ss+rwf*x_ss+rws*S_ss); 

omegab = (1/(1+rb_ss)-beta*del_ss)/(beta^2*(1-del_ss)^2*b_ss);       // Derived from Equation ()

bpr_ss = (1/(1+rd_ss)-1)*d_ss+3*(1/(1+rb_ss)-del_ss)*b_ss
        +(alp_ss+dpf*(1-alp_ss)-1/(1+rl_ss))*x_ss
        + 3*(del_ss+dpb*(1-del_ss)-1/(1+rb_ss))*l_ss
        -3*(omegab/2)*((1-del_ss)*b_ss)^2+Gamma_ss;  

xi     = nu*bpr_ss/F_ss;                                             // Derived from Equation ()

svp_ss = (1/bpr_ss)/(1-nu*(beta*del_ss+beta^2*dpb*(1-del_ss)-1/(1+rb_ss)))/(k*rwb*(1-beta*(1-xi)));
svp_ss = (1/bpr_ss)/(1-nu*(beta*alp_ss+beta^2*dpf*(1-alp_ss)-1/(1+rl_ss)))/(k*rwf*(1-beta*(1-xi)));

svf_ss = (svp_ss-(1/bpr_ss))/nu;                                     // Derived from Equation ()
dF     = svf_ss*(1-beta*(1-xi));                                     // Derived from Equation ()
db     = svp_ss*(b_ss-beta*omegab*(1-del_ss)*b_ss^2);                // Derived from Equation ()

/*----------------------------------------------------------------------------------------
Parameters specific to the government
----------------------------------------------------------------------------------------*/

T_ss = dpf*(1-alp_ss)*x_ss+dpb*(1-del_ss)*l_ss-xi*F_ss;
//T_ss   = 0.00402879/4;

/*----------------------------------------------------------------------------------------
Parameters specific to households
----------------------------------------------------------------------------------------*/

c_ss = w_ss*n_ss+(1-1/(1+rd_ss))*d_ss+l_ss-xi*F_ss-dpf*(1-alp_ss)*x_ss-dpb*(1-del_ss)*l_ss;

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Model equations: Total = 139
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

model;

/*----------------------------------------------------------------------------------------
Household equations: 4 households x 3 equations = 12 equations
----------------------------------------------------------------------------------------*/

@#for region in ["A","B","C","D"]
    T_@{region}*y_@{region}+c_@{region}+d_@{region}/(1+rd_@{region}) 
        = w_@{region}*n_@{region}+d_@{region}(-1);
    (1/c_@{region})*(1/(1+rd_@{region})) 
        = beta/c_@{region}(+1) - chi*(d_@{region}/(1+rd_@{region})-d_ss/(1+rd_ss));
    n_@{region} = n_ss;
@#endfor

/*----------------------------------------------------------------------------------------
Bank equations: 4 banks x 15 equations = 60 equations
----------------------------------------------------------------------------------------*/

// A

svp_A*b_AB(-1)  = beta*svp_A(+1)*omegab*(1-del_AB)*(b_AB(-1))^2 + db;                                                    // Equation (10): FOC wrt del_AB
svp_A*b_AC(-1)  = beta*svp_A(+1)*omegab*(1-del_AC)*(b_AC(-1))^2 + db;                                                    // Equation (10): FOC wrt del_AC
svp_A*b_AD(-1)  = beta*svp_A(+1)*omegab*(1-del_AD)*(b_AD(-1))^2 + db;                                                    // Equation (10): FOC wrt del_AD

svp_A/(1+rb_AB) = beta*svp_A(+1)*del_BA(+1)+beta^2*svp_A(+2)*dpb*(1-del_BA(+1))-dF*k*rwb*(del_ss/del_BA(+1))^etab;       // Equation (11): FOC wrt l_AB
svp_A/(1+rb_AC) = beta*svp_A(+1)*del_CA(+1)+beta^2*svp_A(+2)*dpb*(1-del_CA(+1))-dF*k*rwb*(del_ss/del_CA(+1))^etab;       // Equation (11): FOC wrt l_AC 
svp_A/(1+rb_AD) = beta*svp_A(+1)*del_DA(+1)+beta^2*svp_A(+2)*dpb*(1-del_DA(+1))-dF*k*rwb*(del_ss/del_DA(+1))^etab;       // Equation (11): FOC wrt l_AD 
 
svp_A/(1+rb_BA) = beta*svp_A(+1)*del_AB(+1)+beta^2*svp_A(+2)*omegab*(1-del_AB(+1))^2*(b_AB);                             // Equation (12): FOC wrt b_AB
svp_A/(1+rb_CA) = beta*svp_A(+1)*del_AC(+1)+beta^2*svp_A(+2)*omegab*(1-del_AC(+1))^2*(b_AC);                             // Equation (12): FOC wrt b_AC
svp_A/(1+rb_DA) = beta*svp_A(+1)*del_AD(+1)+beta^2*svp_A(+2)*omegab*(1-del_AD(+1))^2*(b_AD);                             // Equation (12): FOC wrt b_AD

svp_A/(1+rl_A)  = beta*svp_A(+1)*alp_A(+1)+beta^2*svp_A(+2)*dpf*(1-alp_A(+1))-dF*k*rwf*(alp_ss/alp_A(+1))^etaf;          // Equation (13): FOC wrt x_A
svp_A/(1+rd_A)  = beta*svp_A(+1);                                                                                        // Equation (14): FOC wrt d_A
1/bpr_A         = svp_A - svf_A*nu;                                                                                      // Equation (15): FOC wrt bpr_A
svf_A           = dF + beta*svf_A(+1)*(1-xi);                                                                            // Equation (16): FOC wrt F_A

bpr_A           = d_A/(1+rd_A)+b_AB/(1+rb_BA)+b_AC/(1+rb_CA)+b_AD/(1+rb_DA)
                  -(l_AB/(1+rb_AB)+l_AC/(1+rb_AC)+l_AD/(1+rb_AD))-x_A/(1+rl_A)
                  +del_BA*l_AB(-1)+del_CA*l_AC(-1)+del_DA*l_AD(-1)+alp_A*x_A(-1)
                  -(del_AB*b_AB(-1)+del_AC*b_AC(-1)+del_AD*b_AD(-1))-d_A(-1)
                  +dpb*(1-del_BA(-1))*l_AB(-2)+dpb*(1-del_CA(-1))*l_AC(-2)+dpb*(1-del_DA(-1))*l_AD(-2)
                  +dpf*(1-alp_A(-1))*x_A(-2)
                  -omegab/2*(((1-del_AB(-1))*b_AB(-2))^2+((1-del_AC(-1))*b_AC(-2))^2+((1-del_AD(-1))*b_AD(-2))^2)
                  +Gamma_A*S_ss;   
                                                                                                                        // Equation (17): Bank A Profits
F_A             = (1-xi)*F_A(-1)+nu*bpr_A;                                                                               // Equation (18): LoM of F_A

// B

svp_B*b_BA(-1)  = beta*svp_B(+1)*omegab*(1-del_BA)*(b_BA(-1))^2 + db;                                                    // Equation (10): FOC wrt del_BA
svp_B*b_BC(-1)  = beta*svp_B(+1)*omegab*(1-del_BC)*(b_BC(-1))^2 + db;                                                    // Equation (10): FOC wrt del_BC
svp_B*b_BD(-1)  = beta*svp_B(+1)*omegab*(1-del_BD)*(b_BD(-1))^2 + db;                                                    // Equation (10): FOC wrt del_BD

svp_B/(1+rb_BA) = beta*svp_B(+1)*del_AB(+1)+beta^2*svp_B(+2)*dpb*(1-del_AB(+1))-dF*k*rwb*(del_ss/del_AB(+1))^etab;       // Equation (11): FOC wrt l_BA
svp_B/(1+rb_BC) = beta*svp_B(+1)*del_CB(+1)+beta^2*svp_B(+2)*dpb*(1-del_CB(+1))-dF*k*rwb*(del_ss/del_CB(+1))^etab;       // Equation (11): FOC wrt l_BC 
svp_B/(1+rb_BD) = beta*svp_B(+1)*del_DB(+1)+beta^2*svp_B(+2)*dpb*(1-del_DB(+1))-dF*k*rwb*(del_ss/del_DB(+1))^etab;       // Equation (11): FOC wrt l_BD 
 
svp_B/(1+rb_AB) = beta*svp_B(+1)*del_BA(+1)+beta^2*svp_B(+2)*omegab*(1-del_BA(+1))^2*(b_BA);                             // Equation (12): FOC wrt b_BA
svp_B/(1+rb_CB) = beta*svp_B(+1)*del_BC(+1)+beta^2*svp_B(+2)*omegab*(1-del_BC(+1))^2*(b_BC);                             // Equation (12): FOC wrt b_BC
svp_B/(1+rb_DB) = beta*svp_B(+1)*del_BD(+1)+beta^2*svp_B(+2)*omegab*(1-del_BD(+1))^2*(b_BD);                             // Equation (12): FOC wrt b_BD

svp_B/(1+rl_B)  = beta*svp_B(+1)*alp_B(+1)+beta^2*svp_B(+2)*dpf*(1-alp_B(+1))-dF*k*rwf*(alp_ss/alp_B(+1))^etaf;          // Equation (13): FOC wrt x_B
svp_B/(1+rd_B)  = beta*svp_B(+1);                                                                                        // Equation (14): FOC wrt d_B
1/bpr_B         = svp_B - svf_B*nu;                                                                                      // Equation (15): FOC wrt bpr_B
svf_B           = dF + beta*svf_B(+1)*(1-xi);                                                                            // Equation (16): FOC wrt F_B

bpr_B           = d_B/(1+rd_B)+b_BA/(1+rb_AB)+b_BC/(1+rb_CB)+b_BD/(1+rb_DB)
                  -(l_BA/(1+rb_BA)+l_BC/(1+rb_BC)+l_BD/(1+rb_BD))-x_B/(1+rl_B)
                  +del_AB*l_BA(-1)+del_CB*l_BC(-1)+del_DB*l_BD(-1)+alp_B*x_B(-1)
                  -(del_BA*b_BA(-1)+del_BC*b_BC(-1)+del_BD*b_BD(-1))-d_B(-1)
                  +dpb*(1-del_AB(-1))*l_BA(-2)+dpb*(1-del_CB(-1))*l_BC(-2)+dpb*(1-del_DB(-1))*l_BD(-2)
                  +dpf*(1-alp_B(-1))*x_B(-2)
                  -omegab/2*(((1-del_BA(-1))*b_BA(-2))^2+((1-del_BC(-1))*b_BC(-2))^2+((1-del_BD(-1))*b_BD(-2))^2)
                  +Gamma_B*S_ss;   
                                                                                                                        // Equation (17): Bank B Profits
F_B             = (1-xi)*F_B(-1)+nu*bpr_B;                                                                               // Equation (18): LoM of F_B


// C

svp_C*b_CA(-1)  = beta*svp_C(+1)*omegab*(1-del_CA)*(b_CA(-1))^2 + db;                                                    // Equation (10): FOC wrt del_CA
svp_C*b_CB(-1)  = beta*svp_C(+1)*omegab*(1-del_CB)*(b_CB(-1))^2 + db;                                                    // Equation (10): FOC wrt del_CB
svp_C*b_CD(-1)  = beta*svp_C(+1)*omegab*(1-del_CD)*(b_CD(-1))^2 + db;                                                    // Equation (10): FOC wrt del_CD

svp_C/(1+rb_CA) = beta*svp_C(+1)*del_AC(+1)+beta^2*svp_C(+2)*dpb*(1-del_AC(+1))-dF*k*rwb*(del_ss/del_AC(+1))^etab;       // Equation (11): FOC wrt l_CA
svp_C/(1+rb_CB) = beta*svp_C(+1)*del_BC(+1)+beta^2*svp_C(+2)*dpb*(1-del_BC(+1))-dF*k*rwb*(del_ss/del_BC(+1))^etab;       // Equation (11): FOC wrt l_CB 
svp_C/(1+rb_CD) = beta*svp_C(+1)*del_DC(+1)+beta^2*svp_C(+2)*dpb*(1-del_DC(+1))-dF*k*rwb*(del_ss/del_DC(+1))^etab;       // Equation (11): FOC wrt l_CD 
 
svp_C/(1+rb_AC) = beta*svp_C(+1)*del_CA(+1)+beta^2*svp_C(+2)*omegab*(1-del_CA(+1))^2*(b_CA);                             // Equation (12): FOC wrt b_CA
svp_C/(1+rb_BC) = beta*svp_C(+1)*del_CB(+1)+beta^2*svp_C(+2)*omegab*(1-del_CB(+1))^2*(b_CB);                             // Equation (12): FOC wrt b_CB
svp_C/(1+rb_DC) = beta*svp_C(+1)*del_CD(+1)+beta^2*svp_C(+2)*omegab*(1-del_CD(+1))^2*(b_CD);                             // Equation (12): FOC wrt b_CD

svp_C/(1+rl_C)  = beta*svp_C(+1)*alp_C(+1)+beta^2*svp_C(+2)*dpf*(1-alp_C(+1))-dF*k*rwf*(alp_ss/alp_C(+1))^etaf;          // Equation (13): FOC wrt x_C
svp_C/(1+rd_C)  = beta*svp_C(+1);                                                                                        // Equation (14): FOC wrt d_C
1/bpr_C         = svp_C - svf_C*nu;                                                                                      // Equation (15): FOC wrt bpr_C
svf_C           = dF + beta*svf_B(+1)*(1-xi);                                                                            // Equation (16): FOC wrt F_C

bpr_C           = d_C/(1+rd_C)+b_CA/(1+rb_AC)+b_CB/(1+rb_BC)+b_CD/(1+rb_DC)
                  -(l_CA/(1+rb_CA)+l_CB/(1+rb_CB)+l_CD/(1+rb_CD))-x_C/(1+rl_C)
                  +del_AC*l_CA(-1)+del_BC*l_CB(-1)+del_DC*l_CD(-1)+alp_C*x_C(-1)
                  -(del_CA*b_CA(-1)+del_CB*b_CB(-1)+del_CD*b_CD(-1))-d_C(-1)
                  +dpb*(1-del_AC(-1))*l_CA(-2)+dpb*(1-del_BC(-1))*l_CB(-2)+dpb*(1-del_DC(-1))*l_CD(-2)
                  +dpf*(1-alp_C(-1))*x_C(-2)
                  -omegab/2*(((1-del_CA(-1))*b_CA(-2))^2+((1-del_CB(-1))*b_CB(-2))^2+((1-del_CD(-1))*b_CD(-2))^2)
                  +Gamma_C*S_ss;   
                                                                                                                        // Equation (17): Bank C Profits
F_C             = (1-xi)*F_C(-1)+nu*bpr_C;                                                                               // Equation (18): LoM of F_C

// D

svp_D*b_DA(-1)  = beta*svp_D(+1)*omegab*(1-del_DA)*(b_DA(-1))^2 + db;                                                    // Equation (10): FOC wrt del_DA
svp_D*b_DB(-1)  = beta*svp_D(+1)*omegab*(1-del_DB)*(b_DB(-1))^2 + db;                                                    // Equation (10): FOC wrt del_DB
svp_D*b_DC(-1)  = beta*svp_D(+1)*omegab*(1-del_DC)*(b_DC(-1))^2 + db;                                                    // Equation (10): FOC wrt del_DC

svp_D/(1+rb_DA) = beta*svp_D(+1)*del_AD(+1)+beta^2*svp_D(+2)*dpb*(1-del_AD(+1))-dF*k*rwb*(del_ss/del_AD(+1))^etab;       // Equation (11): FOC wrt l_DA
svp_D/(1+rb_DB) = beta*svp_D(+1)*del_BD(+1)+beta^2*svp_D(+2)*dpb*(1-del_BD(+1))-dF*k*rwb*(del_ss/del_BD(+1))^etab;       // Equation (11): FOC wrt l_DB 
svp_D/(1+rb_DC) = beta*svp_D(+1)*del_CD(+1)+beta^2*svp_D(+2)*dpb*(1-del_CD(+1))-dF*k*rwb*(del_ss/del_CD(+1))^etab;       // Equation (11): FOC wrt l_DC 
 
svp_D/(1+rb_AD) = beta*svp_D(+1)*del_DA(+1)+beta^2*svp_D(+2)*omegab*(1-del_DA(+1))^2*(b_DA);                             // Equation (12): FOC wrt b_DA
svp_D/(1+rb_BD) = beta*svp_D(+1)*del_DB(+1)+beta^2*svp_D(+2)*omegab*(1-del_DB(+1))^2*(b_DB);                             // Equation (12): FOC wrt b_DB
svp_D/(1+rb_CD) = beta*svp_D(+1)*del_DC(+1)+beta^2*svp_D(+2)*omegab*(1-del_DC(+1))^2*(b_DC);                             // Equation (12): FOC wrt b_DC

svp_D/(1+rl_D)  = beta*svp_D(+1)*alp_D(+1)+beta^2*svp_D(+2)*dpf*(1-alp_D(+1))-dF*k*rwf*(alp_ss/alp_D(+1))^etaf;          // Equation (13): FOC wrt x_D
svp_D/(1+rd_D)  = beta*svp_D(+1);                                                                                        // Equation (14): FOC wrt d_D
1/bpr_D         = svp_D - svf_D*nu;                                                                                      // Equation (15): FOC wrt bpr_D
svf_D           = dF + beta*svf_D(+1)*(1-xi);                                                                            // Equation (16): FOC wrt F_D

bpr_D           = d_D/(1+rd_D)+b_DA/(1+rb_AD)+b_DB/(1+rb_BD)+b_DC/(1+rb_CD)
                  -(l_DA/(1+rb_DA)+l_DB/(1+rb_DB)+l_DC/(1+rb_DC))-x_D/(1+rl_D)
                  +del_AD*l_DA(-1)+del_BD*l_DB(-1)+del_CD*l_DC(-1)+alp_D*x_D(-1)
                  -(del_DA*b_DA(-1)+del_DB*b_DB(-1)+del_DC*b_DC(-1))-d_D(-1)
                  +dpb*(1-del_AD(-1))*l_DA(-2)+dpb*(1-del_BD(-1))*l_DB(-2)+dpb*(1-del_CD(-1))*l_DC(-2)
                  +dpf*(1-alp_D(-1))*x_D(-2)
                  -omegab/2*(((1-del_DA(-1))*b_DA(-2))^2+((1-del_DB(-1))*b_DB(-2))^2+((1-del_DC(-1))*b_DC(-2))^2)
                  +Gamma_D*S_ss;   
                                                                                                                        // Equation (17): Bank D Profits
F_D             = (1-xi)*F_D(-1)+nu*bpr_D;                                                                               // Equation (18): LoM of F_C

/*----------------------------------------------------------------------------------------
Firm equations: 4 firms x 7 equations = 28 equations
----------------------------------------------------------------------------------------*/

@#for region in ["A","B","C","D"]

    x_@{region}/(1+rl_@{region}) = k_@{region}-(1-tau)*k_@{region}(-1);                            // Equation (46): LoM of k

    fpr_@{region} = tfp_@{region}*y_@{region}-w_@{region}*n_@{region}
        -alp_@{region}*x_@{region}(-1)-(omegaf/2)*(1-alp_@{region}(-1)*x_@{region}(-2))^2;         // Equation (47): Firm profits

    y_@{region} = tfp_@{region}*k_@{region}^mu*n_@{region}^(1-mu);                                 // Equation (48): Production function

    tfp_@{region}*(1-mu)*(k_@{region}/n_@{region})^mu=w_@{region};                                 // Equation (49): FOC wrt n

    tfp_@{region}*mu*(k_@{region}/n_@{region})^(mu-1) 
        = svk_@{region} - beta*(1-tau)*svk_@{region}(+1);                                          // Equation (50): FOC wrt k

    svk_@{region}/(1+rl_@{region}) 
        = beta*alp_@{region}(+1)+beta^2*omegaf*(1-alp_@{region}(+1))^2*x_@{region};                // Equation (51): FOC wrt x

    x_@{region}(-1) = beta*omegaf*(1-alp_@{region})*(x_@{region}(-1))^2+df;                        // Equation (52): FOC wrt alp

@#endfor

/*----------------------------------------------------------------------------------------
Central bank and government: 19 equations
----------------------------------------------------------------------------------------*/

tot_l = l_AC+l_BA+l_CD+l_DB+l_AB+l_AD+l_BC+l_BD+l_CA+l_CB+l_DA+l_DC;
tot_b = b_CA+b_AB+b_DC+b_BD+b_BA+b_DA+b_CB+b_DB+b_AC+b_BC+b_AD+b_CD;

tot_y = y_A+y_B+y_C+y_D;
tot_x = x_A+x_B+x_C+x_D;
                                
T_A*y_A+xi*F_A(-1) = dpb*(1-del_AB(-1))*b_AB(-2)+dpb*(1-del_AC(-1))*b_AC(-2)+dpb*(1-del_AD(-1))*b_AD(-2)+(1-alp_A(-1))*x_A(-2);
T_B*y_B+xi*F_B(-1) = dpb*(1-del_BA(-1))*b_BA(-2)+dpb*(1-del_BC(-1))*b_BC(-2)+dpb*(1-del_BD(-1))*b_BD(-2)+(1-alp_B(-1))*x_B(-2);
T_C*y_C+xi*F_C(-1) = dpb*(1-del_CA(-1))*b_CA(-2)+dpb*(1-del_CB(-1))*b_CB(-2)+dpb*(1-del_CD(-1))*b_CD(-2)+(1-alp_C(-1))*x_C(-2);
T_D*y_D+xi*F_D(-1) = dpb*(1-del_DA(-1))*b_DA(-2)+dpb*(1-del_DB(-1))*b_DB(-2)+dpb*(1-del_DC(-1))*b_DC(-2)+(1-alp_D(-1))*x_D(-2);


/*----------------------------------------------------------------------------------------
Spread: Interbank lending and borrowing rate
----------------------------------------------------------------------------------------*/

IBspread_A  = rb_AC - rb_BA;        // Symmetric case 1: A as intermediary (b/w C and B) 
IBspread_D  = rb_DB - rb_CD;        // Symmetric case 2: D as intermediary (b/w B and C)

IBspread_B = rb_BA - rb_DB;         // Asymmetric case: B as intermediary (b/w A and D)
IBspread_C = rb_CD - rb_AC;         // Asymmetric case: B as intermediary (b/w D and A)

IBspread_AB = rb_AB - rb_BA;        // Reciprocity (A and B)
IBspread_BD = rb_BD - rb_DB;        // Reciprocity (B and D)

/*----------------------------------------------------------------------------------------
Spread: Bank lending rate and interbank borrowing rate
----------------------------------------------------------------------------------------*/

//Lspread_A  = rl_A - rb_BA;
//Lspread_B1 = rl_B - rb_AB;
//Lspread_B2 = rl_B - rb_CB;

/*----------------------------------------------------------------------------------------
Spread: Interbank lending rate and deposit rate
----------------------------------------------------------------------------------------*/

//Dspread_A  = rb_AB - rd_A;
//Dspread_B1 = rb_BA - rd_B;
//Dspread_B2 = rb_BC - rd_B;

/*----------------------------------------------------------------------------------------
Net interbank exposures
----------------------------------------------------------------------------------------*/

//NetIBexp_A = l_AC - b_AB;           // Symmetric case 1: A as intermediary (b/w C and B) 
//NetIBexp_D = l_DB - b_DC;            // Symmetric case 2: D as intermediary (b/w B and C)

//NetIBexp_B1 = l_BA - b_BD;           // Asymmetric case: B as intermediary (b/w A and D)
//NetIBexp_B2 = l_BD - b_BA;           // Asymmetric case: B as intermediary (b/w D and A)

//NetIBexp_AB = l_AB - b_AB;           // Reciprocity (A and B)
//NetIBexp_BD = l_BD - b_BD;           // Reciprocity (B and D)
                                                  
/*----------------------------------------------------------------------------------------
Closing the model: 20 equations 
----------------------------------------------------------------------------------------*/
// Totals 

tot_bpr = bpr_A+bpr_B+bpr_C+bpr_D;
tot_fpr = fpr_A+fpr_B+fpr_C+fpr_D;

rb_av = (rb_AB+rb_AC+rb_AD+rb_BA+rb_BC+rb_BD+rb_CA+rb_CB+rb_CD+rb_DA+rb_DB+rb_DC)/num_links;
rd_av = (rd_A+rd_B+rd_C+rd_D)/num_banks;
rl_av = (rl_A+rl_B+rl_C+rl_D)/num_banks;

ratespread = rl_av - rd_av;

// Shocks

@#include "liqinjvareqs_comp.mod"
@#include "shockvareqs.mod"
     
end;

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Initialise variables at their steady state
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

initval;

d_A = d_ss; d_B = d_ss; d_C = d_ss; d_D = d_ss; 
c_A = c_ss; c_B = c_ss; c_C = c_ss; c_D = c_ss; 
w_A = w_ss; w_B = w_ss; w_C = w_ss;  w_D = w_ss; 
n_A = n_ss; n_B = n_ss; n_C = n_ss; n_D = n_ss;           

k_A = k_ss; k_B = k_ss; k_C = k_ss; k_D = k_ss;
y_A = y_ss; y_B = y_ss; y_C = y_ss; y_D = y_ss;
x_A = x_ss; x_B = x_ss; x_C = x_ss; x_D = x_ss;
alp_A = alp_ss; alp_B = alp_ss; alp_C = alp_ss; alp_D = alp_ss;
fpr_A = fpr_ss; fpr_B = fpr_ss; fpr_C = fpr_ss; fpr_D = fpr_ss;
svk_A = svk_ss; svk_B = svk_ss; svk_C = svk_ss; svk_D = svk_ss;                   

del_AB = del_ss; del_AC = del_ss; del_AD = del_ss; 
del_BA = del_ss; del_BC = del_ss; del_BD = del_ss;  
del_CA = del_ss; del_CB = del_ss; del_CD = del_ss; 
del_DA = del_ss; del_DB = del_ss; del_DC = del_ss;  

l_AB = l_ss; l_AC = l_ss; l_AD = l_ss;
l_BA = l_ss; l_BC = l_ss; l_BD = l_ss;
l_CA = l_ss; l_CB = l_ss; l_CD = l_ss;
l_DA = l_ss; l_DB = l_ss; l_DC = l_ss;

b_AB = b_ss; b_AC = b_ss; b_AD = b_ss;
b_BA = b_ss; b_BC = b_ss; b_BD = b_ss;
b_CA = b_ss; b_CB = b_ss; b_CD = b_ss;
b_DA = b_ss; b_DB = b_ss; b_DC = b_ss;

       
F_A = F_ss; F_B = F_ss; F_C = F_ss; F_D = F_ss;
bpr_A = bpr_ss; bpr_B = bpr_ss; bpr_C = bpr_ss; bpr_D = bpr_ss;
svp_A = svp_ss; svp_B = svp_ss; svp_C = svp_ss; svp_D = svp_ss;
svf_A = svf_ss; svf_B = svf_ss; svf_C = svf_ss; svf_D  = svf_ss;
                               
rd_A = rd_ss; rd_B = rd_ss; rd_C = rd_ss; rd_D = rd_ss;

rb_AB = rb_ss; rb_AC = rb_ss; rb_AD = rb_ss;
rb_BA = rb_ss; rb_BC = rb_ss; rb_BD = rb_ss;
rb_CA = rb_ss; rb_CB = rb_ss; rb_CD = rb_ss;
rb_DA = rb_ss; rb_DB = rb_ss; rb_DC = rb_ss; 
  
rl_A = rl_ss; rl_B = rl_ss; rl_C = rl_ss; rl_D = rl_ss; 
              

T_A = T_ss; T_B = T_ss; T_C = T_ss; T_D = T_ss;
tot_bpr = 4*bpr_ss; tot_fpr = 4*fpr_ss;
tot_y = 4*y_ss; tot_l = 4*l_ss;  tot_b = 4*b_ss;
tot_x = 4*x_ss; rl_av = rl_ss; rd_av = rd_ss;  rb_av = rb_ss;                    

IBspread_A  = 0; IBspread_D  = 0; IBspread_B =0; IBspread_C = 0; IBspread_AB = 0; IBspread_BD = 0;
ratespread = rl_ss - rd_ss;        
//Lspread_A  = 0; Lspread_B1 = 0; Lspread_B2 = 0; Dspread_A  = 0; Dspread_B1 = 0; Dspread_B2 = 0;
//NetIBexp_A = 0; NetIBexp_D = 0; NetIBexp_B1 = 0; NetIBexp_B2 = 0; NetIBexp_AB = 0; NetIBexp_BD = 0;                                             

eps_tfp   = 0;
eps_Gamma = 0;

@#include "liqinjvarinit_comp.mod"
@#include "shockvarinit.mod"                                                                                         

end;

/*%%%%%%%%%%%%%%%%%%%%%%%%%
Compute the steady state
%%%%%%%%%%%%%%%%%%%%%%%%*/

steady(solve_algo=0);

/*%%%%%%%%%%%%%%%%%%%%%%%%%
Check the BK conditions
%%%%%%%%%%%%%%%%%%%%%%%%%*/

check;

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Define magnitude of deterministic productivity and banking shocks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
shocks;

var eps_Gamma; stderr 0.1;
var eps_tfp;   stderr 0.1;

@#if liqinj == "on"
    var eps_m;     stderr 0.1;
@#endif 

end;

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// produce stochastic simulation for variables defined after the command
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

options_.nograph=1;
stoch_simul(order=1,irf=100,periods=5000) 

IBspread_A IBspread_D IBspread_B IBspread_C IBspread_AB IBspread_BD           //1-6
rl_A rl_B rl_C rl_D                                                           //7-10
rd_A rd_B rd_C rd_D                                                           //11-14
l_AB l_BA l_BD                                                                //15-17  
b_AB b_BA b_BD                                                                //18-20
x_A x_B x_C x_D                                                               //21-24
d_A d_B d_C d_D                                                               //25-28   
del_AB del_BA del_BD                                                          //29-31
alp_A alp_B alp_C alp_D                                                       //32-35
tot_y tot_x ratespread rb_av                                                  //36-39 
@#if liqinj == "on"
     m
@#endif 
rb_AC rb_BA rb_CD;                                                            //40-42

     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preparing Dynare output for use in MATLAB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

shockconfig = 'SC2';
liqinj      = 'on';

if strcmp(shockconfig,'SC1') && strcmp(liqinj,'off')

    irfs.SC1 = oo_.irfs;
    save('IRF_compnet_noLI.mat','-struct','irfs','SC1');
    save('varscompnet_noLI.mat','var_list_');

elseif strcmp(shockconfig,'SC1') &&  strcmp(liqinj,'on')

    irfs.SC1 = oo_.irfs;
    save('IRF_compnet_LI.mat','-struct','irfs','SC1');
    save('varscompnet_LI.mat','var_list_');

elseif strcmp(shockconfig,'SC2') &&  strcmp(liqinj,'off')

    irfs.SC2 = oo_.irfs
    save('IRF_compnet_noLI.mat','-struct','irfs','SC2','-append');

elseif strcmp(shockconfig,'SC2') &&  strcmp(liqinj,'on')

    irfs.SC2 = oo_.irfs
    save('IRF_compnet_LI.mat','-struct','irfs','SC2','-append');

end

