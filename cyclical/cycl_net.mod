/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Model under cyclical network structure  %
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

@#define shockconfig = "SC2"
@#define liqinj      = "on"

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Declaration of variables: Total = 99
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

var d_A d_B d_C d_D c_A c_B c_C c_D w_A w_B w_C w_D n_A n_B n_C n_D            // Household variables (16)

k_A k_B k_C k_D y_A y_B y_C y_D x_A x_B x_C x_D alp_A alp_B alp_C alp_D
fpr_A fpr_B fpr_C fpr_D svk_A svk_B svk_C svk_D                                // Firm variables      (24)

del_AB del_BD del_CA del_DC l_AC l_BA l_CD l_DB b_AB b_BD b_CA b_DC          
F_A F_B F_C F_D  bpr_A bpr_B bpr_C bpr_D
svp_A svp_B svp_C svp_D svf_A svf_B svf_C svf_D                                // Bank variables      (28)

rd_A rd_B rd_C rd_D rb_AC rb_BA rb_CD rb_DB rl_A rl_B rl_C rl_D rd_av rl_av rb_av        // Interest rates      (13)

IBspread_A IBspread_B IBspread_C IBspread_D 
LDspread_A LDspread_B LDspread_C LDspread_D
ratespread
//Lspread_A  Lspread_B  Lspread_C Lspread_D  
//Dspread_A  Dspread_B  Dspread_C Dspread_D
//NetIBexp_A NetIBexp_B NetIBexp_C NetIBexp_D

T_A T_B T_C T_D tot_bpr tot_fpr tot_y tot_x  tot_l tot_b        // Other              (14)

@#include "liqinjvardefs_cycl.mod"
@#include "shockvardefs.mod"                                                    

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Declaration of stochastic shocks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

@#include "liqinjvarexo.mod"


/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Declaration of parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

parameters chi beta dpb dpf rwf rwb rws nu etaf etab tau  mu k ykr num_links num_banks
rhob_Gamma rho_Gamma rho_tfp alp_ss del_ss d_ss n_ss S_ss rb_ss rl_ss rd_ss         // Calibration: Parameters and steady states

omegab db dF xi l_ss b_ss bpr_ss F_ss Gamma_ss svf_ss svp_ss omegaf df
tfp_ss  k_ss x_ss y_ss w_ss svk_ss fpr_ss T_ss c_ss                          // Implied steady state parameters 

ratioc2y ratiok2y ratiot2y

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

rb_ss = 0.01;
rl_ss = 0.005;
rd_ss = 0.005;

beta     = 1/(1+rd_ss);

k_ss = n_ss/(ykr)^(1/(1-mu));
x_ss = tau*k_ss*(1+rl_ss);
S_ss = 1*x_ss;

tfp_ss   = 1;
Gamma_ss = rhob_Gamma*S_ss; 

@#include "liqinjpareqs.mod"
@#include "shockpareqs.mod"

num_links = 4;
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

bpr_ss = (1/(1+rd_ss)-1)*d_ss+(1/(1+rb_ss)-del_ss)*b_ss
        +(alp_ss+dpf*(1-alp_ss)-1/(1+rl_ss))*x_ss
        + (del_ss+dpb*(1-del_ss)-1/(1+rb_ss))*l_ss
        -(omegab/2)*((1-del_ss)*b_ss)^2+Gamma_ss;  

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

c_ss = w_ss*n_ss+(1-1/(1+rd_ss))*d_ss-xi*F_ss-dpf*(1-alp_ss)*x_ss-dpb*(1-del_ss)*l_ss;

/*----------------------------------------------------------------------------------------
Steady-state ratios
----------------------------------------------------------------------------------------*/

ratioc2y = c_ss/y_ss;
ratiok2y = k_ss/y_ss;
ratiot2y = T_ss/y_ss;

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Model equations: Total = 99
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

model;

/*----------------------------------------------------------------------------------------
Household equations: 4 households x 3 equations = 12 equations
----------------------------------------------------------------------------------------*/

@#for region in ["A","B","C","D"]
    T_@{region}*y_@{region}+c_@{region}+d_@{region}/(1+rd_@{region})          
        = w_@{region}*n_@{region}+d_@{region}(-1);                                      // Equation (): Budget constraint
    (1/c_@{region})*(1/(1+rd_@{region})) 
        = beta/c_@{region}(+1) - chi*(d_@{region}/(1+rd_@{region})-d_ss/(1+rd_ss));     // Equation (): Euler equation with deposits FOC plugged in 
    n_@{region} = n_ss;                                                                 // Equation (): Exogenous labour supply
@#endfor

//chi*(d_@{region}/(1+rd_@{region})-d_ss/(1+rd_ss))
/*----------------------------------------------------------------------------------------
Bank equations: 4 banks x 9 equations = 36 equations
----------------------------------------------------------------------------------------*/

// A

svp_A*b_AB(-1)  = beta*svp_A(+1)*omegab*(1-del_AB)*(b_AB(-1))^2 + db;                                                    // Equation (10): FOC wrt del_AB
svp_A/(1+rb_AC) = beta*svp_A(+1)*del_CA(+1)+beta^2*svp_A(+2)*dpb*(1-del_CA(+1))-dF*k*rwb*(del_ss/del_CA(+1))^etab;       // Equation (11): FOC wrt l_AC 
svp_A/(1+rb_BA) = beta*svp_A(+1)*del_AB(+1)+beta^2*svp_A(+2)*omegab*(1-del_AB(+1))^2*(b_AB);                             // Equation (12): FOC wrt b_AB
svp_A/(1+rl_A)  = beta*svp_A(+1)*alp_A(+1)+beta^2*svp_A(+2)*dpf*(1-alp_A(+1))-dF*k*rwf*(alp_ss/alp_A(+1))^etaf;          // Equation (13): FOC wrt x_A
svp_A/(1+rd_A)  = beta*svp_A(+1);                                                                                        // Equation (14): FOC wrt d_A
1/bpr_A         = svp_A - svf_A*nu;                                                                                      // Equation (15): FOC wrt bpr_A
svf_A           = dF + beta*svf_A(+1)*(1-xi);                                                                            // Equation (16): FOC wrt F_A
bpr_A           = d_A/(1+rd_A)+b_AB/(1+rb_BA)-l_AC/(1+rb_AC)-x_A/(1+rl_A)+del_CA*l_AC(-1)+alp_A*x_A(-1)
                    -del_AB*b_AB(-1)-d_A(-1)+dpb*(1-del_CA(-1))*l_AC(-2)+dpf*(1-alp_A(-1))*x_A(-2)
                    -omegab/2*((1-del_AB(-1))*b_AB(-2))^2+Gamma_A*S_ss;                                                  // Equation (17): Bank A Profits
F_A             = (1-xi)*F_A(-1)+nu*bpr_A;                                                                               // Equation (18): LoM of F_A

// B

svp_B*b_BD(-1)  = beta*svp_B(+1)*omegab*(1-del_BD)*(b_BD(-1))^2 + db;                                                    // Equation (19): FOC wrt del_BD
svp_B/(1+rb_BA) = beta*svp_B(+1)*del_AB(+1)+beta^2*svp_B(+2)*dpb*(1-del_AB(+1))-dF*k*rwb*(del_ss/del_AB(+1))^etab;       // Equation (20): FOC wrt l_BA 
svp_B/(1+rb_DB) = beta*svp_B(+1)*del_BD(+1)+beta^2*svp_B(+2)*omegab*(1-del_BD(+1))^2*(b_BD);                             // Equation (21): FOC wrt b_BD
svp_B/(1+rl_B)  = beta*svp_B(+1)*alp_B(+1)+beta^2*svp_B(+2)*dpf*(1-alp_B(+1))-dF*k*rwf*(alp_ss/alp_B(+1))^etaf;          // Equation (22): FOC wrt x_B
svp_B/(1+rd_B)  = beta*svp_B(+1);                                                                                        // Equation (23): FOC wrt d_B
1/bpr_B         = svp_B - svf_B*nu;                                                                                      // Equation (24): FOC wrt bpr_B
svf_B           = dF + beta*svf_B(+1)*(1-xi);                                                                            // Equation (25): FOC wrt F_B
bpr_B           = d_B/(1+rd_B)+b_BD/(1+rb_DB)-l_BA/(1+rb_BA)-x_B/(1+rl_B)+del_AB*l_BA(-1)+alp_B*x_B(-1)
                    -del_BD*b_BD(-1)-d_B(-1)+dpb*(1-del_AB(-1))*l_BA(-2)+dpf*(1-alp_B(-1))*x_B(-2)
                    -omegab/2*((1-del_BD(-1))*b_BD(-2))^2+Gamma_B*S_ss;                                                  // Equation (26): Bank B Profits
F_B             = (1-xi)*F_B(-1)+nu*bpr_B;                                                                               // Equation (27): LoM of F_B

// C

svp_C*b_CA(-1)  = beta*svp_C(+1)*omegab*(1-del_CA)*(b_CA(-1))^2 + db;                                                    // Equation (28): FOC wrt del_CA
svp_C/(1+rb_CD) = beta*svp_C(+1)*del_DC(+1)+beta^2*svp_C(+2)*dpb*(1-del_DC(+1))-dF*k*rwb*(del_ss/del_DC(+1))^etab;       // Equation (29): FOC wrt l_CD 
svp_C/(1+rb_AC) = beta*svp_C(+1)*del_CA(+1)+beta^2*svp_C(+2)*omegab*(1-del_CA(+1))^2*(b_CA);                             // Equation (30): FOC wrt b_CA
svp_C/(1+rl_C)  = beta*svp_C(+1)*alp_C(+1)+beta^2*svp_C(+2)*dpf*(1-alp_C(+1))-dF*k*rwf*(alp_ss/alp_C(+1))^etaf;          // Equation (31): FOC wrt x_C
svp_C/(1+rd_C)  = beta*svp_C(+1);                                                                                        // Equation (32): FOC wrt d_C
1/bpr_C         = svp_C - svf_C*nu;                                                                                      // Equation (33): FOC wrt bpr_C
svf_C           = dF + beta*svf_C(+1)*(1-xi);                                                                            // Equation (34): FOC wrt F_C
bpr_C           = d_C/(1+rd_C)+b_CA/(1+rb_AC)-l_CD/(1+rb_CD)-x_C/(1+rl_C)+del_DC*l_CD(-1)+alp_C*x_C(-1)
                    -del_CA*b_CA(-1)-d_C(-1)+dpb*(1-del_DC(-1))*l_CD(-2)+dpf*(1-alp_C(-1))*x_C(-2)
                    -omegab/2*((1-del_CA(-1))*b_CA(-2))^2+Gamma_C*S_ss;                                                  // Equation (35): Bank A Profits                                                                       // Equation (): FOC wrt F_C
F_C             = (1-xi)*F_C(-1)+nu*bpr_C;                                                                               // Equation (36): LoM of F_C

// D

svp_D*b_DC(-1)  = beta*svp_D(+1)*omegab*(1-del_DC)*(b_DC(-1))^2 + db;                                                    // Equation (37): FOC wrt del_DC
svp_D/(1+rb_DB) = beta*svp_D(+1)*del_BD(+1)+beta^2*svp_D(+2)*dpb*(1-del_BD(+1))-dF*k*rwb*(del_ss/del_BD(+1))^etab;       // Equation (38): FOC wrt l_DB 
svp_D/(1+rb_CD) = beta*svp_D(+1)*del_DC(+1)+beta^2*svp_D(+2)*omegab*(1-del_DC(+1))^2*(b_DC);                             // Equation (39): FOC wrt b_DC
svp_D/(1+rl_D)  = beta*svp_D(+1)*alp_D(+1)+beta^2*svp_D(+2)*dpf*(1-alp_D(+1))-dF*k*rwf*(alp_ss/alp_D(+1))^etaf;          // Equation (40): FOC wrt x_D
svp_D/(1+rd_D)  = beta*svp_D(+1);                                                                                        // Equation (41): FOC wrt d_D
1/bpr_D         = svp_D - svf_D*nu;                                                                                      // Equation (42): FOC wrt bpr_D
svf_D           = dF + beta*svf_D(+1)*(1-xi);                                                                            // Equation (43): FOC wrt F_D
bpr_D            = d_D/(1+rd_D)+b_DC/(1+rb_CD)-l_DB/(1+rb_DB)-x_D/(1+rl_D)+del_BD*l_DB(-1)+alp_D*x_D(-1)
                    -del_DC*b_DC(-1)-d_D(-1)+dpb*(1-del_BD(-1))*l_DB(-2)+dpf*(1-alp_D(-1))*x_D(-2)
                    -(omegab/2)*((1-del_DC(-1))*b_DC(-2))^2+Gamma_D*S_ss;                                                // Equation (44): Bank A Profits
F_D             = (1-xi)*F_D(-1)+nu*bpr_D;                                                                               // Equation (45): LoM of F_D

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
Central bank and government: 11 equations
----------------------------------------------------------------------------------------*/

tot_l = l_AC+l_BA+l_CD+l_DB;
tot_b = b_CA+b_AB+b_DC+b_BD;

tot_y = y_A+y_B+y_C+y_D;
tot_x = x_A+x_B+x_C+x_D;
                                
//T*tot_y+xi*(F_A(-1)+F_B(-1)+F_C(-1)+F_D(-1)) = dpb*((1-del_AB(-1))*l_BA(-2)+(1-del_BD(-1))*l_DB(-2)+(1-del_CA(-1))*l_AC(-2)+(1-del_DC(-1))*l_CD(-2))+
 //dpf*((1-alp_A(-1))*x_A(-2)+(1-alp_B(-1))*x_B(-2)+(1-alp_C(-1))*x_C(-2)+(1-alp_D(-1))*x_D(-2));

T_A*y_A+xi*F_A(-1) = dpb*(1-del_AB(-1))*b_AB(-2)+dpf*(1-alp_A(-1))*x_A(-2);
T_B*y_B+xi*F_B(-1) = dpb*(1-del_BD(-1))*b_BD(-2)+dpf*(1-alp_B(-1))*x_B(-2);
T_C*y_C+xi*F_C(-1) = dpb*(1-del_CA(-1))*b_CA(-2)+dpf*(1-alp_C(-1))*x_C(-2);
T_D*y_D+xi*F_D(-1) = dpb*(1-del_DC(-1))*b_DC(-2)+dpf*(1-alp_D(-1))*x_D(-2);
                                                  
/*----------------------------------------------------------------------------------------
Closing the model: 12 equations 
----------------------------------------------------------------------------------------*/
// Totals 

tot_bpr = bpr_A+bpr_B+bpr_C+bpr_D;
tot_fpr = fpr_A+fpr_B+fpr_C+fpr_D;

rb_av = (rb_AC+rb_BA+rb_CD+rb_DB)/num_links;
rd_av = (rd_A+rd_B+rd_C+rd_D)/num_banks;
rl_av = (rl_A+rl_B+rl_C+rl_D)/num_banks;

ratespread = rl_av - rd_av;

// Shocks

@#include "liqinjvareqs_cycl.mod"

// IB Spreads

IBspread_A = rb_AC - rb_BA;
IBspread_B = rb_BA - rb_DB;
IBspread_C = rb_CD - rb_AC;
IBspread_D = rb_DB - rb_CD;

LDspread_A = x_A - d_A;
LDspread_B = x_B - d_B;
LDspread_C = x_C - d_C;
LDspread_D = x_D - d_C;


@#include "shockvareqs.mod"
     
/*----------------------------------------------------------------------------------------
Log-linearisation (to get output in % deviations)
----------------------------------------------------------------------------------------*/

//@#include "cycl_ll_model.mod"

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

del_AB = del_ss; del_BD = del_ss; del_CA = del_ss; del_DC = del_ss;
l_AC = l_ss; l_BA = l_ss; l_CD = l_ss; l_DB = l_ss;
b_AB = b_ss; b_BD = b_ss; b_CA = b_ss; b_DC = b_ss;          
F_A = F_ss; F_B = F_ss; F_C = F_ss; F_D = F_ss;
bpr_A = bpr_ss; bpr_B = bpr_ss; bpr_C = bpr_ss; bpr_D = bpr_ss;
svp_A = svp_ss; svp_B = svp_ss; svp_C = svp_ss; svp_D = svp_ss;
svf_A = svf_ss; svf_B = svf_ss; svf_C = svf_ss; svf_D  = svf_ss;
                               
rd_A = rd_ss; rd_B = rd_ss; rd_C = rd_ss; rd_D = rd_ss;
rb_AC = rb_ss; rb_BA = rb_ss; rb_CD = rb_ss; rb_DB = rb_ss;
rl_A = rl_ss; rl_B = rl_ss; rl_C = rl_ss; rl_D = rl_ss; 
              
T_A = T_ss; T_B = T_ss; T_C = T_ss; T_D = T_ss;
tot_bpr = 4*bpr_ss; tot_fpr = 4*fpr_ss;
tot_y = 4*y_ss; tot_l = 4*l_ss;  tot_b = 4*b_ss; 
tot_x = 4*x_ss; rl_av = rl_ss; rd_av = rd_ss;  rb_av = rb_ss;      

IBspread_A = 0; IBspread_B = 0; IBspread_C = 0; IBspread_D = 0;
LDspread_A = x_ss - d_ss; LDspread_B = x_ss - d_ss;
LDspread_C = x_ss - d_ss; LDspread_D = x_ss - d_ss;
      
eps_tfp   = 0;
eps_Gamma = 0;

@#include "liqinjvarinit_cycl.mod"
@#include "shockvarinit.mod"  
//@#include "cycl_ll_init.mod"                                                                                       

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
//@#include "cycl_ll_vardefs.mod"

IBspread_A IBspread_B IBspread_C IBspread_D   //1-4
rl_A rl_B rl_C rl_D                           //5-8
rd_A rd_B rd_C rd_D                           //9-12
l_AC l_BA l_CD l_DB                           //13-16     
b_AB b_BD b_CA b_DC                           //17-20
x_A x_B x_C x_D                               //21-24
d_A d_B d_C d_D                               //25-28
del_AB del_BD del_CA del_DC                   //29-32
alp_A alp_B alp_C alp_D                       //33-36
tot_y tot_x ratespread rb_av                  //37-40
@#if liqinj == "on"
     m
@#endif                                       //Liquidity injection variable for results section 2
rb_AC rb_BA rb_CD rb_DB;                      //41-44 - Appendix
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preparing Dynare output for use in MATLAB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

shockconfig = 'SC2';
liqinj      = 'on';

if strcmp(shockconfig,'SC1') && strcmp(liqinj,'off')

    irfs.SC1 = oo_.irfs;
    save('IRF_cyclnet_noLI.mat','-struct','irfs','SC1');
    save('varscyclnet_noLI.mat','var_list_');

elseif strcmp(shockconfig,'SC1') &&  strcmp(liqinj,'on')

    irfs.SC1 = oo_.irfs;
    save('IRF_cyclnet_LI.mat','-struct','irfs','SC1');
    save('varscyclnet_LI.mat','var_list_');

elseif strcmp(shockconfig,'SC2') &&  strcmp(liqinj,'off')

    irfs.SC2 = oo_.irfs
    save('IRF_cyclnet_noLI.mat','-struct','irfs','SC2','-append');

elseif strcmp(shockconfig,'SC2') &&  strcmp(liqinj,'on')

    irfs.SC2 = oo_.irfs
    save('IRF_cyclnet_LI.mat','-struct','irfs','SC2','-append');

end

