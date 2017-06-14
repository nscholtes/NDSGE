@#if shockconfig == "SC1"
    tfp_A = tfp_ss; tfp_B = tfp_ss; tfp_C= tfp_ss; tfp_D = tfp_ss; 
@#else
    @#if shockconfig == "SC2"
        Gamma_A = Gamma_ss; tfp_A = tfp_ss; tfp_B = tfp_ss; tfp_C= tfp_ss; tfp_D = tfp_ss;
    @#endif
@#endif
