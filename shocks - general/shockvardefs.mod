
@#if shockconfig == "SC1"
    tfp_A tfp_B tfp_C tfp_D; 
@#else
    @#if shockconfig == "SC2"
        tfp_A tfp_B tfp_C tfp_D Gamma_A;
    @#endif
@#endif
