
@#if shockconfig == "SC1"
    Gamma_A Gamma_B Gamma_C Gamma_D;
@#else
    @#if shockconfig == "SC2"
    Gamma_B Gamma_C Gamma_D; 
    @#endif
@#endif
