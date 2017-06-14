
@#if  CPshockconfig == "SC1"
    Gamma_A Gamma_B Gamma_C Gamma_D; 
@#else
    @#if CPshockconfig == "SC2"
        Gamma_B Gamma_C Gamma_D;
    @#else
        @#if CPshockconfig == "SC3"
          Gamma_A Gamma_C Gamma_D;
        @#else
            @#if CPshockconfig == "SC4"
                Gamma_A Gamma_B Gamma_D;
            @#endif
        @#endif
    @#endif
@#endif
