@#if  CPshockconfig == "SC1"
    Gamma_A = Gamma_ss; Gamma_B = Gamma_ss; Gamma_C = Gamma_ss; Gamma_D = Gamma_ss;
@#else
    @#if CPshockconfig == "SC2"
        Gamma_B = Gamma_ss; Gamma_C = Gamma_ss; Gamma_D = Gamma_ss;
    @#else
        @#if CPshockconfig == "SC3"
          Gamma_A = Gamma_ss; Gamma_C = Gamma_ss; Gamma_D = Gamma_ss;
        @#else
            @#if CPshockconfig == "SC4"
                Gamma_A = Gamma_ss; Gamma_B = Gamma_ss; Gamma_D = Gamma_ss;
            @#endif
        @#endif
    @#endif
@#endif
