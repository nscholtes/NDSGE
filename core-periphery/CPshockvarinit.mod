@#if CPshockconfig == "SC1"
    tfp_A = tfp_ss; tfp_B= tfp_ss; tfp_C= tfp_ss; tfp_D= tfp_ss; 
@#else
    @#if CPshockconfig == "SC2"
        tfp_A= tfp_ss; tfp_B= tfp_ss; tfp_C= tfp_ss; tfp_D= tfp_ss; Gamma_A= Gamma_ss;
    @#else
        @#if CPshockconfig== "SC3"
          tfp_A= tfp_ss; tfp_B= tfp_ss; tfp_C= tfp_ss; tfp_D= tfp_ss; Gamma_B= Gamma_ss;
        @#else
            @#if CPshockconfig == "SC4"
                tfp_A= tfp_ss; tfp_B= tfp_ss; tfp_C= tfp_ss; tfp_D= tfp_ss; Gamma_C= Gamma_ss;
            @#endif
        @#endif
    @#endif
@#endif
