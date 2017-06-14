
@#if CPshockconfig == "SC1"
    tfp_A tfp_B tfp_C tfp_D; 
@#else
    @#if CPshockconfig == "SC2"
        tfp_A tfp_B tfp_C tfp_D Gamma_A;
    @#else
        @#if CPshockconfig == "SC3"
          tfp_A tfp_B tfp_C tfp_D Gamma_B;
        @#else
            @#if CPshockconfig == "SC4"
                tfp_A tfp_B tfp_C tfp_D Gamma_C;
            @#endif
        @#endif
    @#endif
@#endif
