
@#if CPshockconfig == "SC1"
    tfp_A = tfp_A(-1)^rho_tfp*exp(-eps_tfp);
    tfp_B = tfp_B(-1)^rho_tfp*exp(-eps_tfp);
    tfp_C = tfp_C(-1)^rho_tfp*exp(-eps_tfp);
    tfp_D = tfp_D(-1)^rho_tfp*exp(-eps_tfp);
@#else
    @#if CPshockconfig == "SC2"
        Gamma_A = (1-rho_Gamma)*rhob_Gamma+Gamma_A(-1)*rho_Gamma-eps_Gamma;
        tfp_A = tfp_A(-1)^rho_tfp*exp(-eps_tfp);
        tfp_B = tfp_B(-1)^rho_tfp*exp(-eps_tfp);
        tfp_C = tfp_C(-1)^rho_tfp*exp(-eps_tfp);
        tfp_D = tfp_D(-1)^rho_tfp*exp(-eps_tfp);
    @#else
        @#if CPshockconfig == "SC3"
            Gamma_B = (1-rho_Gamma)*rhob_Gamma+Gamma_A(-1)*rho_Gamma-eps_Gamma;
            tfp_A = tfp_A(-1)^rho_tfp*exp(-eps_tfp);
            tfp_B = tfp_B(-1)^rho_tfp*exp(-eps_tfp);
            tfp_C = tfp_C(-1)^rho_tfp*exp(-eps_tfp);
            tfp_D = tfp_D(-1)^rho_tfp*exp(-eps_tfp);
        @#else
            @#if CPshockconfig == "SC4"
                Gamma_C = (1-rho_Gamma)*rhob_Gamma+Gamma_A(-1)*rho_Gamma-eps_Gamma;
                tfp_A = tfp_A(-1)^rho_tfp*exp(-eps_tfp);
                tfp_B = tfp_B(-1)^rho_tfp*exp(-eps_tfp);
                tfp_C = tfp_C(-1)^rho_tfp*exp(-eps_tfp);
                tfp_D = tfp_D(-1)^rho_tfp*exp(-eps_tfp);
            @#endif
        @#endif
    @#endif
@#endif
