
@#if liqinj == "on"
    m=rho_m*m(-1)+pm*(rb_av-rb_ss)+eps_m;
    
    m_BA = b_AB - l_BA;
    m_AC = b_CA - l_AC;
    m_AD = b_DA - l_AD;
           
    m_BA = m/3;
    m_AC = m/3;
    m_AD = m/3;
@#else
    @#if liqinj == "off"
       b_AB = l_BA; 
       b_CA = l_AC;
       b_DA = l_AD;
    @#endif
@#endif
