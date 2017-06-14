
@#if liqinj == "on"
    m=rho_m*m(-1)+pm*(rb_av-rb_ss)+eps_m;
    
    m_AB = b_BA - l_AB;
    m_CA = b_AC - l_CA;
    m_DA = b_AD - l_DA;
           
    m_AB = m/3;
    m_CA = m/3;
    m_DA = m/3;
@#else
    @#if liqinj == "off"
       b_BA = l_AB; 
       b_AC = l_CA;
       b_AD = l_DA;
    @#endif
@#endif
