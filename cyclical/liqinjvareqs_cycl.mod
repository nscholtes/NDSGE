  
@#if liqinj == "on"
    m=rho_m*m(-1)+pm*(rb_av-rb_ss)+eps_m;

    m_AC = b_CA - l_AC;
    m_BA = b_AB - l_BA;
    m_CD = b_DC - l_CD;
    m_DB = b_BD - l_DB;    
           
    m_AC = m/4;
    m_BA = m/4;
    m_CD = m/4;
    m_DB = m/4; 
@#else
    @#if liqinj == "off"
        b_CA = l_AC;
        b_AB = l_BA;
        b_DC = l_CD;
        b_BD = l_DB;
    @#endif
@#endif
