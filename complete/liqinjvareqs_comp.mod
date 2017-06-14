
@#if liqinj == "on"
    m=rho_m*m(-1)+pm*(rb_av-rb_ss)+eps_m;
    
    m_AC = b_CA - l_AC;
    m_BA = b_AB - l_BA;
    m_CD = b_DC - l_CD;
    m_DB = b_BD - l_DB;
    m_AB = b_BA - l_AB;
    m_AD = b_DA - l_AD;
    m_BC = b_CB - l_BC;
    m_BD = b_DB - l_BD;
    m_CA = b_AC - l_CA;
    m_CB = b_BC - l_CB;
    m_DA = b_AD - l_DA;
    m_DC = b_CD - l_DC; 
           
    m_AC = m/12;
    m_BA = m/12;
    m_CD = m/12;
    m_DB = m/12;
    m_AB = m/12;
    m_AD = m/12;
    m_BC = m/12;
    m_BD = m/12;
    m_CA = m/12;
    m_CB = m/12;
    m_DA = m/12;
    m_DC = m/12; 
@#else
    @#if liqinj == "off"
        b_CA = l_AC;
        b_AB = l_BA;
        b_DC = l_CD;
        b_BD = l_DB;
        b_BA = l_AB;
        b_DA = l_AD;
        b_CB = l_BC;
        b_DB = l_BD;
        b_AC = l_CA;
        b_BC = l_CB;
        b_AD = l_DA;
        b_CD = l_DC; 
    @#endif
@#endif
