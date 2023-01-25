struct Pa
{
	double M_A, m_A, M_B, m_B, n_A, n_B, c_n_A, b_n_A, kap_A, c_n_B, b_n_B, kap_B, R_0_A, R_0_B, f_R_A, f_R_B, r, p_A, p_B;

    struct Dark
    {
    	int on;
    	double alpha, mV;
    } dark;

  	struct PN
  	{
  		int PNc, PN1, PN2, PN25, PN3;
  	} pn;	

  	char print_evolve;
	int initial_only;
};

struct Eqm {
	double a1_A, a2_A, a3_A, a1_B, a2_B, a3_B, orb, r_phy;
};
