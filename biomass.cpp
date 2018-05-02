#ifdef INTRASPECIFIC
        t_s->s_output_field[7] += 0.0673*pow(t_wsg*t_Tree_Height*LV*t_dbh*t_dbh*LH*LH*10000, 0.976);  // this is the allometric equ 4 in Chave et al. 2014 Global Change Biology to compute above ground biomass
#else
        t_s->s_output_field[7] += 0.0673*pow(t_s->s_wsg*t_Tree_Height*LV*t_dbh*t_dbh*LH*LH*10000, 0.976);  // this is the allometric equ 4 in Chave et al. 2014 Global Change Biology to compute above ground biomass
#endif


 