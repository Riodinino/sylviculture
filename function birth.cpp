/*##############################################
 ####	            Tree birth              ####
 ####  called by BirthInit and UpdateTree   ####
 ##############################################*/


void Tree::Birth(Species *S, int nume, int site0, float dbh_measured) {
    t_site = site0;
    t_sp_lab = nume;            /* t_sp_lab is the species label of a site. Can be defined even if the site is empty (cf. persistence function defined in Chave, Am Nat. 2001) */
    t_NPPneg=0.0;
    t_s = S+t_sp_lab;
    t_age = 1;
    t_from_Data = 1;
    t_hurt = 0;
    //t_Tree_Height = H0; // change
    t_hmax = (t_s->s_hmax);
    t_ah = (t_s->s_ah);

    if((t_s->s_dmax)*1.5 > dbh_measured) t_dbh = dbh_measured;  // UNIT ? NINO   // force dbh to be within limits of TROLL specifications
    else {
        t_dbh = (t_s->s_dmax); // Yes but, is intraspecific variability to be taken into account in this part ?  Cause it impacts the individual dbh threshold.
        cout << "Warning: DBH_measured > 1.5*DBH_max for species. DBH set to DBH_max for species \n";
    }

#ifdef INTRASPECIFIC
    int dev_rand = genrand2i()%100000;
    t_intraspecific_multiplier_height = d_intraspecific_height[dev_rand];
    t_intraspecific_multiplier_CR = d_intraspecific_CR[dev_rand];
    t_intraspecific_multiplier_N = d_intraspecific_N[dev_rand];
    t_intraspecific_multiplier_P = d_intraspecific_P[dev_rand];
    t_intraspecific_multiplier_LMA = d_intraspecific_LMA[dev_rand];
    t_intraspecific_multiplier_CD = d_intraspecific_CD[dev_rand];
    t_intraspecific_deviation_wsg = d_intraspecific_wsg[dev_rand];
    t_intraspecific_multiplier_dmax = d_intraspecific_dmax[dev_rand];

    t_Pmass = (t_s->s_Pmass) * t_intraspecific_multiplier_P;
    t_Nmass = (t_s->s_Nmass) * t_intraspecific_multiplier_N;
    t_LMA = (t_s->s_LMA) * t_intraspecific_multiplier_LMA;
    t_wsg = maxf((t_s->s_wsg) + t_intraspecific_deviation_wsg,0.1); // cutoff for normal distribution (could also be set at 0.0)
    t_dmax = (t_s->s_dmax) * t_intraspecific_multiplier_dmax;

#ifdef ALLOM_relative
    t_ah *= t_dmax; /* if the Michaelis Menten parameters are inferred from a relative diameter (t_dbh_rel = t_dbh/t_dmax) instead of the unscaled parameter t_dbh, then the equation h_rel = t_hmax * t_dbh_rel/ (t_dbh_rel + t_ah) is equivalent to t_hmax * t_dbh / (t_dbh + t_ah_rescaled), where ah_rescaled = t_ah * t_dmax. This means that, apart from rescaling of t_ah, everything else can be calculated as usually */
#endif
   // t_dbh=t_ah*H0/(t_hmax*t_intraspecific_multiplier_height-H0);    // change  we make sure that all trees start at the same height, this means adapting the diameter depending on inter- and intraspecific changes

    float SLA=10000.0/t_LMA;
    
    t_Vcmaxm=pow(10.0, minf((-1.56+0.43*log10(t_Nmass*1000.0)+0.37*log10(SLA)), (-0.80+0.45*log10(t_Pmass*1000.0)+0.25*log10(SLA)))); // this is equation 2 in Domingues et al 2010 PCE (coefficients from fig7) which made better fits than equation 1 (without LMA)
    t_Jmaxm=pow(10.0, minf((-1.50+0.41*log10(t_Nmass*1000.0)+0.45*log10(SLA)), (-0.74+0.44*log10(t_Pmass*1000.0)+0.32*log10(SLA)))); // added as a Species member variable 14-04-2015; this is equ 2 in Domingues et al 2010 PCE (coefficients from fig7)
    t_Vcmax=t_Vcmaxm*t_LMA;
    t_Jmax=t_Jmaxm*t_LMA;
    t_Rdark=t_LMA*(8.5341-130.6*t_Nmass-567.0*t_Pmass-0.0137*t_LMA+11.1*t_Vcmaxm+187600.0*t_Nmass*t_Pmass)*0.001; //t_Rdark corresponds to leaf maintenance respiration. From Table 6 in Atkin et al 2015 New phytologist v.2.0
    
    t_Gamma = 38.0*iCair; // s_Gamma at 25°C computed according to von Caemmerer 2000 formula: gamma=Kc*O*0.25/(2*Ko), with Kc=260 microbar, Ko=179mbar and O=210 mbar (the last value is from Farquhar et al 1980, the first two one are from von Caemmerer 2000 table 2.3 page 45). gamma is set to 36.9 on Atkin et al 2015. Could be a global variable. v.2.0
    
    t_leaflifespan = pow(10,(2.040816*(2.579713-log10(SLA)))); //this is the expression from Reich et al. 1997 PNAS (provides probably more realistic estimates for species with high LMA).
    //t_leaflifespan=1.5+pow(10,(7.18+3.03*log10(t_LMA*0.0001)));           //this is the expression from Reich et al 1991 Oecologia (San Carlos Rio Negro).
    //t_leaflifespan=0.5+pow(10,(-2.509+1.71*log10(t_LMA)));    //this is the expression from Wright et al 2004 Nature (leaf economics spectrum).

#ifdef LL_limit
    t_leaflifespan = maxf(t_leaflifespan,3.0);
#endif
    
    t_time_young=minf(t_leaflifespan/3.0,1.0);
    t_time_mature=t_leaflifespan/3.0;
    t_time_old=t_leaflifespan-t_time_mature-t_time_young;
    t_litter = 0.0;
    t_fci = 0.0;

//    t_Crown_Radius = exp(1.9472 + 0.5925*log(t_dbh));   /* this is crown allometry derived from data set provided by Jucker et al. 2016 (Global Change Biology) */

//    t_Crown_Radius  *= t_intraspecific_multiplier_CR;
//   t_Crown_Depth = de0 * t_intraspecific_multiplier_CD;
    
    t_dbh_thresh = t_dmax;  // contrary to version without intraspecific variation, dmax already includes variation, so t_dbh_thresh = t_dmax
 //   t_dbhmature=t_dbh_thresh*0.5; // this correponds to the mean thresholds of tree size to maturity, according to Visser et al. 2016 Functional Ecology (suited to both understory short-statured species, and top canopy large-statured species). NOTE that if we decide to keep it as a fixed species-specific value, this could be defined as a Species calss variable, and computed once in Species::Init. -- v230
    //float u=genrand2();
    //t_dbhmature=maxf(0, -(t_dmax)*0.25*log((1-u)/u)+t_dmax*0.5);    // IM test 02-2017, try to introduce intra-species inter-individual variability in dbhmature, following a sigmoidal repartition function, as in Visser et al. 2016 and Wright et al. 2005
    float hrealmax = t_intraspecific_multiplier_height * 1.5*t_hmax*t_dbh_thresh/(1.5*t_dbh_thresh+t_ah);
    
#else
    //t_dbh=(t_s->s_dbh0);
#ifdef ALLOM_relative
    t_ah *= t_s->s_dmax; /* if the Michaelis Menten parameters are inferred from a relative diameter (t_dbh_rel = t_dbh/t_dmax) instead of the unscaled parameter t_dbh, then the equation h_rel = t_hmax * t_dbh_rel/ (t_dbh_rel + t_ah) is equivalent to t_hmax * t_dbh / (t_dbh + t_ah_rescaled), where ah_rescaled = t_ah * t_dmax. This means that, apart from rescaling of t_ah, everything else can be calculated as usually */
#endif
    //t_dbh=t_ah*H0/(t_hmax-H0);    // we make sure that all trees start at the same height, this means adapting the diameter depending on inter- and intraspecific changes

    t_dbh_thresh = ((t_s->s_dmax)-t_dbh)*flor(1.0+log(genrand2())*0.01)+t_dbh;
 //   t_dbhmature=t_dbh_thresh*0.5; // this correponds to the mean thresholds of tree size to maturity, according to Visser et al. 2016 Functional Ecology (suited to both understory short-statured species, and top canopy large-statured species). NOTE that if we decide to keep it as a fixed species-specific value, this could be defined as a Species calss variable, and computed once in Species::Init. -- v230
    //float u=genrand2();
    //t_dbhmature=maxf(0, -(t_s->s_dmax)*0.25*log((1-u)/u)+t_s->s_dmax*0.5);    // IM test 02-2017, try to introduce intra-species inter-individual variability in dbhmature, following a sigmoidal repartition function, as in Visser et al. 2016 and Wright et al. 2005
    float hrealmax=1.5*t_hmax*t_dbh_thresh/(1.5*t_dbh_thresh+t_ah);
#endif
    
    if(_BASICTREEFALL || _TREEFALL) t_Ct = hrealmax*flor(1.0-vC*sqrt(-log(genrand2())));
    
    
    t_ddbh=0.0;
    t_dbhmature=t_dbh_thresh*0.5; // this correponds to the mean thresholds of tree size to maturity, according to Visser et al. 2016 Functional Ecology (suited to both understory short-statured species, and top canopy large-statured species). NOTE that if we decide to keep it as a fixed species-specific value, this could be defined as a Species calss variable, and computed once in Species::Init. -- v230


#ifdef(INTRASPECIFIC) // NINO
    t_Tree_Height = minf(t_intraspecific_multiplier_height * t_hmax * t_dbh/(t_dbh + t_ah),HEIGHT-1); // todef NINO-OK
#else
    t_Tree_Height = t_hmax * t_dbh/(t_dbh + t_ah);
#endif

t_Crown_Radius = exp(1.9472 + 0.5925*log(t_dbh)); /* this is crown allometry derived from data set provided by Jucker et al. 2016 (Global Change Biology) */
if(t_Tree_Height < 5.0) {t_Crown_Depth = de0 + 0.17 * (t_Tree_Height-H0);}
else {t_Crown_Depth = de0+0.26*(t_Tree_Height-H0) - 0.09 * (5.0-H0);} /* allometry deduced from Piste Saint-Elie dataset (unpublished) */
    
#ifdef INTRASPECIFIC
    t_Crown_Radius *= t_intraspecific_multiplier_CR;
    t_Crown_Depth *= t_intraspecific_multiplier_CD;
#endif

    t_dens=dens;
    t_leafarea=t_dens*PI*t_Crown_Radius*LH*t_Crown_Radius*LH*t_Crown_Depth;
    t_youngLA=0.25*t_leafarea;
    t_matureLA=0.5*t_leafarea;
    t_oldLA=0.25*t_leafarea;
    Fluxh(int(t_Tree_Height)+1);
    tempRday=0.0;

// useless in my case, so unchecked
#ifdef CROWN_SHAPE
    t_Crown_Volume_layer = 0.0;
    
    int crown_base = int(t_Tree_Height - t_Crown_Depth);
    int crown_top = int(t_Tree_Height);
    int crown_center = int(t_Tree_Height-0.5*t_Crown_Depth);
    float crown_extent_top = t_Tree_Height - crown_center;
    float crown_extent_base = crown_center - (t_Tree_Height - t_Crown_Depth);
    
    if(t_Crown_Depth > 2.0){        // only consider cases, where more than two layers are reached by the crown
        /* slopes are calculated as horizontal/vertical slopes (inverse of more intuitive definition vertical/horizontal), reasoning: when shape_crown = 1.0, this slope is still defined */
        /* the slope at the base part of the crown is steeper, meaning that the trapezoid narrows down faster (limit: triangle/cone) */
        t_Crown_Slope_Top = t_Crown_Radius * (1.0 - shape_crown) / crown_extent_top;
        t_Crown_Slope_Bottom = min(2 * t_Crown_Slope_Top, t_Crown_Radius/crown_extent_base);
    } else {
        t_Crown_Slope_Top = 0.0;
        t_Crown_Slope_Bottom = 0.0;
    }
    t_Crown_Volume = 0.0;
    for(int h=crown_base;h<crown_top+1;h++){
        /* calculating the height of the current layer */
        float height_layer = 1.0;
        if(crown_top == crown_base) height_layer = t_Crown_Depth;
        else if(h == crown_top) height_layer = (t_Tree_Height-crown_top);
        else if(h == crown_base) height_layer = (crown_base+1-(t_Tree_Height-t_Crown_Depth));
        
        /* calculating the radius of the current layer depending on the respective slopes */
        float radius_layer = t_Crown_Radius;
        if(crown_top == crown_base){}
        else if(h == crown_base) radius_layer -= t_Crown_Slope_Bottom * crown_extent_base;
        else if(h < crown_center ) radius_layer -= t_Crown_Slope_Bottom * (crown_center-h);
        else radius_layer -= t_Crown_Slope_Top * (h - crown_center);            // for h = crown_center full radius
        
        float crown_area_fl = maxf(0.0,PI * radius_layer * radius_layer);
        t_Crown_Volume += crown_area_fl * height_layer;
    }



    (t_s->s_nbind)++;
    nblivetrees++;
}