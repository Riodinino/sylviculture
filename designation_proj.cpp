void Designate() {

	int site, col, row, sp, sph=0, designated;
	float volume, dbh_min[numespharvestable], min_dbh_min, max_dbh_max=0.0;

	/* getting species vector of minimum harvestable diameter */
	for(sp=1;sp<=numesp;sp++)
		if(S[sp].s_harvestable){
			dbh_min[sph]=S[sp].s_dbhmin; 
			sph++;
			if(S[sp].s_dbhmax > max_dbh_max)
				max_dbh_max = S[sp].s_dbhmax;
		}

	/* getting minimum value of minimum harvestable diameter among species*/
	min_dbh_min = dbh_min[0];
	for(sph=1;sph<numespharvestable;sph++)
		if(dbh_min[sph] < min_dbh_min)
			min_dbh_min = dbh_min[sph];

	/* designating tree, increasing minimum harvestable dbh if needed to br under the objective */
	for(min_dbh_min; min_dbh_min < max_dbh_max; min_dbh_min += 0.1){
		volume=0.0;
		designated=0;
		for(site=0;site<sites;site++){
        	if(T[site].t_age > 0										/*alive tree*/
        		&& S[T[site].t_sp_lab].s_harvestable 					/*harvestable species*/
        		&& T[site].t_dbh >= S[T[site].t_sp_lab].s_dbhmin		/*reached minimum dbh*/
        		&& T[site].t_dbh <= S[T[site].t_sp_lab].s_dbhmax){		/*under maximum dbh*/
        		Tlogging[0][site] = 1;
        		volume += -0.0358 + 8.7634*T[site].t_dbh*T[site].t_dbh; /*volume by ONF-2011 in French Guiana - Center (Kourou)*/
        		designated++;
        	}
        }
        if(volume < designated_volume)
        	break;														/*if the volume is under the objective we can stop */
        else
        	for(sp=1;sp<=numesp;sp++)
        		if(S[sp].s_harvestable)
        			S[sp].s_dbhmin += 0.01;								/*if the volume is greater than the objective we need to derease minimum harvestable diameter for all species */
	}

	cout << designated << " trees have been designated, representing " << volume << " m3." << endl;
    cout << "dbh min is now " << min_dbh_min << endl;
} 
/////////////////////////////////////////////////////////////
		/* Get species vector of interest */

	/* getting species vector of minimum harvestable diameter */
	int interest[numespharvestable], max_interest = 0, min_interest = 10;
	for(sp=1;sp<=numesp;sp++)
		if(S[sp].s_harvestable){
			interest[sph]=S[sp].s_interest; 
			sph++;
			if(S[sp].s_interest > max_interest)
				max_interest = S[sp].s_interest;
			if(S[sp].s_interest < min_interest)
				min_interest = S[sp].s_interest;
		}


	/* Calculate total number of stem and volume */
	int des_high, des_none, des_med; // Number of designated tree by order of interest
	volume=0.0; // Volume of designated trees
	designated=0; // Total number of designated tree
	for(site=0;site<sites;site++){
    	if(T[site].t_age > 0										/*alive tree*/
    		&& S[T[site].t_sp_lab].s_harvestable 					/*harvestable species*/
    		&& T[site].t_dbh >= S[T[site].t_sp_lab].s_dbhmin		/*reached minimum dbh*/
    		&& T[site].t_dbh <= S[T[site].t_sp_lab].s_dbhmax){		/*under maximum dbh*/
    		Tlogging[0][site] = 1; /* TAG : logging */
    		volume += -0.0358 + 8.7634*T[site].t_dbh*T[site].t_dbh; /*volume by ONF-2011 in French Guiana - Center (Kourou)*/
    		designated++;
    		switch(S[T[site].t_sp_lab].s_interest){ // new :categorizes the interest
    			case min_interest:
    			des_low++;
    			break;
    			case max_interest
    			:
    			des_high++;
    			break;
    			default:
    			des_med++;
    		}
    	}
    }

		/* Apply a simple conservation rule to high-interest species */

		/* Discard every species of absolute non-interest */
		for(site=0;site<sites;site++){
        	if(T[site].t_age > 0										/*alive tree*/
        		&& S[T[site].t_sp_lab].s_harvestable 					/*harvestable species*/
        		&& T[site].t_dbh >= S[T[site].t_sp_lab].s_dbhmin		/*reached minimum dbh*/
        		&& T[site].t_dbh <= S[T[site].t_sp_lab].s_dbhmax		/*under maximum dbh*/
        		&& S[T[site].t_sp_lab].s_interest == 0){				/* of no interest */	
        		Tlogging[0][site] = 1;
        		volume -= -0.0358 + 8.7634*T[site].t_dbh*T[site].t_dbh; /*volume by ONF-2011 in French Guiana - Center (Kourou)*/
        		designated--;
        	}

	unsigned short *T_H_interest[2]; // put it at the right place
	unsigned short *T_L_interest[2]; // put it at the right place
	unsigned short *T_M_interest[2];

	int it = 0;
    for(site = 0;site<sites;site++){

    	if(Tlogging[0][site] == 1){
    		Tinterest[1][it] = site;
    		Tinterest[2][it] = S[T[site].t_sp_lab].s_interest;
    		it++; // variable : designated tree number
    	}

    }

	for(des = 0; des < it+1; des++){
		if(Tinterest[2][des] == min_interest){
		Tlogging[0][Tinterest[1][des]] = 0;
		volume -= -0.0358 + 8.7634*T[site].t_dbh*T[site].t_dbh; /*volume by ONF-2011 in French Guiana - Center (Kourou)*/
		}
		//if(designated_volume - 2*cols*NV*rows*NH/10000 < volume < designated_volume + 2*cols*NV*rows*NH/10000){
		//	break
		//}
	}

	 
	if(designated_volume - 2*cols*NV*rows*NH/10000 < volume < designated_volume + 2*cols*NV*rows*NH/10000){

	}

	int(ran_ind)
	/* If the volume is still too high, discard randomly in low interest species */
	if(volume > designated_volume + 2*cols*NV*rows*NH/10000){
		while(volume is not in the interval){
			ran_in = genrand()%(it-1);
		}
	}
	/* Or if the volume is now too low, add randomly in growing interest species */
	else if(volume < designated_volume - 2*cols*NV*rows*NH/10000){
	}
	
genrand()%(it-1);