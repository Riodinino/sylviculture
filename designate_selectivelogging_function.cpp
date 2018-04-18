
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