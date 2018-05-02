void Select() {
// careful with the volume area conversion
// careful with declaring other dynamic tab           
	int site, sp, i, rank, rankmax=0, unselected=0;
	float volume=0.0;

	/* Calculating designated volume */
	for(site=0;site<sites;site++) 
		if(Tlogging[0][site] == 1)
			volume += -0.0358 + 8.7634*T[site].t_dbh*T[site].t_dbh; /*volume by ONF-2011 in French Guiana - Center (Kourou)*/
    
	if(volume <= harvested_volume)
        cout << "All designated trees will be harvested." << endl;

	if(volume > harvested_volume){

		/* determining maximal interest rank (leat valuable species) */
		for(sp=1;sp<=numesp;sp++)
        	if(S[sp].s_interest > rankmax)
        		rankmax = S[sp].s_interest;	

        /* determining determining headcount for each rank */
        int rank_nb[rankmax];
        for(rank=0;rank<rankmax;rank++) rank_nb[rank]=0;
        for(site=0;site<sites;site++)
        	if(S[sp].s_harvestable)
        		rank_nb[S[T[site].t_sp_lab].s_interest]++;

        /* removing tree untill wanted volume is reached starting by highest rank */
        for(rank=rankmax-1;rank>=0;rank--){
        	while(rank_nb[rank]>0){
        		site=floor(genrand2()*sites);
        		if(Tlogging[0][site]==1 && S[T[site].t_sp_lab].s_interest==rank){
        			Tlogging[0][site]=0;
        			rank_nb[rank]--;
        			unselected++;
        			volume -= -0.0358 + 8.7634*T[site].t_dbh*T[site].t_dbh; /*volume by ONF-2011 in French Guiana - Center (Kourou)*/
        			if(volume <= harvested_volume) break;
        		}
        	}
        	if(volume <= harvested_volume) break;
        }

        cout << unselected << " trees have been unselected, volume is now of " << volume << " m3." << endl;
	}
}


