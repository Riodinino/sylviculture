void GapDamages() {

	cout << "###   Selective Logging Long Term Damages   ###" << endl;
	int site, row, col, siteG, rowG, colG;
    float deathrate, gaps_deathrate, gaps_hurt, d, dgaps[sites];

    /*Initialise dgaps to maximum distance to a gaps (null effect)*/
    for(site=0;site<sites;site++) dgaps[site] = rows*rows + cols*cols;

    /*Compute for each tree the distance to the closest gap*/
    for(siteG=0;siteG<sites;siteG++){
        if(Tlogging[2][siteG] == 1){
        	if(T[siteG].t_age > 1){ // New trees could have been recruited over a year
        		rowG = floor(siteG/cols);
				colG = siteG-(rowG*cols);
				for(site=0;site<sites;site++){
					row = floor(site/cols);
        			col = site-(row*cols);
         			d = (row - rowG)*(row - rowG) + (col - colG)*(col - colG);
         			if(d < dgaps[site])	dgaps[site] = d;
         		}
			}
		}
    }

    /*Hurt trees depending on their distance to a gaps following an allometry fitted with Paracou data*/
    for(site=0;site<sites;site++){
       	if(T[site].t_age != 0 && T[site].t_dbh > 0.1){ //tree with dbh<10 have not an increased mortality closed to gaps, on the contrary they'll have a tendency 
       		gaps_deathrate = -4.441 + 0.762*exp(0.064*sqrt(dgaps[site]));
       		gaps_deathrate = exp(gaps_deathrate) / (1 + exp(gaps_deathrate)); // Allometry representing gaps damages
       		deathrate = T[site].t_s->DeathRate(T[site].t_PPFD, T[site].t_dbh, T[site].t_NPPneg);
       		if(gaps_deathrate > deathrate){
       			gaps_hurt = T[site].t_Tree_Height/(2*(gaps_deathrate - deathrate));
       			T[site].t_hurt += gaps_hurt;
       		}
       	}
    }        
}