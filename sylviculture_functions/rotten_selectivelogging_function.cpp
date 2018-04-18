void Rot() {

	int site, rotten=0;
	float protten, volume=0.0;

	/* Calculating selected volume */
	for(site=0;site<sites;site++) 
		if(Tlogging[0][site] == 1)
			volume += -0.0358 + 8.7634*T[site].t_dbh*T[site].t_dbh; /*volume by ONF-2011 in French Guiana - Center (Kourou)*/

	/* evaluates each tree probability to be rotten, and remove it if randomly in the risk to be rotten*/
    for(site=0;site<sites;site++){
       	if(Tlogging[0][site]==1){
       		protten = 1 / (1 + exp(-(-5.151 + 0.042*T[site].t_dbh*100))); /*Probability to be rotten*/
       		if(genrand2() < protten){
       			Tlogging[0][site]=0;
       			rotten++;
           		volume -= -0.0358 + 8.7634*T[site].t_dbh*T[site].t_dbh; /*volume by ONF-2011 in French Guiana - Center (Kourou)*/
       		}
       	}
       } 
    cout << rotten << " trees are rotten, volume is now of " << volume << " m3." << endl;
}
