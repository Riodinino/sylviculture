void SecondaryTracks() {

	    int load[sites], tracks[sites], individuals=0, felt=0;
        int site, col, row, site0, row0, col0, siteT, rowT, colT;
        float d, d0, volume=0.0;

        /* Counting number of felt trees to skid */
        for(site=0;site<sites;site++) felt += Tlogging[0][site];
        
        while(felt > 0){

        	/*Computing loadings and tracks distance for each tree*/
        	for(site0=0;site0<sites;site0++){ 
        		load[site0]=0;
        		tracks[site0]=rows*rows + cols*cols;
        		row0 = floor(site0/cols);
        		col0 = site0-(row0*cols);
        		for(site=0;site<sites;site++){
        			if(Tlogging[0][site]==1 || Tlogging[1][site]==1){ // compute distance if the site is a felt tree or a track
        				row = floor(site/cols);
        				col = site-(row*cols);
        				d = (row - row0)*(row - row0) + (col - col0)*(col - col0);
        				if(Tlogging[0][site]==1 && d <= (30*30)) // site can evacuate the tree if is at a distance smaller than 30 meters
        					load[site0]++;
        				if(Tlogging[1][site]==1 && d < tracks[site0]) // save the track distance if it's closest than the previously saved one
        					tracks[site0]=d;
        			}
        		}
        	}

        	/*Seeking the best place to start the secondary track*/
        	site0=0;
        	for(site=0;site<sites;site++){
        		if(load[site]>load[site0]) // best candidate is the one which can evacuate maximum number of trees
        			site0=site;
        		if(load[site]==load[site0] && tracks[site]<tracks[site0]) // for equal loadings, best candidate is the one with a minimum distance to join an existing track
        			site0=site;
        	}	

        	/*Seeking for the closest track*/
        	row0 = floor(site0/cols);
        	col0 = site0-(row0*cols);
        	d0 = rows*rows+cols*cols;
        	for(site=0;site<sites;site++){
        		if(Tlogging[1][site]==1){ // if it's a track compute distance to the track
        			rowT = floor(site/cols);
        			colT = site-(rowT*cols);
        			d = (row0 - rowT)*(row0 - rowT) + (col0 - colT)*(col0 - colT);
        			if(d<d0){ // if the track is closer than the previously saved one, keep the location
        				siteT=site;
        				d0=d;
        			}
        		}
        	}

        	/*Trace the secondary track*/
        	rowT = floor(siteT/cols);
        	colT = siteT-(rowT*cols);
        	do{
        		do {
            		for(int i=-2;i<=2;i++){ 
        				for(int j=-2;j<=2;j++){
        					site = (col0+i)+(row0+j)*cols;
        					if(site>=0 && site<sites) Tlogging[1][site]=1; //flag the track with a size of 4 meters
        				}
        			}
        			for(site=0;site<sites;site++){ 
        				if(Tlogging[0][site]==1){
        					row = floor(site/cols);
        					col = site-(row*cols);
        					d = (row - row0)*(row - row0) + (col - col0)*(col - col0);
        					if(d <= (33*33)){ //unflag served trees in a radius of 30 meters
        						Tlogging[0][site]=0;
        						felt--;
        					}
        				}
        			}
        			if(col0 > colT) col0--; //move in direction of the closest existing track
        			if(col0 < colT) col0++;
        			if(row0 > rowT) row0--;
        			if(row0 < rowT) row0++;
        		} while(row0 != rowT); //stop when we reach the closest existing track
        	} while(col0 != colT);
        	cout << "A secondary track have been traced, " << felt << " trees still need to be evacuated." << endl; //!LONG! computation, console output to follow advancement
		}

		/* Removing trees on secondary tracks */
        for(site=0;site<sites;site++){ 
        	if(Tlogging[1][site] == 1 && T[site].t_age != 0){
        		row = (site/cols);
        		col = site-(row*cols);
        		volume += -0.0358 + 8.7634*T[site].t_dbh*T[site].t_dbh; /*volume by ONF-2011 in French Guiana - Center (Kourou)*/
        		output[36] << "ST" << "\t" << col << "\t" << row << "\t" << T[site].t_age << "\t" << T[site].t_dbh << "\t" << T[site].t_Tree_Height << "\t" << T[site].t_Crown_Radius << "\t" << T[site].t_Crown_Depth << "\t" << T[site].t_sp_lab << endl;
            	T[site].Death();
            	individuals ++;
        	}
        }
    cout << individuals << " trees have been killed for secondary tracks representing " << volume << " m3." << endl;
}