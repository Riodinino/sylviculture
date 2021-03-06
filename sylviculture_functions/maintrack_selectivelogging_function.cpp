void MainTracks() {

    int site, row, col, individuals=0;
    float volume=0.0;
    
    // Find what is the row_num of the road edge
    int row_num;
    float segment_volume;
    float segment_area;

    for(int row_n=0;row_n<rows;row_n++){
        // Calculate segment dimension
        segment_area = (row_n+1)*ncols; // in m²
        segment_area /= 10000; 
        // Extrapolate volume
        segment_volume = harvested_volume*segment_area;
        //Compare it with harvestable volume

        //if(segment_volume >97 && segment_volume < 103){ could be written like that but less permissive
        if(segment_volume < 101){
        	row_num = row_n;
        	break;
        }
        
    }
   // harvested volume or designated volume to be used ?
        
   
    
    for(row=0;row<row_num;row++){ // Track length
        for(col=((cols/2)-3);col<((cols/2)+3);col++){ // Track width and position (central)
        	site = col+row*cols;
        	Tlogging[1][site] = 1;
        	if(T[site].t_age != 0) {
        		volume += -0.0358 + 8.7634*T[site].t_dbh*T[site].t_dbh; /*volume by ONF-2011 in French Guiana - Center (Kourou)*/
        		output[36] << "MT" << "\t" << col << "\t" << row << "\t" << T[site].t_age << "\t" << T[site].t_dbh << "\t" << T[site].t_Tree_Height << "\t" << T[site].t_Crown_Radius << "\t" << T[site].t_Crown_Depth << "\t" << T[site].t_sp_lab << endl;
            	T[site].Death();
            	individuals++;
        	}
       	}
    }    
    cout << individuals << " trees have been killed for the main track representing " << volume << " m3." << endl;
}