void Fell() {

	int site, row, col, felled=0;
	float volume=0.0;

	/* fell the selected tree not rotten */
    for(site=0;site<sites;site++){
        if(Tlogging[0][site]==1){
        	row = floor(site/cols);
        	col = site-(row*cols);
        	volume += -0.0358 + 8.7634*T[site].t_dbh*T[site].t_dbh; /*volume by ONF-2011 in French Guiana - Center (Kourou)*/
        	output[36] << "L" << "\t" << col << "\t" << row << "\t" << T[site].t_age << "\t" << T[site].t_dbh << "\t" << T[site].t_Tree_Height << "\t" << T[site].t_Crown_Radius << "\t" << T[site].t_Crown_Depth << "\t" << T[site].t_sp_lab << endl;
           	T[site].FellTree();
           	felled ++;
        }
    } 
    cout << felled << " trees have been felled representing " << volume << " m3." << endl;
}




