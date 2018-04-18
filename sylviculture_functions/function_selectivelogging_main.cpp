/*#######################################
 ####    Simulate selective logging   ###
 #######################################*/

void SelectiveLogging() {

    if(iter == disturb_iter) {

    	int site;
    	for(site=0;site<sites;site++){
    		Tlogging[0][site]=0;		// tree felling
    		Tlogging[1][site]=0;		// tracks
    		Tlogging[2][site] = 0;			// gaps
    	}

    	cout << "###   Selective Logging   ###" << endl;
    	Designate();
    	Select();
    	Rot();
    	Fell();
    	MainTracks();
    	SecondaryTracks();
        cout << "### Selective Logging done ###" << endl;
    }

    if(iter == (disturb_iter+iterperyear)) {
    	GapDamages();  
    	int i;
    	for (i=0; i<3; i++) delete [] Tlogging[i];	 // free memory
    }
}