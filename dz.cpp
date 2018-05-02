    
/*##################################
 ####           Treefall         ####
 #### called by UpdateTreefall   ####
 ####################################*/
    

/* NEW in TROLL v.2.4: FallTree() function has become Treefall() function, calculation of angle and treefall outside of function, and damages are now added up from several treefalls */
    
void Tree::Treefall(float angle) {
    
    /* treefall statistics */
    nbTreefall1++;
    if(t_dbh*LH>0.1) nbTreefall10++; 
    if(t_dbh*LH>0.3) nbTreefall30++;
    
    int xx,yy;
    int row0,col0,h_int, r_int;

    float h_true = t_Tree_Height*LV;
    h_int = int(h_true*NH);
    row0=t_site/cols;
    col0=t_site%cols;
    
    /* update of Thurt field at the site of the tree, for consistency */
    Thurt[0][t_site+sites] = int(t_Tree_Height);
    
    /* fallen stem destructs other trees */
    for(int h=1;h<h_int;h++) {                      // loop on the fallen stem (horizontally)
        xx=int(flor(col0+h*cos(angle)));          // get projection in col (= xx) direction, where xx is absolute location
        if(xx<cols){
            yy=   int(row0+h*sin(angle));         // get projection in row (= yy) direction, where yy is absolute location
            Thurt[0][xx+(yy+rows)*cols] += int(t_Tree_Height);
            // Thurt[0] where the stem fell, calculation: xx+(yy+rows)*cols= xx + yy*cols + rows*cols = xx + yy*cols + sites / NEW in v.2.4: addition of damage instead of setting equal in order to account for cumulative damage (several treefalls hitting the same site)
        }
    }
    
    /* fallen crown destructs other trees, less damaging than stem */
    xx=col0+int((h_true*NH-t_Crown_Radius)*cos(angle));
    yy=row0+int((h_true*NH-t_Crown_Radius)*sin(angle));
    r_int = int(t_Crown_Radius);
    for(int col=max(0,xx-r_int);col<min(cols,xx+r_int+1);col++) { // loop on the fallen crown (horizontally)
        for(int row=yy-r_int;row<yy+r_int+1;row++) {
            if((col-xx)*(col-xx)+(row-yy)*(row-yy)<r_int*r_int) Thurt[0][col+(row+rows)*cols] += int((t_Tree_Height-t_Crown_Radius*NV*LH)*0.5); // less severe damage than stem / NEW in v.2.4: addition of damage instead of setting equal in order to account for cumulative damage (several treefalls hitting the same site)
        }
    }
    /* v.2.4.0: outputs have been moved to Death() function */
    Death();
}

// Trees are saved in output[36]

