#define __STDC_FORMAT_MACROS

#include "alignment_functions.h"


static int64_t PAM[5][5]={4,-3,-3,-3,-3,-3,4,-3,-3,-3,-3,-3,4,-3,-3,-3,-3,-3,4,-3,-3,-3,-3,-3,-3};
int valOfNucl(char c){
	if(c=='A') return 0;
	if(c=='C') return 1;
	if(c=='G') return 2;
	if(c=='T') return 3;
    if(c=='-') return 4;
	return 4;
}

int64_t compare_letters(char a, char b){
    if(a != 'N' && a != '-') return (a == b) ? POINT : -POINT;
    return -POINT;
}

uint64_t get_max_length_of_sequences(char ** a, uint64_t n_x){
    uint64_t max_len = 0;
    uint64_t len;
    uint64_t i;
    for(i=0;i<n_x;i++){
        len = strlen(a[i]);
        if(len > max_len) max_len = len;
    }
    return max_len;
}

uint64_t get_id_of_longest_sequence(char ** a, uint64_t n_x){
    uint64_t max_len = 0;
    uint64_t max_id = 0;
    uint64_t len;
    uint64_t i;
    for(i=0;i<n_x;i++){
        len = strlen(a[i]);
        if(len > max_len) { max_len = len; max_id = i; }
    }
    return max_id;
}

void pile_chars_up(char ** list_seq_X, char ** list_seq_Y, char * list_piled_X, char * list_piled_Y, uint64_t nx, uint64_t ny, uint64_t posx, uint64_t posy){
    uint64_t i;
    for(i=0;i<nx;i++) list_piled_X[i] = list_seq_X[i][posx];
    for(i=0;i<ny;i++) list_piled_Y[i] = list_seq_Y[i][posy];
}

int64_t compare_piled_up_chars(char * list_piled_X, char * list_piled_Y, uint64_t nx, uint64_t ny){

    if(nx == 1 && ny == 1) return compare_letters(list_piled_X[0], list_piled_Y[0]);

    uint64_t i;
    int64_t acum = 0;
    uint64_t types_X[6] = {0,0,0,0,0,0}; // A, C, G, T, N, -
    uint64_t types_Y[6] = {0,0,0,0,0,0};
    for(i=0;i<nx;i++){
        switch(list_piled_X[i]){
            case 'A': types_X[0]++; break;
            case 'C': types_X[1]++; break;
            case 'G': types_X[2]++; break;
            case 'T': types_X[3]++; break;
            case 'N': types_X[4]++; break;
            case '-': types_X[5]++; break;
        }
    }
    for(i=0;i<ny;i++){
        switch(list_piled_Y[i]){
            case 'A': types_Y[0]++; break;
            case 'C': types_Y[1]++; break;
            case 'G': types_Y[2]++; break;
            case 'T': types_Y[3]++; break;
            case 'N': types_Y[4]++; break;
            case '-': types_Y[5]++; break;
        }
    }
    for(i=0;i<4;i++){
        if(types_X[i] > 0 && types_Y[i] > 0) acum += (types_X[i] + types_Y[i])*POINT;
        types_X[i] -= types_Y[i];
        acum += labs(types_X[i]) * (-POINT);
    }
    for(i=4;i<6;i++){
        acum += (types_X[i] + types_Y[i])*(-POINT);
    }

    return acum;
}

void build_unanchored_alignment(int64_t * cell_path_y, char ** seq_piles_up_X, char ** seq_piles_up_Y, uint64_t n_x, uint64_t n_y){

    uint64_t l_A = get_max_length_of_sequences(seq_piles_up_X, n_x), l_B = get_max_length_of_sequences(seq_piles_up_Y, n_y);
    
    Point p0 = {0,0};
    Point p1 = {l_A/4, l_B/4};
    Point p2 = {l_A - l_A/4, l_B - l_B/4};
    Point p3 = {l_A, l_B};
    calculate_y_cell_path(p0, p1, p2, p3, cell_path_y);
}



void calculate_y_cell_path(Point p0, Point p1, Point p2, Point p3, int64_t * y_points){
    
    //Calculate lines between points
    uint64_t i;

    #ifdef VERBOSE
    printf("Built on\n");
    printf("(%"PRIu64", %"PRIu64")\n", p0.x, p0.y);
    printf("(%"PRIu64", %"PRIu64")\n", p1.x, p1.y);
    printf("(%"PRIu64", %"PRIu64")\n", p2.x, p2.y);
    printf("(%"PRIu64", %"PRIu64")\n", p3.x, p3.y);
    #endif

    long double deltax, deltay, deltaerr, error;
    uint64_t y;

    //P0 to P1
    deltax = p1.x - p0.x;
    deltay = p1.y - p0.y;
    if(deltax != 0) deltaerr = fabsl(deltay/deltax); else deltaerr = 0;
    //printf("Deltas  x: %Le y: %Le Error: %Le\n", deltax, deltay, deltaerr);
    error = deltaerr - 0.5;
    y = p0.y;

    for(i=p0.x;i<p1.x;i++){
        y_points[i] = (int64_t) y;
        error = error + deltaerr;
        if(error >= 0.5){
            y++;
            error = error - 1;
        }
    }

    //P1 to P2

    deltax = p2.x - p1.x;
    deltay = p2.y - p1.y;
    if(deltax != 0) deltaerr = fabsl(deltay/deltax); else deltaerr = 0;
    //printf("Deltas  x: %Le y: %Le Error: %Le\n", deltax, deltay, deltaerr);
    error = deltaerr - 0.5;
    y = p1.y;

    for(i=p1.x;i<p2.x;i++){
        y_points[i] = (int64_t) y;
        error = error + deltaerr;
        if(error >= 0.5){
            y++;
            error = error - 1;
        }
    }
    
    //P2 to P3

    deltax = p3.x - p2.x;
    deltay = p3.y - p2.y;
    if(deltax != 0) deltaerr = fabsl(deltay/deltax); else deltaerr = 0;
    //printf("Deltas  x: %Le y: %Le Error: %Le\n", deltax, deltay, deltaerr);
    error = deltaerr - 0.5;
    y = p2.y;

    for(i=p2.x;i<p3.x;i++){
        y_points[i] = (int64_t) y;
        error = error + deltaerr;
        if(error >= 0.5){
            y++;
            error = error - 1;
        }
    }
    /*
    #ifdef VERBOSE
    for(i=0;i<p3.x;i++){
        printf("%"PRIu64" -> ", y_points[i]);
        if(i % 50 == 0) getchar();
    }
    
    #endif
    */

}
struct best_cell NW(char ** X, uint64_t Xstart, uint64_t Xend, char ** Y, uint64_t Ystart, uint64_t Yend, int64_t iGap, int64_t eGap, struct cell_f ** table, struct positioned_cell * mc, int show, int64_t * cell_path_y, long double * window, uint64_t * current_window_size, uint64_t nx, uint64_t ny){
    

    uint64_t i, j, j_prime;
    int64_t scoreDiagonal = INT64_MIN, scoreLeft = INT64_MIN, scoreRight = INT64_MIN, score = INT64_MIN, delta_dif_1_0, delta_dif_2_1, limit_left, limit_right, j_right_prime = 1, j_left_prime = 1, j_diag_prime = 1;

    struct best_cell bc;
    bc.c.score = INT64_MIN;
    bc.c.xpos = 0; bc.c.ypos = 0;

    // To pile up chars 
    char list_piled_X[nx];
    char list_piled_Y[ny];
    
    //The window size will be a +-15% of the square root of the product of lengths
    int64_t window_size = MIN(MAX_WINDOW_SIZE/2, (uint64_t) (*window * sqrtl((long double) Xend * (long double) Yend)));
    //printf("xlen: %"PRIu64", ylen: %"PRIu64" w-size: %"PRId64"\n", Xend, Yend, window_size);    
    *current_window_size = (uint64_t) window_size;

    //The limits to the window
    limit_left = 0;
    limit_right = 2*window_size + 1;
    if(limit_right > MAX_WINDOW_SIZE) limit_right = MAX_WINDOW_SIZE;
    
    struct positioned_cell mf;
    mf.score = INT64_MIN;
    

    //First row. iCounter serves as counter from zero
    //printf("..0%%");
    //Zero will always be
    pile_chars_up(X, Y, list_piled_X, list_piled_Y, nx, ny, 0, 0);
    table[0][0].score = compare_piled_up_chars(list_piled_X, list_piled_Y, nx, ny);
    //table[0][0].score = compare_letters(X[0], Y[0]);
    mc[0].score = table[0][0].score;
    mc[0].xpos = 0;
    mc[0].ypos = 0;

	for(i=1;i<Yend;i++){
        //table[0][i].score = (X[0] == Y[i]) ? POINT : -POINT;
        if(i < (uint64_t)(cell_path_y[0] + window_size)){
            pile_chars_up(X, Y, list_piled_X, list_piled_Y, nx, ny, 0, i);
            table[0][i].score = compare_piled_up_chars(list_piled_X, list_piled_Y, nx, ny) + iGap + (i-1)*eGap;
            //table[0][i].score = compare_letters(X[0], Y[i]) + iGap + (i-1)*eGap;
        } 
        //table[Xstart][i].xfrom = Xstart;
        //table[Xstart][i].yfrom = i;
        //Set every column max
        
        mc[i].score = compare_piled_up_chars(list_piled_X, list_piled_Y, nx, ny) + iGap + (i-1)*eGap;
        /*
        #ifdef VERBOSE
        printf("%02"PRId64" ", mc[i].score);
        #endif
        */
        mc[i].xpos = 0;
        mc[i].ypos = i;

    }
    /*
    #ifdef VERBOSE
    printf("\n");
    #endif
    */
    //Set row max
    mf.score = table[0][0].score;
    mf.xpos = 0;
    mf.ypos = 0;
    //Init j
    j = MAX(1,(cell_path_y[1] - window_size));

    //Go through full matrix
	for(i=1;i<Xend;i++){
        //Fill first rowcell
        if(cell_path_y[i-1]+window_size < cell_path_y[i]) return bc;
        //Conversion for the j-coordinate
        j_prime = 1;

        //table[i][0].score = (X[i] == Y[0]) ? POINT : -POINT;
        if(cell_path_y[i] - window_size <= 0){
            pile_chars_up(X, Y, list_piled_X, list_piled_Y, nx, ny, i, 0);
            table[i][0].score = compare_piled_up_chars(list_piled_X, list_piled_Y, nx, ny) + iGap + (i-1)*eGap;
            //table[i][0].score = compare_letters(X[i], Y[0]) + iGap + (i-1)*eGap;
            mf.score = table[i][0].score;
        }else{
            pile_chars_up(X, Y, list_piled_X, list_piled_Y, nx, ny, i, 0);
            mf.score = compare_piled_up_chars(list_piled_X, list_piled_Y, nx, ny) + iGap + (i-1)*eGap;
        }

        mf.xpos = i-1;
        mf.ypos = 0;

        delta_dif_1_0 = MAX(1, (cell_path_y[i] - window_size)) - MAX(1,(cell_path_y[i-1] - window_size)); //j-1
        if(i>1) delta_dif_2_1 = MAX(1, (cell_path_y[i-1] - window_size)) - MAX(1, (cell_path_y[i-2] - window_size)); //j-2

        /*
        #ifdef VERBOSE 
        printf("D1_0: %"PRId64" D2_1: %"PRId64"\n", delta_dif_1_0, delta_dif_2_1);
        #endif

        #ifdef VERBOSE
        printf("%02"PRId64" ", mf.score);
        #endif
        */
        //printf("Check on i: (%"PRIu64") from - to (%"PRIu64", %"PRIu64")\n", i, 0L, Xend);
        //printf("I will go from %"PRIu64" to %"PRIu64"\n", (uint64_t) MAX(1,(cell_path_y[i] - window_size)), (uint64_t) MIN(Yend,(cell_path_y[i] + window_size)));
        //getchar();

        /*
        #ifdef VERBOSE
        int64_t r;
        for(r=0;r<MAX(0,(cell_path_y[i] - window_size)); r++){
            printf("  ");
        }
        #endif
        */

        

        

        for(j=(uint64_t)MAX(1,((int64_t)cell_path_y[i] - window_size));j<(uint64_t)MIN((int64_t)Yend,(cell_path_y[i] + window_size)) && j_prime < (uint64_t) limit_right;j++){
            //printf("Doing on : (%"PRIu64",%"PRIu64" and jprime=%"PRIu64"\n", i,j,j_prime);
            //Check if max in row has changed
            //if(j > MAX(1, cell_path_y[i-1] - window_size +1) && mf.score <= table[i][j-2].score){

            //Calculate the real j position in the windowed table
            j_left_prime = ((int64_t)j_prime - (2 - delta_dif_1_0));
            //j_diag_prime = ((int64_t)j_prime - (1 - delta_dif_1_0));
            j_diag_prime = ((int64_t)j_prime - (1 - delta_dif_1_0));
            if(i > 1){
                j_right_prime = ((int64_t)j_prime - (1 - (delta_dif_1_0 + delta_dif_2_1)));
            }

            if(j > (uint64_t) MAX(1, cell_path_y[i-1] - window_size +1) && j < (uint64_t) MIN( (int64_t) Yend,(cell_path_y[i-1] + window_size)) && j_left_prime < limit_right && table[i-1][j_left_prime].score >= mf.score){
                //mf.score = table[i-1][j-2].score;
                mf.score = table[i-1][j_left_prime].score;
                mf.xpos = i-1;
                mf.ypos = j-2;
                if(table[i-1][j_left_prime].score == INT64_MIN){ printf("A: mf.x\t%"PRIu64"\tmf.y\t%"PRIu64"\ts%"PRId64"\n", mf.xpos, mf.ypos, mf.score); printf("@[%"PRIu64", %"PRIu64"] with j_prime: %"PRIu64", wsize: %"PRIu64", cp[i-1]=%"PRId64", cp[i]=%"PRId64"\n", i, j, j_prime, 2*window_size, cell_path_y[i-1], cell_path_y[i]); getchar(); }
                
            }
            //printf("RowMax: %"PRId64"@(%"PRIu64", %"PRIu64")\t", mf.score, mf.xpos, mf.ypos);
            
            //score = (X[i] == Y[j]) ? POINT : -POINT;
            pile_chars_up(X, Y, list_piled_X, list_piled_Y, nx, ny, i, j);
            score = compare_piled_up_chars(list_piled_X, list_piled_Y, nx, ny);
            //score = compare_letters(X[i], Y[j]);

            if(i > 1 && j >= 1 && j-1 >= (uint64_t) MAX(1,(cell_path_y[i-2] - window_size)) && j-1 < (uint64_t) MIN((int64_t)Yend,(cell_path_y[i-2] + window_size)) && j_right_prime >= limit_left && j_right_prime < limit_right && table[i-2][j_right_prime].score >= mc[j-1].score ){
                //mc[j-1].score = table[i-2][j-(1+j_prime)].score;
                //Should be the j_prime we had at cell_path_y
                //MAX(1,(cell_path_y[i] - window_size));j<MIN(Yend,(cell_path_y[i] + window_size))
                
                mc[j-1].score = table[i-2][j_right_prime].score;
                mc[j-1].xpos = i-2;
                mc[j-1].ypos = j-1;

                if(table[i-2][j_right_prime].score == INT64_MIN){ printf("A: j-1\t%"PRIu64"\tmc.xpos\t%"PRIu64"\ts%"PRId64"\n", j-1, mc[j-1].xpos, mc[j-1].score); printf("@[%"PRIu64", %"PRIu64"] with j_prime: %"PRIu64", wsize: %"PRIu64", cp[i-1]=%"PRId64", cp[i]=%"PRId64"\n", i, j, j_prime, 2*window_size, cell_path_y[i-1], cell_path_y[i]); getchar(); }
    
            }


            //Precondition: Upper row needs to reach up to diagonal
            //if((cell_path_y[i-1]+window_size) >= j-1){
            if(j-1 >= MAX(0, (cell_path_y[i-1]-window_size)) && (uint64_t) (cell_path_y[i-1]+window_size) >= j-1 && j_diag_prime >= limit_left && j_diag_prime < limit_right && j_diag_prime < cell_path_y[i-1]+window_size){
                //scoreDiagonal = table[i-1][j-1].score + score;
                //printf("prevdiag: %"PRId64"\n", table[i-1][j_diag_prime].score);
                scoreDiagonal = table[i-1][j_diag_prime].score + score;                
                if(table[i-1][j_diag_prime].score == INT64_MIN){ printf("A: i-1\t%"PRIu64"\tj_diag\t%"PRIu64"\ts%"PRId64"\n", i-1, j_diag_prime, table[i-1][j_diag_prime].score); printf("@[%"PRIu64", %"PRIu64"] with j_prime: %"PRIu64", wsize: %"PRIu64", cp[i-1]=%"PRId64", cp[i]=%"PRId64"\n", i, j, j_prime, 2*window_size, cell_path_y[i-1], cell_path_y[i]); getchar(); }
            }else{
                scoreDiagonal = INT64_MIN;
            }
            
            if(i>=1 && j>1){
                scoreLeft = mf.score + iGap + (j - (mf.ypos+2))*eGap + score;
                
                if(mf.score == INT64_MIN){ printf("A: mf.x\t%"PRIu64"\tmf.y\t%"PRIu64"\ts%"PRId64"\n", mf.xpos, mf.ypos, mf.score); printf("@[%"PRIu64", %"PRIu64"] with j_prime: %"PRIu64", wsize: %"PRIu64", cp[i-1]=%"PRId64", cp[i]=%"PRId64"\n", i, j, j_prime, 2*window_size, cell_path_y[i-1], cell_path_y[i]); getchar(); }
            }else{
                scoreLeft = INT64_MIN;
            }

            if(j>=1 && i>1){
                scoreRight = mc[j-1].score + iGap + (i - (mc[j-1].xpos+2))*eGap + score;
                //if(scoreRight == -12) printf("MC: %"PRId64", From: %"PRIu64", %"PRIu64"->", mc[j-1].score, mc[j-1].xpos, mc[j-1].ypos);
                
                if(mc[j-1].score == INT64_MIN){ printf("A: j-1\t%"PRIu64"\tmc.xpos\t%"PRIu64"\ts%"PRId64"\n", j-1, mc[j-1].xpos, mc[j-1].score); printf("@[%"PRIu64", %"PRIu64"] with j_prime: %"PRIu64", wsize: %"PRIu64", cp[i-1]=%"PRId64", cp[i]=%"PRId64"\n", i, j, j_prime, 2*window_size, cell_path_y[i-1], cell_path_y[i]); getchar(); }
            }else{
                scoreRight = INT64_MIN;
            }
            
            //Choose maximum
            /*
            #ifdef VERBOSE
            printf("The game starts at %"PRId64"\n", MAX(0, cell_path_y[i] - window_size));
            printf("from %c %c and I get to %"PRIu64" while j=%"PRIu64"\n", X[i], Y[j], j_prime, j);
            printf("j_prime: %"PRId64"\n", j_prime);
            printf("j_diag_prime: %"PRId64" limits[%"PRId64", %"PRId64"]\n", j_diag_prime, limit_left, limit_right);
            printf("Score DIAG: %"PRId64"; LEFT: %"PRId64"; RIGHT: %"PRId64"\n", scoreDiagonal, scoreLeft, scoreRight);
            printf("currmf: %"PRId64" mc: %"PRId64"\n", mf.score, mc[j-1].score);
            #endif
            */
            
            if(scoreDiagonal >= scoreLeft && scoreDiagonal >= scoreRight){
                //Diagonal
                
                //fprintf(stdout, "The JPRIME: %"PRId64" actual pos: %"PRIu64"\n", j_prime, j); getchar();
                table[i][j_prime].score = scoreDiagonal;
                table[i][j_prime].xfrom = i-1;
                table[i][j_prime].yfrom = j-1;
                
                                
            }else if(scoreRight > scoreLeft){
                table[i][j_prime].score = scoreRight;
                table[i][j_prime].xfrom = mc[j-1].xpos;
                table[i][j_prime].yfrom = mc[j-1].ypos;
                
            }else{
                //printf("Scores %"PRId64", %"PRId64", %"PRId64"\n", scoreDiagonal, scoreLeft, scoreRight);
                table[i][j_prime].score = scoreLeft;
                table[i][j_prime].xfrom = mf.xpos;
                table[i][j_prime].yfrom = mf.ypos;
            }
        
        
            //check if column max has changed
            //New condition: check if you filled i-2, j-1
            
            if(i == Xend-1 || j == Yend-1){

                if(i == Xend-1 && j != Yend-1){
            		table[i][j_prime].score = table[i][j_prime].score + iGap + (Yend - j)*eGap;
            	}else if(j == Yend-1 && i != Xend-1){
            		table[i][j_prime].score = table[i][j_prime].score + iGap + (Xend - i)*eGap;
            	}
                //Check for best cell
                if(table[i][j_prime].score >= bc.c.score){ 
                    
                    /*
                    if(i == 798 && j == 1052){ // yields 799, 2497
                        printf("in position @ jprime= %"PRIu64" cellpaths [i-1, i] are %"PRId64", %"PRId64"\n", j_prime, cell_path_y[i-1], cell_path_y[i]);
                        printf("Scores %"PRId64", %"PRId64", %"PRId64"\n", scoreDiagonal, scoreLeft, scoreRight);
                        printf("score comes from %"PRIu64", %"PRIu64",\n", mc[j-1].xpos, mc[j-1].ypos);
                        printf("IDlengths: %"PRIu64", %"PRIu64"\n", Xend, Yend);
                        
                        //exit(-1);
                    }
                    */
                    
                    bc.c.score = table[i][j_prime].score; bc.c.xpos = i; bc.c.ypos = j; bc.j_prime = j_prime; 
                }
                //bc.c.score = table[i][j_prime].score; bc.c.xpos = i; bc.c.ypos = j; bc.j_prime = j_prime;
            }
            #ifdef VERBOSE
            //printf("Put score: %"PRId64"\n\n", table[i][j_prime].score);
            /*
            printf("(%"PRId64")%02"PRId64" ", j_diag_prime, table[i][j_prime].score); //printf("->(%"PRIu64", %"PRIu64")", i, j); printf("[%c %c]", X[i], Y[j]);
            */
            //if(scoreDiagonal >= scoreLeft && scoreDiagonal >= scoreRight) printf("*\t");
            //else if(scoreRight > scoreLeft) printf("{\t"); else printf("}\t");
            //getchar();
            #endif
            j_prime++;
        }
        /*
        #ifdef VERBOSE
        printf("\n");
        getchar();
        #endif
        */
    }
        
    return bc;
}

void backtrackingNW(char * X, uint64_t Xstart, uint64_t Xend, char * Y, uint64_t Ystart, uint64_t Yend, struct cell_f ** table, char * rec_X, char * rec_Y, struct best_cell * bc, uint64_t * ret_head_x, uint64_t * ret_head_y, int64_t * cell_path_y, uint64_t window_size, BasicAlignment * ba){
    uint64_t curr_x, curr_y, prev_x, prev_y, head_x, head_y, limit_x, limit_y;
    int64_t k, j_prime, delta_diff = 0;

    limit_x = 2*MAX(Xend, Yend);
    limit_y = limit_x;

    head_x = limit_x;
    head_y = limit_y;

    curr_x = bc->c.xpos;
    curr_y = bc->c.ypos;
    #ifdef VERBOSE
    printf("Optimum : %"PRIu64", %"PRIu64"\n", curr_x, curr_y);
    #endif
    prev_x = curr_x;
    prev_y = curr_y;
    int show = 0;
    printf("Optimum : %"PRIu64", %"PRIu64"\n", curr_x, curr_y);
    printf("Xend: %"PRIu64", YEnd: %"PRIu64"\n", Xend-1, Yend-1);

    for(k=(int64_t)Xend-1; k>(int64_t)curr_x; k--) rec_X[head_x--] = '-';
    for(k=(int64_t)Yend-1; k>(int64_t)curr_y; k--) rec_Y[head_y--] = '-';

    j_prime = bc->j_prime;
    unsigned char first_track = 1;
    
    while(curr_x > 0 && curr_y > 0){

        
        if(first_track == 0){
            delta_diff = MAX(1, cell_path_y[prev_x] - (int64_t) window_size) - MAX(1, cell_path_y[curr_x] - (int64_t)window_size); //j-1
            j_prime = MAX(0, j_prime - (int64_t)(prev_y - curr_y) + (int64_t) delta_diff);

            prev_x = curr_x;
            prev_y = curr_y;
            
            
            #ifdef VERBOSE
            /*
            //printf("Jprime: %"PRId64" :DELTADIF:%"PRId64"\n", j_prime, delta_diff);
            printf("[%c %c]", X[prev_x], Y[prev_y]);
            printf("(%"PRIu64", %"PRIu64") ::: \n", curr_x, curr_y);
            //printf("(%"PRIu64", %"PRIu64") ::: \n", prev_x, prev_y);
            //printf("cellp Prev: %"PRId64" Post: %"PRId64"\n", cell_path_y[prev_x], cell_path_y[curr_x]);
            //printf("the difs? %"PRId64" the other: %"PRId64"\n", MAX(0, cell_path_y[prev_x] - (int64_t) window_size), MAX(0, cell_path_y[curr_x] - (int64_t)window_size));
            getchar();
            */
            #endif

        }

        curr_x = table[prev_x][j_prime].xfrom;
        curr_y = table[prev_x][j_prime].yfrom;
        first_track = 0;
        

        

        if((curr_x == (prev_x - 1)) && (curr_y == (prev_y -1))){
            //Diagonal case
            //printf("DIAG\n");
            if(head_x == 0 || head_y == 0) goto exit_point;
            rec_X[head_x--] = (char) X[prev_x];
            rec_Y[head_y--] = (char) Y[prev_y];
            ba->length++;
        }else if((prev_x - curr_x) > (prev_y - curr_y)){
            //Gap in X
            //printf("Gap X\n");
            if(head_x == 0 || head_y == 0) goto exit_point;
            if(bc->c.xpos != prev_x && bc->c.ypos != prev_y){
                rec_Y[head_y--] = Y[prev_y];
                rec_X[head_x--] = X[prev_x];
            }else{
                rec_Y[head_y--] = '-';
                rec_X[head_x--] = X[prev_x];
            }
            ba->length++;
            
            for(k=(int64_t)prev_x-1;k>(int64_t)curr_x;k--){
                if(head_x == 0 || head_y == 0) goto exit_point;
                #ifdef VERBOSE 
                if(head_x == 0 || head_y == 0){
                    printf("%"PRIu64" %"PRIu64" and prevs are %"PRIu64" %"PRIu64"\n", head_x, head_y, prev_x, prev_y);
                    printf("origin is %"PRIu64", %"PRIu64"\n", bc->c.xpos, bc->c.ypos);
                    uint64_t z;
                    for(z=head_x;z<limit_x;z++){
                        fprintf(stdout, "%c", (char) rec_X[z]);
                    }
                    printf("\n");
                    for(z=head_y;z<limit_y;z++){
                        fprintf(stdout, "%c", (char) rec_Y[z]);
                    }
                    getchar();
                }
                #endif
                rec_Y[head_y--] = '-';
                rec_X[head_x--] = (char) X[k];
                ba->length++;
                ba->egaps++;
            }
            ba->igaps += 1;
            ba->egaps--;
        }else{
            //Gap in Y
            //printf("GAP Y\n");
            //10, 0, 401, 281
            if(head_x == 0 || head_y == 0) goto exit_point;
            if(bc->c.xpos != prev_x && bc->c.ypos != prev_y){
                rec_Y[head_y--] = Y[prev_y];
                rec_X[head_x--] = X[prev_x];
            }else{
                rec_Y[head_y--] = Y[prev_y];
                rec_X[head_x--] = '-';
            }
            ba->length++;

            for(k=(int64_t)prev_y-1;k>(int64_t)curr_y;k--){
                if(head_x == 0 || head_y == 0) goto exit_point;
                #ifdef VERBOSE 
                if(head_x == 0 || head_y == 0){
                    printf("%"PRIu64" %"PRIu64" and prevs are %"PRIu64" %"PRIu64"\n", head_x, head_y, prev_x, prev_y);
                    printf("origin is %"PRIu64", %"PRIu64"\n", bc->c.xpos, bc->c.ypos);
                    uint64_t z;
                    for(z=head_x;z<limit_x;z++){
                        fprintf(stdout, "%c", (char) rec_X[z]);
                    }
                    printf("\n");
                    for(z=head_y;z<limit_y;z++){
                        fprintf(stdout, "%c", (char) rec_Y[z]);
                    }
                    getchar();
                }
                #endif
                rec_X[head_x--] = '-';
                rec_Y[head_y--] = (char) Y[k];
                ba->length++;
                ba->egaps++;
            }
            
            ba->igaps += 1;
            ba->egaps--;
        }
        
    }
    if(curr_x == 0 && curr_y == 0 && (curr_x == (prev_x - 1)) && (curr_y == (prev_y -1))){
        rec_X[head_x--] = (char) X[curr_x];
        rec_Y[head_y--] = (char) Y[curr_y];
        ba->length++;
    }
    
    exit_point:
    if(show == 1)fprintf(stdout, "%"PRIu64", %"PRIu64"\n", head_x, head_y);
    uint64_t huecos_x = 0, huecos_y = 0;
    k=(int64_t)curr_x-1;
    while(k>=0){ if(head_x == 0) break; rec_X[head_x--] = '-'; huecos_x++;  k--; }
    k=(int64_t)curr_y-1;
    while(k>=0){ if(head_y == 0) break; rec_Y[head_y--] = '-'; huecos_y++; k--; }
    
    if(show == 1)fprintf(stdout, "%"PRIu64", %"PRIu64"\n", head_x, head_y);

    if(huecos_x >= huecos_y){
        while(huecos_x > 0) { if(head_y == 0) break; rec_Y[head_y--] = ' '; huecos_x--;}
    }else{
        while(huecos_y > 0) { if(head_x == 0) break; rec_X[head_x--] = ' '; huecos_y--;}
    }

    if(show == 1){
        fprintf(stdout, "%"PRIu64", %"PRIu64"\n", head_x, head_y);
        fprintf(stdout, "%"PRIu64", %"PRIu64"\n", 2*Xend, 2*Yend);
        uint64_t k;
        for(k=head_x;k<limit_x;k++){
            fprintf(stdout, "%c", (char) rec_X[k]);
        }
        printf("\n");
        for(k=head_y;k<limit_y;k++){
            fprintf(stdout, "%c", (char) rec_Y[k]);
        }
        printf("\n");
        getchar();
    }

    *ret_head_x = head_x;
    *ret_head_y = head_y;
    #ifdef VERBOSE
    printf("hx hy: %"PRIu64", %"PRIu64"\n", head_x, head_y);
    #endif
}





/*
Calculates NW table with two rows and stores a cellpath of scores, identities, gaps and starting and ending positions
*/

struct cell NWscore2rows(char * X, uint64_t Xstart, uint64_t Xend, char * Y, uint64_t Ystart, uint64_t Yend, int64_t iGap, int64_t eGap, struct cell * mc, struct cell * f0, struct cell * f1){
    
    uint64_t i, j, k, iCounter=0, jCounter, currEgapR = 0, currEgapG = 0;
	int64_t scoreDiagonal = INT64_MIN, scoreLeft = INT64_MIN, scoreRight = INT64_MIN, score = INT64_MIN;
	struct cell * faux;
	
    struct cell mf;
    

    if(mc == NULL || f0 == NULL || f1 == NULL){
    	printf("Could not allocate memory.\n");
    	exit(-1);
    }
    
    
    mf.score = INT64_MIN;
    

    //First row. iCounter serves as counter from zero
    //printf("..0%%");

    f0[iCounter].score = PAM[valOfNucl(X[Xstart])][valOfNucl(Y[Ystart])];
    f0[iCounter].igaps = 0;
    f0[iCounter].egaps = 0;
    f0[iCounter].ident = (((f0[iCounter].score) > (0)) ? (1) : (0));
    f0[iCounter].xs = Xstart;
    f0[iCounter].ys = Ystart;
    f0[iCounter].xe = Xstart;
    f0[iCounter].ye = Ystart;
    mc[iCounter] = f0[iCounter];
    //printf("    %03"PRId64" ", f0[iCounter].score);
    iCounter++;


    printf("Align between %"PRIu64", %"PRIu64"\n", Ystart+1, Yend);
    for(i=Ystart+1;i<Yend;i++){

        f0[iCounter].score = PAM[valOfNucl(X[Xstart])][valOfNucl(Y[i])] + iGap + ((i-Ystart)-1)*eGap;
    	f0[iCounter].igaps = 1;
    	f0[iCounter].egaps = i - Ystart;
    	f0[iCounter].ident = (((f0[iCounter].score) > (0)) ? (1) : (0));
    	f0[iCounter].xs = Xstart;
    	f0[iCounter].ys = i;
    	f0[iCounter].xe = Xstart;
    	f0[iCounter].ye = i;

    	//Set every column max
    	mc[iCounter] = f0[iCounter];
    	//goodPrint(f0[iCounter].score);
        //printf("%03"PRId64" ", f0[iCounter].score);
    	
    	iCounter++;
	}
	
	//Set row max
	mf = f0[0];
    //printf("\n");
	
	iCounter=1;

	//Go through full matrix with 2 rows
	for(i=Xstart+1;i<Xend;i++){
		//Fill first rowcell
		//printf("ROW: %"PRId64" ||",i);
		//printf("%.2f \n", ((float)i/Xend));
		
        f1[0].score = PAM[valOfNucl(X[i])][valOfNucl(Y[Ystart])] + iGap + ((i-Xstart)-1)*eGap;
		f1[0].xs = i;
		f1[0].ys = Ystart;
		f1[0].ident = (((PAM[valOfNucl(X[i])][valOfNucl(Y[Ystart])]) > (0)) ? (1) : (0));
        f1[0].igaps = 1;
        f1[0].egaps = i - (Xstart + 1);		
		f1[0].xe = i;
		f1[0].ye = Ystart;

		mf = f0[0];

		//goodPrint(f1[0].score);
        //printf("%03"PRId64" ", f1[0].score);
		jCounter=1;
		for(j=Ystart+1;j<Yend;j++){
		//for(j=redir[i][0]+1;j<redir[i][1];j++){
			//Check if max in row has changed
			if(jCounter > 1 && mf.score <= f0[jCounter-2].score){
				mf = f0[jCounter-2];
				mf.xe = i-1;
				mf.ye = j-2;
			}
			
			score = PAM[valOfNucl(X[i])][valOfNucl(Y[j])];
			scoreDiagonal = f0[jCounter-1].score + score;
			if(jCounter>1 && iCounter>=1){
				scoreLeft = mf.score + iGap + (j - (mf.ye+2))*eGap + score;
				currEgapR = (j - (mf.ye+2));
				}else{
					scoreLeft = INT64_MIN;
				}
				
			if(iCounter>1 && jCounter>=1){
				scoreRight = mc[jCounter-1].score + iGap + (i - (mc[j-1].xe+2))*eGap + score;
				currEgapG =  (i - (mc[j-1].xe+2));
				}else{
					scoreRight = INT64_MIN;
				}
			
			//Choose maximum
			//f1[jCounter] = max(max(scoreDiagonal,scoreLeft),scoreRight);
			
			if(scoreDiagonal >= MAX(scoreLeft, scoreRight)){
				//Diagonal
				f1[jCounter] = f0[jCounter-1];
				f1[jCounter].score = scoreDiagonal;
				if(PAM[valOfNucl(X[i])][valOfNucl(Y[j])] > 0) f1[jCounter].ident += 1;
								
			}else if(scoreRight >= scoreLeft){
				//Gap in genome
				f1[jCounter] = mc[jCounter-1];
				f1[jCounter].score = scoreRight;
				f1[jCounter].igaps += 1;
				f1[jCounter].egaps += currEgapG;
				
			}else{
				//Gap in read
				f1[jCounter] = mf;
				f1[jCounter].score = scoreLeft;
				f1[jCounter].igaps += 1;
				f1[jCounter].egaps += currEgapR;
			}

            if(iCounter>1 && jCounter>=1 && f0[jCounter].score >= mc[jCounter-1].score){
                mc[jCounter-1].score = f0[jCounter].score;
                mc[jCounter-1].xe = i-2;
                mc[jCounter-1].ye = j-1;
            }




			//Update movement
			f1[jCounter].xe = i;
			f1[jCounter].ye = j;
			//goodPrint(f1[jCounter].score);
            //printf("%03"PRId64" ", f1[jCounter].score);
			jCounter++;
		}
        //printf("\n");
        /*
		kCounter=0;
		for(k=Ystart;k<Yend;k++){
			//Update column maximum at j
			if(mc[kCounter].score <= f0[kCounter].score){
				mc[kCounter] = f0[kCounter];
				mc[kCounter].xe = i-1;
				mc[kCounter].ye = kCounter;
			}

			kCounter++;
		}
        */
		//Switch rows
		
		faux = f0;
		f0 = f1;
		f1 = faux;
		
		iCounter++;
	}



    int64_t bestScore=INT64_MIN, bestId = 0;
    for(k=Ystart;k<Yend;k++){
    	//printf("start (%"PRIu64",%"PRIu64") end (%"PRIu64",%"PRIu64") score [%"PRId64"] gaps [%"PRIu64"] ident [%"PRIu64"]\n", f0[k].xs, f0[k].ys, f0[k].xe, f0[k].ye, f0[k].score, f0[k].gaps, f0[k].ident);
    	if(f0[k].score >= bestScore){
    		bestScore = f0[k].score;
    		//printf("start (%"PRIu64",%"PRIu64") end (%"PRIu64",%"PRIu64") score [%"PRId64"] gaps [%"PRIu64"] ident [%"PRIu64"]\n", f0[k].xs, f0[k].ys, f0[k].xe, f0[k].ye, f0[k].score, f0[k].gaps, f0[k].ident);
    		bestId = k;
    	}
    }
    return f0[bestId];
}








uint64_t build_multiple_alignment(char ** reconstruct_X, char ** reconstruct_Y, char ** my_x, char ** my_y, uint64_t nx, uint64_t ny, struct cell_f ** table, struct positioned_cell * mc, char * writing_buffer_alignment, uint64_t xlen, uint64_t ylen, int64_t * cell_path_y, long double * window, int64_t iGap, int64_t eGap, BasicAlignment * ba, char * aux){
 

    //Do some printing of alignments here
    uint64_t i, j, k, curr_window_size;
    
    struct best_cell bc = NW(my_x, 0, xlen, my_y, 0, ylen, iGap, eGap, table, mc, 0, cell_path_y, window, &curr_window_size, nx, ny);

    uint64_t max_len = MAX(xlen, ylen);
    uint64_t pos_longest_X = get_id_of_longest_sequence(my_x, nx);
    uint64_t pos_longest_Y = get_id_of_longest_sequence(my_y, ny);
    char * reference;
    if(strlen(my_x[pos_longest_X]) > strlen(my_y[pos_longest_Y])) reference = my_x[pos_longest_X]; else reference = my_y[pos_longest_Y];

    printf("Max len for round %"PRIu64"\n", max_len);
    for(k=0;k<nx;k++){
        backtrackingNW(my_x[k], 0, xlen, reference, 0, max_len, table, reconstruct_X[k], aux, &bc, &i, &j, cell_path_y, curr_window_size, ba);
        i++; j++;
        //fprintf(stdout, "X:%"PRIu64" -> %s\n", k, &reconstruct_X[k][i]);
        memcpy(&my_x[k][0], &reconstruct_X[k][i], strlen(&reconstruct_X[k][i]));
    }

    uint64_t final_y_len;
    for(k=0;k<ny;k++){
        backtrackingNW(reference, 0, max_len, my_y[k], 0, ylen, table, aux, reconstruct_Y[k], &bc, &i, &j, cell_path_y, curr_window_size, ba);
        i++; j++;    
        //fprintf(stdout, "Y:%"PRIu64" -> %s\n", k, &reconstruct_Y[k][j]);
        final_y_len = strlen(&reconstruct_Y[k][j]);
        memcpy(&my_y[k][0], &reconstruct_Y[k][j], final_y_len);
        
    }
    return final_y_len;
}