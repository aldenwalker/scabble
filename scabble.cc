/**************************************************************
*                                                             *
* scabble - calculates the unit scl ball in a 2-plane spanned *
*           by two chains                                     *
*                                                             *  
*                                                             *
*   Copyright Danny Calegari and                              *
*    Alden Walker 2010                                        *
*                                                             *
*                                                             *
*  Released under the GPL license                             *
*                                                             *
**************************************************************/
  

#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <vector>
#include <ctype.h>

#include "word.h"
#include "rational.h"
#include "scabble.h"
#include "lp.h"



using namespace std;


/******************************************************************************/
/*  generate a list of all possible arcs                                      */
/******************************************************************************/
void generate_arcs(vector<arc>* arc_list, vector<string>& w, int WORD) {
  int g,h,i,j;
  int a,b;
  arc tempArc;
  arc_list->resize(0);
  for(g=0;g<WORD;g++){  
    for(i=0;i< (int) w[g].length();i++){
      a=(int) w[g][i];    
      for(h=g;h<WORD;h++){
        if(h==g){
          for(j=i+1;j< (int) w[h].length();j++){
            b=(int) w[h][j];
            if((32+a-b)%64==0){    // are letters inverse?
              tempArc.first = i;
              tempArc.last = j;
              tempArc.first_word=g;
              tempArc.last_word=h;
              arc_list->push_back(tempArc);
              tempArc.last = i;
              tempArc.first = j;
              tempArc.first_word = h;
              tempArc.last_word = g;
              arc_list->push_back(tempArc);
            };    
          };
        } else {
          for(j=0;j< (int) w[h].length();j++){
            b=(int) w[h][j];
            if((32+a-b)%64==0){    // are letters inverse?
              tempArc.first = i;
              tempArc.last = j;
              tempArc.first_word=g;
              tempArc.last_word=h;
              arc_list->push_back(tempArc);
              tempArc.last = i;
              tempArc.first = j;
              tempArc.first_word = h;
              tempArc.last_word = g;
              arc_list->push_back(tempArc);
            };              
          };
        };
      };      
    };
  };

};


/******************************************************************************/
/* generate the polygon list                                                  */
/******************************************************************************/
void generate_polygons(vector<string> w,
                       vector<polygon> &polygon_list,  
                       vector<arc>  &arc_list, 
                       int maxjun) {
                       
  int i,j,k,size,add,close,appeared;
  polygon testpoly;
  testpoly.arc.resize(maxjun);
  int arc_list_length = arc_list.size();

  vector<int> word_length(w.size());
  
  for (i=0; i<(int)w.size(); i++) {
    word_length[i] = w[i].size();
  }

  polygon_list.resize(0);

  for(i=0;i<arc_list_length-1;i++){  // i is initial index of hypothetical polygon
    testpoly.arc[0]=i;
    j=(i+1);
    size=1;
    while(size>0){

      add=1;    // haven't yet added
      close=1;  // haven't yet closed
      
      // does arc j glue up to end of last arc of testpoly?
      if(arc_list[j].first_word == arc_list[testpoly.arc[size-1]].last_word && 
        (arc_list[j].first - arc_list[testpoly.arc[size-1]].last-1)%word_length[arc_list[j].first_word]==0){    
        appeared=1;
        for(k=0;k<size;k++){
          if(j==testpoly.arc[k]){    //has arc j already appeared in testpoly?
            appeared=0;
          }
        }
      
        if(appeared==1){          // if j has not appeared
          testpoly.arc[size]=j;    // then add it to testpoly
          add=0;                  // note that we have added an arc
          
          // does it close up?
          if(arc_list[j].last_word == arc_list[testpoly.arc[0]].first_word && 
            (arc_list[testpoly.arc[0]].first-arc_list[j].last-1)%word_length[arc_list[j].last_word]==0){    
            testpoly.size = size+1;
            polygon_list.push_back(testpoly);
            close=0;        // note that we have closed a polygon
          }
        }
      }
  
  
      if(add==0 && close==1){
        size++;
        j=i+1;
        if(size>=maxjun){
          j=arc_list_length;
        }
      } else {
        j++;    // increment j
      }
      while(j>=arc_list_length){
        size--;
        j=testpoly.arc[size]+1;
      }
    }
  }

}




/*****************************************************************************
This function finds builds the polyhedron which is the preimage under delta
of the line segment between (a,b) and (c,d) and minimizes chi/2 over this
polyhedron.  This point (scaled) gives a vertex of the unit ball.  The function
then returns this vertex.
in: the chains, the arcs, polygons, and weights
   the constraints in is a list of the basic constraints (that everything 
   is positive and the arc equality constraints), since there's no need to 
   regenerate all that.
out: a rational 2d vector
note: this method could be generalized to higher dimensions -- just restrict to
      the affine hyperplane containing the input points
      
note 2: The input points should have scl = 1; in this case, if we find the resulting
        min scl is 1, this function will return the first point -- in this way, we know
        if scl is linear between the two input points
******************************************************************************/
vector<rational> min_point_on_line(vector<vector<string> >&chains, 
                                   vector<arc>& arc_list,
                                   vector<polygon>& polygon_list,
                                   vector<int>& weights,
                                   RatMat* constraints,
                                   vector<int>& equalityType,
                                   vector<vector<rational> >& points,
                                   scallop_lp_solver solver,
                                   int VERBOSE) {
  int i,j,k;
  int firstWordNumber;
  rational coef;
  int numChains = (int)chains.size();
  mpq_t entry;
  mpq_init(entry);
  
  //find the hyperplane which the points lie on
  //I'm going to specialize this for 2d for the moment
  vector<rational> hyperplaneNormal(2, rational(0,1));
  vector<rational> hyperplaneParallel(2, rational(0,1));
  rational hyperplaneValue;          //this is the value both points have
  vector<rational> parallelVal(2, rational(0,1));             //this is the value the parallel function has
  rational temp;
  for (i=0; i<2; i++) {
    hyperplaneParallel[i] = points[0][i] - points[1][i];  
  }
  for (i=0; i<2; i++) {
    parallelVal[i] =  hyperplaneParallel[0]*points[i][0] + hyperplaneParallel[1]*points[i][1] ;
  }
  //the normal is perpendicular
  hyperplaneNormal[0] = -hyperplaneParallel[1];
  hyperplaneNormal[1] = hyperplaneParallel[0];
  
  hyperplaneValue = (hyperplaneNormal[0]*points[0][0]) + (hyperplaneNormal[1]*points[0][1]);
  
  if (VERBOSE==1){
    cout << "Input points: (" <<  points[0][0] << ", " << points[0][1] << "); (" << points[1][0] << ", " << points[1][1] << ")\n";
    cout << "Hyperplane parallel: " <<  hyperplaneParallel[0] << ", " << hyperplaneParallel[1] << "\n";
    cout << "Hyperplane normal: " <<  hyperplaneNormal[0] << ", " << hyperplaneNormal[1] << "\n";
    cout << "Hyperplane value: " << hyperplaneValue << "\n";
  }
  
  //Now, the new equations we need to put in are (normal) . (delta_g, delta_h) == hyperplaneValue
  //and min(parallelvals) <= (parallel).(delta_g, delta_h) <= max(parallelvals)
  
  //note this is 3 new constraints; we will add these rows to constraints, and
  //get rid of them later
  RatMat_change_num_rows(constraints, constraints->nR+3);
  
  //first, the equality constraint
  for (i=0; i<(int)polygon_list.size(); i++) {  //for each poly
 
    coef = rational(0,1);                      //this will be the coefficient
 
    for (j=0; j<polygon_list[i].size; j++) {  //go through the edges
      firstWordNumber = 0;                    //this is the word index of the first word in the chain
      for (k=0; k<numChains; k++) {
        if (arc_list[polygon_list[i].arc[j]].first_word == firstWordNumber &&
            arc_list[polygon_list[i].arc[j]].first == 0) {                        //if we're looking at the first word in a chain, first letter
          coef = (coef) +  (hyperplaneNormal[k]/weights[firstWordNumber]) ;
        }
        firstWordNumber+=chains[k].size();
      }
    }
    coef.get_mpq(entry);
    RatMat_set(constraints, constraints->nR-3, i, entry);
  }
  hyperplaneValue.get_mpq(entry);
  RatMat_set(constraints, constraints->nR-3, constraints->nC-1, entry);
  equalityType.push_back(0);
  
  
  //next, the two inequality constraints
  if (parallelVal[1] < parallelVal[0]) {
    temp = parallelVal[0];
    parallelVal[0] = parallelVal[1];
    parallelVal[1] = temp;
  }

  for (i=0; i<(int)polygon_list.size(); i++) {  //for each poly
 
    coef = rational(0,1);                      //this will be the coefficient
 
    for (j=0; j<polygon_list[i].size; j++) {  //go through the edges
      firstWordNumber = 0;                    //this is the word index of the first word in the chain
      for (k=0; k<numChains; k++) {
        if (arc_list[polygon_list[i].arc[j]].first_word == firstWordNumber &&
            arc_list[polygon_list[i].arc[j]].first == 0) {                        //if we're looking at the first word in a chain, first letter
          coef = coef +  (hyperplaneParallel[k]/weights[firstWordNumber]) ;
        }
        firstWordNumber+=chains[k].size();
      }
    }
    
    coef.get_mpq(entry);
    RatMat_set(constraints, constraints->nR-2, i, entry);
    (-coef).get_mpq(entry);
    RatMat_set(constraints, constraints->nR-1, i, entry);
  }
  
  parallelVal[1].get_mpq(entry);
  RatMat_set(constraints, constraints->nR-2, constraints->nC-1, entry);  //this is x <= b in  a<=x<=b
  
  (-parallelVal[0]).get_mpq(entry);
  RatMat_set(constraints, constraints->nR-1, constraints->nC-1, entry); //this is the -a <= -x
  equalityType.push_back(-1);
  equalityType.push_back(-1);
  
  //Now the constraints matrix and equality type are all set up -- we just need
  //to do the linear programming
  vector<rational> solutionVector(constraints->nC-1, rational());
  rational scl;
  
  if (VERBOSE==1)
    cout << "About to do linear programming\n";
  
  linear_program_from_ratmat( polygon_list,solutionVector, scl, constraints, equalityType, solver, VERBOSE ) ;
  //cout << "About to reset mem functions\n";
  //mp_set_memory_functions(0, 0, 0);


  //Now, the answers *should* be decimal versions of exact rational numbers
  //so we should be able to convert them with rat approx:
  //rational ratScl = approxRat(scl);
  
  if (VERBOSE==1)
    cout << "Done -- got scl = " << scl.get_d() << " = " << scl << "\n";
  
  vector<rational> coefs(2, rational(0,1)); //this is the 2d point under delta
  
  if (scl == rational(1,1)) { //if scl=1, return the first point
    coefs[0] = points[0][0];
    coefs[1] = points[0][1];
  } else {
    //vector<rational> ratSolnVec(constraints->nC-1, rational());
    //for (i=0; i<constraints->nC-1; i++) {
    //  ratSolnVec[i] = approxRat(solutionVector[i]);
    //}
    
    //we must apply the delta function to get the point in g,h space
    for (i=0; i<(int)polygon_list.size(); i++) {  //for each poly
      for (j=0; j<polygon_list[i].size; j++) {  //go through the edges
        firstWordNumber = 0;                    //this is the word index of the first word in the chain
        for (k=0; k<numChains; k++) {
          if (arc_list[polygon_list[i].arc[j]].first_word == firstWordNumber &&
              arc_list[polygon_list[i].arc[j]].first == 0) {                        //if we're looking at the first word in a chain, first letter
            coefs[k] = coefs[k] + solutionVector[i]/weights[firstWordNumber];  //add this many copies of the chain
          }
          firstWordNumber+=chains[k].size();
        }
      }
    }

    //scale it so that the returned point has scl = 1
    for (i=0; i<2; i++) {
      coefs[i] = coefs[i] / scl;
    }
  }
  
  if (VERBOSE==1)
    cout << "new point with scl = 1 at: " << coefs[0] << ", " << coefs[1] << "\n";
  
  
  //remove the three new rows that we added
  RatMat_change_num_rows(constraints, constraints->nR-3);
  equalityType.resize(constraints->nR);
  
  return coefs;
}














//find the scl of a particular point 
rational point_scl(vector<vector<string> >&chains, 
                   vector<arc>& arc_list,
                   vector<polygon>& polygon_list,
                   vector<int>& weights,
                   RatMat* constraints,
                   vector<int>& equalityType,
                   vector<rational>& point, 
                   scallop_lp_solver solver, 
                   int VERBOSE) {
  int i,j,k;
  int firstWordNumber;
  int numChains = (int)chains.size();
  vector<rational> coef(2, rational());
  //note this is 2 new constraints; we will add these rows to constraints, and
  //get rid of them later
  RatMat_change_num_rows(constraints, constraints->nR+2);
  
  //we just have to check that delta_g = point[0] and delta_h=point[1]
  for (i=0; i<(int)polygon_list.size(); i++) {  //for each poly
 
    coef[0] = rational(0,1);                      //this will be the coefficient
    coef[1] = rational(0,1);
    
    for (j=0; j<polygon_list[i].size; j++) {  //go through the edges
      firstWordNumber = 0;                    //this is the word index of the first word in the chain
      for (k=0; k<numChains; k++) {
        if (arc_list[polygon_list[i].arc[j]].first_word == firstWordNumber &&
            arc_list[polygon_list[i].arc[j]].first == 0) {                        //if we're looking at the first word in a chain, first letter
          coef[k] = coef[k] + (rational(1,1)/weights[firstWordNumber]);
        }
        firstWordNumber+=chains[k].size();
      }
    }
    
    RatMat_set_int(constraints, constraints->nR-2, i, coef[0].n(), coef[0].d());
    RatMat_set_int(constraints, constraints->nR-1, i, coef[1].n(), coef[1].d());
  }
  RatMat_set_int(constraints, constraints->nR-2, constraints->nC-1, point[0].n(), point[0].d());
  RatMat_set_int(constraints, constraints->nR-1, constraints->nC-1, point[1].n(), point[1].d());
  equalityType.push_back(0);
  equalityType.push_back(0);
  
   //Now the constraints matrix and equality type are all set up -- we just need
  //to do the linear programming
  vector<rational> solutionVector(constraints->nC-1, rational());
  rational scl;
  linear_program_from_ratmat( polygon_list,solutionVector, scl, constraints, equalityType, solver, VERBOSE ) ;
      //cout << "About to reset mem functions\n";
    //mp_set_memory_functions(0, 0, 0);
  //Now, the answers *should* be decimal versions of exact rational numbers
  //so we should be able to convert them with rat approx:
  //rational ratScl = approxRat(scl);
  
  RatMat_change_num_rows(constraints, constraints->nR-2);
  equalityType.resize(constraints->nR);
  
  return scl;
}






// this function creates the basic constraints (arc constraints, chains)
void create_constraint_matrix(vector<vector<string> >& chains, vector<arc>& arc_list,
                                               vector<polygon>& polygon_list,
                                               vector<int>& weights,
                                               RatMat* constraints,
                                               vector<int>& equalityType) {
  
  int arc_list_length = arc_list.size();
  int polygon_list_length = polygon_list.size();
  int numChains = chains.size();
  int totalNumWords = 0;
  int offset;  // as we build the matrix, we use this for sanity
  int i,j,k,l;
  int n;
  int myWordNumber, firstWordNumber;
  //vector<int> equalityType; //-1 is <=, 0 is ==
  int coefFirstWord;
  int coefMyWord;
  
  for (i=0; i<numChains; i++) {
    totalNumWords += chains[i].size();
  }
  
  
  /*
    build the matrix of constraints for the unit ball:
    each column is a polygon, saying how much of it we want
    rows come from:
     - every arc, saying that that arc must appear the same number of times
       as positive and negative
     - the chains, saying that the boundary decomposes as ag + bh
     The last column is the b vector, as in Ax <= b is the polyhedron
  */
  RatMat_init(constraints,  arc_list_length/2 + totalNumWords - numChains, 
                            polygon_list_length + 1);
  equalityType.resize(arc_list_length/2 + totalNumWords - numChains);
  
  RatMat_reset_int(constraints, 0);
  
  /*  this section removed since I believe glpk will do it for us in lp.cc
  //each polygon must appear >= 0 times:
  for (i=0; i<polygon_list_length; i++) {
    RatMat_set_int(constraints, i,i, -1,1);
    RatMat_set_int(constraints, i, polygon_list_length, 0,1);
    equalityType[i] = -1;
  }
  */
  
  offset = 0;
  
  //each arc must appear the same number of times coming and going
  for (i=0; i<arc_list_length/2; i++) {  //for each arc
    for (j=0; j<polygon_list_length; j++) {  //for each polygon
      //count the number of times the arc appears in the polygon, with signs
      n=0;
      for (k=0; k<polygon_list[j].size; k++) {  // for each arc in the polygon
        if(polygon_list[j].arc[k]==2*i){
           n=n+1;
        }
        if(polygon_list[j].arc[k]==2*i+1){  //the same arc, but backwards
          n=n-1;
        }
      }
      // n is now the coefficient in the offset+i row and jth column
      RatMat_set_int(constraints, offset+i, j, n, 1);
    }
    equalityType[offset+i] = 0;
    // note that these rows are == 0, which the last entry in the jth row already is
  }
  
  offset += arc_list_length/2;
  //the offset is now polygon_list_length + arc_list_length/2
  
  //each chain must appear as one -- that if we have aw + bv, they must appear
  //as caw + cbv = c(aw+bv) (if the weights are 1, then just the same # of times)
  
  firstWordNumber = 0; //note this is a running total through all words
  for (i=0; i<numChains; i++) { //for each chain
    for (j=1; j<(int)chains[i].size(); j++) {  //for each word past the first one
      myWordNumber = firstWordNumber + j;
      for (k=0; k<polygon_list_length; k++) { //for each polygon
        //the number of times our word appears, divided by its weight is the same
        //as the number of times that the first word appears, divided by its weight
        coefFirstWord = 0;
        coefMyWord = 0;
        for (l=0; l<polygon_list[k].size; l++) {
          if (arc_list[polygon_list[k].arc[l]].first_word == myWordNumber &&
              arc_list[polygon_list[k].arc[l]].first == 0 ) { //for each copy of the word, this is 1
            coefMyWord ++;
          }
          else if (arc_list[polygon_list[k].arc[l]].first_word == firstWordNumber &&
                   arc_list[polygon_list[k].arc[l]].first == 0) {
            coefFirstWord++;
          }
        }
        //the entry to put in is mycoef/myweight - firstcoef/firstweight
        RatMat_set_int(constraints, offset, k, 
               weights[firstWordNumber]*coefMyWord - weights[myWordNumber]*coefFirstWord,
                              weights[myWordNumber]*weights[firstWordNumber]);
      }
      equalityType[offset] = 0;
      //again, the right hand side is 0, so we don't need to set it
      offset++;
    }
    firstWordNumber += chains[i].size();
  }  

}





vector<vector<rational> > ball_in_positive_quadrant(vector<vector<string> >& chains,
                                                     vector<string> wordList,
                                                     vector<int> weights,
                                                     int maxjun, scallop_lp_solver solver,
                                                     int VERBOSE) {
  int i,j,k;
  //initialize lrs:
  //scalbal_lrs_init();
  
  vector<arc> arc_list(0);
  vector<polygon> polygon_list(0);
  
  //generate the arcs!
  generate_arcs(&arc_list, wordList, wordList.size());
  int arc_list_length = arc_list.size();
  
  if(VERBOSE==1){
    cout << "generated arcs\n";
    if (VERBOSE >1) {
      cout << arc_list_length << " arcs (start letter, end letter, start word, end word) \n";
      for(i=0;i<arc_list_length;i++){
	      cout << "arc " << i << " : ";
	      cout << arc_list[i].first << " " << arc_list[i].last << " " << arc_list[i].first_word << " " << arc_list[i].last_word << "\n";
  
      }
    }
  };
  
  //generate the polygons!
  generate_polygons(wordList, polygon_list, arc_list, maxjun); 
  int polygon_list_length = polygon_list.size();
  
  if(VERBOSE==1){
    cout << "generated polygons\n";
    if (VERBOSE>1) {
      cout << polygon_list_length << " polygons (cyclic list of arcs) \n";
      for(i=0;i<polygon_list_length;i++){
	      cout << "polygon " << i << " : ";
	      for(j=0;j<polygon_list[i].size;j++){
	        cout << polygon_list[i].arc[j] << " ";
	      };
	      cout << '\n';
      }
    }
  }
  
  
  
  //create the constraint matrix
  RatMat* constraints = new RatMat[1];
  vector<int> equalityType;
  
  create_constraint_matrix(chains, arc_list, polygon_list, weights, constraints,
                                                                 equalityType);
  if (VERBOSE==1) {
    //cout << "Constraint Matrix:\n";
    //RatMat_print(constraints, 1);       
  }                                           
  
  //pick the two basic points (1,0) and (0,1)
  vector< vector<rational> > points;
  points.resize(2);
  for (i=0; i<2; i++) {
    points[i].resize(2);
    points[i][0] = rational(1-i,1);
    points[i][1] = rational(i,1);
  }
  
  //find their scls and scale
  rational scl;
  for (i=0; i<2; i++) {
    scl = point_scl(chains, arc_list, polygon_list, weights, constraints, equalityType, points[i], solver, VERBOSE);
    if (VERBOSE) {
      cout << "Found scl of chain "<< i << ": " << scl << "\n";
    }
    for (j=0; j<2; j++) {
      points[i][j] = points[i][j]/scl;
    }
  }
  
  //thus we now have two points on the scl ball; go into the main loop now
  vector<rational> newPoint(2, rational());
  vector< vector<rational> > pointPair;
  pointPair.resize(2);
  for (i=0; i<2; i++) {
    pointPair[i].resize(2);
  }  
  i=0;
  while (i < (int)points.size()-1) { //note that points.size() will change!
    //take the next two points and find a scl vertex between them
    for (j=0; j<2; j++) {
      for (k=0; k<2; k++) {
        pointPair[j][k] = points[i+j][k];
      }
    }
    newPoint = min_point_on_line(chains, arc_list, polygon_list, weights,
                                                         constraints,
                                                         equalityType,
                                                         pointPair, solver, VERBOSE);
    if (newPoint[0] == points[i][0] && newPoint[1] == points[i][1]) { //if scl is linear
      if (VERBOSE==1) {
        cout << "scl appears to be linear between these points\n";
      }
      i++; //don't add any new points
    } else {
      points.insert(points.begin()+i+1, newPoint); // this causes the point to be inserted after i
      if (VERBOSE==1) {
        cout << "I added a new point: " << newPoint[0] << ", " << newPoint[1] << "\n";
        cout << "Now there are " << points.size() << " points, and i = " << i << "\n"; 
      }
    }
  }
  
  delete[] constraints;
  
  //now points should hold all the points in the polyhedron quadrant
  return points;
}
                                               




void draw_ball(string fileName, vector<vector<string> >& chains, 
                                vector<int>& weights,
                                vector<vector<rational> >& points) {
  int i,j;
  fstream outfile;
  outfile.open(fileName.c_str(), fstream::out);
  
  outfile << "\%!PS-Adobe-2.0 EPSF-2.0\n";
  outfile << "\%\%BoundingBox: 0 0 1024 1024\n\n";      
  outfile << "/Arial findfont\n";
  outfile << "20 scalefont setfont\n";
  outfile << "1 setlinejoin\n";
  outfile << "2 setlinewidth\n";
  outfile << "512 512 translate\n";
  
  //scale so that the farthest point is 3/4 of the way out
  // (in the L^\infty norm)
  rational maxPoint;
  for (i=0; i<(int)points.size(); i++) {
    for (j=0; j<(int)points[i].size(); j++) {
      if (maxPoint < points[i][j]) {
        maxPoint = points[i][j];
      }
    }
  }
  //512 is the amount to the edge, so "1" should be equal to
  //(3/4)*512/maxpoint
  rational scale = rational(5,8)*(rational(512,1)/maxPoint);
  
  
  //draw the grid
  outfile << "1 setlinewidth\n";
  outfile << "0.75 setgray\n";
  double x,y;
  x = scale.get_d()*floor(-1.2*maxPoint.get_d());
  while (x < scale.get_d()*ceil(1.2*maxPoint.get_d())) {
    if (fabs(x) < 0.0001) {
      outfile << "0 setgray\n";
      outfile << x << " -512 moveto\n";
      outfile << x << " 512 lineto\nstroke\n";
      outfile << "0.75 setgray\n";
    } else {
      outfile << x << " -512 moveto\n";
      outfile << x << " 512 lineto\nstroke\n";
    }
    x += scale.get_d();   
  }
  y = scale.get_d()*floor(-1.2*maxPoint.get_d());
  while (y < scale.get_d()*ceil(1.2*maxPoint.get_d())) {
    if (fabs(y) < 0.0001) {
      outfile << "0 setgray\n";
      outfile << "-512 " << y << " moveto\n";
      outfile << "512 " << y << " lineto\nstroke\n";
      outfile << "0.75 setgray\n";
    } else {
      outfile << "-512 " << y << " moveto\n";
      outfile << "512 " << y << " lineto\nstroke\n";
    }
    y+= scale.get_d();      
  } 
  
  outfile << "0 setgray\n";
  outfile << "2 setlinewidth\n";
  
  //draw the polygon
  vector<double> coords(2, 0.0);
  coords[0] = (scale*points[0][0]).get_d();
  coords[1] = (scale*points[0][1]).get_d();
  outfile << coords[0] << " " << coords[1] << " moveto\n";
  for (i=1; i<(int)points.size(); i++) {
    coords[0] = (scale*points[i][0]).get_d();
    coords[1] = (scale*points[i][1]).get_d();
    outfile << coords[0] << " " << coords[1] << " lineto\n";
  }
  outfile << "closepath\nstroke\n";
  
  //and the dots, plus their labels
  double pointLen;
  double pointScaler;
  double xLength, yLength; //these are for the axes labels
  outfile << "5 setlinewidth\n";
  for (i=0; i<(int)points.size(); i++) {
    if (points[i][1] == rational(0,1) && rational(0,1) < points[i][0]) {
      xLength = points[i][0].get_d();
    } else if (points[i][0] == rational(0,1) && rational(0,1) < points[i][1]) {
      yLength = points[i][1].get_d();
    }
    coords[0] = (scale*points[i][0]).get_d();
    coords[1] = (scale*points[i][1]).get_d();
    outfile << coords[0] << " " << coords[1] << " moveto\n";
    outfile << coords[0] << " " << coords[1] << " 2  0 360 arc\n"; 
    outfile << "stroke\n";
    //label
    pointLen = sqrt( (points[i][0]*points[i][0]).get_d() + (points[i][1]*points[i][1]).get_d() );
    pointScaler = (pointLen+0.3)/pointLen;
    outfile << pointScaler*coords[0] << " " << pointScaler*coords[1] << " moveto\n";
    outfile << "((" << points[i][0] << ", " << points[i][1] << ")) dup stringwidth pop 2 div neg 0 rmoveto show\n";
  }
  outfile << "2 setlinewidth\n";
  
  //now the chain labels
  stringstream tempNum;
  tempNum.str("");
  tempNum << weights[0];
  string outString = tempNum.str() + chains[0][0];
  for (i=1; i<(int)chains[0].size(); i++) {
    tempNum.str("");
    tempNum << weights[i];
    outString += " + " + tempNum.str() + chains[0][i];
  }
  outfile << (scale.get_d()*(xLength+0.6)) << " 0 moveto\n";
  outfile << "(" << outString << ") show\n";
  
  tempNum.str("");
  tempNum << weights[chains[0].size()];
  outString = tempNum.str() + chains[1][0];
  for (i=1; i<(int)chains[1].size(); i++) {
    tempNum.str("");
    tempNum << weights[chains[0].size() + i];
    outString += " + " + tempNum.str() + chains[1][i];
  }
  outfile << "0 " << (scale.get_d()*(yLength+0.5)) << " moveto\n";  
  outfile << "(" << outString << ") dup stringwidth pop 2 div neg 0 rmoveto show\n";

  outfile << "\%eof\n";
  
  outfile.close();
}


/*****************************************************************************/
/* create the unit ball in 2 dimensions                                      */
/*****************************************************************************/
void create_print_unit_ball(vector<vector<string> >& chains,
                            vector<int>& weights,
                            string fileName,
                            int maxjun, 
                            scallop_lp_solver solver,
                            int VERBOSE) {
  int i,j;
  
  //first we do it with the regular words
  vector<string> wordList(0);
  for (i=0; i<(int)chains.size(); i++) {
    for (j=0; j<(int)chains[i].size(); j++) {
      wordList.push_back(chains[i][j]);
    }
  } 
  vector< vector<rational> > points1 = ball_in_positive_quadrant(chains, wordList, weights, maxjun, solver, VERBOSE);
  
  //then we swap the second chain
  wordList.resize(0);
  for (i=0; i<(int)chains.size(); i++) {
    for (j=0; j<(int)chains[i].size(); j++) {
      if (i==1) {
        wordList.push_back(inverse(chains[i][j]));
        chains[i][j] = inverse(chains[i][j]);
      } else {
        wordList.push_back(chains[i][j]);
      }
    }
  } 
  vector< vector<rational> > points2 = ball_in_positive_quadrant(chains, wordList, weights, maxjun, solver, VERBOSE);  

  for (i=0; i<(int)chains.size(); i++) {
    for (j=0; j<(int)chains[i].size(); j++) {
      chains[i][j] = inverse(chains[i][j]);
    }
  }

  //collect all the points
  vector< vector<rational> > allPoints;
  allPoints.resize(2*(points1.size() + points2.size()));
  
  for (i=0; i<(int)points1.size(); i++) {
    allPoints[i].resize(2);
    for (j=0; j<2; j++) {
      allPoints[i][j] = points1[i][j];
    }
  }
  for (i=0; i<(int)points2.size(); i++) {  //we do this in reverse order so that they go ccw
    allPoints[points1.size()+i].resize(2);
    allPoints[points1.size()+i][0] = -points2[points2.size()-i-1][0]; //thus we get the left quadrant
    allPoints[points1.size()+i][1] = points2[points2.size()-i-1][1];
  }
  
  //now we reflect all these points through the origin (take negatives)
  for (i=0; i<(int)points1.size()+(int)points2.size(); i++) {
    allPoints[points1.size() + points2.size() + i].resize(2);
    allPoints[points1.size() + points2.size() + i][0] = -allPoints[i][0];
    allPoints[points1.size() + points2.size() + i][1] = -allPoints[i][1];
  }
  
  //now allPoints has all the points on the unit ball
  //for (i=0; i<(int)allPoints.size(); i++) {
  //  cout << "(" <<allPoints[i][0] << ", " << allPoints[i][1] << ")\n";
  //}
  
  //remove the duplicates
  i=0;
  while (i<(int)allPoints.size()-1) {
    if (allPoints[i][0] == allPoints[i+1][0] && 
        allPoints[i][1] == allPoints[i+1][1]) {
      allPoints.erase(allPoints.begin()+i);
      continue;
    }
    i++;
  }
  //remove the last one, since it should be a dup
  allPoints.erase(allPoints.begin()+allPoints.size()-1);
  
  //now allPoints has all the points on the unit ball
  cout << "Unit ball: {";
  for (i=0; i<(int)allPoints.size(); i++) {
    cout << "(" <<allPoints[i][0] << ", " << allPoints[i][1] << "), ";
  }
  cout << "}\n";
  
  //draw the ball
  draw_ball(fileName, chains, weights, allPoints);
  
  
}

/*****************************************************************************/
/* create the unit ball in 3 dimensions!!                                    */
/*****************************************************************************/
void create_print_unit_ball_3D(vector<vector<string> >& chains,
                              vector<int>& weights,
                              string fileName,
                              int maxjun, 
                              scallop_lp_solver solver,
                              int VERBOSE) {
  //in 3 dimensions, we can't just list all the points in order, so we will 
  //produce a list a points (vertices), plus a list of triples of indices, which
  //will be the triangles in a triangulation of the polyhedron
  
  //we need to do all four top orthants separately; then we'll reflect them
  //the chains to take the inverses will be the first two chains, and which 
  //to invert will be given by the binary digits of "orthant"
  //1 = make inverse
  
  vector<vector<vector<rational> > > orthantVertices(4); //these are the actual vertices
  vector<vector<vector<int> > > orthantTriangles(4); //these are the triangles
  
  //it's orthantTriangles[i][j][k], i=which orthant, j=which triangle, k=index in triangle
  
  int orthant;
  int i,j;
  
  for (orthant=0; orthant<4; orthant++) {
    //invert as described above
    for (i=0; i<2; i++) {
      if ( (orthant>>i)&1 == 1 ) {
        //invert this chain
        for (j=0; j<chains[i].size(); j++) {
          chains[i][j] = inverse(chains[i][j]);
        }
      }
    }
    //compute the orthant
    ball_in_positive_orthant(chains, weights, maxjun, solver, VERBOSE,
                             orthantVertices[orthant],
                             orthantTriangles[orthant]);
    //now fix the orthant vertices, and invert the chains back
    //note we need to take the negative of the vertex coordinates which are
    //highlighted 1 in "orthant"
    for (i=0; i<2; i++) {
      if ( (orthant>>i)&1 == 0) {
        continue;
      }
      //invert the chain
      for (j=0; j<chains[i].size(); j++) {
        chains[i][j] = inverse(chains[i][j]);
      }
      //invert the vertices
      for (j=0; j<orthantVertices[orthant].size(); j++) {
        orthantVertices[orthant][j][i] *= -1;
      }
    }
  }
  
  
  //*******
  //now build a big list of all the vertices and triangles -- note that
  //we need to be careful about reindexing the vertices (change the labels 
  //in the triangles)
  
  vector<vector<rational> > allVertices(0);
  vector<vector<int> > allTriangles(0);
  int currentOffset = 0;
  
  
  
  
}

int main(int argc, char* argv[]){


  int i,j,k;
  int VERBOSE = 0;
  int overrideMaxjun = 0;
  int overriddenMaxjun = -1;
  
  scallop_lp_solver solver = EXLP;
  
  init_lp();
  
  string drawFile = "";

  if(argc < 2 || strcmp(argv[1],"-h")==0){
    cout << "usage: scabble [-v] [-glpk] [-mn] fileName [i_1]w_1 [i_1]w_2 . . [i_n]w_n , [j_1]v_1 [j_1]v_2 . . [j_m]v_m  \n";
    cout << "\t\tproduces the unit ball in the plane spanned by the chains, output to fileName\n";
    cout << "-mn overrides the max number of sides for the polygons\n";
    return(0);
  };
  
  while (argv[1][0] == '-') {
    switch (argv[1][1]) {
      case 'v':
        VERBOSE=1;
        break;
      case 'g':
        solver = GLPK_DOUBLE;
        break;
      case 'm':
        overrideMaxjun = 1;
        overriddenMaxjun = atoi(&(argv[1][2]));
        break;
      default:
        cout << "Unrecognized option\n";
        exit(1);
        break;
    }
    for (i=1; i<argc-1; i++) {
      argv[i] = argv[i+1];
    }
    argc--;
  }
  


  //NOTE: it would be really awesome if we knew how big to make the polygon
  //list, since right now it probably does tons of single-element reallocs, 
  //which is just horrible

  vector<vector<string> > chains;
  vector<string> wordList;
  vector<int> weights;
  //int arc_list_length;
  //int polygon_list_length;
  //vector<arc> arc_list(0);
  //vector<polygon> polygon_list(0);
  vector<char> lettersSeen(0);
  int maxjun;
  char c;

  int chainsStart = 2;

  //cout << "starting...\n";
  //fflush(stdout);

  //probably unnecessary
  chains.resize(1);
  chains[0].resize(0);


  //preprocess the words to remove the weights,
  //and add them to the chain lists
  int weightLen;
  int whichChain = 0;
  int wordInd;
  for(i=chainsStart;i<argc;i++){    
    
    if (VERBOSE == 1) {
      cout << "word: " << argv[i] << "\n";
    }
  
    if (string(argv[i]) == ",") {
      whichChain++;
      chains.resize(chains.size()+1);
      chains[whichChain].resize(0);
      continue;
    }
    
    wordList.push_back(argv[i]);
    wordInd = (int)wordList.size()-1;
    
    weightLen = 0;
    j=0;
    while (isdigit(wordList[wordInd][j])) {
      j++;
    }
    if (j>0) {
      weights.push_back(atoi(wordList[wordInd].substr(0,j).c_str()));
    } else {
      weights.push_back(1);
    }
    wordList[wordInd] = wordList[wordInd].substr(j,wordList[wordInd].length()-j);
    
    //make sure the word is cyclically reduced
    cyc_red(wordList[wordInd]);
    
    if (VERBOSE == 1) {
      cout << "cyclically reduced word: " << wordList[wordInd] << "\n";
    }
    
    chains[whichChain].push_back(wordList[wordInd]);
  }
  
  drawFile = string(argv[1]);
  
  //find the rank (count the letters)
  for (i=0; i<(int)wordList.size(); i++) {
    for (j=0; j<(int)wordList[i].size(); j++) {
      c = tolower(wordList[i][j]);
      for (k=0; k<(int)lettersSeen.size(); k++) { 
        if (lettersSeen[k] == c) {
          break;
        }
      }
      if (k==(int)lettersSeen.size()) { //did we not break?
        lettersSeen.push_back(c);
      }
    }
  }
  
  if (overrideMaxjun == 0) {
    maxjun = 2*lettersSeen.size(); //i.e. maxjun is twice the rank
  } else {
    maxjun = overriddenMaxjun;
    cout << "Overridding polygon sides to " << maxjun << "\n";
  }
  
  
  if (VERBOSE == 1) {
    cout << "Number of chains: " << chains.size() << "\n";
    
    cout << "Rank: " << lettersSeen.size() << "\n";
  
    cout << "chains:\n";
    for (i=0; i< (int)chains.size(); i++) {
      for (j=0; j<(int)chains[i].size(); j++) {
        cout << chains[i][j] << " ";
      }
      cout << "\n";
    }
    
    cout << "Word list: ";
    for (i=0; i<(int)wordList.size(); i++) {
      cout << wordList[i] << " ";
    }
    cout << "\n";
    
    cout << "weights: \n";
    for (i=0; i<(int)weights.size(); i++) {
      cout << weights[i] << " ";
    }
    cout << "\n";
  }
       
  
  
  //create, draw the unit ball
  if (chains.size() == 2) {
    create_print_unit_ball(chains, weights, drawFile, maxjun, solver, VERBOSE);
  } else if (chains.size() == 3) {
    create_print_unit_ball_3D(chains, weights, drawFile, maxjun, solver, VERBOSE);
  } else {
    cout << "You have too many chains or something\n";
  }

  
  return 0;
}
