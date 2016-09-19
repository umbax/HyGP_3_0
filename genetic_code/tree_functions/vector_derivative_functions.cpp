// Copyright 2016 Dr Umberto Armani
//  
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//  
//      http://www.apache.org/licenses/LICENSE-2.0
//  
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.


// vector operations

// definition of int_rand function (not recognized by C++)
int int_rand(int max) 
{
	// remember: initialize the random number generator!!
	// the returned number is within the range [0,max-1]
	if (max == 0) {
		cerr << "\n\n vector_derivative_functions.cpp/int_rand : ERROR! max == 0!" << endl;
		exit (-1);
	}

	int v = rand() % max;
	return v;
}

void basic_stat_analysis(double* af,double* out,int n_el)
{
	int COMMENT = 0;

	int n = n_el;
	double ndbl = (double)n;

	if (COMMENT)
		cout << "\nbasic_stat_analysis";

	// min and max
	out[2] = *max_element(af, af+n);
	out[0] = *min_element(af, af+n);

	//mean
	out[1] = 0.;
	for (int i=0; i<n; i++)
		out[1] += af[i];
	out[1]=out[1]/ndbl;

	//unbiased estimate of the population variance
	out[3] = 0.;
	for (int i=0; i<n; i++)
		out[3]+= pow(af[i]-out[1],2.0);
	out[3]=out[3]/(ndbl-1.);

	if (COMMENT) {
		cout << "\nMax :  = "<<out[2];
		cout << "\nMin : = "<<out[0];
		cout << "\nMean = " << out[1];
		cout << "\nVariance = " << out[3];
	}

}

// function that compute the vectorial subtraction x-y
// INPUT:
// vector x, vector y, vector dimension, address of array where to store result components
// OUTPUT
// no output (components of vector x - y are stored in s[])
void vect_sub(double* x, double* y, int dim, double* s)
{
	for (int k=0; k<dim; k++)
		s[k] =  x[k] - y[k];
}

// function that compute the vectorial summation x+y
// INPUT:
// vector x, vector y, vector dimension, address of array where to store result components
// OUTPUT
// no output (components of vector x - y are stored in s[])
void vect_sum(double* x, double* y, int dim, double* s)
{
	for (int k=0; k<dim; k++)
		s[k] =  x[k] + y[k];
}


// function that compute the division between a vector and a scalar
// INPUT:
// vector x, scalar l, vector dimension, address of array where to store result components
// OUTPUT
// no output (components of vector x/l are stored in s[])
void vect_div_scal(double* x, double l, int dim, double* s)
{
	for (int k=0; k<dim; k++)
		s[k] =  x[k]/l;
}


// function that compute the multiplication between a vector and a scalar
// INPUT:
// vector x, scalar l, vector dimension, address of array where to store result components
// OUTPUT
// no output (components of vector x*l are stored in s[])
void vect_mult_scal(double* x, double l, int dim, double* s)
{
	for (int k=0; k<dim; k++)
		s[k] =  x[k]*l;
}


double point_dist(double* x, double* y, int dim)
{
	int COMMENT = 0;
	double dist =0.;

	// dot product <x-y,x-y>
	for (int k=0; k<dim; k++) 
		// define distance (here euclidean norm)
		dist = dist + (x[k]-y[k])*(x[k]-y[k]);
	
	dist = sqrt(dist);
	if (COMMENT) {
		cout << "\nPoint 1: ( " << x[0] << " , " << x[1] << " , " << x[2] << ")" << endl;
		cout << "Point 2: ( " << y[0] << " , " << y[1] << " , " << y[2] << ")" << endl;
		cout << "Distance = " << dist << endl;
	}
	
	return dist;
}

// function to compute the magnitude of a vector
double vect_magn(double *x, int dim)
{
	double magn;
	for (int k=0; k<dim; k++) 
		// define distance (here euclidean norm)
		magn = magn + (x[k])*(x[k]);
	
	magn = sqrt(magn);
	
	return magn;
}



// utility function for sorting the population
int point_comp(const void *a, const void *b)
{
	// recast the arguments
	const double* el1 = *(const double **)a;
	const double* el2 = *(const double **)b;	
	
	int ind = 0;
	//cout << "\nel1[" << ind << "] = " << el1[ind];
	//cout << "\nel2[" << ind << "] = " << el2[ind] << endl;

	// do the comparison
   if (el1[ind]>el2[ind])
		return 1;
	if (el1[ind]==el2[ind])
		return 0;
	if (el1[ind]<el2[ind])
		return -1;
}

// function that returns the row numbers of a desired number of closest neighbours to a selected point of the set
// INPUT:
// address of the array (matrix) containing the points
// number of points 
// number of dimensions (no of variables) of each point
// address of the selected point (now it has to belong to the original set of points)
// size of the patch (no of points the patch is made of - includes the selected point)
// address of the matrix p_neighbours
// OUTPUT:
// no output (distance ([][0]) and no of row ([][1]) of the neighbours (the number of rows is patch_size) is stored in p_neighbours)
void get_neighbours(double** points, int n_points, int dim, double *central_p, int patch_size, double** p_neighbours)
{
	// dynamic allocation of dist, array that will store index and distance of each point
	double **dist = new double* [n_points];
	for 	(int k=0; k<n_points; k++)
		dist[k] = new double [2];

	for (int k=0; k<n_points; k++)  {
		// compute distances of each point from central point (stored in dist[k][0])
		dist[k][0] = point_dist(central_p, &(points[k][0]), dim);
		// copy the corresponding number of the point (no of row in the matrix) 
		dist[k][1] = k;
	}
	// sort array "dist", not points!!
	qsort(dist,n_points, sizeof(dist[n_points-1]) ,point_comp);    //(dim+1)*sizeof(double)

	// update neighbours
	for (int k=0; k<patch_size; k++) {
		p_neighbours[k][0] = dist[k][0];
		p_neighbours[k][1] = dist[k][1];
	}
	
	// free memory: dist won't be used again
	for (int k=0; k< n_points; k++)
		delete dist[k];
	delete [] dist;	
}

// function to compute finite differences
// INPUT:
// value of the function at the initial point
// value of the function at the final point
// magnitude of the vector from initial to final point 
// OUTPUT:
// estimate of the derivative
double fin_differences(double f_initial, double f_final, double ds)
{
	double result;
	//cout << "\nf_initial = " << f_initial;
	//cout << "\nf_final = " << f_final;
	result = (f_final-f_initial)/ds;
	return result;
}


// function that find the closest neighbour and return the value of the estimated derivative using input data
// INPUT:
// address of the matrix containing the data (Z1 Z2  ... ZN TARGET  - row i)
// OUTPUT:
// value of the derivative 
void get_deriv_to_closest_neighbour(double **points, int n_points, int nvar, double** fun_der, int n_der)
{
	// matrix fun_der (nx3):
	// column 0: value of the incremental ratio in initial point using initial and final points (it's just a rough approximation of the derivative!)
	// column 1: row number of selected or initial point (where the derivative is computed)
	// column 2: row number of the final point in sample matrix
	for (int row=0; row<n_der; row++) {

		// select a central point (among the points in the list)
		int n_sel_point = row; //int_rand(n_points);   //=k : this way the validating point are selected
		double *central_p = &(points[n_sel_point][0]);	
		if (nvar==2) cout << "\nPoint selected : row " << n_sel_point << " => (" << central_p[0] << ", " << central_p[1] << ")"  << endl;  
		// update fun_der with the coordinates of the selected point
		fun_der[row][1] = n_sel_point;
		
		// select the size of the patch (no of points)
		int patch_size = 2;  //selected point included
		double** p_neighbours = new double*[patch_size];
		for 	(int k=0; k<patch_size; k++)
			p_neighbours[k] = new double [2];

		// call the function "find_neighbours" to reordinate the point array according to distance
		get_neighbours(points, n_points, nvar, central_p, patch_size, p_neighbours);
	
		// update third entry of fun_der: no of row of the final point
		fun_der[row][2] = p_neighbours[1][1];


		double** neighbours = new double*[patch_size];
		for 	(int k=0; k<patch_size; k++)
			neighbours[k] = new double [nvar];

		// update neighbours
		int i;
		for (int k=0; k<patch_size; k++) {
			i = (int)(p_neighbours[k][1]);
			for (int j=0; j<nvar+1; j++) {
				neighbours[k][j] = points[i][j]; 
			}
		}

		// show results
		cout << "\nPatch size " << patch_size << " . Neighbours, , including the selected point:" << endl;
		for (int k=0; k<patch_size; k++) {	
			cout << "\nPoint " << k << " : (";
			for (int j=0; j<nvar; j++) {
				cout << neighbours[k][j] << " , ";
			}
			cout << ") output : " <<  neighbours[k][nvar];
			cout << "\ndistance " <<p_neighbours[k][0];
			cout << " row number: " << p_neighbours[k][1] << endl;		
		}

		// define the difference vector connecting the two points in the patch
		double* s = new double [nvar]; 
		double s_magn;

		// subtract vectors
		vect_sub(&(neighbours[1][0]),&(neighbours[0][0]), nvar, s); 
		cout << "Vector s = ( " << s[0] << ", " << s[1] << ")" << endl;

		// magnitude of difference vector s
		s_magn = point_dist(&(neighbours[1][0]),&(neighbours[0][0]),nvar);
		cout << "Vector s magnitude = " << s_magn << endl;

	
		// estimate derivative
		fun_der[row][0] = fin_differences(neighbours[0][nvar],neighbours[1][nvar], s_magn);
		cout << "\nEstimated derivative = " << fun_der[row][0] << endl;

		// free memory taken by neighbours_p and neighbours
		for (int k=0; k<patch_size; k++) {
			delete [] p_neighbours[k];
			delete [] neighbours[k];
		}
		delete [] p_neighbours;	
		delete [] neighbours;
//// end cycle
	}
}




