/*
*  R : A Computer Language for Statistical Data Analysis
*  Copyright (C) 2004   The R Development Core Team.
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with this program; if not, a copy is available at
*  http://www.r-project.org/Licenses/
*/

#include <R.h>
//#include "modreg.h" /* for declarations for registration */
/*
x - data matrix
pn - # data points
pp - # dimensions
cen - centers (array)
pk - k
cl - # data points??
pmaxiter - max iters
nc - number centers
wss - 
*/
void kmeans_Rvq(double *x, int *pn, int *pp, double *cen, int *pk, int *cl,
int *pmaxiter, int *nc, double *wss)
{
	int n = *pn, k = *pk, p = *pp, maxiter = *pmaxiter;
	int iter, i, j, c, it, inew = 0;
	double best, dd, tmp;
	Rboolean updated;
	//new variables
	int loc, **indexsorted, *isused, *ixLoc, count;
	double  **distances;
	

	indexsorted = (int **) malloc(k*sizeof(int *));
	isused = (int *) malloc(n*sizeof(int));
	ixLoc = (int *) malloc(k*sizeof(int));
	distances = (double **) malloc(k*sizeof(double *));

	for(i=0;i<k;i++){
		distances[i] = (double *)malloc(n * sizeof(double));
		indexsorted[i] = (int *)malloc(n * sizeof(int));
	}

	for(i = 0; i < n; i++) cl[i] = -1;
	for(iter = 0; iter < maxiter; iter++) {
		updated = FALSE;
		//find distance to all points from each center
		for(i = 0; i < n; i++) {
			for(j = 0; j < k; j++) {
				dd = 0.0;
				for(c = 0; c < p; c++) {
					tmp = x[i+n*c] - cen[j+k*c];
					dd += tmp * tmp;
				}
				distances[j][i] = dd;
			}
			isused[i] = 0;			
		}
		/* for each center, sort by distance to point */
		for(i = 0; i < k; i++) {
			for(j = 0; j < n; j++) {
				indexsorted[i][j] = j;

			}
			rsort_with_index(distances[i],indexsorted[i], n);

		}
		for(j = 0; j<k ; j++) {ixLoc[j] = 0;}
		count = 0;
		while(count < n){
			//printf("i: %u\n", i);
			for(j = 0; j < k; j++) {
				if(count >=100) break;
				//	printf("ixloc: %u\n", ixLoc[j]);
				loc = indexsorted[j][ixLoc[j]];
				//	printf("loc: %u\n", loc);

				while(isused[loc] == 1){
					ixLoc[j] += 1;
					loc = indexsorted[j][ixLoc[j]];
				}
				if(cl[loc] != j+1){
					cl[loc] = j+1;
					updated = TRUE;
				}
				isused[loc] = 1;
				count++;
			}
		}
		if(!updated) break;
		/* update each centre */
		for(j = 0; j < k*p; j++) cen[j] = 0.0;
		for(j = 0; j < k; j++) nc[j] = 0;
		for(i = 0; i < n; i++) {
			it = cl[i] - 1; 
			nc[it]++;
			for(c = 0; c < p; c++) cen[it+c*k] += x[i+c*n];
		}
		for(j = 0; j < k*p; j++) cen[j] /= nc[j % k];
	}

	// set each item to it's actual closest point
	for(i = 0; i < n; i++) {
		/* find nearest centre for each point */
		best = R_PosInf;
		for(j = 0; j < k; j++) {
			dd = 0.0;
			for(c = 0; c < p; c++) {
				tmp = x[i+n*c] - cen[j+k*c];
				dd += tmp * tmp;
			}
			if(dd < best) {
				best = dd;
				inew = j+1;
			}
		}
		cl[i] = inew;

	}
	for(i=0;i<k;i++){
		free(distances[i]);
		free(indexsorted[i]);
	}
	free(indexsorted);
	free(isused);
	free(ixLoc);
	free(distances);
	*pmaxiter = iter + 1;
	for(j = 0; j < k; j++) wss[j] = 0.0;
	for(i = 0; i < n; i++) {
		it = cl[i] - 1;
		for(c = 0; c < p; c++) {
			tmp = x[i+n*c] - cen[it+k*c];
			wss[it] += tmp * tmp;
		}
	}
}