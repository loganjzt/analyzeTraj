// Modified on Oct 11 for 3 dimension Density
// Feb 27, 2018
// 		besides rho_x,y,z , report rho(r) assume there is cylinder symmetry in z direction
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <cstring>
#include <fstream>
#include <ctime>
#include <malloc.h>
#include </usr/include/xdrfile/xdrfile_trr.h> // xdr include file 
#include </usr/include/xdrfile/xdrfile_xtc.h> // xdr include file 
#include "/home/zhitongj/Git/xtc2rho/getZs.cpp"


    using namespace std;


int main(int argc, char **argv){
	
	// input variables
	int firstframe = 500;				
	if(argc >= 3) firstframe = atoi(argv[2]);
	int lastframe = -1;		// ps
	if(argc >= 4) lastframe = atoi(argv[3]);
	int natomPmol = 1;
	if(argc >= 5) natomPmol = atoi(argv[4]);
	double binSize = 0.1;
	if(argc >= 6) binSize = atof(argv[5]);

	printf("# first frame = %d \n",firstframe);
	printf("# last frame = %d \n",lastframe);
	printf("# atoms per molecules= %d \n",natomPmol);

	// XTC variables
	XDRFILE *xd;		// the xtc file

	int natoms;	// number of total atoms
	int step;		// the step counter for the simulation
	float time;		// simulation time

	matrix box;		// box coordinates in a 3x3 matrix
	rvec *coor;		// atom coordinates in a 2-D matrix
	rvec *vel;		// atom coordinates in a 2-D matrix
	rvec *force;		// atom coordinates in a 2-D matrix
	float lambda;
//	float prec;	
	
	//other variables
	int i,j,k;
	int nframe = 0;
	
	// output related
	char f[40];
	FILE *data;

	// density related
	int xbin,ybin,zbin;
	int xi,yi,zi;
	
	int rbin;
	double dx,dy,dz,dr;

	double zs,xs,ys;

	double ***rho;			//rho[xbin][ybin][zbin]
	double *rhor;

	/*read xtc files*/
	read_trr_natoms(argv[1],&natoms);

	coor = (rvec *)malloc(natoms*sizeof(rvec));
	vel = (rvec *)malloc(natoms*sizeof(rvec));
	force = (rvec *)malloc(natoms*sizeof(rvec));

	//open xtc file and loop through each frame
	xd=xdrfile_open(argv[1],"r");
	nframe = 0;
	while( ! read_trr(xd, natoms, &step, &time, &lambda, box, coor, vel, force)){
		if(coor == 0){
			printf("Insufficient memory to load .trr file. \n");
			return 1;
		}
		if(step == 0){
			xbin = int(box[0][0]/binSize);
			ybin = int(box[1][1]/binSize);
			zbin = int(box[2][2]/binSize);
			rbin = min(xbin/2,ybin/2);

			printf("# of atoms = %d\n# Xbin = %d\n# Ybin = %d\n# Zbin = %d\n",natoms,xbin,ybin,zbin);

			rho = new double** [xbin];
			rhor = new double [rbin];
			for(i=0;i<xbin;i++){
				rho[i] = new double* [ybin];
				for(j=0;j<ybin;j++) rho[i][j] = new double [zbin];
			}
			for(i=0;i<xbin;i++){
				for(j=0;j<ybin;j++){
					for(k=0;k<zbin;k++) rho[i][j][k] = 0.0;
				}
			}
			for(i=0;i<rbin;i++) rhor[i] = 0.0;
			printf("# natoms / grid = %.2f\n",double(natoms)/double(zbin)/double(xbin)/double(ybin));

		}

	    if(time >= firstframe && ( time <= lastframe || lastframe < firstframe ) ){
				
			if(int(time*10)%1000 == 0) printf("time = %.2f\n",time);

			xs = getZs(coor,box,binSize,natoms,natomPmol,'x');	
			ys = getZs(coor,box,binSize,natoms,natomPmol,'y');
			zs = getZs(coor,box,binSize,natoms,natomPmol,'z');

			for(i=0;i<natoms/natomPmol;i++){
				if(coor[i*natomPmol][0] - xs >= 0) xi = int((coor[i*natomPmol][0] - xs )/binSize);
				else xi = int((coor[i*natomPmol][0] - xs + box[0][0])/binSize);
				if(xi == xbin) xi = 0;
			
				if(coor[i*natomPmol][1] - ys >= 0) yi = int((coor[i*natomPmol][1] - ys )/binSize);
				else yi = int((coor[i*natomPmol][1] - ys + box[1][1])/binSize);
				if(yi == ybin) yi = 0;

				if(coor[i*natomPmol][2] - zs >= 0) zi = int((coor[i*natomPmol][2] - zs )/binSize);
				else zi = int((coor[i*natomPmol][2] - zs + box[2][2])/binSize);
				if(zi == zbin) zi = 0;

				rho[xi][yi][zi] += 1.0;

				// rho r
				if(coor[i*natomPmol][0] - xs >= 0) dx = coor[i*natomPmol][0] - xs - box[0][0]/2.0;
				else dx = coor[i*natomPmol][0] - xs + box[0][0] - box[0][0]/2.0;

				if(coor[i*natomPmol][1] - ys >= 0) dy = coor[i*natomPmol][1] - ys - box[1][1]/2.0;
				else dy = coor[i*natomPmol][1] - ys + box[1][1] - box[1][1]/2.0;
	
				if(coor[i*natomPmol][2] - zs >= 0) dz = coor[i*natomPmol][2] - zs - box[2][2]/2.0;
				else dz = coor[i*natomPmol][2] - zs + box[2][2] - box[2][2]/2.0;

				if(abs(dx) > box[0][0]/2.0 || abs(dy) > box[1][1]/2.0){
					printf("# Error\n");
					return 1;
				}
				dr = sqrt(dx*dx+dy*dy);
				
				if(dr < box[0][0]/2 && dr < box[1][1]/2 )
					rhor[int(dr/binSize)] += 1.0;
			}
			nframe++;
	    }
	}

	printf("# Finish reading .trr file...\n");
	printf("# Counted frame = %d...\n",nframe);

	for(i=0;i<xbin;i++){
		for(j=0;j<ybin;j++){
			for(k=0;k<zbin;k++) rho[i][j][k] = rho[i][j][k] / double(nframe) / binSize / binSize / binSize;
		}
	}
	for(i=0;i<rbin;i++){
		rhor[i] = rhor[i] /(M_PI * (i+1)*binSize* (i+1)*binSize * box[2][2] - M_PI *i*binSize * i*binSize *box[2][2]) / double(nframe) ;
	}

	double rhoMin,rhoMax;	// get the max and min value of rho

	rhoMin = rho[0][0][0];
	rhoMax = rho[0][0][0];
	sprintf(f,"rho_3d.dat");
	data = fopen(f,"w");
	fprintf(data,"#x \t y \t  z \t Rho3d\n");

	for(i=0;i<xbin;i++){
		for(j=0;j<ybin;j++){
			for(k=0;k<zbin;k++){
				fprintf(data,"%f\t%f\t%f\t%e\n",i*binSize,j*binSize,k*binSize,rho[i][j][k]);
				if(rhoMin > rho[i][j][k]) rhoMin = rho[i][j][k];
				if(rhoMax < rho[i][j][k]) rhoMax = rho[i][j][k];
			}
		}
	}

	fprintf(data,"# rhoMax = %f\n# rhoMin = %f\n",rhoMax,rhoMin);

	//	get the average max and min

	sprintf(f,"rhor_2d.dat");
	data = fopen(f,"w");
	fprintf(data,"#r \t Rhor \n");

	rhoMax = 0.0; 
	rhoMin = 0.0;

	for(i=0;i<rbin;i++){
		fprintf(data,"%f\t%e\n",i*binSize,rhor[i]);
		if(rhoMin > rhor[i]) rhoMin = rhor[i];
		if(rhoMax < rhor[i]) rhoMax = rhor[i];
	}

	fprintf(data,"# rhoMax = %f\n# rhoMin = %f\n",rhoMax,rhoMin);

	fclose(data);

	return 0;
}


