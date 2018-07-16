// calculate the pressure tensor
// read trr or xtc files
// output Pxx(z), Pyy(z) and Pzz(z)
//
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <cstring>
#include <fstream>
#include <ctime>
#include <malloc.h>
#include </usr/include/xdrfile/xdrfile_trr.h> // xdr include file 
#include </usr/include/xdrfile/xdrfile_xtc.h> // xdr include file 

    using namespace std;

 // Old get Zs 
 /*
double getZs(rvec *coor,matrix box,double binSize,int natoms,int natomPmol,char xyz){

	int alpha;
	if(xyz == 'x')  alpha = 0;
	else if(xyz == 'y') alpha = 1;
	else if(xyz == 'z') alpha = 2;
	else{
		printf("getZs input error!\n");
		return 1.0;
	}
	int i,j;
	double zs;
	int zbin = (box[alpha][alpha]/binSize+1);

	double comZ = 0.0;
	for(i=0;i<natoms/natomPmol;i++) comZ += coor[i*natomPmol][alpha]/double(natoms/natomPmol);

	// calculate zs
	if(abs(comZ - box[alpha][alpha]/2.0) > binSize*0.5 ){
		double indicator1,indicator2;
		double integral;

		// get rho'(z)
		double *rhoZtmp = new double [zbin];
		for(i=0;i<zbin;i++) rhoZtmp[i] = 0.0;

		for(i=0;i<natoms/natomPmol;i++){
			rhoZtmp[int(coor[i*natomPmol][alpha]/binSize)] += 1.0/double(natoms/natomPmol)/binSize;
		}
	
		indicator2 = 0.0;	
		for(i=1;i< ( ( zbin-1) *10.0 );i++){
			zs = i*binSize/10.0;
			integral = 0.0;
			for(j = 0 ; (j+1)*binSize < zs ; j++)  integral += (rhoZtmp[j]+rhoZtmp[j+1])/2.0*binSize;
			integral += rhoZtmp[j]*(zs-j*binSize) ;
			indicator1 = integral - zs / box[alpha][alpha] + comZ / box[alpha][alpha] - 0.50;
			if(  indicator2*indicator1 <= 0 && (zs-box[alpha][alpha]/2.0)*(comZ - box[alpha][alpha]/2.0) <= 0 && rhoZtmp[int(zs/binSize)] < 1.2/double(zbin)/binSize  ){
				zs = binSize/10.0 / (abs(indicator1)+abs(indicator2))*abs(indicator2) + (i-1)*binSize/10.0;
				delete [] rhoZtmp;
				return zs;
			}
			indicator2 = indicator1;
		}

		if(indicator2*indicator1 > 0 ){
			printf("# Error when remove c.o.m...%f.(%f)..\n",comZ,box[alpha][alpha]/2.0);
			FILE *data;
			data = fopen("getZS.log","w");
			for(int zi=0;zi<zbin;zi++) fprintf(data,"%f\t%f\n",zi*binSize,rhoZtmp[zi]);

			indicator2 = 0.0;	
			for(int i = 1;i< ( zbin*10.0 );i++){
				zs = i*binSize/10.0;
				integral = 0.0;
				for(j = 0 ; (j+1)*binSize < zs ; j++) integral += (rhoZtmp[j]+rhoZtmp[j+1])/2.0*binSize;
				integral += rhoZtmp[j]*(zs-j*binSize) ;
				indicator1 = integral - zs / box[alpha][alpha] + comZ / box[alpha][alpha] - 0.50;
				zs = binSize/10.0 / (abs(indicator1)+abs(indicator2))*abs(indicator2) + (i-1)*binSize/10.0;
				fprintf(data,"%f\t%f\t%f\n",zs,indicator1,indicator2);
				indicator2 = indicator1;	
			}
			zs = box[alpha][alpha]+1.0;
			return zs;
		}
	}
	else{
		zs = 0.0;
		return zs;
	}
}
*/

double getZs(rvec *coor,matrix box,double binSize,int natoms,int natomPmol,char xyz){

	int alpha;
	double tmp1;
	double dz;
	int i,j;
	double zs = 0.0;

	if(xyz == 'x')  alpha = 0;
	else if(xyz == 'y') alpha = 1;
	else if(xyz == 'z') alpha = 2;
	else{
		printf("getZs input error!\n");
		return 1.0;
	}

	double comZ = 0.0;
	for(i=0;i<natoms/natomPmol;i++){
		comZ += coor[i*natomPmol][alpha]/double(natoms/natomPmol);
	}

	// calculate zs
	if( comZ - box[alpha][alpha]/2.0 > binSize ){
		dz = comZ - box[alpha][alpha]/2.0;
		for(i=1;i< box[alpha][alpha]/2.0/binSize*10.0 ; i++){
			zs = i*binSize/10.0 + dz; 
			tmp1 = 0.0;
			for(j=0;j<natoms/natomPmol;j++){
				if(coor[j*natomPmol][alpha] >= zs) tmp1 += ( coor[j*natomPmol][alpha] - zs ) /double(natoms/natomPmol);
				else tmp1 += (coor[j*natomPmol][alpha] - zs + box[alpha][alpha])/double(natoms/natomPmol);
			}
			printf("%f\t%f\n",zs,tmp1);
			if(tmp1 <= box[alpha][alpha]/2.0) return zs;
		}
		printf("%f\t%f\t%f\t%f\n",comZ,tmp1,box[2][2]/2.0,zs);
		printf("# Error 102\n");
	}
	else if( comZ - box[alpha][alpha]/2.0 < -1.0*binSize ){
		dz = comZ - box[alpha][alpha]/2.0 + box[alpha][alpha];

		for(i=1;i< 10.0*box[alpha][alpha]/2.0/binSize;i++){
			zs = i*binSize/10.0+dz;
			tmp1 = 0.0;
			for(j=0;j<natoms/natomPmol;j++){
				if(coor[j*natomPmol][alpha] >= zs) tmp1 += (coor[j*natomPmol][alpha]-zs)/double(natoms/natomPmol);
				else tmp1 += (coor[j*natomPmol][alpha] - zs + box[alpha][alpha])/double(natoms/natomPmol);
			}
			if(tmp1 >= box[alpha][alpha]/2.0) return zs;
		}
		printf("%f\t%f\t%f\t%f\n",comZ,tmp1,box[2][2]/2.0,zs);
		printf("# Error 102\n");

	}
	else zs = 0.0;

	return zs;
}


int main(int argc, char **argv){
	
	// XTC variables
	XDRFILE *xd;		// the xtc file

	int natoms;	// number of total atoms
	int step;		// the step counter for the simulation
	float time;		// simulation time

	matrix box;		// box coordinates in a 3x3 matrix
	rvec *coor;		// atom coordinates in a 2-D matrix
	rvec *vel;		// atom coordinates in a 2-D matrix
	rvec *force;		// atom coordinates in a 2-D matrix
//	float prec;	
	float lambda;
	

	int firstframe = 500;
	if(argc >= 3) firstframe = atoi(argv[2]);
	printf("# first frame = %d \n",firstframe);

	int lastframe = 2000;			// ps
	if(argc >= 4) lastframe = atoi(argv[3]);
	printf("# last frame = %d \n",lastframe);

	int natomPmol = 1;
	if(argc >= 5) natomPmol = atoi(argv[4]);
	printf("# atoms per molecules= %d \n",natomPmol);

	double binSize = 0.1;
	if(argc >= 6) binSize = atof(argv[5]);
	printf("# binSize = %f \n",binSize);

	double temperature = 148.0;
	if(argc >= 7) temperature = atof(argv[6]);
	printf("# temperature = %f \n",temperature);

	double rCut = 1.0;
	if(argc >= 8) rCut = atof(argv[7]);
	printf("# cut-off = %f \n",rCut);

	printf("#--------------end of input parametars#\n");

	//other variables
	
	double zs;

//	double epsilon = 1.23047;
//	double sigma = 0.373;
//	double beta = 1000.0/8.314/temperature;
	int nframe = lastframe*10+1;
	int zbin;	// divide box into zbin slabs.
	
	double vSlab;
	int nCount = 0;	// count number of frames

	double **pzz;
	double **pxx;
	double **pyy;
		pzz = (double **)malloc(sizeof(double *)*nframe);	
		pxx = (double **)malloc(sizeof(double *)*nframe);	
		pyy = (double **)malloc(sizeof(double *)*nframe);	

//	double pzzTime,pxxTime,pyyTime;

	//char f[40];
//	FILE *data;
	int i,j;

	double comZ,comZold;

	//read xtc files
	read_trr_natoms(argv[1],&natoms);

	coor = (rvec *)malloc(natoms*sizeof(rvec));
	vel = (rvec *)malloc(natoms*sizeof(rvec));
	force = (rvec *)malloc(natoms*sizeof(rvec));

	//open xtc file and loop through each frame
	xd=xdrfile_open(argv[1],"r");
	while( ! read_trr(xd, natoms, &step, &time, &lambda, box, coor, vel, force)){
		if(nCount > nframe) printf("# not enough memory \n");
		if(coor == 0){
			printf("Insufficient memory to load .trr file. \n");
			return 1;
		}
		if(step == 0){
			zbin = int(box[2][2]/binSize);
			vSlab = binSize*box[0][0]*box[1][1];
			for(i=0;i<nframe;i++){
				pzz[i] = (double *)malloc(sizeof(double)*zbin);
				pxx[i] = (double *)malloc(sizeof(double)*zbin);
				pyy[i] = (double *)malloc(sizeof(double)*zbin);
				for(j=0;j<zbin;j++){
					pzz[i][j] = 0.0;
					pxx[i][j] = 0.0;
					pyy[i][j] = 0.0;
				}
			}
			printf("# zbin = %d\n",zbin);
		}

	    if(time >= firstframe && ( time <= lastframe || lastframe < firstframe ) ){
			
			zs = getZs(coor,box,binSize,natoms,natomPmol,'z');		// get comZ in each configuration
			if(zs > box[2][2] || zs < 0.0 ) return 1;
			comZ = 0.0;
			for(i=0;i<natoms;i++){
				if(coor[i][2] >= zs) comZ += (coor[i][2]-zs)/double(natoms);
				else comZ += (coor[i][2]-zs+box[2][2])/double(natoms);
			}
		//	printf("%f\t%f\t%f\n",time,comZ,zs);
			if(abs(comZ-box[2][2]/2.0)> binSize){
				printf("#error!\n");
				comZold = 0.0;
				for(i=0;i<natoms;i++) comZold += coor[i][2]/double(natoms);

				printf("%f\t%f\t%f\t%f\t%f\n",time,comZold,comZ,box[2][2]/2.0,zs);

				double *rhoZtmp = new double [zbin];
				for(i=0;i<zbin;i++) rhoZtmp[i] = 0.0;

				for(i=0;i<natoms/natomPmol;i++)
					rhoZtmp[int(coor[i*natomPmol][2]/binSize)] += 1.0/binSize;
				for(i=0;i<zbin;i++) printf("%f\t%f\n",i*binSize,rhoZtmp[i]);
				return 1;
			}
			nCount++;
    	}
	}
	
	printf("# number of frame counted = %d\n",nCount);	

	return 0;
}

