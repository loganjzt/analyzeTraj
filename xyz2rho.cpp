// Modified on Oct 11 for 3 dimension
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <cstring>
#include <fstream>
#include <ctime>
#include <new>
#include <malloc.h>
#include </usr/include/xdrfile/xdrfile_trr.h> // xdr include file 
#include </usr/include/xdrfile/xdrfile_xtc.h> // xdr include file 

    using namespace std;

int read_xyz(FILE *f, double **coor);
double getZs(double **coor,double** box,double binSize,int natoms,int natomPmol,char xyz);

int main(int argc, char **argv){
	
	int natomPmol = 1;
	double binSize = 0.1;

	//other variables
	
	char f[40];
	int i,j,k;
	
	//these are average values	
	int xbin,ybin,zbin;
	double zs,xs,ys;

	double **box = new double *[3];
	box[0] = new double [3];
	box[1] = new double [3];
	box[2] = new double [3];
	
	box[0][0] = 10.0;
	box[1][1] = 10.0;
	box[2][2] = 10.0;

	xbin = int(box[0][0]/binSize);
	ybin = int(box[1][1]/binSize);
	zbin = int(box[2][2]/binSize);

	FILE *data = fopen("output.xyz","r");
	int natoms;
	fscanf(data,"%d",&natoms);
	printf("natoms = %d\n",natoms);
	rewind(data);

	double **coor = new double* [natoms];
	for(int i = 0 ; i < natoms ; ++i ) 
		coor[i] = new double [3];
	if(coor == NULL){
		printf("error in memory allocation\n");
		return 1;
	}


	printf("# of atoms = %d\n# Xbin = %d\n# Ybin = %d\n# Zbin = %d\n",natoms,xbin,ybin,zbin);

	double *rhoX = new double [xbin];
	double *rhoY = new double [ybin];
	double *rhoZ = new double [zbin];

	for(i=0;i<xbin;i++) rhoX[i] = 0.0;
	for(i=0;i<ybin;i++) rhoY[i] = 0.0;
	for(i=0;i<zbin;i++) rhoZ[i] = 0.0;


	//open xtc file and loop through each frame
	k = 0;

	while( ! read_xyz(data,coor )){
		if(coor == 0){
			printf("Insufficient memory to load .trr file. \n");
			return 1;
		}

		// X direction
		xs = getZs(coor,box,binSize,natoms,natomPmol,'x');
		for(i=0;i<natoms/natomPmol;i++){
			if(coor[i*natomPmol][0] - xs >= 0) rhoX[int((coor[i*natomPmol][0] - xs )/binSize)] += 1.0;
			else rhoX[int((coor[i*natomPmol][0] - xs + box[0][0])/binSize)] += 1.0;
		}

		// Y direction
		ys = getZs(coor,box,binSize,natoms,natomPmol,'y');
		for(i=0;i<natoms/natomPmol;i++){
			if(coor[i*natomPmol][1] - ys >= 0) rhoX[int((coor[i*natomPmol][1] - ys )/binSize)] += 1.0;
			else rhoX[int((coor[i*natomPmol][1] - xs + box[1][1])/binSize)] += 1.0;
		}

		// Z direction
		zs = getZs(coor,box,binSize,natoms,natomPmol,'z');
		for(i=0;i<natoms/natomPmol;i++){
			if(coor[i*natomPmol][2] - zs >= 0) rhoX[int((coor[i*natomPmol][2] - zs )/binSize)] += 1.0;
			else rhoX[int((coor[i*natomPmol][2] - zs + box[2][2])/binSize)] += 1.0;
		}

		k++;
	
	}

	for(i=0;i<xbin;i++) rhoX[i] = rhoX[i]/double(k)/box[2][2]/box[1][1]/binSize;
	for(i=0;i<ybin;i++) rhoY[i] = rhoY[i]/double(k)/box[0][0]/box[2][2]/binSize;
	for(i=0;i<zbin;i++) rhoZ[i] = rhoZ[i]/double(k)/box[0][0]/box[1][1]/binSize;
	
	printf("# Finish reading .trr file...\n");
	printf("# Counted frame = %d...\n",k);

	double nCheck;	// integral of rhoZ to check the area under the curve	
	double rhoMin,rhoMax;	// get the max and min value of rho

	rhoMin = rhoX[0];
	rhoMax = rhoX[0];
	sprintf(f,"rho_x.dat");
	data = fopen(f,"w");
	fprintf(data,"# x\tRhox\tint Rhox\n");

	nCheck = 0.0;
	for(i=0;i<xbin;i++){
		if(i<xbin-1) nCheck += (rhoX[i] + rhoX[i+1])/2.0*double(binSize);
		fprintf(data,"%.3f\t%e\t%e\n",i*binSize,rhoX[i],nCheck);
		if(rhoMin > rhoX[i]) rhoMin = rhoX[i];
		if(rhoMax < rhoX[i] ) rhoMax = rhoX[i];
	}
	fprintf(data,"# rhoMax = %f\n# rhoMin = %f\n",rhoMax,rhoMin);

	rhoMin = rhoY[0];
	rhoMax = rhoY[0];
	sprintf(f,"rho_y.dat");
	data = fopen(f,"w");
	fprintf(data,"# y\tRhoy\tint Rhoy\n");

	nCheck = 0.0;
	for(i=0;i<ybin;i++){
		if(i<ybin-1) nCheck += (rhoY[i] + rhoY[i+1])/2.0*double(binSize);
		fprintf(data,"%.3f\t%e\t%e\n",i*binSize,rhoY[i],nCheck);
		if(rhoMin > rhoY[i]) rhoMin = rhoY[i];
		if(rhoMax < rhoY[i] ) rhoMax = rhoY[i];
	}
	fprintf(data,"# rhoMax = %f\n# rhoMin = %f\n",rhoMax,rhoMin);

	sprintf(f,"rho_z.dat");
	data = fopen(f,"w");
	fprintf(data,"# z\tRhoz\tint Rhoz\n");

	nCheck = 0.0;
	rhoMin = rhoZ[0];
	rhoMax = rhoZ[0];
	for(i=0;i<zbin;i++){
		if(i<zbin-1) nCheck += (rhoZ[i] + rhoZ[i+1])/2.0*double(binSize);
		fprintf(data,"%.3f\t%e\t%e\n",i*binSize,rhoZ[i],nCheck);
		if(rhoMin > rhoZ[i]) rhoMin = rhoZ[i];
		if(rhoMax < rhoZ[i] ) rhoMax = rhoZ[i];
	}
	fprintf(data,"# rhoMax = %f\n# rhoMin = %f\n",rhoMax,rhoMin);

	/*
	//	get the average max and min
 	//average from max to 0.99 max is the average of max
 	//average from min to 1.01 min is the average of min
 
	double rhoMaxAverage = 0.0;
	double rhoMinAverage = 0.0;
	int count1 = 0;
	int count2 = 0;

	for(i=0;i<zbin;i++){
		if(rhoZ[i] >= 0.99* rhoMax){
			rhoMaxAverage += rhoZ[i];
			count1++;
		}
		if(rhoZ[i] <= 1.01 * rhoMin){
			rhoMinAverage += rhoZ[i];
			count2++;
		}
	}

	fprintf(data,"# <rhoMax> = %f\n# <rhoMin> = %f\n",rhoMaxAverage/double(count1),rhoMinAverage/double(count2));
	fprintf(data,"# average from max - 0.99 max, min - 1.01 min\n");
	*/

	return 0;
}

/* func for reading xyz files */

int read_xyz(FILE *f, double **coor){
	int i,length = 0;
	char label = 'C';
	char buffer[100];
	fscanf(f,"%d",&length);
	//printf("natoms = %d\n",length);

	if(feof(f)) return 1; // end of file
	
	label = fgetc(f);
	
	label = fgetc(f);
	
	if(label != '\n'){
		fseek(f,-1,SEEK_CUR);
		fgets(buffer,100,f);
	}
	for(i=0;i<length;i++){
		fscanf(f, "%c %lf %lf %lf\n", &label, &coor[i][0], &coor[i][1], &coor[i][2]);
		//printf("%f %f %f\n", coor[i][0], coor[i][1], coor[i][2]);
	}
	return 0;
}

/* test read_xyz */
/*
int main(){
	int natoms = 0;
	FILE *f = fopen("output.xyz","r");
	fscanf(f,"%d",&natoms);
	printf("natoms = %d\n",natoms);
	rewind(f);
	
	double **coor = new double* [natoms];
	for(int i = 0 ; i < natoms ; ++i ) 
		coor[i] = new double [3];
	if(coor == NULL){
		printf("error in memory allocation\n");
		return 1;
	}
	printf("No error in memory allocation\n");

	while(!read_xyz(f,coor)){
		for(int i=0;i<natoms;i++) printf("%d\t%f\n",i+1,coor[i][0]);
	}

	delete []coor;

	return 0;
}
*/

double getZs(double **coor,double** box,double binSize,int natoms,int natomPmol,char xyz){

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
	if(abs(comZ - box[alpha][alpha]/2.0) > binSize ){
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
	else zs = 0.0;
	return 0.0;
}


