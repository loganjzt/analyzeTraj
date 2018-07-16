/* Oct 11, 2017 
 * modified for 3 dimensions (averaged over the rest two directions) 
 * -- zhitong
 */

/* Jul 14, 2018
 * Check input format and work for both .xtc and .trr files
 * use new instead of malloc
 * -- zhitong
 */ 		

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <cstring>
#include <fstream>
#include <ctime>
#include </usr/include/xdrfile/xdrfile_trr.h> // xdr include file 
#include </usr/include/xdrfile/xdrfile_xtc.h> // xdr include file 

#include "/home/zhitongj/Git/traj2rho/getZs.cpp"	// move the c.o.m of the system to the center of the box


    using namespace std;


int main(int argc, char **argv){
	
	int ifTrr;		// index for input files type. 1 = .trr, 0 = .xtc
	
	if( argv[1][int(sizeof(argv[1])/sizeof(char)) - 1] == 'r' ){
		ifTrr = 1;
		printf("# trr file\n");
	}
	else if( argv[1][int(sizeof(argv[1])/sizeof(char)) - 1] == 'c' ){
		ifTrr = 0;
		printf("# xtc file\n");
	}
	else{
		printf("wrong input files type\n");
		return 0;
	}

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

	/*read xtc files*/
	if(ifTrr == 1) read_trr_natoms(argv[1],&natoms);
	else read_xtc_natoms(argv[1],&natoms);

	rvec *vel;		// atom coordinates in a 2-D matrix
	rvec *force;		// atom coordinates in a 2-D matrix
	float lambda;
	vel = new rvec [natoms];
	force = new rvec [natoms];

	float prec;	
	coor = new rvec [natoms];

	//other variables
	int i,k;
	
	// output related
	char f[40];
	FILE *data;

	// density related
	int xbin,ybin,zbin;
	double zs,xs,ys;
	double *rhoZ,*rhoX,*rhoY;

	//open xtc file and loop through each frame
	xd=xdrfile_open(argv[1],"r");
	k = 0;

	if(ifTrr == 1){ // input is trr file
		while( ! read_trr(xd, natoms, &step, &time, &lambda, box, coor, vel, force)){
			if(coor == 0){
				printf("Insufficient memory to load .trr file. \n");
				return 1;
			}
			if(step == 0){
				xbin = int(box[0][0]/binSize);
				ybin = int(box[1][1]/binSize);
				zbin = int(box[2][2]/binSize);
	
				printf("# of atoms = %d\n# Xbin = %d\n# Ybin = %d\n# Zbin = %d\n",natoms,xbin,ybin,zbin);
	
				rhoX = new double [xbin];
				rhoY = new double [ybin];
				rhoZ = new double [zbin];
	
				for(i=0;i<xbin;i++) rhoX[i] = 0.0;
				for(i=0;i<ybin;i++) rhoY[i] = 0.0;
				for(i=0;i<zbin;i++) rhoZ[i] = 0.0;
				printf("# natoms / zbin = %d\n",natoms/zbin);
	
			}
	
		    if(time >= firstframe && ( time <= lastframe || lastframe < firstframe ) ){
					
				if(int(time*10)%5000 == 0) printf("time = %.2f\n",time);		// print out time each 0.5ns
	
				// X direction
				xs = getZs(coor,box,binSize,natoms,natomPmol,'x');	
				for(i=0;i<natoms/natomPmol;i++){
					if(coor[i*natomPmol][0] - xs >= 0) rhoX[int((coor[i*natomPmol][0] - xs )/binSize)] += 1.0;
					else rhoX[int((coor[i*natomPmol][0] - xs + box[0][0])/binSize)] += 1.0;
				}
	
				// Y direction
				ys = getZs(coor,box,binSize,natoms,natomPmol,'y');
				for(i=0;i<natoms/natomPmol;i++){
					if(coor[i*natomPmol][1] - ys >= 0) rhoY[int((coor[i*natomPmol][1] - ys )/binSize)] += 1.0;
					else rhoY[int((coor[i*natomPmol][1] - ys + box[1][1])/binSize)] += 1.0;
				}
	
				// Z direction
				zs = getZs(coor,box,binSize,natoms,natomPmol,'z');		// get comZ in each configuration
				for(i=0;i<natoms/natomPmol;i++){
					if(coor[i*natomPmol][2] - zs >= 0) rhoZ[int((coor[i*natomPmol][2] - zs )/binSize)] += 1.0;
					else rhoZ[int((coor[i*natomPmol][2] - zs + box[2][2])/binSize)] += 1.0;
				}
	
				k++;
		    }
		}
	}	
	else { // input is xtc file
		while(! read_xtc(xd, natoms, &step, &time, box, coor, &prec) ){
			if(coor == 0){
				printf("Insufficient memory to load .xtc file. \n");
				return 1;
			}
			if(step == 0){
				xbin = int(box[0][0]/binSize);
				ybin = int(box[1][1]/binSize);
				zbin = int(box[2][2]/binSize);
	
				printf("# of atoms = %d\n# Xbin = %d\n# Ybin = %d\n# Zbin = %d\n",natoms,xbin,ybin,zbin);
	
				rhoX = new double [xbin];
				rhoY = new double [ybin];
				rhoZ = new double [zbin];
	
				for(i=0;i<xbin;i++) rhoX[i] = 0.0;
				for(i=0;i<ybin;i++) rhoY[i] = 0.0;
				for(i=0;i<zbin;i++) rhoZ[i] = 0.0;
				printf("# natoms / zbin = %d\n",natoms/zbin);
	
			}
	
		    if(time >= firstframe && ( time <= lastframe || lastframe < firstframe ) ){
					
				if(int(time*10)%5000 == 0) printf("time = %.2f\n",time);		// print out time each 0.5ns
	
				// X direction
				xs = getZs(coor,box,binSize,natoms,natomPmol,'x');	
				for(i=0;i<natoms/natomPmol;i++){
					if(coor[i*natomPmol][0] - xs >= 0) rhoX[int((coor[i*natomPmol][0] - xs )/binSize)] += 1.0;
					else rhoX[int((coor[i*natomPmol][0] - xs + box[0][0])/binSize)] += 1.0;
				}
	
				// Y direction
				ys = getZs(coor,box,binSize,natoms,natomPmol,'y');
				for(i=0;i<natoms/natomPmol;i++){
					if(coor[i*natomPmol][1] - ys >= 0) rhoY[int((coor[i*natomPmol][1] - ys )/binSize)] += 1.0;
					else rhoY[int((coor[i*natomPmol][1] - ys + box[1][1])/binSize)] += 1.0;
				}
	
				// Z direction
				zs = getZs(coor,box,binSize,natoms,natomPmol,'z');		// get comZ in each configuration
				for(i=0;i<natoms/natomPmol;i++){
					if(coor[i*natomPmol][2] - zs >= 0) rhoZ[int((coor[i*natomPmol][2] - zs )/binSize)] += 1.0;
					else rhoZ[int((coor[i*natomPmol][2] - zs + box[2][2])/binSize)] += 1.0;
				}
	
				k++;
		    }
		}
	}

	printf("# Finish reading .trr file...\n");
	printf("# Counted frame = %d...\n",k);

	for(i=0;i<xbin;i++) rhoX[i] = rhoX[i]/double(k)/box[2][2]/box[1][1]/binSize;
	for(i=0;i<ybin;i++) rhoY[i] = rhoY[i]/double(k)/box[0][0]/box[2][2]/binSize;
	for(i=0;i<zbin;i++) rhoZ[i] = rhoZ[i]/double(k)/box[0][0]/box[1][1]/binSize;
	

	double nCheck;			// integral of rhoZ to check the area under the curve	
	double rhoMin,rhoMax;	// get the max and min value of rho

	// output rho (x)
	sprintf(f,"rho_x.dat");
	data = fopen(f,"w");
	fprintf(data,"# x\tRhox\tint Rhox\t");
	fprintf(data,"(# Counted frame = %d...)\n",k);

	nCheck = 0.0;
	rhoMin = rhoX[0];
	rhoMax = rhoX[0];
	for(i=0;i<xbin;i++){
		if(i<xbin-1) nCheck += (rhoX[i] + rhoX[i+1])/2.0*double(binSize);
		fprintf(data,"%.3f\t%e\t%e\n",i*binSize,rhoX[i],nCheck);
		if(rhoMin > rhoX[i]) rhoMin = rhoX[i];
		if(rhoMax < rhoX[i] ) rhoMax = rhoX[i];
	}
	fprintf(data,"# rhoMax = %f\n# rhoMin = %f\n",rhoMax,rhoMin);

	// output rho(y)
	sprintf(f,"rho_y.dat");
	data = fopen(f,"w");
	fprintf(data,"# y\tRhoy\tint Rhoy\n");

	nCheck = 0.0;
	rhoMin = rhoY[0];
	rhoMax = rhoY[0];
	for(i=0;i<ybin;i++){
		if(i<ybin-1) nCheck += (rhoY[i] + rhoY[i+1])/2.0*double(binSize);
		fprintf(data,"%.3f\t%e\t%e\n",i*binSize,rhoY[i],nCheck);
		if(rhoMin > rhoY[i]) rhoMin = rhoY[i];
		if(rhoMax < rhoY[i] ) rhoMax = rhoY[i];
	}
	fprintf(data,"# rhoMax = %f\n# rhoMin = %f\n",rhoMax,rhoMin);


	// output rho(x)
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

	//	get the average max and min
 //	average from max to 0.99 max is the average of max
 //	average from min to 1.01 min is the average of min
 
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

	fclose(data);

	return 0;
}

