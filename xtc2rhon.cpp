/* Jul 16, 2018
 * Check input format and work for both .xtc and .trr files
 * use new instead of malloc
 * output time series to files traj2rhon_nx_ny_nz.dat
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

#include "/home/zhitongj/Git/analyzeTraj/getZs.cpp"	// move the c.o.m of the system to the center of the box

#define PI 3.14159265358979323846

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
		return 1;
	}

	// check input format
	if( argc < 8 || (argc - 8) % 3 != 0 ){
		printf("# Illegal input...\n");
		printf("#    input format : traj files name, firstframe, lastframe, natomPmol, nx ny nz...\n");
		printf("#      example: xtc2rhon.out traj.xtc 0 -1 1 0 0 1 0 0 2\n");
		return 1;
	}


	// input variables
	int firstframe = 500;				
	if(argc >= 3) firstframe = atoi(argv[2]);
	int lastframe = -1;		// ps
	if(argc >= 4) lastframe = atoi(argv[3]);
	int natomPmol = 1;
	if(argc >= 5) natomPmol = atoi(argv[4]);


	int nWave =  int( ( argc - 5 ) / 3 );
	int **wn = new int* [nWave];
	for(int i = 0 ;i < nWave;i++){
		wn[i] = new int [3];
		wn[i][0] = atoi(argv[5+i*3]);
		wn[i][1] = atoi(argv[6+i*3]);
		wn[i][2] = atoi(argv[7+i*3]);
	}

	printf("# first frame = %d \n",firstframe);
	printf("# last frame = %d \n",lastframe);
	printf("# atoms per molecules= %d \n",natomPmol);
	printf("# wn:nx   ny   nz \n");
	for(int i = 0 ;i < nWave;i++){
		printf("#    %2d  %2d  %2d \n",wn[i][0],wn[i][1],wn[i][2]);
	}
	printf("#---- end of input ----\n");


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

	rvec *vel;			// atom coordinates in a 2-D matrix
	rvec *force;		// atom coordinates in a 2-D matrix
	float lambda;
	vel = new rvec [natoms];
	force = new rvec [natoms];

	float prec;	
	coor = new rvec [natoms];


	// other variables
	int i;

	double xs,ys,zs;
	double xi,yi,zi;
	double *rhonRe = new double [nWave];
	double *rhonIm = new double [nWave];

	// output related
	char f[40];
	FILE *data;
	
	
	// loop over all wave number and open the output files
	for(i=0;i<nWave;i++){
		sprintf(f,"traj2rhok_nx%d_ny%d_nz%d.dat",wn[i][0],wn[i][1],wn[i][2]);
		data = fopen(f,"w");
		fprintf(data, "# time rhon_Re\trhon_Im   wn = ( %d %d %d ) \n",wn[i][0],wn[i][1],wn[i][2]);
		fclose(data);
	}

	//open xtc file and loop through each frame
	xd=xdrfile_open(argv[1],"r");
	int nCount = 0;

	if (ifTrr == 1) { // input is trr file
		while( ! read_trr(xd, natoms, &step, &time, &lambda, box, coor, vel, force)){
			if(coor == 0){
				printf("Insufficient memory to load .trr file. \n");
				return 1;
			}
			if(step == 0){
			}

		    if(time >= firstframe && ( time <= lastframe || lastframe < firstframe ) ){
				if(step%5000 == 0 ) printf("time = %f\n",time);

				// initialize rhon s
				for(i=0;i<nWave;i++){
					rhonRe[i] = 0.0;
					rhonIm[i] = 0.0;
				}

				// move the c.o.m of the system to the center of the box in each direction
				xs = getZs(coor,box,0.1,natoms,natomPmol,'x');	
				ys = getZs(coor,box,0.1,natoms,natomPmol,'y');
				zs = getZs(coor,box,0.1,natoms,natomPmol,'z');

				// Calculate new rho_n Real and Imaginary components
				for(int ni=0;ni<natoms/natomPmol;ni++){
					if(coor[ni*natomPmol][0] - xs >= 0) xi = coor[ni*natomPmol][0] - xs;
					else xi = coor[ni*natomPmol][0] - xs + box[0][0];
					if(xi == box[0][0]) xi = 0;
				
					if(coor[ni*natomPmol][1] - ys >= 0) yi = coor[ni*natomPmol][1] - ys;
					else yi = coor[ni*natomPmol][1] - ys + box[1][1];
					if(yi == box[1][1]) yi = 0;

					if(coor[ni*natomPmol][2] - zs >= 0) zi = coor[ni*natomPmol][2] - zs;
					else zi = coor[ni*natomPmol][2] - zs + box[2][2];
					if(zi == box[2][2]) zi = 0;

					for(i=0;i<nWave;i++){
						rhonRe[i] += cos(2.0*wn[i][0]*PI*xi/box[0][0]+2.0*wn[i][1]*PI*yi/box[1][1] + 2.0*wn[i][2]*PI*zi/box[2][2]);
						rhonIm[i] += -sin(2.0*wn[i][0]*PI*xi/box[0][0]+2.0*wn[i][1]*PI*yi/box[1][1] + 2.0*wn[i][2]*PI*zi/box[2][2]);
					}
				}

				// output
				for(i=0;i<nWave;i++){
					sprintf(f,"traj2rhok_nx%d_ny%d_nz%d.dat",wn[i][0],wn[i][1],wn[i][2]);
					data = fopen(f,"a");
					fprintf(data, "%.2f\t%e\t%e\n",time,rhonRe[i],rhonIm[i]);
					fclose(data);
				}
				nCount++;
		    }
		}
		
	}
	else { // input is xtc file
		while(! read_xtc(xd, natoms, &step, &time, box, coor, &prec) ){
			if(coor == 0){
				printf("Insufficient memory to load .trr file. \n");
				return 1;
			}
			if(step == 0){
			}

		    if(time >= firstframe && ( time <= lastframe || lastframe < firstframe ) ){
				if(step%5000 == 0 ) printf("time = %f\n",time);

				// initialize rhon s
				for(i=0;i<nWave;i++){
					rhonRe[i] = 0.0;
					rhonIm[i] = 0.0;
				}

				// move the c.o.m of the system to the center of the box in each direction
				xs = getZs(coor,box,0.1,natoms,natomPmol,'x');	
				ys = getZs(coor,box,0.1,natoms,natomPmol,'y');
				zs = getZs(coor,box,0.1,natoms,natomPmol,'z');

				// Calculate new rho_n Real and Imaginary components
				for(int ni=0;ni<natoms/natomPmol;ni++){
					if(coor[ni*natomPmol][0] - xs >= 0) xi = coor[ni*natomPmol][0] - xs;
					else xi = coor[ni*natomPmol][0] - xs + box[0][0];
					if(xi == box[0][0]) xi = 0;
				
					if(coor[ni*natomPmol][1] - ys >= 0) yi = coor[ni*natomPmol][1] - ys;
					else yi = coor[ni*natomPmol][1] - ys + box[1][1];
					if(yi == box[1][1]) yi = 0;

					if(coor[ni*natomPmol][2] - zs >= 0) zi = coor[ni*natomPmol][2] - zs;
					else zi = coor[ni*natomPmol][2] - zs + box[2][2];
					if(zi == box[2][2]) zi = 0;

					for(i=0;i<nWave;i++){
						rhonRe[i] += cos(2.0*wn[i][0]*PI*xi/box[0][0]+2.0*wn[i][1]*PI*yi/box[1][1] + 2.0*wn[i][2]*PI*zi/box[2][2]);
						rhonIm[i] += -sin(2.0*wn[i][0]*PI*xi/box[0][0]+2.0*wn[i][1]*PI*yi/box[1][1] + 2.0*wn[i][2]*PI*zi/box[2][2]);
					}
				}

				// output
				for(i=0;i<nWave;i++){
					sprintf(f,"traj2rhok_nx%d_ny%d_nz%d.dat",wn[i][0],wn[i][1],wn[i][2]);
					data = fopen(f,"a");
					fprintf(data, "%.2f\t%e\t%e\n",time,rhonRe[i],rhonIm[i]);
					fclose(data);
				}
				nCount++;
		    }
		}
		
	}

	printf("# Finish reading .trr file...\n");
	printf("# Counted frame = %d...\n",nCount);

	return 0;
}
