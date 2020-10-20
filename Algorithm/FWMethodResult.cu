#include "FWMethodResult.cuh"
#include "cuda_implementation/cudaErrorCheck.cuh"
#include <cstdio>
#include <cmath>
#include <cstring>

#define bufSize 1024
#define TARGET 100000

static double m0 = 1.67261e-27;
static double q = 1.60219e-19;
static double c = 2.99793e8;
static double T0 = m0 * c * c/(q*1e9);

void FWMethodResult::runAlgorithm(ParamsCarrier *singleTone){
    int tem1, tem2, tem3;
    double Rig, p1AU, Tkin, w0, w, p100AU, r, sumac;
    double binb[25], binw[25], binc[25];
    double spe1e2[100] = {0}, spe1e2N[100] = {0};
    double spe4e2[400] = {0}, spe4e2N[400] = {0};
    double spe1e3[1000] = {0},spe1e3N[1000] = {0};
    double spelog[24] = {0};
    for(int i = 0 ; i < 25 ; i++){
        binb[i] = exp((((i+1)/5.0)-3.0)*log(10));
    }
    for(int i = 0 ; i < 24 ; i++){
        binw[i] = binb[i+1]-binb[i];
        binc[i] = (binb[i+1]+binb[i])/2.0;
    }
    FILE *inputFile = fopen("log.dat","r");
    int numberOfIterations = charcount(inputFile) - 1 ;
    int targetArray[] = {numberOfIterations };
    FILE *outputFile; 
    for(int j = 0; j < 1;  ++j){
        rewind(inputFile);
        int targetLine = targetArray[j];
        int actualLine = 0 ;
        for(int i = 0 ; i < numberOfIterations / targetArray[j] ; i++){
            actualLine = 0;
            while (targetLine != actualLine){
                ++actualLine;
                int reader = fscanf(inputFile, " %lf  %lf  %lf  %lf %lf \n", &p100AU, &p1AU,&r,&w0,&sumac);
                if(reader == -1){
			return; 
		}
		Rig = p1AU*c/q;
                Tkin = sqrt((T0*T0*q*q*1e9*1e9) + (q*q*Rig*Rig)) - (T0*q*1e9);
                Tkin = Tkin/(q*1e9);
                w0 = (m0*m0*c*c*c*c) + (p100AU*p100AU*c*c);
                w0 = exp(-1.85*log(w0))/p100AU;
                w0 = w0/1e45;
                w = w0*p1AU*p1AU*exp(sumac);
                if (r<0.3) {
                    w = 0.0;
                }
                for(int i = 0; i < 24; i++) {
                    if (Tkin>binb[i] && Tkin<binb[i+1]) {
                        spelog[i] += w;
                    }
                }
                tem1 = (int)trunc(Tkin);
                if (tem1>0 && tem1<101){
                    spe1e2[tem1 - 1] += w;
                    spe1e2N[tem1 - 1] += 1;
                }
                tem2 = (int)trunc((Tkin+0.05)*10);
                if (tem2>0 && tem2<1001) {
                    spe1e3[tem2 - 1] += w;
                    spe1e3N[tem2 - 1] += 1;
                }
                tem3 = (int)trunc((Tkin+0.125)*4);
                if (tem3>0 && tem3<401){
                    spe4e2[tem3 - 1] += w;
                    spe4e2N[tem3 - 1] += 1;
                }
            }
        }
    }
    fclose(inputFile);
	int isCsv = singleTone->getInt("csv",0); 
    if(isCsv){
    	outputFile = fopen("output_1e3bin.csv", "w");
	fprintf(outputFile, "%s,%s,%s\n", "Tkin", "spe1e3N", "spe1e3");
    }else{
	outputFile = fopen("output_1e3bin.dat", "w"); 
    }	
    for(int i = 0 ; i < 1000 ; i++){
	if(isCsv){
    		fprintf(outputFile, "%3.2f,%.14E,%.14E\n", ((float)(i + 1) / 10), spe1e3N[i], spe1e3[i]);
    	}else{
    		fprintf(outputFile, "%3.2f %.14E %.14E\n", ((float)(i + 1) / 10), spe1e3N[i], spe1e3[i]);
    	}
    }
    fclose(outputFile);
    if(isCsv){
    	outputFile = fopen("output_1e2bin.csv", "w");
	fprintf(outputFile, "%s,%s,%s\n", "Tkin", "spe1e2N", "spe1e2");
    }else{
    	outputFile = fopen("output_1e2bin.dat", "w");
    }	
    for(int i = 0 ; i < 100 ; i++){
	if(isCsv){
    		fprintf(outputFile, "%d,%.14E,%.14E\n", i + 1,spe1e2N[i], spe1e2[i]);
    	}else{
        	fprintf(outputFile, "%d %.14E %.14E\n", i + 1,spe1e2N[i], spe1e2[i]);
    	}
    }
    fclose(outputFile);
	if(isCsv){    
	outputFile = fopen("output_logbin.csv", "w");
	fprintf(outputFile, "%s,%s\n", "binc","spelog/binw");
    }else{
    	outputFile = fopen("output_logbin.dat", "w");
    }	
    for(int i = 0 ; i < 24 ; i++){
	if(isCsv){
 		fprintf(outputFile, "%.14E,%.14E\n", binc[i],spelog[i]/binw[i]);
    	}else{
 		fprintf(outputFile, "%.14E %.14E\n", binc[i],spelog[i]/binw[i]);
    	}
    }
    fclose(outputFile);
    if(isCsv){    
    	outputFile = fopen("output_4e2bin.csv", "w");
	fprintf(outputFile, "%s,%s,%s\n", "Tkin", "spe4e2N", "spe4e2");
    }else{
    	outputFile = fopen("output_4e2bin.dat", "w");
    }	
    for(int i = 0 ; i < 400 ; i++){
	if(isCsv){
        	fprintf(outputFile, "%3.4f,%.14E,%.14E\n",((float)(i + 1) / 4),spe4e2N[i], spe4e2[i]);
    	}else{
        	fprintf(outputFile, "%3.4f %.14E %.14E\n",((float)(i + 1) / 4),spe4e2N[i], spe4e2[i]);
    	}
    }
    fclose(outputFile);
}

int FWMethodResult::charcount( FILE *const fin ) {
    char buf[bufSize];
    int count = 0;
    while (fgets(buf, sizeof(buf), fin) != NULL) {
        ++count;
    }
    return count;
}
