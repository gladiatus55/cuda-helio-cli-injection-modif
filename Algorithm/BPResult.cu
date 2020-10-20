#include "BPResult.cuh"
#include "cuda_implementation/cudaErrorCheck.cuh"
#include <cstdio>
#include <cmath>
#include <cstring>

#define bufSize 1024
#define TARGET 100000

static double m0 = 1.67261e-27;
static double q = 1.60219e-19;
static double c = 2.99793e8;
static double T0w = m0 * c * c;

void BPResult::runAlgorithm(ParamsCarrier *singleTone){
    int tem1, tem2, tem3;
    double w, Rig, p1AU, Tkin, r,p, Tkinw, Rig1AU, Tkininj;
    double binb[300], binw[300], binc[300];
    double spe1e2[100] = {0}, spe1e2N[100] = {0};
    double spe4e2[400] = {0}, spe4e2N[400] = {0};
    double spe1e3[1000] = {0},spe1e3N[1000] = {0};
    double spelog[300] = {0}, spelogN[300] = {0};
    int NT = 80;
    double Tmin = 0.01;
    double Tmax = 200.0;
    double dlT = log10(Tmax/Tmin)/NT;
    double X = log10(Tmax);
    double larg = exp((dlT/2.0)*log(10.0))-exp((-1.0*dlT/2.0)*log(10.0));
    for(int i = 0 ; i < NT ; i++){
        spelog[i] = 0;
        double tem = X - (i*dlT);
        binb[i] = exp(tem*log(10.0));
        binc[i] = sqrt(binb[i]* exp((tem+dlT)*log(10.0)));
        binw[i] = larg*binc[i];
    }
    FILE *inputFile = fopen("log.dat","r");
    int numberOfIterations = charcount(inputFile) - 1;
    int targetArray[] = {numberOfIterations };
    FILE *outputFile;
    for(int j = 0; j < 1;  ++j){
        rewind(inputFile);
        int targetLine = targetArray[j];
        int actualLine = 0;
        for(int i = 0 ; i < numberOfIterations / targetArray[j] ; i++){
            actualLine = 0;
            while (targetLine != actualLine){
                ++actualLine;
                int reader = fscanf(inputFile, "%lf %lf %lf %lf\n", &Tkininj,&Tkin,&r,&w);
                if(reader == -1){
			return; 
		}
		Tkinw = Tkininj*1e9*q;
                Rig = sqrt(Tkinw*(Tkinw+(2*T0w)))/q;
                p = Rig*q/c;
                Tkinw = Tkin*1e9*q;
                Rig1AU = sqrt(Tkinw*(Tkinw+(2*T0w)))/q;
                p1AU = Rig1AU*q/c;
                w = (m0*m0*c*c*c*c) + (p*p*c*c);
                w = exp(-1.85*log(w))/p;
                w = w*p1AU*p1AU;
                for(int ii = 0 ; ii < 299 ; ii++){
                    if ((Tkin>binb[ii+1])&&(Tkin<binb[ii])) {
                        spelog[ii+1] = spelog[ii+1] + w;
                        spelogN[ii + 1] = spelogN[ii + 1] + 1;
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
	fprintf(outputFile, "%s,%s,%s\n", "binc", "spelogN", "spelog/binw");
    }else{
    	outputFile = fopen("output_logbin.dat", "w");
    }	
    for(int i = 2 ; i < 80 ; i++){
	if(isCsv){
    		fprintf(outputFile, "%3.2f,%.14E,%.14E\n", binc[i], spelogN[i],1e-6*spelog[i]/binw[i]);
    	}else{
        	fprintf(outputFile, "%3.2f %.14E %.14E\n", binc[i], spelogN[i],1e-6*spelog[i]/binw[i]);
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

int BPResult::charcount( FILE *const fin ) {
    char buf[bufSize];
    int count = 0;
    while (fgets(buf, sizeof(buf), fin) != NULL) {
        ++count;
    }
    return count;
}
