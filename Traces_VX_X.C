//
//  converter_V1_0.C
//  
//
//  Created by Ponce, Francisco on 1/8/16.
//
//

#define MAX_TO_READ 80000000
#define N_detectors 4
#define N_HIST 30000


#include <iostream>
#include <fstream>
#include <stdio.h>
#include <snprintf.h>
#include <math.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TMinuit.h>
#include <TMath.h>


#include <TGraphErrors.h>
#include <TSpectrum.h>
#include <TCanvas.h>

//*/ Standard structure
struct info{
    char path[500]; // This is the main file path location
    char name[500]; // This is the main file name without extension
    char ext[10]; // This is the main file extension
    char filename[500]; // This is a blank space for creating a file name when needed.
    int current; // Needs to be initialized to 0 when created this marks the current location in the data being processed.
    int last; // Needs to be initialized to 0 when created this is the end of the data for the file.
    int counter; // Needs to be initialized to 0 when created general number tracker?
    int channel;
    int event[N_HIST + 1];
    
    TH1F *histogram[N_HIST + 1];
    TGraph *graph[N_HIST + 1];
    double trace[N_HIST + 1][4096];
    double xaxis[4096];
};

//*/ TEST FILES: Produce files with < 32768 bins
char tfile[500] = "111112g_Be7_12p1A_XIA_p108d";
char tfile00[500] = "111112-d_b14-p228b";
char tfile01[500] = "111112-d_b14-p228c";
char tfile02[500] = "111112-d_b14-p228e";
char tfile03[500] = "111112-d_b14-p229a";
char tfile04[500] = "111112-d_b14-p217f";
char tfile05[500] = "111112-h_b14-p281a";
char *folder[13] = {tfile00, tfile01, tfile02, tfile03, tfile04, tfile05};
char fpath[500] = "/Users/ponce10/Root_Folder/Data/233U_to_229Th/Low_Activity_Four_Pinholes/160506-160603_STJ-d_Traces_XIA_233U/";
//char tpath[500] = "/Users/ponce10/Root_Folder/Data/Code_Testing/Traces/";
char tpath[500] = "/Applications/Data-Processed/Livermore/Test/";
char text[500] = ".bin";

// End File names and extentions
info *outgamma[4];

// Start functions

/* **************************************************************************************************** */
double_t func_00(double_t *x, double_t *par){
    /* **************************************************************************************************** */
    /*/
     This function is a superposition of an exponential decay and gaussian function.
     y(x, par) = par[0]/(2pi*par[2]**2)**0.5*exp(-(E-par[1])**2/(2*par[2]**2)) + par[3]*exp(-E/par[4])
     
     Number of parameters: 5
     /*/
    double x0 = x[0];
    double_t f = par[3]*TMath::Exp(-1.0*x0/par[4]);
    f += par[0]/pow(2*TMath::Pi()*par[2]*par[2], 0.5)*TMath::Gaus(*x, par[1], par[2]);
    return f;
}

/* **************************************************************************************************** */
int xiaload(char* filename, info *detector[N_detectors]){
    /* **************************************************************************************************** */
    
    /*
     This function will convert the xia binary list mode data to ascii format with tags.
     */
    
    // Initialize Function Variables:
    int num_read = 0;
    int num_headers = 0;
    int det_num = -1;
    int ecount = 0;
    int exit = 0;
    int vnum = 0;
    int tcheck[4] = {0, 0, 0, 0};
    int loc = 0;
    unsigned char data[16384];
    unsigned char buffer[512];
    unsigned char a,b;
    
    info *det;
    FILE *in;
    // End
    
    // Open binary XIA file for parsing
    in = fopen(filename,"rb");
    if (in == NULL) printf("Failed to open file: %s\nExiting\n",filename);
    // Loop through the file and copy content to storage
    else {
        std::cout << "Loading data:\n" << filename;
        std::cout << "\nTotal Counts:\n";
        while (!feof(in) && (num_read < MAX_TO_READ)) {
            fscanf(in,"%c%c",&a,&b);
            if((a == 0x55) && (b == 0xAA)) {//HEADER STARTS WITH A 16 BIT TOKEN Ox55AA, skip it
                //fread(buffer,1,508,in); //read in header data
                for(int i = 0; i < 508; ++i) fread(&buffer[i], 1, 1, in);
                if(ecount == 0){
                    vnum = buffer[0] + 256*buffer[1];
                    printf("Buffer Header Size: %d\n", vnum);
                    vnum = buffer[2] + 256*buffer[3];
                    printf("Mapping Mode: %d\n", vnum);
                    vnum = buffer[4] + 256*buffer[5];
                    printf("Run Number: %d\n", vnum);
                    vnum = buffer[6] + 256*buffer[7] + 256*256*buffer[8] + 256*256*256*buffer[9];
                    printf("Sequential Buffer Number: %d\n", vnum);
                    vnum = buffer[10] + 256*buffer[11];
                    printf("Buffer ID: %d\n", vnum);
                    vnum = buffer[12] + 256*buffer[13];
                    printf("Number of events in buffer: %d\n", vnum);
                    vnum = buffer[14] + 256*buffer[15] + 256*256*buffer[16] + 256*256*256*buffer[17];
                    printf("Starting Event Number: %d\n", vnum);
                    vnum = buffer[18] + 256*buffer[19];
                    printf("Module #: %d\n", vnum);
                    vnum = buffer[20] + 256*buffer[21];
                    printf("AMREV: %d\n", vnum);
                    vnum = buffer[22] + 256*buffer[23];
                    printf("CODEREV: %d\n", vnum);
                    vnum = buffer[24] + 256*buffer[25];
                    printf("CODEVAR: %d\n", vnum);
                    vnum = buffer[26] + 256*buffer[27];
                    printf("SYSREV: %d\n", vnum);
                    vnum = buffer[28] + 256*buffer[29];
                    printf("SYSVAR: %d\n", vnum);
                    vnum = buffer[30] + 256*buffer[31];
                    printf("FIPPIREV: %d\n", vnum);
                    vnum = buffer[32] + 256*buffer[33];
                    printf("FIPPIVAR: %d\n", vnum);
                    vnum = 0;
                    for (int i = 34; i < 46; i++) vnum += buffer[i];
                    printf("Reserved 0: %d\n", vnum);
                    vnum = buffer[46] + 256*buffer[47] + 256*256*buffer[48] + 256*256*256*buffer[49];
                    printf("Number of words in buffer: %d\n", vnum);
                    vnum = buffer[50] + 256*buffer[51];
                    printf("Bank0ChEna: %d\n", vnum);
                    vnum = buffer[52] + 256*buffer[53];
                    printf("Bank1ChEna: %d\n", vnum);
                    vnum = buffer[54] + 256*buffer[55];
                    printf("Bank2ChEna: %d\n", vnum);
                    vnum = buffer[56] + 256*buffer[57];
                    printf("Bank3ChEna: %d\n", vnum);
                    vnum = buffer[58] + 256*buffer[59];
                    printf("Bank0TraceEna: %d\n", vnum);
                    vnum = buffer[60] + 256*buffer[61];
                    printf("Bank1TraceEna: %d\n", vnum);
                    vnum = buffer[62] + 256*buffer[63];
                    printf("Bank2TraceEna: %d\n", vnum);
                    vnum = buffer[64] + 256*buffer[65];
                    printf("Bank3TraceEna: %d\n", vnum);
                    vnum = 0;
                    for (int i = 66; i < 124; i++) vnum += buffer[i];
                    printf("Reserved 0: %d\n", vnum);
                    vnum = buffer[124] + 256*buffer[125];
                    printf("List Mode Variant: %d\n", vnum);
                    vnum = buffer[126] + 256*buffer[127];
                    printf("Words per Event: %d\n", vnum);
                    vnum = buffer[128] + 256*buffer[129];
                    printf("Number of event in buffer: %d\n", vnum);
                    vnum = 0;
                    for (int i = 130; i < 508; i++) vnum += buffer[i];
                    printf("Reserved 0: %d\n", vnum);
                }
                ++num_headers;
            }
            /*/
            else if ((a == 0xEE) && (b == 0xEE)) {//Event starts with a token too, OxEEEE
                ecount = 0;
                while(true){
                    fscanf(in,"%c%c",&a,&b);
                    if ((a == 0xEE) && (b == 0xEE)) break;
                    if (ecount > 16000) break;
                    ecount++;
                }
                printf("Number of words in a info is half of %d\n", ecount);
                break;
            };
            /*/// Used to count the number of words in a single info;
            else if ((a == 0xEE) && (b == 0xEE)) {//Event starts with a token too, OxEEEE
                ++num_read;
                for(int i = 0; i < 8198; ++i) fread(&data[i],1,1,in);
                det_num = data[0] + 256*data[1]; // Channel ID
                ecount++;
                for (int i = 0; i < N_detectors; i++){
                    if (detector[i]->channel == det_num){
                        det = detector[i];
                        if(det_num == 0 && det->current != 0){
                            loc = 0;
                            for(int j = 0; j < 4; j++) loc += tcheck[j];
                            if (loc != 4) {
                                for (int j = 0; j < 4; j++) if (tcheck[j]){
                                    detector[j]->current--;
                                    detector[j]->last--;
                                    detector[j]->histogram[detector[j]->current]->Reset();
                                }
                            }
                            else for (int j = 0; j < 4; j++) tcheck[j] = 0;
                        }
                        tcheck[det_num/2] = 1;
                        
                        if (det->current < N_HIST){
                            det->event[det->current] = 0;
                            for (int i = 2; i < 6; i++) det->event[det->current] += data[i]*pow(256, i - 2); // Event Number
                            for (int i = 3; i < 4099; i++) det->histogram[det->current]->Fill(i - 2.5, data[2*i] + 256*data[2*i + 1]);
                            //for (int i = 3; i < 4099; i++) det->histogram[det->current]->SetPoint(det->histogram[det->current]->GetN(), i - 3, data[2*i] + 256*data[2*i + 1]);
                            for (int i = 3; i < 4099; i++) det->trace[det->current][i - 3] = data[2*i] + 256*data[2*i + 1];
                            det->current++;
                            det->last++;
                            
                        }
                    }
                    if (detector[i]->current == N_HIST) exit++;
                }
                std::cout << ecount << "\r";
                if (N_detectors <= exit) break;
                else exit = 0;
            };
            //*/
        };
        fclose(in);
        printf("\n");
        exit = 1;
    }
    return exit;
}

/* **************************************************************************************************** */
int tsave(info *det){
/* **************************************************************************************************** */

    //Initialize Variables
    char str[500];
    
    ofstream file;
    //End
    
    
    snprintf(str, sizeof(str), "%s%s_Traces_Ch_%d.txt", fpath, det->name, det->channel);
    file.open(str);
    if(file.fail()) printf("Failed to open file \n%s\n", str);
    else {
        //printf("Saving traces under\n%s\n", str);
        file << "000_xaxis\t";
        for (int i = 0; i < N_HIST + 1; i++) file << i << "_Trace" << (i == N_HIST ? "\n" : "\t");
        for (int i = 0; i < 4096; i++){
            file << det->xaxis[i] << "\t";
            for (int j = 0; j < N_HIST + 1; j++) file << det->trace[j][i] << (j == N_HIST ? "\n" : "\t");
        }
        file.close();
    }
    
    
    
    return 1;
}

/* **************************************************************************************************** */
int plot(info *pointer, int length = 0, int gap = 0, float dtime = 0.4){
/* **************************************************************************************************** */
    
    TGraph *graph = new TGraph(4096);
    TGraph *filter = new TGraph(4096);
    float pre = 0;
    float post = 0;
    float base = 0;
    length = length*1.0/dtime;
    gap = gap*1.0/dtime;
    
    for (int i = 250; i < 750; i++) base += pointer->trace[pointer->current][i]*1.0/500;
    for (int i = 0; i < 4096; i++) {
        graph->SetPoint(graph->GetN(), dtime*i, pointer->trace[pointer->current][i] - base);
    }
    //graph->Draw();
    
    if (length) {
        for (int i = (int)(length); i < 4096 - (int)(length + gap); i++) {
            pre = 0;
            post = 0;
            for (int j = i - length; j < i; j++) pre += pointer->trace[pointer->current][j]*1.0/length;
            for (int j = i + gap; j < i + length + gap; j++) post += pointer->trace[pointer->current][j]*1.0/length;
            filter->SetPoint(filter->GetN(), dtime*i, post - pre);
        }
        filter->Draw("same");
    }
    std::cout << "Pulse Trace: " << pointer->current << "\n";
    pointer->current++;
    return 1;
}

/* **************************************************************************************************** */
double vpp(double trace[4096], int len = 16, int gap = 2){
/* **************************************************************************************************** */
    
    // Initialize Function Variables
    float base = 0;
    double full[4096];
    double pre = 0;
    double post = 0;
    // End
    for (int i = 0; i < 4096; i++) full[i] = 0;
    for (int i = (int)(len); i < 4096 - (int)(len + gap); i++) {
        pre = 0;
        post = 0;
        for (int j = i - len; j < i; j++) pre += trace[j]*1.0/len;
        for (int j = i + gap; j < i + len + gap; j++) post += trace[j]*1.0/len;
        full[i] = post - pre;
    }
    //printf("Vpp: ");
    //for (int i = len; i < len + 5; i++) printf("%f\t", full[i]);
    //printf("\n");
    return TMath::MaxElement(4096, full);
}

/* **************************************************************************************************** */
double area(double trace[4096], int length = 4096){
    /* **************************************************************************************************** */
    
    // Initialize Function Variables
    double base = 0;
    double full[4096];
    // End
    base = 0;
    for (int j = 50; j < 550; j++) base += trace[j]*1.0/500;
    for (int j = 0; j < length; j++) full[j] = trace[j] - base;
    base = 0.4*(full[0] + full[4095])/2;
    for (int j = 1; j < 4095; j++) base += 0.4*full[j];
    //printf("Area: ");
    //for (int i = 0; i < 5; i++) printf("%f\t", full[i]);
    //printf("\n");
    return base;
}

/* **************************************************************************************************** */
int tracecount(char* filename){
/* **************************************************************************************************** */
    
    /*
     This function will convert the xia binary list mode data to ascii format with tags.
     */
    // Initialize Function Variables:
    int num_read = 0;
    int num_headers = 0;
    int det_num[2*N_detectors];
    int exit = 0;
    unsigned char data[2*8192];
    unsigned char a,b;
    FILE *in;
    TH1I *hist = new TH1I(filename, filename, 8, 0, 7);
    // End
    
    for (int i = 0; i < 2*N_detectors; i++) det_num[i] = 0;
    
    // Open binary XIA file for parsing
    in = fopen(filename,"rb");
    if (in == NULL) printf("Failed to open file: %s\nExiting\n",filename);
    // Loop through the file and copy content to storage
    else {
        std::cout << "Loading data:\n" << filename;
        std::cout << "\nTotal Counts:\n";
        while (!feof(in) && (num_read < MAX_TO_READ)) {
            fscanf(in,"%c%c",&a,&b);
            if((a == 0x55) && (b == 0xAA)) {//HEADER STARTS WITH A 16 BIT TOKEN Ox55AA, skip it
                fread(data,1,508,in); //read in header data
                ++num_headers;
            }
            else if ((a == 0xEE) && (b == 0xEE)) {//Event starts with a token too, OxEEEE
                ++num_read;
                for(int i = 0; i < 8198; ++i) fread(&data[i],1,1,in);
                det_num[data[0]]++; // Determine channel of event
                hist->Fill(data[0]);
                std::cout << num_read << "\r";
            }
        }
    }
    fclose(in);
    printf("\n");
    exit = 1;
    hist->Draw();
    return exit;
}

/* **************************************************************************************************** */
int calculate(char* filename, info *detector[N_detectors]){
    /* **************************************************************************************************** */
    /*
     This function will convert the xia binary list mode data to ascii format with tags.
     */
    
    for (int i = 0; i < 4; i++) detector[i]->channel = 2*i;
    
    // Initialize Function Variables:
    int ecount = 0;
    int tcheck[4] = {0, 0, 0, 0};
    int loc = 0;
    
    int num_read = 0;
    int num_headers = 0;
    int det_num = -1;
    int exit = 0;
    int counter = 0;
    char str[500];
    unsigned char data[2*8192];
    unsigned char a,b;
    
    ofstream outfile;
    info *det;
    FILE *in;
    // End
    
    for (int i = 0; i < N_detectors; i++) {
        detector[i]->current = 0;
        detector[i]->last = 0;
    }
    sprintf(str, "%s%s.bin", tpath, filename);
    // Open binary XIA file for parsing
    in = fopen(str,"rb");
    if (in == NULL) printf("Failed to open file: %s\nExiting\n", str);
    // Loop through the file and copy content to storage
    else {
        std::cout << "Loading data:\n" << str;
        sprintf(str, "%s%s%s", fpath, filename, "_Area.txt");
        outfile.open(str);
        outfile << "STJ_0-Area\tSTJ_0-Volt\tSTJ_2-Area\tSTJ_2-Volt\t";
        outfile << "STJ_4-Area\tSTJ_4-Volt\tSTJ_6-Area\tSTJ_6-Volt\n";
        std::cout << "Saving to:\n" << str;
        std::cout << "\nTotal Counts:\n";
        while (!feof(in) && (num_read < MAX_TO_READ)) {
            fscanf(in,"%c%c",&a,&b);
            if((a == 0x55) && (b == 0xAA)) {//HEADER STARTS WITH A 16 BIT TOKEN Ox55AA, skip it
                fread(data,1,508,in); //read in header data
                ++num_headers;
            }
            else if ((a == 0xEE) && (b == 0xEE)) {//Event starts with a token too, OxEEEE
                ++num_read;
                for(int i = 0; i < 8198; ++i) fread(&data[i],1,1,in);
                det_num = data[0] + 256*data[1]; // Channel ID
                ecount++;
                if (ecount%10000 == 0) printf("%d\r", ecount);
                for (int i = 0; i < N_detectors; i++) if (detector[i]->channel == det_num){
                    det = detector[i];
                    tcheck[det_num/2] = 1;
                    for (int i = 3; i < 4099; i++) det->trace[0][i - 3] = 0; // Zero out trace
                    for (int i = 3; i < 4099; i++) det->trace[0][i - 3] = data[2*i] + 256*data[2*i + 1]; // Copy over trace
                    det->current++;
                    det->last++;
                    
                    if (det->current < N_HIST + 1){
                        det->event[det->current] = 0;
                        for (int i = 2; i < 6; i++) det->event[det->current] += data[i]*pow(256, i - 2); // Event Number
                        for (int i = 3; i < 4099; i++) det->trace[det->current][i - 3] = 0; // Zero out trace
                        for (int i = 3; i < 4099; i++) det->trace[det->current][i - 3] = data[2*i] + 256*data[2*i + 1]; // Copy over trace
                    } // Save a few traces for later
                }
                loc = 0;
                for(int j = 0; j < 4; j++) loc += tcheck[j];
                if (loc == 4){
                    for (int i = 0; i < N_detectors; i++){
                        det = detector[i];
                        outfile << area(det->trace[0]) << "\t";
                        outfile << vpp(det->trace[0]) << (i == N_detectors - 1 ? "\n" : "\t");
                    }
                    for (int j = 0; j < 4; j++) tcheck[j] = 0;
                }
            }
            if(ecount == 400000) break;
        }
        fclose(in);
        outfile.close();
        exit = 1;
    }
    return exit;
}

/* **************************************************************************************************** */
int process(int type = 0,  char file[500] = tfile, int randomvalue = 0, char path[500] = tpath){
/* **************************************************************************************************** */
    /* This function is intended to perform the start to finish main function of processing either
       CAEN or XIA data.
    */
    if(type == 0){
        info *data[N_detectors];
        for (int i = 0; i < N_detectors; i++) {
            data[i] = new info;
            data[i]->current = 0;
            data[i]->last = 0;
            data[i]->channel = 2*i;
            sprintf(data[i]->filename, "%s%s%s", path, file, text);
            sprintf(data[i]->path, "%s", path);
            sprintf(data[i]->name, "%s", file);
            sprintf(data[i]->ext, "%s", text);
            outgamma[i] = data[i];
            
            char str[500];
            for (int j = 0 ; j < N_HIST; j++){
                data[i]->graph[j] = new TGraph[4096];
                snprintf(str, sizeof(str), "Det_%d-info_%d", data[i]->channel, j);
                data[i]->histogram[j] = new TH1F(str, str, 4096, 0, 4096);
            }
        }
        xiaload(data[0]->filename, data);
        
        float hold = 0;
        int length = 16;
        int gap = 2;
        for (int i = 0; i < N_detectors; i++){
            data[i]->last = 2;
            if (0 <= randomvalue) for (int j = 0; j < 4096; j++) {
                data[i]->xaxis[j] = 1*j;
                data[i]->trace[0][j] = data[i]->trace[randomvalue][j];
                data[i]->trace[1][j] = 0;
                if (length <= j && j < 4096 - length - gap) {
                    hold = 0;
                    for (int k = j - length; k < j; k++) hold -= data[i]->trace[randomvalue][k]*1.0/length;
                    for (int k = j + gap; k < j + length + gap; k++) hold += data[i]->trace[randomvalue][k]*1.0/length;
                    data[i]->trace[1][j] = hold;
                }
            }
            tsave(data[i]);
        }
        for (int i = 0; i < N_detectors; i++) data[i]->current = 0;
        
    }
    else if(type == 1){
        info *data[N_detectors];
        for (int i = 0; i < N_detectors; i++) {
            data[i] = new info;
            data[i]->current = 0;
            data[i]->last = 0;
            data[i]->channel = 2*i;
            sprintf(data[i]->filename, "%s%s%s", path, file, text);
            sprintf(data[i]->path, "%s", path);
            sprintf(data[i]->name, "%s", file);
            sprintf(data[i]->ext, "%s", text);
            outgamma[i] = data[i];
            
            char str[500];
            for (int j = 0 ; j < N_HIST; j++){
                data[i]->graph[j] = new TGraph[4096];
                snprintf(str, sizeof(str), "Det_%d-info_%d", data[i]->channel, j);
                data[i]->histogram[j] = new TH1F(str, str, 4096, 0, 4096);
            }
        }
        calculate(data[0]->name, data);
    }
    else if(type == 2){
        info *data[N_detectors];
        for (int i = 0; i < N_detectors; i++) {
            data[i] = new info;
            data[i]->current = 0;
            data[i]->last = 0;
            data[i]->channel = 2*i;
            sprintf(data[i]->filename, "%s%s%s", path, file, text);
            sprintf(data[i]->path, "%s", path);
            sprintf(data[i]->name, "%s", file);
            sprintf(data[i]->ext, "%s", text);
            outgamma[i] = data[i];
        }
        xiaload(data[0]->filename, data);
    }
    return type;
}


