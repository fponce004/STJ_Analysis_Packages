//
//  converter_VX_X.C
//  
//
//  Created by Ponce, Francisco on 10/11/17.
//
//

#define MAX_TO_READ 200000000
//#define MAX_TO_READ 2
#define N_Events 80000000
//#define N_Events 6
#define N_detectors 8
#define NUM_COLUMN 1440
//#define NUM_COLUMN 10000
#define NUM_PEAKS 100
#define NPV 100
#define CAEN 8

#define XIA 2
//#define tBin 32768
#define tBin 8192
#define FWHM 1.0
#define WAVELENGTH 3.5
#define fitmodel false // A value of true implies linear and code will adjust, a comment of false will implement quadratic
#define standard true // true means use polynomial fits to calibration false means use 3.5 = peak difference


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
#include <TVector.h>
#include <TMatrix.h>
#include "TRandom.h"
#include "TRandom3.h"
//#include <TBenchmark.h> // Not actually needed...


struct info{
    int current; // Needs to be initialized to 0 when created this marks the current location in the data being processed.
    int last; // Needs to be initialized to 0 when created this is the end of the data for the file.
    int cslice; // Needs to be initialized to 0 when created this marks the current time slice in the hist column.
    int endslice; // Needs to be initialized to 0 when created this marks the last time slice in the hist
    int id[N_Events]; // Use to mark the type of event
    int indexvalue; // Use to indicate file number...
    int numbin; // Use to indicate the number of bins being used in the hist array...
    char path[500]; // This is the main file path location
    char name[500]; // This is the main file name without extension
    char ext[10]; // This is the main file extension
    char filename[500]; // This is a blank space for creating a file name when needed.
    short channel; // Needs to be initialized to detector channel number.
    
    
    bool veto[N_Events]; // To be used to exclude points for veto
    bool high[N_Events]; // To be used to exclude points for HE coincidence
    bool recoil[N_Events]; // To be used to determine recoil candidacy
    float bin[N_Events];
    short tag[N_Events];
    //double energy[N_Events]; // Used later when the bin values are converted over to energy
    long long time[N_Events];
    
    TH1D *histogram[NUM_COLUMN];
    
    // This stores the time intervals for the histogram cuts
    long long thyme[NUM_COLUMN + 1];
    // Histograms of the data parsed into time cuts
    double hist[NUM_COLUMN][tBin];
    // This is an array for storing the xaxis values of a histogram
    float xaxis[tBin];
    // Gaussian parameter fits : {{fit_number, energy, peak, dpeak, sigma, dsigma, area, darea}_0, ...}
    double param[NUM_PEAKS*NUM_COLUMN][8];
    // Calibration fits : {{fit_number, slope, dslope, offset, doffset, quad, dquad}_0, ... }
    double calib[NUM_COLUMN][7];
    // Output files the columns are: Energy, CS-NV, CS-Coin, CS-10ms, CS-50ms, CS-100ms, S-10ms, S-100ms, S-200ms
    double output[9][5];
    // This will contain the fits for the various bin and background values the 8th is for red_chi
    double fit[250][17];
};

struct alpha_matrix{
    float matrix[tBin + 1][tBin + 1];
    float newt[tBin + 1];
    float old[tBin + 1];
};

// Start File names and extentions


// FILTER FILES: Produce files with < 32768 bins
char ffile[500] = "Slow_b14-p208b"; // Need to hardline to 0, 1, 3
char fpath[500] = "/Users/ponce10/Root_Folder/Data/Code_Testing/FILTER/";
char fext[500] = ".bin";

// These are the files for the bias vs responsivity test...
char xfile[500] = "181027_111112-h_b14-p299b";
char xfile00[500] = "01_Bias_100_b14-p226a";
char xfile01[500] = "02_Bias_020_b14-p226a";
char xfile02[500] = "03_Bias_040_b14-p226a";
char xpath[500] = "/Applications/Data-Processed/Livermore/181026_No_Source_STJ_Test/";
char xext[500] = ".bin";
char *xfolder[20] = {xfile00, xfile01, xfile02};
//*/

// These are the files for the bias vs responsivity test...
//char tfile[500] = "STJ_Run_Vs_Bias";
char tfile[500] = "111112g_Be7_12p1A_XIA_p108d";
char tfile00[500] = "111112d_U234_12p2A_p132b";
char tfile01[500] = "111112d_U234_12p1A_p133c";
char tfile02[500] = "111112d_U234_12p2A_p134e";
char tfile03[500] = "111112g_Be7_12p1A_XIA_p112b";
char tfile04[500] = "111112g_Be7_12p1A_XIA_p112d";
char tfile05[500] = "111112g_Be7_12p1A_XIA_p113b";
char tfile06[500] = "111112g_Be7_12p1A_XIA_p113e";
char tfile07[500] = "111112g_Be7_12p1A_XIA_p114e";
char tfile08[500] = "111112g_Be7_12p1A_XIA_p115b";
char tfile09[500] = "111112g_Be7_12p1A_XIA_p115e";
char tfile10[500] = "111112g_Be7_12p1A_XIA_p115g";
char tfile11[500] = "111112g_Be7_12p05A_XIA_p110a";
char tfile12[500] = "111112g_Be7_12p05A_XIA_p111b";
char tfile13[500] = "111112g_Be7_12p15A_XIA_p110e";
char tfile14[500] = "Ortec";
char tfile15[500] = "";
char tfile16[500] = "";
char tfile17[500] = "";
//char tpath[500] = "/Applications/Data-Processed/Livermore/181125/";
//char tpath[500] = "/Applications/Data-Processed/Livermore/190601/";
//char tpath[500] = "/Applications/Data-Processed/Livermore/190623/";
//char tpath[500] = "/Applications/Data-Processed/Livermore/190626/";
//char tpath[500] = "/Applications/Data-Processed/Livermore/191017/";
//char tpath[500] = "/Applications/Data-Processed/Livermore/191206/";
//char tpath[500] = "/Applications/Data-Processed/Livermore/Test/";
char tpath[500] = "/Applications/Data-Processed/Livermore/Ortec/";
char text[500] = ".bin";
char *tfolder[20] = {tfile00, tfile01, tfile02, tfile03, tfile04, tfile05,
                     tfile06, tfile07, tfile08, tfile09, tfile10, tfile11,
                     tfile12, tfile13, tfile14, tfile15, tfile16, tfile17};
//*/

// These are the files for the bias vs responsivity test...
char cfile[500] = "01_Bias_100_b14-p226a";
char cfile00[500] = "01_Bias_100_b14-p226a";
char cfile01[500] = "02_Bias_020_b14-p226a";
char cfile02[500] = "03_Bias_040_b14-p226a";
char cpath[500] = "/Users/ponce10/Root_Folder/Data/Code_Testing/Responsivity_Vs_Bias/";
char cext[500] = ".bin";
char *cfolder[20] = {cfile00, cfile01, cfile02};
//*/


char ufile[500] = "111112-h_b14-p278e";
char ufile01[500] = "111112-a_208um";
char ufile02[500] = "ADC1_Laser_600Hz_x500_b13-p248_merge";
char upath00[500] = "/Users/ponce10/Root_Folder/Data/Uncertainty_Partition_Time/";
char upath01[500] = "/Users/ponce10/Root_Folder/Data/Uncertainty_Reconstructed_Energy_Bin_Size/";
char upath02[500] = "/Users/ponce10/Root_Folder/Data/Uncertainty_Fit_Range_And_Background/";
char upath03[500] = "/Users/ponce10/Root_Folder/Data/Monte_Carlo_Simulation/";
char upath04[500] = "/Users/ponce10/Root_Folder/Data/Code_Testing/Traces/";
char upath05[500] = "/Volumes/PONCE10/";
char *upath[20] = {upath00, upath01, upath02, upath03};
char uext[500] = ".bin";

char ch_ext[500] = "_Ch_";
// End File names and extentions
char *folder[20] = {upath00, upath01, upath02, upath03};

char testfile00[500] = "Data_E";
char ctpath[500] = "/Users/ponce10/Root_Folder/Data/Code_Testing/TEST";
char *testfile = testfile00;

char mcpath[500] = "/Users/ponce10/Root_Folder/Data/";
char mcfile00[500] = "Monte_Carlo_Simulation";


info *troll;
info *bridge;
TH1D *thist[8];
TH2I *globalhist[N_detectors];
TH1D *dhost;
TH1I *ihost;
// Start functions

/* **************************************************************************************************** */
template <class type>
void display(type* item, int start = 0, int end = 10){
    /* **************************************************************************************************** */
    for(int i = start; i < end; i++) std::cout << "Point " << i << "\t" << item[i] << "\n";
}

/* **************************************************************************************************** */
template <class type>
type round_decimal(type value, int dec_place = 0){
    /* **************************************************************************************************** */
    return 1.0*round(value*pow(10, dec_place))/pow(10, dec_place);
}

/* **************************************************************************************************** */
int xiaload(char* filename, info *detector[N_detectors], double start = 0, double stop = 0){
/* **************************************************************************************************** */

    /*
     This function will convert the xia binary list mode data to ascii format with tags.
     */
    
    // Initialize Function Variables:
    int num_read = 0;
    int energy = 0;
    int num_headers = 0;
    int det_num = -1;
    int events = 0;
    int tag = 0;
    int exit = 0;
    long long time_now;
    long long time_ago[4];
    long long prevtime[N_detectors];
    unsigned char data[1024];
    unsigned char a,b;
    
    info *det;
    FILE *in;
    // End
    
    start = 0; //Just so warning is used impliment later
    // Open binary XIA file for parsing
    in = fopen(filename,"rb");
    if (in == NULL) printf("Failed to open file: %s\nExiting\n",filename);
    // Loop through the file and copy content to storage
    else {
        std::cout << "Loading data:\n" << filename;
        std::cout << "\nTotal Counts:\n";
        while (!feof(in) && (num_read < MAX_TO_READ)) {
            fscanf(in,"%c%c",&a,&b);
            if((a == 0x55) && (b == 0xAA)) {
                fread(data,1,508,in); //read in header data
                ++num_headers;
            }//HEADER STARTS WITH A 16 BIT TOKEN Ox55AA, skip it
            else if ((a == 0xEE) && (b == 0xEE)) {
                ++num_read;
                for(int i = 0; i < 14; ++i) fread(&data[i],1,1,in);
                // Readout single event
                time_now = 0;
                for (int i = 6; i < 14; i++) time_now += (long long)data[i]*pow(256, i - 6); // Determine time of event
                energy = data[4] + 256*data[5]; // Determine energy of event
                det_num = data[0]; // Determine channel of event
                tag = data[2];// Determine tag of event
                events++;
                for (int i = 0; i < N_detectors; i++) if (detector[i]->channel == det_num){
                    det = detector[i];
                    if (false || (det->current == 0 && time_now*20 > 0)|| (det->time[det->current - 1] != time_now*20 && time_now*20 > 0)){
                        det->tag[det->current] = tag;
                        det->bin[det->current] = energy;
                        det->time[det->current] = time_now*20; // The 20 is the cycle unit for XIA files
                        det->current++;
                        det->last++;
                    }// change to true in order to get all duplicate time events...
                    else if (det->current != 0 && det->time[det->current - 1] == time_now*20){
                        det->current--;
                        det->last--;
                    }// Remove Duplicate timestamp events...
                    if (N_Events <= det->last) exit = 1;
                }
                if (events%1000000 == 0) printf("%.2fM\r", events/pow(10, 6));
                if (stop != 0 && stop < time_now*20/(3.6*pow(10, 12))) break;
                if (exit) {
                    printf("A detector has reached its count limit...\n");
                    break;
                }
            };//Event starts with a token too, OxEEEE
        };
        fclose(in);
        printf("\nCounts in each Channel:\n");
        for (int i = 0; i < N_detectors; i++) {
            std::cout << "Channel " << detector[i]->channel << " : " << detector[i]->last << "\n";
        }
        exit = 1;
    }
    return exit;
}

/* **************************************************************************************************** */
int caenload(info* pointer, int type = 0, int start = 0, int stop = 100){
/* **************************************************************************************************** */
    /*
     This script will load both xia-converted files and caen files unprocessed and processed alike based on
     the type value. The start and stop entries are to allow only sections of the data to be loaded. The
     time entry is intended to allow files to be synchronized and will be added at a later time and should
     be in units of minutes.
     
     type       =>  This sets the load method from pure CAEN file .dat with header to prior processed file
     rebin      =>  (int) Sets the rebin if any of the loaded data
     start      =>  (int) Time in hours to start load data
     stop       =>  (int) Time in hours to stop loading data
     
     */
    
    
    // Initialize Function Variables:
    ifstream file;
    string trash;
    int czeros = 0;
    int counter = 0;
    int trash_bin = 0;
    int trash_tag = 0;
    int exit = 0;
    long long trash_time = 0;
    // End
    
    file.open(pointer -> filename);
    std::cout << "\nLoading data:\n" << pointer -> filename;
    if(file.fail()) std::cout << "\nUnable to load file.";
    else{
        std::cout << "\nTotal events:\n";
        if (type == -1) {
            for (int i = 0; i < 4; i++) file >> trash;
            for (int i = 0; i < start; i++) file >> trash_time >> trash_bin;
            
            while(file >> pointer->time[counter] >> pointer->bin[counter]){
                pointer->time[counter] = 10*pointer->time[counter]; // *10 to get to ns
                if (pointer -> bin[counter] <= 0) {czeros++; continue;}
                else if (pointer->time[counter] < start*3.6*pow(10, 12)) continue;
                else if (stop*3.6*pow(10, 12) <= pointer->time[counter]) {counter++; break;}
                counter++;
                std::cout << counter << "\r";
            }
        }// Load a CAEN File which has a 4 line header and contains the spurrious bin 0's
        else if (type == 0){
            
            while(file >> trash) if (trash == "Tag") break;
            for (int i = 0; i < start; i++) file >> trash_time >> trash_bin >> trash_tag;
            while(file >> pointer -> time[counter] >> pointer -> bin[counter] >> pointer -> tag[counter]){
                pointer ->tag[counter] = 0;
                counter++;
                std::cout << counter << "\r";
                if (stop*3.6*pow(10, 12) <= pointer->time[counter - 1]) break;
            }
        }// Load an ascii file with three columns in the order time, bin, tag and reset tag to 0
        else if (type == 1){
            while(file >> trash) if (trash == "Tag") break;
            for (int i = 0; i < start; i++) file >> trash_time >> trash_bin >> trash_tag;
            while(file >> pointer -> time[counter] >> pointer -> bin[counter] >> pointer -> tag[counter]){
                counter++;
                std::cout << counter << "\r";
                if (stop*3.6*pow(10, 12) <= pointer->time[counter - 1]) break;
            }
        }// Load an ascii file with three columns in the order time, bin, tag
        else if (type == 2){
            while(file >> trash) if (trash == "Bin") break;
            for (int i = 0; i < start; i++) file >> trash_time >> trash_bin;
            while(file >> pointer -> time[counter] >> pointer -> bin[counter]){
                counter++;
                std::cout << counter << "\r";
                if (stop*3.6*pow(10, 12) <= pointer->time[counter - 1]) break;
            }
        }// Load an ascii file with two columns in the order time, bin
        std::cout << "\n";
        pointer->last = counter;
        file.close();
        exit = 1;
        if (czeros) printf("Number of skipped zeros %d\n", czeros);
    }
    return exit;
}

/* **************************************************************************************************** */
int file_load(info* pointer, int type = 0){
/* **************************************************************************************************** */
    /*
     This function will save data contained in an info struct with the current tagging upto the last
     data entry marked by last.
     INPUTS:
     pointer    => Location to load Data, should already contain filename and path in pointer->filename
     type       => Indicates the information to be saved
     value      => bool, true -> save signal, false -> save laser
     
     */
    // Initialize Function Variables:
    ifstream file;
    file.open(pointer->filename);
    int exit = 0;
    int counter = 0;
    int calibrow = (fitmodel ? 5: 7);
    float value = 0;
    string trash;
    // End
    
    if(file.fail()) std::cout << "\nUnable to open file...\n" << pointer->filename << "\n";
    else{
        std::cout << "\nLoaded data from:\t" << pointer->filename << "\n";
        pointer->current = 0;
        std::cout << "\nTotal events:\n";
        if (type == 0) {
            while(file >> trash) if (trash == "Tag") break;
            while(file >> pointer->time[pointer->current] >> pointer->bin[pointer->current] >> pointer->tag[pointer->current]){
                std::cout << ++pointer->current << "\r";
            }
        }
        else if (type == 2){
            while(file >> trash) if (trash == "ID") break;
            while(file >> pointer->time[pointer->current] >> pointer->bin[pointer->current] >> pointer->tag[pointer->current]
                  >> pointer->recoil[pointer->current] >> pointer->id[pointer->current]){
               std::cout << ++pointer->current << "\r";            }
        } // Loads energy calibrated file with IDs...
        else if (type == 3){
            pointer->numbin = 0;
            while(file >> trash) if (trash == "Partitions:") {file >> pointer->endslice; break;} // Clears header...
            for (int i = 0; i < pointer->endslice + 1; i++) file >> trash; // Clears out the column titles..
            counter = 0;
            pointer->current = 0;
            file >> pointer->xaxis[pointer->current];
            while(file >> pointer->hist[counter][pointer->current]){
                counter++;
                if (counter == pointer->endslice) {
                    pointer->current++;
                    pointer->numbin++;
                    file >> pointer->xaxis[pointer->current];
                    counter = 0;
                }
            }
            pointer->cslice = pointer->endslice - 1;
            printf("Loaded window file with %d slices\n", pointer->cslice + 1);
        } // Loads channel window histograms file...
        /*/
        else if (type == -3){
            for (int i = 0; i < pointer->endslice + 2; i++) {
                file >> trash;// Clears out the column titles...
                std::cout << trash << "\n";
            }
            counter = 0;
            pointer->current = 0;
            file >> pointer->xaxis[pointer->current];
            while(file >> pointer->hist[counter][pointer->current]){
                counter++;
                if (counter == pointer->endslice) {
                    pointer->current++;
                    file >> pointer->xaxis[pointer->current];
                    counter = 0;
                }
            }
            pointer->cslice = pointer->endslice - 1;
        } // Loads energy window histograms file.../*/
        else if (type == 4){
            if (fitmodel){ while(file >> trash) if(trash == "dOffset" || trash == "004_dOffset") {file >> trash; break;}} // Clears out the column titles...
            else {while(file >> trash) if(trash == "006_dQuad") {file >> trash; break;}}// Clears out the column titles...
            counter = 0;
            pointer->current = 0;
            std::cout << "Start readout\n";
            while(file >> pointer->calib[pointer->current][counter%calibrow]){
                if (counter == 0 && pointer->current == 0 &&
                    pointer->calib[pointer->current][counter%calibrow] != 0){
                    pointer->current = pointer->calib[pointer->current][counter%calibrow];
                    pointer->calib[pointer->current][counter%calibrow] = pointer->calib[0][counter%calibrow];
                    pointer->calib[0][counter%calibrow] = 0;
                }
                counter++;
                std::cout << counter << "\r";
                if(!(counter%calibrow)) pointer->current++;
            }
            pointer->cslice = pointer->current - 1;
            printf("Loaded calibration file with %d slices\n", pointer->cslice + 1);
        } // Loads calibration file...
        else if (type == 5){
            //{fit_number, energy, peak, dpeak, sigma, dsigma, area, darea}
            for (int j = 0; j < NUM_PEAKS*NUM_COLUMN; j++) pointer->param[j][0] = -1;
            while(file >> trash) if (trash == "dA" || trash == "007_dA") {file >> trash; break;} // Clears out the column titles...
            counter = 0;
            pointer->current = 0;
            pointer->cslice = 0;
            while(file >> pointer->param[NUM_PEAKS*pointer->cslice + pointer->current][counter%8]){
                if(!(counter%8) &&
                   pointer->cslice != pointer->param[NUM_PEAKS*pointer->cslice + pointer->current][0]) {
                    pointer->cslice = pointer->param[NUM_PEAKS*pointer->cslice + pointer->current][0];
                    pointer->param[NUM_PEAKS*(pointer->cslice - 1) + pointer->current][0] = -1;
                    pointer->current = 0; // With the ++ later this makes this value 0 for the next read...
                    pointer->param[NUM_PEAKS*pointer->cslice][0] = pointer->cslice;
                }
                counter++;
                if(!(counter%8)) pointer->current++;
            }
            printf("Loaded parameter file with %d slices\n", pointer->cslice + 1);
        } // Loads parameter file...
        file.close();
        pointer->last = pointer->current;
        pointer->current = 0;
        pointer->endslice = pointer->cslice + 1;
        pointer->cslice = 0;
        exit = 1;
    }
    return exit;
}

/* **************************************************************************************************** */
int file_save(info* pointer, char* new_name, int type = 0, int start = 0, int stop = N_Events){
/* **************************************************************************************************** */
    /*
     This function will save data contained in an info struct with the current tagging upto the last
     data entry marked by last.
     INPUTS:
     pointer    => Data to be saved
     new_name   => New file name for saved data
     type       => Indicates the information to be saved
     value      => bool, true -> save signal, false -> save laser
     
     */
    // Initialize Function Variables:
    ofstream file;
    file.open(new_name);
    int exit = 0;
    // End
    
    if(file.fail()) std::cout << "\nUnable to open save file...\n" << new_name << "\n";
    else{
        //file << "Time (ns)\tBin\tTag\n";
        if (type == 0) {
            file << "Time (ns)\tBin\tTag\n";
            std::cout << "\nTotal events:\n";
            for (int counter = start; counter < pointer->last; counter++) {
                file << pointer->time[counter] << "\t"  << pointer->bin[counter] << "\t"  << pointer->tag[counter] << "\n";
                std::cout << counter + 1 << "\r";
                if (counter >= stop) break;
            }
        } // Saves uncalibrated ascii file, old
        else if (type == 1){
            file << "Time (ns)\tBin\tVeto\tHigh\tRecoil\n";
            std::cout << "\nTotal events:\n";
            for (int counter = start; counter < pointer->last; counter++) {
                if (!pointer->tag[counter]){
                    file << pointer->time[counter] << "\t";
                    file << pointer->bin[counter] << "\t";
                    file << pointer->veto[counter] << "\t";
                    file << pointer->high[counter] << "\t";
                    file << pointer->recoil[counter] << "\n";
                    std::cout << counter + 1 << "\r";
                    if (counter >= stop) break;
                }
            }
        } // Saves calibrated ascii file, old
        else if (type == 2){
            file << "Time (ns)\tBin\tTag\tRecoil\tID\n";
            std::cout << "\nTotal events:\n";
            for (int counter = start; counter < pointer->last; counter++) {
                file << pointer->time[counter] << "\t";
                file << pointer->bin[counter] << "\t";
                file << pointer->tag[counter] << "\t";
                file << pointer->recoil[counter] << "\t";
                file << pointer->id[counter] << "\n";
                std::cout << counter + 1 << "\r";
                if (counter >= stop) break;
            }
        } // Saves an uncalibrated ascii file
        else if (type == 3){
            file << "Time(ns)\tBin\n";
            std::cout << "\nTotal events:\n";
            for (int counter = start; counter < pointer->last; counter++) {
                file << pointer->time[counter] << "\t";
                file << pointer->bin[counter] << "\n";
                std::cout << counter + 1 << "\r";
                if (counter >= stop) break;
            }
        } // In Progress save all events
        file.close();
        std::cout << "\nSuccessfully saved data to:\t" << new_name << "\n";
        exit = 1;
    }
    return exit;
}

/* **************************************************************************************************** */
int fit_save(info *pointer, bool append = false, short sopt = 0, bool domain = false,
             bool verbose = true){
/* **************************************************************************************************** */
    /*
        This function will save the calibration or fit data stored in the master calibration or fit 
        array to disk.
        Parameters:
        pointer ->  Master file containing all the fit parameters and calibration values
        append  ->  bool, indicates whether this save is appending to the end of the file, 
                    default false which over writes the prior file completely.
        sopt    ->  short, options for saving,
                    0 -> save both parameters and calibrations
                    1 -> save parameters only
                    2 -> save calibration only
        domain  ->  bool, indicates the units for the fit as either eV or ch, default false for ch
        verbose ->  Print out the location and name of file being saved, default true
     
     */
    // Initialize Function Variables
    ofstream par, cal;
    int exit = 0;
    int paramrow = 8;
    int calibrow = (fitmodel ? 5 : 7);
    char str[500];
    char word[500];
    char minstr[10];
    char indexstr[10];
    // End
    
    snprintf(word, sizeof(word), "");
    snprintf(indexstr, sizeof(indexstr), "0");
    if(pointer->indexvalue < 0){
        snprintf(word, sizeof(word), "_%d-%d", -pointer->indexvalue, -pointer->indexvalue + 21);
        if (sopt == 2) pointer->indexvalue = 3;
        if (sopt == 1) pointer->indexvalue = 20;
        
    }
    
    if (sopt == 0 || sopt == 1){
        if (10 <= pointer->indexvalue) snprintf(indexstr, sizeof(indexstr), "");
        snprintf(str, sizeof(str), "%s%s_Ch_%d_%s%d_Parameters%s.txt", pointer->path, pointer->name,
                 pointer->channel, indexstr, pointer->indexvalue, word);
        if (append) par.open(str, ios::out| ios::app);
        else par.open(str);
        if (verbose) printf("\nSaving parameters to %s\n", str);
    }
    if (sopt == 0 || sopt == 2){
        if (sopt == 0 &&  9 <= pointer->indexvalue) snprintf(indexstr, sizeof(indexstr), "");
        if (sopt == 2 && 10 <= pointer->indexvalue) snprintf(indexstr, sizeof(indexstr), "");
        snprintf(str, sizeof(str), "%s%s_Ch_%d_%s%d_Calibration%s.txt", pointer->path, pointer->name,
                 pointer->channel, indexstr, pointer->indexvalue + (sopt == 2 ? 0 : 1), word);
        if (append) cal.open(str, ios::out| ios::app);
        else cal.open(str);
        if (verbose) printf("\nSaving calibrations to %s\n", str);
    }
    
    /*/
        The general outline for saving a file is read each column from a row then move to the next row.
     
        Parameters stored as:
        {{fit_number, energy, peak, dpeak, sigma, dsigma, area, darea}_0, ...}
     
        Calibration stored as:
        {{fit_number, slope, dslope, offset, doffset}_0, ... }
     /*/
    
    
    if (par.fail()||cal.fail()) {
        std::cout << "Could not open save location for the fitted data...\n";
        exit = 0;
    }
    else if (append){
        if (sopt == 0 || sopt == 1){
            for (int i = pointer->cslice*NUM_PEAKS; i < (pointer->cslice + 1)*NUM_PEAKS; i++) {
                for (int j = 0; j < paramrow; j++) {
                    if (pointer->param[i][0] != -1) par << pointer->param[i][j] << (j == paramrow - 1 ? "\n" : "\t");
                }
            }
            par.close();
        }
        if (sopt == 0 || sopt == 2){
            for (int j = 0; j < calibrow; j++) {
                if (pointer->calib[pointer->cslice][0] != -1) {
                    cal << pointer->calib[pointer->cslice][j] << (j == calibrow - 1 ? "\n" : "\t");
                }
            }
            cal.close();
        }
        exit = 1;
    }
    else {
        if (sopt == 0 || sopt == 1){
            snprintf(minstr, sizeof(minstr), "%s", (domain ? "eV" : "Ch"));
            snprintf(str, sizeof(str), "000_Time Slice\t001_Energy\t");
            snprintf(str, sizeof(str), "%s002_Centroid (%s)\t003_dC (%s)\t", str, minstr, minstr);
            snprintf(str, sizeof(str), "%s004_Sigma (%s)\t005_dS (%s)\t", str, minstr, minstr);
            snprintf(str, sizeof(str), "%s006_Area (cts/%s)\t007_dA (cts/%s)\n", str, minstr, minstr);
            par << str;
            for (int i = 0; i < NUM_COLUMN*NUM_PEAKS; i++) {
                for (int j = 0; j < paramrow; j++) {
                    if (pointer->param[i][0] != -1) par << pointer->param[i][j] << (j == paramrow - 1 ? "\n" : "\t");
                }
            }
            par.close();
        }
        if (sopt == 0 || sopt == 2){
            cal << "000_Time Slice\t001_Slope (eV/Ch)\t002_dSlope (eV/Ch)\t003_Offset (eV)\t004_dOffset (eV)";
            if(fitmodel) cal << "\n";
            else cal << "\t005_Quad (eV/Ch^2)\t006_dQuad (eV/Ch^2)\n";
            for (int i = 0; i < pointer->endslice; i++) {
                for (int j = 0; j < calibrow; j++) {
                    if (pointer->calib[i][0] != -1) cal << pointer->calib[i][j] << (j == calibrow - 1 ? "\n" : "\t");
                }
            }
            cal.close();
        }
        exit = 1;
    }
    return exit;
}

/* **************************************************************************************************** */
int hist_save(info *pointer, int histrow = tBin, int histcol = 0, bool verbose = true) {
/* **************************************************************************************************** */
    /*/
        This function saves the histograms stored in the pointer->hist arrays. Does nothing to
        histograms stored in pointer->histogram.
    /*/
    // Initialize Function Variables
    ofstream file;
    int exit = 0;
    
    // End
    sprintf(pointer->filename, "%s%s_Histograms_Channel_%d.txt", pointer->path,
            pointer->name, pointer->channel);
    file.open(pointer->filename);
    if(file.fail()) std::cout << "\nUnable to open save file...\n" << pointer->filename << "\n";
    else{
        if (verbose) std::cout << "Saving to: " << pointer->filename << "\n";
        file << "Energy (eV)\t";
        for (int i = 0; i < histcol; i++) file << "Histogram_" << i << (i == histcol - 1 ? "\n" : "\t");
        for (int i = 0; i < histrow; i++){
            file << pointer->xaxis[i] << "\t";
            for (int j = 0; j < histcol; j++) file << pointer->hist[j][i] << (j == histcol - 1 ? "\n" : "\t");
        }
        file.close();
        exit = 1;
    }
    return exit;
}

/* **************************************************************************************************** */
int xia2ascii(char folder[500], char file[500], char ext[500], char chan[500]){
/* **************************************************************************************************** */
    /*
     This function is intended to convert the XIA Binary file into an ascii file with the bin, time and
     tag information in three columns.
     */
    // Initialize Function Variables
    info *detector[N_detectors];
    int exit = 0;
    for (int i = 0; i < N_detectors; i++) {
        detector[i] = new info;
        detector[i]->current = 0;
        detector[i]->last = 0;
        detector[i]->channel = i;
    }
    // End
    
    std::cout << "Loading data from file\n";
    sprintf(detector[0]->filename, "%s%s%s", folder, file, ext);
    exit = xiaload(detector[0]->filename, detector, 0, 24);
    for(int i = 0; i < N_detectors; i++) {
        if (exit == 1 && detector[i]->last > 0) {
            std::cout << "Detector number:\t" << i << "\tLast:\t" << detector[i]->last << "\n";
            //sprintf(detector[i]->filename, "%s%s%s%s%s%d%s", folder, file, "/", file, chan, i, ".txt");
            snprintf(detector[i]->filename, sizeof(detector[i]->filename), "%s%s/%s_Ch_%d_ASCII.txt", folder, file, file, detector[i]->channel);
            exit = file_save(detector[i], detector[i]->filename);
        }
    }
    return exit;
}

/* **************************************************************************************************** */
int trigger_period(info* pointer){
/* **************************************************************************************************** */
    /*
     This function determines an average period of the laser repetition rate from the trigger file
     */
    return int(1.0*(pointer->time[pointer->last - 1] - pointer->time[0])/(pointer->last - 1));
}

/* **************************************************************************************************** */
double trigger_offset(info* pointer, info* trig, int period = 0){
/* **************************************************************************************************** */
    /*
     This function determines the average time offset between the trigger channel and the data channel
     by subtracting all the data events that occur in the interval of the trigger and looking for a
     peak which indicates a periodic event with some offset from the trigger event.
    */
    // Initialize Function Variables
    int counter = 0;
    int max = 0;
    int exit = 0;
    bool trouble = false;
    double output = 0;
    long long next = 0;
    TF1 *fun;
    TH1I *hist;
    //TSpectrum *peak_off = new TSpectrum(5);
    // End
    
    std::cout << "\nCalculating the trigger-data offset...\n";
    if (period == 0) period = trigger_period(trig);
    if(trouble) printf("Assuming a period of %d from values %lld and %lld\n", period, trig->time[0], trig->time[trig->last - 1]);
    pointer->current = 0;
    trig->current = 0;
    hist = new TH1I("Hist", "T-Difference", period, 0, 2*period);
    hist->Reset();
    
    while (pointer->time[pointer->current] < trig->time[trig->current]) pointer->current++; // advance until pointer is ahead of trigger
    
    if(trouble) printf("Tagging events, starting at pointer->current of %d out of %d\n", pointer->current, pointer->last);
    //You want to scan through the pointer values and adjust the trigger values only when necessary to bound the pointer values.
    for (int i = pointer->current; i < pointer->last; i++){
        while (0 < trig->current && pointer->time[i] < trig->time[trig->current]) trig->current--;
        while (trig->current < trig->last && trig->time[trig->current + 1] < pointer->time[i]) trig->current++;
        if (trig->time[trig->current] <= pointer->time[i] && pointer->time[i] < trig->time[trig->current + 1]){
            hist->Fill(pointer->time[i] - trig->time[trig->current]);
        }
        if(i < pointer->last - 1 && trig->time[trig->last - 1] < pointer->time[i + 1]) break;
    }
    max = hist->GetBinCenter(hist->GetMaximumBin());
    fun = new TF1("fun", "gaus", max - 2500, max + 2500);
    hist->Fit(fun, "OQR");
    output = fun->GetParameter(1);
    if(!trouble) hist->Delete();
    if(trouble) ihost = hist;
    return output;
}

/* **************************************************************************************************** */
int trigger_tag(info* pointer, info* trig,  int sig = 4000, double offset = -1){
/* **************************************************************************************************** */
    /*/
     This function will tag the laser events in the data file based on the laser trigger file
     pointer    =>  Detector channel
     trig       =>  Laser trigger channel
     sig        =>  (int) Maximum time difference between events to be considered coincidence in us.
     
     Needs Improvement...
    /*/
    // Initialize Function Variables:
    int counter = 0;
    //int trouble = 0;
    //int marker = 0;
    
    // End
    offset = (offset < 0 ? trigger_offset(pointer, trig) : offset);
    
    /*/
     Time match the two files so that the current data time stamp is after the current trigger time
     stamp. 
     Note: 
     -100       => Skipped element in the tag column
     -1         => Event in coincidence with the laser trigger file to within sigma and the offset.
     /*/
    pointer->current = 0;
    trig->current = 0;
    for(int i = 0; i < pointer->last; i++) pointer->tag[i] = 0;
    for(int i = 0; i < trig->last; i++) trig->tag[i] = 0;
    
    printf("Signal first time point: %lld\n", pointer->time[pointer->current]);
    printf("Tag first time point: %lld\n", trig->time[trig->current]);
    
    while (pointer->time[pointer->current] < trig->time[trig->current])  {
        pointer->tag[pointer->current] = -100;
        pointer->current++;
    }
    while (trig->time[trig->current + 1] < pointer->time[pointer->current])  trig->current++;
    printf("Signal starting time: %lld\n", pointer->time[pointer->current]);
    printf("Tag starting time: %lld\n", trig->time[trig->current]);
    
    printf("\nStart tagging points in %s\n", pointer->filename);
    printf("Using points in %s\n",trig->filename);
    printf("Using an offset value of %0.2f\n", offset);
    printf("Starting point in signal is %d out of %d\n", pointer->current, pointer->last);
    printf("Starting point in trigger is %d out of %d\n", trig->current, trig->last);
    printf("Tag laser events\n");
    for (int i = pointer->current; i < pointer->last; i++) {
        while (0 < trig->current &&
               pointer->time[i] < trig->time[trig->current]) trig->current--;
        while (trig->current < trig->last &&
               trig->time[trig->current + 1] < pointer->time[i]) trig->current++;
        
        if (trig->time[trig->current] + offset - sig < pointer->time[i] &&
            pointer->time[i] < trig->time[trig->current] + offset + sig &&
            pointer->tag[i] != -100) {
            pointer->tag[i] = 1;
            trig->tag[trig->current] = 1;
            counter++;
        } // Tag laser event
        else if (pointer->tag[i] != -100) pointer -> tag[i] = 0; // Tag signal event
        if (trig->time[trig->last - 1] < pointer->time[i]){
            pointer->current = i;
            break;
        }
        else if (i == pointer->last - 1) pointer->current = pointer->last;
    }
    printf("\nSkipping points %d through %d\n", pointer->current, pointer->last);
    while (pointer->current < pointer->last) {
        pointer->tag[pointer->current] = -100;
        pointer->current++;
    }
    printf("Total number of tagged events %d out of %d data points and %d laser tags\n",
           counter, pointer->last, trig->last);
    return (0 < counter ? int(offset) : 0);
}

/* **************************************************************************************************** */
template<class type>
void sort_list(type *unsorted, int num){
/* **************************************************************************************************** */
    /*
     This function will sort an array into ascending order lowest value first in the event of ties the
     first entry found remains first.
     */
    // Initalize Function Variables:
    type sorted[NUM_PEAKS];
    int sort_val[NUM_PEAKS];
    
    for (int i = 0; i < NUM_PEAKS; i++) sort_val[i] = 0;
    // End
    
    for (int i = 0; i < num; i++) for (int j = i + 1; j < num; j++) {
        unsorted[i] <= unsorted[j] ? sort_val[j]++: sort_val[i]++;
    }
    for (int i = 0; i < num; i++) sorted[sort_val[i]] = unsorted[i];
    for (int i = 0; i < num; i++) unsorted[i] = sorted[i];
}

/* **************************************************************************************************** */
template<class type>
type sum_array(type *column, int num){
/* **************************************************************************************************** */
    /*
     This function will sum an array
     */
    // Initalize Function Variables:
    type val = 0;
    // End
    for (int i = 0; i < num; i++) val += column[i];
    return val;
}

/* **************************************************************************************************** */
template<class type>
type findmax(type *column, int num){
/* **************************************************************************************************** */
    /*
     This function will locate the maximum in a list.
     */
    // Initalize Function Variables:
    type val = column[0];
    // End
    for (int i = 1; i < num; i++) val = (val < column[i] ? column[i] : val);
    return val;
}

/* **************************************************************************************************** */
template<class type>
type findmin(type *column, int num){
/* **************************************************************************************************** */
    /*
     This function will locate the maximum in a list.
     */
    // Initalize Function Variables:
    type val = column[0];
    // End
    for (int i = 1; i < num; i++) val = (column[i] < val ? column[i] : val);
    return val;
}

/* **************************************************************************************************** */
template<class type1, class type2>
int lesser(type1 *column, type2 find, int num){
/* **************************************************************************************************** */
    /*
     This function will locate where the last incident (or passed) of a value occurs in a list
     */
    // Initalize Function Variables:
    int val = 0;
    // End
    for (int i = 0; i < num; i++) {
        //std::cout << find << " < " << column[i] << "\n"; //Trouble Shooting
        if (column[i]< find ) val = i;
        else break;
    }
    return val;
}

/* **************************************************************************************************** */
template<class type>
int where(type *column, type find, int num){
    /* **************************************************************************************************** */
    /*
     This function will locate where the last incident of a value occurs in a list
     */
    // Initalize Function Variables:
    int val = 0;
    // End
    for (int i = 0; i < num; i++) val = (find == column[i] ? i : val);
    return val;
}

/* **************************************************************************************************** */
template<class type1, class type2>
bool element(type1 *column, type2 find, int num){
/* **************************************************************************************************** */
    /*
     This function will determine if find is in column
     */
    // Initalize Function Variables:
    bool val = false;
    // End
    for (int i = 0; i < num; i++) if (find == column[i]) {val = true; break;}
    return val;
}

/* **************************************************************************************************** */
template<class type1, class type2>
bool dot(type1 *column, type2 *find, int num){
    /* **************************************************************************************************** */
    /*
     This function will multiply the first entry by the second and return it in the first...
     */
    for (int i = 0; i < num; i++) column[i] *= find[i];
    return true;
}

/* **************************************************************************************************** */
void idkey(short *line, int value, int num){
/* **************************************************************************************************** */
    /*
     Takes the ID generated in marker and breaks it down to the individual a - i components
     ID Values:
     _ _ _ _ _ _ _ _ _
     A B C D E F G H I
     A: Event type
     -1     =>  Laser, tagged from the beginning
     0      =>  Anticoincidence, all values
     i      =>  Coincidence, i detectors
     
     B - E: Detector Bin Value (Channel 0, 2, 4, 6)
     0      =>  0
     1      =>  0 < # < 7000
     2      =>  7000 < ##
     
     F - I: Detector Timing (Channel 0, 2, 4, 6)
     0      => Not Applicable
     i      => i-th timing
     */
    // Initialize Function Variables
    int storage = value;
    //
    for (int k = 0; k < num; k++) {
        line[k] = storage/int(pow(10, num - k - 1));
        storage -= line[k]*int(pow(10, num - k - 1));
    }
}

/* **************************************************************************************************** */
template<class type>
int removept(type array, int line, int start, int finish){
/* **************************************************************************************************** */
    /*/
        This function is intended to remove a data point in pointer given by line and collapse the 
        array so that the removed point no longer shows up on the list
        array       => Any 1-D array
        line        => Parameter line to be removed
        start       => Start location in array
        finish      => End location in array
        If the array does not contain period entries the start = 0 and finish is length of array
    /*/
    for (int i = start; i < finish; i++) {
        if (i - start < line) continue;
        else if (i < finish - 1) array[i] = array[i + 1];
        else if (i == finish - 1) array[i] = 0;
    }
    return 1;
}

/* **************************************************************************************************** */
int removept(info* pointer, int line, int start, int finish){
/* **************************************************************************************************** */
    /*/
        This function is intended to remove a data point in pointer given by line and collapse the
        array so that the removed point no longer shows up on the list
        pointer     => Standard data structure
        line        => Parameter line to be removed
        start       => Start of slice
        finish      => start of next slice
    /*/
    for (int i = start; i < finish; i++) {
        if (i - start < line) continue;
        else if (i < finish - 1) for (int j = 0; j < 8; j++) pointer->param[i][j] = pointer->param[i + 1][j];
        else if (i == finish - 1) {
            pointer->param[i][0] = -1; // The -1 is my standard empty line marker
            for (int j = 1; j < 8; j++) pointer->param[i][j] = 0;
        }
    }
    return 1;
}

/* **************************************************************************************************** */
int addpt(info* pointer, double array[8], int line, int start, int finish){
/* **************************************************************************************************** */
    /*/
        This function is intended to add a data point in pointer given by line and collapse the
        array so that the removed point no longer shows up on the list
        pointer     => Standard data structure
        array       => This is the array to insert at line
        line        => Parameter line to be removed
        start       => Start of slice
        finish      => start of next slice
     /*/
    for (int i = start; i < finish + 1; i++) {
        if (i - start < line) continue;
        else if (i - start == line) for(int j = 0; j < 8; j++) pointer->param[i][j] = array[j];
        else if (line < i - start) for (int j = 0; j < 8; j++) pointer->param[i + 1][j] = pointer->param[i][j];
    }
    return 1;
}

/* **************************************************************************************************** */
double_t func_00(double_t *x, double_t *par){
/* **************************************************************************************************** */
    /*/
     This function models the calibration as N-distinct Gaussian peaks
    /*/
    double_t f = 0;
    double_t A_n = 0;
    double_t S_n = 0;
    double_t C_n = 0;
    for (int i = 0; i < par[0]; i++)  {
        A_n = par[3*i + 1];
        S_n = par[3*i + 3];
        C_n = par[3*i + 2];
        f += A_n/pow(2*TMath::Pi()*S_n*S_n, 0.5)*TMath::Gaus(*x, C_n, S_n);
    }
    return f;
}

/* **************************************************************************************************** */
double_t func_01(double_t *x, double_t *par){
/* **************************************************************************************************** */
    /*/
     This function models the calibration as one Gaussian at integer locations but with distinct sigma 
     and amplitudes set by poisson distribution. Input parameters are of the form:
     fit is conducted in energy space so [m] = eV/ch and [b] = eV
     p = {n-initial (inclusive), n-final (exclusive), m, b, A_0, S_0, A_1, S_1, ...}
    /*/
    double_t f = 0;
    double_t A_n = 0;
    double_t S_n = 0;
    double_t m = par[2];
    double_t b = par[3];
    double_t mu = par[4];
    double_t x0 = x[0];
    for (int i = 0; i < par[1] - par[0]; i++)  {
        S_n = par[2*i + 5];
        f += (pow(mu, i + 1)*TMath::Exp(-(i + 1))/TMath::Factorial(i+1))/pow(2*TMath::Pi()*S_n*S_n, 0.5)*TMath::Gaus(m*x0 + b, i + par[0], S_n);
    }
    return f;
}

/* **************************************************************************************************** */
double_t bkgfunc_00(double_t *x, double_t *par){
    /* **************************************************************************************************** */
    /*/
     This function is a superposition of an exponential decay and gaussian function.
     y(x, par) = par[0]/(2pi*par[2]**2)**0.5*exp(-(E-par[1])**2/(2*par[2]**2)) + par[3]/par[4]*exp(-E/par[4])
     
     par[3] represents the area of the exponential decay...
     Number of parameters: 5
     /*/
    double x0 = x[0];
    double_t f = par[3]/par[4]*TMath::Exp(-1.0*x0/par[4]);
    f += par[0]/pow(2*TMath::Pi()*par[2]*par[2], 0.5)*TMath::Gaus(*x, par[1], par[2]);
    return f;
}

/* **************************************************************************************************** */
double_t bkgfunc_01(double_t *x, double_t *par){
    /* **************************************************************************************************** */
    /*/
     This function is a superposition of two exponential decays and gaussian function.
     y(x, par) = par[0]/(2pi*par[2]**2)**0.5*exp(-(E-par[1])**2/(2*par[2]**2))
     + par[3]*exp(-E/par[4] + par[5]*exp(-E/par[6])
     
     par[3] represents the area of the exponential decay...
     par[5] represents the area of the exponential decay...
     
     Number of parameters: 7
     /*/
    double x0 = x[0];
    double_t f = par[3]/par[4]*TMath::Exp(- 1.0*x0/par[4]) + par[5]/par[6]*TMath::Exp(- 1.0*x0/par[6]);
    f += par[0]/pow(2*TMath::Pi()*par[2]*par[2], 0.5)*TMath::Gaus(*x, par[1], par[2]);
    return f;
}

/* **************************************************************************************************** */
double_t bkgfunc_02(double_t *x, double_t *par){
    /* **************************************************************************************************** */
    /*/
     This function is a superposition of a line and gaussian function.
     y(x, par) = par[0]/(2pi*par[2]**2)**0.5*exp(-(E-par[1])**2/(2*par[2]**2)) + par[4]*x + par[3]
     
     Number of parameters: 5
     /*/
    double x0 = x[0];
    double_t f = par[3] + par[4]*x0;
    f += par[0]/pow(2*TMath::Pi()*par[2]*par[2], 0.5)*TMath::Gaus(*x, par[1], par[2]);
    return f;
}

/* **************************************************************************************************** */
double_t bkgfunc_03(double_t *x, double_t *par){
    /* **************************************************************************************************** */
    /*/
     This function is a superposition of a quadratic and gaussian function.
     y(x, par) = par[0]/(2pi*par[2]**2)**0.5*exp(-(E-par[1])**2/(2*par[2]**2))
     + par[5]*x**2 + par[4]*x + par[3]
     
     Number of parameters: 6
     /*/
    double x0 = x[0];
    double_t f = par[3] + par[4]*x0 + par[5]*x0*x0;
    f += par[0]/pow(2*TMath::Pi()*par[2]*par[2], 0.5)*TMath::Gaus(*x, par[1], par[2]);
    return f;
}

/* **************************************************************************************************** */
double_t bkgfunc_04(double_t *x, double_t *par){
    /* **************************************************************************************************** */
    /*/
     This function is a superposition of an inverse and gaussian function.
     y(x, par) = par[0]/(2pi*par[2]**2)**0.5*exp(-(E-par[1])**2/(2*par[2]**2)) + par[3]/(x + par[4])
     
     Number of parameters: 5
     /*/
    double x0 = x[0];
    double_t f =  par[3]/(par[4] + x0);
    f += par[0]/pow(2*TMath::Pi()*par[2]*par[2], 0.5)*TMath::Gaus(*x, par[1], par[2]);
    return f;
}

/* **************************************************************************************************** */
template<class type1, class type2>
double_t rchi(type1 *xo, type2 *xe, int num, int dof){
    /* **************************************************************************************************** */
    /*/
     This is the reduced chi square function:
     X2 = sum(xo_i - xe_i)
     Input Parameters:
     xo     => Observed data
     xe     => Expected value (evaluated model)
     num    => Number of points in arrray
     dof    => Number of parameters in model
     
     Output:
     f      => Reduced chi square value
     /*/
    double f = 0;
    for(int i = 0; i < num; i++) f += (xo[i] - xe[i])*(xo[i] - xe[i])/(num - dof - 1.0);
    return f;
}

/* **************************************************************************************************** */
template<class type1, class type2>
void stat(type1 *output, type2 *input, type2 *weights, int num){
/* **************************************************************************************************** */
    /*/
     This will calculate the weighted average and weighted standard deviation for the input array and
     return the values inside of the output array.
     
     Uses the equations for weighted mean and standard deviation based off the NIST website paper
     /*/
    
    output[0] = 0;
    output[1] = 0;
    for (int i = 0; i < num; i++) output[0] += weights[i]*input[i];
    output[0] = output[0]*1.0/sum_array(weights, num);
    
    for (int i = 0; i < num; i++) output[1] += num*weights[i]*pow(input[i] - output[0], 2)/(num - 1);
    output[1] = output[1]/sum_array(weights, num);
    output[1] = pow(output[1], 0.5);
}

/* **************************************************************************************************** */
template<class type1, class type2>
void derivative(type1 *output, type2 *input, type2 *weights, int num, int der = 1, double space = 1){
/* **************************************************************************************************** */
    /*/
     This will calculate the n-th order derivative of the input array
    
     Uses the equations for O(h^4) central difference method for numerical calculations,
     1 -> 4 are from the book Numerical Methods Using Matlab and may be found online as well,
     5+     are derived using the methods to get 1 -> 4 and those were recreated to confirm conformity
     
     Coefficients always read in order from after center to before
     ..., C_2, C_1, C_0, C_-1, C_-2, ...
     starting with the lead term
     
     conum refers to the number of coefficient to the left excluding 0.
     
     This function assumes that the weights will be given in the form of 1/error
    /*/
    // Initialize Variables]
    double coef[20];
    int conum = 0;
    bool track = false;
    int denom = 1;
    // End
    
    for (int i = 0; i < 20; i++) coef[i] = 0;
    for (int i = 0; i < 2; i++) output[i] = 0;
    if (der == 1){
        conum = 2;
        denom = 12;
        coef[0] = -1;   // C_2
        coef[1] =  8;   // C_1
        coef[2] =  0;   // C_0
    }
    if (der == 2){
        conum = 2;
        denom = 12;
        coef[0] =  -1;
        coef[1] =  16;
        coef[2] = -30;
    }
    if (der == 3){
        conum = 3;
        denom = 8;
        coef[0] = -1;   // C_3
        coef[1] =  8;   // C_2
        coef[2] = 13;   // C_1
        coef[3] =  0;   // C_0
    }
    if (der == 4){
        conum = 3;
        denom = 6;
        coef[0] =  -1;
        coef[1] =  12;
        coef[2] = -39;
        coef[3] =  56;
    }
    if (der == 5){
        conum = 4;
        denom = 6;
        coef[0] =  -1;   // C_4
        coef[1] =   9;   // C_3
        coef[2] = -26;   // C_2
        coef[3] =  29;   // C_1
        coef[4] =   0;   // C_0
    }
    if (der == 6){
        conum = 4;
        denom = 4;
        coef[0] =   -1;
        coef[1] =   12;
        coef[2] =  -52;
        coef[3] =  126;
        coef[4] = -150;
    }
    denom *= pow(space, der);
    if (coef[conum] == 0) track = true;
    for (int i = 0; i < conum; i++) coef[conum + i + 1] = pow(-1, (track ? 1 : 0))*coef[conum - i - 1];
    for (int i = 0; i < num; i++) output[0] += coef[i]*input[i]/denom;
    for (int i = 0; i < num; i++) output[1] += pow(coef[i]/weights[i], 2);
    output[1] = pow(output[1], 0.5)/denom;
}

/* **************************************************************************************************** */
template<class type1, class type2>
void nl_offset(type1 *output, type2 *input, type2 *weights, int num, int der_order = 2, double space = 1){
    /* **************************************************************************************************** */
    /*/
     This will calculate the n-th order derivative non-linear offset.
     
     The numerical derivatives stored in derarr are 
     derarr = [D1, dD1, D2, dD2, D3, dD3, ...]
     
     /*/
    // Initialize Variables]
    double store1 = 0;
    double store2 = 0;
    double derarr[50];
    double derval[2];
    // End
    
    for (int i = 0; i < 2; i++) output[i] = 0;
    for (int i = 0; i < 50; i++) derarr[i] = 0;
    for (int i = 1; i < der_order; i++){
        derivative(derval, input, weights, num, i, space);
        derarr[2*i - 2] = derval[0];
        derarr[2*i - 1] = derval[1];
    }
    store1 = 0;
    for(int n = 2; n < der_order + 1; n++){
        store2 = 0;
        for(int m = n; m < der_order + 1; m++) store2 += pow(-space, m - n)*derarr[2*m - 2]/TMath::Factorial(m - n);
        store1 += (n - 1)*store2/TMath::Factorial(n);
    }
    output[0] = store1/derarr[0];
    
}

/* **************************************************************************************************** */
int background(info *pointer, TH1D *hist, double pvalue[7], int type = 0, float start = 0,
               float end = 0, int num = 0, bool trouble = false){
/* **************************************************************************************************** */
    /*/
     This function will fit the histogram with a gaussian function and background function:
     The fit functions are:
     bkgfunc_00    =>  Gaussian + Exponential Decay
     bkgfunc_01    =>  Gaussian + Exponential Decay + Exponential Decay
     bkgfunc_02    =>  Gaussian + Line
     bkgfunc_03    =>  Gaussian + Quadratic
     bkgfunc_04    =>  Gaussian + Inverse
     Input Parameters:
     hist   => pointer to histogram containing data to be fitted
     type   => indicates the function to be used
     param  => array of initial guess values for the data and specified background not all values used
     start  => start point of energy range for fit (inclusive)
     end    => end point of energy range for fit (inclusive)
     num    => number of bins in histogram (not needed?)
     
     Note that whatever the errorbars are on hist is what will be used this function does not
     alter the error values or recalculate them.
     /*/
    // Initialize Variables
    int counts[8192];
    int left = 0;
    int right = 0;
    int mid = 0;
    int pnum = 0;
    int exit = 0;
    char str[500];
    char word[500];
    float energy[8192];
    double *fitparam;
    const double *fiterror;
    TF1 *fun = new TF1("Single_exp", bkgfunc_00, start, end, pnum);
    ofstream outfile;
    // End
    
    if (!start) start = 76.5 - 5;
    if (!end) end = 76.5 + 5;
    for (int i = 0; i < num; i++) {
        energy[i] = hist->GetXaxis()->GetBinCenter(i + 1);
        counts[i] = hist->GetBinContent(i + 1);
    }
    left = lesser(energy, start, num);
    right = lesser(energy, end, num);
    
    if (trouble){
        printf("Left Bound: %d\tRight Bound: %d\n", left, right);
        printf("Left Energy: %.2f\tRight Energy: %.2f\n", energy[left], energy[right]);
    }
    
    // Select background
    if (0 <= type && type <= 4) {
        fun->Delete();
        if (type == 0) {
            if (!pvalue[4]) {
                pvalue[4] = TMath::Log(1.0*counts[right]/counts[left])/(energy[right] - energy[left]);
            }
            if (!pvalue[3]) pvalue[3] = counts[left]*TMath::Exp(energy[left]/pvalue[4]);
            
            pnum = 5;
            fun = new TF1("Single_exp", bkgfunc_00, start, end, pnum);
            fun->SetParameter(3, pvalue[3]);
            fun->SetParLimits(3, 0, pow(10, 9));
            fun->SetParameter(4, pvalue[4]);
            fun->SetParLimits(4, 0, pow(10, 4));
        } // Single exp
        else if (type == 1) {
            mid = lesser(energy, 73, num);
            if (!pvalue[4]) pvalue[4] = TMath::Log(1.0*counts[mid]/counts[left])/(energy[mid] - energy[left]);
            if (!pvalue[3]) pvalue[3] = counts[left]*TMath::Exp(energy[left]/pvalue[4]);
            
            mid = lesser(energy, 80, num);
            if (!pvalue[6]) pvalue[6] = TMath::Log(1.0*counts[right]/counts[mid])/(energy[right] - energy[mid]);
            if (!pvalue[5]) pvalue[5] = counts[right]*TMath::Exp(energy[right]/pvalue[6]);
            
            pnum = 7;
            fun = new TF1("Double_exp", bkgfunc_01, start, end, pnum);
            fun->SetParameter(3, pvalue[3]);
            fun->SetParLimits(3, 0, pow(10, 7));
            fun->SetParameter(4, pvalue[4]);
            fun->SetParLimits(4, 1, pow(10, 3));
            fun->SetParameter(5, pvalue[5]);
            fun->SetParLimits(5, -10, pow(10, 10)); // Should allow value to approach 0
            fun->SetParameter(6, pvalue[6]);
            fun->SetParLimits(6, 1, pow(10, 4));
        } // Double exp
        else if (type == 2) {
            if (!pvalue[4]) pvalue[4] = (counts[right] - counts[left])/(energy[right] - energy[left]);
            if (!pvalue[3]) pvalue[3] = counts[left] - energy[left]*pvalue[4];
            
            pnum = 5;
            fun = new TF1("Linear", bkgfunc_02, start, end, pnum);
            fun->SetParameter(3, pvalue[3]);
            fun->SetParameter(4, pvalue[4]);
        } // Linear
        else if (type == 3) {
            mid = lesser(energy, 73, num);
            if (!pvalue[5]) {
                pvalue[5] = (counts[mid] - counts[left])/(energy[mid] - energy[left]);
                pvalue[5] -= (counts[right] - counts[mid])/(energy[right] - energy[mid]);
                pvalue[5] = pvalue[5]*1.0/(energy[left] - energy[right]);
                
            }
            if (!pvalue[4]) {
                pvalue[4] = (counts[mid] - counts[left])*1.0/(energy[mid] - energy[left]);
                pvalue[4] -= pvalue[5]*(energy[mid] - energy[left]);
            }
            if (!pvalue[3]) {
                pvalue[3] = counts[left] - pvalue[4]*energy[left] - pvalue[5]*energy[left]*energy[left];
            }
            //pvalue[5] = 1;
            
            pnum = 6;
            fun = new TF1("Quadratic", bkgfunc_03, start, end, pnum);
            fun->SetParameter(3, pvalue[3]);
            fun->SetParameter(4, pvalue[4]);
            fun->SetParameter(5, pvalue[5]);
        } // Quadratic
        else if (type == 4) {
            if (!pvalue[3]) {
                pvalue[3] = (energy[right] - energy[left]);
                pvalue[3] *= counts[left]*counts[right]*1.0/(counts[left] - counts[right]);
            }
            if (!pvalue[4]) pvalue[4] = pvalue[3]/counts[0] - energy[0];
            
            pnum = 5;
            fun = new TF1("Inverse", bkgfunc_04, start, end, pnum);
            fun->SetParameter(3, pvalue[3]);
            fun->SetParameter(4, pvalue[4]);
        } // Inverse
        
        fun->SetParameter(0, pvalue[0]);
        fun->SetParLimits(0, 0, 100000);
        fun->SetParameter(1, pvalue[1]);
        if (pvalue[1] < 100) fun->SetParLimits(1, 70, 80);
        //fun->SetParLimits(1, 2250, 2750);
        fun->SetParameter(2, pvalue[2]);
        fun->SetParLimits(2, 0, 10);
        //fun->SetParLimits(2, 10, 100);
        
        /*/ Trouble Shooting
         std::cout << "Start: " << start << "\tEnd: " << end << "\n";
         std::cout << "Left: " << left << "\tRight: " << right << "\n";
         std::cout << "counts_left: " << counts[left] << "\n";
         if(mid) std::cout << "counts_mid: " << counts[mid] << "\n";
         std::cout << "counts_right: " << counts[right] << "\n";
         std::cout << "energy_left: " << energy[left] << "\n";
         if (mid) std::cout << "energy_mid: " << energy[mid] << "\n";
         std::cout << "energy_right: " << energy[right] << "\n";
         for (int i = 0; i < pnum; i++) std::cout << "p_" << i << " : " << pvalue[i] << "\t";
         std::cout << "\n\n";
         //*/
        
        hist->Fit(fun, "LEMR");// "LLEMR" for low statistics bin contents < 100
        fitparam = fun->GetParameters();
        fiterror = fun->GetParErrors();
        pointer->fit[pointer->cslice][0] = type;
        pointer->fit[pointer->cslice][1] = (end - start)/2.0;
        pointer->fit[pointer->cslice][2] = 1.0*fun->GetChisquare()/fun->GetNDF();
        
        
        for (int i = 0; i < pnum; i++) {
            pointer->fit[pointer->cslice][2*i + 3] = fitparam[i];
            pointer->fit[pointer->cslice][2*i + 4] = fiterror[i];
        }
        for (int i = pnum; i < 7; i++){
            pointer->fit[pointer->cslice][2*i + 3] = 0;
            pointer->fit[pointer->cslice][2*i + 4] = 0;
        }
        
        // Save fit values...
        snprintf(word, sizeof(word), "_22_Pu_Fit");
        if(pointer->indexvalue < 0){
            snprintf(word, sizeof(word), "_%d-%d", -pointer->indexvalue, -pointer->indexvalue + 21);
        }
        snprintf(str, sizeof(str), "%sPu_Fits/%s_Ch_%d%s.txt", pointer->path, pointer->name, pointer->channel, word);
        printf("Saving to: %s\n", str);
        outfile.open(str);
        snprintf(str, sizeof(str), "000_bkg\t001_Range_%.2f\t002_rchi\t", (end + start)/2.0);
        outfile << str;
        for (int i = 0; i < 7; i++){
            snprintf(str, sizeof(str), "00%d_p%d\t00%d_dp%d%s", 2*i + 3, i, 2*i + 4, i, (i == 6 ? "\n" : "\t"));
            outfile << str;
        }
        for (int i = 0; i <= pointer->cslice; i++) {
            for (int j = 0; j < 17; j++) outfile << pointer->fit[i][j] << (j == 16 ? "\n" : "\t");
        }
        outfile.close();
        fun->Delete();
        if (pointer->fit[pointer->cslice][3] <= pointer->fit[pointer->cslice][4]){
            printf("Area failure, %0.4f +/- %0.4f\n", pointer->fit[pointer->cslice][3], pointer->fit[pointer->cslice][4]);
            exit = 0;
        } // Area Check
        /*/else if (pointer->fit[pointer->cslice][9] >= pow(10, 6)){//<= pointer->fit[pointer->cslice][10]){
            printf("p3 failure, %0.4f +/- %0.4f\n", pointer->fit[pointer->cslice][9], pointer->fit[pointer->cslice][10]);
            exit = 0;
        } //*/// p3 Check
        /*/else if (pointer->fit[pointer->cslice][11] <= pointer->fit[pointer->cslice][12]){
            printf("p4 failure, %0.4f +/- %0.4f\n", pointer->fit[pointer->cslice][11], pointer->fit[pointer->cslice][12]);
            exit = 0;
        }//*/// p4 Check
        else if (78 <= pointer->fit[pointer->cslice][5] || pointer->fit[pointer->cslice][5]<= 74){
            printf("Centroid failure, 74 <= %0.4f <= 78\n", pointer->fit[pointer->cslice][5]);
            exit = 0;
        } // Centroid check
        /*/else if (pointer->fit[pointer->cslice][6] <= 0.05 || 0.11 <= pointer->fit[pointer->cslice][6]){
            printf("dC failure, 0.05 <= %0.4f\n", pointer->fit[pointer->cslice][6]);
            exit = 0;
            
        }/*/// Check dC is greater than 0.05 eV ???
        else if (pointer->fit[pointer->cslice][7] <= 0.5 || 1.5 <= pointer->fit[pointer->cslice][7]){
            printf("Sigma failure, 0.5 <= %0.4f <= 1.5\n", pointer->fit[pointer->cslice][7]);
            exit = 0;
            
        }// Check that sigma is greater than 0.5 eV
        /*/else if (pointer->fit[pointer->cslice][8] <= 0.05){
            printf("dS failure, 0.05 <= %0.4f\n", pointer->fit[pointer->cslice][8]);
            exit = 0;
            
        }/*/// Check that dS is greater than 0.05 eV
        else if (0 < pointer->cslice && pointer->fit[pointer->cslice][5] == pointer->fit[pointer->cslice - 1][5]){
            printf("Cent Identical Fit, %0.4f = %0.4f\n", pointer->fit[pointer->cslice][5], pointer->fit[pointer->cslice - 1][5]);
            exit = 0;
            
        }// Check that prior centroid is not identical to current
        else if (0 < pointer->cslice && pointer->fit[pointer->cslice][7] == pointer->fit[pointer->cslice - 1][7]){
            printf("Sigma Identical Fit, %0.4f = %0.4f\n", pointer->fit[pointer->cslice][7], pointer->fit[pointer->cslice - 1][7]);
            exit = 0;
        }
        else exit = 1;
        //*/
        //exit = 1; // Ignore check...
    }
    printf("Exiting with %d\n", exit);
    return exit;
}

/* **************************************************************************************************** */
void array2hist(info *pointer, TH1D *hist) {
/* **************************************************************************************************** */
    /*
     This function will load the current array in the info file into a specific histogram.
     */
    for(int i = 0; i < hist->GetNbinsX(); i++) hist->Fill(pointer->xaxis[i], pointer->hist[pointer->cslice][i]);
}

/* **************************************************************************************************** */
void array2hist(info *pointer) {
/* **************************************************************************************************** */
    /*
     This function will load the current array in the info file into its corresponding histogram in the 
     info file...
     */
    for(int i = 0; i < pointer->histogram[pointer->cslice]->GetNbinsX(); i++) {
        pointer->histogram[pointer->cslice]->Fill(pointer->xaxis[i], pointer->hist[pointer->cslice][i]);
    }
}

/* **************************************************************************************************** */
template<class type>
void hist2array(type *column, TH1D *hist) {
/* **************************************************************************************************** */
    /*/
     This function will copy over a histogram to an array.
     Histograms start index value at 1.
     /*/
    for(int i = 0; i < hist->GetNbinsX(); i++) column[i] = hist->GetBinContent(i + 1);
}

/* **************************************************************************************************** */
void hist2array(info *pointer) {
/* **************************************************************************************************** */
    /*
     This function will load the current histogram in the info file into its corresponding array in the
     info file...
     */
    for(int i = 0; i < pointer->histogram[pointer->cslice]->GetNbinsX(); i++) {
        pointer->hist[pointer->cslice][i] = pointer->histogram[pointer->cslice]->GetBinContent(i + 1);
        pointer->xaxis[i] = pointer->histogram[pointer->cslice]->GetBinCenter(i + 1);
    }
}

/* **************************************************************************************************** */
int window_save(info* pointer, long long use_time = 15, int numbin = tBin){
/* **************************************************************************************************** */
/*/
    Save window histograms and partition boundaries. Saves histograms stored in pointer->hist[i][...]
    array does nothing to histograms inside pointer->histogram[i].
    
    Will save from the current slice to the end slice...
 
    Will save under pointer->name
 
/*/
    // Initialize Function Variables
    char str[500];
    ofstream outfile;
    // End
    
    snprintf(str, sizeof(str), "%s%s.txt", pointer->path, pointer->name);
    outfile.open(str);
    
    if(outfile.fail()){
        printf("Failed to open file for saving windows: %s\n", str);
    }
    else{
        printf("Saving to\n%s", str);
        outfile << "Time Start (ns):\t" << pointer->thyme[0] << "\n";
        outfile << "Time Finish (ns):\t" << pointer->thyme[pointer->endslice] << "\n";
        outfile << "Window Step (min):\t" << (int)(0.1*use_time/(60*pow(10, 9))) << "\n";
        outfile << "Window Width (min):\t" << (use_time/(60*pow(10, 9))) << "\n";
        outfile << "Partitions:\t" << pointer->endslice << "\n";
        outfile << "x_axis\t";
        for (int i = pointer->cslice; i < pointer->endslice; i++){
            sprintf(str, "Window_%d", i);
            outfile << str << (i == pointer->endslice - 1? "\n" : "\t");
        }
        printf("\nNumber of rows: %d\tStarting Column: %d\tEnding Column: %d\n", numbin, pointer->cslice, pointer->endslice);
        for (int j = 0; j < numbin; j++) for (int i = pointer->cslice; i < pointer->endslice; i++) {
            //printf("row: %d\tcolumn: %d\r", j, i);
            if (i == pointer->cslice) outfile << pointer->xaxis[j] << "\t";
            outfile << pointer->hist[i][j] << (i == pointer->endslice - 1? "\n" : "\t");
        }
        outfile.close();
    }
    printf("Finished...\n");
    return !outfile.fail();

}

/* **************************************************************************************************** */
template<class type>
double avg_diff(type *array, int num, bool trouble = false){
/* **************************************************************************************************** */
    /*/
     This function is intended to work on pointer->current parameters in determining the average
     difference between peak points to determine an estimate for the average peak difference even when
     there are extra or skipped peaks in the data set...
    /*/
    // Initialize Function Variables
    double average = 0;
    double peaks[NUM_PEAKS];
    double difference[NUM_PEAKS];
    
    
    TH1D *hist = new TH1D("Avgdiff", "Avgdiff", tBin, 0, tBin);
    // End
    
    trouble = false;// Here to stop warnings...
    
    for (int i = 1; i < num; i++) hist->Fill(array[i] - array[i - 1]);
    
    
    
    return average;
} // Work in Progress...

/* **************************************************************************************************** */
int monte_carlo(info* pointer, int numcounts = pow(10, 6), int numpeaks = 30, double mu = 10,
                double sigma = 0.8, double lower = 0, double upper = tBin, int numbin = tBin,
                int type = 0, bool trouble = false){
/* **************************************************************************************************** */
    /*/
     This function will generate a monte carlo simulation of the STJ response to the laser for a given 
     set of values in number of peaks, average number of photons, total number of counts and a given set
     of calibration values assuming a photon energy of 3.5eV.
     
     pointer->current should point to the location in pointer with the values of area??? no don't need 
     this because I am using a general mu and number of peaks...
     
     20161130
     The variable statsigsq is used to adjust the sigma of each peak to account for energy resolution
     from statistics.
     statsigsq = E*F*epsilon = 3.5*0.2*1.7delta = 0.7*1.7*0.7*10^-3 = 8.3*10^-4
     
     To introduce a nonlinearity to the monte carlo simulation I have changed the centroids to depend
     on a small quadratic term that will center about 49 eV and decrease by 0.05 eV at 0 and 100 eV.
     nonlinearity = 0.05*(E - 49)^2/45.5^2
     
     20161201
     Time: 900
     I have corrected the nonlinearity to be accounted for by nonconstant binwidths...
     set centnonlin = 0
     
     Time: 1258
     The Area being used was off... I was doing the poisson amplitude not amplitude *sqrt(2pi)*sigma
     changed this to proper values now. This has no effect on the fits just what the expected areas
     should come out to.
     
     20161206
     Time: 1115
     Added Monte Carlo Simulation of actual signal, using sigma as sigma, numpeaks is used as a
     number of minutes that the simulation represents (ie 30 -> 0.5 hrs), mu is used to adjust the
     centroid location. Based on the actual data the incoming rate is 3 counts/0.2eV/hr so this
     value is set at 15 counts/hr * percentage of hour.
     
     20161207
     Time: 2017
     Changed statsigsq: 5*10^-3 
     This takes into account the multiple tunneling F -> F + 1 => 0.2 -> 1.2
     Removed centnonlin as unnecessary since this is achieved through the variable bin size of the
     modeling process...
     
     20170104
     Time: 853
     Added file save of list data for MCS, this will allow me to reload the exact same data and 
     rebin it into variable bins to compare nonlinearity and linearity on identical data not just
     nominally identical data sets...
     
     20170117
     Time: 1244
     In signal simulation changed the tau value to 40 here was 37 before. Additionally changed sigma
     to depend on energy and locked p3 (area of decay) to the area of the gaussian so that the amplitudes
     matched when decay is at 60 eV.
     
     20170119
     Time: 912
     Need to add a piece for simulating the signal, the current time markers work fine for the laser
     however, the signal has significantly fewer events so the periodic times are no good...
     Found an error in prior stuff, rerunning to have proper numbers...
     
     20170525
     Added delta to the code, this enables modeling of the mu dependent energy.
     
    /*/
    // Initialize Function Variable
    TCanvas *canvas;
    char str[500];
    long long thyme = 0;
    double delta = (type < 0 ? TMath::Abs(type) : 0)*mu/22; // Set laser substrate offset as a function of mu for full eV at 77 eV.
    double statsigsq = 5*pow(10, -3);
    double *xaxis = new double[numbin];
    double *yaxis = new double[numbin];
    double *cent = new double[numpeaks];
    TF1 *fun;
    TH1D *quant;
    TH1D *simul;
    TH1D *recon;
    TH1D *hist = new TH1D("FuncForm", "FuncForm", numbin, lower, upper);
    TH1D *monte = pointer->histogram[pointer->cslice];
    TRandom3 *random = new TRandom3(0);
    ofstream outfile;
    // End
    
    // Create function for Monte Carlo
    if (trouble) printf("Setting up Function\n");
    
    if (0 < type){
        sigma = pow(statsigsq*mu/3.5 + sigma*sigma, 0.5);
        fun = new TF1("fun", bkgfunc_00, lower, upper, 5);
        fun->SetParameter(0, numpeaks); // Set gaussian area...
        fun->SetParameter(1, mu); // Set centroid
        fun->SetParameter(2, sigma); // Set sigma
        fun->SetParameter(3, 2*40*pow(2*TMath::Pi()*sigma*sigma, -0.5)*TMath::Exp(60./40)*numpeaks); // Set decay area                    HERE HERE HERE HERE HERE HERE HERE HERE Comment back in!!!!!!!!!!!
        //fun->SetParameter(3, 0);
        fun->SetParameter(4, 40); // Set tau
        numcounts = (1 + 2*40*TMath::Exp(60./40)*pow(2*TMath::Pi()*sigma*sigma, -0.5))*numpeaks/numcounts;
        snprintf(str, sizeof(str), "%sSignal_MCS_List.txt", pointer->path);
    } // Signal Simulations
    else {
        fun = new TF1("fun",func_00, lower, upper, 3*numpeaks + 1);
        fun->FixParameter(0, numpeaks);
        for (int i = 0; i < numpeaks; i++){
            fun->SetParameter(3*i + 1, pow(mu, i + 1)*TMath::Exp(-mu)/TMath::Factorial(i + 1)*pow(2*TMath::Pi()*(statsigsq*((i + 1) + delta/3.5) + sigma*sigma), 0.5)); // Set Area
            fun->SetParameter(3*i + 2, 3.5*(i + 1) + delta); // Set centroid including delta offset...
            fun->SetParameter(3*i + 3, pow(statsigsq*((i + 1) + delta/3.5) + sigma*sigma, 0.5)); // Set Sigma
            snprintf(str, sizeof(str), "%sLaser_MCS_List.txt", pointer->path);
        }
    } // Laser Simulations
    
    // Setup 'list mode' file saving...
    if (pointer->cslice == 0){
        outfile.open(str);
        outfile << "p0\t" << fun->GetParameter(0) << "\n";
        outfile << "p1\t" << fun->GetParameter(1) << "\n";
        outfile << "p2\t" << fun->GetParameter(2) << "\n";
        outfile << "p3\t" << fun->GetParameter(3) << "\n";
        outfile << "p4\t" << fun->GetParameter(4) << "\n";
        outfile << "Numcounts\t" << numcounts << "\n";
        outfile << "000_Time\t001_Energy\n";
    }
    else outfile.open(str, ios::out| ios::app);
    
    thyme = ((long long)(pointer->cslice*numcounts)*10000000); // Assume a 10ms (100 Hz) laser
    // Create histogram for Monte Carlo
    if(trouble) printf("Setting up model histogram\n");
    for (int i = 0; i < numbin; i++) {
        hist->Fill(hist->GetBinCenter(i + 1), fun->Eval(hist->GetBinCenter(i + 1)));
    }
    // Create Monte Carlo
    for(int i = 0; i < numcounts; i++){
        xaxis[0] = random->Rndm();
        hist->GetQuantiles(1, yaxis, xaxis);
        monte->Fill(yaxis[0]);
        outfile << thyme  << "\t";// Assume a 10ms (100 Hz) laser
        outfile << yaxis[0] << "\n";
        thyme+= 10000000;
    }
    
    if (trouble){
        canvas = new TCanvas;
        simul = new TH1D("Sim", "Sim", numbin, lower, upper);
        recon = new TH1D("Recon", "Recon", numbin, lower, upper);
        quant = new TH1D("Quant", "Quant", numbin, lower, upper);
        
        for (int i = 0; i < numbin; i++) xaxis[i] = (i + 1.0)/numbin;
        gRandom->SetSeed();
        //simul->FillRandom("fun", numcounts);
        simul = hist;
        simul->GetQuantiles(numbin, yaxis, xaxis);
        for (int i = 0; i < numbin; i++) quant->Fill(xaxis[i], yaxis[i]);
        for(int i = 0; i < numcounts; i++){
            xaxis[0] = random->Rndm();
            simul->GetQuantiles(1, yaxis, xaxis);
            recon->Fill(yaxis[0]);
        }
        canvas->Divide(1, 2);
        canvas->cd(1);
        hist->Draw("hist");
        simul->Draw("hist same");
        canvas->cd(2);
        recon->Draw("hist");
        monte->Draw("hist same");
        canvas->Update();
        canvas->Update();
    }
    else {
        hist->Delete();
        fun->Delete();
        
    }
    return 1;
}
    
/* **************************************************************************************************** */
int peak_check(info* pointer, int num, float sigma = 5, double lcl = 0, double ucl = 8000,
               int lal = 600, int pspace = 0, bool trouble = false){
/* **************************************************************************************************** */
    /*
     This function is intended to check the current fitted/located peaks to the prior peaks for
     reference and come up with the next guess for peak fitter. 
     
     pointer    ->  data, cslice should point to the current histogram being fitted, the current peaks
                    should be located in param[cslice] and ideally sorted in increasing order.
     num        ->  number of peaks in current slice from TSpectrum or the attempted fit...
     
     */
    
    
    // Initalize Function Variables
    TH1D *hist = pointer->histogram[pointer->cslice];
    
    int exit = 0;
    int counter = 0;
    int pnum[NUM_PEAKS];
    bool track[NUM_PEAKS];
    char str[500];
    float var = 0.25;   // Sets the range a peak may be off from expected to +/- 0.25 eV
    float varoff = 0;   // Used to adjust the center of the allowed range for which the peak lands...
    double avgdiff = 0;
    double cpar[NUM_PEAKS];
    double ppar[NUM_PEAKS];
    double dummyvar = 0;
    double slope = 0;
    double prior = 0;
    double current = 0;
    double offset = 0;
    // End
    
    pspace = 0; // Just here to stop warnings...
    
    // Copy over current and prior peaks, sort current peaks...
    if(trouble) printf("Peak Num\tCurrent\tPrior\n");
    for (int i = 0; i < NUM_PEAKS; i++) {
        pnum[i] = 0;
        cpar[i] = pointer->param[pointer->cslice*NUM_PEAKS + i][2];
        ppar[i] = pointer->param[(pointer->cslice - 1)*NUM_PEAKS + i][2];
        track[i] = false;
        if (trouble && (cpar[i] != 0 || ppar[i] != 0)) printf("%d\t%.2f\t%.2f\n", i, cpar[i], ppar[i]);
    }
    sort_list(cpar, num);
    
    printf("Removing out of range points\n");
    counter = 0;
    for (int i = num - 1; i >= 0; i--){
        if (cpar[i] < lcl || ucl < cpar[i]){
            if(counter == 0) printf("lcl: %.2f\tucl: %.2f\nNum\tCent\n", lcl, ucl);
            printf("%d\t%.2f\n", i, cpar[i]);
            removept(cpar, i, 0, NUM_PEAKS);
            counter++;
        }
    }// Remove peaks out of centroid range
    if(counter == 0) printf("No peaks removed...\n");
    num -= counter;
    
    /*/
     Parameter fist: param => {{fit_number, energy, peak, dpeak, sigma, dsigma, area, darea}_0, ...}
     Calibration fits : calib => {{fit_number, slope, dslope, offset, doffset}_0, ... }
     
     Allow current peaks to differ from prior peak calibration by +/- 0.5 eV
     /*/
    // Take prior calibration fit values as approximates to check current peaks...
    slope = pointer->calib[pointer->cslice - 1][1];
    offset = pointer->calib[pointer->cslice - 1][3];
    offset = (offset < 0 ? -offset: offset);
    avgdiff = 3.5/slope;
    // Mark points within 0.25 of avgdiff of prior peak
    for (int i = 0; i < NUM_PEAKS; i++){
        if (ppar[i] == 0) break;
        for (int j = 0; j < num; j++) {
            if (ppar[i] - 0.25*avgdiff <= cpar[j] && cpar[j] <= ppar[i] + 0.25*avgdiff) track[j] = true;
        }
    }
    
    
    for (int i = 0; i < num; i++){
        dummyvar = cpar[i]*slope - 3.5*round_decimal((cpar[i]*slope - (i == 0 ? 0 : varoff))/3.5, 0); // Amount the current peak is off of expectation in eV...
        if(dummyvar < -3.5 || 3.5 < dummyvar) {
            printf("Error in peak_check(...)\n"); // Peak is out of range for rounding, should never happend...
            printf("Number\tdummyvar\tpeak\tslope\tvaroff\n");
            printf("%d\t%.4f\t%.4f\t%.4f\t%.4f\n", i, dummyvar, cpar[i], slope, varoff);
        }
        if (-(var + offset) <= dummyvar - varoff && dummyvar - varoff <= var + offset) {
        //if (-(var + offset) <= dummyvar && dummyvar <= var + offset) {
            track[i] = true;
            pnum[i] = round_decimal((cpar[i]*slope - varoff)/3.5, 0); // records peak number...nothing now maybe later...
            varoff = dummyvar;
        }
        if(trouble){
            printf("peak_check, dummyvar: %.4f, cpar[i]: %.4f, ", dummyvar, cpar[i]);
            printf("slope: %.4f, abs(offset): %.4f, ", slope, offset);
            printf("varoff: %.4f, expected: %.4f\t", varoff, cpar[i]*slope);
            printf("track: %s\n", (track[i] ? "true": "false"));
        }
    } // Mark good peaks...
    
    for (int i = 0; i < NUM_PEAKS; i++){
        if (pointer->param[(pointer->cslice - 1)*NUM_PEAKS + i][0] == -1) break;
        dummyvar = pointer->param[(pointer->cslice - 1)*NUM_PEAKS + i][2];
        for (int j = 0; j < num; j++) if (!track[i] && dummyvar - 0.2*avgdiff <= cpar[i] && cpar[i] <= dummyvar + 0.2*avgdiff) track[i] = true;
    }
    
    
    
    if(trouble) {
        printf("Average peak distance: %.2f\nPeaks before cleaning\nNum\tPeak\tTrack\n", avgdiff);
        for(int i = 0; i < num; i++) printf("%d\t%.2f\t%s\n", i, cpar[i], (track[i] ? "true": "false"));
    }
    
    printf("Removing none peak points\nNum\tCentroid\n");
    counter = 0;
    for (int i = num - 1; i >= 0; i--){
        if (track[i]) continue;
        printf("%d\t%.2f\n", i, cpar[i]);
        removept(cpar, i, 0, NUM_PEAKS);
        counter++;
    }// Remove none peaks...
    if(counter == 0) printf("No peaks removed...\n");
    num -= counter;
    if(num <= 0){
        printf("Error all peaks have been remove...Reset to all peaks\n");
        for (int i = 0; i < counter; i++) {
            cpar[i] = pointer->param[pointer->cslice*NUM_PEAKS + i][2];
            track[i] = true;
            num = counter;
        }
    } // Removal error, repopulate centroids...
    sort_list(cpar, num);
    
    printf("Removing out of range points\n");
    counter = 0;
    for (int i = num - 1; i >= 0; i--){
        if (cpar[i] < lcl || ucl < cpar[i]){
            if(counter == 0) printf("lcl: %.2f\tucl: %.2f\nNum\tCent\n", lcl, ucl);
            printf("%d\t%.2f\n", i, cpar[i]);
            removept(cpar, i, 0, NUM_PEAKS);
            counter++;
        }
    }// Remove peaks out of centroid range
    if(counter == 0) printf("No peaks removed...\n");
    num -= counter;
    
    counter = 0;
    printf("Combining double points, 0.25*avgdiff: %f\n", 0.25*avgdiff);
    for (int i = num - 1; i >= 0; i--){
        if (0 <= cpar[i] - cpar[i - 1] && cpar[i] - cpar[i - 1] <= 0.25*avgdiff){
            if(counter == 0) printf("cpar[i - 1]\tcpar[i]\n");
            printf("%.2f\t%.2f\n", cpar[i - 1], cpar[i]);
            cpar[i - 1] = (cpar[i - 1] + cpar[i])/2.0;
            cpar[i] = 0;
            counter++;
        }
        if ( cpar[i] == 0) removept(cpar, i, 0, NUM_PEAKS);
    }// Combine repeated peaks
    if(counter == 0) printf("No peaks removed...\n");
    num -= counter;
    
    printf("Missing points being added.\n");
    counter = 0;
    for (int i = 1; i < num; i++){
        prior = (cpar[i - 1]*slope)/3.5;
        current = (cpar[i]*slope)/3.5;
        if (1 < round_decimal(current - prior, 0)){
            if(trouble){
                printf("\t\tPrior\tCurrent\n");
                printf("Peak Number: %.2f\t%.2f\n", prior, current);
                printf("Centroid: %.2f\t%.2f\n", cpar[i - 1], cpar[i]);
                printf("Slope: %.2f\n", slope);
                printf("Number already in list of peaks? %s\n", (element(pnum, round_decimal(prior, 0) + 1, NUM_PEAKS) ? "Yes" : "No"));
            }
            cpar[num] = (round_decimal(prior, 0) + 1)*3.5/slope - offset/slope;
            if (!element(pnum, round_decimal(prior, 0) + 1, NUM_PEAKS) && cpar[i - 1] < cpar[num] && cpar[num] < cpar[i]){
                num++;
                sort_list(cpar, num);
                if(trouble){
                    printf("Sorted peaks:\n");
                    for(int i = 0; i < num; i++) printf("%.2f\n", cpar[i]);
                }
            }
            else {
                printf("ERROR new value is not between prior and current points: %.2f\n", cpar[num]);
                printf("\t\tPrior\tCurrent\n");
                printf("Peak Number: %.2f\t%.2f\n", prior, current);
                printf("Centroid: %.2f\t%.2f\n", cpar[i - 1], cpar[i]);
                printf("Slope: %.2f\n", slope);
                printf("Number already in list of peaks? %s\n\n", (element(pnum, round_decimal(prior, 0) + 1, NUM_PEAKS) ? "Yes" : "No"));
            }
        }
        else if(current <= prior) {
            printf("\n\n\nERROR, current centroid is less than or equal to prior centroid...\n");
            printf("Prior: %.2f\tCurrent: %.2f\n\n\n", prior, current);
        }
    }// Add skipped peaks
    if(counter == 0) printf("No peaks removed...\n");
    
    counter = 0;
    printf("Combining double points, again...\n");
    for (int i = num - 1; i >= 0; i--){
        if (0 <= cpar[i] - cpar[i - 1] && cpar[i] - cpar[i - 1] <= 0.25*avgdiff){
            if(counter == 0) printf("cpar[i - 1]\tcpar[i]\n");
            printf("%.2f\t%.2f\n", cpar[i - 1], cpar[i]);
            cpar[i - 1] = (cpar[i - 1] + cpar[i])/2.0;
            cpar[i] = 0;
            counter++;
        }
        if ( cpar[i] == 0) removept(cpar, i, 0, NUM_PEAKS);
    }// Combine repeated peaks
    if(counter == 0) printf("No peaks removed...\n");
    num -= counter;
    
    if(trouble) {
        printf("Peaks after cleaning\nNum\tPeak\n");
        for(int i = 0; i < num; i++) printf("%d\t%.2f\n", i, cpar[i]);
    }
    // Wrap up centroid checks...
    
    dummyvar = 0;
    for (int i = 0; i < num; i++) pointer->param[pointer->cslice*NUM_PEAKS + i][2] = cpar[i];
    for (int i = 0; i < num; i++) {
        dummyvar = hist->GetXaxis()->FindBin(pointer->param[pointer->cslice*NUM_PEAKS + i][2]);
        dummyvar = (hist->GetBinContent(dummyvar - 2) + hist->GetBinContent(dummyvar - 1) +
                    hist->GetBinContent(dummyvar) + hist->GetBinContent(dummyvar + 1) +
                    hist->GetBinContent(dummyvar + 2))*sigma/5.0;
        pointer->param[pointer->cslice*NUM_PEAKS + i][4] = (sigma < 0 ? -sigma : sigma);
        pointer->param[pointer->cslice*NUM_PEAKS + i][6] = pow(2*TMath::Pi(), 0.5)*(dummyvar < 0 ? -dummyvar : dummyvar);
    } // Set initial parameter guesses
    
    if(trouble) {
        printf("Preclean parameters\nNum\tPeak\tSigma\tArea\n");
        for(int i = 0; i < num; i++) printf("%d\t%.2f\t%.2f\t%.2f\n", i, pointer->param[pointer->cslice*NUM_PEAKS + i][2],
                                            pointer->param[pointer->cslice*NUM_PEAKS + i][4],
                                            pointer->param[pointer->cslice*NUM_PEAKS + i][6]);
    }
    
    counter = num - 1;
    while (pointer->param[NUM_PEAKS*pointer->cslice + counter][6] < lal){
        printf("Removing point: %d\tArea: %.2f\n", counter, pointer->param[NUM_PEAKS*pointer->cslice + counter][6]);
        removept(pointer, counter, NUM_PEAKS*pointer->cslice, NUM_PEAKS*(pointer->cslice + 1));
        counter--;
        if (counter < 0) break;
    } // Removes peaks from the right until first peak above threshold
    num = counter + 1;
    counter = 0;
    while (pointer->param[NUM_PEAKS*pointer->cslice][6] < lal){
        printf("peak_check, Removing point: %d\tArea: %.2f\n", counter, pointer->param[NUM_PEAKS*pointer->cslice][6]);
        removept(pointer, 0, NUM_PEAKS*pointer->cslice, NUM_PEAKS*(pointer->cslice + 1));
        printf("Successfully removed...\n");
        counter++;
        if (counter >= num) break;
    } // Removes peaks from the left until first peak above threshold
    num -= counter;
    counter = NUM_PEAKS*pointer->cslice;
    for (int i = 1; i < num - 1; i++){
        if(pointer->param[counter + i][6] - pointer->param[counter + i - 1][6] < 0 &&
           pointer->param[counter + i][6] - pointer->param[counter + i + 1][6] < 0){
            pointer->param[counter + i][6] = pointer->param[counter + i + 1][6];
            pointer->param[counter + i][6] += pointer->param[counter + i - 1][6];
            pointer->param[counter + i][6] = pointer->param[counter + i][6]/2.0;
        }
    }// No local minima...
    
    for (int i = num; i < NUM_PEAKS; i++) {
        pointer->param[pointer->cslice*NUM_PEAKS + i][2] = 0;
        pointer->param[pointer->cslice*NUM_PEAKS + i][4] = 0;
        pointer->param[pointer->cslice*NUM_PEAKS + i][6] = 0;
    }
    if(trouble) {
        printf("Input parameters\nNum\tPeak\tSigma\tArea\n");
        for(int i = 0; i < num; i++) printf("%d\t%.2f\t%.2f\t%.2f\n", i, pointer->param[pointer->cslice*NUM_PEAKS + i][2],
                                            pointer->param[pointer->cslice*NUM_PEAKS + i][4],
                                            pointer->param[pointer->cslice*NUM_PEAKS + i][6]);
    }
    exit = num;
    return exit;
}

/* **************************************************************************************************** */
int parse_time(info* pointer, long long use_time = 15, float val = 5, bool setval = true,
               int numbin = tBin, int lower = 0, int upper = tBin){
/* **************************************************************************************************** */
    /* 
        This function is intended to parse the data inside of pointer into arrays that can be readily
        saved and moved into and out of histograms for fitting purposes.
        use_time    => This is the time window in minutes (default is 15 minutes).
        val         => This is the fixed window size in minutes (default is 5 minutes)
        setval      => This switches the code to do a fixed window size of val over the use_time frame
        numbin      => Number of bins in histogram to be generated, default is tBin
        lower       => Lower limit of histogram, default is 0 (inclusive)
        upper       => Upper limit of histogram, default is tBin (exclusive)
     
     New addition, this function will now save the generated histograms into pointer->hist in order to
     save the generated histograms for future processing without the need to reload the entire data
     file again.
     */
    // Initialize Function Variables
    TH1D *host[50];
    TH1D *hist;
    char titles[50][500];
    char str[500];
    int partition = 0;
    int numpart = 10;
    ofstream outfile;
    // End
    
    if(use_time > 0){
        printf("Parsing data for channel: %d\n", pointer->channel);
        pointer->current = 0;
        pointer->cslice = 0;
        pointer->endslice = 0;
    }
    else use_time = 15;
    // Advance current position on pointer to laser event and store the time.
    while(pointer->current < pointer->last && !pointer->tag[pointer->current]) pointer->current++;
    if (pointer->last <= pointer->current) {
        printf("\n\n\nERROR ERROR ERROR\n\n\nNo laser events found!!!!\n\n\nERROR ERROR ERROR\n\n\n");
    }
    pointer->thyme[partition] = pointer->time[pointer->current];
    
    /*/
     Use fixed window size set by val
     /*/
    if (setval) {
        use_time = (long long)(val*60*pow(10, 9)); // Convert to nanoseconds
        pointer->numbin = pointer->histogram[pointer->endslice]->GetNbinsX();
        while(pointer->current < pointer->last){
            hist = pointer->histogram[pointer->endslice];
            hist->Reset();
            while(pointer->time[pointer->current] < pointer->thyme[partition] + use_time){
                if (pointer->tag[pointer->current]) hist->Fill(pointer->bin[pointer->current]);
                pointer->current++;
                if(pointer->last <= pointer->current) break;
            }
            if (pointer->endslice == 0) {
                for (int i = 0; i < tBin; i++) pointer->xaxis[i] = 0;
                for (int i = 0; i < pointer->numbin; i++) pointer->xaxis[i] = hist->GetBinCenter(i + 1);
            }// Clear xaxis if this is first slice and copy over
            hist2array(pointer->hist[pointer->endslice], hist);
            if (pointer->current == pointer->last) {
                pointer->current--;
                if (0.9*use_time <= pointer->time[pointer->current] - pointer->thyme[partition]){
                    pointer->endslice++;
                    partition++;
                    pointer->thyme[partition] = pointer->time[pointer->current];
                }
                break;
            }
            else if (pointer->current < pointer->last) {
                pointer->endslice++;
                partition++;
                pointer->thyme[partition] = pointer->time[pointer->current];
            }
        }
        partition++; // Make sure the last partition value is read out...
    } // For set partition range...
    else{
        /*
         Generates the first numpart partitions from the data adding up to a total of use_time and logs the
         enclosing time partitions into the time partition array called thyme.
         */
        use_time = use_time*60*pow(10, 9); // Convert to nanoseconds
        sprintf(str, "%s%s_Part_Hist_Ch_%d.txt", pointer->path, pointer->name, pointer->channel);
        outfile.open(str);
        for (int i = 0; i < numpart; i++) {
            partition++;
            sprintf(titles[i], "Partition_%d", i);
            host[i] = new TH1D(titles[i], titles[i], numbin, lower, upper);
            
            while(pointer->current < pointer->last &&
                  pointer->time[pointer->current] <= pointer->thyme[partition - 1] + 0.1*use_time){
                if (pointer->tag[pointer->current]) host[i]->Fill(pointer->bin[pointer->current]);
                pointer->current++;
            }
            pointer->thyme[partition] = pointer->time[pointer->current];
            outfile << "Partition_" << partition - 1 << "-" << partition << "\t";
            for (int j = 0; j < numbin; j++) outfile << host[i]->GetBinContent(j + 1) << (j == numbin - 1 ? "\n" : "\t");
        }
        
        /*                                      Cyclic
         The main histogram is cleared and reconstructed from the existing 10 partitions and written to
         the histogram array of the pointer. The earliest partition is then overwritten with the next
         partition of size use_counts/10. The partition count is then advanced along with the end
         position of the histogram array. Finally the enclosing time is added to the time partition
         array called thyme.
         
         The individual partition histograms are saved in rows not columns.
         */
        if (!pointer->cslice) for (int i = 0; i < tBin; i++) pointer->xaxis[i] = 0;// Clear xaxis if this is first slice
        for (int i = 0; i < numbin; i++) pointer->xaxis[i] = host[0]->GetBinCenter(i + 1);// Copy over xaxis
        pointer->endslice = 0;
        while(pointer->current < pointer->last){
            hist = pointer->histogram[pointer->endslice];
            hist->Reset();
            for (int i = 0; i < numpart; i++) hist->Add(host[i]);
            hist2array(pointer->hist[pointer->endslice], hist);
            host[partition%numpart]->Reset();
            while(pointer->current < pointer->last &&
                  pointer->time[pointer->current] <= pointer->thyme[partition] + 0.1*use_time){
                if (pointer->tag[pointer->current]) host[partition%numpart]->Fill(pointer->bin[pointer->current]);
                pointer->current++;
            }
            if (pointer->current < pointer->last) {
                outfile << "Partition_" << partition << "-" << partition + 1 << "\t";
                for (int j = 0; j < numbin; j++) outfile << host[partition%numpart]->GetBinContent(j + 1) << (j == numbin - 1 ? "\n" : "\t");
                pointer->endslice++;
                partition++;
                pointer->thyme[partition] = pointer->time[pointer->current];
            }
            pointer->numbin = hist->GetNbinsX();
        }
        for (int i = 0; i < numpart; i++) host[i]->Delete();
        outfile.close();
    } // For moving window average...
    pointer->cslice = 0; // Set so that all histograms are saved...
    window_save(pointer, use_time, pointer->numbin);
    //sprintf(str, "%s%s_Part_Bounds_Ch_%d.txt", pointer->path, pointer->name, pointer->channel);
    sprintf(str, "%s%s_Part_Bounds.txt", pointer->path, pointer->name);
    outfile.open(str);
    outfile << "Boundary\tTime(ns)\n";
    for (int i = 0; i < partition; i++) outfile << i << "\t" << pointer->thyme[i] << "\n";
    outfile.close();
    return 1;
}

/* **************************************************************************************************** */
int ngaussfit(info* pointer,  int num = 0, double ranleft = 0, double ranright = 0, float minsep = 0,
              float maxsep = 0, bool override = false, bool usexaxis = false){
/* **************************************************************************************************** */
    /*
     This function will perform an N-Gaussian fit to the data over the selected range and given
     variable parameters.
     ranleft    =>  This is the lower bound of the fitting range.
     ranright   =>  This is the upper bound of the fitting range.
     minsep     =>  This is the minimum separation for two consecutive peaks. (In Channels)
     maxsep     =>  This is the maximum separation for two consecutive peaks. (In Channels)
    */
    
    // Initalize Function Variables
    TSpectrum *peak_off = new TSpectrum(NUM_PEAKS);
    TH1D *hist;
    TF1 *fun;
    int exit = 1;
    int order[NUM_PEAKS];
    int counter = 0;
    bool trouble = false; // Used to activate trouble shooting features...
    float perror = 0.5;
    double avgdiff = 0;
    double comp = 0;
    double centroid[NUM_PEAKS];
    double sigma[NUM_PEAKS];
    double area[NUM_PEAKS];
    double error[NUM_PEAKS];
    double check = 0;
    double pcent = 0;
    double psigma = 0;
    double parea = 0;
    double avgsig = 0;
    double *fitparam;
    const double *fiterror;
    // End
    
    //std::cout << "Fitting superposition of N-Gaussians...\n"; // Trouble shooting
    // Prepare unaltered input parameters to match data
    if (!num) while (pointer->param[pointer->cslice*NUM_PEAKS + num][2]) num++;
    if (!minsep || !maxsep || !ranleft || !ranright) {
        for (int i = 0; i < num; i++) centroid[i] = pointer->param[pointer->cslice*NUM_PEAKS + i][2];
        for (int i = num; i < NUM_PEAKS; i++) centroid[i] = 0;
        if (!minsep) minsep = 0.25*(centroid[num-1] - centroid[0])/(num - 1);
        if (!maxsep) maxsep = 1.75*(centroid[num-1] - centroid[0])/(num - 1);
        if (!ranleft) ranleft = centroid[0] - minsep;
        if (!ranright) ranright = centroid[num - 1] + minsep;
    }
    if(0 < pointer->cslice) avgdiff = 3.5/pointer->calib[pointer->cslice - 1][1];
    else {
        avgdiff += pointer->param[pointer->cslice*NUM_PEAKS + num - 1][2]/2.0/(num - 1);
        avgdiff -= pointer->param[pointer->cslice*NUM_PEAKS][2]/2.0/(num - 1);
    }
    if (avgdiff <= 0){
        printf("Average difference between peaks is %.2f\nUsing default value of 70\n", avgdiff);
        avgdiff = 150; //-avgdiff;
    }
    if (trouble) {
        printf("Starting input N-Gaussian parameters\tCurrent slice: %d\n", pointer->cslice);
        printf("Range is (%.2f, %.2f)\n", ranleft, ranright);
        printf("Sigma limits: 0.1 to %d\n", 100);
        printf("usexaxis value: %d\n", usexaxis);
        printf("Initial values:\nPeak\tCentroid\tSigma\tArea\n");
        for (int i = 0; i < num; i++) {
            printf("%d\t%.2f\t%.2f\t%.2f\n", i, pointer->param[pointer->cslice*NUM_PEAKS + i][2],
                   pointer->param[pointer->cslice*NUM_PEAKS + i][4],
                   pointer->param[pointer->cslice*NUM_PEAKS + i][6]);
        }
    }
    
    /*
        Copy over histogram slice data and setup N-Gaussian function with parameters:
                {fit_number, energy, peak, dpeak, sigma, dsigma, area, darea}
        Note that for TF1 the optional entries are:
                "R" : use the set function range
                "Q" : Minimal printing
                "O" : Do not plot
    */
    hist = pointer->histogram[pointer->cslice]; // Point to the histogram being fitted...
    //hist->Sumw2(kFALSE);
    //hist->Sumw2();
    
    if (trouble){
        printf("Peak\tCentroid\tSigma\tArea\n");
        for (int i = 0; i < num; i++){
            std::cout << pointer->cslice*NUM_PEAKS + i << "\t";
            std::cout << pointer->param[pointer->cslice*NUM_PEAKS + i][2] << "\t";
            std::cout << pointer->param[pointer->cslice*NUM_PEAKS + i][4] << "\t";
            std::cout << pointer->param[pointer->cslice*NUM_PEAKS + i][6] << "\n";
        }
    }
    
    // Set initial fit parameters and fit limits
    fun = new TF1("fun",func_00, ranleft, ranright, 3*num + 1);
    fun->FixParameter(0, num);
    for (int i = 0; i < num; i++){
        pcent = pointer->param[pointer->cslice*NUM_PEAKS + i][2]; // Centroid
        psigma = pointer->param[pointer->cslice*NUM_PEAKS + i][4];// Sigma
        fun->SetParameter(3*i + 1, pointer->param[pointer->cslice*NUM_PEAKS + i][6]); // Set Area
        fun->SetParameter(3*i + 2, pcent); // Set centroid
        fun->SetParameter(3*i + 3, psigma); // Set sigma
        if(!override){
            fun->SetParLimits(3*i + 3, 0.1, 100); // Limit sigma values to positive < 100
            fun->SetParLimits(3*i + 2, pcent - (usexaxis ? 0.5 : avgdiff), pcent + (usexaxis ? 0.5 : avgdiff)); // Limit centroid to +/- average distance between peaks
            fun->SetParLimits(3*i + 1, 1, pow(10, 8)); // Limit area values to positive
            //fun->SetParLimits(3*i + 1, 1, pow(10, 4)); // Limit area values to positive
        }
    }
    hist->Fit(fun, "QRLEM");
    fitparam = fun->GetParameters();
    fiterror = fun->GetParErrors();
    
    /*/
     Determine order of fitted peaks. Check for double counting or skipped peaks, if no
     error copy over 'good' peaks ignore bad ones. First step is to sort the fitted peaks in
     order of increasing centroid value. I am looking for peaks that are to closely spaced
     or to far apart indicating a missing peak. For the sigma, I am checking that the fit
     did not converge on the limits
    /*/
    for (int i = 0; i < NUM_PEAKS; i++) order[i] = centroid[i] = sigma[i] = area[i] = 0;
    for (int i = 0; i < num; i++) for (int j = i + 1; j < num; j++) {
        fitparam[3*i + 2] <= fitparam[3*j + 2] ? order[j]++ : order[i]++;
    }
    // Check centroid values...
    for (int i = 0; i < num; i++) centroid[order[i]] = fitparam[3*i + 2];
    for (int i = 0; i < num; i++) {
        pcent = pointer->param[pointer->cslice*NUM_PEAKS + order[i]][2];
        if (0 < i && centroid[i] - centroid[i - 1] <= 0.75*avgdiff) {
            printf("Proximity alert\ni - 1: %.2f\ti: %.2f\t", centroid[i - 1], centroid[i]);
            printf("Difference: %.2f\t<= 0.75*Avgdiff: %.2f\n", centroid[i] - centroid[i - 1], 0.75*avgdiff);
            exit = 0;
            centroid[i - 1] = 0;
        }
        else if(0 < i && 2.5*avgdiff <= centroid[i] - centroid[i - 1]){
            printf("Skipped alert\ni: %.2f\ti + 1: %.2f\t", centroid[i], centroid[i + 1]);
            printf("Difference: %.2f\t2.5*Avgdiff: %.2f\n", centroid[i] - centroid[i - 1], 2*avgdiff);
            exit = 0;
            centroid[i - 1] = 0;
        }
        else if(round_decimal(centroid[i], 1) == round_decimal(pcent - (usexaxis ? 0.5 : avgdiff), 1) ||
                round_decimal(centroid[i], 1) == round_decimal(pcent + (usexaxis ? 0.5 : avgdiff), 1)){
            printf("Fitted centroid has reached limits\nLower: %.2f\tCent:  %.2f\tUpper: %.2f\n",
                   round_decimal(pcent - (usexaxis ? 0.5 : avgdiff), 1), round_decimal(centroid[i], 1),
                   round_decimal(pcent + (usexaxis ? 0.5 : avgdiff), 1));
            exit = 0;
            centroid[i] = 0;
        }
        else if (0 < i) centroid[i - 1] = 1;
        if (i == num - 1 && !centroid[i - 1]) centroid[i] = 0;
        else if ( i == num - 1) centroid[i] = 1;
    }
    // Check sigma values
    for (int i = 0; i < num; i++) sigma[order[i]] = fitparam[3*i + 3];
    for (int i = 0; i < num; i++) {
        psigma = pointer->param[pointer->cslice*NUM_PEAKS + order[i]][4];
        if(round_decimal(sigma[i], 4) == 100 || round_decimal(sigma[i], 4) == 0.1){
            printf("Sigma value is at rails\nsigma: %.2f\t for initial value: %.2f\n", sigma[i], psigma);
            sigma[i] = 0;
            exit = 0;
        }
        else if (sigma[i] <= 0) {
            printf("Sigma value is negative\n");
            sigma[i] = 0;
            exit = 0;
        }
        else sigma[i] = 1;
    }
    // Check area values
    for (int i = 0; i < num; i++) area[order[i]] = fitparam[3*i + 1];
    for (int i = 0; i < num; i++){
        if (pow(10, 9) < fitparam[3*i + 1]){
            printf("ith Area (%d) values reaching upper limit: %.2f\n", order[i], fitparam[3*i + 1]);
            area[order[i]] = 0;
            exit = 0;
        }
        else if (area[i] <= 0) {
            printf("Area value is negative\n");
            area[i] = 0;
            exit = 0;
        }
        else area[order[i]] = 1;
    }
    // Check error values
    for (int i = 0; i < NUM_PEAKS; i++) error[i] = 0;
    for (int i = 0; i < num ; i++){
        /*/ Centroid error check dC = sig/sqrt(BinWidth*Area)
        if (fitparam[3*i + 3]*5.0/pow(hist->GetBinWidth(1)*fitparam[3*i + 1], 0.5) < fiterror[3*i + 2]) {
            printf("i: %d\t Centroid Uncertainty issue\ndC: %.2f\tsig: %.2f\tsqrt(area): %.2f\n", order[i], fiterror[3*i + 2], fitparam[3*i + 3], pow(fitparam[3*i + 1], 0.5));
            printf("5*dC_Theory: %.2f\n", fitparam[3*i + 3]*5.0/pow(fitparam[3*i + 1], 0.5));
            centroid[order[i]] = 0;
            exit = 0;
        }/*/
        // Area error check dA/A < 50%
        if (perror <= fiterror[3*i + 1]/fitparam[3*i + 1]) {
            printf("i: %d\t Area uncertainty to high\nArea: %.2f\tdA: %.2f\n", order[i], fitparam[3*i + 1], fiterror[3*i + 1]);
            area[order[i]] = 0;
            exit = 0;
        }
        // Sigma error check dS/S < 50%
        if (perror <= fiterror[3*i + 3]/fitparam[3*i + 3]) {
            printf("i: %d\t Sigma uncertainty to high\nSigma: %.2f\tdS: %.2f\n", order[i], fitparam[3*i + 3], fiterror[3*i + 3]);
            sigma[order[i]] = 0;
            exit = 0;
        }
        /*/ Centroid error check |dS/dC - 1| < 0.5 (Based on plots)
        if (fiterror[3*i + 3]/fiterror[3*i + 2] < 0.5 || 1.5 < fiterror[3*i + 3]/fiterror[3*i + 2]){
            printf("i: %d\t dS/dC: %.2f\ndS: %.2f\tdC: %.2f\n", order[i], fiterror[3*i + 3]/fiterror[3*i + 2], fiterror[3*i + 3], fiterror[3*i + 2]);
            centroid[order[i]] = 0;
            sigma[order[i]] = 0;
            exit = 0;
        }//*/
    }
    
    if (trouble){
        printf("Slot location: %d\n", NUM_PEAKS*pointer->cslice);
        printf("Peak\tCentroid\tSigma\tArea\n");
        for (int i = 0; i < num; i++) {
            printf("%d\t%.2f\t%.2f\t%.2f\n", order[i], fitparam[3*i + 2], fitparam[3*i + 3], fitparam[3*i + 1]);
        }
        printf("Check Values\nOrder\tCent\tSig\tArea\n");
        for (int i = 0; i < num; i++) printf("%d\t%.2f\t%.2f\t%.2f\n", order[i], centroid[i], sigma[i], area[i]);
    }
    if (exit == 0) {
        printf("Failed to fit N-Gaussian...\tCurrent slice: %d out of %d\n", pointer->cslice, pointer->endslice);
        printf("Range is (%.2f, %.2f)\n", ranleft, ranright);
        printf("Sigma limits: 0.1 to %d\n", 100);
        printf("usexaxis value: %d\n", usexaxis);
        printf("Initial values:\nPeak\tCentroid\tSigma\tArea\n");
        for (int i = 0; i < num; i++) {
            printf("%d\t%.2f\t%.2f\t%.2f\n", i, pointer->param[pointer->cslice*NUM_PEAKS + i][2],
                   pointer->param[pointer->cslice*NUM_PEAKS + i][4],
                   pointer->param[pointer->cslice*NUM_PEAKS + i][6]);
        }
        printf("Fit values: \n");
        for (int i = 0; i < num; i++) {
            std::cout << order[i] << "\t" << fitparam[3*i + 2] << "\t";
            std::cout << fitparam[3*i + 3] << "\t" << fitparam[3*i + 1] << "\n";
        }
    }
    for (int i = 0; i < num; i++) {
        if ((centroid[order[i]] && sigma[order[i]] && area[order[i]]) || override) {
            pointer->param[NUM_PEAKS*pointer->cslice + order[i]][0] = pointer->cslice;   // Slice number for fits
            pointer->param[NUM_PEAKS*pointer->cslice + order[i]][1] = order[i];          // Peak Order, this will eventually become energy
            pointer->param[NUM_PEAKS*pointer->cslice + order[i]][2] = fitparam[3*i + 2]; // Copy Centroid value (Ch)
            pointer->param[NUM_PEAKS*pointer->cslice + order[i]][3] = fiterror[3*i + 2]; // Copy Centroid uncertainty
            pointer->param[NUM_PEAKS*pointer->cslice + order[i]][4] = fitparam[3*i + 3]; // Copy Sigma value (Ch)
            pointer->param[NUM_PEAKS*pointer->cslice + order[i]][5] = fiterror[3*i + 3]; // Copy Sigma uncertainty
            pointer->param[NUM_PEAKS*pointer->cslice + order[i]][6] = fitparam[3*i + 1]; // Copy Area value (Ch*Counts)
            pointer->param[NUM_PEAKS*pointer->cslice + order[i]][7] = fiterror[3*i + 1]; // Copy Area uncertainty
        }
        else {
            pointer->param[NUM_PEAKS*pointer->cslice + order[i]][0] = pointer->cslice;   // Slice number for fits
            pointer->param[NUM_PEAKS*pointer->cslice + order[i]][1] = order[i];          // Peak Order, this will eventually become energy
            pointer->param[NUM_PEAKS*pointer->cslice + order[i]][3] = 0;                 // Copy Centroid uncertainty
            pointer->param[NUM_PEAKS*pointer->cslice + order[i]][5] = 0;                 // Copy Sigma uncertainty
            pointer->param[NUM_PEAKS*pointer->cslice + order[i]][7] = 0;                 // Copy Area uncertainty
        }
    }
    for( int i = 0; i < num; i++){
        if (round_decimal(pointer->param[NUM_PEAKS*pointer->cslice + i][4], 5) == 100) {
            printf("Potential Sigma Limit\n");
            exit = 0;
        }
        if (0.5*pow(10, 8) <= pointer->param[NUM_PEAKS*pointer->cslice + i][6]) {
            printf("Potential Area Limit\n");
            exit = 0;
        }
    }
    
    if(trouble){
        std::cout << "First partition fitted: \n";
        std::cout << "Peak\tCentroid\tSigma\tArea\n";
        for (int i = 0; i < num; i++) {
            std::cout << order[i] << "\t" << fitparam[3*i + 2] << "\t";
            std::cout << fitparam[3*i + 3] << "\t" << fitparam[3*i + 1] << "\n";
        }
        std::cout << "Copied over: \n";
        std::cout << "Peak\tCentroid\tSigma\tArea\n";
        for (int i = 0; i < num; i++) {
            std::cout << order[i] << "\t";
            std::cout << pointer->param[NUM_PEAKS*pointer->cslice + order[i]][2] << "\t";
            std::cout << pointer->param[NUM_PEAKS*pointer->cslice + order[i]][4] << "\t";
            std::cout << pointer->param[NUM_PEAKS*pointer->cslice + order[i]][6] << "\n";
        }
    }
    for (int i = num; i < NUM_PEAKS; i++) pointer->param[NUM_PEAKS*pointer->cslice + i][0] = -1;
    if(false){
        printf("Fit was %s\n\n\n", (exit ? "Successful" : "Unsuccessful"));
        exit = 0;
    }
    if (trouble) {
        printf("Exit input N-Gaussian parameters\tCurrent slice: %d\n", pointer->cslice);
        printf("Range is (%.2f, %.2f)\n", ranleft, ranright);
        printf("Sigma limits: 0.1 to %d\n", 100);
        printf("usexaxis value: %d\n", usexaxis);
        printf("Initial values:\nPeak\tCentroid\tSigma\tArea\n");
        for (int i = 0; i < num; i++) {
            printf("%d\t%.2f\t%.2f\t%.2f\n", i, pointer->param[pointer->cslice*NUM_PEAKS + i][2],
                   pointer->param[pointer->cslice*NUM_PEAKS + i][4],
                   pointer->param[pointer->cslice*NUM_PEAKS + i][6]);
        }
    }
    if (override) exit = 1;
    return exit;
}

/* **************************************************************************************************** */
int calibration(info* pointer, int lal = 800, bool specialize = false, bool trouble = false,
                float lstart = 63){
/* **************************************************************************************************** */
    /*
     This function will take the data in pointer and perform a fit to the individual arrays in the
     histogram array. This function explicitly does not use the single (3.5 eV) photon peak to
     calculate the calibration of the energy. Additionally, beyond using the peaks to get the proper
     energy number only centroids with greater than 25% of the maximum centroid area are used in the
     final calibration fit.
     
     domain  =>  bool, indicates if the input values are in eV or bin, default is false for bin
     
     Currently, sigma is not being used and is commented out...
     
     
     */
    // Initialize Function Variables:
    int exit = 0;
    int cur = pointer->cslice;
    int counter = 0;
    int pshift = 0; // This will keep track of the total number of peaks shifted
    int attempt = 0;
    int mainnum = 0;
    bool usepeak[2*NUM_PEAKS];
    bool domain = true;
    float unit = 1.0;
    //float lstart = 70;
    double slope = 0;
    double offset = 0;
    double area[2*NUM_PEAKS];
    double energy[2*NUM_PEAKS];
    double cent[2*NUM_PEAKS];
    double dcent[2*NUM_PEAKS];
    double *out;
    double prior_energy[2*NUM_PEAKS];
    double prior_cent[2*NUM_PEAKS];
    const double *outerr;
    
    TGraphErrors *linefit;
    TF1 *fun;
    // End
    
    if (trouble){
        std::cout << "Starting calibration...\n";
        std::cout << "Block read from pointer->param\n";
        std::cout << "Fit#\tEnergy\tPeak\tdPeak\tSigma\tdSigma\tArea\tdArea\n";
        for (int i = cur*NUM_PEAKS; i < (cur + 1)*NUM_PEAKS; i++) {
            if ( pointer->param[i][0] == -1){
                for (int j = 0; j < 8; j++) std::cout << pointer->param[i][j] << (j == 7? "\n" : "\t");
                break;
            }
            for (int j = 0; j < 8; j++) std::cout << pointer->param[i][j] << (j == 7? "\n" : "\t");
        }
    } // Trouble Shooting
    
    /*/
     This copies over the parameter values from pointer->param
     Gaussian parameter fits : {{fit_number, energy, peak, dpeak, sigma, dsigma, area, darea}_0, ...}
     Note that fit_number is the partition number and that within the partition there are NUM_PEAKS
     slots dedicated to peaks found. If slots are not used the fit_number is set to -1 to indicate end
     of data for the assigned space. As such the if statement will break from copying the parameters
     once it reaches the empty statement of -1 or the section has reached the assigned max number.
     Then remove the two endpoint fits as these will be prone to the most error.
     /*/
    for (int i = 0; i < 2*NUM_PEAKS; i++) usepeak[i] = false;
    
    // Change i -> counter dependence in filling arrays and scan through entire allocated data block...
    counter = 0;
    for (int i = cur*NUM_PEAKS; i < (cur + 1)*NUM_PEAKS; i++){
        if (pointer->param[i][0] == -1) break;
        energy[counter] = pointer->param[i][1];
        cent[counter] = pointer->param[i][2];
        dcent[counter] = pointer->param[i][3];
        area[counter] = pointer->param[i][6];
        usepeak[counter] = true;
        counter++;
    }
    mainnum = counter;
    for (int i = 0; i < mainnum; i++) if (area[i] + pow(area[i], 0.5) < lal) usepeak[i] = false;
    
    if(trouble){
        printf("Number of peaks originally copied over %d\n", mainnum);
        printf("Values in array for calibration:\nEnergy\tCentroid\tdCent\tArea\tUse\n");
        for (int i = 0; i < mainnum; i++){
            printf("%.2f\t%.2f\t%.2f\t%.2f\t%d\n", energy[i], cent[i], dcent[i], area[i], usepeak[i]);
        }
    }
    
    
    // Determine if the energy array units
    for (int i = 0; i < mainnum; i++) if(0 != (energy[i]/3.5 - round_decimal(energy[i]/3.5, 0))){
        printf("Not in Energy units: %.2f\n", energy[i]);
        unit = 3.5;
        domain = false;
        break;
    }
    
    /*/
     Setup initial slope and off set guess in units of ch/# and ch originally, I have now updated it to
     give units of #/ch and # so no additional calculations are needed later.
     /*/
    if(true || !domain){                                                // true, is here as a temporary thing when 'finalized' and not running individuals trouble shooting commands  this should not be needed...
        slope = (energy[counter - 1] - energy[0])/(cent[counter - 1] - cent[0]);
        offset = energy[0] - cent[0]*slope;
    }
    else{
        slope = pointer->calib[cur][1];
        offset = pointer->calib[cur][3];
    }
    if(trouble) printf("Using slope of %f and offset of %f\n", slope, offset);
    
    pshift = 0;
    while(offset + pshift < -1) pshift++;
    while (1 < offset + pshift) pshift--;
    if(!domain){
        if (0.5 < offset + pshift) pshift--;
        else if(offset + pshift < -0.5) pshift++;
    }
    if(trouble){
        printf("Returned a pshift value of %d and using a unit value of %f\n", pshift, unit);
    }
    
    /*/
     The offset is actually one of two possible options as the offset will be between two points
    /*/
    
    
    for (int i = 0; i <  mainnum; i++) {
        if (!domain) energy[i] = unit*(energy[i] + pshift);
        pointer->param[i + cur*NUM_PEAKS][1] = energy[i];
    }
    
    if (specialize){
        printf("\n\n\nWARNING WARNING WARNING\n");
        printf("Using peaks [%0.2f, %0.2f] eV only\n", lstart, lstart + 21);
        printf("WARNING WARNING WARNING\n\n\n");
        if (trouble) {
            printf("Peaks Available:\nEnergy\tCentroid\n");
            for (int i = 0; i < mainnum; i++) if(usepeak[i]) printf("%.2f\t%.2f\n", energy[i], cent[i]);
        }
        for (int i = 0; i < mainnum; i++) if (energy[i] < lstart|| lstart + 21 < energy[i]) usepeak[i] = false;        // Use to limit calibration to peaks between 61 and 85 eV
        printf("Peaks being used:\nEnergy\tCentroid\tArea\tUse\n");
        counter = 0;
        for (int i = 0; i < mainnum; i++) if(usepeak[i]){
            counter++;
            printf("%.2f\t%.2f\t%.2f\t", energy[i], cent[i], area[i]);
            if(energy[i] < lstart + 3.5 || lstart + 17.5 < energy[i]) usepeak[i] = false;
            printf("%s\n", usepeak[i] ? "Yes" : "No");
        }
    }// Limit calibration to 63.5 - 81.5 eV range...
    else {
        usepeak[0] = false;
        usepeak[mainnum - 1] = false;
    }// Exclude end points from calibration
    
    //printf("Initial guesses:\nSlope: %.2f\tOffset: %.2f\tpshift: %d\n", slope, offset, pshift);
    
    /*/
     Set up to perform linear fit to choosen centroid and energy values. Note that the function fun
     parameter entry values are fun->SetParameters(p0 = offset, p1 = slope). Additionally, beyond the first
     calibration values all subsequent ones use the prior values adjusted for the current peak offset.
     /*/
    attempt = 0;
    fun = new TF1("fun2", "pol1");
    linefit = new TGraphErrors(NUM_PEAKS);
    if(trouble) printf("Peaks Loading to fitter\nEnergy\tCent\n");
    for (int i = 0; i < mainnum; i++) if(usepeak[i]){
        if(trouble) printf("%f\t%f\n", energy[i], cent[i]);
        if (standard){
            linefit->SetPoint(linefit->GetN(), cent[i], energy[i]);
            linefit->SetPointError(linefit->GetN() - 1, dcent[i], 0);
        }
        else{
            if(i != 0 && usepeak[i - 1]){
                linefit->SetPoint(linefit->GetN(), cent[i] - cent[i - 1], 3.5);
                linefit->SetPointError(linefit->GetN() - 1, pow(dcent[i]*dcent[i] + dcent[i-1]*dcent[i-1], 0.5), 0);
            }
        }
    }
    if ((specialize && counter == 7) || !specialize){
        if (trouble) printf("Attempt %d, using slope %f and offset %f\n", attempt, slope, offset);
        while (attempt < 2){
            if (0 < attempt && 0 < cur && pointer->calib[cur - 1][1] != 0) {
                slope = pointer->calib[cur - 1][1]/unit;
                offset = pointer->calib[cur - 1][3]/unit;
            }
            fun->SetParameters(unit*(offset + pshift), unit*slope);
            if (!standard) fun->FixParameter(0, 0);
            linefit->Fit(fun, "Q");
            if(round_decimal(pow(fun->GetParameter(1), 2), 5) <= 3.5/2) break;
            attempt++;
        }
        /*/
         Calibration fits : {{fit_number, slope, dslope, offset, doffset}_0, ... }
         Save fit calibration to calibration array and delete the functioin and graph
         /*/
        //printf("Output values:\nSlope: %.2f\tOffset: %.2f\n", fun->GetParameter(0), fun->GetParameter(1));
        out = fun->GetParameters();
        outerr = fun->GetParErrors();
        pointer->calib[cur][0] = pointer->param[cur*NUM_PEAKS][0];
        pointer->calib[cur][1] = out[1];
        pointer->calib[cur][2] = outerr[1];
        pointer->calib[cur][3] = out[0];
        pointer->calib[cur][4] = outerr[0];
        if (trouble) {
            printf("Successfully performed calibration\nslope: %f\toffset: %f\n", out[1], out[0]);
        }
        exit = 1;
    }
    else{
        pointer->calib[cur][0] = pointer->param[cur*NUM_PEAKS][0];
        pointer->calib[cur][1] = 0;
        pointer->calib[cur][2] = 0;
        pointer->calib[cur][3] = 0;
        pointer->calib[cur][4] = 0;
        printf("Failed to perform calibration\ncounter: %d", counter);
        exit = 0;
    }
    
    if (!fitmodel && exit){
        printf("p1: %f\tp0: %f\n", out[1], out[0]);
        printf("Linear Test: Expected %f\t using %f\t get %f\n", pointer->param[cur*NUM_PEAKS][1], pointer->param[cur*NUM_PEAKS][2], pointer->param[cur*NUM_PEAKS][2]*out[1] + out[0]);
        printf("p2 guess: %f\n", (pointer->param[cur*NUM_PEAKS][1] - pointer->param[cur*NUM_PEAKS][2]*out[1] - out[0])/(pointer->param[cur*NUM_PEAKS][2]*pointer->param[cur*NUM_PEAKS][2]));
        delete fun;
        fun = new TF1("fun3", "pol2");
        //if (pointer->cslice == 0) fun->SetParameters(pointer->calib[cur][3], pointer->calib[cur][1], 0.001); // Original Line 191024
        if (pointer->cslice == 0) fun->SetParameters(pointer->calib[cur][3], pointer->calib[cur][1],
                                                     (pointer->param[cur*NUM_PEAKS][1] - pointer->param[cur*NUM_PEAKS][2]*out[1]
                                                      - out[0])/(pointer->param[cur*NUM_PEAKS][2]*pointer->param[cur*NUM_PEAKS][2]));
        else fun->SetParameters(pointer->calib[cur - 1][3], pointer->calib[cur - 1][1], pointer->calib[cur - 1][5]);
        //fun->FixParameter(0, 0); // Fixes offset to 0
        //fun->SetParLimits(0, -2*pointer->calib[cur][3], 2*pointer->calib[cur][3]); // Non-zero offset
        //fun->SetParLimits(1, 0.75*pointer->calib[cur][1], 1.25*pointer->calib[cur][1]);
        //fun->SetParLimits(2, -1, 1);
        linefit->Fit(fun, "Q");
        out = fun->GetParameters();
        outerr = fun->GetParErrors();
        pointer->calib[cur][0] = pointer->param[cur*NUM_PEAKS][0];
        pointer->calib[cur][1] = out[1];                                // Linear Term
        pointer->calib[cur][2] = outerr[1];                             // Linear Uncertainty Term
        pointer->calib[cur][3] = out[0];                                // Offset Term
        pointer->calib[cur][4] = outerr[0];                             // Offset Uncertainty Term
        pointer->calib[cur][5] = out[2];                                // Quadratic Term
        pointer->calib[cur][6] = outerr[2];                             // Quadratic Uncertainty Term
        printf("p2: %f\tp1: %f\tp0: %f\n", out[2], out[1], out[0]);
        printf("Quadratic Test: Expected %f\t using %f\t get %f\n", pointer->param[cur*NUM_PEAKS][1], pointer->param[cur*NUM_PEAKS][2], pointer->param[cur*NUM_PEAKS][2]*pointer->param[cur*NUM_PEAKS][2]*out[2] + pointer->param[cur*NUM_PEAKS][2]*out[1] + out[0]);
        exit = 1;
    }
    fun->Delete();
    linefit->Delete();
    return exit;
}

/* **************************************************************************************************** */
int peakfinder(info* pointer,  bool *override, double *sig, double *amp, bool usexaxis = false,
               int lcl = 100, int ucl = 600, int lal = 600, int pspace = 13, bool useenergy = false,
               float fracleft = 2, float fracright = 2, bool trouble = false){
/* **************************************************************************************************** */
    /* 
        This function will locate the peaks for a single histogram. Note that as of now nothing in this
        function attempts to quantify if a peak is in fact a peak although a small attempt is made at
        making sure peaks are not double counted or skipped.
        
        poff_amp    =>  This is the amp fraction value for TSpectrum to do an initial peak find, it is
                        the fraction of the largest amplitude that sets a bound on the acceptable lowest
                        amplitude a located peak may have. For the highest gain setting a value of 0.035 
                        is needed, at the moderate gain setting a value of 0.005 is used.
        poff_sig    =>  This is the sig value for TSpectrum to do an initial peak find, it is the
                        expected width of the sigma for a located peak. For the highest gain setting a
                        value of 25 is needed, at the moderate gain setting a value of 1 is used.
        fracleft    =>  This refers to the denominator value of the fraction of the average peak distance
                        to be used to determine the left hand bound for the fitting range. For the highest
                        gain settings a value of 2 is best for moderate gain a value of 10 is better.
        fracright   =>  This refers to the denominator value of the fraction of the average peak distance
                        to be used to determine the right hand bound for the fitting range. For the 
                        highest gain settings a value of 2 is best no change for moderate gain noted.
        usexaxis    =>  Change load histogram to use the currently stored x-axis
        lcl         =>  Set a lower centroid limit, centroids below this are discarded
        ucl         =>  Set an upper centroid limit, centroids above this are discarded
        lal         =>  Set a minimum threshold on peak area, areas below this are discarded
        pspace      =>  Set a minimum spacing between peaks to be considered double counting
     */
    
    // Initalize Function Variables
    TSpectrum *peak_off = new TSpectrum(NUM_PEAKS);
    TH1D *hist;
    TCanvas *canvas = new TCanvas("Main", "main", 1000, 1000);
    
    int exit = 0;
    int num = 0;
    int counter = 0;
    int attempt = 0;
    int histlowerlim = lcl;
    int histupperlim = ucl;
    int mnum = 0; // master number of peaks found from peak_off;
    bool autorun = true;
    bool tracker = false;
    short priortrack[NUM_PEAKS];
    char str[500];
    float energy_sigma = 0.5; // Sets the sigma value to use when doing energy domain fits
    double comp = 0;
    double avgdiff = 0;
    double parhold[NUM_PEAKS];
    double temphold[NUM_PEAKS];
    double sigmahold[NUM_PEAKS];
    double centhold[NUM_PEAKS];
    double priorcent[NUM_PEAKS];
    double openvar = 0;
    double poff_sig = *sig;
    double poff_amp = *amp;
    double dummyvar = 0;
    double slope = 0;
    double offset = 0;
    double ranleft = 0;
    double ranright = 0;
    // End
    /*
     Indicate the start of peak finding and identify partition being processed. If no peaks are found
     the program will give the user an opportunaty to change the sig and amp values, note that a
     comman failure mode is that the program fail to load the histogram properly...
    */
    std::cout << "Looking for peaks in " << pointer->cslice << "...\n";
    hist = pointer->histogram[pointer->cslice];
    hist->Sumw2(kFALSE);
    hist->Sumw2();
    num = peak_off->Search(hist, poff_sig, "", poff_amp);
    if (trouble) printf("Found %d peaks in partition %d\n", num, pointer->cslice);
    while (num < 1) {
        canvas = new TCanvas("peak_off", "peak_off", 1000, 1000);
        hist->Draw("same hist");
        canvas->Update();
        canvas->Update();
        std::cout << "The sigma and amplitude parameters appear to be no good please enter new ones...\n";
        std::cout << "Current poff_sig: " << poff_sig << "\t Current poff_amp: " << poff_amp << "\n";
        std::cout << "Enter new sigma: ";
        while (true) {
            std::cin >> poff_sig;
            if (std::cin.fail()) {
                std::cin.clear();
                std::cin.ignore();
                std::cout << "Invalid input please enter one of the available options...\n";
            }
            else {
                *sig = poff_sig;
                *amp = poff_amp;
                break;
            }
        } // Error catching for sigma
        std::cout << "Enter new amplitude: ";
        while (true) {
            std::cin >> poff_amp;
            if (std::cin.fail()) {
                std::cin.clear();
                std::cin.ignore();
                std::cout << "Invalid input please enter one of the available options...\n";
            }
            else break;
        } // Error catching for amplitude
        std::cout << "Rerunning search with poff_sig: " << poff_sig << "\t poff_amp: " << poff_amp << "\n";
        num = peak_off->Search(hist, poff_sig, "", poff_amp);
        std::cout << "After rerunning search with poff_sig: " << poff_sig << "\t poff_amp: " << poff_amp << "\n";
    } // Catches situations where peaks where not found
    pointer->current = 0;
    for (int i = 0; i < NUM_PEAKS; i++) { parhold[i] = 0; temphold[i] = 0; centhold[i] = 0;}
    for (int i = 0; i < num; i++) centhold[i] = peak_off->GetPositionX()[i];
    if (trouble) {
        printf("Printing all the peaks found by ->Search\n");
        for (int i = 0; i < num; i++) printf("%d\t%0.2f\n", i, centhold[i]);
    }
    
    sort_list(centhold, num);
    mnum = num; // set the master number of peaks
    
    /*/
     Set parameters for fitting. The centroid is set by the peak location found by TSpectrum, the initial
     sigma value (comp) is currently set at 5 bins or 0.75eV. This is set to an average value based on
     prior runs supplied value later. The initial area is take as Area = sqrt(2Pi)*sigma*Amplitude where 
     amplitude is the histogram bin value averaged across the peak location.
     
     Peaks below and above centroid limits or below area limits are removed along with duplicate peaks as
     set by pspace.
     
     Note:
     - At the end both ranleft and ranright mark the range for the fit to be performed over
     
     /*/
    //comp = (useenergy ? 0.75 : 5);
    comp = (useenergy ? 0.75 : pspace/3.0);
    //avgdiff = 1.0*(centhold[num - 1] - centhold[0])/(num - 1);
    counter = 0;
    attempt = num - 1;
    while (counter < num && centhold[counter] < lcl) counter++;
    while (attempt > 0 && centhold[attempt] > ucl) attempt--;
    avgdiff = 1.0*(centhold[attempt] - centhold[counter])/(attempt - counter);
    //printf("%f, %d, %d, %f, %f\n", avgdiff, attempt, counter, centhold[attempt], centhold[counter]);
    attempt = 0;
    counter = 0;
    
    //gBenchmark->Start("test"); //Added to try fix plotting issue
    
    if(pointer->cslice != 0){
        dummyvar = 0;
        counter = 0;
        for (int i = 0; i < NUM_PEAKS; i++){
            if (pointer->param[(pointer->cslice - 1)*NUM_PEAKS + i][0] == -1) break;
            dummyvar += pointer->param[(pointer->cslice - 1)*NUM_PEAKS + i][4];
            counter++;
        }
        if(dummyvar != 0 && counter != 0) {
            dummyvar = dummyvar/counter;
            if (dummyvar < 100) {
                printf("Prior comp: %0.2f\n", comp);
                printf("Prior avgdiff: %0.2f\n", avgdiff);
                comp = dummyvar;
                avgdiff = 1.1*3.5/pointer->calib[pointer->cslice - 1][1];
                printf("New comp: %0.2f\n", comp);
                printf("New avgdiff: %0.2f\n", avgdiff);
            }
            else printf("Prior average sigma is to large...%.2f\n", dummyvar);
        }
    } // Use prior fit to get an average value for sigma and peak spacing...
    else {
        //avgdiff = 13;
        printf("Using an avgdiff of %f and comp of %f\n", avgdiff, comp);
        printf("Initial centroids\n");
        //comp = 0.5*avgdiff;
        for (int i = 0; i < num; i++){
            printf("%d\t%0.2f\n", i, centhold[i]);
        }
    }
    for (int i = 0; i < NUM_PEAKS; i++) {
        pointer->param[pointer->cslice*NUM_PEAKS + i][2] = 0;
        pointer->param[pointer->cslice*NUM_PEAKS + i][4] = 0;
        pointer->param[pointer->cslice*NUM_PEAKS + i][6] = 0;
    }// Ensure that the allocated memory is zeroed
    
    // Perform initial guess...
    if(trouble) printf("Number of peaks before cleaning %d\n", num);
    
    for (int i = 0; i < num; i++) pointer->param[pointer->cslice*NUM_PEAKS + i][2] = centhold[i];
    if (pointer->cslice != 0) num = peak_check(pointer, num, comp, lcl, ucl, lal, pspace, trouble);
    else{
        for (int i = 0; i < num; i++) {
            dummyvar = hist->GetXaxis()->FindBin(pointer->param[pointer->cslice*NUM_PEAKS + i][2]);
            dummyvar = (hist->GetBinContent(dummyvar - 2) + hist->GetBinContent(dummyvar - 1) +
                        hist->GetBinContent(dummyvar) + hist->GetBinContent(dummyvar + 1) +
                        hist->GetBinContent(dummyvar + 2))*(usexaxis ? energy_sigma : comp)/5.0;
            pointer->param[pointer->cslice*NUM_PEAKS + i][4] = (usexaxis ? energy_sigma : comp);
            pointer->param[pointer->cslice*NUM_PEAKS + i][6] = pow(2*TMath::Pi(), 0.5)*dummyvar;
        }// Set initial parameter guesses
        counter = 0;
        for (int i = num - 1; i >= 0; i--){
            if (pointer->param[NUM_PEAKS*pointer->cslice + i][2] < lcl ||
                pointer->param[NUM_PEAKS*pointer->cslice + i][2] > ucl) {
                printf("%d\t%.2f\n", i, pointer->param[pointer->cslice*NUM_PEAKS + i][2]);
                removept(pointer, i, NUM_PEAKS*pointer->cslice, NUM_PEAKS*(pointer->cslice + 1));
                counter++;
            }
        }// Remove peaks with centroids outside of set ranges...
        num -= counter;
        counter = num - 1;
        while (pointer->param[NUM_PEAKS*pointer->cslice + counter][6] < lal){
            printf("High, Removing point: %d\tArea: %.2f\n", counter, pointer->param[NUM_PEAKS*pointer->cslice + counter][6]);
            removept(pointer, counter, NUM_PEAKS*pointer->cslice, NUM_PEAKS*(pointer->cslice + 1));
            printf("Succesfully removed\n");
            counter--;
            if (counter < 0) break;
        } // Removes peaks from the right until first peak above threshold
        num = counter + 1;
        counter = 0;
        while (pointer->param[NUM_PEAKS*pointer->cslice][6] < lal){
            printf("Low, Removing point: %d\tArea: %.2f\n", counter, pointer->param[NUM_PEAKS*pointer->cslice][6]);
            removept(pointer, 0, NUM_PEAKS*pointer->cslice, NUM_PEAKS*(pointer->cslice + 1));
            printf("Succesfully removed\n");
            counter++;
            if (counter >= num) break;
        } // Removes peaks from the left until first peak above threshold
        num -= counter;
        counter = NUM_PEAKS*pointer->cslice;
        for (int i = 1; i < num - 1; i++){
            if(pointer->param[counter + i][6] - pointer->param[counter + i - 1][6] < 0 &&
               pointer->param[counter + i][6] - pointer->param[counter + i + 1][6] < 0){
                pointer->param[counter + i][6] = pointer->param[counter + i + 1][6];
                pointer->param[counter + i][6] += pointer->param[counter + i - 1][6];
                pointer->param[counter + i][6] = pointer->param[counter + i][6]/2.0;
            }
        }// No local minima...
    }
    ranright = (int)(pointer->param[pointer->cslice*NUM_PEAKS + num - 1][2] + 1.0*avgdiff/fracright);
    ranleft = (int)(pointer->param[pointer->cslice*NUM_PEAKS][2] - 1.0*avgdiff/fracleft);
    if (ranleft < 0) ranleft = 0;
    if (ranright > tBin) ranright = tBin - 1;

    if(trouble){
        printf("Number of peaks: %d\nRange: (%.2f, %.2f)\n", num, ranleft, ranright);
        printf("Initial Values:\nPeak\tSigma\tArea\n");
        for(int i = 0; i < num; i++) {
            printf("%.2f\t%.2f\t%.2f\n",
                   pointer->param[pointer->cslice*NUM_PEAKS + i][2],
                   pointer->param[pointer->cslice*NUM_PEAKS + i][4],
                   pointer->param[pointer->cslice*NUM_PEAKS + i][6]);
        }
        
    }

    /*/
     Here the full N-Gaussian function is fit to the slice, if it fails the first time then it will
     copy over the successful fits and try again with changes to the ranges.
     /*/
    //autorun = false;
    attempt = 0;
    while (true) {
        if((autorun && ngaussfit(pointer,  num, ranleft, ranright, 0, 0, *override, useenergy)) || attempt == -100) {
            exit = 1;
            for(int i = 0; i < num; i++) if(round_decimal(pointer->param[pointer->cslice*NUM_PEAKS + i][4], 3) == 100) exit = 0;
            if (exit || attempt == -100 || *override){
                if (attempt == 0) break;
                //else if (0 < attempt && attempt <= 5 && ngaussfit(pointer,  num, ranleft, ranright, 0, 0, false, useenergy)) {
                else if (0 < attempt && attempt <= 5) {
                    std::cout << "Exiting on attempt: " << attempt << "\n";
                    std::cout << "Final N-Gaussian:\t" << pointer->cslice;
                    std::cout << "\nFit values:\nLocation\tCentroid\tdCent\tSigma\tArea\n";
                    for (int i = 0; i < num; i++) {
                        std::cout << i << "\t";
                        std::cout << pointer->param[NUM_PEAKS*pointer->cslice + i][2] << "\t";
                        std::cout << pointer->param[NUM_PEAKS*pointer->cslice + i][3] << "\t";
                        std::cout << pointer->param[NUM_PEAKS*pointer->cslice + i][4] << "\t";
                        std::cout << pointer->param[NUM_PEAKS*pointer->cslice + i][6] << "\n";
                    }
                    std::cout << "Lower Bound\tUpper Bound\n";
                    std::cout << ranleft << "\t" << ranright << "\n";
                    break;
                }
                else if ((attempt == -100 && ngaussfit(pointer, num, ranleft, ranright, 0, 0, true, useenergy))) {
                    std::cout << "Exiting with override!\n";
                    std::cout << "Final N-Gaussian:\t" << pointer->cslice;
                    std::cout << "\nFit values:\nLocation\tCentroid\tdCent\tSigma\tArea\n";
                    for (int i = 0; i < num; i++) {
                        std::cout << i << "\t";
                        std::cout << pointer->param[NUM_PEAKS*pointer->cslice + i][2] << "\t";
                        std::cout << pointer->param[NUM_PEAKS*pointer->cslice + i][3] << "\t";
                        std::cout << pointer->param[NUM_PEAKS*pointer->cslice + i][4] << "\t";
                        std::cout << pointer->param[NUM_PEAKS*pointer->cslice + i][6] << "\n";
                    }
                    std::cout << "Lower Bound\tUpper Bound\n";
                    std::cout << ranleft << "\t" << ranright << "\n";
                    break;
                }
                else if (attempt == -100){
                    printf("\n\n\nOverride Failed...\n\n\n");
                    attempt = 5;
                }
            }
        } // Exit criterion
        else if (!autorun){
            std::cout << "Enter next attempt option, enter (-1)exit\n";
            std::cout << "Options: (5)display current values, (6)ranleft, (7)ranright, ";
            std::cout << "(8)centroid, (9)sigma, (10)area, (11)histbound, (12)removept, ";
            std::cout << "(13)toggle autorun, (14)remove peaks with small areas, ";
            std::cout << "(15)addpt, (16)full value reset, (17)Force Auto, ";
            std::cout << "(18)Plot spectrum(may not work): ";
            while (true) {
                std::cin >> attempt;
                if (cin.fail()) {
                    std::cin.clear();
                    std::cin.ignore();
                    std::cout << "Invalid input please enter one of the available options...\n";
                }
                else break;
            }
            if (attempt != -100 && attempt < 0) break;
        }
        printf("Starting Attempt: %d\nOverride: %s\n", attempt, (*override ? "true" : "false"));
        if (5 == attempt || 20 < attempt){
            std::cout << "Require user input...\n";
            std::cout << "Displaying current parameters, setting auto to false...\n";
            autorun = false;
            printf("Current slice: %d out of %d\tNumber of peaks: %d\n", pointer->cslice, pointer->endslice, num);
            printf("ranleft: %.2f\tranright: %.2f\n", ranleft, ranright);
            printf("Average peak distance: %.2f\n", avgdiff);
            std::cout << "Current values:\nLocation\tCentroid\tdCent\tSigma\tArea\n";
            for (int i = 0; i < num; i++) {
                std::cout << i << "\t";
                std::cout << pointer->param[NUM_PEAKS*pointer->cslice + i][2] << "\t";
                std::cout << pointer->param[NUM_PEAKS*pointer->cslice + i][3] << "\t";
                std::cout << pointer->param[NUM_PEAKS*pointer->cslice + i][4] << "\t";
                std::cout << pointer->param[NUM_PEAKS*pointer->cslice + i][6] << "\n";
            }
        }
        else if (0 <= attempt && attempt < 5 && (useenergy || usexaxis)){
            if (pointer->cslice != 0) num = peak_check(pointer, num, comp, lcl, ucl, lal, pspace, trouble);
            else attempt = 4;
            if (attempt == 0){
                ranleft = floor(pointer->param[NUM_PEAKS*pointer->cslice][2] - 1);
                ranright = ceil(pointer->param[NUM_PEAKS*pointer->cslice + num - 1][2] + 1);
            }
            else if (attempt == 1){
                ranleft = floor(pointer->param[NUM_PEAKS*pointer->cslice][2] - 1.5);
                ranright = ceil(pointer->param[NUM_PEAKS*pointer->cslice + num - 1][2] + 1);
            }
            else if (attempt == 2){
                ranleft = floor(pointer->param[NUM_PEAKS*pointer->cslice][2] - 1.5);
                ranright = ceil(pointer->param[NUM_PEAKS*pointer->cslice + num - 1][2] + 1.5);
            }
            else if (attempt == 3){
                ranleft = floor(pointer->param[NUM_PEAKS*pointer->cslice][2] - 1);
                ranright = ceil(pointer->param[NUM_PEAKS*pointer->cslice + num - 1][2] + 1.5);
            }
            else if (attempt == 4){
                if(round_decimal(pointer->param[NUM_PEAKS*pointer->cslice + num - 1][4], 3) == 1.2){
                    removept(pointer, num - 1, NUM_PEAKS*pointer->cslice, NUM_PEAKS*(pointer->cslice + 1));
                    num--;
                }
                if(round_decimal(pointer->param[NUM_PEAKS*pointer->cslice][4], 3) == 1.2){
                    removept(pointer, 0, NUM_PEAKS*pointer->cslice, NUM_PEAKS*(pointer->cslice + 1));
                    num--;
                }
                ranleft = floor(pointer->param[NUM_PEAKS*pointer->cslice][2] - 1.25);
                ranright = ceil(pointer->param[NUM_PEAKS*pointer->cslice + num - 1][2] + 1.25);
                if (*override) attempt = -101;
            }
            attempt++;
        } // To be used when partition is in energy domain, needed????
        else if (0 <= attempt && attempt < 5){
            if (pointer->cslice != 0) num = peak_check(pointer, num, comp, lcl, ucl, lal, pspace, trouble);
            else attempt = 4;
            if (attempt == 0){
                ranleft = floor(pointer->param[NUM_PEAKS*pointer->cslice][2] - comp);
                ranright = ceil(pointer->param[NUM_PEAKS*pointer->cslice + num - 1][2] + comp);
            }
            else if (attempt == 1){
                ranleft = floor(pointer->param[NUM_PEAKS*pointer->cslice][2] - 1.5*comp);
                ranright = ceil(pointer->param[NUM_PEAKS*pointer->cslice + num - 1][2] + 1.5*comp);
            }
            else if (attempt == 2){
                ranleft = floor(pointer->param[NUM_PEAKS*pointer->cslice][2] - 0.3*avgdiff);
                ranright = ceil(pointer->param[NUM_PEAKS*pointer->cslice + num - 1][2] + 0.3*avgdiff);
            }
            else if (attempt == 3){
                ranleft = floor(pointer->param[NUM_PEAKS*pointer->cslice][2] - 0.4*avgdiff);
                ranright = ceil(pointer->param[NUM_PEAKS*pointer->cslice + num - 1][2] + 0.4*avgdiff);
            }
            else if (attempt == 4){
                ranleft = floor(pointer->param[NUM_PEAKS*pointer->cslice][2] - 0.5*avgdiff); // Was 15 before
                ranright = ceil(pointer->param[NUM_PEAKS*pointer->cslice + num - 1][2] + 0.5*avgdiff); // was 15 before
            }
            attempt++;
        } // Bypass 0 -> 4 below for simple change in range of fit...
        else if (attempt ==  0){
            num = 0;
            printf("Reset Values:\nPeak\tCentroid\tSigma\tArea\n");
            for (int i = 0; i < mnum; i++){
                dummyvar = hist->GetXaxis()->FindBin(centhold[i]);
                dummyvar = hist->GetBinContent(dummyvar)*(usexaxis ? energy_sigma : comp);
                pointer->param[pointer->cslice*NUM_PEAKS + i][2] = centhold[i];
                pointer->param[pointer->cslice*NUM_PEAKS + i][4] = (usexaxis ? energy_sigma : comp);
                pointer->param[pointer->cslice*NUM_PEAKS + i][6] = pow(2*TMath::Pi(), 0.5)*dummyvar;
                printf("%d\t%.2f\t%.2f\t%.2f\n", i, centhold[i], (usexaxis ? energy_sigma : comp), pow(2*TMath::Pi(), 0.5)*dummyvar);
                num++;
            } // Reset parameters to initial attempt values...
            printf("Removing points out of bounds...\nPeak\tCentroid\n");
            for (int i = mnum - 1; i >= 0; i--){
                if (pointer->param[NUM_PEAKS*pointer->cslice + i][2] < lcl ||
                    pointer->param[NUM_PEAKS*pointer->cslice + i][2] > ucl) {
                    printf("%d\t%.2f\n", i, pointer->param[pointer->cslice*NUM_PEAKS + i][2]);
                    removept(pointer, i, NUM_PEAKS*pointer->cslice, NUM_PEAKS*(pointer->cslice + 1));
                    num--;
                }
            }// Remove peaks with centroids outside of set ranges...
            printf("Removing points below area limits...\nPeak\tCentroid\tArea\n");
            if (true){
                counter = 0;
                while (pointer->param[NUM_PEAKS*pointer->cslice + num - 1 - counter][6] < lal){
                    printf("%.2f\t%.2f\n", pointer->param[pointer->cslice*NUM_PEAKS + num - 1 - counter][2],
                           pointer->param[pointer->cslice*NUM_PEAKS + num - 1 - counter][6]);
                    removept(pointer, num - 1 - counter, NUM_PEAKS*pointer->cslice, NUM_PEAKS*(pointer->cslice + 1));
                    counter++;
                } // Removes peaks from the right until first peak above threshold
                num -= counter;
                counter = 0;
                while (pointer->param[NUM_PEAKS*pointer->cslice + counter][6] < lal){
                    printf("%.2f\t%.2f\n", pointer->param[pointer->cslice*NUM_PEAKS + counter][2],
                           pointer->param[pointer->cslice*NUM_PEAKS + counter][6]);
                    removept(pointer, 0, NUM_PEAKS*pointer->cslice, NUM_PEAKS*(pointer->cslice + 1));
                    counter++;
                } // Removes peaks from the left until first peak above threshold
                num -= counter;
                counter = NUM_PEAKS*pointer->cslice;
                for (int i = 1; i < num - 1; i++){
                    if(pointer->param[counter + i][6] - pointer->param[counter + i - 1][6] < 0 &&
                       pointer->param[counter + i][6] - pointer->param[counter + i + 1][6] < 0){
                        pointer->param[counter + i][6] = pointer->param[counter + i + 1][6];
                        pointer->param[counter + i][6] += pointer->param[counter + i - 1][6];
                        pointer->param[counter + i][6] = pointer->param[counter + i][6]/2.0;
                    }
                }// No local minima...
            }// Remove peaks with small areas...
            printf("Set of parameters looking for in prior runs...\nPeak\tCentroid\n");
            for (int i = 0; i < num; i++){
                printf("%d\t%.2f\n", i, pointer->param[pointer->cslice*NUM_PEAKS + i][2]);
            }
            if (pointer->cslice){
                std::cout << "Attempting Prior fit values\n";
                std::cout << "Peak#\tCentroid\tSigma\tArea\tSlice\n";
                /*/
                 The i-th index is going to scan through the values in the current partitions alotted section in reverse direction.
                 The j-th index indicates which prior partition section (up to 5 partitions away) is currently being scanned.
                 The k-th index is the value in the j-th partition being compared to the i-th value.
                 tracker is used to break out from the double for loops once an acceptable canidate has been found.
                 /*/
                counter = 0;
                for (int i = NUM_PEAKS*pointer->cslice + num - 1; i >= NUM_PEAKS*pointer->cslice; i--) {
                    if (pointer->param[i][2] == 0) continue;
                    tracker = false;
                    for (int j = pointer->cslice - 1; j >= (pointer->cslice > 5 ? pointer->cslice - 5: 0); j--) {
                        for (int k = NUM_PEAKS*j; k < NUM_PEAKS*(j + 1); k++) {
                            if (pointer->param[k][2] - comp <= pointer->param[i][2] &&
                                pointer->param[i][2] <= pointer->param[k][2] + comp){
                                if (0.8*lal <= pointer->param[k][6]){
                                    printf("Replaceing line:\nCent: %.2f\tSigma: %.2f\tArea: %.2f\n", pointer->param[i][2], pointer->param[i][4], pointer->param[i][6]);
                                    printf("New line:\nCent: %.2f\tSigma: %.2f\tArea: %.2f\n", pointer->param[k][2], pointer->param[k][4], pointer->param[k][6]);
                                    pointer->param[i][2] = pointer->param[k][2];
                                    pointer->param[i][4] = pointer->param[k][4];
                                    pointer->param[i][6] = pointer->param[k][6];
                                    tracker = true;
                                    break;
                                }
                                else if (j == (pointer->cslice > 5 ? pointer->cslice - 5: 0)){
                                    printf("Removing line:\nCent: %.2f\tSigma: %.2f\tArea: %.2f\n", pointer->param[i][2], pointer->param[i][4], pointer->param[i][6]);
                                    removept(pointer, i - NUM_PEAKS*pointer->cslice, NUM_PEAKS*pointer->cslice, NUM_PEAKS*(pointer->cslice + 1));
                                    counter++;
                                    tracker = true;
                                    break;
                                }
                            }
                            else if (j == (pointer->cslice > 5 ? pointer->cslice - 5: 0) &&
                                     k == NUM_PEAKS*(j + 1) - 1 && i != NUM_PEAKS*pointer->cslice){
                                printf("Changing entry values for line:\nCent: %.2f\tSigma: %.2f\tArea: %.2f\n", pointer->param[i][2], pointer->param[i][4], pointer->param[i][6]);
                                printf("New entry values for line:\nCent: %.2f\tSigma: %.2f\tArea: %.2f\n", pointer->param[i][2], pointer->param[i-1][4], pointer->param[i-1][6]);
                                pointer->param[i][4] = pointer->param[i - 1][4];
                                pointer->param[i][6] = pointer->param[i - 1][6];
                            }
                        }
                        if (tracker) break;
                    }
                }
                num -= counter;
            }
            ranleft = (int)(pointer->param[pointer->cslice*NUM_PEAKS][2] - 0.3*avgdiff);
            ranright = (int)(pointer->param[pointer->cslice*NUM_PEAKS + num - 1][2] + 0.3*avgdiff);
            attempt = 0;
            attempt++;
        }   // Max left side bound (Edited for Pu to remove unknown low energy peak...)/*/
        else if (attempt ==  1){
            num = 0;
            for (int i = 0; i < mnum; i++){
                dummyvar = hist->GetXaxis()->FindBin(centhold[i]);
                dummyvar = hist->GetBinContent(dummyvar)*(usexaxis ? energy_sigma : comp);
                pointer->param[pointer->cslice*NUM_PEAKS + i][2] = centhold[i];
                pointer->param[pointer->cslice*NUM_PEAKS + i][4] = (usexaxis ? energy_sigma : comp);
                pointer->param[pointer->cslice*NUM_PEAKS + i][6] = pow(2*TMath::Pi(), 0.5)*dummyvar;
                num++;
            }// Reset parameters to initial attempt values...
            for (int i = num - 1; i >= 0; i--){
                if (pointer->param[NUM_PEAKS*pointer->cslice + i][2] < lcl ||
                    pointer->param[NUM_PEAKS*pointer->cslice + i][2] > ucl) {
                    removept(pointer, i, NUM_PEAKS*pointer->cslice, NUM_PEAKS*(pointer->cslice + 1));
                    num--;
                }
            }// Remove peaks with centroids outside of set ranges...
            if (true){
                counter = 0;
                while (pointer->param[NUM_PEAKS*pointer->cslice + num - 1 - counter][6] < lal){
                    removept(pointer, num - 1 - counter, NUM_PEAKS*pointer->cslice, NUM_PEAKS*(pointer->cslice + 1));
                    counter++;
                } // Removes peaks from the right until first peak above threshold
                num -= counter;
                counter = 0;
                while (pointer->param[NUM_PEAKS*pointer->cslice + counter][6] < lal){
                    removept(pointer, 0, NUM_PEAKS*pointer->cslice, NUM_PEAKS*(pointer->cslice + 1));
                    counter++;
                } // Removes peaks from the left until first peak above threshold
                num -= counter;
                counter = NUM_PEAKS*pointer->cslice;
                for (int i = 1; i < num - 1; i++){
                    if(pointer->param[counter + i][6] - pointer->param[counter + i - 1][6] < 0 &&
                       pointer->param[counter + i][6] - pointer->param[counter + i + 1][6] < 0){
                        pointer->param[counter + i][6] = pointer->param[counter + i + 1][6];
                        pointer->param[counter + i][6] += pointer->param[counter + i - 1][6];
                        pointer->param[counter + i][6] = pointer->param[counter + i][6]/2.0;
                    }
                }// No local minima...
            }// Remove peaks with small areas...
            if (pointer->cslice){
                std::cout << "Attempting Prior fit values\n";
                std::cout << "Peak#\tCentroid\tSigma\tArea\tSlice\n";
                /*/
                 The i-th index is going to scan through the values in the current partitions alotted section in reverse direction.
                 The j-th index indicates which prior partition section (up to 5 partitions away) is currently being scanned.
                 The k-th index is the value in the j-th partition being compared to the i-th value.
                 tracker is used to break out from the double for loops once an acceptable canidate has been found.
                 /*/
                counter = 0;
                for (int i = NUM_PEAKS*pointer->cslice + num - 1; i >= NUM_PEAKS*pointer->cslice; i--) {
                    if (pointer->param[i][2] == 0) continue;
                    tracker = false;
                    for (int j = pointer->cslice - 1; j >= (pointer->cslice > 5 ? pointer->cslice - 5: 0); j--) {
                        for (int k = NUM_PEAKS*j; k < NUM_PEAKS*(j + 1); k++) {
                            if (pointer->param[k][2] - comp <= pointer->param[i][2] &&
                                pointer->param[i][2] <= pointer->param[k][2] + comp){
                                if (0.8*lal <= pointer->param[k][6]){
                                    pointer->param[i][2] = pointer->param[k][2];
                                    pointer->param[i][4] = pointer->param[k][4];
                                    pointer->param[i][6] = pointer->param[k][6];
                                    tracker = true;
                                    break;
                                }
                                else if (j == (pointer->cslice > 5 ? pointer->cslice - 5: 0)){
                                    removept(pointer, i - NUM_PEAKS*pointer->cslice, NUM_PEAKS*pointer->cslice, NUM_PEAKS*(pointer->cslice + 1));
                                    counter++;
                                    tracker = true;
                                    break;
                                }
                            }
                            else if (j == (pointer->cslice > 5 ? pointer->cslice - 5: 0) &&
                                     k == NUM_PEAKS*(j + 1) - 1 && i != NUM_PEAKS*pointer->cslice){
                                pointer->param[i][4] = pointer->param[i - 1][4];
                                pointer->param[i][6] = pointer->param[i - 1][6];
                            }
                        }
                        if (tracker) break;
                    }
                }
                num -= counter;
            }
            ranleft = (int)(pointer->param[pointer->cslice*NUM_PEAKS][2] - 13);
            ranright = (int)(pointer->param[pointer->cslice*NUM_PEAKS + num - 1][2] + 13);
            attempt = 1;
            attempt++;
        }   // Update fit values to reflect prior fits for that peak.
        else if (attempt ==  2){
            if (true){
                counter = 0;
                while (pointer->param[NUM_PEAKS*pointer->cslice + num - 1 - counter][6] < lal){
                    removept(pointer, num - 1 - counter, NUM_PEAKS*pointer->cslice, NUM_PEAKS*(pointer->cslice + 1));
                    counter++;
                } // Removes peaks from the right until first peak above threshold
                num -= counter;
                counter = 0;
                while (pointer->param[NUM_PEAKS*pointer->cslice + counter][6] < lal){
                    removept(pointer, 0, NUM_PEAKS*pointer->cslice, NUM_PEAKS*(pointer->cslice + 1));
                    counter++;
                } // Removes peaks from the left until first peak above threshold
                num -= counter;
                counter = NUM_PEAKS*pointer->cslice;
                for (int i = 1; i < num - 1; i++){
                    if(pointer->param[counter + i][6] - pointer->param[counter + i - 1][6] < 0 &&
                       pointer->param[counter + i][6] - pointer->param[counter + i + 1][6] < 0){
                        pointer->param[counter + i][6] = pointer->param[counter + i + 1][6];
                        pointer->param[counter + i][6] += pointer->param[counter + i - 1][6];
                        pointer->param[counter + i][6] = pointer->param[counter + i][6]/2.0;
                    }
                }// No local minima...
            }// Remove peaks with small areas...
            counter = 0;
            for (int i = num - 1; i > 0; i--){
                if (0 < pointer->param[NUM_PEAKS*pointer->cslice + i][2] - pointer->param[NUM_PEAKS*pointer->cslice + i - 1][2] &&
                    pointer->param[NUM_PEAKS*pointer->cslice + i][2] - pointer->param[NUM_PEAKS*pointer->cslice + i - 1][2] < pspace) { // This was set to 75 for the Dec 2015 pu239
                    removept(pointer, i, NUM_PEAKS*pointer->cslice, NUM_PEAKS*(pointer->cslice + 1));
                    counter--;
                }
            }// This should remove duplicate peaks...
            num -= counter;
            ranleft = (int)(pointer->param[pointer->cslice*NUM_PEAKS][2] - 1.2*comp);
            ranright = (int)(pointer->param[pointer->cslice*NUM_PEAKS + num - 1][2] + 1.2*comp);
            attempt = 2;
            attempt++;
        }   // Send to remove all subpar peaks
        else if (attempt ==  3){
            ranleft = (int)(pointer->param[pointer->cslice*NUM_PEAKS][2] - comp);
            ranright = (int)(pointer->param[NUM_PEAKS*pointer->cslice + num - 1][2] + comp);
            attempt = 3;
            attempt++;
        }   // Decrease fit range...
        else if (attempt ==  4){
            printf("Performing full reset...\n");
            num = 0;
            for (int i = 0; i < mnum; i++){
                dummyvar = hist->GetXaxis()->FindBin(centhold[i]);
                dummyvar = hist->GetBinContent(dummyvar)*(usexaxis ? energy_sigma : comp);
                pointer->param[pointer->cslice*NUM_PEAKS + i][2] = centhold[i];
                pointer->param[pointer->cslice*NUM_PEAKS + i][4] = (usexaxis ? energy_sigma : comp);
                pointer->param[pointer->cslice*NUM_PEAKS + i][6] = pow(2*TMath::Pi(), 0.5)*dummyvar;
                num++;
            }// Reset parameters to initial attempt values...
            for (int i = mnum - 1; i >= 0; i--){
                if (pointer->param[NUM_PEAKS*pointer->cslice + i][2] < lcl ||
                    pointer->param[NUM_PEAKS*pointer->cslice + i][2] > ucl) {
                    removept(pointer, i, NUM_PEAKS*pointer->cslice, NUM_PEAKS*(pointer->cslice + 1));
                    num--;
                }
            }// Remove peaks with centroids outside of set ranges...
            if (true){
                counter = 0;
                while (pointer->param[NUM_PEAKS*pointer->cslice + num - 1 - counter][6] < lal){
                    removept(pointer, num - 1 - counter, NUM_PEAKS*pointer->cslice, NUM_PEAKS*(pointer->cslice + 1));
                    counter++;
                } // Removes peaks from the right until first peak above threshold
                num -= counter;
                counter = 0;
                while (pointer->param[NUM_PEAKS*pointer->cslice + counter][6] < lal){
                    removept(pointer, 0, NUM_PEAKS*pointer->cslice, NUM_PEAKS*(pointer->cslice + 1));
                    counter++;
                } // Removes peaks from the left until first peak above threshold
                num -= counter;
                counter = NUM_PEAKS*pointer->cslice;
                for (int i = 1; i < num - 1; i++){
                    if(pointer->param[counter + i][6] - pointer->param[counter + i - 1][6] < 0 &&
                       pointer->param[counter + i][6] - pointer->param[counter + i + 1][6] < 0){
                        pointer->param[counter + i][6] = pointer->param[counter + i + 1][6];
                        pointer->param[counter + i][6] += pointer->param[counter + i - 1][6];
                        pointer->param[counter + i][6] = pointer->param[counter + i][6]/2.0;
                    }
                }// No local minima...
            }// Remove peaks with small areas...
            if (pointer->cslice){
                std::cout << "Attempting Prior fit values\n";
                std::cout << "Peak#\tCentroid\tSigma\tArea\tSlice\n";
                /*/
                 The i-th index is going to scan through the values in the current partitions alotted section in reverse direction.
                 The j-th index indicates which prior partition section (up to 5 partitions away) is currently being scanned.
                 The k-th index is the value in the j-th partition being compared to the i-th value.
                 tracker is used to break out from the double for loops once an acceptable canidate has been found.
                 /*/
                counter = 0;
                for (int i = NUM_PEAKS*pointer->cslice + num - 1; i >= NUM_PEAKS*pointer->cslice; i--) {
                    if (pointer->param[i][2] == 0) continue;
                    tracker = false;
                    for (int j = pointer->cslice - 1; j >= (pointer->cslice > 5 ? pointer->cslice - 5: 0); j--) {
                        for (int k = NUM_PEAKS*j; k < NUM_PEAKS*(j + 1); k++) {
                            if (pointer->param[k][2] - comp <= pointer->param[i][2] &&
                                pointer->param[i][2] <= pointer->param[k][2] + comp){
                                if (lal < pointer->param[k][6]){
                                    pointer->param[i][2] = pointer->param[k][2];
                                    pointer->param[i][4] = pointer->param[k][4];
                                    pointer->param[i][6] = pointer->param[k][6];
                                    tracker = true;
                                    break;
                                }
                                else if (j == (pointer->cslice > 5 ? pointer->cslice - 5: 0)){
                                    removept(pointer, i - NUM_PEAKS*pointer->cslice, NUM_PEAKS*pointer->cslice, NUM_PEAKS*(pointer->cslice + 1));
                                    counter++;
                                    tracker = true;
                                    break;
                                }
                            }
                            else if (j == (pointer->cslice > 5 ? pointer->cslice - 5: 0) &&
                                     k == NUM_PEAKS*(j + 1) - 1 && i != NUM_PEAKS*pointer->cslice){
                                pointer->param[i][4] = pointer->param[i - 1][4];
                                pointer->param[i][6] = pointer->param[i - 1][6];
                            }
                        }
                        if (tracker) break;
                    }
                }
                num -= counter;
            }
            ranleft = (int)(pointer->param[pointer->cslice*NUM_PEAKS][2] - 15);
            ranright = (int)(pointer->param[pointer->cslice*NUM_PEAKS + num - 1][2] + 15);
            attempt = 4;
            if (*override) attempt = -101;
            attempt++;
        }   // Full reset
        else if (attempt ==  6){
            std::cout << "\nEnter new lower bound: ";
            while (true) {
                std::cin >> dummyvar;
                if (cin.fail()) {
                    std::cin.clear();
                    std::cin.ignore();
                    printf("Invalid input please enter one of the available options...\n");
                }
                else if (dummyvar < ranright) {
                    ranleft = dummyvar;
                    break;
                }
                else printf("Entry was above ranright...\n");
            }
            attempt = 5;
        }   // Change ranleft
        else if (attempt ==  7){
            std::cout << "\nEnter new upper bound: ";
            while (true) {
                std::cin >> dummyvar;
                if (cin.fail()) {
                    std::cin.clear();
                    std::cin.ignore();
                    printf("Invalid input please enter one of the available options...\n");
                }
                else if (ranleft < dummyvar) {
                    ranright = dummyvar;
                    break;
                }
                else printf("Entry was below ranleft...\n");
            }
            
            attempt = 5;
        }   // Change ranright
        else if (attempt ==  8){
            std::cout << "\nEnter line you wish to change: ";
            while (true) {
                std::cin >> attempt;
                if (cin.fail()) {
                    std::cin.clear();
                    std::cin.ignore();
                    std::cout << "Invalid input please enter one of the available options...\n";
                }
                else break;
            }
            std::cout << "\nEnter new centroid value: ";
            if (0 <= attempt && attempt <= num) while (true) {
                std::cin >> pointer->param[NUM_PEAKS*pointer->cslice + attempt][2];
                if (cin.fail()) {
                    std::cin.clear();
                    std::cin.ignore();
                    std::cout << "Invalid input please enter one of the available options...\n";
                }
                else break;
            }
            attempt = 5;
        }   // Change Centroid
        else if (attempt ==  9){
            std::cout << "\nEnter line you wish to change: ";
            while (true) {
                std::cin >> attempt;
                if (cin.fail()) {
                    std::cin.clear();
                    std::cin.ignore();
                    std::cout << "Invalid input please enter one of the available options...\n";
                }
                else break;
            }
            std::cout << "\nEnter new sigma value: ";
            if (0 <= attempt && attempt <= num) while (true) {
                std::cin >> pointer->param[NUM_PEAKS*pointer->cslice + attempt][4];
                if (cin.fail()) {
                    std::cin.clear();
                    std::cin.ignore();
                    std::cout << "Invalid input please enter one of the available options...\n";
                }
                else break;
            }
                
                
            attempt = 5;
        }   // Change Sigma
        else if (attempt == 10){
            std::cout << "\nEnter line you wish to change: ";
            while (true) {
                std::cin >> attempt;
                if (cin.fail()) {
                    std::cin.clear();
                    std::cin.ignore();
                    std::cout << "Invalid input please enter one of the available options...\n";
                }
                else break;
            }
            std::cout << "\nEnter new Area value: ";
            if (0 <= attempt && attempt <= num) while (true) {
                std::cin >> pointer->param[NUM_PEAKS*pointer->cslice + attempt][6];
                if (cin.fail()) {
                    std::cin.clear();
                    std::cin.ignore();
                    std::cout << "Invalid input please enter one of the available options...\n";
                }
                else break;
            }
            attempt = 5;
        }  // Change Area
        else if (attempt == 11){
            std::cout << "Enter new lower bound as a negative and upper bound as positive: ";
            while (true) {
                std::cin >> attempt;
                if (cin.fail()) {
                    std::cin.clear();
                    std::cin.ignore();
                    std::cout << "Invalid input please enter one of the available options...\n";
                }
                else break;
            }
            if (attempt <= 0 && -attempt < histupperlim) histlowerlim = -attempt;
            else if (histlowerlim < attempt) histupperlim = attempt;
            attempt = 5;
        }  // Change histogram bounds
        else if (attempt == 12){
            std::cout << "\nWARNING YOU ARE ABOUT TO REMOVE A POINT FROM THE FIT!!!!!\n";
            std::cout << "Please enter the peak that you wish to remove from fitting (Enter -1 to abort): ";
            while (true) {
                std::cin >> attempt;
                if (cin.fail()) {
                    std::cin.clear();
                    std::cin.ignore();
                    std::cout << "Invalid input please enter one of the available options...\n";
                }
                else break;
            }
            if (0 <= attempt && attempt < num) {
                removept(pointer, attempt, NUM_PEAKS*pointer->cslice, NUM_PEAKS*(pointer->cslice + 1));
                num = num - 1;
            }
            else std::cout << "\nNo peak was removed\n";
            attempt = 5;
        }  // Remove point
        else if (attempt == 13){
            printf("Toogling auto from %s to %s\n", (autorun ? "true" : "false"), (!autorun ? "true" : "false"));
            autorun = !autorun;
            attempt = 5;
        }  // Toggle auto
        else if (attempt == 14){
            for (int i = num - 1; i >= 0; i--) if (pointer->param[NUM_PEAKS*pointer->cslice + i][6] < lal) {
                removept(pointer, i, NUM_PEAKS*pointer->cslice, NUM_PEAKS*(pointer->cslice + 1));
                num--;
            } // Remove peaks with small areas...
            for (int i = num - 1; i > 0; i--){
                if (0 < pointer->param[NUM_PEAKS*pointer->cslice + i][2] - pointer->param[NUM_PEAKS*pointer->cslice + i - 1][2] &&
                    pointer->param[NUM_PEAKS*pointer->cslice + i][2] - pointer->param[NUM_PEAKS*pointer->cslice + i - 1][2] < pspace) { // This was set to 75 for the Dec 2015 pu239
                    removept(pointer, i, NUM_PEAKS*pointer->cslice, NUM_PEAKS*(pointer->cslice + 1));
                    num--;
                }
            }// This should remove duplicate peaks...
            ranleft = (int)(pointer->param[pointer->cslice*NUM_PEAKS][2] - 1.2*comp);
            ranright = (int)(pointer->param[pointer->cslice*NUM_PEAKS + num - 1][2] + 1.2*comp);
            attempt = 5;
        }  // Remove points w/ small area
        else if (attempt == 15){
            std::cout << "Entering new point\n";
            std::cout << "Enter new centroid value: ";
            std::cin >> pointer->param[NUM_PEAKS*pointer->cslice + num][2];
            pointer->param[NUM_PEAKS*pointer->cslice + num][4] = pointer->param[NUM_PEAKS*pointer->cslice + num - 1][4];
            pointer->param[NUM_PEAKS*pointer->cslice + num][6] = pointer->param[NUM_PEAKS*pointer->cslice + num - 1][6];
            num++;
            attempt = 5;
        }  // Add point
        else if (attempt == 16){
            std::cout << "Performing full reset...\n";
            num = 0;
            printf("Reset Values:\nPeak\tCentroid\tSigma\tArea\n");
            for (int i = 0; i < mnum; i++){
                dummyvar = hist->GetXaxis()->FindBin(centhold[i]);
                dummyvar = hist->GetBinContent(dummyvar)*(usexaxis ? energy_sigma : comp);
                pointer->param[pointer->cslice*NUM_PEAKS + i][2] = centhold[i];
                pointer->param[pointer->cslice*NUM_PEAKS + i][4] = (usexaxis ? energy_sigma : comp);
                pointer->param[pointer->cslice*NUM_PEAKS + i][6] = pow(2*TMath::Pi(), 0.5)*dummyvar;
                printf("%d\t%.2f\t%.2f\t%.2f\n", i, centhold[i], (usexaxis ? energy_sigma : comp), pow(2*TMath::Pi(), 0.5)*dummyvar);
                num++;
            } // Reset parameters to initial attempt values...
            printf("Removing points out of bounds...\nPeak\tCentroid\n");
            for (int i = mnum - 1; i >= 0; i--){
                if (pointer->param[NUM_PEAKS*pointer->cslice + i][2] < lcl ||
                    pointer->param[NUM_PEAKS*pointer->cslice + i][2] > ucl) {
                    printf("%d\t%.2f\n", i, pointer->param[pointer->cslice*NUM_PEAKS + i][2]);
                    removept(pointer, i, NUM_PEAKS*pointer->cslice, NUM_PEAKS*(pointer->cslice + 1));
                    num--;
                }
            }// Remove peaks with centroids outside of set ranges...
            ranleft = (int)(pointer->param[pointer->cslice*NUM_PEAKS][2] - 10);
            ranright = (int)(pointer->param[pointer->cslice*NUM_PEAKS + num - 1][2] + 10);
            attempt = 5;
        }  // Full reset
        else if (attempt == 17){
            printf("Confirm that you wish to exit user edit mode and automate the rest... 0 is exit no change 1 is switch override on\n");
            while (true) {
                std::cin >> attempt;
                if (cin.fail()) {
                    std::cin.clear();
                    std::cin.ignore();
                    std::cout << "Invalid input please enter one of the available options...\n";
                }
                else break;
            }
            if (attempt) {
                autorun = true;
                *override = autorun;
            }
            attempt = 5;
        } // Set to forced auto
        else if (attempt == 18){
            sprintf(str, "%s%s%s%d%s%d%s", "File: ", pointer->name, " STJ # ", pointer->channel, " Slice: ", pointer->cslice, " N-Gaussian");
            canvas = new TCanvas(str, str, 1000, 1000);
            hist->Draw("hist");
            hist->GetXaxis()->SetRangeUser(histlowerlim, histupperlim);
            canvas->Modified();
            canvas->Update();
            hist->Draw("same");
            //gBenchmark->Show("test");
        }// Plot data...
        else if (attempt > -100) attempt++;
        else if (attempt < -100) attempt = 5;
    }
    if (trouble){
        std::cout << "Copy over data and exiting...\n";
        std::cout << "Peaks from fitter for slice: " << pointer->cslice << "\n";
        std::cout << "Location\tCentroid\tdCent\tSigma\tArea\n";
        for (int i = 0; i < num; i++) {
            std::cout << i << "\t";
            std::cout << pointer->param[NUM_PEAKS*pointer->cslice + i][2] << "\t";
            std::cout << pointer->param[NUM_PEAKS*pointer->cslice + i][3] << "\t";
            std::cout << pointer->param[NUM_PEAKS*pointer->cslice + i][4] << "\t";
            std::cout << pointer->param[NUM_PEAKS*pointer->cslice + i][6] << "\n";
        }
    }
    exit = calibration(pointer, lal, false, trouble); // Run initial calibration on peaks for use in next partition
    return exit;
}

/* **************************************************************************************************** */
int energy(info*pointer, int val = 10, int numslice = 1000){
/* **************************************************************************************************** */
    /*
        This function will convert the bin value to energy for all entries in pointer. Because the data
        parsed as a moving average of 10 windows the first and last 5 thyme values will not be converted
        these entries will come in as a 0 instead to distinctly mark them.
     
        update:
        val is to be used in connection with the val of the partitioning function to get the right center
        for the window...
     
        numslice indicates the total number of expected output slices for example if there are 20
        partitions and the window is 5 then the first (and last) 4 will not have enough calibration points
        either to the left (or right) and need to be 'dropped'. As such there are only 22 output
        partition calibrations for partitions [4, 26) but not for [0, 4) or [26, 30). In this example
        val -> 5 and numslice = 30 - 2*(5 - 1) = 22.
     */
    // Initialize Function Variables
    int exit = 0;
    int counter = 0;
    int segment = val;
    char str[500];
    double aslope[100000];
    double daslo[100000];
    double aoffset[100000];
    double daoff[100000];
    double slope[100];
    double dslo[100];
    double offset[100];
    double doff[100];
    ofstream outfile;
    // End
    
    
    std::cout << "Starting energy Calibration\n";
    /*/
     Calculate the weighted average for the slope and offset of the individual partitions, again looking
     only at partitions with maximum windows. Equations: (for offset m -> b)
     m-bar = sum(m_i/dm_i**2)/sum(1/dm_i**2)
     dm-bar = (sum(1/dm_i**2))**(-1/2)
     
     Calibration fits : {{fit_number, slope, dslope, offset, doffset}_0, ... }
    /*/
    for (int i = 0; i < numslice; i++) {
        for (int j = i; j < i + segment; j++) {
            slope[j - i] = pointer->calib[j][1];
            dslo[j - i] = 1.0/pow(pointer->calib[j][2], 2);
            offset[j - i] = pointer->calib[j][3];
            doff[j - i] = 1.0/pow(pointer->calib[j][4], 2);
        }
        aslope[i] = 0;
        for (int j = 0; j < segment; j++) aslope[i] += slope[j]*dslo[j];
        aslope[i] = aslope[i]/sum_array(dslo, segment + 1);
        daslo[i] = pow(sum_array(dslo, segment + 1), -0.5);
        
        aoffset[i] = 0;
        for (int j = 0; j < segment; j++) aoffset[i] += offset[j]*doff[j];
        aoffset[i] = aoffset[i]/sum_array(doff, segment + 1);
        daoff[i] = pow(sum_array(doff, segment + 1), -0.5);
        /*/ Trouble Shooting
        if (i <= 20) {
            std::cout << "\n";
            std::cout << pow(sum_array(dslo, segment + 1), -0.5) << "\n";
            std::cout << daslo[i] << "\n";
            std::cout << pow(sum_array(doff, segment + 1), -0.5) << "\n";
            std::cout << &daoff[i] << "\n";
            std::cout << "\nSlope-bar: " << aslope[i] << " +/- " << daslo[i];
            std::cout << "\tOffset-Bar: " << aoffset[i] << " +/- " << daoff[i] << "\n";
            std::cout << "Number\tm\tdm\tb\tdb\n";
            for (int k = 0; k <= segment; k++) {
                std::cout << k + i << "\t" << slope[k] << "\t"  << dslo[k] << "\t" << offset[k] << "\t";
                std::cout << doff[k] << "\n";
            }
        }
        //*/
        
    }
    /*/
     Save the calibrations values used in each partition including the partition bounds.
     >> /Path/<File>_Partition_Calibration_Values.txt
    /*/
    snprintf(str, sizeof(str), "%s%s_Ch_%d_09_Partition_Calibration_Values.txt", pointer->path,
             pointer->name, pointer->channel);
    outfile.open(str);
    outfile << "000_Time Slice\t001_Slope (eV/Ch)\t002_dSlope (eV/Ch)\t003_Offset (eV)\t004_dOffset (eV)\n";
    for (int i = 0; i < numslice; i++) {
        outfile << i + segment - 1 << "\t";
        outfile << aslope[i] << "\t" << daslo[i] << "\t";
        outfile << aoffset[i] << "\t" << daoff[i] << "\n";
    }
    outfile.close();
    return exit;
}

/* **************************************************************************************************** */
int xveto(info* pointer, info* secondary, int type = 0, double dt = 1, double range = 10000){
/* **************************************************************************************************** */
    /*/ This function is intended to veto events based on the type number provided
        Input parameters:
        pointer -> This is the main file in which vetoing is to occur
        type    -> Indicates the type of vetoing that will be performed
        dt      -> Time after a specific event (HE for example) to veto in milliseconds
        range   -> Time +/- in which two events in different detectors are considered the same in ns
     
        A veto value of:
        true    -> Use this point
        false   -> Do not use this point
    /*/
    // Initialize Function Variables
    int exit = 1;
    int counter = 0;
    int one, two = 0;
    int location = 0;
    int high[2][1000000];
    bool action = false;
    // End
    
    //std::cout << "Starting xveto type: " << type << "\n";
    dt = dt*pow(10, 6); // Convert to ns
    if (type == -4){ for (int i = 0; i < pointer->last; i++) pointer->recoil[i] = false;
    }   // Set all recoil False
    else if (type == -3){ for (int i = 0; i < pointer->last; i++) pointer->high[i] = false;
    }   // Set all high False
    else if (type == -2) {
        for (int i = 0; i < pointer->last; i++){
            pointer->veto[i] = false;
            pointer->high[i] = false;
            pointer->recoil[i] = false;
        }
    }  // Set all veto False
    else if (type == -1){ // Unveto all events
        for (int i = 0; i < pointer->last; i++){
            pointer->veto[i] = true;
            pointer->high[i] = true;
            pointer->recoil[i] = true;
        }
    }   // Set all veto True
    else if (type == 0){
        for (int i = 0; i < pointer->last; i++) if (pointer->tag[i] == 1) {
            pointer->veto[i] = false;
            pointer->id[i] = -1;
        }
    }    // Set tag == 1 (Laser) veto false and id to -1
    else if (type == 1){
        /*/ In this veto I assume that there exist two time sets, t_i -> pointer and tau_j ->secondary
         I will use this notation to outline what I am doing below. /*/
        counter = 0;
        secondary->current = 1;
        for (int i = 0; i < pointer->last; i++){
            /*/ The idea here is to get a bounding time:
                                    tau_j-1 <= t_i <= tau_j
             The first if-statement will skip all tagged events
             The first while-loop will make sure that tau_j-1 is less than t_i
             The second while-loop will make sure that tau_j is greater than t_i
             The issue is that a situation may arise where j = 0 or j_max is either a lower or upper
             bound this is dealt with by adding the additional criterion to each of the last two
             if-statements to check if the value they would check is beyond the data set bounds./*/
            if(pointer->tag[i] == 1) continue;
            while (1 <= secondary->current &&
                   pointer->time[i] < secondary->time[secondary->current - 1]) secondary->current--;
            while (secondary->current < secondary->last &&
                   secondary->time[secondary->current] < pointer->time[i]) secondary->current++;
            if (secondary->current < secondary->last  &&
                secondary->tag[secondary->current] != 1 &&
                (8000 < secondary->bin[secondary->current] || 8000 < pointer->bin[i]) &&
                pointer->time[i] - range <= secondary->time[secondary->current] &&
                secondary->time[secondary->current] <= pointer->time[i] + range) {
                pointer->high[i] = false;
                secondary->high[secondary->current] = false;
                counter++;
            }
            else if (0 < secondary->current &&
                     secondary->tag[secondary->current - 1] != 1 &&
                     (8000 < secondary->bin[secondary->current - 1] || 8000 < pointer->bin[i]) &&
                     pointer->time[i] - range <= secondary->time[secondary->current - 1] &&
                     secondary->time[secondary->current - 1] <= pointer->time[i] + range) {
                pointer->high[i] = false;
                secondary->high[secondary->current - 1] = false;
                counter++;
            }
        }
        //std::cout << "\nPoints vetoed: " << counter << " out of " << pointer->last << "\n";
    }    // Veto all events with a coincidence event in another detector with at least one of the events greater than 8000
    else if (type == 2){
        for (int i = 0; i < pointer->last; i++) pointer->id[i] = 1;
    }    // Set ID array true (1) for all entries
    else if (type == 3){
        for (int i = 0; i < pointer->last; i++) pointer->id[i] = 0;
    }    // Set ID array false (0) for all entries
    else if (type == 4){
        /*/ In this veto I assume that there exist two time sets, t_i -> pointer and tau_j ->secondary
            I will use this notation to outline what I am doing below. /*/
        counter = 0;
        secondary->current = 1;
        for (int i = 0; i < pointer->last; i++){
            /*/ The idea here is to get a bounding time:
                        tau_j-1 <= t_i <= tau_j
                The first if-statement will skip all uncalibrated points and tagged events
                The first while-loop will make sure that tau_j-1 is less than t_i
                The second while-loop will make sure that tau_j is greater than t_i
                The issue is that a situation may arise where j = 0 or j_max is either a lower or upper
                bound this is dealt with by adding the additional criterion to each of the last two 
                if-statements to check if the value they would check is beyond the data set bounds./*/
            if(pointer->tag[i] == 1) continue;
            while (1 <= secondary->current &&
                   pointer->time[i] < secondary->time[secondary->current - 1]) secondary->current--;
            while (secondary->current < secondary->last &&
                   secondary->time[secondary->current] < pointer->time[i]) secondary->current++;
            if (secondary->current < secondary->last  && secondary->tag[secondary->current] != 1 &&
                pointer->time[i] - range <= secondary->time[secondary->current] &&
                secondary->time[secondary->current] <= pointer->time[i] + range) {
                pointer->veto[i] = false;
                secondary->veto[secondary->current] = false;
                counter++;
                }
            else if (0 < secondary->current && secondary->tag[secondary->current - 1] != 1 &&
                pointer->time[i] - range <= secondary->time[secondary->current - 1] &&
                secondary->time[secondary->current - 1] <= pointer->time[i] + range) {
                pointer->veto[i] = false;
                secondary->veto[secondary->current - 1] = false;
                counter++;
                }
        }
        //std::cout << "\nPoints vetoed: " << counter << " out of " << pointer->last << "\n";
    }    // Veto all events with a coincidence event in another detector regardless of bin value (Fully debugged)
    else if (type == 5){
        pointer->current = 0;
        for (int i = 0; i < pointer->last; i++) {
            /*/ The 8000 is arbitrary and is intended to reset the vetoing if there is a high energy event
                in the event that it was not vetoed in coincidence but in close proximity to a coincidence
                event/*/
            if (pointer->tag[i] == 1) continue; // All laser events are ignored
            else if (pointer->veto[i] == false ) {
                action = true;
                location = i;
            }
            else if (action &&  pointer->bin[i] < 8000 && pointer->time[i] <= pointer->time[location] + dt) {
                pointer->veto[i] = false;
            }
            else action = false;
        }
    }    // Veto all events after time dt of a coincidence event (use after type 4)
    else if (type == 6){
        /*/ This will passively 'enable' events after high energy events 8000 < bin (arbitrarily set)
            note that it requires that the prior process be type 4 or type 5. /*/
        pointer->current = 0;
        for(int i = 0; i < pointer->last; i++){
            if (pointer->tag[i] == 1) continue; // All laser events are ignored
            else if (pointer->veto[i] == false) action = false;
            else if (8000 <= pointer->bin[i]) {
                pointer->current = i;
                action = true;
            }
            else if (action && pointer->time[i] <= pointer->time[pointer->current] + dt){
                pointer->veto[i] = true;
            }
            else pointer->veto[i] = false;
        }
    }    // Veto all events except for time dt after none coincidence HE event (use after type 4 or 5)
    return exit;
}

/* **************************************************************************************************** */
int diagnostics(info *pointer, int type = 0, bool plot = false){
/* **************************************************************************************************** */
    /*
     This function is intended to
     */
    // Initialize Function Variables
    int exit = 0;
    int duplicate = 0;
    int order = 0;
    int recoil = 0;
    int alpha = 0;
    int coin = 0;
    int dt = 5;
    int time = 0;
    int num = 0;
    char str[500];
    float ebin = 0.3;
    TH1D *hist[3];
    TH2I *therm;
    ofstream outfile;
    // End
    
    sprintf(str, "%s/%s/%s_Diagnositic_%d.txt", pointer->path, pointer->name, pointer->name, type);
    outfile.open(str);
    if (outfile.fail()) {
        outfile.close();
        std::cout << "Failed to open file: " << str << "\n";
        std::cout << "Using default location...\n";
        sprintf(str, "/Users/ponce10/Root_Folder/Default_Path/%s_Channel_%d_Diagnositic_%d.txt", pointer->name, pointer->channel, type);
        outfile.open(str);
    }
    for (int i = 0; i < 3; i++){
        sprintf(str, "Histogram_%d", i);
        hist[i] = new TH1D(str, str, tBin, 0, tBin);
    }
    time = floor((pointer->time[pointer->last - 1] - pointer->time[0])/(3.6*pow(10,12)));
    num = time/12; // This calculates the maximum into 5 minute intervals
    time = int(num/5.0);
    therm = new TH2I("Thermal", "Thermal", num, 0, time, tBin, 0, tBin);
    
    if (type == 0){
        for (int i = 0; i < pointer->last; i++){
            if (pointer->high[i] == false) alpha++;
            if (pointer->veto[i] == false && pointer->high[i]) coin++;
            if (pointer->veto[i] && 8000 < pointer->bin[i]) recoil++;
            if (i != 0 && pointer->time[i] - pointer->time[i-1] == 0) duplicate++;
            if (i != 0 && pointer->time[i] - pointer->time[i-1] < 0) order++;
        }
        std::cout << "File Diagnostics: ";
        std::cout << "Status:\tValue\n";
        std::cout << "Total time (hrs)\t" << (pointer->time[pointer->last - 1] - pointer->time[0])/(3.6*pow(10, 12)) << "\n";
        std::cout << "Start (ns)\t" << pointer->time[0] << "\n";
        std::cout << "Finish (ns)\t" << pointer->time[pointer->last - 1] << "\n";
        std::cout << "Events\t" << pointer->last << "\n";
        std::cout << "Recoils\t" << recoil << "\n";
        std::cout << "Alphas\t" << alpha << "\n";
        std::cout << "LE Coin\t" << coin << "\n";
        
        outfile << "Status:\tValue\n";
        outfile << "Total time (hrs)\t" << (pointer->time[pointer->last - 1] - pointer->time[0])/(3.6*pow(10, 12)) << "\n";
        outfile << "Start (ns)\t" << pointer->time[0] << "\n";
        outfile << "Finish (ns)\t" << pointer->time[pointer->last - 1] << "\n";
        outfile << "Events\t" << pointer->last << "\n";
        outfile << "Recoils\t" << recoil << "\n";
        outfile << "Alphas\t" << alpha << "\n";
        outfile << "LE Coin\t" << coin << "\n";
        exit = 1;
    } // Marker diagnositic values
    else if(type == 1){
        for (int i = 0; i < 3; i++) hist[i]->Reset();
        for (int i = 0; i < pointer->last; i++) {
            hist[0]->Fill(pointer->bin[i]);
            if (pointer->id[i] == -1) hist[1]->Fill(pointer->bin[i]);
            else hist[2]->Fill(pointer->bin[i]);
            therm->Fill(pointer->time[i], pointer->bin[i]);
        }
        outfile << "000_Center\t001_Full\t002_Laser\t003_Signal\n";
        for (int i = 0; i < tBin; i++){
            outfile << hist[0]->GetXaxis()->GetBinCenter(i + 1) << "\t";
            outfile << hist[0]->GetBinContent(i + 1) << "\t";
            outfile << hist[1]->GetBinContent(i + 1) << "\t";
            outfile << hist[2]->GetBinContent(i + 1) << "\n";
        }
        outfile.close();
        sprintf(str, "%s/%s/%s_Channel_%d_Diagnositic_%d-2.txt", pointer->path, pointer->name, pointer->name, pointer->channel, type);
        outfile.open(str);
        if (outfile.fail()) {
            outfile.close();
            std::cout << "Failed to open file: " << pointer->filename << "\n";
            std::cout << "Using default location...\n";
            sprintf(str, "/Users/ponce10/Root_Folder/Default_Path/%s_Channel_%d_Diagnositic_%d-2.txt", pointer->name, pointer->channel, type);
            outfile.open(str);
        }
        outfile << "000_Time\t001_UCenter\t002_Uncalibrated\n";
        for (int i = 0; i < tBin; i++) for (int j = 0; j < num; j++){
            outfile << therm->GetXaxis()->GetBinCenter(j + 1) << "\t";
            outfile << therm->GetYaxis()->GetBinCenter(i + 1) << "\t";
            outfile << therm->GetBinContent(i + 1, j + 1) << "\n";
        }
        outfile.close();
    } // Uncalibrated file separation and heat map
    /*/ Remove later 20161021
    else if(type == 2){
        for(int i = 0; i < 3; i++) hist[i]->Delete();
        therm->Delete();
        for (int i = 0; i < 3; i++){
            sprintf(str, "Histogram_%d", i);
            hist[i] = new TH1D(str, str, int(108/ebin), 4, 112);
        }
        therm = new TH2I("Thermal", "Thermal", num, 0, time, int(108/ebin), 4, 112);
        for (int i = 0; i < 3; i++) hist[i]->Reset();
        
        for (int i = 0; i < pointer->last; i++) {
            hist[0]->Fill(pointer->energy[i]);
            if (pointer->id[i] == -1) hist[1]->Fill(pointer->energy[i]);
            else hist[2]->Fill(pointer->energy[i]);
            therm->Fill(pointer->time[i], pointer->energy[i]);
        }
        outfile << "000_Center\t001_Full\t002_Laser\t003_Signal\n";
        for (int i = 0; i < int(108/ebin); i++){
            outfile << hist[0]->GetXaxis()->GetBinCenter(i + 1) << "\t";
            outfile << hist[0]->GetBinContent(i + 1) << "\t";
            outfile << hist[1]->GetBinContent(i + 1) << "\t";
            outfile << hist[2]->GetBinContent(i + 1) << "\n";
        }
        outfile.close();
     
        sprintf(str, "%s/%s/%s_Channel_%d_Diagnositic_%d-2.txt", pointer->path, pointer->name, pointer->name, pointer->channel, type);
        outfile.open(str);
        if (outfile.fail()) {
            outfile.close();
            std::cout << "Failed to open file: " << pointer->filename << "\n";
            std::cout << "Using default location...\n";
            sprintf(str, "/Users/ponce10/Root_Folder/Default_Path/%s_Channel_%d_Diagnositic_%d-2.txt", pointer->name, pointer->channel, type);
            outfile.open(str);
        }
        outfile << "000_Time\t001_UCenter\t002_Uncalibrated\n";
        for (int i = 0; i < int(108/ebin); i++) for (int j = 0; j < num; j++){
            outfile << therm->GetXaxis()->GetBinCenter(j + 1) << "\t";
            outfile << therm->GetYaxis()->GetBinCenter(i + 1) << "\t";
            outfile << therm->GetBinContent(i + 1, j + 1) << "\n";
        }
        outfile.close();
    } // Calibrated file separation and heat map
    /*/
    if (plot){
        therm->Draw();
        globalhist[0] = therm;
    }
    else{
        for(int i = 0; i < 3; i++) hist[i]->Delete();
        therm->Delete();
    }
    return exit;
}

/* **************************************************************************************************** */
int HEvents(info* pointer[N_detectors], long long dt = 1){
/* **************************************************************************************************** */
    /*/ This function is intended to produce a figure/data file of the signal behavior before and after
        a HE event to get an understanding of the signal behavior.
            dt => Is the time before/after a high energy event in ms
     
    /*/
    // Initialize Function Variables:
    int exit = 0;
    int initial = 1;
    int check = 0;
    int final = 1;
    int rec = 0;
    int sub = 0;
    int counter = 0;
    TGraph *hesub[4];
    TGraph *herec[4];
    TCanvas *canvas[4];
    TH1D *forward;
    TH1D *backward;
    TH1D *diff = new TH1D("Diff", "Diff", 8192, 0, 100);
    ofstream frec, fsub;
    // End
    
    dt = dt*pow(10, 6);
    
    
    // This will set all vetoes to true, set the laser false on all channels and veto coincidence events
    // where at least one channel was over bin value 8000.
    for (int i = 0; i < 4; i++){
        xveto(pointer[i], pointer[i], -1);
        xveto(pointer[i], pointer[i], 0);
    }
    for (int i = 0; i < 3; i++) for (int j = i + 1; j < 4; j++) {
        xveto(pointer[i], pointer[j], 1);
    }
    
    /*/ At this point the high energy coincidence events would be vetoed as false. Next I need to
     get all time before and after the veto ignoring the laser. From here on we are looking at each 
     detector independently./*/
    for (int chan = 0; chan < N_detectors; chan++) {
        if (chan != 3) continue;
        
        // Setup the canvas and save files
        sprintf(pointer[chan]->filename, "%s%d", "Channel", pointer[chan]->channel);
        canvas[chan] = new TCanvas(pointer[chan]->filename,  pointer[chan]->filename, 1000, 1000);
        canvas[chan]->Divide(1, 2, 0, 0);
        rec = sub = 0;
        hesub[chan] = new TGraph(10000000);
        herec[chan] = new TGraph(10000000);
        forward = new TH1D("Forward", "Forward", 8192, 0, 100);
        backward = new TH1D("Backward", "Backward", 8192, 0, 100);
        
        
        sprintf(pointer[chan]->filename, "%s%s%s%d%s", pointer[chan]->path, pointer[chan]->name, "_HE_Rec_Channel", pointer[chan]->channel, ".txt");
        frec.open(pointer[chan]->filename);
        //frec << "Delta\tBin\n";
        frec << "Delta\tTime\tBin\tEnergy\n";
        sprintf(pointer[chan]->filename, "%s%s%s%d%s", pointer[chan]->path, pointer[chan]->name, "_HE_Sub_Channel", pointer[chan]->channel, ".txt");
        fsub.open(pointer[chan]->filename);
        //fsub << "Delta\tBin\n";
        fsub << "Delta\tTime\tBin\tEnergy\n";
        
        // Scan through current detector data
        pointer[chan]->current = 0;
        for (int i = 0; i < pointer[chan]->last; i++){
            if (pointer[chan]->tag[i] == 1) continue; // skip laser criterion or an otherwise usable point
            else if (pointer[chan]->veto[i] && pointer[chan]->bin[i] < 8000) continue;
            else if (pointer[chan]->veto[i] && 8000 < pointer[chan]->bin[i]){
                initial = 0;
                check = 1;
                final = 1;
                while (0 < i - check && pointer[chan]->time[i] <= pointer[chan]->time[i - check] + dt + 3*pow(10, 7)) {
                    if (pointer[chan]->tag[i - check] != 1 && 8000 < pointer[chan]->bin[i - check]) {
                        //(8000 < pointer[chan]->bin[i - check] || pointer[chan]->veto[i - check] == false)) {
                        initial = 0;
                        break;
                    } // Point is no good because there is a nearby high energy event.
                    if (pointer[chan]->tag[i - check] != 1 && pointer[chan]->time[i] <= pointer[chan]->time[i - check] + dt) {
                        initial = check;
                    }
                    check++;
                } // Lower bound
                while (i + final < pointer[chan]->last && pointer[chan]->time[i + final] - dt <= pointer[chan]->time[i]) {
                    if (pointer[chan]->tag[i + final] != 1 && 8000 < pointer[chan]->bin[i + final]){
                        //(8000 < pointer[chan]->bin[i + final] || pointer[chan]->veto[i + final] == false)){
                        final = 0;
                        break;
                    } // Point is no good because there is a nearby high energy event.
                    final++;
                } // Upper bound
                if (initial != 0 && final != 0) {
                    rec++;
                    for (int j = i - initial; j < i + final; j++) if (pointer[chan]->tag[j] != 1 && pointer[chan]->veto[j]){
                        frec << pointer[chan]->time[j] - pointer[chan]->time[i] << "\t" << pointer[chan]->time[j] << "\t";
                        //frec << pointer[chan]->bin[j] << "\t" << pointer[chan]->energy[j] << "\n";
                        //if (pointer[chan]->time[j] - pointer[chan]->time[i] < -3*pow(10, 7)*0 ) backward->Fill(pointer[chan]->energy[j]);
                        //else if (3*pow(10, 7)*0 < pointer[chan]->time[j] - pointer[chan]->time[i]) forward->Fill(pointer[chan]->energy[j]);
                        //frec << pointer[chan]->time[j] - pointer[chan]->time[i] << "\t" << pointer[chan]->bin[j] << "\n";
                        if (rec < 1000){
                            herec[chan]->SetPoint(herec[chan]->GetN(), (double)(pointer[chan]->time[j] - pointer[chan]->time[i]), (double)pointer[chan]->bin[j]);
                        }
                    }
                    counter++;
                }
            }
            else{/*/
                initial = 0;
                check = 1;
                final = 1;
                sub++;
                while (0 < i - check && pointer[chan]->time[i] <= pointer[chan]->time[i - check] + dt + 3*pow(10, 7)) {
                    if (pointer[chan]->tag[i - check] != 1 && 8000 < pointer[chan]->bin[i - check]) {
                        //(8000 < pointer[chan]->bin[i - check] || pointer[chan]->veto[i - check] == false)) {
                        initial = 0;
                        break;
                    }
                    if (pointer[chan]->tag[i - check] != 1 && pointer[chan]->time[i] <= pointer[chan]->time[i - check] + dt) {
                        initial = check;
                    }
                    check++;
                } // Lower bound
                while (i + final < pointer[chan]->last && pointer[chan]->time[i + final] - dt <= pointer[chan]->time[i]) {
                    if (pointer[chan]->tag[i + final] != 1 && 8000 < pointer[chan]->bin[i + final]){
                       //(8000 < pointer[chan]->bin[i + final] || pointer[chan]->veto[i + final] == false)){
                        final = 0;
                        break;
                    }
                    final++;
                } // Upper bound
                if (initial != 0 && final != 0) {
                    for (int j = i - initial; j < i + final; j++) if (pointer[chan]->tag[j] != 1 ){
                        fsub << pointer[chan]->time[j] - pointer[chan]->time[i] << "\t" << pointer[chan]->bin[j] << "\n";
                        if (sub <1000) {
                            hesub[chan]->SetPoint(hesub[chan]->GetN(), (double)(pointer[chan]->time[j] - pointer[chan]->time[i]), (double)pointer[chan]->bin[j]);
                        }
                    
                    }
                  /*/ continue;}
        }
        std::cout << "Number of Sub: " << sub << "\tNumber of Rec: " << rec << "\r";
        std::cout << "\nNumber of Points: " << counter <<  "\n";
        frec.close();
        fsub.close();
        canvas[chan]->cd(1);
        forward->Draw();
        backward->Draw("same");
        //hesub[chan]->Draw();
        canvas[chan]->cd(2);
        diff->Add(backward, forward, -1);
        diff->Draw();
        //herec[chan]->Draw();
        
        sprintf(pointer[chan]->filename, "%s%s%s%d%s", pointer[chan]->path, pointer[chan]->name, "_HE_Rec_Sub_Channel", pointer[chan]->channel, "_Hist.txt");
        frec.open(pointer[chan]->filename);
        frec << "Forward\tBackward\n";
        for (int i = 0; i < 8192; i++) frec << forward->GetBinContent(i+1) << "\t" << backward->GetBinContent(i+1) << "\n";
        frec.close();
        
    }

    return exit;
}

/* **************************************************************************************************** */
int ptcheck(info *point, info* check, int range = 10){
/* **************************************************************************************************** */
    /*/ This function is intended to check if the event at point->current is within range of an event
        in check->time.
        range   -> This is the range to look at in microseconds
        point   -> container for event to check
        check   -> container of events to check against
    /*/
    // Initialize Function Variables
    int exit = 0;
    int counter = 0;
    long long timeago = point->time[point->current];
    // End
    
    range = range*pow(10, 3); // Convert to ns
    
    /*/ The idea here is to get a bounding time:
                tau_j-1 <= t_i <= tau_j
        The first if-statement will skip all uncalibrated points and tagged events
        The first while-loop will make sure that tau_j-1 is less than t_i
        The second while-loop will make sure that tau_j is greater than t_i
        The issue is that a situation may arise where j = 0 or j_max is either a lower or upper
        bound this is dealt with by adding the additional criterion to each of the last two
        if-statements to check if the value they would check is beyond the data set bounds./*/
    while (1 <= check->current &&
           timeago < check->time[check->current - 1]) check->current--;
    while (check->current < check->last &&
           check->time[check->current] < timeago) check->current++;
    if (check->current < check->last  &&
        timeago - range <= check->time[check->current] &&
        check->time[check->current] <= timeago + range) {
        exit = 1;
    }
    else if (0 < check->current &&
             timeago - range <= check->time[check->current - 1] &&
             check->time[check->current - 1] <= timeago + range) {
        check->current = check->current - 1;
        exit = 1;
    }
    return exit;
}

/* **************************************************************************************************** */
int marker(info *pointer[N_detectors], int range = 10, bool markrecoil = true, bool trouble = false){
/* **************************************************************************************************** */
    /*/
        This function is will mark the ID column of the loaded file. 
     Inputs:
        data[N] =>  N-detector data
        range    =>  Coincidence gate in us
     
     Setup to perform coincidence and recoil marking, note that the coincidence marking 'should'
     be the same regardless of choosen detector. Coincidence file is generated recoding the bin
     and time for each of the four detectors.
     Define Recoil Time Condition (units are in ms):
                                xaxis = stj2 - stj0
                                yaxis = stj4 - stj6
                                    pval[0]         pval[1]
                            STJ 0: x(4.5, 5.5)      y(-3, 3)        =>      0
                            STJ 2: x(-5.5, -4.5)    y(-3, 3)        =>      1
                            STJ 4: x(-3, 3)         y(4.5, 5.5)     =>      2
                            STJ 6: x(-3, 3)         y(-5.5, -4.5)   =>      3
     
     This section identifies recoil events as coincidence events in which there is a 0 in one detector
     and non-zero values in the remaining detectors. The recoils are marked on the high column as true
     and events after the recoil (non-laser, non-coincidence) are saved to the file
     <file>_Recoil_10ms_Channel_#.txt
     Final step is to record the fully processed energy calibrated data as
     <file>_Energy_Calibrated_Channel_#.txt
     the columns are:    Time (ns)   Bin     Energy (eV)     Veto    High
     which will allow for different time cuts without reprocessing the data.
     
     Check that points are in coincidence and none of the points are laser events, mark
     veto column false to indicate a coincidence event.
     Note:
                                            7000 <= ##
                                            0 < # < 7000
     
     ID Values:
     _ _ _ _ _ _ _ _ _
     A B C D E F G H I
     A: Event type
     -1     =>  Laser, tagged from the beginning
     0      =>  Anticoincidence, all values
     i      =>  Coincidence, i detectors
     
     B - E: Detector Bin Value (Channel 0, 2, 4, 6)
     0      =>  0
     1      =>  0 < # < 7000
     2      =>  7000 < ##
     
     F - I: Detector Timing (Channel 0, 2, 4, 6)
     0      => Not Applicable
     i      => i-th timing
    /*/
    
    // Initialize Function Variables
    int exit = 1;
    int channel = 0; // Choice of stj is irrelevent
    int counter = 0;
    int idval = 0;
    int cdet = 0;
    int *(pval[2]);
    int valone[2] = {4500, 5500};
    int valtwo[2] = {-3000, 3000};
    int alphacutofflimit = 7000;
    int lasermark = 1;
    int tcount[13]; // Trouble Shooting
    int numdet = 0;
    int actdet[N_detectors];
    char str[500];
    short values[N_detectors];
    long long lineup[N_detectors];
    long long t20 = 0;
    long long t64 = 0;
    ofstream outfile;
    // End
    // Zero out storage arrays and determine number of active detectors.
    for (int i = 0; i < N_detectors; i++) {lineup[i] = 0; values[i] = 0; actdet[i] = 0;}
    for (int i = 0; i < N_detectors; i++) if (0 < pointer[i]->last) {numdet++; actdet[i] = 1;}
    //                      Set veto, recoil and high column default values
    for (int i = 0; i < N_detectors; i++) if (actdet[i]){
        xveto(pointer[i], pointer[i], -1);    // Set all veto, recoil and high to true (this means use all points)
        xveto(pointer[i], pointer[i], -3);    // Set all high to false
        xveto(pointer[i], pointer[i], -4);    // Set all recoil to false
        xveto(pointer[i], pointer[i], 3);     // Set all ID values to false (0)
        xveto(pointer[i], pointer[i], 0);     // Set all vetos for laser events false and ID to -1
        pointer[i]->current = 0;
    }
    
    // Look for close proximity points below range value in active detectors
    for (int i = 0; i < N_detectors; i++) if (actdet[i]){
        counter = 0;
        for (int j = 1; j < pointer[i]->last; j++) if (pointer[i]->time[j] - pointer[i]->time[j - 1] < range*pow(10, 3)) {
            counter++;
            if (trouble) printf("\nEntry Value\tTime\n%d\t%lld\n%d\t%lld\n", j - 1, pointer[i]->time[j - 1], j, pointer[i]->time[j]);
        }
        if (counter) printf("STJ %d has %d points within gate range of one another...\n", 2*i, counter);
    }
    if (trouble) printf("Finish searching for proximity events...\n");
    
    // Scan through all the data points and increase id value by 10**(channel_#/2)
    /*/
     This section will run through each available detector and identify all coincidence events, first line.
     At the start of each detector check all detector starting points are reset to zero, second line.
     Next the data in the i-th detector is checked, third line, with the tracking value for the data set to the
     counter j. Next, the remaining detectors i + 1 -> end are checked against the current detector note that
     you do not need to check BA once AB is checked because each detector has a dedicated location to identify
     a coincidence event. Laser events are ignored for coincidence, this is being adjusted...
    /*/// Mark Coincidence events for each detector
    for (int i = 0; i < N_detectors; i++) if (actdet[i]){
        for (int k = 0; k < N_detectors; k++) pointer[k]->current = 0; // Reset pointer current to 0 for each detector
        for (int j = 0; j < pointer[i]->last; j++){
            pointer[i]->current = j;
            for (int k = i + 1; k < N_detectors; k++) if (actdet[k]){
                if (!pointer[i]->tag[pointer[i]->current] && ptcheck(pointer[i], pointer[k]) && !pointer[k]->tag[pointer[k]->current]) {
                    if (!pointer[i]->id[pointer[i]->current]) pointer[i]->id[pointer[i]->current]+= pow(10, i);
                    if (!pointer[k]->id[pointer[k]->current]) pointer[k]->id[pointer[k]->current]+= pow(10, k);
                    pointer[i]->id[pointer[i]->current]+= pow(10, k);
                    pointer[k]->id[pointer[k]->current]+= pow(10, i);
                } // Two events are in coincidence and neither of them is a laser event
                if(ptcheck(pointer[i], pointer[k]) && (pointer[i]->tag[pointer[i]->current] || pointer[k]->tag[pointer[k]->current])) {
                    if (!pointer[i]->tag[pointer[i]->current]) for (int m = i + 1; m < k; m++) {
                        if (ptcheck(pointer[i], pointer[m])) {
                            pointer[m]->id[pointer[m]->current] = -1;
                            pointer[m]->tag[pointer[m]->current] = 1;
                        }
                    }
                    pointer[i]->id[pointer[i]->current] = -1;
                    pointer[k]->id[pointer[k]->current] = -1;
                    pointer[i]->tag[pointer[i]->current] = 1;
                    pointer[k]->tag[pointer[k]->current] = 1;
                } // Two events are in coincidence and at least one of them is a laser event mark them as laser instead of coincidence
            }
        }
    }
    if(trouble) printf("Finish marking coincidence events...\n");
    
    
    /*/ Trouble Shooting
    std::cout << "4-Detector Coincidence\n";
    TH1D *thist[4];
    for (int i = 0; i < N_detectors; i++)if (actdet[i]) {
        counter = 0;
        sprintf(str, "%s%d", "STJ ", 2*i);
        thist[i] = new TH1D(str, str, 3000, 0, 1500);
        for (int j = 0; j < pointer[i]->last; j++){
            thist[i]->Fill(pointer[i]->id[j]);
            if (pointer[i]->id[j] == 1111) counter++;
        }
        std::cout << "Total coincidence events in Detector " << 2*i << ": " << counter << "\n";
        std::cout << "Channel Max Value: " << findmax(pointer[i]->id, pointer[i]->last) << "\n";
    }
    TCanvas *canvas = new TCanvas("Trouble", "Trouble", 1000, 1000);
    TGraph *plot = new TGraph(10000000);
    TGraph *plot2 = new TGraph(10000000);
    
    //*/// Trouble Shooting
    
    // Create file of 4-detector coincidence
    sprintf(str, "%s%s_Coinicidence.txt", pointer[0]->path, pointer[0]->name);
    outfile.open(str);
    if(outfile.fail()) printf("Failed to open save location: %s\n", str);
    for (int i = 0; i < N_detectors; i++) if(actdet[i]) {
        outfile << pointer[i]->channel << "_Time (ns)\t";
        outfile << pointer[i]->channel << "_Bin" << (i == N_detectors - 1 ? "\n": "\t");
        
    }
    //outfile << "0_Time(ns)\t0_Bin\t2_Time(ns)\t2_Bin\t4_Time(ns)\t4_Bin\t6_Time(ns)\t6_Bin\n";
    
    // This section will mark the recoil column for each individual detector
    if (trouble) for (int i = 0; i < 13; i++) tcount[i] = 0;
    if (markrecoil) for (int i = 0; i < N_detectors; i++) {
        if (trouble) for (int i = 0; i < 5; i++) tcount[i] = 0;
        for (int j = 0; j < pointer[i]->last; j++){
            if (pointer[i]->id[j] < 0) continue; // Skip laser events
            /*/else if (pointer[i]->id[j] == 0 && pointer[i]->bin[j] == 0) {
                pointer[i]->recoil[j] = true;
            } // Anticoincidence recoil criterion/*/
            else if (pointer[i]->id[j] == 1111){
                pointer[i]->current = j;
                for (int k = 0; k < N_detectors; k++) if (k != i) ptcheck(pointer[i], pointer[k]);
                for (int k = 0;  k < N_detectors; k++) lineup[k] = pointer[k]->time[pointer[k]->current];
                if (i == 0) for (int m = 0; m < N_detectors; m++){
                    outfile << pointer[m]->time[pointer[m]->current] << "\t";
                    outfile << pointer[m]->bin[pointer[m]->current] << (m == N_detectors - 1 ? "\n" : "\t");
                }
                if (findmax(lineup, N_detectors) - findmin(lineup, N_detectors) < range*pow(10, 3)) {
                    for (int k = 0; k < N_detectors; k++) {
                        // Tag recoils and some alphas from (0, #, #, #), (0, 0, #, #), (0, 0, 0, #), and (0, 0, 0, 0)
                        if ((!pointer[k]->bin[pointer[k]->current] &&
                             0 < pointer[(k + 1)%4]->bin[pointer[(k + 1)%4]->current] &&
                             pointer[(k + 1)%4]->bin[pointer[(k + 1)%4]->current] < alphacutofflimit &&
                             0 < pointer[(k + 2)%4]->bin[pointer[(k + 2)%4]->current] &&
                             pointer[(k + 2)%4]->bin[pointer[(k + 2)%4]->current] < alphacutofflimit &&
                             0 < pointer[(k + 3)%4]->bin[pointer[(k + 3)%4]->current] &&
                             pointer[(k + 3)%4]->bin[pointer[(k + 3)%4]->current] < alphacutofflimit) ||
                            (!pointer[k]->bin[pointer[k]->current] &&
                             !pointer[(k + 1)%4]->bin[pointer[(k + 1)%4]->current] &&
                             0 < pointer[(k + 2)%4]->bin[pointer[(k + 2)%4]->current] &&
                             pointer[(k + 2)%4]->bin[pointer[(k + 2)%4]->current] < alphacutofflimit &&
                             0 < pointer[(k + 3)%4]->bin[pointer[(k + 3)%4]->current] &&
                             pointer[(k + 3)%4]->bin[pointer[(k + 3)%4]->current] < alphacutofflimit) ||
                            (!pointer[k]->bin[pointer[k]->current] &&
                             !pointer[(k + 1)%4]->bin[pointer[(k + 1)%4]->current] &&
                             !pointer[(k + 2)%4]->bin[pointer[(k + 2)%4]->current] &&
                             0 < pointer[(k + 3)%4]->bin[pointer[(k + 3)%4]->current] &&
                             pointer[(k + 3)%4]->bin[pointer[(k + 3)%4]->current] < alphacutofflimit) ||
                            (!pointer[k]->bin[pointer[k]->current] &&
                             !pointer[(k + 1)%4]->bin[pointer[(k + 1)%4]->current] &&
                             !pointer[(k + 2)%4]->bin[pointer[(k + 2)%4]->current] &&
                             !pointer[(k + 3)%4]->bin[pointer[(k + 3)%4]->current])
                            ){
                            if (trouble){
                                tcount[0]++; // total number of events satisfying initial recoil criterion
                                if (!pointer[(k + 3)%4]->bin[pointer[(k + 3)%4]->current]) tcount[4]++; // total number of events satisfying initial recoil criterion (0, 0, 0, 0)
                                else if (!pointer[(k + 2)%4]->bin[pointer[(k + 2)%4]->current]) tcount[3]++; // total number of events satisfying initial recoil criterion (0, 0, 0, #)
                                else if (!pointer[(k + 1)%4]->bin[pointer[(k + 1)%4]->current]) tcount[2]++; // total number of events satisfying initial recoil criterion (0, 0, #, #)
                                else if (!pointer[(k + 0)%4]->bin[pointer[(k + 0)%4]->current]) tcount[1]++; // total number of events satisfying initial recoil criterion (0, #, #, #)
                            }
                            if (k == 0 || k == 1) {pval[0] = valone; pval[1] = valtwo;} // STJ 0 and 2 x{4.6, 5.5} and y{-3, 3};
                            else {pval[0] = valtwo; pval[1] = valone;} // STJ 4 and 6 x{-3, 3} and y{4.6, 5.5}
                            t20 = pointer[1]->time[pointer[1]->current] - pointer[0]->time[pointer[0]->current];
                            t64 = pointer[3]->time[pointer[3]->current] - pointer[2]->time[pointer[2]->current];
                            
                            //if (i == 0) plot->SetPoint(plot->GetN(), t20,  t64); // Trouble Shooting
                            
                            cdet = where(lineup, findmin(lineup, N_detectors), N_detectors);
                            if ((cdet == 0 && (pval[0][0] <=  t20 &&  t20 <= pval[0][1] && pval[1][0] <=  t64 &&  t64 <= pval[1][1])) ||
                                (cdet == 1 && (pval[0][0] <= -t20 && -t20 <= pval[0][1] && pval[1][0] <=  t64 &&  t64 <= pval[1][1])) ||
                                (cdet == 2 && (pval[0][0] <=  t20 &&  t20 <= pval[0][1] && pval[1][0] <=  t64 &&  t64 <= pval[1][1])) ||
                                (cdet == 3 && (pval[0][0] <=  t20 &&  t20 <= pval[0][1] && pval[1][0] <= -t64 && -t64 <= pval[1][1]))){
                                for (int m = 0; m < 4; m++) pointer[m]->recoil[pointer[m]->current] = true;
                                if (trouble) {
                                    if (!pointer[(k + 3)%4]->bin[pointer[(k + 3)%4]->current]) tcount[8]++; // total number of events satisfying initial recoil criterion (0, 0, 0, 0)
                                    else if (!pointer[(k + 2)%4]->bin[pointer[(k + 2)%4]->current]) tcount[7]++; // total number of events satisfying initial recoil criterion (0, 0, 0, #)
                                    else if (!pointer[(k + 1)%4]->bin[pointer[(k + 1)%4]->current]) tcount[6]++; // total number of events satisfying initial recoil criterion (0, 0, #, #)
                                    else if (!pointer[(k + 0)%4]->bin[pointer[(k + 0)%4]->current]) tcount[5]++; // total number of events satisfying initial recoil criterion (0, #, #, #)
                                    tcount[cdet + 9]++;
                                }
                                // Calculates the number of tagged recoils
                                //plot2->SetPoint(plot2->GetN(), t20, t64); // Trouble Shooting
                                } // Recoil condition
                            //else if (i == 0 && k == 1) std::cout << pval[0][0] << " <= " <<  -t20 << " <= " << pval[0][1]  << "\t" << pval[1][0] << " <= "  << t64 << " <= " << pval[1][1] << "\n"; // Trouble
                            }
                        
                    }
                }
            } // Four detector coincidence recoil criterion
            /*/else {
                pointer[i]->current = j;
                for (int k = 0;  k < N_detectors; k++) values[k] = lineup[k] = 0;
                idkey(values, pointer[i]->id[pointer[i]->current], N_detectors);
                
                
                for (int k = 0; k < N_detectors; k++) if (k != i && values[k]) ptcheck(pointer[i], pointer[k]);
                if ((!pointer[i]->bin[pointer[i]->current] || alphacutofflimit <= pointer[i]->bin[pointer[i]->current])){
                    for (int k = 0; k < N_detectors; k++) if (k != i){
                        if (!values[k]) continue;
                        else if (values[k] && 0 < pointer[k]->bin[pointer[k]->current] &&
                                 pointer[k]->bin[pointer[k]->current] < alphacutofflimit){
                            pointer[i]->recoil[pointer[i]->current] = true;
                        }
                        else {
                            pointer[i]->recoil[pointer[i]->current] = false;
                            break;
                        }
                    }
                }
            } // Two or three detector coincidence recoil criterion/*/
        }
        if (trouble){
            std::cout << "\nNumber of Candidates: " << tcount[0] << "\n";
            std::cout << "Number of Possible Recoils (0, #, #, #): " << tcount[1]<< "\n";
            std::cout << "Number of Possible Recoils (0, 0, #, #): " << tcount[2] << "\n";
            std::cout << "Number of Possible Recoils (0, 0, 0, #): : " << tcount[3] << "\n";
            std::cout << "Number of Possible Recoils (0, 0, 0, 0): : " << tcount[4] << "\n";
            std::cout << "Number of Accepted  Recoils (0, #, #, #): " << tcount[5]<< "\n";
            std::cout << "Number of Accepted  Recoils (0, 0, #, #): " << tcount[6] << "\n";
            std::cout << "Number of Accepted  Recoils (0, 0, 0, #): : " << tcount[7] << "\n";
            std::cout << "Number of Accepted  Recoils (0, 0, 0, 0): : " << tcount[8] << "\n";
            std::cout << "Accepted Zeros STJ 0: " << tcount[9] << "\n";
            std::cout << "Accepted Zeros STJ 2: " << tcount[10] << "\n";
            std::cout << "Accepted Zeros STJ 4: " << tcount[11] << "\n";
            std::cout << "Accepted Zeros STJ 6: " << tcount[12] << "\n";
        }
        for(int i = 0; i < 13; i++) tcount[i] = 0;
        outfile.close();
    }
    if (trouble) printf("Finished marking recoils and untagging laser events...\n");

    for (int i = 0; i < N_detectors; i++) if (actdet[i]){
        for (int k = 0; k < N_detectors; k++) pointer[k]->current = 0;
        for (int j = 0; j < pointer[i]->last; j++) if (0 < pointer[i]->id[j]) {
            pointer[i]->current = j;
            idval = pow(10, 8);
            for (int k = 0; k < N_detectors; k++) if (actdet[k]){
                values[k] = pointer[i]->id[pointer[i]->current]/pow(10, 3 - k);
                pointer[i]->id[pointer[i]->current] -= values[k]*pow(10, 3 - k);
                if (k != i && values[k]) {
                    ptcheck(pointer[i], pointer[k], range);
                    if (pointer[k]->id[pointer[k]->current] != -1) idval += pow(10, 8);
                    else values[k] = 0;
                }
            }
            
            for (int k = 0; k < N_detectors; k++){
                if (values[k]) lineup[k] = pointer[k]->time[pointer[k]->current];
                else lineup[k] = 0;
            }
            sort_list(lineup, N_detectors);
            counter = 0;
            for (int k = 0; k < N_detectors; k++) if (!lineup[k]) counter++;
            for (int k = 0; k < N_detectors; k++) if (actdet[k])for (int m = 0; m < N_detectors; m++)if (actdet[m]){
                if (values[k] && pointer[k]->time[pointer[k]->current] == lineup[m]) {
                    idval += (m + 1 - counter)*pow(10, 3 - k);
                    // Remove laser events for primary event from calibration
                    lasermark = 0;
                    if (m + 1 - counter == 1) {
                        while (pointer[k]->current + lasermark < pointer[k]->last && pointer[k]->time[pointer[k]->current + lasermark] - pointer[k]->time[pointer[k]->current] <= 20*pow(10, 6)) {
                            if (pointer[k]->tag[pointer[k]->current + lasermark]) pointer[k]->tag[pointer[k]->current + lasermark] = false;
                            lasermark++;
                        }
                    }
                }
            } // Enter the time category for each detector
            
            for (int k = 0; k < N_detectors; k++) if (actdet[k]){
                if (values[k]) {
                    if (!pointer[k]->bin[pointer[k]->current]) idval += 0*pow(10, 7 - k);
                    else if (7000 <= pointer[k]->bin[pointer[k]->current]) idval += 2*pow(10, 7 - k);
                    else idval += 1*pow(10, 7 - k);
                }
                else idval += 9*pow(10, 7 - k);
            } // Enter the bin category for each detector
            for (int k = 0; k < N_detectors; k++) if (actdet[k]) if (values[k]) pointer[k]->id[pointer[k]->current] = idval;
        /*/
            if (pointer[0]->time[pointer[0]->current] == 539992905000){
                std::cout << "This is the trouble point Finish\n";
                std::cout << "Time: " << pointer[0]->time[pointer[0]->current] << "\n";
                std::cout << "Bin: " << pointer[0]->bin[pointer[0]->current] << "\n";
                std::cout << "Tag: " << pointer[0]->tag[pointer[0]->current] << "\n";
                std::cout << "Veto: " << pointer[0]->veto[pointer[0]->current] << "\n";
                std::cout << "High: " << pointer[0]->high[pointer[0]->current] << "\n";
                std::cout << "id: " << pointer[0]->id[pointer[0]->current] << "\n";
            }/*/
        }
        
    } // Enter ID information and remove laser events from calibration...
    if (trouble) printf("Recorded ID information..\n");
    
    //for (int i = 0; i < N_detectors; i++) thist[i]->Draw("same"); // Trouble Shooting
    //plot->Draw(); // Trouble Shooting
    //plot2->Draw("Same"); // Trouble Shooting
    return exit;
}

/* **************************************************************************************************** */
double alpha_ij(float old_i0, float old_i1, float new_j0, float new_j1){
    /* **************************************************************************************************** */
    /*
     This function is intended to compute the matrix element alpha_ij for the conversion from one
     histogram -> old to a new histogram -> new. Note that the resulting values for the new histogram
     will in all likelihood be no integar values as a result of this conversion...
     */
    float numerator = (old_i1 < new_j1 ? old_i1 : new_j1) - (old_i0 > new_j0 ? old_i0 : new_j0);
    numerator = (numerator > 0 ? numerator : 0);
    return numerator/(old_i1 - old_i0);
    
}

/* **************************************************************************************************** */
int revise_energy(info *pointer, int tval = tBin, int numbin = tBin, double lower = 0,
                     double upper = tBin, bool trouble = false){
/* **************************************************************************************************** */
    /*
     This function is intended to convert the intensities (counts) from one histogram -> old to another
     histogram -> new. Note that the resulting intensity values for the new histogram will in all
     likelihood be non-integar values as a result of this conversion... The old histogram must be
     stored in an actual TH1D histogram under pointer->histogram[pointer->cslice] while the new 
     histogram will overwrite/replace the entries in the array pointer->hist[pointer->cslice][...]
     
     pointer    ->  Data, cslice should point to the histogram for conversion
     numbin     ->  Number of bins in the new histogram
     tval       ->  Number of bins in the old histogram
     lower      ->  Lower bound, inclusive
     upper      ->  Upper bound, exclusive
     
     OUTPUT
     To pointer->hist[pointer->current]
     
     Notes:
     This function does not change pointer->current
     This function updates pointer->xaxis to the new histograms center bin values...
     
     
     Tested and working 20160910
     */
    // Initialize Function Variables
    alpha_matrix *alpha = new alpha_matrix;
    int current = pointer->cslice;
    int cnew = 0;   // Current location on new range
    int cold = 0;   // Current location on old range
    int lold = 0;   // Last location on old range with overlap with new range
    int tcount = 0; //trouble
    // End
    
    // If calibration is marked as skip by m = b = 0 then the returned hist array is empty
    if (pointer->calib[current][1] == 0 && pointer->calib[current][3] == 0){
        for(int i = 0; i < tBin; i++) pointer->hist[current][i] = 0;
    }
    else{
        if (trouble){
            printf("Current partition: %d\n", current);
            printf("slope: %.2f\toffset: %.2f\n", pointer->calib[current][1], pointer->calib[current][3]);
        }
        // Generate histogram bounds
        for (int i = 0; i < tval + 1; i++){
            alpha->old[i] = pointer->histogram[current]->GetBinLowEdge(i + 1);
            if (fitmodel) alpha->old[i] = pointer->calib[current][1]*alpha->old[i] + pointer->calib[current][3];                        //                              Linear Calibration
            else alpha->old[i] = pointer->calib[current][1]*alpha->old[i] + pointer->calib[current][3] +  pointer->calib[current][5]*alpha->old[i]*alpha->old[i] ;    //                              Quadratic Calibration
        }
        for (int i = 0; i < numbin + 1; i++) alpha->newt[i] = lower + i*(upper - lower)*1.0/numbin;
        while(cold + 1 < tval && alpha->old[cold + 1] <= alpha->newt[cnew]) cold++;
        while(lold + 1 < tval && alpha->old[lold] < alpha->newt[numbin]) lold++;
        if(trouble){
            printf("Old x-axis edges: \n");
            for(int i = 0; i < (true ? tval + 1 : 10); i++) printf("%.2f\t", alpha->old[i]);
            printf("\nNew x-axis edges: \n");
            for(int i = 0; i < (true ? numbin + 1 : 10); i++) printf("%.2f\t", alpha->newt[i]);
            printf("\ncold: %d\tlold: %d\n", cold, lold);
            printf("Starting points old -> %.2f\t\tnew -> %.2f\n", alpha->old[cold], alpha->newt[cnew]);
            printf("End points old -> %.2f\t\tnew -> %.2f\n", alpha->old[lold], alpha->newt[numbin]);
        }
        for (int j = 0; j < numbin; j++) for (int i = cold; i < lold; i++){
            alpha->matrix[j][i - cold] = alpha_ij(alpha->old[i], alpha->old[i + 1], alpha->newt[j], alpha->newt[j + 1]);
        } // Calculate transformation matrix, alpha_ij
        if (trouble){
            tcount = 0;
            printf("\n\nAlpha Matrix:\n");
            for (int j = 0; j < numbin; j++) for (int i = cold; i < lold; i++){
                printf("%.2f\t", alpha->matrix[j][i - cold]);
                if(i == lold - 1){
                    if (tcount + cold < lold) printf("\t\t%.2f\n", pointer->histogram[current]->GetBinContent(tcount + cold + 1));
                    else printf("\n");
                    tcount++;
                }
            }
        }
        for(int i = 0; i < tBin; i++) pointer->hist[current][i] = 0; // Clear contents in hist array
        if(trouble) tcount = 0;
        for (int j = 0; j < numbin; j++) for(int i = cold; i < lold; i++) {
            pointer->hist[current][j] += alpha->matrix[j][i - cold]*pointer->histogram[current]->GetBinContent(i + 1);
            if (trouble) {
                printf("%.2f*%.2f", alpha->matrix[j][i - cold], pointer->histogram[current]->GetBinContent(i + 1));
                (i == lold - 1 ? printf(" = %.2f\n", pointer->hist[current][j]) : printf(" + "));
            }
        } // Populate hist array using values in histogram and transformation matrix
        if (trouble){
            printf("Old x-axis centers: \n");
            for(int i = cold; i < lold; i++) std::cout << (alpha->old[i + 1] + alpha->old[i])/2.0 << "\t";
            printf("\nOld y-axis: \n");
            for(int i = cold; i < lold; i++) std::cout << pointer->histogram[current]->GetBinContent(i + 1) << "\t";
            printf("\nNew x-axis centers: \n");
            for(int i = 0; i < numbin; i++) std::cout << (alpha->newt[i + 1] + alpha->newt[i])/2.0 << "\t";
            printf("\nNew y-axis: \n");
            for(int i = 0; i < numbin; i++) std::cout << pointer->hist[current][i] << "\t";
        }
        for (int i = 0; i < tval; i++) pointer->xaxis[i] = 0;
        for (int i = 0; i < numbin; i++) pointer->xaxis[i] = (alpha->newt[i] + alpha->newt[i + 1])/2.0; // Changes x-axis to new axis
    }
    delete alpha;
    return 1;
}

/* **************************************************************************************************** */
void residual(info *pointer, bool domain = false){
/* **************************************************************************************************** */
    /*/
        This function calculates the residuals for the current slice directed by pointer->cslice
     
        domain    => bool, toogle whether the fitted centroid unit is in eV or channel number,
                     default is false which refers to channel domain.
     
     20161208
     Residuals:         Fitted - Theory
     changed to =>      Theory - Fitted
     
    /*/
    // Initialize Function Variables
    char str[500];
    char word[500];
    double energy = 0;
    double fit = 0;
    double errfit = 0;
    double slope = pointer->calib[pointer->cslice][1];
    double errsl = pointer->calib[pointer->cslice][2];
    double offset = pointer->calib[pointer->cslice][3];
    double erroff = pointer->calib[pointer->cslice][4];
    double quad = pointer->calib[pointer->cslice][5];
    double errquad = pointer->calib[pointer->cslice][6];
    
    ofstream outfile;
    // End
    snprintf(word, sizeof(word), "");
    if(pointer->indexvalue < 0){
        snprintf(word, sizeof(word), "_%d-%d", -pointer->indexvalue, -pointer->indexvalue + 21);
        pointer->indexvalue = 2;
    }
    if (pointer->indexvalue < 10) snprintf(str, sizeof(str), "%s%s_Ch_%d_0%d_Residuals%s.txt", pointer->path, pointer->name, pointer->channel, pointer->indexvalue, word);
    else snprintf(str, sizeof(str), "%s%s_Ch_%d_%d_Residuals%s.txt", pointer->path, pointer->name, pointer->channel, pointer->indexvalue, word);
    
    if (!pointer->cslice) outfile.open(str);
    else outfile.open(str, ios::app);
    if(outfile.fail()){
        outfile.close();
        std::cout << "Failed to open file: " << str << "\n";
        std::cout << "Using default location...\n";
        snprintf(str, sizeof(str), "");
        sprintf(str, "/Users/ponce10/Root_Folder/Default_Path/%s_Ch_%d_02_Residuals.txt", pointer->name, pointer->channel);
        outfile.open(str);
    }
    if(!outfile.fail()){
        if (!pointer->cslice) outfile << "000_Slice\t001_Energy\t002_Residual\t003_Uncertainty\n";
        for (int i = 0; i < NUM_PEAKS; i++){
            if(pointer->param[NUM_PEAKS*pointer->cslice + i][0] == -1) continue;
            energy = pointer->param[NUM_PEAKS*pointer->cslice + i][1];
            fit = pointer->param[NUM_PEAKS*pointer->cslice + i][2];
            errfit = pointer->param[NUM_PEAKS*pointer->cslice + i][3];
            outfile << pointer->cslice << "\t" << energy << "\t";
            if (fitmodel) outfile << energy - (domain ? fit : slope*fit + offset) << "\t"; // linear
            else outfile << energy - (domain ? fit : slope*fit + quad*fit*fit) << "\t"; // quadratic without offset
            outfile << (domain ? errfit : pow(fit*fit*errsl*errsl + erroff*erroff + errfit*errfit*slope*slope, 0.5));
            outfile << "\n";
        }
    }
    outfile.close();
}

/* **************************************************************************************************** */
int driftcor(info *pointer, int start = 0, int end = 0, int standnum = 0, int val4 = 0, bool trouble = false){
/* **************************************************************************************************** */
    /*
     Unlike its predecessor calibration(...) this function will only correct for the drift of the
     detector and keep the x-axis in the original uncalibrated (but shifted) energy axis or bin 
     axis.
     pointer    =>  This is the detector whose gain drifted and needs correction,
                    the data must already be processed and contain the energy column in eV and peaks
     
     trouble    =>  As always this will be used to trouble shoot this function
     
     The idea on speak and sperr is that this is the standard for converting all partitions into. The idea
     here is to do a linear fit between the standard partition and the other partitions with the slope
     being the gain difference and the offset the difference between partitions for the laser specific
     offset term which appears to differ systematically from beginning to end of run. Should not be the 
     case but looks like it is the case.
     
     NOTE: Since both the gain and laser offset term are saved at this stage the user must hard code the
     removal of the offset term when handling the signal half so that only the gain is taken.
    */
    // Initialize Function Variables
    int spt = 0;
    int tpt = 0;
    int gpt = 0;
    int numarray[NUM_COLUMN];
    float enarray[NUM_COLUMN];
    double speak[NUM_PEAKS];
    double sperr[NUM_PEAKS];
    double sener[NUM_PEAKS];
    double tpeak[NUM_PEAKS];
    double tperr[NUM_PEAKS];
    double tener[NUM_PEAKS];
    double gain;
    double offset;
    double *out;
    const double *outerr;
    TGraphErrors *linefit;
    TF1 *fun;
    // End
    
    val4 = 0; // Place holders for possible entries...
    
    for (int i = 0; i < NUM_COLUMN; i++) {
        numarray[i] = enarray[i] = 0;
        speak[i] = sperr[i] = sener[i] = 0;
        tpeak[i] = tperr[i] = tener[i] = 0;
    }// Zero out arrays
    for (int i = 0; i < NUM_COLUMN; i++){
        for (int j = 0; j < NUM_PEAKS; j++) {
            if(j == 0) enarray[i] = pointer->param[NUM_PEAKS*i][1];
            if(pointer->param[NUM_PEAKS*i + j][0] != -1) numarray[i]++;
            else break;
        }
    }// Count the number of peaks in each partition and the starting energy
    if (false && trouble){
        printf("Partition\t Energy\tNumber\n");
        for(int i = 0; i < pointer->endslice; i++) printf("%d\t%0.2f\t%d\n", i, enarray[i], numarray[i]);
    }// Print out list of values from above
    
    //Temporary need to make this more advanced for the data without long range intensities
    spt = int(enarray[standnum]/3.5);
    for (int i = spt + 1; i < spt + numarray[standnum] - 1; i++){
        sener[i] = 1;
        speak[i] = pointer->param[standnum*NUM_PEAKS + i - spt][2];
        sperr[i] = pointer->param[standnum*NUM_PEAKS + i - spt][3];
    }// Skip the first and last peaks...
    if (trouble){
        //printf("Block read from param\n");
        //for (int i = 0; i < 2*NUM_PEAKS; i++) printf("%d\t%0.2f\t%0.2f\n", i, pointer->param[i][0], pointer->param[i][1]);
        printf("Using this as the standard\n");
        printf("Peak\tdError\tEnergy\tUse\n");
        for(int i = 0; i < NUM_PEAKS; i++) if(speak[i] != 0) printf("%0.4f\t%0.4f\t%0.4f\t%0.4f\n", speak[i], sperr[i], 3.5*i, sener[i]);
    }// Print out list of values from above
    pointer->cslice = 0;
    
    for (int i = 0; i < pointer->endslice; i++){
        for (int j = 0; j < NUM_PEAKS; j++) tener[j] = tpeak[j] = tperr[j] = 0;
        tpt =int(enarray[i]/3.5);
        for (int j = tpt + 1; j < tpt + numarray[i] - 1; j++){
            tener[j] = 1;
            tpeak[j] = pointer->param[NUM_PEAKS*i + j - tpt][2];
            tperr[j] = pointer->param[NUM_PEAKS*i + j - tpt][3];
        }// Skip the first and last peaks...
        if (i == 0 && trouble){
            printf("On partition number: %d\n", i);
            printf("Using this as the test\n");
            printf("Delta\tdDelta\tEnergy\tUse\n");
            for(int j = 0; j < NUM_PEAKS; j++) if(tpeak[j] != 0) printf("%0.4f\t%0.4f\t%0.4f\t%0.4f\n", tpeak[j], tperr[j], 3.5*j, tener[j]);
        }// Print out list of values from above
        /*/
         The index values on j will exclude the endpoints on the standard to which the gain matching is performed
         The test value endpoints will likewise be excluded...
        /*/
        tpt = gpt = 0;
        dot(tener, sener, NUM_PEAKS); // Will return tener with 1's in overlapping energy positions.
        fun = new TF1("fun2", "pol1");
        linefit = new TGraphErrors(NUM_PEAKS);
        for (int j = 0; j < NUM_PEAKS; j++) if (tener[j]){
            if (tpt == 0) tpt = j;
            if (i == standnum && trouble) printf("%0.4f\t%0.4f\n", tpeak[j], speak[j]);
            linefit->SetPoint(linefit->GetN(), tpeak[j], speak[j]);
            linefit->SetPointError(linefit->GetN() - 1, tperr[j], sperr[j]);
            gpt++;
        }
        if (gpt <= 0) {
            printf("ERROR ERROR ERROR\nThere was no overlap of the two arrays no gain matching possible\n");
            pointer->calib[i][0] = i;
            pointer->calib[i][1] = 0;
            pointer->calib[i][2] = 0;
            pointer->calib[i][3] = 0;
            pointer->calib[i][4] = 0;
        }
        else {
            if (0 < i && pointer->calib[i - 1][1] != 0) {
                gain = pointer->calib[i - 1][1];
                offset = pointer->calib[i - 1][3];
            }
            else{
                gain = (speak[tpt + 1] - speak[tpt])/(tpeak[tpt + 1] - tpeak[tpt]);
                offset = speak[tpt] - gain*tpeak[tpt];
            }
            fun->SetParameters(offset, gain);
            linefit->Fit(fun, "Q");
            out = fun->GetParameters();
            outerr = fun->GetParErrors();
            pointer->calib[i][0] = i;
            pointer->calib[i][1] = out[1];
            pointer->calib[i][2] = outerr[1];
            pointer->calib[i][3] = out[0];
            pointer->calib[i][4] = outerr[0];
        }
        if (i == 0) printf("Number of bins in original, new, and end point is %d\n", pointer->numbin);
        revise_energy(pointer, pointer->numbin, pointer->numbin, start, end, false);
        pointer->cslice++;
    }// This is where the gain correction values will be generated...
    printf("Drift Correction Calculation Completed...\n");
    return 1;
}

/* **************************************************************************************************** */
int process(int type = 0,  float sig = 0, char file[500] = tfile, double randomvalue = 0, int val1 = 0){
/* **************************************************************************************************** */
    /* This function is intended to perform the start to finish main function of processing either 
       CAEN or XIA data.
    */
    if(type == 1){
        /*/                                     NOTES (20171011)
         This piece is intended to play with the sigma and amplitude of the peak search function
         to determine good values.
         /*/
        // Function Variables:
        int counter = 0;
        int attempt = 0;
        int exit = 0;
        int num = 0;
        int rebin = val1;
        char str[500];
        char *path = tpath;
        char *ext = text;
        bool trouble = true;
        double binstart = 0;
        double binend = 0;
        double poff_sig = 1;
        double poff_amp = 0.05;
        
        info *data = new info;
        TSpectrum *peak_off = new TSpectrum(NUM_PEAKS);
        TH1D *hist;
        TCanvas *canvas;
        // End
        
        data->current = 0;
        data->last = 0;
        data->endslice = 0;
        data->channel = int(2*randomvalue);
        snprintf(data->path, sizeof(str), "%s%s/", path, file);
        snprintf(data->name, sizeof(str), "%s", file);
        snprintf(data->ext, sizeof(str), "%s", ext);
        
        // Load to buffer array and fill histograms
        snprintf(data->filename, sizeof(str), "%s%s/%s_Ch_%d_00_Window_Hist.txt", path, file, file, data->channel);
        file_load(data, 3); // Windows file load channel domain
        binstart = (3*data->xaxis[0] - data->xaxis[1])/2;
        binend = data->xaxis[data->last - 1] + (data->xaxis[1] - data->xaxis[0])/2;
        printf("Number of bins: %d\nStart: %.2f\nEnd: %.2f\n", data->last, binstart, binend);
        printf("Number of slices: %d\n", data->endslice);
        for (int i = 0; i < data->endslice; i++){
            snprintf(str, sizeof(str), "Partition_%d_Channel_%d", i, data->channel);
            data->histogram[i] = new TH1D(str, str, data->last, binstart, binend);
            for (int j = 0; j < data->last; j++) data->histogram[i]->Fill(data->xaxis[j], data->hist[i][j]);
            data->histogram[i]->Rebin(rebin);
            data->histogram[i]->Sumw2(kFALSE);
            data->histogram[i]->Sumw2();
        }
        hist = data->histogram[data->cslice];
        std::cout << "Looking for peaks in " << data->cslice << "...\n";
        hist = data->histogram[data->cslice];
        sprintf(str, "%s%s%s%d", "File: ", data->name, " Slice: ", data->cslice);
        canvas = new TCanvas(str, str, 1000, 1000);
        while(true){
            num = peak_off->Search(hist, poff_sig, "", poff_amp);
            hist->Draw("hist");
            hist->GetXaxis()->SetRangeUser(binstart, binend);
            canvas->Update();
            canvas->Update();
            for(int i = 0; i < num; i++) printf("%d\t%0.2f\n", i, peak_off->GetPositionX()[i]);
            std::cout << "Please enter peak you want to rerun\n";
            while (true) {
                std::cin >> attempt;
                if (cin.fail()) {
                    std::cin.clear();
                    std::cin.ignore();
                    std::cout << "Invalid input please enter one of the available options...\n";
                }
                else if (attempt < 0){
                    exit = 1;
                    break;
                }
                else {
                    data->cslice = attempt;
                    break;
                }
            }
            if (exit == 1) break;
        }
        dhost = hist;
    } // Produce heat map of laser data in XIA file
    else if(type ==  2){
        info *data[N_detectors];
        for (int i = 0; i < N_detectors; i++) {
            data[i] = new info;
            data[i]->current = 0;
            data[i]->last = 0;
            data[i]->endslice = 0;
            data[i]->channel = 2*i; //edit this line to reduce number of files needed.
            sprintf(data[i]->filename, "%s%s%s", tpath, file, text);
            sprintf(data[i]->path, "%s", tpath);
            sprintf(data[i]->name, "%s", file);
            sprintf(data[i]->ext, "%s", text);
            /*/
            sprintf(data[i]->filename, "%s%s%s", xpath, xfile, xext);
            sprintf(data[i]->path, "%s", xpath);
            sprintf(data[i]->name, "%s", xfile);
            sprintf(data[i]->ext, "%s", xext);
            /*/
        }
        xiaload(data[0]->filename, data, 1);
        double sigma_value = 40;
        double area_value = 0.005;
        bool override = false;
        
        //*/
        int channel = 3;
        parse_time(data[channel], 30); // Doing 30 minute partitions...
        data[channel]->cslice = 0;
        for (int i = 0; i < NUM_PEAKS*NUM_COLUMN; i++) data[channel]->param[i][0] = -1;
        for (int i = 0; i < NUM_COLUMN; i++) data[channel]->calib[i][0] = -1;
        std::cout << "Total slices: " << data[channel]->endslice << "\n";
        for (int i = data[channel]->cslice; i < data[channel]->endslice; i++) {
            peakfinder(data[channel], &override, &sigma_value, &area_value);
            calibration(data[channel]);
            data[channel]->cslice++;
        }
        //*/
        //fit_save(data[channel], false);
        energy(data[channel]);
        
        xveto(data[channel], data[channel], -1);
        xveto(data[channel], data[channel], 0);
        
        for (int i = 0; i < 4; i++) if (channel != i) {
            std::cout << "Veto Channel: " << data[channel]->channel << " with Channel: " << data[i]->channel << "\n";
            xveto(data[channel], data[i], 4);
            xveto(data[channel], data[i], 1);
        }
        
        sprintf(data[channel]->filename, "%s%s%s%d%s", tpath, file, "_Fitted_CH_", 2*channel, ".txt");
        //file_save(data[channel], data[channel]->filename, 1);
        
        
        
        //status(data[channel], 1);
        //status(data[channel], 2);
        
        /*/ energy bin no longer true, remove later, 20161021
        std::cout << "Time in file: " << (data[channel]->time[data[channel]->last - 1] - data[channel]->time[0])/(3.6*pow(10, 12)) << " hrs\n";
        
        TH2I *hist = new TH2I("Uncalibrated", "Uncalibrated", 8192, 0, 50, 1440, 0, 24);
        ofstream outsavefile;
        sprintf(data[channel]->filename, "%s%s%s", tpath, file, "_Laser_Only.txt");
        outsavefile.open(data[channel]->filename);
        outsavefile << "0_Time (Hrs)\t1_Bin\t2_Energy (eV)\n";
        for (int i = 0; i < data[channel]->last; i++) if (data[channel]->tag[i] == 1 && data[channel]->energy[i] != 0){
            hist->Fill(data[channel]->energy[i], data[channel]->time[i]/(3.6*pow(10, 12)));
            outsavefile << data[channel]->time[i]/(3.6*pow(10, 12)) << "\t";
            outsavefile << data[channel]->bin[i] << "\t";
            outsavefile << data[channel]->energy[i] << "\n";
        }
        outsavefile.close();
        hist->Draw("colz");
        //*/
        
        
    }// Type 2 is reserved for XIA files
    else if(type ==  3){
        /*/
         This section is intended to process a run using the CAEN system, the file will be calibrated,
         the coincidences with high energy events will be marked, and the laser removed from the
         output file. Both a calibration and a peak file will be generated.
        /*/
        // Function Variables:
        int loadtime = 30;              // In hours
        int parsetime = 5;              // In minutes
        int rebinval = 8;                 // Sets values closer to those for XIA
        int counter = 0;
        bool singular = true;           // Use false for moving window average...
        char str[500];
        char *path = cpath;
        char *ext = cext;
        
        TH1I *histo = new TH1I("Trouble", "Trouble", 1000000, 0, 20000000);
        info *data[N_detectors];
        // End
        
        //                      Load Data to buffer for processing
        for (int i = 0; i < N_detectors; i++) {
            data[i] = new info;
            data[i]->current = 0;
            data[i]->last = 0;
            data[i]->endslice = 0;
            data[i]->channel = i; // CAEN system uses 0, 1, 2, 3, ...
            sprintf(data[i]->filename, "%s%s_ch00%d%s", path, file, i, ext);
            sprintf(data[i]->path, "%s", path);
            sprintf(data[i]->name, "%s_ch00%d", file, i);
            sprintf(data[i]->ext, "%s", ext);
            caenload(data[i], -1, 0, loadtime);
            data[i]->current = 0;
        }
        trigger_tag(data[0], data[3]);
        data[3]->current = 0;
        trigger_tag(data[1], data[3]);
        data[3]->last = 0;
        
        for (int i = 0; i < N_detectors; i++) {
            if (i > 1) continue;
            sprintf(data[i]->path, "%s%s/", path, file);
            sprintf(data[i]->name, "%s_Ch_%d_00_Window_Hist", file, 2*i);
            std::cout << "Start:\tFinish\n";
            std::cout << data[i]->time[0]/pow(10, 6) << " us\t";
            std::cout << data[i]->time[data[i]->last - 1]/pow(10, 9) << " s\n";
            for (int j = 0; j < NUM_COLUMN; j++){
                sprintf(str, "Partition_%d_Channel_%d", j, data[i]->channel);
                data[i]->histogram[j] = new TH1D(str, str, tBin, 0, tBin);
                data[i]->histogram[j]->Rebin(rebinval);
            } // Creates all the histograms...
            data[i]->current = 0;
            data[i]->cslice = 0;
            parse_time(data[i], 20, parsetime, singular); // Use false for moving window average...
        }
        for (int i = 0; i < N_detectors; i++) {
            if (i > 1) continue;
            sprintf(data[i]->path, "%s%s/", path, file);
            sprintf(data[i]->name, "%s_Ch_%d_00_Signal_Hist", file, 2*i);
            printf("Working on channel %d\n", data[i]->channel);
            // Flip the signal and laser markers inorder to use the parse function
            for (int j = 0; j < data[i]->last; j++){
                if(data[i]->tag[j] == 0) data[i]->tag[j] = 1; // Signal file is composed of only non-coincidence events
                else if (data[i]->tag[j] == -100) continue;
                else data[i]->tag[j] = 0;
            }
            for (int j = 0; j < NUM_COLUMN; j++) data[i]->histogram[j]->Reset(); // Clear prior
            data[i]->current = 0;
            data[i]->cslice = 0;
            parse_time(data[i], 20, parsetime, singular); // Use false for moving window average...
            printf("Cleaning up histograms\n\n");
            for (int j = 0; j < NUM_COLUMN; j++) data[i]->histogram[j]->Delete(); // Clean up
        }
        printf("Cleaning up\n");
        for (int i = 0; i < N_detectors; i++) delete data[i];
    }// Trouble Shooting:              Load CAEN file and generate signal and laser window files...
    else if(type ==  4){
        /*
            This section is intended to process a run using the XIA system, the file will be calibrated, 
            the coincidences with high energy events will be marked, and the laser removed from the 
            output file. Both a calibration and a peak file will be generated.
        */
        // Function Variables:
        char str[500];
        char *path = tpath;
        char *ext = text;
        bool singular = true;          // Use false for moving window average...
        bool override = false;
        double sigma_value = 5;
        double area_value = 0.005;
        int lal = 200;
        int lcl = 100;
        int ucl = 500;
        int loadtime = 30;               // In hours
        int parsetime = 2;              // In minutes
        int numbin = int(216/0.2);
        double lower = 4;
        double upper = 220;
        info *data[N_detectors];
        // End
        
        for (int i = 0; i < N_detectors; i++) {
            data[i] = new info;
            data[i]->current = 0;
            data[i]->last = 0;
            data[i]->endslice = 0;
            data[i]->channel = 2*i;
            sprintf(data[i]->filename, "%s%s%s", path, file, ext);
            sprintf(data[i]->path, "%s", path);
            sprintf(data[i]->name, "%s", file);
            sprintf(data[i]->ext, "%s", ext);
        } // Load Data to buffer for processing
        xiaload(data[0]->filename, data, 0, loadtime); // Use Full Time Available
        marker(data, 10); // Mark all coincidence events
        
        
        //                      Perform energy calibration on each channel
        for (int i = 0; i < N_detectors; i++) {
            if (i <= 1) continue; // Use to prevent calibration of channel 0 and 2
            sprintf(data[i]->path, "%s%s/", path, file);
            std::cout << "Start:\tFinish\n";
            std::cout << data[i]->time[0]/pow(10, 6) << " us\t";
            std::cout << data[i]->time[data[i]->last - 1]/pow(10, 9) << " s\n";
            diagnostics(data[i], 1);
            for (int j = 0; j < tBin; j++){
                sprintf(str, "Partition_%d_Channel_%d", j, data[i]->channel);
                data[i]->histogram[j] = new TH1D(str, str, tBin, 0, tBin);
            } // Creates all the histograms...
            
            parse_time(data[i], 20, parsetime, singular); // Use false for moving window average...
            std::cout << "Total slices: " << data[i]->endslice << "\n";
            data[i]->cslice = 0;
            for (int j = 0; j < NUM_PEAKS*NUM_COLUMN; j++) data[i]->param[j][0] = -1;
            for (int j = 0; j < NUM_COLUMN; j++) data[i]->calib[j][0] = -1;
            std::cout << "Starting at " << data[i]->cslice << " ending at " << data[i]->endslice << "\n";
            for (int j = data[i]->cslice; j < data[i]->endslice; j++) {
                peakfinder(data[i], &override, &sigma_value, &area_value, false, lcl, ucl, lal);
                calibration(data[i]);
                revise_energy(data[i], tBin, numbin, lower, upper);
                data[i]->cslice++;
            }
            fit_save(data[i], false);
            //energy(data[i], 0, true); // Set to false for moving window...
            diagnostics(data[i], 2);
        }
        
        // Save calibrated data with all events events
        for (int j = 0; j < 4; j++) {
            sprintf(str, "%s%s%s%d%s", path, file, "_Energy_Calibrated_Channel_", 2*j, ".txt");
            file_save(data[j], str, 2); // Save all data including uncalibrated points
        }
        
        // Perform Fit over each partition...
        for (int i = 0; i < N_detectors; i++) {
            if (i <= 1) continue; // Use to prevent calibration of channel...
            //sprintf(data[i]->path, "%s%s/", path, file);
            sprintf(data[i]->name, "%s_Energy_Domain", file);
            std::cout << "Reprocessing data in Energy domain\n";
            std::cout << "Start:\tFinish\n";
            std::cout << data[i]->time[0]/pow(10, 6) << " us\t";
            std::cout << data[i]->time[data[i]->last - 1]/pow(10, 9) << " s\n";
            data[i]->cslice = 0;
            //data[i]->endslice = 0;
            for(int j = 0; j < tBin; j++) data[i]->xaxis[j] = 0;
            for(int j = 0; j < numbin; j++) data[i]->xaxis[j] = lower + (j + 0.5)*(upper - lower)*1.0/numbin;
            hist_save(data[i], numbin, data[i]->endslice);
            for (int j = 0; j < tBin; j++) data[i]->histogram[j]->Delete();  //                                            Remove or comment out when doing energy fit comparison...
                
            /*/
            for (int j = 0; j < data[i]->endslice; j++){
                data[i]->histogram[j]->Delete();
                sprintf(str, "Partition_%d_Channel_%d", j, data[i]->channel);
                data[i]->histogram[j] = new TH1D(str, str, numbin, lower, upper);
                array2hist(data[i], numbin);
                data[i]->cslice++;
            } // Creates all the histograms...
            parse_time(data[i], 0, parsetime, singular, numbin, lower, upper, true);// Use first bool as false for moving window average...
            std::cout << "Total slices: " << data[i]->endslice << "\n";
            sigma_value = 1;
            area_value = 0.005;
            lcl = 10;
            ucl = 100;
            lal = 100;
            data[i]->cslice = (singular ? 0 : 10);
            std::cout << "Starting at " << data[i]->cslice << " ending at " << data[i]->endslice << "\n";
            for (int j = data[i]->cslice; j < data[i]->endslice - (singular ? 0 : 10); j++) {
                peakfinder(data[i], &sigma_value, &area_value, true, lcl, ucl, lal, 1, true);
                calibration(data[i]);
                data[i]->cslice++;
            }
            //*/
            //fit_save(data[i], false);
        }
        for (int i = 0; i < N_detectors; i++) delete data[i];
    }// XIA Full Recoil processing
    else if(type ==  5){
        char str[500];
        char *path = tpath;
        char *ext = text;
        //                      Load Data to buffer for processing
        info *data[N_detectors];
        for (int i = 0; i < N_detectors; i++) {
            data[i] = new info;
            data[i]->current = 0;
            data[i]->last = 0;
            data[i]->endslice = 0;
            data[i]->channel = 2*i; //edit this line to reduce number of files needed.
            sprintf(data[i]->filename, "%s%s%s", path, file, ext);
            sprintf(data[i]->path, "%s", path);
            sprintf(data[i]->name, "%s", file);
            sprintf(data[i]->ext, "%s", ext);
        }
        xiaload(data[0]->filename, data); // Use upto first 21 hrs only
        marker(data, 10); // Mark all coincidence events
        
        
        TH2I *histo[N_detectors];
        for (int i = 0; i < N_detectors; i++) {
            sprintf(str, "%s%d", data[i]->name, 2*i);
            histo[i] = new TH2I(str, str, 144, 0, 12, 1024, 0, 1023);
            globalhist[i] = histo[i];
            data[i]->current = 0;
            while (data[i]->current < data[i]->last) {
                if (!data[i]->tag[data[i]->current]) histo[i]->Fill(data[i]->time[data[i]->current]/(60*pow(10, 9)), data[i]->bin[data[i]->current]);
                data[i]->current++;
            }
        }
    }// XIA Heat Map
    else if(type ==  6){
        char str[500];
        char *path = xpath;
        char *ext = xext;
        TH1D *hist[4];
        TH1D *anti[4];
        ofstream outfile;
        
        //                      Load Data to buffer for processing
        info *data[N_detectors];
        for (int i = 0; i < N_detectors; i++) {
            data[i] = new info;
            data[i]->current = 0;
            data[i]->last = 0;
            data[i]->endslice = 0;
            data[i]->channel = 2*i; //edit this line to reduce number of files needed.
            sprintf(data[i]->filename, "%s%s%s", path, file, ext);
            sprintf(data[i]->path, "%s", path);
            sprintf(data[i]->name, "%s", file);
            sprintf(data[i]->ext, "%s", ext);
        }
        xiaload(data[0]->filename, data, 1, 0); // Use upto first 21 hrs only
        for (int i = 0; i < N_detectors; i++) {
            sprintf(str, "%s%s%s%d%s", path, file, "_Histogram_Ch_", 2*i, ".txt");
            hist[i] = new TH1D(str, str, 1024, 0, 1024);
            sprintf(str, "%s%s%s%d%s", path, file, "_Histogranti_Ch_", 2*i, ".txt");
            anti[i] = new TH1D(str, str, 1024, 0, 1024);
            thist[2*i] = hist[i];
            thist[2*i + 1] = anti[i];
            for (int j = 0; j < data[i]->last; j++) {
                if(data[i]->tag[j]) hist[i]->Fill(data[i]->bin[j]);
                else anti[i]->Fill(data[i]->bin[j]);
            }
        }
        sprintf(str, "%s%s%s", path, file, "_Histogram.txt");
        outfile.open(str);
        outfile << "STJ_0\tSTJ_2\tSTJ_4\tSTJ_6\n";
        for (int j = 0; j < 1024; j++) for (int i = 0; i < N_detectors; i++) {
            outfile << hist[i]->GetBinContent(j + 1) << (i == 3 ? "\n" : "\t");
        }
        outfile.close();
        /*/
        for (int i = 0; i < 4; i++) {
            hist[i]->Delete();
            anti[i]->Delete();
        }
        /*/
        
    }// Bias Vs Responsivity
    else if(type ==  7){
        char str[500];
        char *path = tpath;
        char *ext = text;
        int counter = 0;
        
        info *data[N_detectors];
        for (int i = 0; i < N_detectors; i++) {
            data[i] = new info;
            data[i]->current = 0;
            data[i]->last = 0;
            data[i]->endslice = 0;
            data[i]->channel = 2*i;
            sprintf(data[i]->filename, "%s%s%s", tpath, file, text);
            sprintf(data[i]->path, "%s", tpath);
            sprintf(data[i]->name, "%s", file);
            sprintf(data[i]->ext, "%s", text);
            //caenload(data[i], -1, 1);
        }
        xiaload(data[0]->filename, data, 1);
        marker(data, 10); // Mark all coincidence events
        for (int i = 0; i < N_detectors; i++) {
            sprintf(data[i]->filename, "%s%s_Ch_%d.txt", tpath, file, 2*i);
            file_save(data[i], data[i]->filename, 2);
        }
        
    }// XIA Reprocess data to mark recoil events...
    else if(type ==  8){
        
        // Function Variables:
        int channel = sig;
        int counter = 0;
        int ctwo = 0;
        int window = 10;
        int recounter = 0;
        int acounter = 0;
        int *(pval[2]);
        int valone[2] = {4600, 5500};
        int valtwo[2] = {-3000, 3000};
        char str[500];
        bool gate = false;
        long long lineup[4];
        long long timeago = 0;
        long long timerange = 10*pow(10, 6); // This is a 10 ms time window...
        long long t20 = 0;
        long long t64 = 0;
        TH1D *recoil = new TH1D("Recoil", "Recoil", 512, 0, 3);
        ofstream recoilfile;
        ofstream outfile;
        // End
        
        info *data[N_detectors];
        for (int i = 0; i < N_detectors; i++) {
            data[i] = new info;
            data[i]->current = 0;
            data[i]->last = 0;
            data[i]->endslice = 0;
            data[i]->channel = 2*i; //edit this line to reduce number of files needed.
            sprintf(data[i]->filename, "%s%s%s", tpath, file, text);
            sprintf(data[i]->path, "%s", tpath);
            sprintf(data[i]->name, "%s", file);
            sprintf(data[i]->ext, "%s", text);
        }
        xiaload(data[0]->filename, data, 0, 5); // Use upto first 10 hrs only
        
        channel = sig;
        for (int i = 0; i < 4; i++) {
            xveto(data[i], data[i], -1);    // Set all veto, recoil and high to true (this means use all points)
            xveto(data[i], data[i], -3);    // Set all high to false
            xveto(data[i], data[i], -4);    // Set all recoil to false
            xveto(data[i], data[i], 0);     // Set all vetos for laser events false
            data[i]->current = 0;
        }
        
        sprintf(str, "%s%s%s", tpath, file, "_Coincidence.txt");
        outfile.open(str);
        outfile << "0_Bin\t0_Time\t2_Bin\t2_Time\t4_Bin\t4_Time\t6_Bin\t6_Time\n";
        for (int i = 0; i < data[channel]->last; i++) {
            data[channel]->current = i;
            std::cout << i << "\r";
            /*/
             Check that points are in coincidence and none of the points are laser events, mark
             veto column false to indicate a coincidence event.
             Note:
             8000 <= ##
             0 < # < 8000
             /*/
            if (!data[channel]->tag[i] &&
                ptcheck(data[channel], data[(channel + 1)%4]) &&
                !data[(channel + 1)%4]->tag[data[(channel + 1)%4]->current] &&
                ptcheck(data[channel], data[(channel + 2)%4]) &&
                !data[(channel + 2)%4]->tag[data[(channel + 2)%4]->current] &&
                ptcheck(data[channel], data[(channel + 3)%4]) &&
                !data[(channel + 3)%4]->tag[data[(channel + 3)%4]->current]) {
                for (int j = 0; j < 4; j++) lineup[j] = data[j]->time[data[j]->current];
                if (findmax(lineup, 4) - findmin(lineup, 4) < window*pow(10, 3)) {
                    for (int j = 0; j < 4; j++) {
                        data[j]->veto[data[j]->current] = false;
                        outfile << data[j]->bin[data[j]->current] << "\t";
                        outfile << data[j]->time[data[j]->current] << (j == 3 ? "\n" : "\t");
                        counter++;
                    } // Coincidence condition
                    // Tag alpha events specifically know as (0, 0, 0, 0) and (##, ##, ##, ##)
                    if ((8000 <= data[0]->bin[data[0]->current] &&
                         8000 <= data[1]->bin[data[1]->current] &&
                         8000 <= data[2]->bin[data[2]->current] &&
                         8000 <= data[3]->bin[data[3]->current]) ||
                        (!data[0]->bin[data[0]->current] &&
                         !data[1]->bin[data[1]->current] &&
                         !data[2]->bin[data[2]->current] &&
                         !data[3]->bin[data[3]->current])) {
                            for (int j = 0; j < 4; j++)  data[j]->high[data[j]->current] = true; // Primary Alpha condition
                            ctwo++;
                        }
                    // Tag recoils and some alphas from (0, #, #, #) and (##, #, #, #)
                    else for (int j = 0; j < 4; j++) if((!data[j]->bin[data[j]->current] ||  8000 <= data[j]->bin[data[j]->current]) &&
                            (data[(j + 1)%4]->bin[data[(j + 1)%4]->current] &&
                             data[(j + 1)%4]->bin[data[(j + 1)%4]->current] < 8000 &&
                             data[(j + 2)%4]->bin[data[(j + 2)%4]->current] &&
                             data[(j + 2)%4]->bin[data[(j + 2)%4]->current] < 8000 &&
                             data[(j + 3)%4]->bin[data[(j + 3)%4]->current] &&
                             data[(j + 3)%4]->bin[data[(j + 3)%4]->current] < 8000)) {
                                if (data[0]->recoil[data[0]->current] || data[0]->high[data[0]->current]) {
                                    std::cout << "Detector: " << j << "\t";
                                    std::cout << "STJ 0:\t" << data[0]->bin[data[0]->current] << "\t";
                                    std::cout << "STJ 2:\t" << data[1]->bin[data[1]->current] << "\t";
                                    std::cout << "STJ 4:\t" << data[2]->bin[data[2]->current] << "\t";
                                    std::cout << "STJ 6:\t" << data[3]->bin[data[3]->current] << "\n";
                                }
                                if (j == 0 || j == 1) {pval[0] = valone; pval[1] = valtwo;} // STJ 0 and 2 x{4.6, 5.5} and y{-3, 3};
                                else {pval[0] = valtwo; pval[1] = valone;} // STJ 4 and 6 x{-3, 3} and y{4.6, 5.5}
                                t20 = data[1]->time[data[1]->current] - data[0]->time[data[0]->current];
                                t64 = data[2]->time[data[2]->current] - data[3]->time[data[3]->current];
                                if ((j == 0 && (pval[0][0] <=  t20 &&  t20 <= pval[0][1] && pval[1][0] <=  t64 &&  t64 <= pval[1][1])) ||
                                    (j == 1 && (pval[0][0] <= -t20 && -t20 <= pval[0][1] && pval[1][0] <=  t64 &&  t64 <= pval[1][1])) ||
                                    (j == 2 && (pval[0][0] <=  t20 &&  t20 <= pval[0][1] && pval[1][0] <=  t64 &&  t64 <= pval[1][1])) ||
                                    (j == 3 && (pval[0][0] <=  t20 &&  t20 <= pval[0][1] && pval[1][0] <= -t64 && -t64 <= pval[1][1]))){
                                    recounter++;
                                    for (int k = 0; k < 4; k++) data[k]->recoil[data[k]->current] = true;
                                } // Recoil condition
                                else {
                                    ctwo++;
                                    for (int k = 0; k < 4; k++) data[k]->high[data[k]->current] = true;
                                }// Secondary Alpha condition
                                
                                //std::cout << "x0 < t20 < x1:\t" << pval[0][0] << " <\t" << t20 << "\t< " << pval[0][1] << "\t";
                                //std::cout << "y0 < t64 < y1:\t" << pval[1][0] << " <\t" << t64 << "\t< " << pval[1][1] << "\t\t";
                                //std::cout << j << "\t" << data[0]->recoil[data[0]->current] << "\n";
                                //break;
                            }
                }
                if (data[channel]->high[data[channel]->current]) acounter++;
            }
        }
        outfile.close();
        std::cout << "\nNumber of coincidence Events: " << counter << "\n";
        std::cout << "\nNumber of Recoil Events: " << recounter << "\n";
        std::cout << "\nNumber of Alpha Events: " << acounter << "\n";
        std::cout << "\nNumber of ctwo: " << ctwo << "\n";
        for (int j = 0; j < 4; j++) {
            sprintf(str, "%s%s%s%d%s", tpath, file, "_Energy_Calibrated_Channel_", 2*j, ".txt");
            file_save(data[j], str, 1);
        }
    }// Recoil, alpha, and coin marker tester
    else if(type ==  9){
        /*
         This section will save the data from the xia into ascii
         
         */
        // Function Variables:
        char str[500];
        char *path = tpath;
        char *ext = text;
        double loadtime = 20;                                  // In hours
        float parsetime = (randomvalue == 0 ? 5 : randomvalue); // minutes
        //int det_channels[32] = {0, 2, 4, 6};
        //int det_channels[32] = {6, 8, 12, 14};
        int markerrange = 10;                                   // Range to mark coincidence in microsecond...
        int det_channels[32] = {0, 2, 4, 6, 8, 10, 12, 14};
        info *data[N_detectors];
        // End
        
        printf("\n\n\nNOTICE NOTICE NOTICE\n\n\n");
        printf("Partition size: %.2f min, ", parsetime);
        printf("Load time: %0.2f hrs, ", loadtime);
        printf("Coincidence range: %d us", markerrange);
        printf("\n\n\nNOTICE NOTICE NOTICE\n\n\n");
        
        for (int i = 0; i < N_detectors; i++) {
            data[i] = new info;
            data[i]->current = 0;
            data[i]->last = 0;
            data[i]->endslice = 0;
            data[i]->channel = det_channels[i];
            sprintf(data[i]->filename, "%s%s%s", path, file, ext);
            sprintf(data[i]->path, "%s", path);
            sprintf(data[i]->name, "%s", file);
            sprintf(data[i]->ext, "%s", ext);
        } // Load Data to buffer for processing
        xiaload(data[0]->filename, data, 0, loadtime); // Use Full Time Available
        
        for (int i = 0; i < N_detectors; i++){
            //info* pointer, char* new_name, int type = 0, int start = 0, int stop = N_Events
            //str to generate file name
            snprintf(str, sizeof(str), "%s%s/%s_Ch_%d_000_ASCII.txt", path, file, file, data[i]->channel); // Window file name..
            file_save(data[i], str, 0);
        }
        
        
        
        
        
        
    } // Trouble Shooting:              To save XIA to ascii
    else if(type == 10){
    } // Trouble Shooting:              Empty
    else if(type == 11){
        /*
         This section will run a fit on the fully energy calibrated file to confirm position of peaks
         and acts as a sanity check for the calibration...
         
         Loads the calibrated window file for processing, cannot use individual list mode points due 
         to addition of histogram error...
         */
        // Function Variables:
        char str[500];
        char *path = tpath;
        char *ext = text;
        bool override = false;
        double sigma_value = 1;
        double area_value = 0.005;
        int lal = 100;
        int lcl = 6;
        int ucl = 100;
        info *data = new info;
        // End
        
        data->current = 0;
        data->last = 0;
        data->endslice = int(randomvalue);
        data->channel = 4;
        sprintf(data->filename, "%s%s_Ch_%d_05_Energy_Domain_Partition_Hist.txt", path, file, data->channel);
        //sprintf(data->filename, "%s%s_Energy_Domain_Partition_Energy_Domain_Partition_Energy_Domain_Moving_Window_Hist_Ch_%d.txt", path, file, data->channel);
        
        sprintf(data->path, "%s", path);
        sprintf(data->name, "%s_Full_Laser_Fit", file);
        sprintf(data->ext, "%s", ext);
        file_load(data, 3);
        data->cslice = 1;
        data->histogram[0] = new TH1D("Full", "Full", int(216/0.2), 4, 220);
        printf("Number of line in array: %d\n", data->last);
        for (int i = 0; i < data->endslice; i++) {
            sprintf(str, "Partition_%d_Channel_%hd", i, data->channel);
            data->histogram[data->cslice] = new TH1D(str, str, int(216/0.2), 4, 220);
            array2hist(data);
            data->histogram[0]->Add(data->histogram[data->cslice]);
            data->cslice++;
        }
        
        //*/
        for (int j = 0; j < NUM_PEAKS*NUM_COLUMN; j++) data->param[j][0] = -1;
        for (int j = 0; j < NUM_COLUMN; j++) data->calib[j][0] = -1;
        data->cslice = 0;
        data->endslice = 1;
        for (int i = 0; i < data->endslice; i++){
            peakfinder(data, &override, &sigma_value, &area_value, true, lcl, ucl, lal, 1, false);
            calibration(data);
            data->cslice++;
        }
        fit_save(data, false);
        //*/
    } // Process Full Laser Energy Domain...
    else if(type == 12){
        /*
         This section will test if all of the laser in a file has been properly identified...
         */
        // Function Variables:
        char str[500];
        char *path = tpath;
        char *ext = text;
        int loadtime = 30;                  // In hours
        float parsetime = 6;                // In minutes ... Had this on 6.49 for awhile do to broad laser test file processing was great...Instead of 2 minutes going to run 5 minutes...
        int det_channels[32] = {0, 2, 4, 6};
        int markerrange = 10;               // Range to mark coincidence in microsecond...
        int rebinval = 1;
        int range = 10;
        int tcounter = 0;
        int lcounter = 0;
        long long ptime = 0;
        long long lentime = 100; // This is in msec
        //int det_channels[32] = {0, 2, 4, 6, 8, 10, 12, 14};
        info *data[N_detectors];
        
        TH1D *las[N_detectors];
        TH1D *sig[N_detectors];
        TH1D *hist[N_detectors];
        
        // End
        
        printf("\n\n\nNOTICE NOTICE NOTICE\n\n\n");
        printf("Partition size: %.2f min, ", parsetime);
        printf("Load time: %d hrs, ", loadtime);
        printf("Rebin Value: %d, ", rebinval);
        printf("Coincidence range: %d us", markerrange);
        printf("\n\n\nNOTICE NOTICE NOTICE\n\n\n");
        
        for (int i = 0; i < N_detectors; i++) {
            data[i] = new info;
            data[i]->current = 0;
            data[i]->last = 0;
            data[i]->endslice = 0;
            data[i]->channel = det_channels[i];
            sprintf(data[i]->filename, "%s%s%s", path, file, ext);
            sprintf(data[i]->path, "%s", path);
            sprintf(data[i]->name, "%s_Ch_%d", file, 2*i);
            sprintf(data[i]->ext, "%s", ext);
        } // Load Data to buffer for processing
        /*/
         A straight load from the XIA file will mark all lasers that fell in range of the tag
         window with a 1 in ->tag. Next I would like to use marker to go into more detail but 
         need to work that out more here... For that purpose I'll be using ptcheck to determine
         coincidence events. Keep in mind that not all detectors may see the laser there is a
         none zero possibility of observing no photons. So need to step confirmation... First
         locate all points that are in coincidence. Then determine which of these points fall 
         into the period of laser using known laser points to reset to zero periodically...
        /*/
        xiaload(data[0]->filename, data, 0, loadtime); // Use Full Time Available
        //marker(data, markerrange);
        for (int i = 0; i < N_detectors; i++){
            data[i]->current = 0;
            while(data[i]->current < data[i]->last){
                for (int j = 0; j < N_detectors; j++) if (j != i && ptcheck(data[i], data[j], range)){
                    data[i]->id[data[i]->current]++;
                    data[j]->id[data[j]->current]++;
                }// This will increase the id count for coincidence events including laser events...
                data[i]->current++;
            }
        }
        
        // Create output histograms...
        for (int i = 0; i < N_detectors; i++) {
            snprintf(str, sizeof(str), "No_ID_%d", data[i]->channel);
            las[i] = new TH1D(str, str, 8192, 0, 10.24); // Units are in milliseconds need to divide the inputs values to get correct units... -> 10^-6
            snprintf(str, sizeof(str), "ID_%d", data[i]->channel);
            sig[i] = new TH1D(str, str, 8192, 0, 10.24);
            snprintf(str, sizeof(str), "signal_%d", data[i]->channel);
            hist[i] = new TH1D(str, str, 8192, 0, 8192);
        }
        
        for (int i = 0; i < N_detectors; i++){
            data[i]->current = 0;
            for (int j = 0; j < data[i]->last - 1; j++){
                if (data[i]->time[j + 1] == data[i]->time[j] && data[i]->bin[j] == 0) break;
                data[i]->current++;
            }
            
            printf("Detector %d\tStarting at %d of %d\nTag\tBin\tTime\n", data[i]->channel, data[i]->current, data[i]->last);
            //for (int j = (data[i]->current < 10 ? 0 : data[i]->current - 10); j < data[i]->current + 10; j++) printf("%d\t%d\t%lld\n", data[i]->tag[j], data[i]->bin[j], data[i]->time[j]);
            for (int j = (data[i]->current < 10 ? 0 : data[i]->current - 10); j < data[i]->current + 10; j++) printf("%d\t%f\t%lld\n", data[i]->tag[j], data[i]->bin[j], data[i]->time[j]);
        }
        
        
        //
        for (int i = 0; i < N_detectors; i++) {
            las[i]->Reset();
            sig[i]->Reset();
            ptime = 0;
            data[i]->current = 0;
            while (data[i]->tag[data[i]->current] != 1) data[i]->current++;
            ptime = data[i]->time[data[i]->current];
            printf("Starting at %d out of %d for detector %d\n", data[i]->current, data[i]->last, data[i]->channel);
            printf("Number of out of order points\n0");
            tcounter = 0;
            lcounter = 0;
            for (int j = 1; j < data[i]->last; j++){
                if(data[i]->tag[j] != 1) hist[i]->Fill(data[i]->bin[j]);
                //if(j != data[i]->last - 1 && data[i]->tag[j] != 1 && (data[i]->time[j] - ptime)*pow(10, -6) < lentime) {
                las[i]->Fill((data[i]->time[j] - data[i]->time[j - 1])*pow(10, -6));
                //    tcounter++;
                //}
                if(data[i]->time[j] < data[i]->time[j - 1]){
                    tcounter++;
                    printf("%d\r", tcounter);
                }
                if(data[i]->tag[j] != 1) lcounter++;
                if(2 <= data[i]->id[j] || 2 <= data[i]->id[j - 1]){// && (data[i]->time[j] - ptime)*pow(10, -6) < lentime) {
                    sig[i]->Fill((data[i]->time[j] - data[i]->time[j - 1])*pow(10, -6));
                }
                else if (data[i]->tag[j] == 1) ptime = data[i]->time[j];
                else if ((data[i]->time[j] - ptime)*pow(10, -6) > lentime) {
                    while (data[i]->tag[j] != 1) j++;
                }
            }
            printf("\n");
            printf("Number of none laser points: %d\n", lcounter);
            tcounter = 0;
            for(int j = 0; j < las[i]->GetNbinsX(); j++){
                data[i]->xaxis[j] = las[i]->GetBinCenter(j + 1);
                data[i]->hist[0][j] = las[i]->GetBinContent(j + 1);
                data[i]->hist[1][j] = sig[i]->GetBinContent(j + 1);
            }
            data[i]->endslice = 2;
            //window_save(data[i], 15, las[i]->GetNbinsX());
            hist[i]->Rebin(8);
            if (i == 0) hist[i]->Draw();
            else hist[i]->Draw("same");
            
            //delete data[i];
        }
        //*/
        
        
        
    } // Trouble Shooting:              Post-Processing Laser Tagging...
    else if(type == 13){
        /*
         This section will load the binary xia file and create an ascii file.
         */
        // Function Variables:
        info *detector[N_detectors];
        char *path = tpath;
        int exit = 0;
        // End
        
        for (int i = 0; i < N_detectors; i++) {
            detector[i] = new info;
            detector[i]->current = 0;
            detector[i]->last = 0;
        }
        
        if (sig == 0) sig = 1; // default 1 hour conversion...
        printf("Loading data from file\n");
        sprintf(detector[0]->filename, "%s%s%s", path, file, ".bin");
        exit = xiaload(detector[0]->filename, detector, 0, sig);
        for(int i = 0; i < N_detectors; i++) {
            if (exit == 1 && detector[i]->last > 0) {
                printf("Detector number:\t%d\tLast:\t%d\n", i, detector[i]->last);
                sprintf(detector[i]->filename, "%s%s/%s_Ch_%d.txt", path, file, file, i);
                exit = file_save(detector[i], detector[i]->filename);
            }
        }
    } // Trouble Shooting:              Create ascii file of sig time...
    else if(type == 14){
        info *data = new info;
        TH1D *hist = new TH1D("Main", "Main", 11, -0.45, 1.2);
        hist->Fill(0.00, 15);
        hist->Fill(0.15, 30);
        hist->Fill(0.30, 45);
        hist->Fill(0.45, 30);
        hist->Fill(0.60, 15);
        data->endslice = 1;
        data->cslice = 0;
        data->histogram[0] = hist;
        hist->Draw();
        data->calib[0][1] = 1;
        data->calib[0][3] = 0;
        revise_energy(data, 11, 4, 1./6, 2./3, true);
        //*/
        printf("\n\n\n");
        revise_energy(data, 11, 4, 0, 0.8, true);
        printf("\n\n\n");
        revise_energy(data, 11, 3, 0, 0.75, true);
        printf("\n\n\n");
        revise_energy(data, 11, 8, -0.15, 1.05, true);
        printf("\n\n\n");
        revise_energy(data, 11, 12, -0.2, 1, true);
        //*/
    } // Trouble Shooting:              revise_energy test
    else if(type == 15){
        /*
         This section will load the binary xia file and create the window/partition histograms for 
         both the laser and signal generating two window/partition histogram files for quick us later.
         
         Note that a parse time of 5 minutes has been used with the exception of the Dec data were the laser
         has a 6 minute period...
         
         randomvalue    =>      is used to change the parsetime for partition size
         sig            =>      is use to set the rebin value of the histogram
         
         */
        // Function Variables:
        char str[500];
        char *path = tpath;
        char *ext = text;
        bool singular = true;                                   // Use false for moving window average...
        double loadtime = 20;                                  // In hours
        float parsetime = (randomvalue == 0 ? 5 : randomvalue); // minutes
        //int det_channels[32] = {0, 2, 4, 6};
        //int det_channels[32] = {6, 8, 12, 14};
        int markerrange = 10;                                   // Range to mark coincidence in microsecond...
        int rebinval = (sig <= 0 ? 1 : sig);
        int det_channels[32] = {0, 2, 4, 6, 8, 10, 12, 14};
        info *data[N_detectors];
        // End
        
        printf("\n\n\nNOTICE NOTICE NOTICE\n\n\n");
        printf("Partition size: %.2f min, ", parsetime);
        printf("Load time: %0.2f hrs, ", loadtime);
        printf("Rebin Value: %d, ", rebinval);
        printf("Coincidence range: %d us", markerrange);
        printf("\n\n\nNOTICE NOTICE NOTICE\n\n\n");
        
        for (int i = 0; i < N_detectors; i++) {
            data[i] = new info;
            data[i]->current = 0;
            data[i]->last = 0;
            data[i]->endslice = 0;
            data[i]->channel = det_channels[i];
            sprintf(data[i]->filename, "%s%s%s", path, file, ext);
            sprintf(data[i]->path, "%s", path);
            sprintf(data[i]->name, "%s", file);
            sprintf(data[i]->ext, "%s", ext);
        } // Load Data to buffer for processing
        xiaload(data[0]->filename, data, 0, loadtime); // Use Full Time Available
        marker(data, markerrange); // Mark all coincidence events
        
        // Parse data into individual files...
        for (int i = 0; i < N_detectors; i++) {
            sprintf(data[i]->path, "%s%s/", path, file);
            sprintf(data[i]->name, "%s_Ch_%d_00_Window_Hist", file, 2*i);
            std::cout << "Start:\tFinish\n";
            std::cout << data[i]->time[0]/pow(10, 6) << " us\t";
            std::cout << data[i]->time[data[i]->last - 1]/pow(10, 9) << " s\n";
            for (int j = 0; j < NUM_COLUMN; j++){
                sprintf(str, "Partition_%d_Channel_%d", j, data[i]->channel);
                data[i]->histogram[j] = new TH1D(str, str, tBin, 0, tBin);
                data[i]->histogram[j]->Rebin(rebinval);
            } // Creates all the histograms...
            data[i]->numbin = data[i]->histogram[0]->GetNbinsX();
            printf("Number of bins in each histogram for detector %d is %d\n", data[i]->channel, data[i]->numbin);
            data[i]->current = 0;
            data[i]->cslice = 0;
            parse_time(data[i], 20, parsetime, singular); // Use false for moving window average...
        }
        for (int i = 0; i < N_detectors; i++) {
            sprintf(data[i]->path, "%s%s/", path, file);
            sprintf(data[i]->name, "%s_Ch_%d_00_Signal_Hist", file, 2*i);
            // Temporarily flip the signal and laser markers inorder to use the parse function
            for (int j = 0; j < data[i]->last; j++){
                if(data[i]->id[j] == 0) data[i]->tag[j] = 1; // Signal file is composed of only non-coincidence events
                else data[i]->tag[j] = 0;
            }
            for (int j = 0; j < NUM_COLUMN; j++) data[i]->histogram[j]->Reset(); // Clear prior
            data[i]->current = 0;
            data[i]->cslice = 0;
            parse_time(data[i], 20, parsetime, singular); // Use false for moving window average...
            // Clean Up...
            //chost = data[i]->histogram[0];
            for (int j = 0; j < NUM_COLUMN; j++) data[i]->histogram[j]->Delete();
            delete data[i];
        }
    } // Trouble Shooting:              Load binary file and generate signal and laser window files...
    else if(type == 16){
        /*/                                     NOTES (20160914)
            This section is a trouble shooting piece for the peakfinder algorithm allowing quick tests
            on ideas such as lower and upper centroid limitations, peak finding values and so on.
            This is final and should require no additional changes unless other sections of the code have
            been changed and features have been altered or removed.
         
            pointer->last from the window load will indicate the number of bins...
         
            Output is the parameter file...
            
            The values of lal, lcl, and ucl may need to be adjusted from experiment to experiment depending 
            on acquisition settings.
            
            20161213
            For the partition width uncertainty quantification the lal was set to 100 for all peaks and the 
            range was set to between 300 and 800 bins which correspondes to roughly to 40 - 100 eV
         
            Use sig to toogle trouble....when user starts from anywhere other then 0...
         /*/
        // Function Variables:
        char str[500];
        char *path = tpath;
        char *ext = text;
        bool override = false; //(sig < 0 ? false : true);
        bool trouble = false; //(sig == 0 ? false : true);
        double sigma_value = 5; // Change back to 15 for ORTEC files
        double area_value = 0.001;// 0.001
        int rebin = int(val1 = 0 ? 1 : TMath::Abs(val1));
        int lal = 1000;
        int lcl = 10; //40; //20;
        int ucl = 3000; //8190; //150;
        int pspace = 15; //40;//10
        int counter = 0;
        double binstart = 0;
        double binend = 0;
        info *data = new info;
        // End
        
        data->current = 0;
        data->last = 0;
        data->endslice = 0;
        data->channel = int(2*randomvalue);
        snprintf(data->path, sizeof(str), "%s%s/", path, file);
        snprintf(data->name, sizeof(str), "%s", file);
        snprintf(data->ext, sizeof(str), "%s", ext);
        
        // Load to buffer array and fill histograms
        snprintf(data->filename, sizeof(str), "%s%s/%s_Ch_%d_00_Window_Hist.txt", path, file, file, data->channel); // Window file name...
        file_load(data, 3); // Windows file load channel domain
        binstart = (3*data->xaxis[0] - data->xaxis[1])/2;
        binend = data->xaxis[data->last - 1] + (data->xaxis[1] - data->xaxis[0])/2;
        printf("\nNumber of bins: %d\nStart: %.2f\nEnd: %.2f\n", data->last, binstart, binend);
        printf("Number of slices: %d\n", data->endslice);
        for (int j = 0; j < NUM_PEAKS*NUM_COLUMN; j++) data->param[j][0] = -1;
        for (int i = 0; i < data->endslice; i++){
            snprintf(str, sizeof(str), "Partition_%d_Channel_%d", i, data->channel);
            data->histogram[i] = new TH1D(str, str, data->last, binstart, binend); // Calibration invariant
            for (int j = 0; j < data->last; j++){
                data->histogram[i]->Fill(data->xaxis[j], data->hist[i][j]);
            }
            data->histogram[i]->Rebin(rebin);
            data->histogram[i]->Sumw2(kFALSE);
            data->histogram[i]->Sumw2();
        }
        data->cslice = TMath::Abs(sig);
        if (data->cslice != 0){
            counter = data->endslice;
            sprintf(data->filename, "%s%s/%s_Ch_%d_01_Parameters.txt", path, file, file, data->channel);
            sprintf(data->path, "%s%s/", path, file);
            sprintf(data->name, "%s", file);
            sprintf(data->ext, ".txt");
            file_load(data, 5);
            //This change was needed because of changes in calibration that looks back at older calibration calculations to calculate the next term
            data->cslice = 0; //TMath::Abs(sig) - 1;
            while (data->cslice < TMath::Abs(sig)){
                calibration(data, lal, false, trouble);
                data->cslice++;
            }
            //data->cslice = TMath::Abs(sig);
            data->endslice = counter;
            counter = 0;
            printf("Total slices: %d\n", data->endslice);
        }
        data->indexvalue = 1;
        for (int j = TMath::Abs(sig); j < data->endslice; j++) {
            peakfinder(data, &override, &sigma_value, &area_value, false, lcl, ucl, lal, pspace, false, 2, 2, trouble);
            data->cslice++;
            if (sig < 0) break;
            fit_save(data, false, 1);
        }
        printf("\nNumber of bins: %d\nStart: %.2f\nEnd: %.2f\n", data->last, binstart, binend);
        printf("Number of slices: %d\n", data->endslice);
        data->indexvalue = 1;
        fit_save(data, false, 1);
        // Clean up
        for (int i = 0; i < data->endslice; i++) data->histogram[i]->Delete();
        delete data;
        
    } // Trouble Shooting:              Load partition/window and run peakfinder(...)
    else if(type == 17){
        /*/                                     NOTES (20160914)
         This section is a trouble shooting piece for the calibration algorithm allowing quick tests
         on ideas such as area or sigma limitations on the calibration outputs. This is final and
         should require no additional changes unless other sections of the code have been changed
         and features have been altered or removed.
        
         Output: calibration and residual file...
         
         The lal variable may be adjusted depending on the needs of data. Based on a plot of dC vs Area the ideal 
         limit is peaks with areas above 800 counts but 700 could be used if a lower limited is needed... The choice 
         here is somewhat arbitrary
         
        /*/
        
        // Initialize Variables
        int counter = 0;
        int lal = 200;
        float lstart = (sig == 0 ? 66.5 : -sig);
        bool trouble = false;
        bool specialize = false;
        char str[500];
        char *path = tpath;
        info *data = new info;
        // End
        
        // General data object setup
        data->current = 0;
        data->last = 0;
        data->endslice = 0;
        data->channel = 2*int(randomvalue);
        sprintf(data->filename, "%s%s/%s_Ch_%d_01_Parameters.txt", path, file, file, data->channel);
        sprintf(data->path, "%s%s/", path, file);
        sprintf(data->name, "%s", file);
        sprintf(data->ext, ".txt");
        file_load(data, 5); // Parameter file load
        data->cslice = 0;
        for (int i = 0; i < data->endslice; i++){
            calibration(data, lal, specialize, trouble, lstart);
            data->indexvalue = (sig == 0 ? 2 : sig);
            residual(data);
            data->cslice++;
        }
        data->indexvalue = (sig == 0 ? 3 : sig);
        fit_save(data, false, 2);
        printf("Performing cleanup\n\n\n");
        delete data;
    } // Trouble Shooting:              Load Parameter file and run calibration(...) and residual(...)
    else if(type == 18){
        /*/                                     NOTES (20161012)
         This section is a trouble shooting piece for the creation of the revised_energy partitions and eventually (if needed)
         the construction of the moving window averaging (currently commented out).
         
         Output is the energy domain histograms (decimal valued) betweeen the energy set by lower and upper and number of bins
         set by numbin, default is 4 -> 220 eV with 216 bins at 0.2 eV width. This range was choosen because it allows the use of bins
         sizes : 0.1, 0.2, 0.3, 0.4, 0.5. The range could be reduced 4 -> 112 eV with 108 bins and still achieve the same 
         versatility, the additional range is choosen to make sure nothing is occuring in this range.
         
         20161213
         Updated the numbin value to be set by int values of sig with default size of 0.2 eV
         
         /*/
        
        // Function Variables:
        int counter = 0;
        double binstart = 0;
        double binend = 0;
        int numbin = 2*2160/(sig <= 0 ? 2 : sig); //Use 1, 2, 3, 4, 5 for 0.1, 0.2, 0.3, 0.4, 0.5 eV/bin
        //int numpart = 30; // Moving window average is this value at 2 minutes each...
        char str[500];
        char word[500];
        char *path = tpath;
        char *ext = text;
        double lower = 0;
        double upper = 2*216 + lower;
        info *data = new info;
        TH1D *hist = new TH1D("Main", "Main", numbin, lower, upper);
        // End
        
        printf("Converting to Number of bins: %d\nStart: %.2f\nEnd: %.2f\n", numbin, lower, upper);
        
        // Set info general values...
        data->current = 0;
        data->last = 0;
        data->endslice = 0;
        data->channel = int(2*randomvalue);
        snprintf(data->path, sizeof(data->path), "%s%s/", path, file);
        snprintf(data->name, sizeof(data->name), "%s", file);
        snprintf(data->ext, sizeof(data->ext), ".txt");
        
        // Load histograms
        snprintf(data->filename, sizeof(data->filename), "%s%s/%s_Ch_%d_00_Window_Hist.txt", path, file, file, data->channel); // Window file name...
        file_load(data, 3);
        binstart = (3*data->xaxis[0] - data->xaxis[1])/2;
        binend = data->xaxis[data->last - 1] + (data->xaxis[1] - data->xaxis[0])/2;
        for (int j = 0; j < NUM_PEAKS*NUM_COLUMN; j++) data->param[j][0] = -1;
        for (int i = 0; i < data->endslice; i++){
            snprintf(str, sizeof(str), "Partition_%d_Channel_%d", i, data->channel);
            data->histogram[i] = new TH1D(str, str, data->last, binstart, binend); // Calibration invariant
            if (i == 0) {
                printf("\nNumber of bins: %d\nStart: %.2f\nEnd: %.2f\n", data->last, binstart, binend);
                printf("Number of slices: %d\n", data->endslice);
            }
            for (int j = 0; j < data->last; j++){
                data->histogram[i]->Fill(data->xaxis[j], data->hist[i][j]);
            }
        } // Generate BIN domain histograms
        
        // Load calibration file
        snprintf(str, sizeof(word), "");
        if(sig < 0) {
            snprintf(word, sizeof(word), "_%d-%d", int(-sig), int(-sig) + 21);
            snprintf(data->filename, sizeof(data->filename), "%s%s/%s_Ch_%d_03_Calibration%s.txt", path, file, file, data->channel, word);
        }
        else snprintf(data->filename, sizeof(data->filename), "%s%s/%s_Ch_%d_03_Calibration.txt", path, file, file, data->channel);
        
        file_load(data, 4);
        
        // Run revised energy to transform partition bin basis into energy basis...
        printf("Starting conversion to Number of bins: %d\n", numbin);
        data->cslice = 0;
        for (int j = 0; j < data->endslice; j++) {
            revise_energy(data, tBin, numbin, lower, upper);
            data->cslice++;
        }
        printf("End conversion to Number of bins: %d\n", numbin);
        // Save energy domain partitions under new file name...
        snprintf(str, sizeof(word), "");
        if(sig < 0) {
            snprintf(word, sizeof(word), "_%d-%d", int(-sig), int(-sig) + 21);
            snprintf(data->name, sizeof(data->name), "%s_Ch_%d_04_Energy_Domain_Partition_Hist%s", file, data->channel, word);
        }
        else snprintf(data->name, sizeof(data->name), "%s_Ch_%d_04_Energy_Domain_Partition_Hist", file, data->channel);
        for(int j = 0; j < tBin; j++) data->xaxis[j] = 0;
        for(int j = 0; j < numbin; j++) data->xaxis[j] = lower + (j + 0.5)*(upper - lower)*1.0/numbin;
        data->cslice = 0;
        printf("Number of columns In: %d\n", numbin);
        window_save(data, 0, numbin);
        /*/
        // Recreate partitions in energy domain
        for (int i = 0; i < data->endslice; i++){
            snprintf(str, sizeof(str), "Partition_%d_Channel_%d", i, data->channel);
            data->histogram[i]->Delete();
            data->histogram[i] = new TH1D(str, str, numbin, lower, upper);
            for (int j = 0; j < data->last; j++){
                data->histogram[i]->Fill(data->xaxis[j], data->hist[i][j]);
            }
        }
        
        // Create moving window energy domain histograms the range (numpart - 1) deals with right side inclusive...
        snprintf(data->name, sizeof(data->name), "%s_Ch_%d_05_Energy_Domain_Moving_Hist", file, data->channel); // Prevents overwritting of prior files...
        for(int i = 0; i < data->endslice - (numpart - 1); i++){
            hist->Reset();
            for(int j = 0; j < numpart; j++) hist->Add(data->histogram[i + j]);
            for(int j = 0; j < tBin; j++) data->hist[i][j] = 0;
            for(int j = 0; j < numbin; j++) data->hist[i][j] = hist->GetBinContent(j + 1);
        }
        for(int i = 0; i < numbin; i++) data->xaxis[i] = hist->GetBinCenter(i + 1);
        data->cslice = 0;
        window_save(data, numbin, data->endslice - (numpart - 1)); // Save moving window in energy domain...
        /*/
        // Cleanup...
        for (int i = 0; i < data->endslice; i++) data->histogram[i]->Delete();
        delete data;
        hist->Delete();
    } // Trouble Shooting:              Load partition/window and calibration files and run revise_eneregy(...)
    else if(type == 19){
        /*/                                     NOTES (20161012)
         This section is 'identical' to type 16 with the exception that it now loads the revised_energy partition file
         (this may later be adjusted to the moving window average).
         
         output is the fitted peaks in the energy domain. 
         
         In a perfect world the residuals obtained from the fits here would match perfectly with the residuals calculated
         in type 17. Residuals must currently be calculated after the fact as it is currently not scripted...
         
         /*/
        
        // Initialize Function Variables
        int lal = 500;
        int lcl = 2;
        int ucl = 200;
        int counter = 0;
        bool override = true;
        double binstart = 0;
        double binend = 0;
        int pspace = 1;
        //int numpart = 30;  // This is the number of partitions comprising the window generated by the type section above...
        char str[500];
        char word[500];
        char *path = tpath;
        char *ext = text;
        double sigma_value = 1;
        double area_value = 0.001;
        info *data = new info;
        // End
        
        data->current = 0;
        data->last = 0;
        data->endslice = 0;
        data->channel = int(2*randomvalue);
        snprintf(data->path, sizeof(data->path), "%s%s/", path, file);
        snprintf(data->name, sizeof(data->name), "%s", file);
        snprintf(data->ext, sizeof(data->ext), ".txt");
        
        // Load histograms
        //snprintf(data->filename, sizeof(data->filename), "%s%s/%s_Ch_%d_05_Energy_Domain_Moving_Hist.txt", path, file, file, data->channel);
        snprintf(str, sizeof(word), "");
        if(sig < 0) {
            snprintf(word, sizeof(word), "_%d-%d", int(-sig), int(-sig) + 21);
            snprintf(data->filename, sizeof(data->filename), "%s%s/%s_Ch_%d_04_Energy_Domain_Partition_Hist%s.txt", path, file, file, data->channel, word);
        }
        else snprintf(data->filename, sizeof(data->filename), "%s%s/%s_Ch_%d_04_Energy_Domain_Partition_Hist.txt", path, file, file, data->channel);
        file_load(data, 3);
        binstart = (3*data->xaxis[0] - data->xaxis[1])/2;
        binend = data->xaxis[data->last - 1] + (data->xaxis[1] - data->xaxis[0])/2;
        data->numbin = data->last;
        printf("\nNumber of bins: %d\nStart: %.2f\nEnd: %.2f\n", data->last, binstart, binend);
        printf("Number of slices: %d\n", data->endslice);
        for (int j = 0; j < NUM_PEAKS*NUM_COLUMN; j++) data->param[j][0] = -1;
        
        data->histogram[0] = new TH1D("Full", "Full", data->last, binstart, binend); // Calibration invariant
        for (int i = 0; i < data->endslice; i++){
            for (int j = 0; j < data->last; j++) data->histogram[0]->Fill(data->xaxis[j], data->hist[i][j]);
        }
        data->histogram[0]->Sumw2(kFALSE);
        data->histogram[0]->Sumw2();
        data->cslice = 0;
        data->endslice = 1;
        hist2array(data);snprintf(data->name, sizeof(data->name), "%s_Ch_%d_05_Full_Laser_Hist", file, data->channel);
        window_save(data, 0, data->numbin);
        snprintf(data->name, sizeof(data->name), "%s", file);
        for (int j = 0; j < data->endslice; j++) {
            peakfinder(data, &override, &sigma_value, &area_value, true, lcl, ucl, lal, pspace, true);
            data->cslice++;
        }
        printf("\nNumber of bins: %d\nStart: %.2f\nEnd: %.2f\n", data->last, binstart, binend);
        printf("Number of slices: %d\n", data->endslice);
        data->indexvalue = 6;
        fit_save(data, false, 1);
        for (int i = 0; i < data->endslice; i++) data->histogram[i]->Delete();
    } // Trouble Shooting:              Load window file (energy domain) run peakfinder(...)
    else if(type == 20){
        // Initialize Variables
        int counter = 0;
        int numpart = 30;
        char str[500];
        char *path = tpath;
        info *data = new info;
        // End
        
        // General data object setup
        data->current = 0;
        data->last = 0;
        data->endslice = 0;
        data->channel = 2*int(randomvalue);
        sprintf(data->filename, "%s%s/%s_Ch_%d_06_Parameters.txt", path, file, file, data->channel);
        sprintf(data->path, "%s%s/", path, file);
        sprintf(data->name, "%s_Energy_Domain", file);
        sprintf(data->ext, ".txt");
        file_load(data, 5); // Parameter file load
        data->cslice = 0;
        data->indexvalue = 7;
        for (int i = 0; i < data->endslice; i++){
            calibration(data);
            residual(data, true);
            data->cslice++;
        }
        data->indexvalue = 8;
        fit_save(data, false, 2);
        delete data;
    } // Trouble Shooting:              Load Parameter file (moving -> energy) and rerun calibration
    else if(type == 21){
        // Function Variables:
        int counter = 0;
        int eslot = 0;
        int lal = 200;
        int lcl = 200;
        int ucl = 2000;
        int pspace = 10;
        int fbin = 4*2160/(sig <= 0 ? 2 : sig);// The multiplicative factor here matches fupp for 0.2 eV/bin
        int range = 15;
        int backtype = 0;
        int standnum = 0; //50; // This is the partition to use as the standard when calculating the gain matching
        int attempt = 0;
        bool override = false; //(sig < 0 ? false : true);
        bool trouble = false;
        bool specialize = false;
        bool automate = true;
        char str[500];
        char *path = tpath;
        char *ext = text;
        char fname[5][100] = {"Single_Exp", "Double_Exp", "Linear", "Quadratic", "Inverse"};
        float lstart = (0 <= sig ? 66.5 : -sig);
        float center = 76.5;
        double pvalue[7] = {1000, 76.5, 0.9, 20000, 40, 0, 0}; //{1000, 76.8, 0.9, 20000, 40, 0, 0}; //
        double sigma_value = 3;
        double area_value = 0.00005;
        double binstart = 0;
        double binend = 0;
        double flow = 1;
        double fupp = 4*216 + flow; // The multiplicative factor here matches fbin for 0.2 eV/bin
        double gain = 0;
        double dgai = 0;
        double offs = 0;
        double doff = 0;
        double poin = 0;
        double dpoi = 0;
        double centroid[NUM_PEAKS][NUM_COLUMN + 1]; // The additional term is used to keep track of position in filling and total counts
        double delta[NUM_PEAKS][NUM_COLUMN];
        double output[2];
        
        info *data = new info;
        TH1D *laser;
        TH1D *signal;
        TCanvas *canvas;
        // End
        
        if (true){
            data->current = 0;
            data->last = 0;
            data->endslice = 0;
            data->channel = 2*int(randomvalue);
            snprintf(data->filename, sizeof(data->filename), "%s%s/%s_Ch_%d_01_Parameters.txt",
                     path, file, file, data->channel);
            sprintf(data->path, "%s%s/", path, file);
            sprintf(data->name, "%s", file);
            sprintf(data->ext, ".txt");
            file_load(data, 5); // Parameter file load
            data->cslice = 0;
            snprintf(data->filename, sizeof(data->filename), "%s%s/%s_Ch_%d_00_Window_Hist.txt",
                     path, file, file, data->channel);
            file_load(data, 3); // Window file load
            
            binstart = (3*data->xaxis[0] - data->xaxis[1])/2;
            binend = data->xaxis[data->last - 1] + (data->xaxis[1] - data->xaxis[0])/2;
            printf("\nNumber of bins: %d\nStart: %.2f\nEnd: %.2f\n", data->numbin, binstart, binend);
            printf("Number of slices: %d\n", data->endslice);
            for (int i = 0; i < data->endslice; i++){
                snprintf(str, sizeof(str), "Partition_%d_Channel_%d", i, data->channel);
                data->histogram[i] = new TH1D(str, str, data->numbin, binstart, binend); // Calibration invariant
                for (int j = 0; j < data->numbin; j++){
                    data->histogram[i]->Fill(data->xaxis[j], data->hist[i][j]);
                }
                data->histogram[i]->Sumw2(kFALSE);
                data->histogram[i]->Sumw2();
            }// Generate BIN domain histograms
            eslot = data->endslice;
        } // Initial setup
        
        if(true){
            driftcor(data, binstart, binend, standnum, 0, trouble);
            data->indexvalue = 2;
            fit_save(data, false, 2);// Save Drift correction gain values
        } // Drift Correction
        else{
            snprintf(data->filename, sizeof(data->filename), "%s%s/%s_Ch_%d_02_Calibration.txt",
                     path, file, file, data->channel);
            file_load(data, 4); // Window file load
        } // OR Load Prior Calculation

        if (true){
            for(int j = 0; j < data->numbin; j++) data->xaxis[j] = 0;
            for(int j = 0; j < data->numbin; j++) data->xaxis[j] = binstart + (j + 0.5)*(binend - binstart)*1.0/data->numbin;
            data->cslice = 0;
            snprintf(data->name, sizeof(data->name), "%s_Ch_%d_03_Laser_Gain_Matched_Partitions", file, data->channel);
            window_save(data, 0, data->numbin); // save gain matched partitions
            for (int i = 0; i < data->endslice; i++) data->histogram[i]->Delete(); //Clean up
            laser = new TH1D("Laser", "Laser", data->numbin, binstart, binend);
            for (int i = 0; i < data->endslice; i++) for (int j = 0; j < data->last; j++) laser->Fill(data->xaxis[j], data->hist[i][j]);
            data->cslice = 0;
            data->endslice = 1;
            data->histogram[0] = laser;
            hist2array(data);
            snprintf(data->name, sizeof(data->name), "%s_Ch_%d_04_Full_Laser", file, data->channel);
            window_save(data, 0, data->numbin);// Save full laser file histogram
        } // Apply gain matching, save new partitions, sum file, and drift correction terms
        else{
            for (int i = 0; i < data->endslice; i++) data->histogram[i]->Delete(); //Clean up
            laser = new TH1D("Laser", "Laser", data->numbin, binstart, binend);
        }
        
        if (false){
            printf("Loading signal file...\n");
            snprintf(data->filename, sizeof(data->filename), "%s%s/%s_Ch_%d_00_Signal_Hist.txt", path, file, file, data->channel);
            file_load(data, 3); // Signal file load channel domain
            binstart = (3*data->xaxis[0] - data->xaxis[1])/2;
            binend = data->xaxis[data->last - 1] + (data->xaxis[1] - data->xaxis[0])/2;
            printf("\nNumber of bins: %d\nStart: %.2f\nEnd: %.2f\n", data->last, binstart, binend);
            printf("Number of slices: %d\n", data->endslice);
            for (int i = 0; i < data->endslice; i++){
                snprintf(str, sizeof(str), "Partition_%d_Channel_%d", i, data->channel);
                data->histogram[i] = new TH1D(str, str, data->last, binstart, binend); // Calibration invariant
                for (int j = 0; j < data->last; j++) data->histogram[i]->Fill(data->xaxis[j], data->hist[i][j]);
            } // Generate BIN domain histograms
            data->cslice = 0;
            //for (int i = 0; i < data->endslice; i++) data->calib[i][3] = 0; // Clear offset term for laser only (?)
            for (int i = 0; i < data->endslice; i++){
                revise_energy(data, data->numbin, data->numbin, binstart, binend, trouble);
                data->cslice++;
            } // Convert signal to gain corrected basis
            
            for(int j = 0; j < data->numbin; j++) data->xaxis[j] = 0;
            for(int j = 0; j < data->numbin; j++) data->xaxis[j] = binstart + (j + 0.5)*(binend - binstart)*1.0/data->numbin;
            data->cslice = 0;
            snprintf(data->name, sizeof(data->name), "%s_Ch_%d_05_Signal_Gain_Matched_Partitions", file, data->channel);
            window_save(data, 0, data->numbin); // save gain matched partitions
            
            for (int i = 0; i < data->endslice; i++) data->histogram[i]->Delete(); //Clean up
            
            signal = new TH1D("Signal", "Signal", data->numbin, binstart, binend);
            for (int i = 0; i < data->endslice; i++) for (int j = 0; j < data->last; j++) signal->Fill(data->xaxis[j], data->hist[i][j]);
            data->cslice = 0;
            data->endslice = 1;
            data->histogram[0] = signal;
            hist2array(data);
            snprintf(data->name, sizeof(data->name), "%s_Ch_%d_06_Full_Signal", file, data->channel);
            window_save(data, 0, data->numbin);// Save full laser file histogram
            snprintf(data->name, sizeof(data->name),  "%s", file); // Resets base name for other files...
        } // Load signal file apply gain matching, save new partitions and sum file
        else{
            signal = new TH1D("Signal", "Signal", data->numbin, binstart, binend);
        }
        
        if (true){
            data->cslice = 0;
            data->endslice = 3;
            data->histogram[0] = laser;
            data->histogram[1] = signal;
            data->histogram[2] = signal;
            for (int i = 0; i < data->endslice; i++){
                hist2array(data);
                data->cslice++;
            }
            laser->Sumw2(kFALSE);
            laser->Sumw2();
            //laser->Scale(1, "width");
            signal->Sumw2(kFALSE);
            signal->Sumw2();
            //signal->Scale(1, "width");
        } // Set up data object to contain laser, signal, signal data in slots 0, 1 and 2, respectively.
        else{
            data->cslice = 0;
            snprintf(data->filename, sizeof(data->filename), "%s%s/%s_Ch_%d_04_Full_Laser.txt", path, file, file, data->channel);
            file_load(data, 3); // Window file load
            data->cslice = 0;
            array2hist(data, laser);
            snprintf(data->filename, sizeof(data->filename), "%s%s/%s_Ch_%d_06_Full_Signal.txt", path, file, file, data->channel);
            file_load(data, 3); // Window file load
            data->cslice = 0;
            array2hist(data, signal);
            hist2array(data->hist[0], laser);
            hist2array(data->hist[1], signal);
            hist2array(data->hist[2], signal);
            data->histogram[0] = laser;
            data->histogram[1] = signal;
            data->histogram[2] = signal;
            data->endslice = 3;
            laser->Sumw2(kFALSE);
            laser->Sumw2();
            //laser->Scale(1, "width");
            signal->Sumw2(kFALSE);
            signal->Sumw2();
            //signal->Scale(1, "width");
        } // OR Load Prior Data
        
        if (false){
            data->cslice = 0;
            for (int j = 0; j < NUM_PEAKS*NUM_COLUMN; j++) data->param[j][0] = -1; // 'Clear' loaded parameters
            peakfinder(data, &override, &sigma_value, &area_value, false, lcl, ucl, lal, pspace, false, 2, 2, trouble);
            calibration(data, lal, specialize, trouble, lstart);
            for (int i = 0; i < 7; i++) data->calib[1][i] = data->calib[0][i]; // Copy calibration to slot 1 with offset term
            for (int i = 0; i < 7; i++) data->calib[2][i] = (i < 3 ? data->calib[0][i] : 0); // Copy calibration to slot 2 no offset term
            data->indexvalue = 7;
            fit_save(data, false, 0); // Save parameters and calibration
            data->indexvalue = 9;
            data->endslice = 1;
            residual(data); // Save residuals
            data->endslice = 3;
        } // Perform final fit to sum laser, save parameters, perform calibration, save calibration and residuals
        else{
            data->cslice = 0;
            data->endslice = eslot;
            lal = 0;
            sprintf(data->name, "%s", file);
            for (int i = 0; i < NUM_PEAKS; i++) for(int j = 0; j < NUM_COLUMN + 1; j++) centroid[i][j] = 0;
            for(int i = 0; i < data->endslice; i++){
                gain = data->calib[i][1];
                dgai = data->calib[i][2];
                offs = data->calib[i][3];
                doff = data->calib[i][4];
                counter = 0;
                while (counter < NUM_PEAKS){
                    eslot = int(data->param[NUM_PEAKS*i + counter][1]/3.5);
                    poin = data->param[NUM_PEAKS*i + counter][2];
                    dpoi = data->param[NUM_PEAKS*i + counter][3];
                    data->param[NUM_PEAKS*i + counter][2] = gain*poin + offs;
                    data->param[NUM_PEAKS*i + counter][3] = pow(dgai*dgai*poin*poin + gain*gain*dpoi*dpoi + doff*doff, 0.5);
                    data->param[NUM_PEAKS*i + counter][4] = 0;
                    data->param[NUM_PEAKS*i + counter][5] = 0;
                    data->param[NUM_PEAKS*i + counter][6] = 0;
                    data->param[NUM_PEAKS*i + counter][7] = 0;
                    if (counter != 0 && data->param[NUM_PEAKS*i + counter + 1][0] != -1){
                        centroid[eslot][int(centroid[eslot][NUM_COLUMN])] = data->param[NUM_PEAKS*i + counter][2];
                        delta[eslot][int(centroid[eslot][NUM_COLUMN])] = 1.0/data->param[NUM_PEAKS*i + counter][3];
                        centroid[eslot][NUM_COLUMN]++;
                    }
                    counter++;
                    if (data->param[NUM_PEAKS*i + counter][0] == -1) break;
                }
            }// Transform the original centroids to the gain matched values.
            data->indexvalue = 11;
            fit_save(data, false, 1);
            for(int i = 0; i < NUM_PEAKS*NUM_COLUMN; i++) data->param[i][0] = -1;
            printf("Starting Averages Caluculations\n");
            counter = 0;
            for(int i = 0; i < NUM_PEAKS; i++) if (10 <= centroid[i][NUM_COLUMN]){
                printf("On Peak energy %0.2f with %d peaks\n", 3.5*i, int(centroid[i][NUM_COLUMN]));
                stat(output, centroid[i], delta[i], centroid[i][NUM_COLUMN]);
                data->param[counter][0] = 0;
                data->param[counter][1] = 3.5*i;
                data->param[counter][2] = output[0];
                data->param[counter][3] = output[1]*pow(centroid[i][NUM_COLUMN], -0.5); // Measurement Uncertainty
                data->param[counter][4] = output[1]; // Standard Deviation
                counter++;
            }
            data->cslice = 0;
            data->endslice = 3;
            calibration(data, lal, specialize, trouble, lstart);
            for (int i = 0; i < 7; i++) data->calib[1][i] = data->calib[0][i]; // Copy calibration to slot 1 with offset term
            for (int i = 0; i < 7; i++) data->calib[2][i] = (i < 3 ? data->calib[0][i] : 0); // Copy calibration to slot 2 no offset term
            data->indexvalue = 12;
            fit_save(data, false, 0); // Save parameters and calibration
            data->indexvalue = 14;
            data->endslice = 1;
            residual(data); // Save residuals
            data->endslice = 3;
        } // Perform Calibration using weighted average values of the gain matched parameters
        
        if (true){
            data->cslice = 0;
            for (int i = 0; i < data->endslice; i++){
                data->cslice = 0;
                snprintf(data->name, sizeof(data->name), "%s_Ch_%d_2%d_Full_Laser", file, data->channel, i);
                window_save(data, 0, data->numbin);// Save full laser file histogram
                
                data->cslice = i;
                revise_energy(data, data->numbin, fbin, flow, fupp, trouble);
                data->cslice++;
            }
            
            
            
            
            data->cslice = 0;
            snprintf(data->name, sizeof(data->name), "%s_Ch_%d_10_Calibrated_Data", file, data->channel);
            window_save(data, 0, fbin); // save gain matched partitions
            snprintf(data->name, sizeof(data->name),  "%s", file); // Resets base name for other files...
            signal->Delete(); // Clean up
            laser->Delete(); // Clean up
        } // Perform calibration transformation on laser and signal then save
        
        if (false){
            printf("Performing Decay Fit\n");
            data->cslice = 0;
            data->endslice = 2;
            data->histogram[0] = new TH1D("With Offset", "With Offset", fbin, flow, fupp);
            data->histogram[1] = new TH1D("Without Offset", "Without Offset", fbin, flow, fupp);
            snprintf(data->path, sizeof(data->path),  "%s", path); // Resets base name for other files...
            for (int i = 0; i < data->endslice; i++){
                signal = data->histogram[i];
                data->cslice = i + 1;
                array2hist(data, signal);
                data->cslice = i;
                hist2array(data);
                signal->Sumw2(kFALSE);
                signal->Sumw2();
                //signal->Scale(1, "width");
                automate = background(data, signal, pvalue, backtype, center - range, center + range, data->last, trouble);
                while (trouble || !automate) {
                    if (0 <= attempt && attempt < 7) background(data, signal, pvalue, backtype, center - range, center + range, data->last, trouble);
                    printf("Current background fit is %s\n", fname[backtype]);
                    printf("To Change a parameter enter the paramter line value (-1 to exit): ");
                    while (true) {
                        std::cin >> attempt;
                        if (cin.fail()) {
                            std::cin.clear();
                            std::cin.ignore();
                            std::cout << "Invalid input please enter one of the available options...\n";
                        }
                        else break;
                    }
                    if (attempt == -1) {
                        automate = true;
                        break;
                    }
                    else if (0 <= attempt && attempt < 7) while (true) {
                        printf("\n\t\t\tCurrent Parameters\n\nFunction: %s\n", fname[backtype]);
                        printf("Binning: %.2f\nRange: +/- %d\n", (binend - binstart)*1.0/data->last, range);
                        printf("Set values\n");
                        for (int k = 0; k < 7; k++) printf("p%d\t%f\n", k, pvalue[k]);
                        printf("Please enter new value: ");
                        std::cin >> pvalue[attempt];
                        printf("\n\n");
                        if (cin.fail()) {
                            std::cin.clear();
                            std::cin.ignore();
                            std::cout << "Invalid input please enter one of the available options...\n";
                        }
                        else break;
                    }
                    else{
                        printf("Prior centroid was: %f\n", data->fit[data->cslice - 1][5]);
                        printf("Current centroid is : %f\n", data->fit[data->cslice][5]);
                        canvas = new TCanvas(data->name, data->name, 1000, 1000);
                        signal->GetXaxis()->SetRangeUser(center - range, center + range);
                        signal->Draw();
                        canvas->Update();
                        canvas->Update();
                    }
                }
            }
            for (int i = 0; i < data->endslice; i++) data->histogram[i]->Delete();
        } // Perform signal fit before exiting...TEMPORARY Suspension of lines...
    
        if (false){
            printf("Performing Decay Fit in uncalibrated space...\n");
            data->cslice = 0;
            data->endslice = 1;
            data->histogram[0] = signal;
            hist2array(data);
            snprintf(data->path, sizeof(data->path),  "%s", path); // Resets base name for other files...
            background(data, signal, pvalue, backtype, center - range, center + range, data->last, trouble);
            for (int i = 0; i < data->endslice; i++) data->histogram[i]->Delete();
        } // Perform signal fit before exiting...

        
        printf("Performing cleanup\n\n\n");
        delete data;
    } // Trouble Shooting:              Full XIA drift correction run after initial peak finding.
    else if(type == 22){
        /*/                                     NOTES (20160921)
         This section is a trouble shooting piece demonstrating how to load the window/parameter file. The histograms are loaded into
         the data->hist[...] array and the corresponding histograms are generated and populate data->histogram[i][...].
        /*/
        // Function Variables:
        char str[500];
        char *path = tpath;
        char *ext = text;
        int counter = 0;
        double binstart = 0;
        double binend = 0;
        info *data = new info;
        // End
        
        // General data file setup
        data->current = 0;
        data->last = 0;
        data->endslice = 0;
        data->channel = int(2*randomvalue);
        snprintf(data->path, sizeof(str), "%s%s/", path, file);
        snprintf(data->name, sizeof(str), "%s", file);
        snprintf(data->ext, sizeof(str), ".txt");
        
        // Load partition/window file
        snprintf(data->filename, sizeof(data->filename), "%s%s/%s_Ch_%d_00_Window_Hist.txt", path, file, file, data->channel);
        file_load(data, 3);
        
        // Generate histograms
        binstart = (3*data->xaxis[0] - data->xaxis[1])/2;
        binend = data->xaxis[data->last - 1] + (data->xaxis[1] - data->xaxis[0])/2;
        for (int j = 0; j < NUM_PEAKS*NUM_COLUMN; j++) data->param[j][0] = -1;
        for (int i = 0; i < data->endslice; i++){
            snprintf(str, sizeof(str), "Partition_%d_Channel_%d", i, data->channel);
            data->histogram[i] = new TH1D(str, str, data->last, binstart, binend);
            for (int j = 0; j < data->last; j++) data->histogram[i]->Fill(data->xaxis[j], data->hist[i][j]);
        }
        
        // Load parameter file
        snprintf(data->filename, sizeof(data->filename), "%s%s/%s_Ch_%d_01_Parameters.txt", path, file, file, data->channel);
        file_load(data, 5);
        
        // Load calibration file
        snprintf(data->filename, sizeof(data->filename), "%s%s/%s_Ch_%d_03_Calibration.txt", path, file, file, data->channel);
        file_load(data, 4);
        
        
        // Clean up -> Delete histograms and Data file
        for (int i = 0; i < data->endslice; i++) data->histogram[i]->Delete();
        delete data;
    
    } // Trouble Shooting:              Example on loading a partition/window, parameter, and calibration file...
    else if(type == 23){
        // Initialize Function Variables
        int lal = 20;
        //int lcl = 5;
        int lcl = 6;
        int ucl = 500;
        int counter = 0;
        double binstart = 0;
        double binend = 0;
        int pspace = 2;
        char str[500];
        char word[500];
        char *path = tpath;
        char *ext = text;
        double sigma_value = 1;
        double area_value = 0.0005;
        info *data = new info;
        int numbin = 0;
        bool override = false;
        bool nonlineardet = false;
        double xaxis[tBin + 1];
        double noncent = 75;
        double nonline = 0.25/(50*50);//(noncent*noncent); Changed on 170123 to make the range 50 eV regardless of center
        TH1D *hist;
        troll = data;
        // End
        printf("Starting Full Laser fitting\n\n\n");
        printf("Detector is %s\n\n\n", nonlineardet ? "Non-Linear" : "Linear");
        data->current = 0;
        data->last = 0;
        data->endslice = 0;
        data->channel = int(2*randomvalue);
        snprintf(data->path, sizeof(data->path), "%s%s/", path, file);
        snprintf(data->name, sizeof(data->name), "%s", file);
        snprintf(data->ext, sizeof(data->ext), ".txt");
        
        // Load to buffer array and fill histograms
        snprintf(str, sizeof(word), "");
        if(sig < 0) snprintf(word, sizeof(word), "_%d-%d", int(-sig), int(-sig) + 21);
        
        snprintf(data->filename, sizeof(data->filename), "%s%s/%s_Ch_%d_04_Energy_Domain_Partition_Hist%s.txt", path, file, file, data->channel, word);
        file_load(data, 3);
        binstart = (3*data->xaxis[0] - data->xaxis[1])/2;
        binend = data->xaxis[data->last - 1] + (data->xaxis[1] - data->xaxis[0])/2;
        numbin = data->last;
        
        printf("\nNumber of bins: %d\nStart: %.2f\nEnd: %.2f\n", numbin, binstart, binend);
        printf("Number of slices: %d\n", data->endslice);
        printf("X-axis Values\n x_0: %.2f\tx_1: %.2f\tx_f: %.2f\n", data->xaxis[0], data->xaxis[1], data->xaxis[data->last - 1]);
        for (int j = 0; j < NUM_PEAKS*NUM_COLUMN; j++) data->param[j][0] = -1;
        for (int i = 0; i < data->endslice; i++){
            snprintf(str, sizeof(str), "Partition_%d_Channel_%d", i, data->channel);
            data->histogram[i] = new TH1D(str, str, numbin, binstart, binend); // Calibration invariant
            if (nonlineardet){
                data->histogram[i]->Delete();
                numbin = 8192;
                data->numbin = 8192;
                for (int i = 0; i <= tBin + 1; i++) xaxis[i] = 0.1*i + nonline*pow(0.1*i - noncent, 2);
                data->histogram[i] = new TH1D(str, str, tBin, xaxis);
            }// For Variable bin widths...
            for (int j = 0; j < numbin; j++) data->histogram[i]->Fill(data->xaxis[j], data->hist[i][j]);
            data->histogram[i]->Sumw2(kFALSE);
            data->histogram[i]->Sumw2();
        }
        snprintf(str, sizeof(str), "Partition_Sum_Channel_%d", data->channel);
        hist = new TH1D(str, str, numbin, binstart, binend);
        if (nonlineardet){
            hist->Delete();
            hist = new TH1D(str, str, tBin, xaxis);
        }// For Variable bin widths...
        for (int i = 0; i < data->endslice; i++) hist->Add(data->histogram[i]);
        hist->Sumw2(kFALSE);
        hist->Sumw2();
        if (nonlineardet) hist->Scale(1, "width");// For Variable bin widths...
        
        
        data->histogram[0] = hist;
        data->cslice = 0;
        peakfinder(data, &override, &sigma_value, &area_value, true, lcl, ucl, lal, pspace, true);
        
        data->indexvalue = sig; // insure different name...
        fit_save(data, false, 1);
        
        // Clean up
        for (int i = 1; i < data->endslice; i++) data->histogram[i]->Delete();
        
        data->endslice = 1;
        data->cslice = 0;
        sprintf(data->name, "%s_Ch_%d_020_Energy_Full_Laser_Hist", data->name, data->channel);
        for(int i = 0; i < numbin; i++) data->hist[0][i] = data->histogram[0]->GetBinContent(i + 1);
        window_save(data, 0, numbin);
        // Clean up
        data->histogram[0]->Delete();
        //troll = data;
        //hist->Delete();
    } // Trouble Shooting:              Load Energy partition combine into single histogram perform full laser fit
    else if(type == 24){
        // Initialize Function Variables
        int lal = 20;
        int lcl = 12;
        int ucl = 80;
        int counter = 0;
        double binstart = 0;
        double binend = 0;
        int pspace = 1;
        int numpart = 30;  // This is the number of partitions comprising the window generated by the type section above...
        char str[500];
        char *path = tpath;
        char *ext = text;
        bool override = true;
        double sigma_value = 1;
        double area_value = 0.001;
        info *data = new info;
        int numbin = 1024;
        // End
        
        data->current = 0;
        data->last = 0;
        data->endslice = 0;
        data->channel = int(2*randomvalue);
        snprintf(data->path, sizeof(data->path), "%s%s/", path, file);
        snprintf(data->name, sizeof(data->name), "%s", file);
        snprintf(data->ext, sizeof(data->ext), ".txt");
        
        // Load to buffer array and fill histograms
        //snprintf(data->filename, sizeof(data->filename), "%s%s/%s_Ch_%d_05_Energy_Domain_Moving_Hist.txt", path, file, file, data->channel);
        snprintf(data->filename, sizeof(data->filename), "%s%s/%s_Ch_%d_00_Window_Hist.txt", path, file, file, data->channel);
        file_load(data, 3);
        snprintf(data->filename, sizeof(data->filename), "%s%s/%s_Ch_%d_03_Calibration.txt", path, file, file, data->channel);
        file_load(data, 4); // Calibration file load
        
        
        
        
        
        for (int j = 0; j < NUM_PEAKS*NUM_COLUMN; j++) data->param[j][0] = -1;
        for (int i = 0; i < data->endslice; i++){
            binstart = data->calib[i][3];
            binend = data->calib[i][1]*numbin + data->calib[i][3];
            printf("\nNumber of bins: %d\nStart: %.2f\nEnd: %.2f\n", numbin, binstart, binend);
            printf("Number of slices: %d\n", data->endslice);
            snprintf(str, sizeof(str), "Partition_%d_Channel_%d", i, data->channel);
            data->histogram[i] = new TH1D(str, str, numbin, binstart, binend); // Calibration invariant
            
            for (int j = 0; j < numbin; j++){
                data->histogram[i]->Fill(data->xaxis[j]*data->calib[i][1] + data->calib[i][3], data->hist[i][j]);
            }
        }
        //chost = data->histogram[0];
        //*/
        for (int j = 0; j < data->endslice; j++) {
            peakfinder(data, &override, &sigma_value, &area_value, true, lcl, ucl, lal, pspace, true);
            data->cslice++;
        }
        data->indexvalue = 20; // insure different name...
        fit_save(data, false, 1);
        for (int i = 0; i < data->endslice; i++) data->histogram[i]->Delete();
        //*/
    } // Trouble Shooting:              Load window file (energy domain) and calibration files (independent histogram x-axis) run peakfinder(...)  => Perform fits of reconstructed energy
    else if(type == 25){
        /*/
            20161213
            Updated the numbin value to depend on sig with default value of 0.2eV
        /*/
        // Function Variables:
        int counter = 0;
        bool sansoff = true; // True means do not use offset calculated in calibration for calculating the revised energy
        double binstart = 0;
        double binend = 0;
        int numbin = 2*2160/(sig <= 0 ? 2 : sig);  //Use 1, 2, 3, 4, 5 for 0.1, 0.2, 0.3, 0.4, 0.5 eV/bin
        //int numpart = 30; // Moving window average is this value at 2 minutes each...
        char str[500];
        char word[500];
        char *path = tpath;
        char *ext = text;
        double lower = 0;
        double upper = 2*216 + lower;
        info *data = new info;
        TH1D *hist = new TH1D("Main", "Main", numbin, lower, upper);
        // End
        //printf("sig is %d, numbin is %d, numbin float is %d, difference is %d\n", sig, numbin, 2160/(sig <= 0 ? 2 : sig), 2160/(sig <= 0 ? 2 : sig) - numbin);
        printf("Converting to Number of bins: %d\nStart: %.2f\nEnd: %.2f\n", numbin, lower, upper);
        
        data->current = 0;
        data->last = 0;
        data->endslice = 0;
        data->channel = int(2*randomvalue);
        snprintf(data->filename, sizeof(data->filename), "%s%s/%s_Ch_%d_00_Signal_Hist.txt", path, file, file, data->channel); // Window file name...
        snprintf(data->path, sizeof(data->path), "%s%s/", path, file);
        snprintf(data->name, sizeof(data->name), "%s", file);
        snprintf(data->ext, sizeof(data->ext), ".txt");
        
        // Load to buffer array and fill histograms
        file_load(data, 3); // Signal file load channel domain
        binstart = (3*data->xaxis[0] - data->xaxis[1])/2;
        binend = data->xaxis[data->last - 1] + (data->xaxis[1] - data->xaxis[0])/2;
        printf("\nNumber of bins: %d\nStart: %.2f\nEnd: %.2f\n", data->last, binstart, binend);
        printf("Number of slices: %d\n", data->endslice);
        
        for (int j = 0; j < NUM_PEAKS*NUM_COLUMN; j++) data->param[j][0] = -1;
        for (int i = 0; i < data->endslice; i++){
            snprintf(str, sizeof(str), "Partition_%d_Channel_%d", i, data->channel);
            data->histogram[i] = new TH1D(str, str, data->last, binstart, binend); // Calibration invariant
            for (int j = 0; j < data->last; j++){
                data->histogram[i]->Fill(data->xaxis[j], data->hist[i][j]);
            }
        } // Generate BIN domain histograms
        
        snprintf(str, sizeof(word), "");
        if(sig < 0) snprintf(word, sizeof(word), "_%d-%d", int(-sig), int(-sig) + 21);
        
        snprintf(data->filename, sizeof(data->filename), "%s%s/%s_Ch_%d_03_Calibration%s.txt", path, file, file, data->channel, word);
        file_load(data, 4); // Calibration file load
        
        // Clear calibration offset
        if (sansoff) for (int i = 0; i < data->last; i++) data->calib[i][3] = 0;
        
        
        // Run revised energy to transform partition basis into energy basis...
        printf("Start converting to Number of bins: %d\n", numbin);
        data->cslice = 0;
        for (int j = 0; j < data->endslice; j++) {
            revise_energy(data, tBin, numbin, lower, upper);
            data->cslice++;
        }
         printf("Start converting to Number of bins: %d\n", numbin);
        // Over write x-axis and create single histogram of all the data...
        for(int j = 0; j < tBin; j++) data->xaxis[j] = 0;
        for(int j = 0; j < numbin; j++) data->xaxis[j] = lower + (j + 0.5)*(upper - lower)*1.0/numbin;
        for (int i = 0; i < numbin; i++) for (int j = 0; j < data->endslice; j++) hist->Fill(data->xaxis[i], data->hist[j][i]);
        snprintf(data->name, sizeof(data->name), "%s_Ch_%d_20_Energy_Partition_Signal_Hist", file, data->channel);
        data->cslice = 0;
        printf("Number of columns In: %d\n", numbin);
        window_save(data, 0, numbin);
        
        for (int i = 0; i < numbin; i++) data->hist[0][i] = hist->GetBinContent(i + 1);
        
        // Save energy domain partitions and clean up
        
        snprintf(str, sizeof(word), "");
        if(sig < 0) snprintf(word, sizeof(word), "_%d-%d", int(-sig), int(-sig) + 21);
        
        snprintf(data->name, sizeof(data->name), "%s_Ch_%d_20_Energy_Full_Signal_Hist%s", file, data->channel, word); // Prevents overwritting of prior files...
        data->cslice = 0;
        for (int i = 1; i < data->endslice; i++) data->histogram[i]->Delete();
        data->endslice = 1;
        printf("Number of columns In: %d\n", numbin);
        window_save(data, 0, numbin); // Save partition in energy domain...
        data->histogram[0]->Delete();
        delete data;
        hist->Delete();
    } // Trouble Shooting:              Load signal file and calibration files and run revised_energy(...)
    else if(type == 26){
        // Initialize Variables
        int exit = 0;
        int range = 25; // Should be 25
        int attempt = 0;
        int ftype = 0; // Matches order for fname
        bool trouble = (randomvalue <= 0 ? false : true);
        bool automate = false;
        char str[500];
        char word[500];
        char *path = tpath;
        char fname[5][100] = {"Single_Exp", "Double_Exp", "Linear", "Quadratic", "Inverse"};
        float center = 76.8;
        double pvalue[7] = {50, 76.8, 0.7, 20000, 40, 0, 0};
        double val = 0;
        double binstart = 0;
        double binend = 0;
        TH1D *hist;
        TCanvas *canvas;
        info *data = new info;
        // End
        
        // General info object setup...
        randomvalue = TMath::Abs(randomvalue);
        data->current = 0;
        data->last = 0;
        sprintf(data->path, "%s", path);
        sprintf(data->name, "%s", file);
        sprintf(data->ext, "%s", text);
        data->channel = int(2*randomvalue);
        // Load signal data in revise_energy(...) space and reconstruct histogram...
        
        snprintf(str, sizeof(word), "");
        if(sig < 0) snprintf(word, sizeof(word), "_%d-%d", int(-sig), int(-sig) + 21);
        
        if (sig < 0) sprintf(data->filename, "%s%s/%s_Ch_%d_20_Energy_Full_Signal_Hist%s.txt", path, file, file, data->channel, word);
        else sprintf(data->filename, "%s%s/%s_Ch_%d_20_Energy_Full_Signal_Hist.txt", path, file, file, data->channel);
        if (0 < sig) sprintf(data->filename, "%s%s/%s_Ch_%d_10_Calibrated_Data.txt", path, file, file, data->channel);
        
        printf("Opening file: %s", data->filename);
        file_load(data, 3);
        //data->endslice = 1;// Toggle to process enmass...
        binstart = (3*data->xaxis[0] - data->xaxis[1])/2;
        binend = data->xaxis[data->last - 1] + (data->xaxis[1] - data->xaxis[0])/2;
        printf("\nNumber of bins: %d\nStart: %.2f\nEnd: %.2f\n", data->last, binstart, binend);
        printf("First bin: %f\nLast: %.2f\n", data->xaxis[0], data->xaxis[data->last - 1]);
        printf("Number of slices: %d\n", data->endslice);
        for (int i = 0; i < data->endslice; i++){
            snprintf(str, sizeof(str), "Full_Signal_Channel_%d_%d", data->channel, i);
            data->histogram[i] = new TH1D(str, str, data->last, binstart, binend);
            for (int j = 0; j < data->last; j++) data->histogram[i]->Fill(data->xaxis[j], data->hist[i][j]);
            data->histogram[i]->Sumw2(kFALSE);
            data->histogram[i]->Sumw2();
            //data->histogram[i]->Scale(1, "width");
        }
        hist = data->histogram[int(0 <= sig ? sig : 0)];
        // Perform background fit to data...
        snprintf(str, sizeof(str), "%s%s%d", data->name, "_Channel_", data->channel);
        printf("Run initial background\n");
        data->cslice = 0;
        for (int iter = (0 <= sig ? sig : 0) ; iter < data->endslice; iter++){
            //if (iter == 0) hist = data->histogram[iter];
            //else hist->Add(data->histogram[iter]);
            hist = data->histogram[iter];
            for (int i = 0; i < 5; i++){
                if (false) {
                    for (int k = 3; k < 7; k++) pvalue[k] = 0;
                    pvalue[0] = 150;
                    pvalue[1] = 76.8;
                    pvalue[2] = 0.9;
                }
                for (int j = 0; j < 21; j++){
                    automate = true;
                    if (!trouble) {
                        range = 10 + 2*j;
                        ftype = i;
                    }
                    data->indexvalue = TMath::Abs(sig);
                    automate = background(data, hist, pvalue, ftype, center - range, center + range, data->last, trouble);
                    automate = background(data, hist, pvalue, ftype, center - range, center + range, data->last, trouble);
                    //automate = false;
                    if (trouble) printf("Starting loop\n");
                    while (!automate) {
                        if (0 <= attempt && attempt < 7) background(data, hist, pvalue, ftype, center - range, center + range, data->last, trouble);
                        printf("Current background fit is %s\n", fname[i]);
                        printf("To Change a parameter enter the paramter line value (-1 to exit): ");
                        while (true) {
                            std::cin >> attempt;
                            if (cin.fail()) {
                                std::cin.clear();
                                std::cin.ignore();
                                std::cout << "Invalid input please enter one of the available options...\n";
                            }
                            else break;
                        }
                        if (attempt == -1) {
                            automate = true;
                            break;
                        }
                        else if (0 <= attempt && attempt < 7) while (true) {
                            printf("\n\t\t\tCurrent Parameters\n\nFunction: %s\n", fname[i]);
                            printf("Binning: %.2f\nRange: +/- %d\n", (binend - binstart)*1.0/data->last, range);
                            printf("Set values\n");
                            for (int k = 0; k < 7; k++) printf("p%d\t%f\n", k, pvalue[k]);
                            printf("Please enter new value: ");
                            std::cin >> pvalue[attempt];
                            printf("\n\n");
                            if (cin.fail()) {
                                std::cin.clear();
                                std::cin.ignore();
                                std::cout << "Invalid input please enter one of the available options...\n";
                            }
                            else break;
                        }
                        else{
                            printf("Prior centroid was: %f\n", data->fit[data->cslice - 1][5]);
                            printf("Current centroid is : %f\n", data->fit[data->cslice][5]);
                            canvas = new TCanvas(data->name, data->name, 1000, 1000);
                            hist->GetXaxis()->SetRangeUser(center - range, center + range);
                            hist->Draw();
                            canvas->Update();
                            canvas->Update();
                        }
                    }
                    if (automate && 0 < j) for (int k = 0; k < 7; k++) pvalue[k] = data->fit[data->cslice][2*k + 3];
                    if (trouble) i = j = 9999;
                    else data->cslice++;
                }
            }
            data->cslice++; // Toggle to process enmass...
        }// Toggle to process enmass...
        if (trouble) thist[0] = hist;
        else hist->Delete();
    } // Trouble Shooting:              Load revised energy signal file and run fit with background
    else if(type == 27){
        
        // Initialize Function Variables
        int counter = 0;
        int lal = 0;
        char str[500];
        char *path = tpath;
        char *ext = text;
        bool trouble = false;
        bool specialize = true;
        double lstart = 66.5;
        double correction = 0;
        double dcorr = 0;
        double *out;
        const double *outerr;
        
        info *data = new info;
        TGraphErrors *linefit;
        TF1 *fun;
        // End
        
        data->current = 0;
        data->last = 0;
        data->endslice = 0;
        data->channel = int(2*randomvalue);
        snprintf(data->path, sizeof(str), "%s%s/", path, file);
        snprintf(data->name, sizeof(str), "%s", file);
        snprintf(data->ext, sizeof(str), "%s", ext);
        
        // Load to gain matched info
        snprintf(data->filename, sizeof(str), "%s%s/%s_Ch_%d_12_Parameters.txt", path, file, file, data->channel);
        file_load(data, 5);
        snprintf(data->filename, sizeof(str), "%s%s/%s_Ch_%d_02_Calibration.txt", path, file, file, data->channel);
        file_load(data, 4);
        
        fun = new TF1("fun2", "pol1");
        linefit = new TGraphErrors(data->endslice);
        
        for (int i = 0; i < data->endslice; i++) {
            linefit->SetPoint(linefit->GetN(), data->calib[i][1], data->calib[i][3]);
            linefit->SetPointError(linefit->GetN() - 1, data->calib[i][2], data->calib[i][4]);
        }
        fun->SetParameters(5, -5);
        linefit->Fit(fun, "Q");
        out = fun->GetParameters();
        outerr = fun->GetParErrors();
        
        if(true){
            data->cslice = 0;
            data->endslice = 0;
            snprintf(data->filename, sizeof(str), "%s%s/%s_Ch_%d_13_Calibration.txt", path, file, file, data->channel);
            file_load(data, 4);
            
            data->cslice = 3; //data->endslice;
            data->endslice = 4; //++;
            data->calib[data->cslice][0] = data->cslice;
            data->calib[data->cslice][1] = out[1];
            data->calib[data->cslice][2] = outerr[1];
            data->calib[data->cslice][3] = out[0];
            data->calib[data->cslice][4] = outerr[0];
            //correction = (TMath::Abs(out[1]) + TMath::Abs(out[0]))/2;
            //dcorr = pow(outerr[1]*outerr[1] + outerr[0]*outerr[0], 0.5)/2;
            correction = out[0];
            dcorr = outerr[0];
            data->calib[data->cslice - 1][3] = data->calib[0][3] - data->calib[0][1]*correction;
            data->calib[data->cslice - 1][4] = pow(data->calib[0][4]*data->calib[0][4] +
                                                   data->calib[0][2]*data->calib[0][2]*correction*correction +
                                                   dcorr*dcorr*data->calib[0][1]*data->calib[0][1], 0.5);
            data->indexvalue = 13;
            fit_save(data, false, 2);// Save Drift correction gain values
        }
        else{
            for (int i = 0; i < NUM_PEAKS; i++){
                printf("%0.4f\n", data->param[i][2]);
                data->param[i][2] = data->param[i][2] - TMath::Abs(out[1]);
                data->param[i][3] = pow(data->param[i][3]*data->param[i][3] + outerr[1]*outerr[1], 0.5);
                if (data->param[i][0] == -1) break;
            }
            data->indexvalue = 15;
            fit_save(data, false, 1);
            data->cslice = 0;
            data->endslice = 1;
            calibration(data, lal, specialize, trouble, lstart);
            data->indexvalue = 16;
            fit_save(data, false, 2);
        }
        
        
        delete data;
    } // Trouble Shooting:              Run offset correction from gain value
    else if(type == 28){
        /*/                                     NOTES (20161130)
         This section will be used to generate Monte Carlo simulation of the STJ laser response. I assume
         a poisson distribution to the peak amplitudes and a sigma that is scaling with the sqrt of energy
         for the peak.
         
         20161208
         I have run three of these prior to today, the first was a linear model with the sigma varing as
         a function of energy using the statistical uncertainty eFE added in quadrature to a theoretical 
         electronic noise of 0.8. The mu was kept constant at 20 photons, the bin sizes where 0.1 eV and 
         the range was [0, 220) eV. A partition was given 30k counts which at 100Hz represents 5 minutes 
         of laser acquisition. A total of 50 peaks were used in the model function and a total of 100 
         simulations were run, this represents making the exact same measurement 100 time, the idea is 
         to get an understanding of how well my code will reconstruct the true values of the data each 
         time.
         
         The second MCS introduced a nonlinearity in the centroids centered at 49eV and decreasing by
         0.05 eV at 0 and 100 eV:
                                        -0.05*(x - 49)^2/49^2
         All other parameters were identical to before at this point I believe I also changed over to
         8192 bins so the range is from [0, 819.2) eV
         
         The third MCS introduced a nonlinearity through the bin size of the histogram. The bin edges of
         the simulation the bins edges were assigned using 
                                    bin = E - 0.05*(x - 50)/50^2
         at this point I believe I also changed over to 8192 bins so the range is from [-0.05, 807..)eV
         
         And now for what I am doing now...
         I have adjusted the sigma to take into account the multiple tunneling F = 0.2 -> F = 1.2 and
         adjusted the input areas so that the calulated amplitudes will correspond to the poisson 
         distribution amplitude values. 
                                P(mu, N) -> P(mu, N)*sqrt(2*pi*sigma^2) = Area
                        0.64 + 8.3*10^-4 -> 0.64 + 5*10^-3*n            = sigma^2
         still want to use 30k counts and introduce a moving mu value now to compare to actual data...
         For this purpose I will use the mu value as calculated from a real data set as the area 
         weighted position of peak position for file p242a, generated on python...
         
         Generated a constant mu = 20 MCS and variable mu MCS (based on p242a) at 30k counts, linear, 
         over an energy range of [0, 819.2) with a bin width or 0.1eV and no offset.
         Folders:   Cons_mu_Linear_2
                    Var_mu_Linear
         
         Additionally generated a variable mu MCS (based on p242a) at 30k counts, quadratic, centered
         at 73.5 eV and a +0.05 off set at 113.5 and 33.5 eV.
         Folder:    Var_mu_Nonlinear
         
         20170525
         sig now selects between 1: signal simulation and <= 0: laser with the negative terms returning
         the eV offset of the laser substrate interaction at 77 eV (i.e. -1 => 1 eV offset of laser energy)
         
        /*/
        // Function Variables
        int numcounts = 30000;
        int numpeaks = 50;
        int simtype = (sig <= 0 ? sig : 1); // See notes
        int counter = 0;
        char str[500];
        char *path = mcpath;
        char *ext = text;
        bool trouble = false;
        bool nonlinear = false;// (Obsolete) set to false until removed from code
        bool varmu = false;
        double mu = 20;
        double lower = 0; // Changed this on 170117 was 0
        double upper = 819.2; // Changed this on 170117 was 819.2
        double sigma = 0.8;
        double xaxis[tBin + 1];
        double noncent = 73.5;
        double nonline = 0.05/(noncent*noncent);
        double vmu[tBin];
        info *data = new info;
        
        ifstream file;
        string trash;
        // End

        data->channel = 4;
        data->current = 0;
        data->cslice = 0;
        data->endslice = 100;
        data->numbin = tBin;
        snprintf(data->path, sizeof(str), "%s", mcpath);
        snprintf(data->name, sizeof(str), "Monte_Carlo_Simulation_Ch_%d_00_Window_Hist", data->channel);
        snprintf(data->ext, sizeof(str), "%s", ext);
        
        if(0 < simtype){
            numpeaks = 16600;
            mu = 76.8;
            sigma = 0.8;
            numcounts = data->endslice;
        } // Set signal parameters for modeling
        if(varmu){
            printf("Using time dependent mu for MCS generation\n");
            mu = 0;
            snprintf(str, sizeof(str), "%sMonte_Carlo_Simulation/Var_mu.txt", mcpath);
            file.open(str);
            if(file.fail()){
                printf("Failed to open file:\n%s\n", str);
            }
            file >> trash;
            file >> trash;
            while(file >> trash >> vmu[counter]) counter++;
            data->endslice = counter;
            printf("First 10 mu values...\n");
            for (int i = 0; i < 10; i ++) printf("%d\t%f\n", i, vmu[i]);
        } // Upload time dependent mu values
        else{
            printf("Using constant mu for MCS generation\n");
            for (int i = 0; i < data->endslice; i++) vmu[i] = mu;
        } // Set constant mu values
        if (nonlinear) {
            for (int i = 0; i < int(sizeof(xaxis)/sizeof(xaxis[0])); i++){
                xaxis[i] = (upper - lower)*i/data->numbin + pow((upper - lower)*i/data->numbin - noncent, 2)*nonline;
            }
            for(int i = 0; i < data->endslice; i++){
                snprintf(str, sizeof(str), "Window_%d", i);
                data->histogram[i] = new TH1D(str, str, data->numbin, xaxis);
                monte_carlo(data, numcounts, numpeaks, vmu[i], sigma, lower, upper, data->numbin, 0, trouble);
                hist2array(data);
                data->cslice++;
            }
        }// (Obsolete, this step is now done by response function in process 30) This simulates a Non-linear detector
        else {
            printf("Running simtype %d\n", simtype);
            for(int i = 0; i < data->endslice; i++){
                snprintf(str, sizeof(str), "Window_%d", i);
                data->histogram[i] = new TH1D(str, str, data->numbin, lower, upper);
                monte_carlo(data, numcounts, numpeaks, vmu[i], sigma, lower, upper, data->numbin, simtype, trouble);
                hist2array(data);
                data->cslice++;
            }
        }// This simulates a linear detector on a partitioned time bases generating a single list output
        
        data->cslice = 0;
        data->indexvalue = 0;
        //window_save(data, 0, data->numbin);// (Obsolete)
        // Clean up
        if (trouble) troll = data;
        else {
            for (int i = 0; i < data->endslice; i++) data->histogram[i]->Delete();
            delete data;
        }
    } // Trouble Shooting:              Monte Carlo Simulation generation Laser with mu value or signal...
    else if(type == 29){
        /*/                                     NOTES (20161206)
         This section will be used to generate Monte Carlo simulation of the STJ signal response. I assume
         a poisson distribution to the peak amplitude and a single exponential decay
         /*/
        // Function Variables
        int numcounts = 30000;
        int numpeaks = 5; // represents minutes here
        double mu = 76.8; // represents centroid here
        char str[500];
        char *path = tpath;
        char *ext = text;
        bool trouble = true;
        bool nonlinear = false;
        double lower = 0;
        double upper = 819.2;
        double sigma = 0.8;
        double xaxis[tBin + 1];
        double noncent = 50;
        double nonline = -0.05/(noncent*noncent);
        info *data = new info;
        // End
        
        data->channel = 4;
        data->current = 0;
        data->cslice = 0;
        data->endslice = 100;
        data->numbin = tBin;
        snprintf(data->path, sizeof(str), "%s", mcpath);
        snprintf(data->name, sizeof(str), "Monte_Carlo_Simulation_Ch_%d_00_Window_Hist", data->channel);
        snprintf(data->ext, sizeof(str), "%s", ext);
        
        if (nonlinear) {
            for (int i = 0; i < int(sizeof(xaxis)/sizeof(xaxis[0])); i++) xaxis[i] = (upper - lower)*i/data->numbin + pow((upper - lower)*i/data->numbin - noncent, 2)*nonline;
            for(int i = 0; i < data->endslice; i++){
                snprintf(str, sizeof(str), "Window_%d", i);
                data->histogram[i] = new TH1D(str, str, tBin, xaxis);
                monte_carlo(data, numcounts, numpeaks, mu, sigma, lower, upper, data->numbin, 1, trouble);
                hist2array(data);
                data->cslice++;
            }
        }// This simulates a Non-linear detector
        else {
            for(int i = 0; i < data->endslice; i++){
                snprintf(str, sizeof(str), "Window_%d", i);
                data->histogram[i] = new TH1D(str, str, data->numbin, lower, upper);
                monte_carlo(data, numcounts, numpeaks, mu, sigma, lower, upper, data->numbin, 1, trouble);
                hist2array(data);
                data->cslice++;
            }
        }// This simulates a linear detector
        
        data->cslice = 0;
        data->indexvalue = 0;
        troll = data;
        window_save(data, 0, data->numbin);
    } // Trouble Shooting:              Monte Carlo Simulation generation signal (in Progress)...
    else if(type == 30){
        /*/                                     NOTES
         20170104
         This piece is intended to load the MCS listmode data and generate a linear and non-linear histogram.
         I set the range and bin size so that the numbers are more 'realistic' in that they are not rounded numbers.
         Time: 1148
         Set lower and upper bounds to -1/6.0 and 1350 this should provide an some offset and variation in the histograms
         
         20170516
         Updated the code to use a response function X(E) to generate histogram, set to work for up to cubic responses
         without the laser-substrate term (eoffset = 0). Need to work out how to use the eoffset term for non linear cases
         works for linear case though.
         
         For operating type 16 on the actual energy data use
            sigma_value = 3;
            area_value = 0.001;
            lal = 150;
            lcl =  10;
            ucl = 200;
            pspace = 2;
         
         For the responseivity values use:
            sigma_value = 5;
            area_value = 0.01;
            lal =  150;
            lcl =  300;
            ucl = 4500;
            pspace = 75;
         
         20170525
         Set the laser-substrate into the MCS list generation of process 28, so the responsivity cannot be set differently
         for the laser-detector and laser-substrates. Should be fine can work out additional details if needed later.
         /*/
        // Function Variables
        int counter = 0;
        int simtype = (sig <= 0 ? 0 : 1);
        int numparts = 100;
        int rebinval = 4;                       // 1 For use of actual energy values, else 4
        char str[500];
        char *path = mcpath;
        char *ext = text;
        bool trouble = false;
        //bool line = false;
        bool singular = true;               // Use false for moving window average...
        bool varmu = false;
        float parsetime = 5; // minutes
        double b[4] = {0, 30, -0.01, 0.0001};   // [Constant, Linear, Quadratic, Cubic] term, linear = 1 For use of actual energy values, else 30, constant = 0
        double lower = 0;                       // 1 For use of actual energy values, else 0
        double upper = tBin;                     // 217 For use of actual energy values, else tBin
        double timeresp = 0;                     // Use this to introduce time dependent responsivities
        double vmu[tBin];
        
        TF1 *fun;
        TF1 *efun;
        info *data = new info;
        ifstream file;
        string trash;
        // End

        data->channel = 4;
        data->current = 0;
        data->cslice = 0;
        data->endslice = 100;
        data->numbin = tBin;                    // 2160/2 For use of actual energy values, else tBin
        snprintf(data->path, sizeof(str), "%s%s/", mcpath, mcfile00);
        //snprintf(data->name, sizeof(str), "Monte_Carlo_Simulation_Ch_%d_00_Window_Hist", data->channel);
        snprintf(data->ext, sizeof(str), "%s", ext);
        
        
        if(varmu){
            printf("Using time dependent mu\n");
            snprintf(str, sizeof(str), "%sMonte_Carlo_Simulation/Var_mu.txt", mcpath);
            file.open(str);
            if(file.fail()){
                printf("Failed to open file:\n%s\n", str);
            }
            file >> trash;
            file >> trash;
            while(file >> trash >> vmu[counter]) counter++;
            data->endslice = counter;
            printf("Total number of mu entries, %d\n", counter);
            file.close();
        }
        
        snprintf(str, sizeof(str), "%sMonte_Carlo_Simulation/Laser_MCS_List.txt", mcpath);
        if (simtype == 1) snprintf(str, sizeof(str), "%sMonte_Carlo_Simulation/Signal_MCS_List.txt", mcpath);
        file.open(str);
        if(file.fail()) printf("Failed to open file:\n%s\n", str);
        file >> trash;
        file >> trash;
        data->current = 0;
        while(file >> data->time[data->current] >> data->bin[data->current]){
            data->tag[data->current] = 1; //(simtype == 0 ? 1 : 0);
            data->current++;
        }
        data->last = data->current;
        data->current = 0;
        printf("Starting time adjustments\n");
        if (simtype != 0) for (int i = 0; i < data->last; i++){
            data->time[i] = (long long)(parsetime*60*pow(10, 9)*(i/(data->last/numparts))) + (long long)(10000000*(i%(data->last/numparts)));
            if (i == data->last - 1) data->time[i] = (long long)(parsetime*60*pow(10, 9)*((i + 1)/(data->last/numparts))) - 1;
        }
        
        // Create Detector Response Function
        sig = (sig == 0 ? 1: TMath::Abs(sig));
        printf("Function Model type: pol%d\n", int(sig));
        snprintf(str, sizeof(str), "pol%d", int(sig));
        fun = new TF1("fun2", str);
        for (int i = 0; i < (sig + 1 < 4 ? sig + 1: 4); i++) fun->SetParameter(i, b[i]);
        
        // Create output histograms/partitions
        for(int i = 0; i < NUM_COLUMN; i++){
            sprintf(str, "Partition_%d_Channel_%d", i, data->channel);
            data->histogram[i] = new TH1D(str, str, data->numbin, lower, upper);
            data->histogram[i]->Rebin(rebinval);
        }
        snprintf(data->path, sizeof(str), "%s%s/", mcpath, mcfile00);
        snprintf(data->name, sizeof(str), "MCS_Response_Fcn_Ch_12_00_Window_Hist");
        counter = 0;
        for (int i = 0; i < data->last; i++){
            /*/if (varmu){
                if (counter - data->time[i]/(5*60*pow(10, 9)) < 0) counter++;
                eoffset = vmu[counter]*1.0/22;
            }
            if (line) data->bin[i] = fun->Eval(data->bin[i] + eoffset) - (eoffset == 0 ? 0 : (varmu ? 0: 0.25))*fun->Eval(data->bin[i])*data->time[i]/(8.64*pow(10, 13));// (Obsolete offset is handled in MCS list generation)
            else data->bin[i] = (1 - 0*0.25*data->time[i]/(8.64*pow(10, 13)))*fun->Eval(data->bin[i]); // Use for non-linear cases... (Obsolete offset is handled in MCS list generation)/*/// Obsolete
            data->bin[i] = (1 - timeresp*0.25*data->time[i]/(8.64*pow(10, 13)))*fun->Eval(data->bin[i]);
        }
        parse_time(data, 20, parsetime, singular);
        //Clean up
        for (int i = 0; i < NUM_COLUMN; i ++) data->histogram[i]->Delete();
        delete data;
    } // Trouble Shooting:              Bin List Monte Carlo Simulation in either linear or non-linear...
    else if(type == 31){
        // Function Variables:
        char str[500];
        char *path = tpath;
        char *ext = text;
        bool override = (sig < 0 ? false : true);
        bool trouble = false; //(sig == 0 ? false : true);
        double sigma_value = 3; // Change back to 15 for ORTEC files
        double area_value = 0.005;
        int folderlength = 27;
        int lal = 500;
        int lcl =  2500; //60;
        int ucl = 7500; //600;
        int pspace = 100;
        int counter = 0;
        double binstart = 0;
        double binend = 0;
        info *host = new info;
        info *data = new info;
        // End
        
        host->current = 0;
        host->last = 0;
        host->endslice = 0;
        host->channel = int(2*randomvalue);
        snprintf(host->path, sizeof(str), "%sFull_Experiment/", path);
        snprintf(host->name, sizeof(str), "Experiment-Laser_Ch_%d", host->channel);
        snprintf(host->ext, sizeof(str), "%s", ext);
        snprintf(host->filename, sizeof(str), "%sFull_Experiment/Experiment-Laser_Ch_%d", path, host->channel);
        
        
        for (int i = 0; i < folderlength; i++){
            if (i == 24) continue;
            data->current = 0;
            data->last = 0;
            data->endslice = 0;
            data->channel = int(2*randomvalue);
            snprintf(data->path, sizeof(str), "%s%s/", path, folder[i]);
            snprintf(data->name, sizeof(str), "%s", folder[i]);
            snprintf(data->ext, sizeof(str), "%s", ext);
            
            // Load to buffer array and fill histograms
            snprintf(data->filename, sizeof(str), "%s%s/%s_Ch_%d_04_Full_Laser.txt", path, folder[i], folder[i], data->channel); // Window file name...
            file_load(data, 3); // Windows file load channel domain
            for(int j = 0; j < data->last; j++) host->hist[host->cslice][j] = data->hist[0][j];
            host->cslice++;
            host->endslice++;
        }
        for(int j = 0; j < data->last; j++) host->xaxis[j] = data->xaxis[j];
        host->cslice = 0;
        window_save(host, 0, data->numbin); // save gain matched partitions
        
        host->cslice = 0;
        host->endslice = 0;
        snprintf(host->name, sizeof(str), "Experiment-Signal_Ch_%d", host->channel);
        for (int i = 0; i < folderlength; i++){
            if (i == 24) continue;
            data->current = 0;
            data->last = 0;
            data->endslice = 0;
            data->channel = int(2*randomvalue);
            snprintf(data->path, sizeof(str), "%s%s/", path, folder[i]);
            snprintf(data->name, sizeof(str), "%s", folder[i]);
            snprintf(data->ext, sizeof(str), "%s", ext);
            
            // Load to buffer array and fill histograms
            snprintf(data->filename, sizeof(str), "%s%s/%s_Ch_%d_06_Full_Signal.txt", path, folder[i], folder[i], data->channel); // Window file name...
            file_load(data, 3); // Windows file load channel domain
            for(int j = 0; j < data->last; j++) host->hist[host->cslice][j] = data->hist[0][j];
            host->cslice++;
            host->endslice++;
        }
        host->cslice = 0;
        for(int j = 0; j < data->last; j++) host->xaxis[j] = data->xaxis[j];
        window_save(host, 0, data->numbin); // save gain matched partitions
        
        
    } // Trouble Shooting:              Merge laser/signal files from a single full experiment...
    else if(type == 32){
        // Function Variables:
        int counter = 0;
        int eslot = 0;
        int standnum = 0; //50; // This is the partition to use as the standard when calculating the gain matching
        int lal = 200;
        bool trouble = false;
        bool specialize = true;
        char str[500];
        char *path = tpath;
        char *ext = text;
        float lstart = (0 <= sig ? 66.5 : -sig);
        double gain = 0;
        double dgai = 0;
        double offs = 0;
        double doff = 0;
        double poin = 0;
        double dpoi = 0;
        double centroid[NUM_PEAKS][NUM_COLUMN + 1]; // The additional term is used to keep track of position in filling and total counts
        double delta[NUM_PEAKS][NUM_COLUMN];
        double output[2];
        
        info *data = new info;
        TH1D *laser;
        TH1D *signal;
        // End
        if(true){
            data->current = 0;
            data->last = 0;
            data->endslice = 0;
            data->channel = 2*int(randomvalue);
            snprintf(data->filename, sizeof(data->filename), "%s%s/%s_Ch_%d_01_Parameters.txt", path, file, file, data->channel);
            sprintf(data->path, "%s%s/", path, file);
            sprintf(data->name, "%s", file);
            sprintf(data->ext, ".txt");
            file_load(data, 5); // Parameter file load
            data->cslice = 0;
            snprintf(data->filename, sizeof(data->filename), "%s%s/%s_Ch_%d_02_Calibration.txt", path, file, file, data->channel);
            file_load(data, 4); // Gain Match file load
        } // Initial setup
        if(true){
            data->cslice = 0;
            for (int i = 0; i < NUM_PEAKS; i++) for(int j = 0; j < NUM_COLUMN + 1; j++) centroid[i][j] = 0;
            for(int i = 0; i < data->endslice; i++){
                gain = data->calib[i][1];
                dgai = data->calib[i][2];
                offs = data->calib[i][3];
                doff = data->calib[i][4];
                counter = 0;
                while (counter < NUM_PEAKS){
                    eslot = int(data->param[NUM_PEAKS*i + counter][1]/3.5);
                    poin = data->param[NUM_PEAKS*i + counter][2];
                    dpoi = data->param[NUM_PEAKS*i + counter][3];
                    data->param[NUM_PEAKS*i + counter][2] = gain*poin + offs;
                    data->param[NUM_PEAKS*i + counter][3] = pow(dgai*dgai*poin*poin + gain*gain*dpoi*dpoi + doff*doff, 0.5);
                    data->param[NUM_PEAKS*i + counter][4] = 0;
                    data->param[NUM_PEAKS*i + counter][5] = 0;
                    data->param[NUM_PEAKS*i + counter][6] = 0;
                    data->param[NUM_PEAKS*i + counter][7] = 0;
                    if (counter != 0 && data->param[NUM_PEAKS*i + counter + 1][0] != -1){
                        centroid[eslot][int(centroid[eslot][NUM_COLUMN])] = data->param[NUM_PEAKS*i + counter][2];
                        delta[eslot][int(centroid[eslot][NUM_COLUMN])] = 1.0/data->param[NUM_PEAKS*i + counter][3];
                        centroid[eslot][NUM_COLUMN]++;
                    }
                    counter++;
                    if (data->param[NUM_PEAKS*i + counter][0] == -1) break;
                }
            }// Transform the original centroids to the gain matched values.
            data->indexvalue = 11;
            fit_save(data, false, 1);
        } // Calculate the gain matched terms and generate the master centroid and weight matrices
        if(true){
            for(int i = 0; i < NUM_PEAKS*NUM_COLUMN; i++) data->param[i][0] = -1;
            printf("Starting Averages Caluculations\n");
            counter = 0;
            for(int i = 0; i < NUM_PEAKS; i++) if (10 <= centroid[i][NUM_COLUMN]){
                printf("On Peak energy %0.2f with %d peaks\n", 3.5*i, int(centroid[i][NUM_COLUMN]));
                stat(output, centroid[i], delta[i], centroid[i][NUM_COLUMN]);
                data->param[counter][0] = 0;
                data->param[counter][1] = 3.5*i;
                data->param[counter][2] = output[0];
                data->param[counter][3] = output[1]*pow(centroid[i][NUM_COLUMN], -0.5); // Measurement Uncertainty
                data->param[counter][4] = output[1]; // Standard Deviation
                counter++;
            }
            data->cslice = 0;
            data->endslice = 3;
            calibration(data, lal, specialize, trouble, lstart);
            for (int i = 0; i < 7; i++) data->calib[1][i] = data->calib[0][i]; // Copy calibration to slot 1 with offset term
            for (int i = 0; i < 7; i++) data->calib[2][i] = (i < 3 ? data->calib[0][i] : 0); // Copy calibration to slot 2 no offset term
            data->indexvalue = 12;
            fit_save(data, false, 0); // Save parameters and calibration
            data->indexvalue = 14;
            data->endslice = 1;
            residual(data); // Save residuals
        } // Save the average weighted values for the centroids, std and calibration, the edge peaks were excluded in the prior step

        
    } // Trouble Shooting:              This will calculate the adjusted gain centroids and new centroids for final calibration...
    else if(type == 33){
        // Function Variables:
        int counter = 0;
        int eslot = 0;
        int standnum = 0;
        bool trouble = false;
        char str[500];
        char *path = tpath;
        char *ext = text;
        double binstart = 0;
        double binend = 0;
        
        info *data = new info;
        TH1D *laser;
        TH1D *signal;
        TCanvas *canvas;
        // End
        
        if (true){
            data->current = 0;
            data->last = 0;
            data->endslice = 0;
            data->channel = 2*int(randomvalue);
            snprintf(data->filename, sizeof(data->filename), "%s%s/%s_Ch_%d_01_Parameters.txt",
                     path, file, file, data->channel);
            sprintf(data->path, "%s%s/", path, file);
            sprintf(data->name, "%s", file);
            sprintf(data->ext, ".txt");
            file_load(data, 5); // Parameter file load
            data->cslice = 0;
            snprintf(data->filename, sizeof(data->filename), "%s%s/%s_Ch_%d_00_Window_Hist.txt",
                     path, file, file, data->channel);
            file_load(data, 3); // Window file load
            
            binstart = (3*data->xaxis[0] - data->xaxis[1])/2;
            binend = data->xaxis[data->last - 1] + (data->xaxis[1] - data->xaxis[0])/2;
            printf("\nNumber of bins: %d\nStart: %.2f\nEnd: %.2f\n", data->numbin, binstart, binend);
            printf("Number of slices: %d\n", data->endslice);
            for (int i = 0; i < data->endslice; i++){
                snprintf(str, sizeof(str), "Partition_%d_Channel_%d", i, data->channel);
                data->histogram[i] = new TH1D(str, str, data->numbin, binstart, binend); // Calibration invariant
                for (int j = 0; j < data->numbin; j++){
                    data->histogram[i]->Fill(data->xaxis[j], data->hist[i][j]);
                }
                data->histogram[i]->Sumw2(kFALSE);
                data->histogram[i]->Sumw2();
            }// Generate BIN domain histograms
            eslot = data->endslice;
        } // Initial setup
        
        if(true){
            for (int i = 0; i < data->endslice; i++){
                //i = 2;
                driftcor(data, binstart, binend, i, 0, trouble);
                data->indexvalue = i;
                fit_save(data, false, 2);// Save Drift correction gain values
                //break;
            }
        } // Drift Correction
        
    } // Troubleshooting:               Process an entire run using each partition as the standard and save the gain matched calibration values.
    else if(type == -100){
        for (int j = 0; j < 2; j++) for (int i = 1; i < 4; i++){
            process(15, 0, folder[i], j + 2);
            process(16, 0, folder[i], j + 2);
            process(17, 0, folder[i], j + 2);
            process(18, 0, folder[i], j + 2); }
    } // Mass Laser processing
    else if(type == -200){
        for (int i = 1; i < 10; i++){
            process(25, 0, folder[i], randomvalue);
        }
        for (int i = 1; i < 10; i++){
            process(26, 0, folder[i], randomvalue);
        }
    } // Mass  processing
    else if(type == -300){
        sprintf(tpath, "/Users/ponce10/Root_Folder/Data/239_Pu_to_235U/High_Activity_Four_Pinholes/160711-160817_STJ-d_LE_XIA_239Pu/");
        process(21, 0, folder[5], 3);
        for (int i = 16; i < 20; i++) process(21, 0, folder[i], 3);

    } // Mass  processing
    else if(type == -400){
        sprintf(tpath, "/Users/ponce10/Root_Folder/Data/239_Pu_to_235U/High_Activity_Four_Pinholes/160711-160817_STJ-d_LE_XIA_239Pu/");
        for (int i = 20; i < 25; i++) process(21, 0, folder[i], 3);
    }
    else if(type == -500){
        snprintf(tpath, sizeof(tpath), "%s", upath[3]);
        char str[500];
        char *path = tpath;
        char *ext = text;
        bool override = true;
        bool trouble = false;
        double sigma_value = 5;
        double area_value = 0.005;
        int lal = 1000;
        int lcl = 200;
        int ucl = 5000;
        int pspace = 75;
        int counter = 0;
        double binstart = 0;
        double binend = 0;
        info *data = new info;
        // End
        
        data->current = 0;
        data->last = 0;
        data->endslice = 0;
        data->channel = int(2*randomvalue);
        snprintf(data->path, sizeof(str), "%s%s/", path, file);
        snprintf(data->name, sizeof(str), "%s", file);
        snprintf(data->ext, sizeof(str), "%s", ext);
        
        // Load to buffer array and fill histograms
        snprintf(data->filename, sizeof(str), "%s%s/%s_Ch_%d_00_Window_Hist.txt", path, file, file, data->channel); // Window file name...
        file_load(data, 3); // Windows file load channel domain
        binstart = (3*data->xaxis[0] - data->xaxis[1])/2;
        binend = data->xaxis[data->last - 1] + (data->xaxis[1] - data->xaxis[0])/2;
        printf("\nNumber of bins: %d\nStart: %.2f\nEnd: %.2f\n", data->last, binstart, binend);
        printf("Number of slices: %d\n", data->endslice);
        for (int j = 0; j < NUM_PEAKS*NUM_COLUMN; j++) data->param[j][0] = -1;
        for (int i = 0; i < data->endslice; i++){
            snprintf(str, sizeof(str), "Partition_%d_Channel_%d", i, data->channel);
            data->histogram[i] = new TH1D(str, str, data->last, binstart, binend); // Calibration invariant
            for (int j = 0; j < data->last; j++){
                data->histogram[i]->Fill(data->xaxis[j], data->hist[i][j]);
            }
            data->histogram[i]->Sumw2(kFALSE);
            data->histogram[i]->Sumw2();
        }
        troll = data;
    } // Mass Signal processing
    std::cout << "\n";
    return type;
}

