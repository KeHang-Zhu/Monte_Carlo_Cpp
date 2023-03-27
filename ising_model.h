/*Hearders*/
#include <stdio.h>
#include <iostream>
#include <sys/time.h>
#include <stdlib.h>
#include <math.h>
#include <numeric>

// -- common parameters and variables ------------------------------

const double tm32 = pow(2.0, -32);
const double eps = pow(1.0, -14); // very small number
bool prtflg;                      // flag for write2file
const int Mxint = 2147483647;     // maximum integer
const int Mnint = -2147483647;    // minimum integer

// -- observe ----------------------------------------------------------

const int NObs = 2;
int Npb; // # samples per block
int Nbk; // the {Nbk}th block is the current one
int Np;  // the pointer

// -- parameters and variables -------------------------------------
const double Pi = 3.1415926536;
const int D = 2; // dimensionality
int Lx, Ly, Vol; // the system size
                 // Lx: lattice size
int *sp;
double t;    // spin spin coupling
double h;    // external magnetic field
double beta; // temperature
double *p_flip;
double *dE_flip;

// -- Lattice and State --------------------------------------------
int nnb;  //  # neighboring vertices
int nnb1; // nnb+1
int nnb2; // nnb/2
int *nnsite;
int *backdir;
double p, dE;
double E;
int MP;
int wormstep;

/*Global parameters for calculating error and correlations*/
const int Max_block = 4000;
const int N_block = 1000;
const double Max_cor = 0.1;
const int N_print = 100;
const int Ntoss = 2000;
const int Nsamp = Ntoss * 10;
const int N_measure = Ntoss / 10;

int itoss;

bool flg_cor;

struct Obs
{
    std::string nam;       // the name of observables
    double vnr;            // current value (integer)
    double val;            // the accumulated value of observables
    double cor;            // correlation
    double erb;            // errobar
    double blk[Max_block]; // blocks for calculating the errorbars
    double a;              //
    double b;              // obser=a*val+b
};

// observables
Obs *Obser;

void initialize();
void allo_aux();
void init_cnf();
void markov();
void worm();
void def_prob();
void def_latt();
double getE();
double getMP();
void init_stat();
void measure();
void set_Obs();
void init_Obs(int i_, std::string nam0_, double a0, double b0);
void merge_blk();
bool cal_cor_dev();
int spin2index(int d_spin, int s);
void printing();
bool check();
void coll_data();
