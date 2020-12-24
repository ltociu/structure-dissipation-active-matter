/* ----------------------------------------------------------------------
 LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
 www.cs.sandia.gov/~sjplimp/lammps.html
 Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories
 
 Copyright (2003) Sandia Corporation.  Under the terms of Contract
 DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
 certain rights in this software.  This software is distributed under
 the GNU General Public License.
 
 See the README file in the top-level LAMMPS directory.
 ------------------------------------------------------------------------- */

// c++_driver = simple example of how an umbrella program
//              can invoke LAMMPS as a library on some subset of procs
// Syntax: c++_driver P in.lammps
//         P = # of procs to run LAMMPS on
//             must be <= # of procs the driver code itself runs on
//         in.lammps = LAMMPS input script
// See README for compilation instructions

#include <ctime>
#include <iostream>
#include <random>
#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "math.h"
#include "string.h"
#include "lammps.h"         // these are LAMMPS include files
#include "input.h"
#include "atom.h"
#include "library.h"
#include "modify.h"
#include "fix.h"
#include "assert.h"
#include "random_park.h"

using namespace LAMMPS_NS;


int main(int narg, char **arg)
{
    // setup MPI and various communicators
    // driver runs on all procs in MPI_COMM_WORLD
    // comm_lammps only has 1st P procs (could be all or any subset)
    
    MPI_Init(&narg,&arg);
    char outputfile[100];
    char inputline[1000];
    char inputline1[1000];
    if (narg == 1) {
        printf("Syntax: c++_driver Proc in.lammps Tau gamma\n");
        exit(1);
    }
    int me,nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD,&me);
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    
    int nprocs_lammps = atoi(arg[1]);
    double tau_x3 = atof(arg[3]);
    double tau = tau_x3 / 3.0;
    double Pe = atof(arg[4]);
    double gamma = atof(arg[5]);   // Rho, N, L are declared in the lammps input file.
     

    // Declare a few more variables.
    double dt = 0.000001 * gamma;
    double rho = RHO;  // This is set by the job submission script. This should be the only instance of the word RHO in this file.
    int t_taus = int(round(25*gamma));  // Production runtime in approximate units of \tau.
    int t_total = t_taus * tau_x3 * 0.4; // Total time in the production run.
    double frac = 0;  // Fraction of inactive particles
   
    if (nprocs_lammps > nprocs) {
        if (me == 0) //if (me==0) is true when the processor has id 0.
            printf("ERROR: LAMMPS cannot use more procs than available\n");
        MPI_Abort(MPI_COMM_WORLD,1);
    }
    
    int lammps;
    if (me < nprocs_lammps) lammps = 1;
    else lammps = MPI_UNDEFINED;
    MPI_Comm comm_lammps;
    MPI_Comm_split(MPI_COMM_WORLD,lammps,0,&comm_lammps);
    
    // open LAMMPS input script
    
    FILE *fp;
    if (me == 0) {
        fp = fopen(arg[2],"r");
        if (fp == NULL) {
            printf("ERROR: Could not open LAMMPS input script\n");
            MPI_Abort(MPI_COMM_WORLD,1);
        }
    }
    

    // run the input script thru LAMMPS one line at a time until end-of-file
    // driver proc 0 reads a line, Bcasts it to all procs
    // (could just send it to proc 0 of comm_lammps and let it Bcast)
    // all LAMMPS procs call input->one() on the line
    
    LAMMPS *lmp;
   
    // Setting random number generator to set initial v_a

    time_t seed1; time_t seed2;
    seed1 = time(NULL);
    seed2 = seed1/2;
    RanPark *ran = new RanPark(lmp,seed1);
    sprintf(inputline,"Seed1 is %li and Seed2 is %li" ,seed1, seed2);

    std::cout << inputline <<std::endl;

    char *arg1 = "-screen";
    char *arg2 = "none";
    char *commargs[2]={arg1,arg2};
    if (lammps == 1) lmp = new LAMMPS(1,commargs,comm_lammps);
    // SEE EARLIER PART OF CODE FOR CONDITIONS WHEN lammps!=1. This happens only when input is faulty.
    
    
    int n,loopi,loopj;
    char line[1024];
    while (1) {
        if (me == 0) {
            if (fgets(line,1024,fp) == NULL) n = 0;
            else n = strlen(line) + 1;
            if (n == 0) fclose(fp);
        }
        MPI_Bcast(&n,1,MPI_INT,0,MPI_COMM_WORLD);
        if (n == 0) break;
        MPI_Bcast(line,n,MPI_CHAR,0,MPI_COMM_WORLD);
        if (lammps == 1) lmp->input->one(line);
    }
    


   // Add lines to the input file and read them. Use this to make sure the Pex, Pey, gamma and tau you set in this file are the same as those in the lammps input. 


   sprintf(inputline,"timestep %f",dt);
   if (lammps==1) lmp->input->one(inputline);

   sprintf(inputline,"fix       1 all bd_ABP 1.0 %f %f %f %f %li", gamma, tau, Pe, frac, seed2);  
   if (lammps==1) lmp->input->one(inputline);

   int natoms = static_cast<int> (lmp->atom->natoms);
   int nlocal= static_cast<int> (lmp->atom->nlocal);
   int npassive=(int)(floor(nlocal*frac));

  //Initialize v_a and noise


  double **v_a = lmp->atom->v_a;
  double **noise=lmp->atom->noise;
  double fd_term = sqrt(2 * dt * 1 / (gamma)); // t_target has to be set to 1 for this.

  std::cout << "Setting v_a values" <<std::endl;
  for (int i = 0; i < npassive; i++) {     
     v_a[i][0] = 0;
     v_a[i][1] = 0;
  }
  for (int i = npassive; i < nlocal ; i++) {   	  
     v_a[i][0] = 0;
     v_a[i][1] = 0;
  }
  for (int i = 0; i < nlocal ; i++) {     
     noise[i][0] = fd_term * ran->gaussian();
     noise[i][1] = fd_term * ran->gaussian();
     noise[i][2] = fd_term * ran->gaussian();  
  }

    loopi=0;

    // Equilibration runs ...

    do{
        loopi+=1;
        lmp->input->one("run 1 post no");
    }
    while (loopi<1);


    sprintf(inputline,"run  %i",(int)(round(0.2*t_total/dt)));
    if (lammps==1) lmp->input->one(inputline);

    /*
    // Code for visualizing trajectories
    // Write restart files every 25 tau
    // sprintf(inputline,"restart         %d aoup.R%.2f.TA%.2f.U%.2f.restart", int(round(25*tau/dt)), rho, Pe, tau);
    // if (lammps == 1) lmp->input->one(inputline);
    // Dump every tau/dt during production run
    sprintf(inputline,"dump            15 all atom %d aoup.R%.2f.U%.2f.TA%.2f.lammpstrj", int(round(0.01*0.125*gamma*tau_x3*0.4/dt)), rho, tau, Pe);
    if (lammps == 1) lmp->input->one(inputline);
    sprintf(inputline,"dump            16 all custom %d aoup.R%.2f.U%.2f.TA%.2f.forces.lammpstrj id fx fy vx vy", int(round(0.01*0.125*gamma*tau_x3*0.4/dt)), rho, tau, Pe);
    if (lammps == 1) lmp->input->one(inputline);
    // Perform equilibration
    sprintf(inputline,"run  %i",(int)(round(gamma*tau_x3*0.4/(2*dt))));
    if (lammps==1) lmp->input->one(inputline);
    // Start the production run
    sprintf(inputline,"run  %i",(int)(round(t_total*tau_x3*0.4/(4*dt))));
    if (lammps == 1) lmp->input->one(inputline);
    // Write a restart file
    sprintf(inputline,"write_restart   aoup.R%.2f.U%.2f.TA%.2f.restart.*", rho, tau, Pe);
    if (lammps == 1) lmp->input->one(inputline);
    // End code for visualizing trajectories
     */

    
    // Code for measuring pair correlation functions
    FILE *fileout2;
    sprintf(outputfile,"gr.AOU.density%.2f.Tau%.2f.Pe%.2f.Gamma%.2f.stats",rho,tau_x3,Pe,gamma);
    fileout2=fopen(outputfile,"w");

    double box_size = 150.0;
    int num_bins = 15000;
    double dr = 0.005;
    double rmax = dr * (double) num_bins;
    double r[num_bins];
    double gr[num_bins];
    for (int i=0; i<num_bins; i++){
        r[i] = dr*i;
        gr[i] = 0.0;
    }


    double **x=lmp->atom->x;

    int counter1=0;

    for (int loopk=0;loopk<(int)(t_taus/0.25);loopk++){

        sprintf(inputline1,"run  %i",(int)(0.25*tau_x3*0.4/dt));
        if (lammps==1) lmp->input->one(inputline1);

        counter1+=1;
        for (int loopi=0;loopi < nlocal;loopi++){
            //for (int i=0; i<=num_bins; i++){
            for (int loopj=0;loopj < nlocal; loopj++){
                double dx=x[loopj][0] - x[loopi][0];
                double dy=x[loopj][1] - x[loopi][1];
                if(dx > (box_size/2.0)) dx-=box_size;
                if(dx < -(box_size/2.0)) dx+=box_size;
                if(dy > (box_size/2.0)) dy-=box_size;
                if(dy < -(box_size/2.0)) dy+=box_size;
                if (loopi != loopj && fmax(fabs(dx), fabs(dy)) <= rmax) {
                    double rij = sqrt(dx*dx + dy*dy);
                    int i = (int) (rij/dr);
                    if (i < num_bins) {
                        gr[i] += 1.0/nlocal;
                    //if (loopi != loopj && r[i]*r[i] <= dx*dx+dy*dy && (r[i]+dr)*(r[i]+dr) > dx*dx+dy*dy){
                    //    gr[i]+=1.0/nlocal;
                    }
                }
            }
        }
    }
    double normalize;
    fprintf(fileout2,"%f\t%f\n", (r[0] + dr/2), 0.0);
    for (int i=1; i<num_bins; i++){
        // normalize = counter1 * M_PI * (2*r[i]*dr + dr*dr) * rho;
        normalize = counter1 * M_PI * (2*r[i]*dr) * rho;
        fprintf(fileout2,"%f\t%f\n", (r[i] + dr/2), (gr[i]/normalize));
    }
    fflush(fileout2);
    // End code for measuring pair correlation functions
    

    if (lammps == 1) delete lmp;

    // close down MPI
    
    MPI_Finalize();
}



