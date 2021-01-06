/* ----------------------------------------------------------------------
i LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
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
//
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
    if (narg < 6) {
        printf("Syntax: c++_driver Proc in.lammps tau Pe gamma\n");
        exit(1);
    }
    int me,nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD,&me);
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    
    int nprocs_lammps = atoi(arg[1]);
    double tau=atof(arg[3]);
    double Pe=atof(arg[4]);
    double gamma=atof(arg[5]);   
    // The arguments do not include rho since all simulations will be done at rho=0.5 for now.
    // Rho, N, L are declared in the lammps input file called AOUP_input.in
     

    // Declare a few more variables.
    double dt=0.0005;
    double rho=0.5;
    int L=25;
    int t_total=20;  // Total time in the production run.
    double frac=0;  // Fraction of inactive particles. Code can be used to run simulations where only a fraction of particles are active.
   
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
    
    
    LAMMPS *lmp;
   
    // Setting random number generator to set initial v_a and to generate noise terms during simulation through the fix_AOUP files

    time_t seed1, seed2;
    seed1 = 635247624; // should be changed every simulation.
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

   sprintf(inputline,"fix       1 all bd_AOUP 1.0 %f %f %f %f %li", gamma, tau, Pe, frac, seed2);  
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

    // Equilibration run, for 50 tau

    do{
        loopi+=1;
        lmp->input->one("run 1 post no");
    }
    while (loopi<1);

    sprintf(inputline,"run  %i",(int)(floor(50*tau/dt)));
    if (lammps==1) lmp->input->one(inputline);

    // To have a more compact and easy-to-read code, three back-to-back simulations of length t_total = 20 will be run. During the first,
    // we will calculate rate of work, averaged every time step. During the second we will calculate g(r), averaged every t = 0.01. During 
    // the third, we will export snapshots of the particles every t = tau.

    // Calculate rate of work \dot{w} = -<f*v_a>, averaged every time step.  The rate of work is calculated over active as well as passive
    // particles, and printed out to a file. The rate of work over passive particles will naturally be zero.

    FILE *fileout;
    sprintf(outputfile,"dw.AOU.density%.2f.Tau%.2f.Pe%.2f.Gamma%.2f",rho,tau,Pe,gamma);
    fileout=fopen(outputfile,"w");

    double dw [2];
    double **f=lmp->atom->f;


    for (int loopk=0;loopk<t_total*tau/dt;loopk++){

        for (int loopi=0;loopi < npassive;loopi++){
            dw[0]+=0;
        }
        for (int loopi=npassive;loopi < nlocal;loopi++){
            dw[1]+=-1/gamma*(f[loopi][0]*v_a[loopi][0]+f[loopi][1]*v_a[loopi][1])/((nlocal-npassive)*t_total*tau/dt);
        }
        lmp->input->one("run 1");
    }
    fprintf(fileout,"%f\t%f\n", dw[0], dw[1]);
    fflush(fileout);

    
    FILE *fileout2;
    sprintf(outputfile,"gr.AOU.density%.2f.Tau%.2f.Pe%.2f.Gamma%.2f",rho,tau,Pe,gamma);
    fileout2=fopen(outputfile,"w");

    // Calculate g(r) for r between 0.5 and 3.5, averaged every 0.01

    double r[300];
    double gr[300];

    for (int i=0; i<300; i++){
        r[i]=0.01*(i+1) + 0.5;
        gr[i]=0;
    }
   

    double **x=lmp->atom->x;
    int counter1=0;

    double dtgr =  0.01; // time interval for capturing g(r) to average it

    for (int loopk=0;loopk<(int)(t_total*tau/dtgr);loopk++){

        sprintf(inputline1,"run  %i",(int)(dtgr/dt));
        if (lammps==1) lmp->input->one(inputline1);

        counter1+=1;

        //The g(r) is calculated between active particles and the passive bath, in general. If all particles are driven,
        //however, it will become a regular g(r)

        for (int loopi=npassive;loopi < nlocal;loopi++){
            for (int i=0; i<300; i++){
                for  (int loopj=0;loopj < nlocal; loopj++){
                    double dx=x[loopj][0] - x[loopi][0];
                    double dy=x[loopj][1] - x[loopi][1];
                    if(dx > L/2) dx-=L;
                    if(dx < -L/2) dx+=L;
                    if(dy > L/2) dy-=L;
                    if(dy < -L) dy+=L;
                    if ( r[i]*r[i] <= dx*dx+dy*dy && r[i+1]*r[i+1] > dx*dx+dy*dy){
                        gr[i]+=1.0/(nlocal-npassive);
                    }
                }
            }
        }
    }
    for (int i=0; i<300; i++){
        fprintf(fileout2,"%f\t%f\n", r[i], gr[i]/(counter1*2*3.14*r[i]*(r[i+1]-r[i])*rho));
    }
    fflush(fileout2);


    //Save coordinates every t = tau
    
    FILE *fileout3;
    sprintf(outputfile,"positions.AOU.density%.2f.Tau%.2f.Pe%.2f.Gamma%.2f",rho,tau,Pe,gamma);
    fileout3=fopen(outputfile,"w");

    int counter2=0;

    for (int loopk=0;loopk<(int)(t_total);loopk++){

        sprintf(inputline1,"run  %i",(int)(tau/dt));
        if (lammps==1) lmp->input->one(inputline1);

        counter2+=1;
        fprintf(fileout3,"Snapshot\t%i\n", counter2);
        for (int loopi=0;loopi < nlocal;loopi++){
            fprintf(fileout3,"%f\t%f\n", x[loopi][0], x[loopi][1]);

        }
    }
    fflush(fileout3);



    if (lammps == 1) delete lmp;

    // close down MPI
    
    MPI_Finalize();
}



