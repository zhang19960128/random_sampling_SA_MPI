#include "sa.h"
#include "simann.h"
#include "atom.h"
#include "readpara.h"
#include <math.h>
#include <mpi.h>
#include <iostream>
#include "output.h"
#include <fstream>
#include <new>
namespace saconst{
 double sa_temp=0.0001;
 double sa_ratio=0.95;
 double sa_max=3;
 double sa_eweight=5.0;
 double sa_fweight=1.0;
 double sa_sweight=0.0;
 double sa_nt=3;
 double sa_ns=3;
 int sa_atom_num;
};
//extern void WriteFiles();
#define NPAR 12
#define NEPS 20
#define EPS 1.0e-6
void SimulatedAnnealing(double (*PenaltyFunc)(double*, box*,int,int,int), 
			box* system,
      double *xacc,
			int    N,//number of parameters to be optimized.
			int    NT,
			int    NS,
			int    sa_max,
			double sa_temp,
			double sa_ratio,
			double *vm,
			double *ub,
			double *lb,
			double *c)
{
    /*c is an array with every element to be 2.0*/
   int j,h,h0,m;
   double penaltyACC, penaltyOPT;
   double penaltySTAR[NEPS];
   int QUIT;
   double vmratio;
   double metropolis;
   double penaltyp;
   double sampling=new double control::database.size();
   double* optimize_temp;
   double ratio,ratiop;
   int nnew,nfcnev,*nacp;
   int iseed1, iseed2;
   int saiter;
   double *xopt, *xp;

   nacp = new int [N];
   xopt = new double [N];
   xp   = new double [N];

   /* Simulated Annealing output */
	 int world_rank;
	 std::fstream safile,optout;
	 MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);
	 if(world_rank==0){
	 	safile.open("sa.out",std::fstream::out);
	 	optout.open("opt.out",std::fstream::out);
	 }
   /* initilization of xopt */
   for(h0=0;h0<N;h0++)
      xopt[h0] = xacc[h0];

   penaltyACC = 1.e50;
   penaltyOPT = 1.e50;
   nnew = nfcnev = 0;
   for(h0=0;h0<N;h0++)
      nacp[h0] = 0;

   for(h0=0;h0<NEPS;h0++)
      penaltySTAR[h0] = 1.0e+20;
   /*                                      *
    * Simulated Annealing input parameters *
    *                                      */
   iseed1 = 1; iseed2 = 2;
   rmarin(&iseed1,&iseed2);
   for(saiter=0;saiter<sa_max;saiter++){
      /*                *
       *  m, j, h loop  *
       *                */
      for(m=0;m<NT;m++){
	 for(j=0;j<NS;j++){
            for(h=0;h<N;h++)
	       {

		  /* x -> x + dx */
		  if(nfcnev != 0)
		     {
			for(h0=0;h0<N;h0++){
			   if(h0==h)
			      {
				 xp[h0] = xacc[h0] + (ranmar()*2-1.0) * vm[h0];
				 if((xp[h0] < lb[h0])||(xp[h0] > ub[h0]))
				    xp[h0] = lb[h0] + (ub[h0]-lb[h0])*ranmar();
			      }
			   else
				 xp[h0] = xacc[h0];
			}
		     }
		  else{
		     for(h0=0;h0<N;h0++)
						xp[h0] = xacc[h0];
			}
			MPI_Bcast(xp,N,MPI_DOUBLE,0,MPI_COMM_WORLD);
			/*calculating penalty*/
			penaltyp=0.0;
		  for(size_t i=0;i<control::ionsize.size();i++){
		    optimize_temp= PenaltyFunc(xp,control::database[i],control::ionsize[i],control::minienergytick[i],i);//Zhenbang
            penaltyp=optimize_temp[0]+penaltyp;
            sampling[i]=optimize_temp[1];
          }
		  nfcnev += 1;
			if(world_rank==0){
		  if (penaltyp < penaltyACC){
				for(h0=0;h0<N;h0++)
			   	xacc[h0] = xp[h0];
				nacp[h] += 1;
				penaltyACC = penaltyp;
				if (penaltyp < penaltyOPT){
			  	 penaltyOPT=penaltyp;
			   	 for(h0=0;h0<N;h0++)
			     	 xopt[h0] = xp[h0];
			       nnew++;
				 	   write_opt_parameter(optout);
		   		   safile<<"  NEW OPTIMUM "<<penaltyOPT<<"\t"<<saiter<<"\t"<<m<<"\t"<<j<<"\t"<<h<<"\t"<<xp[h]<<" "<<"the sampling is: ";
                   for(size_t i=0;i<control::ionsize.size();i++){
                        saiter<<sampling[i]<<" ";
                   }
                   saiter<<std::endl;
						for(size_t i=0;i<control::ionsize.size();i++){
							write_defile(i);
						}
				}
		     }
		   else{
			/* Metropolis */
					metropolis = (penaltyACC - penaltyp)/sa_temp;
					ratio = exprep(metropolis);
					ratiop =  ranmar();
					if (ratiop < ratio){
			   		for(h0=0;h0<N;h0++){
			      		xacc[h0] = xp[h0];
			   		}
			   	penaltyACC = penaltyp;
			   	nacp[h] += 1;
						}
		     }
			}
	  } /* h */
	 } /* j */
	 
	 
	 /* Adjust vm so that approximately half of all evaluations are accepted. */
	 if(world_rank==0){
	 for(h0=0;h0<N;h0++){
	    vmratio = (double)(nacp[h0])/(double)NS;
	    if (vmratio > 0.6){
	       vm[h0] = vm[h0] * (1.0 + c[h0]*(vmratio-0.6)/0.4);
	    }
	    else if (vmratio < 0.4){
	       vm[h0] = vm[h0] / (1.0 + c[h0]*(0.4 - vmratio)/0.4);
	    }
	    if (vm[h0] > (ub[h0]-lb[h0])){
	       vm[h0] = ub[h0] - lb[h0];
	    }
	    nacp[h0] = 0;
	 }
	 }
      } /* m */
	if(world_rank==0){
      for(h0=NEPS-1;h0>0;h0--)
       	  penaltySTAR[h0] = penaltySTAR[h0-1];
      penaltySTAR[0] = penaltyACC;
      for(h0=0;h0<NEPS;h0++)
	 		QUIT = 1;
      for(h0=0;h0<NEPS;h0++){
         if( fabs(penaltyACC-penaltySTAR[h0]) > EPS ) QUIT = 0;
				 std::cout<<penaltySTAR[h0]<<"the QUIT signal is: "<<QUIT<<std::endl;
			}
      /* Quit if penaltyfunction is optimized !! */
      
      sa_temp = sa_temp * sa_ratio;
      
      
      penaltyACC = penaltyOPT;
      for(h0=0;h0<N;h0++)
         xacc[h0] = xopt[h0];
	}
			MPI_Bcast(&QUIT,1,MPI_INT,0,MPI_COMM_WORLD);
			if(world_rank==0){
				std::cout<<saiter<<std::endl;
			}
      if (QUIT == 1) break;
 } /* saiter */
	 if(world_rank==0){
   	safile.close();
		optout.close();
	 }
   delete [] nacp;
   delete [] xopt;
   delete [] xp;
}
