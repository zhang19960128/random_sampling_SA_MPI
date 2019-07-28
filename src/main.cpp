#include "atom.h"
#include <iostream>
#include <fstream>
#include <string>
#include "readion.h"
#include <ctime>
#include "readpara.h"
#include "penalty.h"
#include "sa.h"
#include <mpi.h>
#include <time.h>
int main(int argc,char* argv[]){
	 MPI_Init(NULL,NULL);
	 readPT(argv[1]);
	 int size_box;
	 MPI_Barrier(MPI_COMM_WORLD);
   SimulatedAnnealing(&PenaltyFunc,
			 control::database[0],
			 control::xop,
			 control::paracount_bvv+control::paracount_charge,
			 saconst::sa_nt,
			 saconst::sa_ns,
			 saconst::sa_max,
			 saconst::sa_temp,
			 saconst::sa_ratio,
			 control::vm,
			 control::ub,
			 control::lb,
			 control::c);
		MPI_Finalize();
}
