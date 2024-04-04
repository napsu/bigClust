#*************************************************************************
#*                                                                       *
#*     Makefile for BigClust, the stochastic nonsmooth clustering        *
#*     software by Napsu Karmitsa 2024 (last modified 03.04.2024).       *
#*                                                                       *
#*     The work was financially supported by the Research Council of     *
#*     Finland (Project No. #345804 and #345805).                        *
#*                                                                       *
#*     BigClust software is covered by the MIT license.                  *
#*                                                                       *
#*************************************************************************
# 
FF = gfortran -O3
OPEN = 
OBJ = parameters.o initbigclust.o clusteringmod.o objfun.o subpro.o slmb.o bigclust.o

all:	bigclust


parameters.o: parameters.f95
	$(FF) -c $(OPEN) parameters.f95

initclust.mod: r_precision.mod initbigclust.o initbigclust.f95
	$(FF) -c $(OPEN) initbigclust.f95

initslmb.mod: r_precision.mod initclust.mod initbigclust.o initbigclust.f95
	$(FF) -c $(OPEN) initbigclust.f95

initbigclust.o: r_precision.mod initbigclust.f95
	$(FF) -c $(OPEN) initbigclust.f95

clusteringmod.o: r_precision.mod param.mod initclust.mod clusteringmod.f95
	$(FF) -c $(OPEN) clusteringmod.f95

objfun.o: r_precision.mod initclust.mod objfun.f95
	$(FF) -c $(OPEN) objfun.f95

subpro.o: r_precision.mod param.mod subpro.f95
	$(FF) -c $(OPEN) subpro.f95

slmb.o: r_precision.mod param.mod initclust.mod initslmb.mod exe_time.mod obj_fun.mod subpro.mod slmb.f95
	$(FF) -c $(OPEN) slmb.f95

bigclust.o: r_precision.mod param.mod exe_time.mod initclust.mod initslmb.mod clusteringmod.mod slmb_mod.mod bigclust.f95
	$(FF) -c $(OPEN) bigclust.f95

bigclust: $(OBJ)
	$(FF) -o bigclust $(OPEN) $(OBJ)

clean:	
	rm bigclust $(OBJ) *.mod
	echo Clean done	