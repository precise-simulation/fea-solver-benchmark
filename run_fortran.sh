#!/bin/bash
# FEM Benchmark Fortran test run script.

# Copyright 2013-2022 Precise Simulation, Ltd.


################################################################
#   INPUT PARAMETERS
################################################################

# line=$(cat testrun_param.txt | cut -d " " -f1)
# args=(${line// / })
# x="${args[0]}"
# y="${args[1]}"
# echo $x # 13
# echo $y # 243

# NX="32 64 128 256"

rm -f *.log

NX0=$(cat testrun_param.txt | cut -d " " -f1 | sed 's/\r$//')
NLEV=$(cat testrun_param.txt | cut -d " " -f2 | sed 's/\r$//')
NLOOPS=$(cat testrun_param.txt | cut -d " " -f3 | sed 's/\r$//')
NX=()
for (( i=1; i<=$NLEV; i+=1 )); do
    NX+=" "$(expr $NX0 \* $((2**($i-1))) )
done

# echo "${NX[@]}"
# exit 0
# read
################################################################

cd src_fortran
make

for N in $NX ; do

cat >data/featfem.dat <<END_OF_DATA
======================================================================
File for input data for featfem
======================================================================
unit numbers and file names on /FILES/:
----------------------------------------------------------------------
'data/featfem.out'    CFILE  (name of protocol file)
----------------------------------------------------------------------
Values for /OUTPUT/ :
----------------------------------------------------------------------
0              M
0              MT
0              ICHECK
2              MSHOW  (1=file output,2=1+term. output)
0              IOUT   (write solution to GMV file)
----------------------------------------------------------------------
Values for grid :
----------------------------------------------------------------------
${N}           NX     (number of horizontal cells)
-1             NY     (number of vertical cells)
1D0            LX     (rectangle horizontal length)
1D0            LY     (rectangle vertical length)
2              ISORT  (grid ordering 1=xyz, 2=feat, 3=CM)
----------------------------------------------------------------------
Values for assembly :
----------------------------------------------------------------------
1              IASM   (1=FEAT2D assembly, 2=custom)
4              ICUB   (Cubature formula)
----------------------------------------------------------------------
Values for linear solver :
----------------------------------------------------------------------
6              ISOL   (1=JAC, 2=GS, 3=SOR, 4=4-Col. GS, 5=PCG, 6=MG)
1000           NIT    (maximum number of linear iterations)
0.7D0          OMEGA  (relaxation parameter for ISOL=1,3,4)
1D-6           EPS    (tolerance/stopping criteria)
----------------------------------------------------------------------
Values for multigrid parameters :
----------------------------------------------------------------------
9              NLEV   (number of grid levels)
0              ICYCLE (mg-cycle: 0=F, 1=V, 2=W)
1              IRELMG (monitor relative (or absolute) defect)
2 3              ISM    (smoother: 1=JAC, 2=SOR(GS if RLXSLC=1), 3=SOR 4-Col, 4=SSOR)
1              NPRESM (number of presmoothing steps)
1              NPOSSM (number of postsmoothing steps)
0.8D0          RLXSMF (relaxation for the mg smoother)
1              NSMFAC (factor for pre/postsm. steps on coarser levels)
3              IREST  (restriction op: 1=injection, 2=half w., 3=full w.)
1              ISOLC  (coarse grid solver 1=JAC, 2=SOR(GS if RLXSLC=1))
1.0D-3         DMPSLC (damping of residuals for mg-it.)
500            NSLC   (maximum number of coarse grid solver iterations)
0.8D0          RLXSLC (relaxation for the coarse grid solver)
0              IMTIME (check mg timings)
----------------------------------------------------------------------
END_OF_DATA

  for (( i=1; i<=$NLOOPS; i+=1 )); do
    echo "Iteration: $i"
    ./featfem
  done

done

echo
echo '------   Simulations finished !!!   ------'
echo

mv *.log ../output
cd ..
