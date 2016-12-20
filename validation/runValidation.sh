#!/bin/bash

cd ../build
src/RapidSim.exe ../validation/B2Kee 100000
src/RapidSim.exe ../validation/Bs2Jpsiphi 100000
src/RapidSim.exe ../validation/Bd2D0rho0 100000
src/RapidSim.exe ../validation/Bs2D0Kpi 10000
src/RapidSim.exe ../validation/D02Kpi 100000

cd ../validation
root -b -q -l "compareHistograms.C(\"B2Kee\")"
root -b -q -l "compareHistograms.C(\"Bs2Jpsiphi\")"
root -b -q -l "compareHistograms.C(\"Bd2D0rho0\")"
root -b -q -l "compareHistograms.C(\"Bs2D0Kpi\")"
root -b -q -l "compareHistograms.C(\"D02Kpi\")"
