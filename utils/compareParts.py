#!/usr/bin/env python
from __future__ import print_function

import math

#get number of digits after the decimal
def scale(x):
    max_digits = 14
    int_part = int(abs(x))
    magnitude = 1 if int_part == 0 else int(math.log10(int_part)) + 1
    if magnitude >= max_digits:
        return (magnitude, 0)
    frac_part = abs(x) - int_part
    multiplier = 10 ** (max_digits - magnitude)
    frac_digits = multiplier + int(multiplier * frac_part + 0.5)
    while frac_digits % 10 == 0:
        frac_digits /= 10
    scale = int(math.log10(frac_digits))
    return scale

def compare(title, rs, eg):
    scaleEG = scale(eg)
    scaleRS = scale(rs)

    diff = rs-eg

    #check for zero
    if rs==0. and eg!=0.:
        print(title, "RS is zero", rs, eg)
        return True

    if rs!=0. and eg==0.:
        print(title, "EG is zero", rs, eg)
        return False

    #check for significant difference
    if abs(diff) >= 0.999*10.**-min(scaleRS,scaleEG):
        print(title, "DIFFERENT", rs, eg)
        return True

    #report which is more precise
    if scaleRS > scaleEG:
        print(title, "RS more precise", rs, eg)
        return False
    elif scaleRS < scaleEG:
        print(title, "EG more precise", rs, eg)
        return True
    else:
        print(title, "SAME precision", rs, eg)
        return True

from collections import namedtuple
Particle = namedtuple("Particle","id name mass width charge spin")

partsRS = {}
partsEG = {}

fileRS = open("rapidsimParts.txt")
linesRS = fileRS.readlines()
fileRS.close()

fileEG = open("evtgenParts.txt")
linesEG = fileEG.readlines()
fileEG.close()

for line in linesRS:
    pid = int(line.split()[0])
    name = line.split()[1]
    aname = line.split()[2]
    mass = float(line.split()[3])
    width = float(line.split()[4])
    charge = float(line.split()[5])
    spin = float(line.split()[6])
    partsRS[pid] = Particle(pid, name, mass, width, charge, spin)
    if(aname!="---"): partsRS[-pid] = Particle(-pid, aname, mass, width, -charge, spin)

for line in linesEG:
    pid = int(line.split()[4])
    name = line.split()[3]
    mass = float(line.split()[5])
    width = float(line.split()[6])
    charge = float(line.split()[8])/3.0
    spin = float(line.split()[9])/2.0
    partsEG[pid] = Particle(pid, name, mass, width, charge, spin)

rsKeys = set(partsRS.keys())
egKeys = set(partsEG.keys())

rsOnly = rsKeys.difference(egKeys)
egOnly = egKeys.difference(rsKeys)
both   = rsKeys.intersection(egKeys)

print("\n---------SUMMARY---------\n")
print("RS only:", len(rsOnly), "EG only:", len(egOnly), "both:", len(both), "\n")

print("\n---------DIFFERENCES---------\n")

diffM=0
diffG=0
diffC=0
diffS=0

for part in sorted(both):
    if part < 0 and -part in both:
        continue
    if partsRS[part].mass != partsEG[part].mass or \
       partsRS[part].width != partsEG[part].width or \
       partsRS[part].charge != partsEG[part].charge or \
       partsRS[part].spin != partsEG[part].spin:
        print(part, partsRS[part].name)
        if partsRS[part].mass != partsEG[part].mass:
            if compare("mass:   ", partsRS[part].mass, partsEG[part].mass):
                diffM+=1
        if partsRS[part].width != partsEG[part].width:
            if compare("width:  ", partsRS[part].width, partsEG[part].width):
                diffG+=1
        if partsRS[part].charge != partsEG[part].charge:
            print("charge: ", partsEG[part].charge, partsRS[part].charge)
            diffC+=1
        if partsRS[part].spin != partsEG[part].spin:
            print("spin:   ", partsEG[part].spin, partsRS[part].spin)
            diffS+=1
        print("\n")

print("\n---------DIFF SUMMARY---------")
print("M G C S")
print(diffM, diffG, diffC, diffS)

print "\n---------ONLY in EvtGen---------\n"
for part in sorted(egOnly):
    if part < 0 and -part in egOnly:
        continue
    print(part, partsEG[part].name)
