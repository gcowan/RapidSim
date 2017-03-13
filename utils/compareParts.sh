#!/bin/bash

wget http://svn.cern.ch/guest/evtgen/tags/R01-06-00/evt.pdl
cat evt.pdl | grep "^add  p Particle" > evtgenParts.txt
cat ../config/particles.dat | grep -v "^ID" > rapidsimParts.txt

python compareParts.py
