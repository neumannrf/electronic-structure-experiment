#!/bin/bash


for PW in 700 750 800 850 900 950 1000 1050 1100 1150 1200 1250 1300 1350 1400 1500 1600 1700 1800 1900 2000 2500 3000 3500 4000; do

cd $PW

jbsub -queue x86_24h -proj EMRE -name MgMOF-$PW -mem 256G -cores 1x36 -out log.oe -err log.oe bash run.sh

cd ..

done

