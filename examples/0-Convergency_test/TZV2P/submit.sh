#!/bin/bash


for PW in 600 1600 1700 1800 1900 2000 2500 3000 3500 4000; do

cd $PW

jbsub -queue x86_24h -proj EMRE -name MgMOF-$PW -mem 256G -cores 1x36 -out log.oe -err log.oe bash run.sh

cd ..

done

