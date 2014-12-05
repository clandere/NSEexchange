#!/bin/bash

echo
echo "***************************************"
echo "* Testing model against benchmarks... *"
echo "***************************************"
echo

Rscript test_eta.R

Rscript test_simulated_eta_distbn.R

Rscript test_phi_prediction.R
