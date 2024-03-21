#!/usr/bin/env python3
import sys

import sys
import os
import yaml
import pprint
import time
import datetime
import warnings

import optuna
import logging
import uproot

import pathlib
import matplotlib

matplotlib.use("pdf")
import matplotlib.pyplot as plt
import random
import subprocess
import multiprocessing
import numpy as np
import json
import array
import sys
import argparse
import pandas as pd
import matplotlib.pyplot as plt

from typing import Optional, Union
from pathlib import Path

srcDir = Path(__file__).resolve().parent


def run_ckf(params, names, outDir):
    if len(params) != len(names):
        raise Exception("Length of Params must equal names")

    ckf_script = srcDir / "ml_chain_na60plus.py"
    nevts = "-n 1"
    #indir = "--indir=" + str(srcDir)
    outdir = "--output=" + str(outDir)

    ret = ["python3"]
    ret.append(ckf_script)
    ret.append(nevts)
    #ret.append(indir)
    ret.append(outdir)

    i = 0
    for param in params:
        arg = "--sf_" + names[i] + "=" + str(param)
        ret.append(arg)
        i += 1
    print(ret)
    # Run CKF for the given parameters
    subprocess.call(ret)


class Objective:
    def __init__(self, k_dup, k_time, k_d0):
        self.res = {
            "eff": [],
            "fakerate": [],
            "duplicaterate": [],
            "runtime": [],
            "d0": [],
        }

        self.k_dup = k_dup
        self.k_time = k_time
        self.k_d0 = k_d0

    def __call__(self, trial, ckf_perf=True):
        params = []
        #example for categorical data
        #maxSeedsPerSpM = trial.suggest_categorical("true_false_variable", [True, False])
        #params.append(maxSeedsPerSpM)
        maxSeedsPerSpM = trial.suggest_int("maxSeedsPerSpM", 0, 10)
        params.append(maxSeedsPerSpM)
        cotThetaMax = trial.suggest_float("cotThetaMax", 5.0, 10.0)
        params.append(cotThetaMax)
        sigmaScattering = trial.suggest_float("sigmaScattering", 0.2, 50)
        params.append(sigmaScattering)
        radLengthPerSeed = trial.suggest_float("radLengthPerSeed", 0.001, 0.1)
        params.append(radLengthPerSeed)
        impactMax = trial.suggest_float("impactMax", 0.1, 25)
        params.append(impactMax)
        keys = [
            "maxSeedsPerSpM",
            "cotThetaMax",
            "sigmaScattering",
            "radLengthPerSeed",
            "impactMax",
        ]

        get_tracking_perf(self, ckf_perf, params, keys)
        efficiency = self.res["eff"][-1]
        penalty = (
            self.res["fakerate"][-1]
            #+ self.res["duplicaterate"][-1] / self.k_dup
            #+ self.res["runtime"][-1] / self.k_time
            #+ self.res["d0"][-1] / self.k_d0
        )

        return efficiency #- penalty


def get_tracking_perf(self, ckf_perf, params, keys):
    if ckf_perf:
        outDirName = "Optimization"
        outputfile = srcDir / outDirName / "performance_ckf.root"
        effContName = "particles"
        contName = "tracks"
    else:
        outDirName = "Output_Seeding"
        outputfile = srcDir / outDirName / "performance_seeding.root"
        effContName = "seeds"
        contName = "seeds"

    outputDir = Path(srcDir / outDirName)
    outputDir.mkdir(exist_ok=True)
    run_ckf(params, keys, outputDir)
    rootFile = uproot.open(outputfile)
    self.res["eff"].append(rootFile["eff_" + effContName].member("fElements")[0])
    
    self.res["fakerate"].append(rootFile["fakerate_" + contName].member("fElements")[0])
    self.res["duplicaterate"].append(
        rootFile["duplicaterate_" + contName].member("fElements")[0]
    )

    timingfile = srcDir / outDirName / "timing.tsv"
    timing = pd.read_csv(timingfile, sep="\t")

    if ckf_perf:
        time_ckf = float(
            timing[timing["identifier"].str.match("Algorithm:TrackFindingAlgorithm")][
                "time_perevent_s"
            ]
        )

    time_seeding = float(
        timing[timing["identifier"].str.match("Algorithm:SeedingAlgorithm")][
            "time_perevent_s"
        ]
    )
    if ckf_perf:
        self.res["runtime"].append(time_ckf + time_seeding)
    else:
        self.res["runtime"].append(time_seeding)
    outputfile = srcDir / outDirName / "performance_track_fitter_ckf.root"
    rootFile = uproot.open(outputfile)
    histogram = rootFile["res_d0"]

    # Get bin edges and bin contents
    bin_edges = histogram.edges
    bin_contents = histogram.values

        # Calculate the RMS
    # Note: Exclude the underflow and overflow bins if needed
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    rms = np.sqrt(np.sum(bin_contents * (bin_centers**2)) / np.sum(bin_contents))

    self.res["d0"].apend(rms)

def main(k_dup = 5, k_time = 5):

    # Initializing the objective (score) function
    objective = Objective(k_dup, k_time)

    outputDir = Path("OptunaResults")
    outputDir.mkdir(exist_ok=True)
    # Open and read JSON file
    try:
        with open(outputDir / "results.json",'r') as file:
            start_values = json.load(file)
    except:
        start_values = {
            'maxSeedsPerSpM': 3,
            'cotThetaMax': 5.809222379141632,
            'sigmaScattering': 7.293572849910582,
            'radLengthPerSeed': 0.07827007716884793,
            'impactMax': 10.47494916604592
        }

    
    # Optuna logger
    optuna.logging.get_logger("optuna").addHandler(logging.StreamHandler(sys.stdout))
    study_name = "test_study"

    # creating a new optuna study
    study = optuna.create_study(
        study_name=study_name,
        storage="sqlite:///{}.db".format(study_name),
        direction="maximize",
        load_if_exists=True,
    )

    study.enqueue_trial(start_values)
    # Start Optimization
    study.optimize(objective, n_trials=1)

    # Printout the best trial values
    print("Best Trial until now", flush=True)
    for key, value in study.best_trial.params.items():
        print(f"    {key}: {value}", flush=True)

    with open(outputDir / "results.json", "w") as fp:
        json.dump(study.best_params, fp)

    fig = optuna.visualization.plot_param_importances(study)
    fig.show()
    #plt.savefig("param_importance.png")

    fig = optuna.visualization.plot_optimization_history(study)
    fig.show()
    #plt.savefig("opt_history.png")

if __name__ == "__main__":
    main()