#!/usr/bin/env python3
"""
Based on the Swave software package (https://github.com/songbowang125/Swave).
Originally licensed under the GPL-3.0.

Publication:
Wang, S., Xu, T., Zhang, P. & Ye, K. Population-level structural variant 
characterization using pangenome graphs. Nat Genet (2026). 
https://doi.org/10.1038/s41588-026-02538-6

Modified and refactored for Nextflow integration.
Copyright (c) 2026 Jonah Kapski <Jonah.Kapski@edu.ruhr-uni-bochum.de>
"""

import sys
import os
import logging
import argparse

from src.swave.predict import process_projections_to_predictions


logging.basicConfig(
    level=logging.INFO,
    format='[%(levelname)s] %(message)s',
    stream=sys.stdout
)

SWAVE_MODEL_PATH = "/app/src/swave/predict/LSTM-l1-fc64-bi.pth"


def main():
    parser = argparse.ArgumentParser(description="Reads the projections from a .pkl file, runs the trained Bi-LSTM modelto predict"
                                     "structural variant labels for each snarl, and saves the predictions into a compressed .pkl file.")
    
    parser.add_argument("--projections_pkl", required=True, help="Path to the gzipped projections .pkl file.")
    parser.add_argument("--predictions_out_prefix", required=True, help="Output prefix for the predictions .pkl file.")
    
    # Optional parameters
    parser.add_argument("--model", default=None, help="Path to a custom trained model .pth file. If None, official pretrained model is used.")
    parser.add_argument("--device", default="cpu", choices=["cpu", "gpu"], help="Device to run prediction on (cpu or gpu).")
    parser.add_argument("--cpu_threads", type=int, default=4, help="Number of CPU threads for PyTorch (only used if device is cpu).")
    
    options = parser.parse_args()
    
    if options.model is None:
        model_path = SWAVE_MODEL_PATH
        logging.info(f"No custom model provided. Using official pretrained model from: {model_path}")
    else:
        model_path = options.model
        logging.info(f"Using model from: {model_path}")
    
    if not os.path.exists(model_path):
        logging.error(f"Model file not found: {model_path}")
        sys.exit(1)
    
    process_projections_to_predictions(
            projections_pickle_path=options.projections_pkl,
            model_path=model_path,
            predictions_output_prefix=options.predictions_out_prefix,
            device=options.device,
            num_threads=options.cpu_threads
        )

if __name__ == "__main__":
    main()
