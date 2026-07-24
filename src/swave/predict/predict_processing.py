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
import torch
import torch.utils.data as data
import pickle
import gzip
import logging

from .model import DecoderLSTM, Signal_DataSet, Signal_DataCollection


logging.basicConfig(
    level=logging.INFO,
    format='[%(levelname)s] %(message)s',
    stream=sys.stdout
)


def process_projections_to_predictions(projections_pickle_path, model_path, predictions_output_prefix, device="cpu", num_threads=4):
    """
    Loads precomputed projection matrices from a gzipped pickle file, runs the trained Bi-LSTM model to predict
    structural variant labels for each snarl and writes the predicted labels to a new gzipped pickle file.
    
    Args:
        projections_pickle_path (str): Path to the gzipped pickle file containing the projection matrices
        model_path (str): Path to the trained Bi-LSTM model file
        predictions_output_prefix (str): Prefix for the output gzipped pickle file containing predicted labels
        device (str): Device to run the model on ("cpu" or "cuda")
        num_threads (int): Number of threads to use for data loading and processing
    """
    
    if device == "gpu" and torch.cuda.is_available():
        actual_device = torch.device("cuda")
    else:
        actual_device = torch.device("cpu")
        torch.set_num_threads(num_threads)
    
    logging.info(f"Loading projection matrices from: {projections_pickle_path}")
    
    with gzip.open(projections_pickle_path, 'rb') as pickle_file:
        snarl_projections_dict = pickle.load(pickle_file)
    
    predict_dataset = Signal_DataSet(snarl_projections_dict)
    
    snarl_prediction_res = {}
    
    if len(predict_dataset) == 0:
        logging.warning("No valid snarl projections found for prediction. Output will be empty.")
    else:
        predict_loader = data.DataLoader(
            predict_dataset, 
            batch_size=128, 
            shuffle=False, 
            collate_fn=Signal_DataCollection, 
            pin_memory=True
        )
        
        logging.info(f"Initializing Bi-LSTM network on device: {actual_device}")
        
        model = DecoderLSTM(
            n_input=predict_dataset.data_dim,
            n_hidden=64,
            n_layer=1,
            architecture="LSTM",
            bidirect=True
        ).to(actual_device)
        
        model.load_state_dict(torch.load(model_path, map_location=actual_device))
        model.eval()
        
        logging.info("Starting prediction on projection matrices...")
        
        with torch.no_grad():
            for ids, Xs, _, X_lens in predict_loader:
                Xs = Xs.to(actual_device)
                
                pred = model(Xs, X_lens)
                pred_labels = torch.max(pred, 2)[1].cpu().data.numpy()
                
                for index in range(len(ids)):
                    snarl_prediction_res[ids[index]] = pred_labels[index][:X_lens[index]]
        
    output_pickle_path = predictions_output_prefix + "_predictions.pkl.gz"
    
    with gzip.open(output_pickle_path, 'wb') as output_file:
        pickle.dump(snarl_prediction_res, output_file, protocol=pickle.HIGHEST_PROTOCOL)
    
    logging.info(f"{len(snarl_prediction_res)} model predictions completed and saved.")
