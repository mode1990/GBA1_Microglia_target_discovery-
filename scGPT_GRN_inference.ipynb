"""
Complete scGPT Gene Regulatory Network (GRN) Analysis Pipeline
This script performs end-to-end GRN analysis using pre-trained scGPT models
"""

import os
import sys
import json
import torch
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from anndata import AnnData
import scanpy as sc
import gdown

# Import scGPT modules (assumes torchtext patch is already applied)
sys.path.insert(0, "../")
import scgpt as scg
from scgpt.tasks import GeneEmbedding
from scgpt.tokenizer.gene_tokenizer import GeneVocab
from scgpt.model import TransformerModel
from scgpt.preprocess import Preprocessor
from scgpt.utils import set_seed


class scGPT_GRN_Pipeline:
    """
    Complete pipeline for Gene Regulatory Network analysis using scGPT
    """
    
    def __init__(self, working_dir, model_dir, device='cpu'):
        """
        Initialize the GRN analysis pipeline
        
        Parameters:
        -----------
        working_dir : str
            Path to working directory
        model_dir : str or Path
            Path to pre-trained scGPT model directory
        device : str
            Device to run model on ('cpu' or 'cuda')
        """
        self.working_dir = Path(working_dir)
        self.model_dir = Path(model_dir)
        self.device = torch.device(device)
        
        # Set working directory
        os.chdir(self.working_dir)
        print(f"✓ Working directory set to: {os.getcwd()}")
        
        # Model components
        self.vocab = None
        self.model = None
        self.gene2idx = None
        self.gene_embeddings = None
        self.embed = None
        
        # Special tokens for scGPT
        self.special_tokens = ["<pad>", "<cls>", "<eoc>"]
        self.pad_token = "<pad>"
        self.pad_value = 0
        self.n_input_bins = 51
        
    def download_model(self, gdrive_url):
        """
        Download pre-trained model from Google Drive
        
        Parameters:
        -----------
        gdrive_url : str
            Google Drive folder URL containing the model
        """
        print("Downloading model files from Google Drive...")
        gdown.download_folder(gdrive_url, quiet=False, use_cookies=False)
        print("✓ Model files downloaded successfully")
    
    def load_vocabulary(self):
        """Load and prepare vocabulary from model directory"""
        print("Loading vocabulary...")
        vocab_file = self.model_dir / "vocab.json"
        
        self.vocab = GeneVocab.from_file(vocab_file)
        
        # Add special tokens if not present
        for token in self.special_tokens:
            if token not in self.vocab:
                self.vocab.append_token(token)
        
        self.gene2idx = self.vocab.get_stoi()
        print(f"✓ Vocabulary loaded: {len(self.vocab)} genes")
        
    def load_model(self):
        """Load pre-trained scGPT model with configuration"""
        print("Loading pre-trained scGPT model...")
        
        model_config_file = self.model_dir / "args.json"
        model_file = self.model_dir / "best_model.pt"
        
        # Load model configuration
        with open(model_config_file, "r") as f:
            model_configs = json.load(f)
        
        print(f"Model configuration loaded from {model_config_file}")
        
        # Extract model parameters
        embsize = model_configs["embsize"]
        nhead = model_configs["nheads"]
        d_hid = model_configs["d_hid"]
        nlayers = model_configs["nlayers"]
        n_layers_cls = model_configs["n_layers_cls"]
        
        # Initialize model
        ntokens = len(self.vocab)
        self.model = TransformerModel(
            ntokens,
            embsize,
            nhead,
            d_hid,
            nlayers,
            vocab=self.vocab,
            pad_value=self.pad_value,
            n_input_bins=self.n_input_bins,
        )
        
        # Load model weights
        try:
            self.model.load_state_dict(torch.load(model_file, map_location=self.device))
            print(f"✓ Full model loaded from {model_file}")
        except RuntimeError as e:
            print(f"Full load failed, attempting partial load: {e}")
            model_dict = self.model.state_dict()
            pretrained_dict = torch.load(model_file, map_location=self.device)
            pretrained_dict = {
                k: v for k, v in pretrained_dict.items()
                if k in model_dict and v.shape == model_dict[k].shape
            }
            for k, v in pretrained_dict.items():
                print(f"  Loading params {k} with shape {v.shape}")
            model_dict.update(pretrained_dict)
            self.model.load_state_dict(model_dict)
            print("✓ Partial model loaded successfully")
        
        # Move model to device
        self.model.to(self.device)
        self.model.eval()
        print(f"✓ Model moved to {self.device}")
    
    def extract_gene_embeddings(self, adata):
        """
        Extract gene embeddings from the pre-trained model
        
        Parameters:
        -----------
        adata : AnnData
            Annotated data object containing gene expression data
        
        Returns:
        --------
        dict : Gene embeddings filtered by genes in adata
        """
        print("Extracting gene embeddings...")
        
        gene_ids = np.array([idx for idx in self.gene2idx.values()])
        
        with torch.no_grad():
            embeddings = self.model.encoder(
                torch.tensor(gene_ids, dtype=torch.long).to(self.device)
            )
            embeddings = embeddings.detach().cpu().numpy()
        
        # Filter embeddings for genes present in the dataset
        hvg_genes = set(adata.var.index.tolist())
        self.gene_embeddings = {
            gene: embeddings[i] 
            for i, gene in enumerate(self.gene2idx.keys()) 
            if gene in hvg_genes
        }
        
        print(f"✓ Gene embeddings extracted for {len(self.gene_embeddings)} genes")
        return self.gene_embeddings
    
    def identify_gene_programs(self, resolution=10, min_genes=5):
        """
        Identify gene programs using Louvain clustering
        
        Parameters:
        -----------
        resolution : float
            Louvain clustering resolution parameter
        min_genes : int
            Minimum number of genes required per program
        
        Returns:
        --------
        dict : Gene programs (metagenes) with their constituent genes
        """
        print(f"Identifying gene programs (resolution={resolution})...")
        
        # Create GeneEmbedding object
        self.embed = GeneEmbedding(self.gene_embeddings)
        
        # Perform Louvain clustering
        gdata = self.embed.get_adata(resolution=resolution)
        metagenes = self.embed.get_metagenes(gdata)
        
        # Filter programs by minimum gene count
        self.gene_programs = {
            mg: genes for mg, genes in metagenes.items() 
            if len(genes) >= min_genes
        }
        
        print(f"✓ Identified {len(self.gene_programs)} gene programs")
        for mg, genes in list(self.gene_programs.items())[:5]:
            print(f"  {mg}: {len(genes)} genes")
        
        return self.gene_programs
    
    def score_gene_programs(self, adata):
        """
        Score gene programs in the dataset
        
        Parameters:
        -----------
        adata : AnnData
            Annotated data object
        """
        print("Scoring gene programs in dataset...")
        self.embed.score_metagenes(adata, self.gene_programs)
        print(f"✓ Gene program scores added to adata.obs")
    
    def visualize_gene_programs(self, adata, group_by, output_file=None, dpi=600):
        """
        Visualize gene program activation across groups
        
        Parameters:
        -----------
        adata : AnnData
            Annotated data object with scored gene programs
        group_by : str
            Column name in adata.obs to group by (e.g., 'Genotype', 'celltype')
        output_file : str, optional
            Path to save the figure (PDF format recommended)
        dpi : int
            Resolution for saved figure
        """
        print(f"Visualizing gene programs by {group_by}...")
        
        sns.set(font_scale=0.35)
        
        # Create visualization
        self.embed.plot_metagenes_scores(adata, self.gene_programs, group_by)
        plt.tight_layout()
        
        # Save figure if output path provided
        if output_file:
            plt.savefig(output_file, format='pdf', dpi=dpi, bbox_inches='tight')
            print(f"✓ Figure saved to {output_file}")
        
        plt.show()
    
    def run_complete_pipeline(self, adata, group_by='Genotype', 
                            resolution=10, min_genes=5, 
                            output_file='gene_programs.pdf'):
        """
        Run the complete GRN analysis pipeline
        
        Parameters:
        -----------
        adata : AnnData
            Annotated data object containing gene expression data
        group_by : str
            Column name for grouping in visualization
        resolution : float
            Louvain clustering resolution
        min_genes : int
            Minimum genes per program
        output_file : str
            Output file path for visualization
        """
        print("\n" + "="*60)
        print("Starting Complete scGPT GRN Analysis Pipeline")
        print("="*60 + "\n")
        
        # Step 1: Load vocabulary
        self.load_vocabulary()
        
        # Step 2: Load model
        self.load_model()
        
        # Step 3: Extract gene embeddings
        self.extract_gene_embeddings(adata)
        
        # Step 4: Identify gene programs
        self.identify_gene_programs(resolution=resolution, min_genes=min_genes)
        
        # Step 5: Score gene programs
        self.score_gene_programs(adata)
        
        # Step 6: Visualize results
        self.visualize_gene_programs(adata, group_by=group_by, output_file=output_file)
        
        print("\n" + "="*60)
        print("Pipeline completed successfully!")
        print("="*60 + "\n")
        
        return self.gene_programs


# =============================================================================
# Example Usage
# =============================================================================

if __name__ == "__main__":
    
    # Configuration
    WORKING_DIR = "/data/outputs/cellranger_outputs/Analysis_by_Mo_Multi/"
    MODEL_DIR = "/data/outputs/cellranger_outputs/Analysis_by_Mo_Multi/scGPT_brain"
    GDRIVE_URL = 'https://drive.google.com/drive/folders/1vf1ijfQSk7rGdDGpBntR5bi5g6gNt-Gx'
    
    # Initialize pipeline
    pipeline = scGPT_GRN_Pipeline(
        working_dir=WORKING_DIR,
        model_dir=MODEL_DIR,
        device='cpu'  # Use 'cuda' if GPU available
    )
    
    # Optional: Download model if needed
    # pipeline.download_model(GDRIVE_URL)
    
    # Load your data (example - replace with your actual data loading)
    # adata = sc.read_h5ad("your_microglia_data.h5ad")
    
    # Run complete pipeline
    # gene_programs = pipeline.run_complete_pipeline(
    #     adata=adata,
    #     group_by='Genotype',  # or 'celltype', 'condition', etc.
    #     resolution=10,
    #     min_genes=5,
    #     output_file='gene_programs_microglia.pdf'
    # )
    
    print("Pipeline class initialized. Load your AnnData and run:")
    print("gene_programs = pipeline.run_complete_pipeline(adata)")
