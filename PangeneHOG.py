#!/usr/bin/env python3
import os
import subprocess
import argparse
import pandas as pd
from typing import List, Dict, Optional
from Bio import SeqIO
import logging

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

class PangenePipeline:
    def __init__(self, config: Optional[Dict] = None):
        """
        Initialize the pangene pipeline.
        
        Args:
            config: Dictionary containing configuration parameters for clustering and tools.
        """
        self.working_dir = './pangene_results'
        self.cluster_dir = os.path.join(self.working_dir, 'clustering')
        os.makedirs(self.working_dir, exist_ok=True)
        os.makedirs(self.cluster_dir, exist_ok=True)
        
        # Default configuration
        self.config = {
            'diamond_path': 'diamond',
            'seqtk_path': 'seqtk',
            'miniprot_cmd': ['miniprot', '--outs=0.97', '--no-cs', '-Iut16'],
            'pangene_cmd': ['pangene'],
            'clustering': {
                'identity': 90,
                'coverage': 80,
                'threads': 8,
                'memory': '64G'
            }
        }
        
        if config:
            self.config.update(config)
    
    def clean_fasta_nulls_and_stars(self, input_fasta: str, output_fasta: str) -> None:
        """
        Remove null characters, trailing '*' from sequence lines, and Windows CRLF endings.
        
        Args:
            input_fasta: Path to input FASTA file
            output_fasta: Path to write cleaned output FASTA
        """
        logger.info(f"Cleaning FASTA: {input_fasta}")
        
        with open(input_fasta, 'r') as f_in, open(output_fasta, 'w') as f_out:
            for line in f_in:
                if line.startswith('>'):
                    # Header line - write as is
                    f_out.write(line)
                else:
                    # Sequence line - remove nulls, trailing *, and normalize line endings
                    line = line.replace('\0', '').rstrip('*\r\n') + '\n'
                    f_out.write(line)
        
        logger.info(f"Cleaned FASTA written to: {output_fasta}")
    
    def filter_amino_acids(self, input_fasta: str, output_fasta: str) -> None:
        """
        Filter out non-standard amino acids and whitespace from sequence lines.
        
        Args:
            input_fasta: Path to input FASTA file
            output_fasta: Path to write filtered output FASTA
        """
        logger.info(f"Filtering non-standard amino acids in: {input_fasta}")
        
        valid_aa = set('ACDEFGHIKLMNPQRSTVWY')
        
        with open(input_fasta, 'r') as f_in, open(output_fasta, 'w') as f_out:
            for line in f_in:
                if line.startswith('>'):
                    # Header line - write as is
                    f_out.write(line)
                else:
                    # Filter sequence line to only contain valid amino acids
                    filtered_seq = ''.join(c for c in line.upper() if c in valid_aa)
                    f_out.write(filtered_seq + '\n')
        
        logger.info(f"Filtered FASTA written to: {output_fasta}")
    
    def _run_command(self, cmd: List[str], output_file: Optional[str] = None, **kwargs) -> bool:
        """
        Run a shell command with error handling.
        
        Args:
            cmd: Command to run as a list of strings
            output_file: Optional file path to write command output
            **kwargs: Additional arguments to subprocess.run()
            
        Returns:
            bool: True if command succeeded, False otherwise
        """
        cmd_str = ' '.join(str(c) for c in cmd)
        logger.info(f"Running command: {cmd_str}")
        
        try:
            if output_file:
                with open(output_file, 'w') as f:
                    result = subprocess.run(
                        cmd, 
                        stdout=f, 
                        stderr=subprocess.PIPE,
                        text=True,
                        check=False,
                        **kwargs
                    )
            else:
                result = subprocess.run(
                    cmd, 
                    capture_output=True, 
                    text=True,
                    check=False,
                    **kwargs
                )
                if result.stderr:
                    logger.warning(f"Command stderr: {result.stderr}")
            
            if result.returncode != 0:
                error_msg = result.stderr if result.stderr else "No error message"
                logger.error(f"Command failed with exit code {result.returncode}")
                logger.error(f"Error output: {error_msg}")
                return False
                
            return True
            
        except Exception as e:
            logger.error(f"Error running command: {e}")
            if hasattr(e, 'stderr') and e.stderr:
                if isinstance(e.stderr, bytes):
                    logger.error(f"Error output: {e.stderr.decode('utf-8', errors='replace')}")
                else:
                    logger.error(f"Error output: {e.stderr}")
            return False
    
    def cluster_proteins(self, input_fasta: str, skip_filtering: bool = False) -> Optional[str]:
        """
        Cluster proteins using DIAMOND linclust.
        
        Args:
            input_fasta: Path to input FASTA file
            skip_filtering: If True, skip the FASTA cleaning and AA filtering steps
            
        Returns:
            Path to clustered FASTA file or None if clustering failed
        """
        logger.info("Starting protein clustering...")
        
        if skip_filtering:
            logger.info("Skipping FASTA cleaning and filtering as requested")
            clean2 = input_fasta  # Use the input file directly
        else:
            
            logger.info(f"Cleaning FASTA: {input_fasta}")
            clean1 = os.path.join(self.cluster_dir, 'pan_proteome.clean1.fa')
            self.clean_fasta_nulls_and_stars(input_fasta, clean1)
            
            
            logger.info(f"Filtering non-standard amino acids in: {clean1}")
            clean2 = os.path.join(self.cluster_dir, 'pan_proteome.clean2.fa')
            self.filter_amino_acids(clean1, clean2)
        
        
        cluster_tsv = os.path.join(self.cluster_dir, 'pan_lin_clusters.tsv')
        
        cmd = [
            self.config['diamond_path'], 'linclust',
            '-d', clean2,  # Input FASTA file
            '-o', cluster_tsv,  # Output TSV file
            '-M', str(self.config['clustering']['memory']),  # Memory limit
            '--approx-id', str(int(self.config['clustering']['identity'])),  # Approx identity (0-100)
            '--member-cover', str(int(self.config['clustering']['coverage'])),  # Member coverage (0-100)
            '--threads', str(self.config['clustering']['threads']),  # Number of threads
            '--tmpdir', self.cluster_dir  # Directory for temporary files
        ]
        if not self._run_command(cmd):
            logger.error("DIAMOND linclust failed")
            return None
        
        
        reps_file = os.path.join(self.cluster_dir, 'reps.ids')
        with open(cluster_tsv, 'r') as f_in, open(reps_file, 'w') as f_out:
            seen = set()
            for line in f_in:
                rep_id = line.split('\t')[0]
                if rep_id not in seen:
                    f_out.write(f"{rep_id}\n")
                    seen.add(rep_id)
        
        
        clustered_fasta = os.path.join(self.cluster_dir, 'pan.proteomic.clustered.faa')
        cmd = [
            self.config['seqtk_path'], 'subseq',
            clean2,
            reps_file
        ]
        if not self._run_command(cmd, clustered_fasta):
            logger.error("Failed to extract representative sequences")
            return None
        
        
        db_prefix = os.path.join(self.cluster_dir, 'repsDB')
        cmd = [
            self.config['diamond_path'], 'makedb',
            '--in', clustered_fasta,
            '-d', db_prefix
        ]
        if not self._run_command(cmd):
            logger.error("Failed to create DIAMOND database")
            return None
        
        #All-vs-all BLASTP
        blastp_out = os.path.join(self.cluster_dir, 'reps_vs_reps.m8')
        
        cmd = [
            self.config['diamond_path'], 'blastp',
            '--db', f"{db_prefix}.dmnd",  # DIAMOND database
            '--query', clustered_fasta,  # Query sequences
            '--out', blastp_out,  # Output file
            '--outfmt', '6', 'qseqid', 'sseqid', 'corrected_bitscore',  # Output format
            '--approx-id', str(int(self.config['clustering']['identity'])),  # Approx identity (0-100)
            '--query-cover', str(int(self.config['clustering']['coverage'])),  # Query coverage (0-100)
            '--threads', str(self.config['clustering']['threads']),  # Number of threads
            '--tmpdir', self.cluster_dir  # Directory for temporary files
        ]
        if not self._run_command(cmd):
            logger.error("All-vs-all BLASTP failed")
            return None
        
        logger.info(f"Protein clustering completed. Clustered FASTA: {clustered_fasta}")
        return clustered_fasta

    def prepare_protein_set(self, pan_proteome_file):
        """
        Prepare protein set from pan_proteome.fa file.
        Converts sequences to pangene format.
        """
        output_file = os.path.join(self.working_dir, 'proteins.faa')
        
       
        with open(pan_proteome_file, 'r') as f_in, open(output_file, 'w') as f_out:
            for record in SeqIO.parse(f_in, 'fasta'):
                # Format sequence ID as GeneName:ProteinID
                gene_name = record.id.split('|')[0]  # Assuming gene name is before first |
                protein_id = f"{gene_name}:{record.id}"
                record.id = protein_id
                record.description = ''
                SeqIO.write(record, f_out, 'fasta')

        return output_file

    def run_miniprot(self, genome_file: str, protein_file: str) -> str:
        """
        Run miniprot alignment with configurable command.
        
        Args:
            genome_file: Path to genome FASTA file
            protein_file: Path to protein FASTA file
            
        Returns:
            Path to the output PAF file
        """
        output_file = os.path.join(self.working_dir, f"{os.path.basename(genome_file)}.paf")
        
        #configurable miniprot command and append input files
        cmd = self.config['miniprot_cmd'] + [genome_file, protein_file]
        
        if not self._run_command(cmd, output_file):
            raise RuntimeError(f"Miniprot failed for {genome_file}")
            
        return output_file

    def construct_pangene_graph(self, paf_files: List[str]) -> str:
        """
        Construct pangene graph from PAF files using configurable command.
        
        Args:
            paf_files: List of PAF file paths
            
        Returns:
            Path to the output GFA file
        """
        output_file = os.path.join(self.working_dir, 'graph.gfa')
        
        
        cmd = self.config['pangene_cmd'] + paf_files
        
        if not self._run_command(cmd, output_file):
            raise RuntimeError("Failed to construct pangene graph")
            
        return output_file

    def analyze_graph(self, gfa_file):
        """Analyze the pangene graph and generate gene presence/absence matrix."""
        bubble_file = os.path.join(self.working_dir, 'bubble.txt')
        matrix_file = os.path.join(self.working_dir, 'gene_presence_absence.Rtab')
        
        #pangene.js call to create bubble.txt
        cmd1 = ['pangene.js', 'call', gfa_file]
        try:
            with open(bubble_file, 'w') as f:
                subprocess.run(cmd1, stdout=f, check=True)
        except subprocess.CalledProcessError as e:
            logger.error(f"Error creating bubble.txt: {e}")
            raise
        
        #pangene.js for gfa2matrix
        cmd2 = ['pangene.js', 'gfa2matrix', gfa_file]
        try:
            with open(matrix_file, 'w') as f:
                subprocess.run(cmd2, stdout=f, check=True)
            return matrix_file
        except subprocess.CalledProcessError as e:
            logger.error(f"Error creating gene presence/absence matrix: {e}")
            raise

    def run_pipeline(self, pan_proteome_file: str, genome_files: List[str],
                    skip_clustering: bool = False, skip_filtering: bool = False) -> str:
        """
        Run the complete pangene pipeline with optional protein clustering.
        
        Args:
            pan_proteome_file: Path to pan_proteome.fa file from PanHOG
            genome_files: List of genome FASTA files
            skip_clustering: If True, skip the protein clustering step
            skip_filtering: If True, skip the FASTA cleaning and AA filtering steps
            
        Returns:
            Path to the final output file
        """
        logger.info("Starting pangene pipeline...")
        
        
        if not skip_clustering:
            logger.info("Running protein clustering...")
            clustered_proteome = self.cluster_proteins(pan_proteome_file, skip_filtering=skip_filtering)
            if not clustered_proteome:
                logger.warning("Protein clustering failed, using original proteome")
                clustered_proteome = pan_proteome_file
        else:
            logger.info("Skipping protein clustering as requested")
            clustered_proteome = pan_proteome_file
        
        
        protein_file = self.prepare_protein_set(clustered_proteome)
        logger.info(f"Protein set prepared: {protein_file}")

        
        paf_files = []
        for genome in genome_files:
            logger.info(f"Running miniprot for {os.path.basename(genome)}...")
            paf = self.run_miniprot(genome, protein_file)
            paf_files.append(paf)

        
        logger.info("Constructing pangene graph...")
        gfa_file = self.construct_pangene_graph(paf_files)

        
        logger.info("Analyzing pangene graph...")
        output_file = self.analyze_graph(gfa_file)
        logger.info(f"Analysis complete. Results saved to {output_file}")

        return output_file

def main():
    parser = argparse.ArgumentParser(description='Run pangene pipeline with PanHOG results')
    parser.add_argument('--pan-proteome', required=True, help='Path to pan_proteome.fa file')
    parser.add_argument('--genomes', nargs='+', required=True, help='List of genome FASTA files')
    
    # Clustering options
    parser.add_argument('--skip-clustering', action='store_true', 
                       help='Skip the protein clustering step')
    parser.add_argument('--skip-filtering', action='store_true',
                       help='Skip FASTA cleaning and amino acid filtering steps')
    parser.add_argument('--diamond-path', default='diamond',
                       help='Path to DIAMOND executable')
    parser.add_argument('--seqtk-path', default='seqtk',
                       help='Path to seqtk executable')
    parser.add_argument('--identity', type=float, default=90.0,
                       help='Minimum percent identity for clustering (0-100)')
    parser.add_argument('--coverage', type=float, default=80.0,
                       help='Minimum coverage for clustering (0-100)')
    parser.add_argument('--threads', type=int, default=8,
                       help='Number of threads to use')
    parser.add_argument('--memory', default='64G',
                       help='Memory limit for DIAMOND (e.g., 64G)')
    
    # Advanced: Custom command overrides
    parser.add_argument('--miniprot-cmd', help='Custom miniprot command (overrides other miniprot options)')
    parser.add_argument('--pangene-cmd', help='Custom pangene command (overrides other pangene options)')
    
    args = parser.parse_args()
    
    #prepare conf.
    config = {
        'diamond_path': args.diamond_path,
        'seqtk_path': args.seqtk_path,
        'clustering': {
            'identity': args.identity,
            'coverage': args.coverage,
            'threads': args.threads,
            'memory': args.memory
        }
    }
    
    # Add custom commands if provided
    if args.miniprot_cmd:
        config['miniprot_cmd'] = args.miniprot_cmd.split()
    if args.pangene_cmd:
        config['pangene_cmd'] = args.pangene_cmd.split()
    
    # initialize pipeline with config
    pipeline = PangenePipeline(config=config)
    
    # Run
    try:
        output_file = pipeline.run_pipeline(
            args.pan_proteome, 
            args.genomes,
            skip_clustering=args.skip_clustering,
            skip_filtering=args.skip_filtering
        )
        logger.info(f"Pipeline completed successfully. Results saved to {output_file}")
    except Exception as e:
        logger.error(f"Pipeline failed: {e}", exc_info=True)
        raise SystemExit(1)

if __name__ == "__main__":
    main()

