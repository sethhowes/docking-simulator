import os
import subprocess

from fastapi import FastAPI
from ray import serve

import boto3
from botocore import UNSIGNED
from botocore.config import Config

from rdkit import Chem
from rdkit.Chem import AllChem, SDWriter
import selfies as sf

from data import DockingRequest, DockingResponse, TargetBox

app = FastAPI()


@serve.deployment(num_replicas=2, ray_actor_options={"num_gpus": 1})
@serve.ingress(app)
class UniDock:
    def __init__(self, data_dir="/app/data"):
        self.selfies_path = os.path.join(data_dir, "selfies.txt")
        self.ligand_dir = os.path.join(data_dir, "ligands/")
        self.target_pdb_path = os.path.join(data_dir, "target.pdb")
        self.target_pdbqt_path = os.path.join(data_dir, "target.pdbqt")
        self.output_dir = os.path.join(data_dir, "results/")

    @app.post("/",)
    def run_pipeline(self, request: DockingRequest) -> DockingResponse:
        retrieve_file_from_s3(
            DockingRequest.bucket,
            DockingRequest.selfies_object,
            os.path.join(self.data_dir, "selfies.txt")
            )
        retrieve_file_from_s3(
            DockingRequest.bucket,
            DockingRequest.target_object.
            os.path.join(self.data_dir, "target.pdb")
            )
        convert_selfies_to_sdf_files(self.selfies_filepath, self.ligand_dir)
        convert_pdb_to_pdbqt(self.target_pdb_path, self.target_pdbqt_path)
        dock_molecules(
            self.ligand_dir,
            self.target_pdbqt_path,
            DockingRequest.target_box,
            self.output_dir
            )
        


def retrieve_file_from_s3(bucket_name, object_name, destination_path):
    # Create an S3 client that doesn't require credentials
    s3 = boto3.client('s3', config=Config(signature_version=UNSIGNED))

    # The destination data path in the Docker container
    local_file_path = os.path.join(destination_path, object_name)

    # Download file
    s3.download_file(bucket_name, object_name, destination_path)


def convert_selfies_to_sdf_files(selfies_filepath: str, output_directory: str):
    with open(selfies_filepath, 'r') as file:
        selfies = file.readlines()
    
    # Remove any whitespace and newline characters
    selfies = [selfie.strip() for selfie in selfies]

    # Save molecule as SDF file
    os.makedirs(output_directory)
    
    for i, selfie in enumerate(selfies):
        # Convert SELFIE to SMILES
        smiles = sf.decoder(selfie)
        
        # Convert SMILES to RDKit molecule objects
        mol = Chem.MolFromSmiles(smiles)

        # Generate 3D coordinates
        AllChem.EmbedMolecule(mol)
        AllChem.UFFOptimizeMolecule(mol)
        
        # Save sdf file
        sdf_path = os.path.join(output_directory, f'mol_{i+1}.sdf')
        writer = SDWriter(sdf_path)
        writer.write(mol)
        writer.close()


def convert_pdb_to_pdbqt(input_filepath: str, output_filepath: str):
    subprocess.run(["unidocktools", "proteinprep", "-r", input_filepath, "-o", output_filepath])


def dock_molecules(ligand_directory: str, target_filepath: str, bounding_box: TargetBox, output_directory: str):
    # List ligand files
    ligand_filepaths = os.listdir(ligand_directory)
    
    # Run unidock
    subprocess.run([
        "unidocktools", "unidock_pipeline",
        "-r", target_filepath,
        "-l", *ligand_filepaths,
        "-sd", output_directory,
        "-cx", bounding_box.center_x, "-cy", bounding_box.center_y, "-cz", bounding_box.center_z
        ])
