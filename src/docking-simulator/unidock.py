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


@serve.deployment(num_replicas=2, ray_actor_options={"num_cpus": 0.2, "num_gpus": 0})
@serve.ingress(app)
class UniDock:
    def __init__(self, data_dir="/app/data"):
        self.data_dir = data_dir
        self.selfies_path = os.path.join(data_dir, "selfies.txt")
        self.ligand_dir = os.path.join(data_dir, "ligands/")
        self.target_pdb_path = os.path.join(data_dir, "target.pdb")
        self.target_pdbqt_path = os.path.join(data_dir, "target.pdbqt")
        self.output_dir = os.path.join(data_dir, "results/")

    @app.get("/")
    def hello_world(self):
        return {"Hello": "world!"}

    @app.post("/")
    def get_docking_scores(self, request: DockingRequest) -> DockingResponse:
        self.retrieve_file_from_s3(
            DockingRequest.bucket,
            DockingRequest.selfies_object
            )
        self.retrieve_file_from_s3(
            DockingRequest.bucket,
            DockingRequest.target_object
            )
        self.convert_selfies_to_sdf_files(self.selfies_filepath, self.ligand_dir)
        self.convert_pdb_to_pdbqt(self.target_pdb_path, self.target_pdbqt_path)
        self.dock_molecules(
            self.ligand_dir,
            self.target_pdbqt_path,
            DockingRequest.target_box,
            self.output_dir
            )
        return {"Response": "test"}

    def retrieve_file_from_s3(self, bucket_name, object_name):
        # Create an S3 client that doesn't require credentials
        s3 = boto3.client('s3', config=Config(signature_version=UNSIGNED))

        # Download file
        s3.download_file(bucket_name, object_name, self.selfies_path)

    def convert_selfies_to_sdf_files(self):
        with open(self.selfies_path, 'r') as file:
            selfies = file.readlines()
        
        # Remove any whitespace and newline characters
        selfies = [selfie.strip() for selfie in selfies]

        # Save molecule as SDF file
        os.makedirs(self.ligand_dir)
        
        for i, selfie in enumerate(selfies):
            # Convert SELFIE to SMILES
            smiles = sf.decoder(selfie)
            
            # Convert SMILES to RDKit molecule objects
            mol = Chem.MolFromSmiles(smiles)

            # Generate 3D coordinates
            AllChem.EmbedMolecule(mol)
            AllChem.UFFOptimizeMolecule(mol)
            
            # Save sdf file
            sdf_path = os.path.join(self.ligand_dir, f'mol_{i+1}.sdf')
            writer = SDWriter(sdf_path)
            writer.write(mol)
            writer.close()

    def convert_pdb_to_pdbqt(self):
        subprocess.run(["unidocktools", "proteinprep", "-r", self.target_pdb_path, "-o", self.target_pdbqt_path])

    def dock_molecules(self, bounding_box: TargetBox):
        # List ligand files
        ligand_filepaths = os.listdir(self.ligand_dir)
        
        # Run unidock
        subprocess.run([
            "unidocktools", "unidock_pipeline",
            "-r", self.target_pdbqt_path,
            "-l", *ligand_filepaths,
            "-sd", self.output_dir,
            "-cx", bounding_box.center_x, "-cy", bounding_box.center_y, "-cz", bounding_box.center_z
            ])

docking_app = UniDock.bind()
