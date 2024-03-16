from typing import List
import selfies as sf

def smiles_to_selfies(smiles_list: List[str]):
    """
    Convert a list of SMILES to a list of SELFIES
    
    Parameters:
    - smiles_list: List of SMILES strings.
    
    Returns:
    - List of SELFIES strings.
    """
    selfies_list = [sf.encoder(smiles) for smiles in smiles_list]
    return selfies_list


if __name__=="__main__":
    with open('/home/convexity-research/seth/docking-simulator/small-molecules.smi', 'r') as file:
        # Read in smiles and identifiers
        lines = file.readlines()
        # Remove whitespace
        lines = [line.strip() for line in lines]
        # Remove identifiers
        smiles = [line.split()[0] for line in lines]

    # Convert list of smiles to list of selfies
    selfies = smiles_to_selfies(smiles)
    # Write one selfie per line
    with open('/home/convexity-research/seth/docking-simulator/small-molecules.txt', 'w') as file:
        for selfie in selfies:
            file.write(selfie + '\n')
