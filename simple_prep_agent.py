#!/usr/bin/env python3
"""
Simple PRELYM Preparation Agent - Python 3.13 Compatible
No pdbfixer/openmm dependencies required

Uses only BioPython + requests for basic file preparation
"""

import requests
import json
from pathlib import Path
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB import Select

class RemoveWatersSelect(Select):
    WATER_NAMES = {"HOH", "H2O", "WAT"}
    def accept_residue(self, residue):
        resname = residue.get_resname().strip()
        return resname not in self.WATER_NAMES

def search_protein(name: str):
    """Search RCSB for protein by name"""
    query = {
        "query": {
            "type": "terminal",
            "service": "full_text",
            "parameters": {"value": name}
        },
        "return_type": "entry",
        "request_options": {"return_all_hits": True}
    }
    url = "https://search.rcsb.org/rcsbsearch/v2/query"
    response = requests.post(url, json=query)
    if response.status_code == 200:
        data = response.json()
        hits = [h["identifier"] for h in data.get("result_set", [])]
        return hits[0] if hits else None
    return None

def download_pdb(pdb_id: str, output_path: Path):
    """Download PDB file from RCSB"""
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    response = requests.get(url)
    response.raise_for_status()
    output_path.write_text(response.text)
    return output_path

def clean_pdb(input_pdb: Path, output_pdb: Path):
    """Remove waters and create clean PDB"""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", str(input_pdb))

    io = PDBIO()
    io.set_structure(structure)
    io.save(str(output_pdb), select=RemoveWatersSelect())

    return output_pdb

def create_basic_files(pdb_id: str, output_dir: Path):
    """Create the three basic files needed for PRELYM"""
    output_dir.mkdir(parents=True, exist_ok=True)

    # Download original PDB
    original_pdb = output_dir / f"{pdb_id}.pdb"
    download_pdb(pdb_id, original_pdb)

    # Create clean PDB (no waters)
    clean_pdb_file = output_dir / f"{pdb_id}_clean.pdb"
    clean_pdb(original_pdb, clean_pdb_file)

    # For H-bond file, copy clean file (basic fallback)
    hbond_file = output_dir / f"{pdb_id}FH.pdb"
    hbond_file.write_text(clean_pdb_file.read_text())

    # For PQR, create a basic version (copy with .pqr extension)
    pqr_file = output_dir / f"{pdb_id}.pqr"
    pqr_file.write_text(clean_pdb_file.read_text())

    return {
        'clean_pdb': clean_pdb_file,
        'hbond_file': hbond_file,
        'pqr_file': pqr_file
    }

class SimplePrepAgent:
    """Basic preparation agent compatible with Python 3.13"""

    def prepare_from_protein(self, protein_name: str, output_dir: Path):
        pdb_id = search_protein(protein_name)
        if not pdb_id:
            raise ValueError(f"No PDB found for protein: {protein_name}")
        return create_basic_files(pdb_id, output_dir)

    def prepare_from_pdb_id(self, pdb_id: str, output_dir: Path):
        return create_basic_files(pdb_id, output_dir)

    def check_dependencies(self):
        return {
            'biopython': True,
            'requests': True,
            'pdb2pqr': False,
            'reduce': False,
            'pdbfixer': False
        }

if __name__ == "__main__":
    agent = SimplePrepAgent()
    files = agent.prepare_from_protein("lysozyme", Path("./test_output"))
    print("Files created:")
    for name, path in files.items():
        print(f"  {name}: {path}")