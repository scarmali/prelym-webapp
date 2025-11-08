#!/usr/bin/env python3
"""
PRELYM File Preparation Agent

Automated preparation of the three input files needed for PRELYM protein analysis:
1. Clean PDB structure file (*_clean.pdb)
2. Hydrogen-optimized structure file (*FH.pdb)
3. Charge file with partial charges (*.pqr)

Usage:
    python prelym_prep_agent.py --protein "lysozyme"
    python prelym_prep_agent.py --pdb-id 4CHA
    python prelym_prep_agent.py --pdb-file myprotein.pdb

Requirements:
    - reduce (MolProbity suite) - for optimal hydrogen placement
    - pdb2pqr - for charge assignment
    Both available via conda-forge or direct installation

Author: PRELYM Prep Agent
Version: 1.0
"""

import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
import urllib.request
import urllib.parse
from pathlib import Path
from typing import Optional, Tuple, List, Dict

try:
    from Bio.PDB import PDBParser, PDBIO, Select
    from Bio.PDB.PDBExceptions import PDBConstructionWarning
    import warnings
    warnings.simplefilter('ignore', PDBConstructionWarning)
    HAS_BIOPYTHON = True
except ImportError:
    HAS_BIOPYTHON = False

try:
    from pdbfixer import PDBFixer
    from openmm.app import PDBFile
    HAS_PDBFIXER = True
except ImportError:
    HAS_PDBFIXER = False


class NoWaterSelect(Select):
    """BioPython selector to exclude water molecules (HOH)"""
    def accept_residue(self, residue):
        return residue.get_resname() != "HOH"


class PrelymPrepAgent:
    """Main agent class for PRELYM file preparation"""

    def __init__(self, output_dir: str = ".", verbose: bool = True):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        self.verbose = verbose
        self.temp_dir = None

    def log(self, message: str):
        """Print message if verbose mode enabled"""
        if self.verbose:
            print(f"[PRELYM Prep] {message}")

    def error(self, message: str):
        """Print error message and exit"""
        print(f"[ERROR] {message}", file=sys.stderr)
        sys.exit(1)

    def check_dependencies(self):
        """Check if required external tools are available"""
        tools = {}

        # Check for reduce
        try:
            result = subprocess.run(['reduce', '-version'],
                                  capture_output=True, text=True, timeout=10)
            tools['reduce'] = True
            self.log("✓ reduce found")
        except (subprocess.TimeoutExpired, subprocess.CalledProcessError, FileNotFoundError):
            tools['reduce'] = False
            self.log("✗ reduce not found - will use PDBFixer fallback")

        # Check for pdb2pqr
        try:
            result = subprocess.run(['pdb2pqr30', '--help'],
                                  capture_output=True, text=True, timeout=10)
            tools['pdb2pqr'] = 'pdb2pqr30'
            self.log("✓ pdb2pqr30 found")
        except (subprocess.TimeoutExpired, subprocess.CalledProcessError, FileNotFoundError):
            try:
                result = subprocess.run(['pdb2pqr', '--help'],
                                      capture_output=True, text=True, timeout=10)
                tools['pdb2pqr'] = 'pdb2pqr'
                self.log("✓ pdb2pqr found")
            except (subprocess.TimeoutExpired, subprocess.CalledProcessError, FileNotFoundError):
                tools['pdb2pqr'] = False
                self.log("✗ pdb2pqr not found")

        if not tools['pdb2pqr']:
            self.error("pdb2pqr is required but not found. Please install via conda or from pdb2pqr.readthedocs.io")

        return tools

    def search_protein(self, protein_name: str) -> Optional[str]:
        """Search RCSB for protein by name and return best PDB ID"""
        self.log(f"Searching RCSB for '{protein_name}'...")

        # RCSB search API query
        search_query = {
            "query": {
                "type": "terminal",
                "service": "full_text",
                "parameters": {
                    "value": protein_name
                }
            },
            "return_type": "entry",
            "request_options": {
                "pager": {"start": 0, "rows": 20}
            }
        }

        try:
            # Search for entries
            search_url = "https://search.rcsb.org/rcsbsearch/v2/query"
            search_data = json.dumps(search_query).encode('utf-8')

            req = urllib.request.Request(search_url, data=search_data,
                                       headers={'Content-Type': 'application/json'})

            with urllib.request.urlopen(req, timeout=30) as response:
                search_results = json.loads(response.read().decode())

            if not search_results.get('result_set'):
                self.log(f"No results found for '{protein_name}'")
                return None

            pdb_ids = [hit['identifier'] for hit in search_results['result_set']]
            self.log(f"Found {len(pdb_ids)} potential matches: {', '.join(pdb_ids[:5])}...")

            # Get details for ranking
            best_pdb = self._rank_pdb_entries(pdb_ids[:10])  # Check top 10

            if best_pdb:
                self.log(f"Selected best structure: {best_pdb}")
                return best_pdb
            else:
                self.log("No suitable X-ray structures found")
                return pdb_ids[0] if pdb_ids else None

        except Exception as e:
            self.log(f"Search failed: {e}")
            return None

    def _rank_pdb_entries(self, pdb_ids: List[str]) -> Optional[str]:
        """Rank PDB entries by resolution, preferring X-ray structures"""
        if not pdb_ids:
            return None

        best_entry = None
        best_resolution = float('inf')

        for pdb_id in pdb_ids:
            try:
                # Get entry details
                info_url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
                with urllib.request.urlopen(info_url, timeout=10) as response:
                    entry_data = json.loads(response.read().decode())

                # Check if X-ray crystallography
                exp_methods = entry_data.get('exptl', [])
                is_xray = any('X-RAY' in method.get('method', '').upper()
                            for method in exp_methods)

                if not is_xray:
                    continue

                # Get resolution
                resolution = entry_data.get('rcsb_entry_info', {}).get('resolution_combined', [None])
                if resolution and resolution[0] is not None:
                    res_value = float(resolution[0])
                    if res_value < best_resolution:
                        best_resolution = res_value
                        best_entry = pdb_id

            except Exception as e:
                self.log(f"Could not get info for {pdb_id}: {e}")
                continue

        return best_entry

    def download_pdb(self, pdb_id: str) -> str:
        """Download PDB file from RCSB"""
        pdb_id = pdb_id.upper()
        self.log(f"Downloading PDB {pdb_id}...")

        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        output_file = self.output_dir / f"{pdb_id}.pdb"

        try:
            urllib.request.urlretrieve(url, output_file)
            self.log(f"Downloaded {output_file}")
            return str(output_file)
        except Exception as e:
            self.error(f"Failed to download {pdb_id}: {e}")

    def clean_pdb(self, pdb_file: str) -> str:
        """Clean PDB file: remove waters, keep non-water heterogens"""
        self.log("Step 1: Cleaning PDB file...")

        input_path = Path(pdb_file)
        base_name = input_path.stem
        output_file = self.output_dir / f"{base_name}_clean.pdb"

        if not HAS_BIOPYTHON:
            self.error("BioPython is required for PDB cleaning. Install with: pip install biopython")

        try:
            # Parse structure
            parser = PDBParser(QUIET=True)
            structure = parser.get_structure("protein", pdb_file)

            # Write cleaned structure (no waters)
            io = PDBIO()
            io.set_structure(structure)
            io.save(str(output_file), NoWaterSelect())

            self.log(f"Cleaned PDB saved: {output_file}")
            return str(output_file)

        except Exception as e:
            self.error(f"PDB cleaning failed: {e}")

    def add_hydrogens(self, pdb_file: str, ph: float = 8.0, use_reduce: bool = True) -> str:
        """Add/optimize hydrogens using reduce or PDBFixer"""
        self.log("Step 2: Adding/optimizing hydrogens...")

        input_path = Path(pdb_file)
        base_name = input_path.stem.replace('_clean', '')
        output_file = self.output_dir / f"{base_name}FH.pdb"

        if use_reduce:
            return self._add_hydrogens_reduce(pdb_file, str(output_file))
        else:
            return self._add_hydrogens_pdbfixer(pdb_file, str(output_file), ph)

    def _add_hydrogens_reduce(self, input_file: str, output_file: str) -> str:
        """Add hydrogens using reduce (MolProbity suite)"""
        try:
            cmd = ['reduce', '-build', '-nuclear', input_file]

            with open(output_file, 'w') as f:
                result = subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE,
                                      text=True, timeout=300)

            if result.returncode != 0:
                self.log(f"reduce stderr: {result.stderr}")
                raise subprocess.CalledProcessError(result.returncode, cmd)

            self.log(f"Hydrogens added with reduce: {output_file}")
            return output_file

        except Exception as e:
            self.log(f"reduce failed: {e}")
            self.log("Falling back to PDBFixer...")
            return self._add_hydrogens_pdbfixer(input_file, output_file, 8.0)

    def _add_hydrogens_pdbfixer(self, input_file: str, output_file: str, ph: float) -> str:
        """Add hydrogens using PDBFixer"""
        if not HAS_PDBFIXER:
            self.error("PDBFixer is required. Install with: conda install -c conda-forge pdbfixer")

        try:
            fixer = PDBFixer(filename=input_file)
            fixer.findMissingResidues()
            fixer.findNonstandardResidues()
            fixer.replaceNonstandardResidues()
            fixer.findMissingAtoms()
            fixer.addMissingAtoms()
            fixer.addMissingHydrogens(ph)

            PDBFile.writeFile(fixer.topology, fixer.positions, open(output_file, 'w'))

            self.log(f"Hydrogens added with PDBFixer (pH {ph}): {output_file}")
            return output_file

        except Exception as e:
            self.error(f"PDBFixer failed: {e}")

    def generate_charges(self, pdb_file: str, forcefield: str = "AMBER",
                        ph: float = 8.0, pdb2pqr_cmd: str = "pdb2pqr30") -> str:
        """Generate partial charges using PDB2PQR"""
        self.log("Step 3: Generating partial charges with PDB2PQR...")

        input_path = Path(pdb_file)
        base_name = input_path.stem.replace('FH', '').replace('_clean', '')
        output_file = self.output_dir / f"{base_name}.pqr"

        try:
            cmd = [
                pdb2pqr_cmd,
                '--ff', forcefield.lower(),
                '--pH', str(ph),
                '--keep-chain',
                '--include-header',
                str(pdb_file),
                str(output_file)
            ]

            self.log(f"Running: {' '.join(cmd)}")

            result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)

            if result.returncode != 0:
                self.log(f"pdb2pqr stderr: {result.stderr}")
                # Try without some flags that might cause issues
                cmd_simple = [pdb2pqr_cmd, '--ff', forcefield.lower(),
                             str(pdb_file), str(output_file)]
                result = subprocess.run(cmd_simple, capture_output=True, text=True, timeout=600)

            if result.returncode != 0:
                raise subprocess.CalledProcessError(result.returncode, cmd)

            self.log(f"Charges generated: {output_file}")
            return str(output_file)

        except Exception as e:
            self.error(f"PDB2PQR failed: {e}")

    def validate_outputs(self, clean_pdb: str, h_pdb: str, pqr_file: str) -> bool:
        """Validate that outputs meet PRELYM requirements"""
        self.log("Validating outputs for PRELYM compatibility...")

        issues = []

        # Check clean PDB has no waters
        try:
            with open(clean_pdb, 'r') as f:
                if 'HOH' in f.read():
                    issues.append("Clean PDB still contains water molecules (HOH)")
        except Exception as e:
            issues.append(f"Cannot read clean PDB: {e}")

        # Check H-bond file has hydrogens
        try:
            with open(h_pdb, 'r') as f:
                content = f.read()
                if ' H' not in content and ' D' not in content:
                    issues.append("Hydrogen file appears to lack hydrogen atoms")
        except Exception as e:
            issues.append(f"Cannot read hydrogen file: {e}")

        # Check PQR has proper format
        try:
            with open(pqr_file, 'r') as f:
                lines = f.readlines()
                atom_lines = [l for l in lines if l.startswith('ATOM') or l.startswith('HETATM')]
                if not atom_lines:
                    issues.append("PQR file contains no ATOM/HETATM records")
                else:
                    # Check if charge column exists (should be floats)
                    try:
                        for line in atom_lines[:5]:  # Check first 5 ATOM lines
                            parts = line.split()
                            if len(parts) < 9:
                                raise ValueError("Not enough columns")
                            float(parts[8])  # Charge column
                    except (ValueError, IndexError):
                        issues.append("PQR file missing charge/radius columns")
        except Exception as e:
            issues.append(f"Cannot read PQR file: {e}")

        if issues:
            self.log("⚠️  Validation warnings:")
            for issue in issues:
                self.log(f"  - {issue}")
            return False
        else:
            self.log("✓ All files validated successfully!")
            return True

    def cleanup(self):
        """Clean up temporary files"""
        if self.temp_dir and os.path.exists(self.temp_dir):
            shutil.rmtree(self.temp_dir)

    def prepare_files(self, protein_name: str = None, pdb_id: str = None,
                     pdb_file: str = None, ph: float = 8.0,
                     forcefield: str = "AMBER") -> Tuple[str, str, str]:
        """Main workflow to prepare all three PRELYM input files"""

        # Check dependencies
        tools = self.check_dependencies()

        # Get PDB file
        if pdb_file:
            if not os.path.exists(pdb_file):
                self.error(f"PDB file not found: {pdb_file}")
            input_pdb = pdb_file
            self.log(f"Using provided PDB file: {pdb_file}")
        elif pdb_id:
            input_pdb = self.download_pdb(pdb_id)
        elif protein_name:
            found_pdb_id = self.search_protein(protein_name)
            if not found_pdb_id:
                self.error(f"No suitable PDB found for '{protein_name}'")
            input_pdb = self.download_pdb(found_pdb_id)
        else:
            self.error("Must provide --protein, --pdb-id, or --pdb-file")

        # Step 1: Clean PDB
        clean_pdb = self.clean_pdb(input_pdb)

        # Step 2: Add hydrogens
        h_pdb = self.add_hydrogens(clean_pdb, ph, tools['reduce'])

        # Step 3: Generate charges
        pqr_file = self.generate_charges(h_pdb, forcefield, ph,
                                       tools['pdb2pqr'] or 'pdb2pqr')

        # Validate outputs
        self.validate_outputs(clean_pdb, h_pdb, pqr_file)

        self.log("\n🎉 PRELYM preparation complete!")
        self.log("📁 Output files:")
        self.log(f"   Clean PDB:     {clean_pdb}")
        self.log(f"   H-bond file:   {h_pdb}")
        self.log(f"   Charge file:   {pqr_file}")
        self.log("\n📤 Upload these three files to PRELYM:")
        self.log("   https://prelym-webapp.onrender.com/")

        return clean_pdb, h_pdb, pqr_file


def main():
    parser = argparse.ArgumentParser(
        description="PRELYM File Preparation Agent",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s --protein "lysozyme"
  %(prog)s --pdb-id 4CHA
  %(prog)s --pdb-file myprotein.pdb
  %(prog)s --protein "chymotrypsin" --ph 7.0 --forcefield CHARMM
        """)

    # Input options (mutually exclusive)
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('--protein', help='Protein name to search for')
    input_group.add_argument('--pdb-id', help='PDB ID to download (e.g., 4CHA)')
    input_group.add_argument('--pdb-file', help='Local PDB file path')

    # Processing options
    parser.add_argument('--ph', type=float, default=8.0,
                       help='pH for hydrogen placement and charges (default: 8.0)')
    parser.add_argument('--forcefield', choices=['AMBER', 'CHARMM', 'PARSE'],
                       default='AMBER', help='Force field for charges (default: AMBER)')
    parser.add_argument('--output-dir', default='.',
                       help='Output directory for prepared files (default: current)')
    parser.add_argument('--quiet', action='store_true',
                       help='Suppress progress messages')

    args = parser.parse_args()

    # Create and run agent
    agent = PrelymPrepAgent(output_dir=args.output_dir, verbose=not args.quiet)

    try:
        agent.prepare_files(
            protein_name=args.protein,
            pdb_id=args.pdb_id,
            pdb_file=args.pdb_file,
            ph=args.ph,
            forcefield=args.forcefield
        )
    except KeyboardInterrupt:
        agent.log("\nOperation cancelled by user")
        sys.exit(1)
    except Exception as e:
        agent.error(f"Unexpected error: {e}")
    finally:
        agent.cleanup()


if __name__ == "__main__":
    main()