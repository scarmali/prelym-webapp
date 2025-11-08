# PRELYM Preparation Agent

Automated preparation of the three input files needed for PRELYM protein lysine modification prediction analysis.

## Quick Start

```bash
# Install dependencies
pip install -r prelym_prep_requirements.txt
conda install -c conda-forge reduce pdb2pqr

# Use the agent
python prelym_prep_agent.py --protein "lysozyme"
python prelym_prep_agent.py --pdb-id 4CHA
python prelym_prep_agent.py --pdb-file myprotein.pdb
```

## What It Does

The agent automates the [PRELYM File Preparation Guide](https://prelym-webapp.onrender.com/prepare) workflow:

1. **🔍 Protein Search** (if `--protein` given)
   - Searches RCSB PDB using full-text search
   - Ranks candidates by resolution and selects best X-ray structure
   - Downloads the PDB file automatically

2. **🧹 Step 1: Clean PDB** → `*_clean.pdb`
   - Removes water molecules (HOH)
   - Keeps non-water heterogens
   - Preserves protein structure integrity

3. **⚛️ Step 2: Add Hydrogens** → `*FH.pdb`
   - **Preferred**: Uses `reduce` (MolProbity suite) for optimal hydrogen placement
   - **Fallback**: Uses PDBFixer if reduce not available
   - Optimizes hydrogen bond networks

4. **⚡ Step 3: Generate Charges** → `*.pqr`
   - Uses PDB2PQR with specified force field (default: AMBER)
   - Assigns partial charges and atomic radii
   - Sets protonation states at specified pH (default: 8.0)

5. **✅ Validation**
   - Verifies no waters in clean PDB
   - Confirms hydrogens present in H-bond file
   - Validates PQR format with charge/radius columns

## Installation

### Quick Setup (Conda Recommended)
```bash
# Create environment with all tools
conda create -n prelym-prep python=3.9
conda activate prelym-prep
conda install -c conda-forge biopython pdbfixer openmm reduce pdb2pqr

# Or install with pip + conda tools
pip install -r prelym_prep_requirements.txt
conda install -c conda-forge reduce pdb2pqr
```

### Manual Installation
1. **Python dependencies**: `pip install biopython pdbfixer openmm`
2. **reduce**: Download from [MolProbity](http://molprobity.biochem.duke.edu/)
3. **pdb2pqr**: Download from [pdb2pqr.readthedocs.io](https://pdb2pqr.readthedocs.io/)

## Usage Examples

### Search and Prepare by Protein Name
```bash
python prelym_prep_agent.py --protein "chymotrypsin"
# Searches RCSB, finds best structure, prepares all files
```

### Use Specific PDB ID
```bash
python prelym_prep_agent.py --pdb-id 4CHA
# Downloads 4CHA.pdb and prepares files
```

### Use Local PDB File
```bash
python prelym_prep_agent.py --pdb-file myprotein.pdb
# Uses your existing PDB file
```

### Advanced Options
```bash
python prelym_prep_agent.py --protein "lysozyme" \\
    --ph 7.0 \\
    --forcefield CHARMM \\
    --output-dir ./prepared_files
```

## Output Files

For input `4CHA`, the agent generates:
- `4CHA_clean.pdb` - Clean structure (no waters)
- `4CHAFH.pdb` - Hydrogen-optimized structure
- `4CHA.pqr` - Partial charges and radii

**Ready to upload to PRELYM**: https://prelym-webapp.onrender.com/

## Command Line Options

```
Required (one of):
  --protein NAME      Search RCSB for protein by name
  --pdb-id ID         Download specific PDB ID (e.g., 4CHA)
  --pdb-file PATH     Use local PDB file

Optional:
  --ph FLOAT          pH for protonation (default: 8.0)
  --forcefield FF     AMBER, CHARMM, or PARSE (default: AMBER)
  --output-dir DIR    Output directory (default: current)
  --quiet             Suppress progress messages
```

## Dependencies

### Required Python Packages
- **biopython** - PDB parsing and manipulation
- **pdbfixer** - Missing atom/hydrogen addition (fallback)
- **openmm** - Molecular simulation toolkit

### Required External Tools
- **pdb2pqr** - Partial charge assignment (REQUIRED)
- **reduce** - Optimal hydrogen placement (recommended)

### Installation Check
The agent automatically checks for tools and reports what's available:
```
[PRELYM Prep] ✓ reduce found
[PRELYM Prep] ✓ pdb2pqr30 found
```

## Troubleshooting

### "pdb2pqr not found"
```bash
conda install -c conda-forge pdb2pqr
# or download from pdb2pqr.readthedocs.io
```

### "reduce not found"
Agent will use PDBFixer fallback, but for best results:
```bash
conda install -c conda-forge reduce
```

### "No results found for protein"
Try different search terms:
```bash
python prelym_prep_agent.py --protein "hen egg white lysozyme"
python prelym_prep_agent.py --protein "HEWL"
```

### Large proteins timeout
Increase timeout or use local file:
```bash
python prelym_prep_agent.py --pdb-file large_protein.pdb
```

## Integration with PRELYM

1. **Run the agent** to prepare files
2. **Upload to PRELYM** at https://prelym-webapp.onrender.com/
3. **Set probe radius** and run analysis
4. **Download results** when complete

The agent ensures files meet PRELYM's exact requirements for successful analysis.

## Advanced Usage

### Batch Processing
```bash
for protein in "lysozyme" "chymotrypsin" "trypsin"; do
    python prelym_prep_agent.py --protein "$protein" --output-dir "$protein"_prep
done
```

### Custom pH/Force Field
```bash
# Physiological pH with CHARMM
python prelym_prep_agent.py --pdb-id 1LYZ --ph 7.4 --forcefield CHARMM
```

## Notes

- ⚠️ **PRELYM analysis** still runs on the website - this agent only prepares input files
- 🔬 **Force field choice** affects charge assignment - AMBER works well for most proteins
- 🧪 **pH matters** for protonation states - use physiological pH (7.4) or experimental conditions
- 📊 **Resolution ranking** prioritizes high-quality X-ray structures when searching

---

**Need help?** Check the [PRELYM website](https://prelym-webapp.onrender.com/) for more information about lysine modification prediction.