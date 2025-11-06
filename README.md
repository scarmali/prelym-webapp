# PRELYM Web Application
PRELYM: Automated Prediction of Lysine Modifications on Proteins

A user-friendly web interface for predicting lysine reactivity patterns in protein structures using computational analysis.

## Features

- **Web-based Interface**: Easy-to-use drag-and-drop file upload
- **Real-time Progress Tracking**: Monitor analysis progress with detailed status updates
- **Comprehensive Analysis**:
  - Exposed Surface Area (ESA) calculation using FreeSASA
  - pKa prediction using PROPKA
  - Hydrogen bond analysis using MDTraj
  - Secondary structure determination using DSSP
  - Surface charge calculations
  - Decision tree-based reactivity prediction
- **Results Visualization**: Interactive table preview and Excel download
- **Robust Error Handling**: File validation and detailed error messages

## Installation

### Option 1: Using Conda (Recommended)

1. Clone or download this repository
2. Navigate to the directory:
   ```bash
   cd prelym-webapp
   ```
3. Run the setup script:
   ```bash
   python setup.py
   ```
4. Activate the environment:
   ```bash
   conda activate ATRP
   ```

### Option 2: Using pip

1. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```
2. Install DSSP separately (required for secondary structure analysis):
   ```bash
   # On Ubuntu/Debian:
   sudo apt-get install dssp

   # On macOS with Homebrew:
   brew install dssp

   # Using conda:
   conda install -c speleo3 dssp
   ```

## Usage

1. Start the web application:
   ```bash
   python app.py
   ```

2. Open your web browser and navigate to:
   ```
   http://localhost:5000
   ```

3. Upload your files:
   - **PDB Structure File**: Main protein structure (.pdb)
   - **H-Bond Analysis File**: PDB file for hydrogen bond analysis (.pdb)
   - **Charge File**: File containing charge information (.pqr, .txt)

4. Set the probe radius (default: 1.4 Å)

5. Click "Start Analysis" and monitor progress

6. Download results as Excel file when complete

## File Requirements

### PDB Structure File
- Standard PDB format
- Contains ATOM or HETATM records
- Should include lysine residues for analysis

### H-Bond Analysis File
- PDB format file used for hydrogen bond calculations
- Can be the same as the main PDB file or a processed version

### Charge File
- Contains atomic charge information
- Formats supported: PQR, TXT
- Should include coordinates and partial charges

## Output

The analysis generates an Excel file (`Results.xlsx`) containing:

- **Chain**: Protein chain identifier
- **Residue**: Residue number
- **Amino Group**: Residue type (LYS for lysine, N for N-terminus)
- **ESA**: Exposed Surface Area (Ų)
- **pKa**: Predicted pKa value
- **Secondary Structure**: Helix, Strand, or Coil
- **H-Donor**: Whether residue acts as hydrogen donor (Yes/No)
- **Area Of Lower Charge**: Whether in low charge density area (Yes/No)
- **Interaction**: Final prediction (fast-reacting, slow-reacting, non-reacting)

## Technical Details

The web application preserves the exact functionality of the original PRELYM script:

1. **ESA Calculation**: Uses FreeSASA with configurable probe radius
2. **pKa Prediction**: PROPKA-based calculations via propkatraj
3. **Hydrogen Bonding**: Baker-Hubbard algorithm implementation
4. **Secondary Structure**: DSSP analysis for structure classification
5. **Surface Charges**: Electrostatic calculations
6. **Decision Tree**: Rule-based classification of lysine reactivity

## Dependencies

Core scientific libraries:
- BioPython (1.78)
- MDAnalysis (1.0.0)
- FreeSASA (2.1.0)
- PROPKA/propkatraj (1.1.0)
- pandas, numpy, scipy

Web framework:
- Flask (2.3.3)
- Werkzeug (2.3.7)

External tools:
- DSSP (for secondary structure)

## Troubleshooting

### Common Issues

1. **DSSP not found**: Install DSSP using conda or system package manager
2. **File upload errors**: Check file formats and sizes (max 50MB per file)
3. **Analysis fails**: Ensure PDB files contain valid ATOM records
4. **Dependencies missing**: Run setup.py or install from requirements.txt

### Error Messages

- **"File doesn't appear to be a valid PDB file"**: Check PDB format
- **"Charge file appears to be empty"**: Verify charge file contains data
- **"Probe radius must be between 0.1 and 10.0 Å"**: Adjust probe radius value

## Development

To modify or extend the application:

1. Original analysis functions are in `atrp_predictor.py`
2. Web interface logic is in `app.py`
3. Frontend templates are in `templates/`
4. Static files go in `static/` (if needed)

## License

See LICENSE.txt for details.

## Citation

If you use PRELYM in your research, please cite the original work.
