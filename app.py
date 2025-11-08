#!/usr/bin/env python3
"""
PRELYM Web Application
Flask-based web interface for protein lysine modification prediction
Preserves exact functionality of original atrp_predictor.py script
"""

import os
import uuid
import json
import shutil
from flask import Flask, render_template, request, jsonify, send_file, url_for
from werkzeug.utils import secure_filename
import pandas as pd
from threading import Thread
import time
from datetime import datetime

# Import ALL the original functions to maintain exact functionality
from atrp_predictor import (
    pdb2esa, pdb2pka, pdb2hbond, pdb2ss, pdb2res, pka_final,
    search, pdb2charge, final_data_table, decision_tree
)

# Import the PRELYM preparation agent (optional)
try:
    from prelym_prep_agent import PrelymPrepAgent
    PREP_AGENT_AVAILABLE = True
except ImportError:
    PREP_AGENT_AVAILABLE = False
    print("Warning: PRELYM preparation agent not available (requires pdbfixer/openmm)")

app = Flask(__name__)
app.config['SECRET_KEY'] = 'your-secret-key-here'
app.config['UPLOAD_FOLDER'] = 'uploads'
app.config['RESULTS_FOLDER'] = 'results'
app.config['MAX_CONTENT_LENGTH'] = 50 * 1024 * 1024  # 50MB max file size

# Create directories if they don't exist
os.makedirs(app.config['UPLOAD_FOLDER'], exist_ok=True)
os.makedirs(app.config['RESULTS_FOLDER'], exist_ok=True)

# Store job status with file persistence for production
import pickle
import os

# Store the original app directory for consistent file paths
APP_ROOT_DIR = os.path.abspath(os.path.dirname(__file__))

def load_job_status():
    """Load job status from file"""
    status_file = os.path.join(APP_ROOT_DIR, app.config['RESULTS_FOLDER'], 'job_status.pkl')
    try:
        if os.path.exists(status_file):
            with open(status_file, 'rb') as f:
                loaded_status = pickle.load(f)
                print(f"Loaded {len(loaded_status)} jobs from {status_file}")
                return loaded_status
        else:
            print(f"Job status file not found: {status_file}")
    except Exception as e:
        print(f"Error loading job status: {e}")
        import traceback
        traceback.print_exc()
    return {}

def save_job_status():
    """Save job status to file"""
    global job_status
    status_file = os.path.join(APP_ROOT_DIR, app.config['RESULTS_FOLDER'], 'job_status.pkl')
    try:
        os.makedirs(os.path.dirname(status_file), exist_ok=True)
        with open(status_file, 'wb') as f:
            pickle.dump(job_status, f)
        print(f"Job status saved to {status_file}")
    except Exception as e:
        print(f"Error saving job status: {e}")
        print(f"Current working directory: {os.getcwd()}")
        print(f"Target file: {status_file}")
        print(f"APP_ROOT_DIR: {APP_ROOT_DIR}")
        import traceback
        traceback.print_exc()

# Load existing job status or create new dict
job_status = load_job_status()

ALLOWED_EXTENSIONS = {'pdb', 'txt', 'pqr'}  # Added pqr for charge files

def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

def process_protein_analysis(job_id, pdb_file, hbond_file, charge_file, probe_radius):
    """Background task to process protein analysis - uses original decision_tree function exactly"""
    try:
        job_status[job_id]['status'] = 'processing'
        job_status[job_id]['progress'] = 10
        job_status[job_id]['message'] = 'Starting analysis...'
        save_job_status()

        # Change to job directory to match original script behavior
        original_cwd = os.getcwd()
        job_dir = os.path.dirname(pdb_file)
        os.chdir(job_dir)

        try:
            job_status[job_id]['progress'] = 20
            job_status[job_id]['message'] = 'Calculating Exposed Surface Area (ESA)...'

            job_status[job_id]['progress'] = 30
            job_status[job_id]['message'] = 'Calculating pKa values...'

            job_status[job_id]['progress'] = 40
            job_status[job_id]['message'] = 'Analyzing hydrogen bonds...'

            job_status[job_id]['progress'] = 50
            job_status[job_id]['message'] = 'Determining secondary structure...'

            job_status[job_id]['progress'] = 60
            job_status[job_id]['message'] = 'Calculating surface charges...'

            job_status[job_id]['progress'] = 70
            job_status[job_id]['message'] = 'Building final data table...'

            job_status[job_id]['progress'] = 80
            job_status[job_id]['message'] = 'Running decision tree analysis...'

            # Call the EXACT original function with the EXACT same parameters
            # This preserves all the original logic and output format
            decision_tree(
                os.path.basename(pdb_file),      # filename
                os.path.basename(hbond_file),    # filename1
                os.path.basename(charge_file),   # filename2
                probe_radius                     # probe_radius
            )

            job_status[job_id]['progress'] = 90
            job_status[job_id]['message'] = 'Finalizing results...'

            # The original function creates 'Results.xlsx' in current directory
            results_file = os.path.join(APP_ROOT_DIR, app.config['RESULTS_FOLDER'], f'{job_id}_results.xlsx')
            print(f"Moving Results.xlsx to: {results_file}")
            print(f"Current working directory: {os.getcwd()}")
            print(f"APP_ROOT_DIR: {APP_ROOT_DIR}")

            if os.path.exists('Results.xlsx'):
                shutil.move('Results.xlsx', results_file)
                job_status[job_id]['results_file'] = results_file
                print(f"Successfully moved results file. Stored path: {results_file}")
            else:
                print(f"Results.xlsx not found in current directory: {os.getcwd()}")
                print(f"Directory contents: {os.listdir('.')}")
                raise Exception("Results file was not generated")

            job_status[job_id]['status'] = 'completed'
            job_status[job_id]['progress'] = 100
            job_status[job_id]['message'] = 'Analysis completed successfully!'
            save_job_status()

        finally:
            # Always restore original working directory
            os.chdir(original_cwd)

    except Exception as e:
        os.chdir(original_cwd)  # Ensure we restore directory even on error
        job_status[job_id]['status'] = 'error'
        job_status[job_id]['message'] = f'Error: {str(e)}'
        save_job_status()
        print(f"Error in job {job_id}: {str(e)}")  # For debugging

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/prepare')
def prepare():
    return render_template('prepare.html', prep_agent_available=PREP_AGENT_AVAILABLE)

@app.route('/about')
def about():
    return render_template('about.html')

def validate_file_content(file_path, file_type):
    """Validate file content based on type"""
    try:
        with open(file_path, 'r') as f:
            content = f.read(1000)  # Read first 1000 chars

        if file_type == 'pdb':
            # Check for PDB format markers
            if not any(line.startswith(('ATOM', 'HETATM', 'HEADER')) for line in content.split('\n')):
                return False, "File doesn't appear to be a valid PDB file"

        elif file_type == 'charge':
            # Check for charge file format (should have numeric data)
            lines = content.strip().split('\n')
            if len(lines) < 2:
                return False, "Charge file appears to be empty or too short"

        return True, "Valid"

    except Exception as e:
        return False, f"Error reading file: {str(e)}"

@app.route('/upload', methods=['POST'])
def upload_files():
    try:
        # Check if all required files are present
        if 'pdb_file' not in request.files or 'hbond_file' not in request.files or 'charge_file' not in request.files:
            return jsonify({'error': 'Missing required files. Please upload all three files.'}), 400

        pdb_file = request.files['pdb_file']
        hbond_file = request.files['hbond_file']
        charge_file = request.files['charge_file']

        # Validate probe radius
        try:
            probe_radius = float(request.form.get('probe_radius', 1.4))
            if probe_radius <= 0 or probe_radius > 20:
                return jsonify({'error': 'Probe radius must be between 0.1 and 20.0 Å'}), 400
        except ValueError:
            return jsonify({'error': 'Invalid probe radius value'}), 400

        # Validate files exist
        if not all([pdb_file.filename, hbond_file.filename, charge_file.filename]):
            return jsonify({'error': 'All files must be selected'}), 400

        # Validate file extensions
        if not all([allowed_file(pdb_file.filename), allowed_file(hbond_file.filename), allowed_file(charge_file.filename)]):
            return jsonify({'error': 'Invalid file type. Only PDB, TXT, and PQR files are allowed'}), 400

        # Validate file sizes
        max_size = 50 * 1024 * 1024  # 50MB
        for file_obj in [pdb_file, hbond_file, charge_file]:
            if len(file_obj.read()) > max_size:
                return jsonify({'error': f'File {file_obj.filename} is too large (max 50MB)'}), 400
            file_obj.seek(0)  # Reset file pointer

        # Generate unique job ID
        job_id = str(uuid.uuid4())

        # Create job directory
        job_dir = os.path.join(app.config['UPLOAD_FOLDER'], job_id)
        os.makedirs(job_dir, exist_ok=True)

        # Save files with validation
        pdb_filename = secure_filename(pdb_file.filename)
        hbond_filename = secure_filename(hbond_file.filename)
        charge_filename = secure_filename(charge_file.filename)

        pdb_path = os.path.join(job_dir, pdb_filename)
        hbond_path = os.path.join(job_dir, hbond_filename)
        charge_path = os.path.join(job_dir, charge_filename)

        # Save files
        pdb_file.save(pdb_path)
        hbond_file.save(hbond_path)
        charge_file.save(charge_path)

        # Validate file contents
        validations = [
            (pdb_path, 'pdb', 'PDB structure file'),
            (hbond_path, 'pdb', 'H-bond analysis file'),
            (charge_path, 'charge', 'Charge file')
        ]

        for file_path, file_type, file_desc in validations:
            is_valid, message = validate_file_content(file_path, file_type)
            if not is_valid:
                # Clean up files
                shutil.rmtree(job_dir, ignore_errors=True)
                return jsonify({'error': f'{file_desc}: {message}'}), 400

        # Initialize job status
        job_status[job_id] = {
            'status': 'queued',
            'progress': 0,
            'message': 'Job queued for processing',
            'created_at': datetime.now().isoformat(),
            'files': {
                'pdb': pdb_filename,
                'hbond': hbond_filename,
                'charge': charge_filename
            },
            'probe_radius': probe_radius
        }
        save_job_status()  # Persist to file

        # Start background processing
        thread = Thread(target=process_protein_analysis,
                       args=(job_id, pdb_path, hbond_path, charge_path, probe_radius))
        thread.daemon = True  # Dies when main thread dies
        thread.start()

        return jsonify({'job_id': job_id, 'status': 'Job submitted successfully'})

    except Exception as e:
        return jsonify({'error': f'Unexpected error: {str(e)}'}), 500

@app.route('/status/<job_id>')
def get_job_status(job_id):
    global job_status
    if job_id not in job_status:
        print(f"Job {job_id} not found in job_status. Available jobs: {list(job_status.keys())}")
        # Try to reload job status from file
        job_status = load_job_status()
        if job_id not in job_status:
            print(f"Job {job_id} still not found after reload")
            return jsonify({'error': 'Job not found'}), 404

    return jsonify(job_status[job_id])

@app.route('/results/<job_id>')
def get_results(job_id):
    global job_status
    if job_id not in job_status:
        print(f"Job {job_id} not found in job_status. Available jobs: {list(job_status.keys())}")
        # Try to reload job status from file
        job_status = load_job_status()
        if job_id not in job_status:
            print(f"Job {job_id} still not found after reload")
            return jsonify({'error': 'Job not found'}), 404

    # If job exists but not completed, try reloading status in case it completed recently
    if job_status[job_id]['status'] != 'completed':
        print(f"Job {job_id} status is '{job_status[job_id]['status']}', reloading to check for updates")
        job_status = load_job_status()
        if job_id not in job_status or job_status[job_id]['status'] != 'completed':
            print(f"Job {job_id} still not completed after reload. Current status: {job_status.get(job_id, {}).get('status', 'not found')}")
            return jsonify({'error': 'Job not completed yet'}), 400

    results_file = job_status[job_id].get('results_file')
    if not results_file:
        print(f"No results_file in job_status for {job_id}")
        return jsonify({'error': 'Results file not found'}), 404

    if not os.path.exists(results_file):
        print(f"Results file does not exist: {results_file}")
        print(f"Current working directory: {os.getcwd()}")
        print(f"APP_ROOT_DIR: {APP_ROOT_DIR}")
        # Try to find the file with APP_ROOT_DIR
        alt_results_file = os.path.join(APP_ROOT_DIR, app.config['RESULTS_FOLDER'], f'{job_id}_results.xlsx')
        if os.path.exists(alt_results_file):
            print(f"Found results file at alternate location: {alt_results_file}")
            results_file = alt_results_file
        else:
            print(f"Alternate results file also not found: {alt_results_file}")
            return jsonify({'error': 'Results file not found'}), 404

    return send_file(results_file, as_attachment=True, download_name=f'prelym_results_{job_id}.xlsx')

@app.route('/preview/<job_id>')
def preview_results(job_id):
    global job_status
    if job_id not in job_status:
        print(f"Job {job_id} not found in job_status. Available jobs: {list(job_status.keys())}")
        # Try to reload job status from file
        job_status = load_job_status()
        if job_id not in job_status:
            print(f"Job {job_id} still not found after reload")
            return jsonify({'error': 'Job not found'}), 404

    # If job exists but not completed, try reloading status in case it completed recently
    if job_status[job_id]['status'] != 'completed':
        print(f"Job {job_id} status is '{job_status[job_id]['status']}', reloading to check for updates")
        job_status = load_job_status()
        if job_id not in job_status or job_status[job_id]['status'] != 'completed':
            print(f"Job {job_id} still not completed after reload. Current status: {job_status.get(job_id, {}).get('status', 'not found')}")
            return jsonify({'error': 'Job not completed yet'}), 400

    results_file = job_status[job_id].get('results_file')
    if not results_file:
        print(f"No results_file in job_status for {job_id}")
        return jsonify({'error': 'Results file not found'}), 404

    if not os.path.exists(results_file):
        print(f"Results file does not exist: {results_file}")
        # Try to find the file with APP_ROOT_DIR
        alt_results_file = os.path.join(APP_ROOT_DIR, app.config['RESULTS_FOLDER'], f'{job_id}_results.xlsx')
        if os.path.exists(alt_results_file):
            print(f"Found results file at alternate location: {alt_results_file}")
            results_file = alt_results_file
        else:
            print(f"Alternate results file also not found: {alt_results_file}")
            return jsonify({'error': 'Results file not found'}), 404

    try:
        # Read the Excel file and convert to JSON for preview
        df = pd.read_excel(results_file)
        # Limit to first 100 rows for preview
        preview_df = df.head(100)

        return jsonify({
            'columns': list(preview_df.columns),
            'data': preview_df.to_dict('records'),
            'total_rows': len(df),
            'preview_rows': len(preview_df)
        })
    except Exception as e:
        return jsonify({'error': f'Error reading results: {str(e)}'}), 500

@app.route('/logo.png')
def serve_logo():
    return send_file('LOGO.PNG')

@app.route('/chmyotrypsin/<filename>')
def serve_example_file(filename):
    """Serve example chymotrypsin files for download"""
    try:
        file_path = os.path.join('chmyotrypsin', filename)
        if os.path.exists(file_path):
            return send_file(file_path, as_attachment=True, download_name=filename)
        else:
            return jsonify({'error': 'File not found'}), 404
    except Exception as e:
        return jsonify({'error': str(e)}), 500

# Store preparation job status
prep_job_status = {}

def process_preparation(job_id, input_type, input_value, ph, forcefield, output_dir):
    """Background task to process file preparation using PRELYM preparation agent"""
    try:
        prep_job_status[job_id]['status'] = 'processing'
        prep_job_status[job_id]['progress'] = 10
        prep_job_status[job_id]['message'] = 'Initializing preparation agent...'

        # Initialize the preparation agent
        agent = PrelymPrepAgent()

        prep_job_status[job_id]['progress'] = 20
        prep_job_status[job_id]['message'] = 'Checking system dependencies...'

        # Check dependencies availability
        tools = agent.check_dependencies()
        if not tools['pdb2pqr']:
            raise Exception("pdb2pqr is required but not found. Please install it via conda or from pdb2pqr.readthedocs.io")

        prep_job_status[job_id]['progress'] = 30

        # Get or download PDB file
        pdb_file = None
        if input_type == 'protein':
            prep_job_status[job_id]['message'] = f'Searching for protein: {input_value}'
            pdb_id = agent.search_protein(input_value)
            if not pdb_id:
                raise Exception(f"No PDB structures found for protein: {input_value}")
            prep_job_status[job_id]['message'] = f'Found PDB ID: {pdb_id}. Downloading...'
            pdb_file = agent.download_pdb(pdb_id)
        elif input_type == 'pdb_id':
            prep_job_status[job_id]['message'] = f'Downloading PDB ID: {input_value}'
            pdb_file = agent.download_pdb(input_value)
        elif input_type == 'pdb_file':
            pdb_file = input_value
            prep_job_status[job_id]['message'] = f'Using uploaded file: {os.path.basename(pdb_file)}'

        prep_job_status[job_id]['progress'] = 50
        prep_job_status[job_id]['message'] = 'Cleaning PDB structure...'

        # Step 1: Clean PDB
        clean_pdb = agent.clean_pdb(pdb_file)

        prep_job_status[job_id]['progress'] = 70
        prep_job_status[job_id]['message'] = 'Adding hydrogens...'

        # Step 2: Add hydrogens
        h_pdb = agent.add_hydrogens(clean_pdb, ph=ph)

        prep_job_status[job_id]['progress'] = 85
        prep_job_status[job_id]['message'] = 'Generating charges...'

        # Step 3: Generate charges
        pqr_file = agent.generate_charges(clean_pdb, forcefield=forcefield, ph=ph)

        prep_job_status[job_id]['progress'] = 95
        prep_job_status[job_id]['message'] = 'Validating outputs...'

        # Validate outputs
        if not agent.validate_outputs(clean_pdb, h_pdb, pqr_file):
            raise Exception("Output validation failed. Files may not meet PRELYM requirements.")

        prep_job_status[job_id]['status'] = 'completed'
        prep_job_status[job_id]['progress'] = 100
        prep_job_status[job_id]['message'] = 'File preparation completed successfully!'
        prep_job_status[job_id]['files'] = {
            'clean_pdb': clean_pdb,
            'h_pdb': h_pdb,
            'pqr_file': pqr_file
        }

        # Cleanup temporary files
        agent.cleanup()

    except Exception as e:
        prep_job_status[job_id]['status'] = 'failed'
        prep_job_status[job_id]['message'] = f'Error: {str(e)}'
        print(f"Preparation failed for job {job_id}: {str(e)}")
        import traceback
        traceback.print_exc()

@app.route('/prepare-files', methods=['POST'])
def prepare_files():
    """Start automated file preparation using PRELYM preparation agent"""
    if not PREP_AGENT_AVAILABLE:
        return jsonify({'error': 'PRELYM preparation agent not available (requires pdbfixer/openmm)'}), 503

    try:
        data = request.get_json()

        input_type = data.get('input_type')  # 'protein', 'pdb_id', or 'pdb_file'
        input_value = data.get('input_value')
        ph = float(data.get('ph', 8.0))
        forcefield = data.get('forcefield', 'AMBER')

        if not input_type or not input_value:
            return jsonify({'error': 'input_type and input_value are required'}), 400

        # Generate unique job ID
        job_id = str(uuid.uuid4())

        # Initialize job status
        prep_job_status[job_id] = {
            'status': 'queued',
            'progress': 0,
            'message': 'Job queued for processing...',
            'created_at': datetime.now().isoformat(),
            'input_type': input_type,
            'input_value': input_value,
            'ph': ph,
            'forcefield': forcefield
        }

        # Create output directory
        output_dir = os.path.join(app.config['RESULTS_FOLDER'], job_id)
        os.makedirs(output_dir, exist_ok=True)

        # Start background preparation task
        thread = Thread(target=process_preparation,
                       args=(job_id, input_type, input_value, ph, forcefield, output_dir))
        thread.start()

        return jsonify({
            'job_id': job_id,
            'status': 'queued',
            'message': 'File preparation started'
        })

    except Exception as e:
        return jsonify({'error': str(e)}), 500

@app.route('/prepare-status/<job_id>')
def get_preparation_status(job_id):
    """Get status of file preparation job"""
    if job_id not in prep_job_status:
        return jsonify({'error': 'Job not found'}), 404

    return jsonify(prep_job_status[job_id])

@app.route('/prepare-download/<job_id>/<file_type>')
def download_prepared_file(job_id, file_type):
    """Download prepared files"""
    if job_id not in prep_job_status:
        return jsonify({'error': 'Job not found'}), 404

    job = prep_job_status[job_id]
    if job['status'] != 'completed':
        return jsonify({'error': 'Job not completed'}), 400

    if 'files' not in job:
        return jsonify({'error': 'No files available'}), 404

    file_map = {
        'clean': job['files']['clean_pdb'],
        'hydrogen': job['files']['h_pdb'],
        'charge': job['files']['pqr_file']
    }

    if file_type not in file_map:
        return jsonify({'error': 'Invalid file type'}), 400

    file_path = file_map[file_type]
    if not os.path.exists(file_path):
        return jsonify({'error': 'File not found'}), 404

    return send_file(file_path, as_attachment=True, download_name=os.path.basename(file_path))

if __name__ == '__main__':
    import os
    port = int(os.environ.get('PORT', 5002))
    debug_mode = os.environ.get('FLASK_DEBUG', '0').lower() in ['1', 'true', 'yes']
    app.run(debug=debug_mode, host='0.0.0.0', port=port)