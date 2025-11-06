#!/usr/bin/env python3
"""
Setup script for PRELYM web application
This script helps users set up the web application with proper dependencies
"""

import os
import sys
import subprocess
import platform

def check_conda():
    """Check if conda is available"""
    try:
        subprocess.run(['conda', '--version'], capture_output=True, check=True)
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        return False

def setup_conda_environment():
    """Set up conda environment from environment.yml"""
    print("Setting up conda environment...")
    try:
        # Create environment from yml file
        subprocess.run(['conda', 'env', 'create', '-f', 'environment.yml'], check=True)
        print("✓ Conda environment 'ATRP' created successfully")

        # Install Flask in the environment
        subprocess.run(['conda', 'run', '-n', 'ATRP', 'pip', 'install', 'Flask==2.3.3', 'Werkzeug==2.3.7'], check=True)
        print("✓ Flask installed in conda environment")

        return True
    except subprocess.CalledProcessError as e:
        print(f"✗ Error setting up conda environment: {e}")
        return False

def setup_pip_environment():
    """Set up using pip and requirements.txt"""
    print("Setting up using pip...")
    try:
        subprocess.run([sys.executable, '-m', 'pip', 'install', '-r', 'requirements.txt'], check=True)
        print("✓ Dependencies installed with pip")
        return True
    except subprocess.CalledProcessError as e:
        print(f"✗ Error installing with pip: {e}")
        return False

def check_system_dependencies():
    """Check for system-level dependencies"""
    print("Checking system dependencies...")

    # Check for DSSP
    try:
        subprocess.run(['dssp', '--version'], capture_output=True, check=True)
        print("✓ DSSP found")
    except (subprocess.CalledProcessError, FileNotFoundError):
        print("⚠ DSSP not found. Install via: conda install -c speleo3 dssp")

    return True

def create_directories():
    """Create necessary directories"""
    directories = ['uploads', 'results', 'static', 'templates']
    for directory in directories:
        os.makedirs(directory, exist_ok=True)
        print(f"✓ Created directory: {directory}")

def main():
    print("PRELYM Web Application Setup")
    print("=" * 40)

    # Create directories
    create_directories()

    # Check system dependencies
    check_system_dependencies()

    # Setup Python environment
    if check_conda() and os.path.exists('environment.yml'):
        print("\nUsing conda for setup (recommended)...")
        if setup_conda_environment():
            print("\n✓ Setup complete!")
            print("\nTo run the application:")
            print("1. conda activate ATRP")
            print("2. python app.py")
        else:
            print("\nTrying pip fallback...")
            setup_pip_environment()
    else:
        print("\nConda not available, using pip...")
        setup_pip_environment()

    print("\n" + "=" * 40)
    print("Setup complete! The web application is ready to run.")
    print("\nTo start the server:")
    print("python app.py")
    print("\nThen open: http://localhost:5000")

if __name__ == '__main__':
    main()