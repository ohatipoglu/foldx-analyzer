# macOS Setup Guide for Bioinformatics Analyzer Suite

This guide provides detailed instructions for setting up and running the Bioinformatics Analyzer Suite (FoldX & PLIP) on macOS.

---

## Table of Contents

1. [Prerequisites](#prerequisites)
2. [Step 1: Install Python](#step-1-install-python)
3. [Step 2: Clone the Repository](#step-2-clone-the-repository)
4. [Step 3: Create Virtual Environment](#step-3-create-virtual-environment)
5. [Step 4: Install Dependencies](#step-4-install-dependencies)
6. [Step 5: Install PLIP and OpenBabel](#step-5-install-plip-and-openbabel)
7. [Step 6: Install PyMOL](#step-6-install-pymol)
8. [Step 7: Configure the Application](#step-7-configure-the-application)
9. [Step 8: Run the Application](#step-8-run-the-application)
10. [Step 9: Run Tests](#step-9-run-tests)
11. [Troubleshooting](#troubleshooting)
12. [Quick Reference Commands](#quick-reference-commands)
13. [Anaconda Quick Start Guide](#anaconda-quick-start-guide)
14. [Version Compatibility Matrix](#version-compatibility-matrix)
15. [Uninstallation](#uninstallation)

---

## Prerequisites

Before you begin, ensure you have the following:

- **macOS** version 10.15 (Catalina) or later
- **Administrator access** to install software
- **Internet connection** for downloading packages
- **Terminal access** (built-in Terminal app or iTerm2)
- At least **2GB of free disk space**

---

## Step 1: Install Python

### Option A: Using Homebrew (Recommended)

1. **Install Homebrew** (if not already installed):

   Open Terminal and run:

   ```bash
   /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
   ```

   Follow the on-screen instructions. After installation, add Homebrew to your PATH:

   ```bash
   echo 'eval "$(/opt/homebrew/bin/brew shellenv)"' >> ~/.zprofile
   eval "$(/opt/homebrew/bin/brew shellenv)"
   ```

   > **Note:** On Intel Macs, Homebrew is installed to `/usr/local/bin` instead of `/opt/homebrew/bin`.

2. **Install Python 3.12 or later**:

   ```bash
   brew install python@3.12
   ```

3. **Verify installation**:

   ```bash
   python3 --version
   pip3 --version
   ```

### Option B: Using Official Python Installer

1. Download Python from [python.org](https://www.python.org/downloads/macos/)
2. Run the `.pkg` installer
3. Follow the installation wizard
4. Verify installation:

   ```bash
   python3 --version
   pip3 --version
   ```

---

## Step 2: Clone the Repository

1. **Navigate to your projects directory**:

   ```bash
   cd ~/Projects
   ```

   If the directory doesn't exist, create it:

   ```bash
   mkdir -p ~/Projects
   cd ~/Projects
   ```

2. **Clone the repository**:

   ```bash
   git clone https://github.com/ohatipoglu/foldx-analyzer.git
   cd foldx-analyzer
   ```

3. **Verify the repository structure**:

   ```bash
   ls -la
   ```

   You should see files like `combined_gui.py`, `foldx_core.py`, `pyproject.toml`, etc.

---

## Step 3: Create Virtual Environment

Creating a virtual environment isolates project dependencies from your system Python.

### Option A: Using Python venv

1. **Navigate to the project directory**:

   ```bash
   cd ~/Projects/foldx-analyzer
   ```

2. **Create a virtual environment**:

   ```bash
   python3 -m venv venv
   ```

3. **Activate the virtual environment**:

   For **zsh** (default on macOS Catalina and later):

   ```bash
   source venv/bin/activate
   ```

   For **bash**:

   ```bash
   source venv/bin/activate
   ```

   You should see `(venv)` prefix in your terminal prompt, indicating the virtual environment is active.

4. **Verify virtual environment**:

   ```bash
   which python
   which pip
   ```

   Both should point to paths inside the `venv` directory.

### Option B: Using Conda (Anaconda or Miniconda)

> **Note:** This section applies to both **Anaconda** (full distribution) and **Miniconda** (minimal installation).

1. **Verify conda installation**:

   ```bash
   conda --version
   ```

   If not found, initialize conda:

   ```bash
   # For Anaconda (default install location)
   ~/anaconda3/bin/conda init zsh
   
   # For Miniconda (default install location)
   ~/miniconda3/bin/conda init zsh
   
   # Restart terminal or run:
   source ~/.zshrc
   ```

2. **Navigate to the project directory**:

   ```bash
   cd ~/Projects/foldx-analyzer
   ```

3. **Create a conda environment**:

   ```bash
   conda create -n foldx-env python=3.12
   ```

4. **Activate the conda environment**:

   ```bash
   conda activate foldx-env
   ```

   You should see `(foldx-env)` prefix in your terminal prompt.

5. **Verify environment**:

   ```bash
   which python
   which conda
   ```

> **Tip:** If you already have Anaconda installed, you can skip to Step 4 and use your existing conda installation. Only create a new environment for this project.

---

## Step 4: Install Dependencies

### Option A: Using pip (venv or conda environment)

1. **Upgrade pip**:

   ```bash
   pip install --upgrade pip
   ```

2. **Install project dependencies**:

   ```bash
   pip install -e .
   ```

   This installs all dependencies listed in `pyproject.toml`:
   - pandas
   - numpy
   - matplotlib
   - customtkinter
   - biopython
   - python-docx

3. **Install test dependencies** (optional):

   ```bash
   pip install pytest
   ```

### Option B: Using conda (Anaconda/Miniconda only)

> **Note:** Anaconda comes with many scientific packages pre-installed (numpy, pandas, matplotlib). You may only need to install the missing ones.

1. **Install core dependencies from conda-forge**:

   ```bash
   conda install -c conda-forge pandas numpy matplotlib biopython
   ```

2. **Install remaining packages via pip**:

   ```bash
   pip install customtkinter python-docx
   ```

3. **Install the project in editable mode** (optional):

   ```bash
   pip install -e .
   ```

4. **Install test dependencies** (optional):

   ```bash
   conda install -c conda-forge pytest
   ```

> **Tip:** Using `conda install` for scientific packages is often faster and more stable on macOS, especially for Apple Silicon (M1/M2/M3).

---

## Step 5: Install PLIP and OpenBabel

PLIP (Protein-Ligand Interaction Profiler) requires OpenBabel and has specific installation requirements on macOS.

### Option A: Using Conda (Recommended for Anaconda/Miniconda Users)

> **Note for Anaconda Users:** You already have conda installed. Skip to step 2.

1. **Install Miniconda** (if not already installed and not using Anaconda):

   ```bash
   curl -fsSL https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh -o miniconda.sh
   bash miniconda.sh
   ```

   > **Note:** For Intel Macs, use `Miniconda3-latest-MacOSX-x86_64.sh`

2. **Initialize conda** (if not already initialized):

   ```bash
   # For Anaconda
   ~/anaconda3/bin/conda init zsh
   
   # For Miniconda
   ~/miniconda3/bin/conda init zsh
   
   # Restart terminal or run:
   source ~/.zshrc
   ```

3. **Activate your conda environment**:

   ```bash
   conda activate foldx-env
   ```

4. **Install OpenBabel and PLIP**:

   ```bash
   conda install -c conda-forge openbabel plip
   ```

### Option B: Using Homebrew (Alternative)

1. **Install OpenBabel**:

   ```bash
   brew install open-babel
   ```

2. **Install PLIP via pip**:

   ```bash
   pip install plip
   ```

   > **Warning:** This method may fail due to compilation issues. Conda is preferred.

3. **Verify installation**:

   ```bash
   plip --version
   ```

---

## Step 6: Install PyMOL

PyMOL is required for 3D visualization of protein-ligand interactions.

### Option A: Using Conda (Recommended for Anaconda/Miniconda Users)

> **Note for Anaconda Users:** This is the recommended method as it integrates seamlessly with your existing conda environment.

1. **Ensure your conda environment is active**:

   ```bash
   conda activate foldx-env
   ```

2. **Install PyMOL via conda**:

   ```bash
   conda install -c conda-forge pymol
   ```

### Option B: Using Homebrew

```bash
brew install pymol
```

### Option C: Official PyMOL Installer

1. Download from [PyMOL website](https://pymol.org/2/)
2. Run the installer
3. Follow the installation wizard

### Verify PyMOL Installation

```bash
pymol --version
```

Or launch PyMOL:

```bash
pymol
```

---

## Step 7: Configure the Application

The application creates a `config.ini` file on first run. However, you may want to configure it manually.

1. **Check if config.ini exists**:

   ```bash
   ls -la config.ini
   ```

2. **If it doesn't exist, run the application once** to generate it, or create manually:

   ```bash
   cat > config.ini << EOF
   [Settings]
   pymol_executable_path = pymol

   [Graph]
   figure_dpi = 100
   save_dpi = 300
   stabilizing_threshold = -0.5
   destabilizing_threshold = 0.5

   [PDB]
   neighbor_threshold = 6.0
   disulfide_min_dist = 2.0
   disulfide_max_dist = 2.2
   EOF
   ```

3. **Configure PyMOL path** (if not in system PATH):

   Find PyMOL executable location:

   ```bash
   which pymol
   ```

   Update `config.ini` with the full path:

   ```ini
   [Settings]
   pymol_executable_path = /path/to/pymol
   ```

   Example for conda installation:

   ```ini
   [Settings]
   pymol_executable_path = /Users/yourusername/miniconda3/envs/foldx-env/bin/pymol
   ```

---

## Step 8: Run the Application

1. **Ensure virtual environment is activated**:

   ```bash
   source venv/bin/activate
   ```

   Or if using conda:

   ```bash
   conda activate foldx-env
   ```

2. **Navigate to project directory**:

   ```bash
   cd ~/Projects/foldx-analyzer
   ```

3. **Run the application**:

   ```bash
   python combined_gui.py
   ```

4. **Expected behavior**:

   - A GUI window should open with title "Bioinformatics Analyzer Suite"
   - The window should have a sidebar with navigation options
   - You should see the FoldX Analyzer panel by default

### Troubleshooting Launch Issues

**Error: `ImportError: No module named 'customtkinter'`**

```bash
pip install customtkinter
```

**Error: `ImportError: No module named 'plip'`**

```bash
conda install -c conda-forge plip
```

**Error: `_tkinter.TclError: no display name`**

This occurs when running over SSH without X11 forwarding. Run locally on the Mac.

**Error: `PyMOL executable not found`**

Ensure PyMOL is installed and in PATH:

```bash
which pymol
```

If not found, update `config.ini` with the full path to PyMOL executable.

---

## Step 9: Run Tests

Verify your installation by running the test suite:

1. **Navigate to project directory**:

   ```bash
   cd ~/Projects/foldx-analyzer
   ```

2. **Run pytest**:

   ```bash
   python -m pytest tests/ -v
   ```

3. **Expected output**:

   ```
   ============================= test session starts =============================
   platform darwin -- Python 3.12.x, pytest-8.x.x
   collected 20 items

   tests/test_foldx_core.py::TestDetectCommand::test_determine_command_analysecomplex PASSED
   tests/test_foldx_core.py::TestDetectCommand::test_determine_command_buildmodel PASSED
   ...
   ========================= 20 passed in X.XXs ==================================
   ```

---

## Troubleshooting

### Common Issues and Solutions

#### 1. **Permission Denied Errors**

**Symptom:** `PermissionError: [Errno 13] Permission denied`

**Solution:**

```bash
sudo chown -R $(whoami) ~/Projects/foldx-analyzer
```

Or run with appropriate permissions:

```bash
chmod +x combined_gui.py
```

#### 2. **Certificate Verification Errors**

**Symptom:** `ssl.SSLCertVerificationError` when installing packages

**Solution:**

```bash
/Applications/Python\ 3.12/Install\ Certificates.command
```

Or use:

```bash
pip install --trusted-host pypi.org --trusted-host files.pythonhosted.org <package-name>
```

#### 3. **OpenBabel Not Found**

**Symptom:** `ImportError: No module named 'openbabel'`

**Solution:**

```bash
conda install -c conda-forge openbabel
```

Verify:

```bash
which obabel
```

#### 4. **Tkinter GUI Issues**

**Symptom:** Application crashes or GUI doesn't render properly

**Solution:**

macOS comes with Tcl/Tk, but it may be outdated. Install updated version:

```bash
brew install tcl-tk
```

Then reinstall Python with Homebrew to use the updated Tcl/Tk.

#### 5. **Memory Issues with Large Files**

**Symptom:** Application becomes unresponsive with large `.fxout` files

**Solution:**

- Close other applications to free memory
- Process files in smaller batches
- Increase swap space (advanced users)

#### 6. **PyMOL Visualization Not Working**

**Symptom:** PLIP analysis runs but no visualization appears

**Solution:**

1. Check PyMOL is in PATH:
   ```bash
   which pymol
   ```

2. Test PyMOL independently:
   ```bash
   pymol --version
   ```

3. Update `config.ini` with full PyMOL path

4. Check console logs for error messages

#### 7. **Apple Silicon (M1/M2/M3) Specific Issues**

Some Python packages may have compatibility issues with Apple Silicon.

**Solution:**

Use Rosetta 2 for x86 packages:

```bash
softwareupdate --install-rosetta
```

Or use conda with arm64 packages:

```bash
conda install -c conda-forge --force-reinstall <package-name>
```

#### 8. **matplotlib Backend Issues**

**Symptom:** `UserWarning: matplotlib is currently using agg, which is a non-GUI backend`

**Solution:**

Install PyQt5 or TkAgg backend:

```bash
pip install PyQt5
```

Or set backend in matplotlib config:

```bash
mkdir -p ~/.matplotlib
echo "backend: TkAgg" >> ~/.matplotlib/matplotlibrc
```

---

## Additional Resources

### Documentation

- [Python Documentation](https://docs.python.org/3/)
- [Conda Documentation](https://docs.conda.io/en/latest/)
- [Homebrew Documentation](https://docs.brew.sh/)
- [PLIP Documentation](https://github.com/pharmai/plip)
- [PyMOL Wiki](http://pymol.wikidot.com/)

### Support

- **GitHub Issues:** [Report bugs or request features](https://github.com/ohatipoglu/foldx-analyzer/issues)
- **Discussions:** [Ask questions or share ideas](https://github.com/ohatipoglu/foldx-analyzer/discussions)

### Related Tools

- [FoldX Suite](http://foldxsuite.crg.eu/)
- [Biopython](https://biopython.org/)
- [OpenBabel](http://openbabel.org/)

---

## Quick Reference Commands

```bash
# Activate virtual environment
source venv/bin/activate

# Activate conda environment (Anaconda/Miniconda)
conda activate foldx-env

# Run application
python combined_gui.py

# Run tests
python -m pytest tests/ -v

# Check PyMOL installation
which pymol

# Check PLIP installation
plip --version

# Check Python version
python --version

# List installed packages
pip list

# Update dependencies
pip install --upgrade -r requirements.txt
```

---

## Anaconda Quick Start Guide

> **For macOS users who already have Anaconda installed** - Follow these condensed steps:

### 1. Open Terminal and navigate to project

```bash
cd ~/Projects/foldx-analyzer
```

### 2. Create and activate conda environment

```bash
conda create -n foldx-env python=3.12 -y
conda activate foldx-env
```

### 3. Install scientific packages from conda-forge

```bash
conda install -c conda-forge pandas numpy matplotlib biopython openbabel plip pymol pytest -y
```

### 4. Install remaining packages via pip

```bash
pip install customtkinter python-docx
```

### 5. Install the project

```bash
pip install -e .
```

### 6. Run the application

```bash
python combined_gui.py
```

### 7. (Optional) Run tests

```bash
python -m pytest tests/ -v
```

> **Note:** Your environment name `foldx-env` can be changed to any name you prefer.

---

## Version Compatibility Matrix

| Component | Minimum Version | Recommended Version | Notes |
|-----------|-----------------|---------------------|-------|
| macOS | 10.15 (Catalina) | 12.0 (Monterey) or later | Older versions may have compatibility issues |
| Python | 3.12 | 3.12.x | Python 3.13+ not yet tested |
| PyMOL | 2.4 | 2.5 or later | Required for PLIP visualization |
| PLIP | 1.4.0 | Latest | Requires OpenBabel |
| OpenBabel | 3.1.0 | Latest | Required for PLIP |
| customtkinter | 5.0 | Latest | GUI framework |

---

## Uninstallation

To remove the application and all dependencies:

1. **Deactivate and remove virtual environment**:

   ```bash
   deactivate
   rm -rf ~/Projects/foldx-analyzer/venv
   ```

2. **Remove project directory**:

   ```bash
   rm -rf ~/Projects/foldx-analyzer
   ```

3. **Remove conda environment** (if created):

   ```bash
   conda env remove -n foldx-env
   ```

4. **Remove globally installed packages** (optional):

   ```bash
   pip uninstall foldx-analyzer plip biopython customtkinter
   ```

---

## Changelog

### Version 0.2.0 (Current)

- macOS compatibility improvements
- Enhanced error handling for file operations
- Updated documentation

### Version 0.1.0

- Initial release

---

**Last Updated:** May 3, 2026

**Maintained By:** The FoldX Analyzer Team
