# macOS Setup Guide for Bioinformatics Analyzer Suite (Conda-Only)

This guide provides a streamlined, robust method for setting up the Bioinformatics Analyzer Suite on macOS using only the Conda package manager. This approach avoids common compilation issues by using pre-compiled binaries from the `conda-forge` channel.

---

## Table of Contents

1.  [Prerequisites](#prerequisites)
2.  [Step 1: Install Miniconda](#step-1-install-miniconda)
3.  [Step 2: Configure Conda (Required)](#step-2-configure-conda-required)
4.  [Step 3: Create and Activate the Conda Environment](#step-3-create-and-activate-the-conda-environment)
5.  [Step 4: Install All Dependencies](#step-4-install-all-dependencies)
6.  [Step 5: Install the Application](#step-5-install-the-application)
7.  [Step 6: Configure and Run the Application](#step-6-configure-and-run-the-application)
8.  [Troubleshooting](#troubleshooting)

---

## Prerequisites

Before you begin, ensure you have the following:

-   **macOS**: Version 11.0 (Big Sur) or later.
-   **Architecture**: Apple Silicon (M1/M2/M3) or Intel. The commands are compatible with both.
-   **Terminal**: The built-in Terminal app or a replacement like iTerm2.
-   **Git**: Required for cloning the repository. Can be installed with `xcode-select --install`.
-   **Disk Space**: At least 3GB of free disk space for Miniconda and the required packages.

---

## Step 1: Install Miniconda

Miniconda is a minimal installer for Conda. We will use Homebrew to install it, as it's the simplest method.

1.  **Install Homebrew** (if you don't have it):
    Open your terminal and paste the following command:
    ```bash
    /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
    ```
    Follow the on-screen instructions to complete the installation.

2.  **Install Miniconda using Homebrew**:
    ```bash
    brew install --cask miniconda
    ```

3.  **Initialize Conda for your shell**:
    After the installation, run the following command to set up Conda for use in your terminal.
    ```bash
    conda init "$(basename "${SHELL}")"
    ```
    For example, if you use `zsh` (the default on modern macOS), the command is `conda init zsh`.

4.  **Restart Your Terminal**:
    Close your current terminal window and open a new one. You should now see `(base)` prefixed to your terminal prompt, indicating that Conda is active.

5.  **Verify Installation**:
    ```bash
    conda --version
    ```
    This should print the installed Conda version.

---

## Step 2: Configure Conda (Required)

To comply with Anaconda's new Terms of Service and prevent potential channel errors, you must accept the terms for the main channels.

1.  **Open your Terminal**.

2.  **Run the following commands one by one**:
    ```bash
    conda config --set channel_priority strict
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r
    ```
    This ensures that Conda can access the necessary packages without interruption.

---

## Step 3: Create and Activate the Conda Environment

We will create a dedicated Conda environment for the project to keep its dependencies isolated.

1.  **Clone the Repository** (if you haven't already):
    ```bash
    git clone https://github.com/ohatipoglu/foldx_analyzer.git
    cd foldx_analyzer
    ```

2.  **Create a new Conda environment**:
    This command creates an environment named `foldx-env` with Python 3.12.
    ```bash
    conda create -n foldx-env python=3.12 -y
    ```

3.  **Activate the environment**:
    ```bash
    conda activate foldx-env
    ```
    Your terminal prompt will now be prefixed with `(foldx-env)`. All subsequent commands will execute within this isolated environment.

---

## Step 4: Install All Dependencies

This is the core step. A single command will install all complex scientific packages from the reliable `conda-forge` channel.

1.  **Ensure your `foldx-env` is active**.

2.  **Install bioinformatics and core packages**:
    This command fetches `pandas`, `numpy`, `matplotlib`, `biopython`, `openbabel`, `plip`, and `pymol-open-source` all at once.
    ```bash
    conda install -c conda-forge pandas numpy matplotlib biopython openbabel plip pymol-open-source pytest -y
    ```
    > **Note:** The PyMOL package available on conda-forge is named `pymol-open-source`.

3.  **Install remaining packages with pip**:
    The following packages are best installed using `pip` within the Conda environment.
    ```bash
    pip install customtkinter python-docx
    ```

---

## Step 5: Install the Application

Now that all dependencies are in place, install the Bioinformatics Analyzer Suite itself in "editable" mode. This allows you to get code updates by simply running `git pull`.

```bash
pip install -e .
```

---

## Step 6: Configure and Run the Application

1.  **Run the application for the first time**:
    ```bash
    python foldx_analyzer_core/combined_gui.py
    ```
    The application will launch and automatically create a `config.ini` file in the project's root directory.

2.  **Verify PyMOL Path** (usually not needed with this method):
    The Conda installation automatically places `pymol` in the environment's PATH, so the default `config.ini` setting should work. You can verify the path with:
    ```bash
    which pymol
    ```
    The output should be a path inside your `miniconda3/envs/foldx-env/bin/` directory. If the application ever reports it cannot find PyMOL, you can update the `pymol_executable_path` in `config.ini` with this path.

---

## Troubleshooting

-   **Error: `PyMOL executable not found`**
    -   **Cause**: The application cannot find the PyMOL executable.
    -   **Solution**: Activate the conda environment (`conda activate foldx-env`) and run `which pymol`. Copy the full output path and paste it into the `pymol_executable_path` setting in `config.ini`.

-   **GUI Does Not Render Correctly (e.g., blurry fonts on Retina displays)**
    -   **Cause**: This can sometimes be related to the Tk backend used by `customtkinter`.
    -   **Solution**: This is rare with Conda installs, but if it occurs, ensure your macOS is up to date. No further action is usually required as Conda provides a compatible Tcl/Tk.

-   **SSL Certificate Errors during `pip install`**
    -   **Cause**: Python cannot find the root SSL certificates.
    -   **Solution**: Conda environments are generally immune to this, but if it happens, it can often be fixed by reinstalling `certifi`:
        ```bash
        pip install --upgrade --force-reinstall certifi
        ```
