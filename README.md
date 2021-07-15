## Shared library with functions for analysis and the steering macro

### Combines previous macros for reading different objects from the simulations

## Installation

```bash
$ cd ~/alice
$ git clone https://github.com/MFT-MCHMatching/MFTana.git
$ cd MFTana
$ make
$ make install
```

## Usage

### Step 1: Process MFT simulation 
The first step loads and processes simulated and reconstructed MFT objects 
```bash
root.exe -b -q ~/alice/MFTana/MFTAna.C+
```
and creates `MFTAnaSimTracks.root` with extended MFT standalone analysis objects.

To browse the ROOT file:

```bash
root.exe _rootlogon.C MFTAnaSimTracks.root
```

### Step 2: Analyse `MFTAnaSimTracks.root` and generate histograms
Generates histograms from processed data from step 1.
```bash
root.exe -b -q ~/alice/MFTana/ReadMFTAnaSim.C+
```
Histograms stored on `ReadMFTAnaSim.root`

### Doxygen documentation

https://bovulpes.gitlab.io/mftanasim/files.html