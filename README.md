<div style="text-align:center">
<img src="https://raw.githubusercontent.com/rpoteau/HMO/main/hmo/icons-logos-banner/HMO_Banner.png" alt="HMObanner" width="800"/>
</div>

# Hückel Molecular Orbital (HMO) Viewer & Drawer

**HMO.py** is an interactive application that allows you to draw simple organic molecules, build the corresponding Hückel matrix, calculate π molecular orbitals, and visualize the numerical results as tables and graphics.

## Main Features

- ✏️ **Interactive drawing** of molecules (adding atoms and bonds),
- ⚙️ **Automatic calculation** of π molecular orbitals using the Hückel method,
- 📊 **Display of molecular orbital coefficients**, energies, bond indices, and π atomic charges,
- 🔬 **Calculation of molecular descriptors**:
  - Total π-electron energy,
  - HOMO-LUMO gap,
  - Electronic potential (μ),
  - Chemical hardness (η),
  - Chemical softness (S),
  - Electrophilicity index (ω),
- 💾 **Save/load molecules** in `.hmo` format,
- 📈 **Export of data** to Excel (.xlsx) files,
- 🖼️ **Graphical display of the molecule** in the results window.

## Installation

The required dependencies are listed in the `pyproject.toml` file. If you install the project via `pip install hmo`, all dependencies will be handled automatically.

## Launch

```bash
python HMO.py
```

## Usage

1️⃣ **Drawing:**
- Click to add atoms on the canvas,
- Click and drag between two atoms to create a bond,
- Delete or modify atoms/bonds as needed.

2️⃣ **Run Hückel:**
- Once the molecule is built, click the **Run Hückel** button to perform the calculations.

3️⃣ **Results:**
- Click **Show Numerical Results** to open the detailed window:
  - Molecule displayed with atom numbering,
  - Full table of molecular orbitals and descriptors,
  - Scrollbars are available for larger molecules.

4️⃣ **Saving:**
- Save the molecule using the menu or the dedicated button (in `.hmo` format),
- Export Excel results after running the calculation.

## `.hmo` File Format

A simple text file structured as follows:

```
Nodes:
C· 330.0 180.0
C· 480.0 210.0
...
Bonds:
0 1
1 2
...
```
- Each atom: `Type x y`
- Each bond: `index_atom1 index_atom2`

## Known Limitations

- The method is limited to planar π-conjugated systems,
- Atom names must be recognized by the library (e.g., C·, O:, etc.),
- No graphical export into the Excel file yet (planned for a future version).

## Roadmap

- Improved color and style management,
- Support for user-defined Hückel parameters,
- Optional integration of molecule snapshots in the Excel file.

## Useful Links

- 📚 [Documentation](https://hmo.readthedocs.io/)
- 💻 [GitHub Repository](https://github.com/rpoteau/HMO)
- 📝 [Changelog](https://github.com/rpoteau/HMO/blob/main/CHANGELOG.md)

## Credits

This project was developed by Romuald Poteau (LPCNO, University of Toulouse, CNRS & INSA), with strong contributions of ChatGPT for documentation, the graphical interface, and calculation features.

Thanks to users and testers for their feedback and improvement suggestions.

## Bibliography

- **Erich Hückel**, Quantentheoretische Beiträge zum Benzolproblem. I. Die Elektronenkonfiguration des Benzols und verwandter Verbindungen (**1931**), *Z. Phys.* **70**: 204–286. DOI: [10.1007/BF01339530](https://doi.org/10.1007/BF01339530)

- **Andrew Streitwieser**, Molecular Orbital Theory for Organic Chemists (**1961**), Wiley International Edition, 4th Edition, Wiley.

- **Robert Burns Woodward** & **Roald Hoffmann**, The Conservation of Orbital Symmetry (**1969**) *Angew. Chem., Int. ed. Eng.* **8**: 781-853. DOI: [10.1002/anie.196907811](https://doi.org/10.1002/anie.196907811)

- **Frederic A. Van-Catledge**, A Pariser-Parr-Pople-based set of Hückel molecular orbital parameters (**1980**), *J. Org. Chem.* **45**: 4801–4802. DOI: [10.1021/jo01311a060](https://doi.org/10.1021/jo01311a060)

