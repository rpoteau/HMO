<div style="text-align:center">
<img src="https://raw.githubusercontent.com/rpoteau/HMO/main/hmo/icons-logos-banner/HMO_Banner.png" alt="HMObanner" width="800"/>
</div>

# Changelog

All notable changes to this project will be documented in this file.

## [0.4.1] - 2025-05-13

### Fixed
- Prevented Hückel run on an empty molecule: a warning is now shown if no atoms are defined.
- The "Clear" button now also removes all formal charges, not just the molecular skeleton.
- Closing behavior improved: all related windows (MO viewer, OM diagram, and charge/bond view) are now properly closed when the molecule is cleared.

### Added
- Confirmation prompt when clearing a molecule to prevent accidental deletion.

---

## [0.4.0] - 2025-05-11

### Added
- New `show_charge_bond_view()` method and `descriptors` button for visualizing molecular structures with atom types, formal charges, and π bond indices using matplotlib.
- Added GUI window with export options (PNG, SVG, PDF, EPS) and enhanced styling for atoms and labels.

---

## [0.3.1] - 2025-05-08

### Fixed
- when all MOs are occupied (eg, H2N-CH-NH2 radical), there is no LUMO within HMO. Error when displaying or saving HOMO-LUMO gaps and reactivity descriptors (eta,..) fixed 

---

## [0.3.0] - 2025-05-07

### [Added] Support for formal molecular charges
- Added a "Charge" button to place a non-bonding formal charge on the canvas.
- Charges are displayed as colored disks (red = negative, blue = positive) with the symbol centered.
- Right-click on a charge opens a dialog to edit its value (`+`, `-2`, `1`, etc.).
- Only one formal charge can exist at a time.
- New variable `self.charge_mol` is used in the total π-electron count:
  ```python
  total_pi_electrons -= self.charge_mol
- Imported `simpledialog` from `tkinter`.
- Function `add_plus_if_needed(text)` to display `+x` when appropriate on charge labels.
- Optional coordinates `x`, `y` for `draw()` method in `ChargeNode`, allowing custom positioning.s

### Changed
- Default `charge` parameter in `ChargeNode.__init__` changed from `'-'` to `'-1'`.
- Charge text color is now always set to white, instead of being conditional on its sign.
- Charge drawing now uses custom coordinates (`draw_x`, `draw_y`) instead of `self.x`, `self.y`.
- Attribute `charges` replaced by `formal_charges`.

---

## [0.2.0] - 2025-05-06

### Added
- Full English `README.md` and structured documentation (`docstrings` updated).
- Support for building a standalone Linux application using PyInstaller (with an appropriately configured build script).
- Custom "About" window with logo, author details, and version info; Escape key shortcut and non-resizable design.

### Fixed
- Various bug fixes in the molecule viewer and numerical results window (scrollbars, scaling).

---

## [0.1.5] to [0.1.11] - 2025-05-05

### Added
- `pyproject.toml`, `.gitignore`, `MANIFEST.in` configuration files.
- `push_pyPi.sh` script for automated deployment to PyPI.
- Documentation setup: `index.rst`, `conf.py` for ReadTheDocs.

---

## [0.1.1] to [0.1.4] - 2025-05-05

### Added
- `HuckelParameters` and `HMOViewerParameters` classes.
- Various internal bug fixes and code cleanups.

---

## [0.1.0] - 2025-05-04

### Major change
- Initial merge of `HMOViewer` and `MoleculeDrawer` into a unified project.

---

## MoleculeDrawer [0.0.1] to [0.0.18] & HMOViewer [0.0.1] to [0.0.7] 2025-05-01 - 2025-05-03
    
### Development
- Initial development of `HMOViewer` and `MoleculeDrawer` with iterative improvements alongside ChatGPT.
