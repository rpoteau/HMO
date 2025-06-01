<div style="text-align:center">
<img src="https://raw.githubusercontent.com/rpoteau/HMO/main/hmo/icons-logos-banner/HMO_Banner.png" alt="HMObanner" width="800"/>
</div>

# Changelog

## Versions 0.7.6 - 2025-06-01
### Changed
- vertical flip of the molecule in `show_charge_bond_view()`, to match the orientation used in other graphical representations (e.g., MO diagrams and molecular overview), ensuring visual consistency across views

---

## Versions 0.7.4 & 0.7.5 - 2025-06-01
### Added
- added a **Using the Molecule Builder:** at the top of `HMO.py`, describing drawing rules, atom customization, and charge editing.

### Changed
- recompiled the documentation to reflect the new usage instructions.

---

## Version 0.7.3 - 2025-06-01
### Added
- added detailed docstrings to the following methods for improved documentation and future Sphinx integration: `sanitize_filename()`, `create_toolbar()`, `load_icons()`, `create_button()`, `bind_shortcuts()`, `toggle_eraser()`, `quit_program()`
- added confirmation dialog when quitting the application (`quit_program()` now prompts the user with a Yes/No messagebox before exiting)
- introduced a `__last_update__` variable in `hmo/__init__.py` to store the most recent update date.
- the `show_about()` window now displays this date along with a clickable link to the online changelog on GitHub.
- updated the `push_pypi.sh` script to automatically update the `__last_update__` field with the current date whenever a new version is pushed.

---

## Versions 0.7.1 & 0.7.2 - 2025-06-01
### Changed
- refactored `pyproject.toml` to follow [PEP 621](https://peps.python.org/pep-0621/) and modern setuptools structure:
    - switched to declarative package discovery using `[tool.setuptools.packages.find]`.
    - ensured that all subpackages like hmo.Fonts.OpenSans are included by adding `__init__.py` files.
    - defined package-data for fonts, diagrams, and icons.
- updated build-system section to use setuptools>=61 with build_meta.

---

## Version 0.7.0 - 2025-06-01
### Changed
- multiple matplotlib figures (e.g., MO grid, energy diagram) can now be saved as separate pages in a single PDF file
- `push_pyPi.sh` script also updates the `__version__` variable in `hmo/__init.py__`
- atom and bond coloring in `render_double_panel_charges_bonds()` are now set to lighter tints using the new `lighten_color()` blending utility
- the dimension of the figure returned by `render_double_panel_charges_bonds()` as a two side-by-side panels representing molecular structure and bonding and charge descriptors has been set up to 26 cm x 13 cm
- readability of the code of `render_double_panel_charges_bonds()` improved by defining explicit variables (`widthPlot_cm`, `heightPlot_cm`, `widthPlot`,...)

### Added
- added automatic git commit and tagging in `push_pyPi.sh` when bumping the version in this release script. The code version and the git repository are now synchronized at each release.
- added a utility function `lighten_color()` to adjust color intensity by blending any base color with white, enabling the creation of lighter, fully opaque color variants for graphical elements.
- global charge(s) rendering in `render_double_panel_charges_bonds`: they are now displayed in the sigma skeleton plot with colored circles (red for negative, blue for positive)
- added `df_global_charges` to store user-defined formal charges placed on the molecular canvas. Each entry includes charge value and (x, y) position in grid units.
- integrated `df_global_charges` into `save_dataframe_as_xlsx()` for export to Excel along with other molecular data.
- passed `df_global_charges` as an argument to HMOViewer for further integration and visualization support.
- systematically rendered global formal charges in `draw_skeleton_overview()`, with each charge displayed as a colored circle (red if negative, blue if positive) with a black outline.
- global formal charges are now also drawn in `show_dataframe_in_window()`, overlaying the molecular sketch in the top-left canvas. Same visual rules apply (color-coded circles with centered text).

### Fixed
- fixed a scaling issue in `render_double_panel_charges_bonds()` where the molecule could be incorrectly stretched to fill the panel. Bond lengths are now accurately rendered at the desired scale, unless the molecule is too large, in which case it is proportionally reduced to fit
- prevented a `TypeError` in `render_energy_diagram_matplotlib()` when rendering the energy diagram for molecules with undefined descriptors (e.g., no HOMO-LUMO gap). Descriptors like gap, η, and E_atom/N are now checked for numerical type before formatting. Non-numeric values are safely displayed as N/A.
- Ensured that `df_global_charge` is always created (even empty), avoiding missing attribute errors in export or drawing routines

---

## Versions 0.6.4 & 0.6.5 - 2025-05-22
### Changed
- Graphical documentation now also given in `README.md`.

## Version 0.6.3 - 2025-05-21

### Added
- new binding keys for Save data in a spreadsheet (`save_data`, Ctrl-l), About HMO (`show_about`, Ctrl-h), Export all results to PDF (`export_all_results_to_pdf`, Ctrl-P)

### Fixed
- (Ctrl-S) binding key was used both for saving the speadsheet and the tabular MOs.

---

## Version 0.6.2 - 2025-05-21

### Changed
- Read the Docs documentation update as graphical documentation support: Created detailed Inkscape-ready annotations for main GUI buttons and canvas.

---

## Version 0.6.1 - 2025-05-21

### Changed
- New "Université de Toulouse" logo, that replaces both the "UPS/UT3" and the former "Université de Toulouse" logos.

---

## Version 0.6.0 - 2025-05-20

### Added
New rendering functions, introduced as a foundation for building a unified high-resolution PDF export
- `render_all_OMs_grid_with_dual_lobes()`: exports a consistent grid of MOs with dual lobe shading and energy labels.
- `render_energy_diagram_matplotlib()`: exports a scaled π-energy diagram with fixed physical dimensions (13 cm), spin arrows, Hund's rule, and descriptor section.
- `render_double_panel_charges_bonds()`: generates a double panel (sigma skeleton + bond indices/charges) with auto-scaled coordinates for visual export.

---

## Version 0.5.1 - 2025-05-15

### Changed
- Summary section now selectable and copyable (Text widget instead of Label).

### Added
- Added "Save MOs as a text file" feature: export molecular orbital energies, occupations, and coefficients in blocks of 6 MOs as formatted text.

### Fixed
- Fixed "Close" button visibility on smaller screens by adjusting window layout.
- `self.formal_charges.clear()`, `self.molecule_is_new = True`, `self.charge_mol = 0` introduced in `load_molecule` for a proper initialization
- in `HMO Diagram Viewer` the Save PNG window now appears above the HMO Diagram window

---

## [0.5.0] - 2025-05-13

### Added
- Confirmation prompt when clearing a molecule to prevent accidental deletion.

### Fixed
- Prevented Hückel run on an empty molecule: a warning is now shown if no atoms are defined.
- The "Clear" button now also removes all formal charges, not just the molecular skeleton.
- Closing behavior improved: all related windows (MO viewer, OM diagram, and charge/bond view) are now properly closed when the molecule is cleared.

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
