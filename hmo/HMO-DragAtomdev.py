"""
Hückel Molecule Drawer & Viewer (HMO Tool)

This program provides a complete graphical interface for the interactive construction of
(planar π-conjugated) molecules and the analysis of their molecular orbitals using the Hückel model.
It allows the user to draw the molecular skeleton on a grid, select atom types (C·, N:, O·, etc.),
add bonds, and perform a full Hückel analysis that includes:

- generation of the Hückel matrix,
- solving for orbital energies (eigenvalues),
- calculating the molecular orbital coefficients (eigenvectors),
- determining π-charges on atoms and π-bond orders,
- computing global descriptors: total energy, HOMO-LUMO gap, hardness, chemical potential, etc.

Main features:
--------------
- Manual construction of the molecular skeleton (adding/removing atoms and bonds),
- Saving and loading the sigma skeleton of molecules (.hmo format),
- Running Hückel analysis and directly visualizing results,
- Exporting complete results (MO coefficients, π-charges, bond orders, descriptors) to Excel or PDF,
- Graphical visualization of molecular orbitals (shapes and energies) in a dedicated window.

Display organization (Tkinter-based):
-------------------------------------
The program is built with the Tkinter library, which manages all graphical interfaces:

1. **Molecule Drawer (main window)**:
    - Handled by the `MoleculeDrawer` class, this main window (based on `Tk()`) contains:
    - Left section: a `Canvas` where the user builds the molecule (grid, atoms, bonds),
    - Right section: a toolbar with buttons for main actions (run analysis, save/load molecule,
       access visualizations, export data).

2. **HMO Viewer (dedicated MO visualization window)**:
    - Instantiated by the `HMOViewer` class when the user clicks the visualization button.
    - This secondary window uses `Toplevel` to open independently.
    - It displays:
        - an energy diagram (MO energy levels),
        - graphical representations of occupied and virtual orbitals,
        - a miniature view of the molecule to show its overall topology.

3. **Numerical Results Window**:
    - Managed by the `ResultsViewer` class and also opened via `Toplevel`.
    - Displays:
        - molecular orbital coefficients in a table format,
        - π-charges, π-bond orders,
        - global descriptors (total energy, HOMO-LUMO gap, etc.).

Usage:
------
The user starts by constructing the molecule, then runs the Hückel analysis. After computation,
results can be explored numerically (`ResultsViewer` window) and graphically (`HMOViewer` window),
and exported. The graphical visualization allows easy exploration of molecular orbitals and their energy levels.

Notes:
------
- Most graphical interactions (atom placement, bond drawing) are done with the mouse,
  with keyboard shortcuts for efficiency.
- The program is structured around the main classes: `MoleculeDrawer` (main interface),
  `HMOViewer` (MO visualization), `ResultsViewer` (numerical results), and `Node` (management of individual atoms).
- Hückel parameters and some visulization settings are defined in the `HuckelParameters` class
- Visualization and layout parameters for the HMOViewer interface are provided by the `HMOViewerParameters` class
"""

# ========================================================================================================

import tkinter as tk
from tkinter import filedialog, messagebox, font, ttk, simpledialog
import matplotlib
matplotlib.use("TkAgg") ### mandatory for building a standalone Linux application
import matplotlib.pyplot as plt
from PIL import Image, ImageTk
import numpy as np
import pandas as pd
from pandas.plotting import table
import re
import sys, os
from pathlib import Path

# Helper function (global)
def resource_path(relative_path):
    """Get absolute path to resource, works for dev and PyInstaller."""
    try:
        base_path = sys._MEIPASS
    except AttributeError:
        base_path = os.path.abspath(".")
    return os.path.join(base_path, relative_path)

# ========================================================================================================


class HuckelParameters:
    """
    Container for Huckel model parameters and visualization settings.

    This class centralizes all constants and mappings used for building and analyzing molecules 
    within the Hückel framework. It includes atomic and bond parameters, grid and drawing settings, 
    and atom display options for the GUI.

    Attributes
    ----------
    GRID_SIZE : int
        The size of one grid unit (in pixels) for positioning atoms in the molecular canvas.
    ATOM_RADIUS : int
        The radius (in pixels) of the circle used to draw an atom on the canvas.
    HIGHLIGHT_RADIUS : int
        The radius (in pixels) of the highlight circle used when selecting atoms or bonds.
    ATOM_COLORS : dict
        A dictionary mapping atom type strings (e.g., 'C·', 'O:', 'N+·') to their display color 
        in hex code or standard color names.
    ATOM_OPTIONS : list of str
        The list of atom types available for selection and drawing in the molecule builder.
    Huckel_atomic_parameters : dict
        A dictionary mapping atom types to their Hückel parameters, each of which is a dictionary with:
            - 'alpha_expr': a string expression for the Coulomb integral α (often depends on β),
            - 'Fx': the index of free valence (reactivity index) representing the reactivity 
              of the center,
            - 'n_pi': the number of π-electrons contributed by the atom.
    Huckel_kXY_parameters : dict
        A nested dictionary mapping pairs of atom types to their interaction energies (float),
        representing the interaction between two 2p atomic orbitals (used for the β term in the 
        Hückel matrix).

    Notes
    -----
    - The atom types include common π-conjugated atoms and heteroatoms (e.g., C·, N·, N:, O·, etc.).
    - The ▯ symbol (e.g., 'B▯') indicates that the atom holds a π-vacancy (e.g., boron with an empty 
      p-orbital).
    - A colon (:) after the atom type (e.g., 'O:', 'N:') denotes that the atom provides a lone pair 
      into the π system (typically 2 π-electrons).
    - A dot (·) after the atom type (e.g., 'C·', 'O·', 'N·') indicates the atom provides a single 
      π-electron to the system.
    - The 'Fx' parameter (index of free valence) represents the theoretical reactivity of the atom 
      in the π system and can be used as a qualitative indicator of chemical reactivity.
    - The `Huckel_kXY_parameters` dictionary does not hold mere bond scaling factors, but rather 
      contains **interaction energies between two 2p atomic orbitals**, which define the off-diagonal 
      β terms in the Hückel Hamiltonian.
    - The `alpha_expr` parameter (Coulomb integral) gives the **energy of an electron in a 2p atomic orbital** 
      of a given atom.
    - Typically, α and β are **negative values**. For carbon atoms, it is common to use:
        - α_C ≈ -11.0 eV (Coulomb integral / 2p orbital energy),
        - β_CC ≈ -2.7 eV (interaction between two carbon 2p orbitals).
    - For heteroatoms (atoms other than carbon), the parameters α_X and β_XY are expressed as **functions 
      of α_C and β_C**, enabling a transferability of parameters while keeping carbon as the reference.

    Example
    -------
    >>> # Access a color for carbon
    >>> HuckelParameters.ATOM_COLORS['C·']
    '#909090'

    >>> # Access π-electron count for nitrogen (lone pair donor)
    >>> HuckelParameters.Huckel_atomic_parameters['N:']['n_pi']
    2

    >>> # Access 2p–2p interaction energy between C and O
    >>> HuckelParameters.Huckel_kXY_parameters['C·']['O·']
    1.06
    """

    GRID_SIZE: int = 30
    ATOM_RADIUS: int = 14
    HIGHLIGHT_RADIUS: int = 12

    ATOM_COLORS: dict = {
        'C·':   '#909090',
        'B▯':   '#ffb5b5',
        'N·':   'blue',
        'N:':   'navy',
        'N+·':  'deepskyblue',
        'O+·':  'tomato',
        'O·':   'red',
        'O:':   'darkred',
        'F:':   '#90e050',
        'Si·':  '#f0c8a0',
        'P·':   '#ffae4c',
        'P:':   '#ff8000',
        'S·':   '#ffff30',
        'S:':   '#afaf21',
        'Cl:':  '#afaf21',
        'Me:':  'purple',
        'Br:':  '#ff00ff',
    }

    ATOM_OPTIONS: list = [
        'B▯', 'C·', 'N·', 'N:', 'N+·', 'O·', 'O:',  'O+·', 'F:',
        'Si·', 'P·', 'P:', 'S·', 'S:', 'Cl:', 'Br:',
        'Me:'
    ]

    Huckel_atomic_parameters: dict = {
    'B▯':   {"alpha_expr": "alpha - 0.45 * beta", 'Fx': 1.705, 'n_pi': 0},
    'C·':   {"alpha_expr": "alpha + 0.00 * beta", 'Fx': 1.732, 'n_pi': 1},
    'N·':   {"alpha_expr": "alpha + 0.51 * beta", 'Fx': 1.393, 'n_pi': 1},
    'N:':   {"alpha_expr": "alpha + 1.37 * beta", 'Fx': 1.583, 'n_pi': 2},
    'N+·':  {"alpha_expr": "alpha + 2.00 * beta", 'Fx': None,   'n_pi': 1},
    'O·':   {"alpha_expr": "alpha + 0.97 * beta", 'Fx': 0.909, 'n_pi': 1},
    'O:':   {"alpha_expr": "alpha + 2.09 * beta", 'Fx': 0.942, 'n_pi': 2},
    'O+·':  {"alpha_expr": "alpha + 2.50 * beta", 'Fx': None,   'n_pi': 1},
    'F:':   {"alpha_expr": "alpha + 2.71 * beta", 'Fx': 0.179, 'n_pi': 2},
    'Si·':  {"alpha_expr": "alpha + 0.00 * beta", 'Fx': 1.732, 'n_pi': 1},
    'P·':   {"alpha_expr": "alpha + 0.19 * beta", 'Fx': 1.409, 'n_pi': 1},
    'P:':   {"alpha_expr": "alpha + 0.75 * beta", 'Fx': 1.666, 'n_pi': 2},
    'S·':   {"alpha_expr": "alpha + 0.46 * beta", 'Fx': 0.962, 'n_pi': 1},
    'S:':   {"alpha_expr": "alpha + 1.11 * beta", 'Fx': 1.229, 'n_pi': 2},
    'Cl:':  {"alpha_expr": "alpha + 1.48 * beta", 'Fx': 0.321, 'n_pi': 2},
    'Br:':  {"alpha_expr": "alpha + 1.50 * beta", 'Fx': None,   'n_pi': 2},
    'Me:':  {"alpha_expr": "alpha + 2.00 * beta", 'Fx': None,   'n_pi': 2},
    }

    Huckel_kXY_parameters: dict = {
    'C·':  {'C·':1.00, 'B▯':0.73, 'N·':1.02, 'N:':0.89, 'O·':1.06, 'O:':0.66, 'F:':0.52, 'Si·':0.75, 'P·':0.77, 'P:':0.76, 'S·':0.81, 'S:':0.69, 'Cl:':0.62, 'Me:':0.70, 'Br:':0.30, 'N+·':1.00, 'O+·':1.00},
    'B▯':  {'B▯':0.87, 'N·':0.66, 'N:':0.53, 'O·':0.60, 'O:':0.35, 'F:':0.26, 'Si·':0.57, 'P·':0.53, 'P:':0.54, 'S·':0.51, 'S:':0.44, 'Cl:':0.41, 'Me:':None, 'Br:':None, 'N+·':None, 'O+·':None},
    'N·':  {'N·':1.09, 'N:':0.99, 'O·':1.14, 'O:':0.80, 'F:':0.65, 'Si·':0.72, 'P·':0.78, 'P:':0.81, 'S·':0.83, 'S:':0.78, 'Cl:':0.77, 'Me:':None, 'Br:':None, 'N+·':None, 'O+·':None},
    'N:':  {'N:':0.98, 'O·':1.13, 'O:':0.89, 'F:':0.77, 'Si·':0.43, 'P·':0.55, 'P:':0.64, 'S·':0.68, 'S:':0.73, 'Cl:':0.80, 'Me:':None, 'Br:':None, 'N+·':None, 'O+·':None},
    'O·':  {'O·':1.26, 'O:':1.02, 'F:':0.92, 'Si·':0.65, 'P·':0.75, 'P:':0.82, 'S·':0.84, 'S:':0.85, 'Cl:':0.88, 'Me:':None, 'Br:':None, 'N+·':None, 'O+·':None},
    'O:':  {'O:':0.95, 'F:':0.94, 'Si·':0.24, 'P·':0.31, 'P:':0.39, 'S·':0.43, 'S:':0.54, 'Cl:':0.70, 'Me:':None, 'Br:':None, 'N+·':None, 'O+·':None},
    'F:':  {'F:':1.04, 'Si·':0.17, 'P·':0.21, 'P:':0.22, 'S·':0.28, 'S:':0.32, 'Cl:':0.51, 'Me:':None, 'Br:':None, 'N+·':None, 'O+·':None},
    'Si·': {'Si·':0.64, 'P·':0.62, 'P:':0.52, 'S·':0.61, 'S:':0.40, 'Cl:':0.34, 'Me:':None, 'Br:':None, 'N+·':None, 'O+·':None},
    'P·':  {'P·':0.63, 'P:':0.58, 'S·':0.65, 'S:':0.48, 'Cl:':0.35, 'Me:':None, 'Br:':None, 'N+·':None, 'O+·':None},
    'P:':  {'P:':0.63, 'S·':0.65, 'S:':0.60, 'Cl:':0.55, 'Me:':None, 'Br:':None, 'N+·':None, 'O+·':None},
    'S·':  {'S·':0.68, 'S:':0.58, 'Cl:':0.52, 'Me:':None, 'Br:':None, 'N+·':None, 'O+·':None},
    'S:':  {'S:':0.63, 'Cl:':0.59, 'Me:':None, 'Br:':None, 'N+·':None, 'O+·':None},
    'Cl:': {'Cl:':0.68, 'Me:':None, 'Br:':None, 'N+·':None, 'O+·':None},
    'Me:': {'Me:':None, 'Br:':None, 'N+·':None, 'O+·':None},
    'Br:': {'Br:':None, 'N+·':None, 'O+·':None},
    'N+·': {'N+·':None, 'O+·':None},
    'O+·': {'O+·':None}
    }


# =============================================================================================================================================
class HMOViewerParameters:
    """
    Visualization and layout parameters for the HMOViewer interface.

    This class centralizes all GUI layout settings, dimensions, and rendering 
    parameters used in the molecular orbital (MO) diagram viewer. It includes 
    frame sizes and positions for HOMO/LUMO, diagram coordinates, canvas size, 
    and parameters for drawing molecular orbitals (lobes, bonds).

    Attributes
    ----------
    FRAME_WIDTH : int
        Width (in pixels) of each MO visualization frame (for HOMO/LUMO).
    FRAME_HEIGHT : int
        Height (in pixels) of each MO visualization frame.
    FRAME_LUMO_X : int
        X-coordinate (in pixels) of the LUMO frame's top-left corner.
    FRAME_LUMO_Y : int
        Y-coordinate (in pixels) of the LUMO frame's top-left corner.
    FRAME_HOMO_X : int
        X-coordinate of the HOMO frame's top-left corner (aligned with LUMO frame).
    FRAME_HOMO_Y : int
        Y-coordinate of the HOMO frame's top-left corner (placed below LUMO with spacing).
    BOTTOM_OF_HOMO : int
        Y-coordinate of the bottom of the HOMO frame.
    MARGIN_BOTTOM : int
        Margin (in pixels) below the HOMO frame in the canvas.
    CANVAS_HEIGHT : int
        Total height (in pixels) of the drawing canvas, computed dynamically.
    CANVAS_WIDTH : int
        Total width (in pixels) of the drawing canvas.
    DIAG_X : int
        X-position (in pixels) for the center of the MO energy diagram (left side).
    DIAG_Y_TOP : int
        Top Y-coordinate of the MO energy diagram.
    DIAG_Y_BOTTOM : int
        Bottom Y-coordinate of the MO energy diagram.
    DIAG_WIDTH_UNIT : int
        Width of a unit interval in the diagram (used for scaling purposes).
    SHRINK_FACTOR : float
        Scaling factor applied to lobe sizes to prevent overlap.
    LOBE_OFFSET : int
        Pixel offset for lobe positioning when rendering orbitals.
    TARGET_BOND_PX : int
        Target length (in pixels) for bond drawings between atoms.

    Notes
    -----
    - The canvas is designed to hold both HOMO and LUMO visualizations side by side with appropriate margins.
    - The MO energy diagram is drawn to the left of the MO frames for visual reference.
    - The `SHRINK_FACTOR` and `LOBE_OFFSET` parameters allow fine control of orbital lobe rendering for clarity.
    - The `Path2Imgs` attribute should point to a valid directory containing any additional images used in the visualization (e.g., backgrounds, decorations).

    Example
    -------
    >>> HMOViewerParameters.FRAME_WIDTH
    700

    >>> HMOViewerParameters.Path2Imgs
    PosixPath('DesignMOdiagram')
    """
    
    # Frames for MO visualization
    FRAME_WIDTH = 700
    FRAME_HEIGHT = 450
    
    FRAME_LUMO_X = 350
    FRAME_LUMO_Y = 30
    
    FRAME_HOMO_X = FRAME_LUMO_X
    FRAME_HOMO_Y = FRAME_LUMO_Y + FRAME_HEIGHT + 40  # vertical spacing of 40 px
    
    BOTTOM_OF_HOMO = FRAME_HOMO_Y + FRAME_HEIGHT
    MARGIN_BOTTOM = 50
    CANVAS_HEIGHT = BOTTOM_OF_HOMO + MARGIN_BOTTOM
    
    CANVAS_WIDTH = FRAME_LUMO_X + FRAME_WIDTH + 40 # vertical spacing of 40 px
    
    # MOD Diagram on the left
    DIAG_X = 120
    DIAG_Y_TOP = 50
    DIAG_Y_BOTTOM = 600
    DIAG_WIDTH_UNIT = 30
    
    SHRINK_FACTOR = 0.85
    LOBE_OFFSET = 5
    TARGET_BOND_PX = 60

    # '''
    # Path2Imgs : Path
    # Path object pointing to the directory where design images or assets are stored.
    # '''
    # Path2Imgs = Path("DesignOfMOdiagram")

class HMOViewer:
    """
    GUI application for visualizing Hückel Molecular Orbitals (MOs) and their energy diagrams.
    
    This class provides a graphical interface to render the molecular structure and the associated energy levels of the molecular orbitals (occupied and virtual). It allows the user to explore the molecule's π-system, display orbital diagrams, adjust visualization parameters in real time, and save high-quality images.
    
    Parameters
    ----------
    master : tk.Toplevel or tk.Tk
        The parent Tkinter window that contains the interface.
    df_MOs : pd.DataFrame
        DataFrame containing molecular orbital coefficients for each atom in the system.
    df_atoms : pd.DataFrame
        DataFrame containing atomic positions, indices, and element types.
    df_bonds : pd.DataFrame
        DataFrame specifying bonds between atoms (indices and bond order).
    df_descriptors : pd.DataFrame
        DataFrame containing computed molecular descriptors (e.g., number of π electrons, symmetry).
    project_name : str
        Name of the current project, used in window titles and when saving files.
    
    Attributes
    ----------
    canvas : tk.Canvas
        The canvas widget where the molecule and energy levels are drawn.
    scale_slider : tk.Scale
        Slider widget that allows dynamic adjustment of the overall molecular diagram scale (×0.5 to ×3.0).
    lobe_slider : tk.Scale
        Slider widget that allows dynamic adjustment of orbital lobe sizes (×0.7 to ×1.3).
    homo_lumo_button : tk.Button
        Button that resets the view to the default HOMO and LUMO orbitals.
    skeleton_button : tk.Button
        Button to toggle skeleton-only mode (hide or show atom labels).
    save_button : tk.Button
        Button to save the full diagram as a high-quality PNG image.
    close_button : tk.Button
        Button to close the viewer window.
    show_atom_labels : bool
        Tracks whether atom labels are displayed in the skeleton overview.
    skeleton_items : list
        List of graphical items representing the skeleton; used for easy refresh.
    user_scale_multiplier : float
        Multiplier applied to the default molecule scaling (set by the scale_slider).
    user_lobe_multiplier : float
        Multiplier applied to the orbital lobe scaling (set by the lobe_slider).
    
    Notes
    -----
    - Energy levels can be clicked to view individual MOs.
    - The visualization is interactive: adjusting sliders immediately refreshes the displayed orbitals.
    - Pressing the Escape key closes the viewer.

    Example
    -------
    >>> root = tk.Tk()
    >>> viewer = HMOViewer(
    ...     master=root,
    ...     df_MOs=df_mos,
    ...     df_atoms=df_atoms,
    ...     df_bonds=df_bonds,
    ...     df_descriptors=df_descriptors,
    ...     project_name="Benzene"
    ... )
    >>> root.mainloop()
    
    - The window opens showing the molecular diagram and energy levels.
    - Click the "HOMO/LUMO" button to quickly reset the view.
    - Adjust the "Molecule Scale" and "Lobe Scale" sliders to change the diagram in real-time.
    - Save the diagram using the "Save as PNG" button.
    - Close the window by pressing Escape or clicking "Close".
    """


    def __init__(self, master, df_MOs, df_atoms, df_bonds, df_descriptors, project_name):
        self.master = master
        self.master.title("HMO Diagram Viewer")

        # Stocke les DataFrames directement
        self.df_MOs = df_MOs
        self.df_atoms = df_atoms
        self.df_bonds = df_bonds
        self.df_descriptors = df_descriptors
        self.project_name = project_name

        self.show_atom_labels = True  # État initial : on montre les atomes dans le squelette
        self.skeleton_items = []      # Liste pour stocker les éléments graphiques du squelette

        # Chargement des images
        self.load_images()

        # Initialisation des valeurs
        self.prepare_data()

        # Création du canvas
        self.canvas = tk.Canvas(master, width=HMOViewerParameters.CANVAS_WIDTH, height=HMOViewerParameters.CANVAS_HEIGHT, bg='white')
        self.canvas.pack()

        # Sliders
        self.user_scale_multiplier = 1.0
        self.user_lobe_multiplier = 1.0

        # Shortcuts
        self.master.bind('<Escape>', self.on_escape)

        # === 🟦 Ajout des boutons sous le diagramme ===
        button_frame = tk.Frame(self.master)
        button_frame.pack(pady=10)  # petit espacement vertical
        
        # Bouton HOMO/LUMO
        self.homo_lumo_button = tk.Button(button_frame, text="HOMO/LUMO", command=self.display_default_homo_lumo)
        self.homo_lumo_button.pack(side=tk.LEFT, padx=5)
        
        # === Slider pour la taille globale ===
        self.scale_slider = tk.Scale(button_frame, from_=0.5, to=3.0, resolution=0.1, orient=tk.HORIZONTAL, label="Molecule Scale")
        self.scale_slider.set(1.0)
        self.scale_slider.pack(side=tk.LEFT, padx=5)
        
        # === Slider pour la taille des lobes ===
        self.lobe_slider = tk.Scale(button_frame, from_=0.7, to=1.3, resolution=0.1, orient=tk.HORIZONTAL, label="Lobe Scale")
        self.lobe_slider.set(1.0)
        self.lobe_slider.pack(side=tk.LEFT, padx=5)

        # Bouton Skeleton Only
        self.skeleton_button = tk.Button(button_frame, text="Skeleton Only", command=self.toggle_skeleton)
        self.skeleton_button.pack(side=tk.LEFT, padx=5)
        
        # Bouton Save as PNG
        self.save_button = tk.Button(button_frame, text="Save as PNG", command=self.save_canvas_as_png)
        self.save_button.pack(side=tk.LEFT, padx=5)

        # Bouton Close
        self.close_button = tk.Button(button_frame, text="Close", command=self.master.destroy)
        self.close_button.pack(side=tk.LEFT, padx=5)
        
        self.scale_slider.config(command=lambda event: self.refresh_MOs())
        self.lobe_slider.config(command=lambda event: self.refresh_MOs())

        self.draw_layout()
        self.draw_energy_levels()
        self.display_default_homo_lumo()
        
    def on_escape(self, event):
        # print("Escape key pressed. Closing the app.")
        self.master.destroy()

    def load_images(self):
        e1_png = resource_path(f"DesignOfMOdiagram/1e.png")
        e2_png = resource_path(f"DesignOfMOdiagram/2e.png")
        energy_level_png = resource_path(f"DesignOfMOdiagram/energy_level.png")
        self.img_1e = ImageTk.PhotoImage(Image.open(e1_png).resize((8, 40)))
        self.img_2e = ImageTk.PhotoImage(Image.open(e2_png).resize((23, 41)))
        self.img_energy_level = ImageTk.PhotoImage(Image.open(energy_level_png).resize((40, 5)))

    def prepare_data(self):
    
        def extract_beta_coeff(energy_str):
            pattern = r'(.*?)β'
            match = re.search(pattern, energy_str)
            if match:
                coeff_str = match.group(1).replace('β', '').strip().split()
                if len(coeff_str) == 0:
                    return 1.0
                if coeff_str[0] == '+':
                    return float(coeff_str[1]) if len(coeff_str) > 1 else 1.0
                if coeff_str[0] == '-':
                    return -float(coeff_str[1]) if len(coeff_str) > 1 else -1.0
                return float(coeff_str[0])
            return 0.0
    
        self.energies = []
        self.occupations = []
    
        self.max_coef_global = np.nanmax(np.abs(self.df_MOs.values))
        coords = {row['Atom']: (row['X (grid units)'], row['Y (grid units)']) for _, row in self.df_atoms.iterrows()}
    
        bond_lengths = []
        for _, row in self.df_bonds.iterrows():
            a1, a2 = row['Atom 1'], row['Atom 2']
            x1, y1 = coords[a1]
            x2, y2 = coords[a2]
            dist = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)
            bond_lengths.append(dist)
        self.mean_bond_length = np.mean(bond_lengths)
    
        for col in self.df_MOs.columns:
            energy_expr = col.split('\n')[0].split('=')[-1].replace('α', '').strip()
            energy = extract_beta_coeff(energy_expr)
            self.energies.append(energy)
            occ_raw = col.split('\n')[1].strip()
            occ_num = occ_raw.split('e')[0]  # prend tout avant le 'e'
            self.occupations.append(int(occ_num))
    
        self.energy_groups = {}
        for idx, energy in enumerate(self.energies):
            rounded_e = round(energy, 5)
            self.energy_groups.setdefault(rounded_e, []).append(idx)

    def toggle_skeleton(self):
        """Active/désactive le mode skeleton only."""
        self.show_atom_labels = not self.show_atom_labels
        # print(f"[DEBUG] Skeleton only mode: {self.show_atom_labels}")
        self.refresh_skeleton()
        self.refresh_MOs()

    def refresh_MOs(self):
        """Redessine l'OM actuellement affichée (HOMO & LUMO)."""
        if hasattr(self, 'current_occ_idx'):
            self.display_om(self.current_occ_idx, occupied=True)
        if hasattr(self, 'current_virt_idx'):
            self.display_om(self.current_virt_idx, occupied=False)
    
    def refresh_skeleton(self):
        """Efface et redessine le squelette global."""
        # Effacer l'ancien squelette
        for item in self.skeleton_items:
            self.canvas.delete(item)
        self.skeleton_items = []
    
        # Redessiner à la position sauvegardée
        self.draw_skeleton_overview(
            self.skeleton_x0, self.skeleton_y0,
            self.skeleton_width, self.skeleton_height
        )

    def save_canvas_as_png(self):
        """Sauvegarde le canvas entier en PNG haute qualité."""
        from tkinter import filedialog
        file = filedialog.asksaveasfilename(
            defaultextension=".png",
            filetypes=[("PNG files", "*.png")],
            initialfile=f"{self.project_name}_HMO_diagram.png"
        )
        if file:
            # Sauvegarde du canvas en PostScript (format vectoriel temporaire)
            file_path = Path(file)
            tmp_eps = file_path.with_suffix('.tmp.eps')
            self.canvas.postscript(file=tmp_eps, colormode='color')
            try:
                from PIL import Image
                img = Image.open(tmp_eps)
    
                # 🔥 Améliorer la qualité : ouvrir à haute résolution (~300 DPI)
                # (default is ~72 dpi -> flou, 300 dpi = beaucoup mieux)
                img.load(scale=12)  # facteur qui augmente la résolution
    
                img.save(file, 'png')
                print(f"Canvas saved as {file}")
                # Nettoyer
                tmp_eps.unlink()
            except Exception as e:
                print(f"Erreur de sauvegarde PNG: {e}")

    def draw_layout(self):
        """
        Draws the main layout of the viewer interface.
    
        This includes:
        - The project name at the top.
        - Two frames: one for virtual MOs (LUMO and above) and one for occupied MOs (HOMO and below).
        - Labels for each section.
        - An overview area to display a miniaturized version of the molecular skeleton under the MO diagram.
    
        Notes
        -----
        The overview area is scaled to fit the available space and is drawn using `draw_skeleton_overview`.
        """
        self.canvas.create_text(
            3*HMOViewerParameters.FRAME_LUMO_X/4,
            10,
            text=self.project_name,
            font=('DejaVu Sans', 12, 'bold'),
            fill='#128d85',
            anchor='center'
        )

        self.canvas.create_rectangle(HMOViewerParameters.FRAME_LUMO_X, HMOViewerParameters.FRAME_LUMO_Y, HMOViewerParameters.FRAME_LUMO_X + HMOViewerParameters.FRAME_WIDTH, HMOViewerParameters.FRAME_LUMO_Y + HMOViewerParameters.FRAME_HEIGHT, outline='black')
        self.canvas.create_text(HMOViewerParameters.FRAME_LUMO_X + HMOViewerParameters.FRAME_WIDTH/2, HMOViewerParameters.FRAME_LUMO_Y - 20, text='Virtual MOs', font=('DejaVu Sans', 14, 'bold'))

        self.canvas.create_rectangle(HMOViewerParameters.FRAME_HOMO_X, HMOViewerParameters.FRAME_HOMO_Y, HMOViewerParameters.FRAME_HOMO_X + HMOViewerParameters.FRAME_WIDTH, HMOViewerParameters.FRAME_HOMO_Y + HMOViewerParameters.FRAME_HEIGHT, outline='black')
        self.canvas.create_text(HMOViewerParameters.FRAME_HOMO_X + HMOViewerParameters.FRAME_WIDTH/2, HMOViewerParameters.FRAME_HOMO_Y - 20, text='Occupied MOs', font=('DejaVu Sans', 14, 'bold'))
        descriptors_height = 140
        # print(f"DEBUG y0 = {HMOViewerParameters.CANVAS_HEIGHT - 10=}")
        # print(f"DEBUG width = {HMOViewerParameters.FRAME_LUMO_X - 10=}")
        # print(f"DEBUG height = {HMOViewerParameters.CANVAS_HEIGHT - 10 - HMOViewerParameters.DIAG_Y_BOTTOM - descriptors_height =}")
        self.skeleton_x0 = 10  # ou l'endroit où tu dessines le squelette
        self.skeleton_y0 = HMOViewerParameters.DIAG_Y_BOTTOM + descriptors_height  # sous le diagramme HOMO
        self.skeleton_width = HMOViewerParameters.FRAME_LUMO_X - 10
        self.skeleton_height = HMOViewerParameters.CANVAS_HEIGHT - 10 - HMOViewerParameters.DIAG_Y_BOTTOM - descriptors_height
        self.draw_skeleton_overview(
                                        x0 = self.skeleton_x0,  # ajuste en fonction de la largeur que tu veux
                                        y0 = self.skeleton_y0,  # un peu sous le diagramme
                                        width = self.skeleton_width,
                                        height = self.skeleton_height,
                                    )

    def draw_energy_scale_and_descriptors(self, min_e, max_e, scale, energies_sorted):
        """
        Draws the energy scale, horizontal dotted guide lines, and key descriptors under the MO diagram.
    
        Parameters
        ----------
        min_e : float
            The minimum energy value in the set of molecular orbitals.
        max_e : float
            The maximum energy value.
        scale : float
            The vertical scaling factor to convert energy values to pixel positions.
        energies_sorted : list
            Sorted list of energy values (from high to low).
    
        Notes
        -----
        - Adds a vertical energy scale with an arrow, horizontal guide lines at regular energy steps,
          and labels each line with its corresponding (α + β) expression.
        - Displays summary descriptors like total energy, atomization energy, HOMO-LUMO gap, and hardness.
        """
        # print(f"draw_energy_scale_and_descriptors {max_e=}")
        # print(f"draw_energy_scale_and_descriptors {min_e=}")
        # === 🟦 ÉCHELLE ÉNERGÉTIQUE VERTICALE ===
        start_y = HMOViewerParameters.DIAG_Y_BOTTOM - (0.3) * scale
        end_y = HMOViewerParameters.DIAG_Y_BOTTOM - ( (max_e - 0.2) - min_e ) * scale

        # Prendre un des groupes (le premier trouvé)
        for energy in energies_sorted:
            group = self.energy_groups[energy]
            n_degenerate = len(group)
            total_width = n_degenerate * (self.img_energy_level.width() + 10)
            start_x = HMOViewerParameters.DIAG_X + (self.img_energy_level.width() // 2) - total_width // 2 + (self.img_energy_level.width() // 2)
            # x_pos du premier dans le groupe :
            center_x = start_x
            break
            
        x_pos_scale = center_x
    
        self.canvas.create_line(
            x_pos_scale, start_y, x_pos_scale, end_y,
            fill='blue', width=2, arrow=tk.LAST, state="disabled"
        )
    
        # === 🟦 POINTILLÉS HORIZONTAUX ===
        max_deg = max(len(group) for group in self.energy_groups.values())
        line_length = (max_deg * (self.img_energy_level.width() + 10)) * 1.1
    
        step = 0.5
        current = step * np.floor(min_e / step)
        while current >= max_e:
            y_pos = HMOViewerParameters.DIAG_Y_BOTTOM - ( current - min_e ) * scale
        
            self.canvas.create_line(
                center_x - line_length/2, y_pos, center_x + line_length/2, y_pos,
                fill='blue', dash=(4, 2), state="disabled"
            )
        
            if np.isclose(current, 0.0):
                label = 'α'
            elif current > 0:
                label = f"α + {abs(current):.1f} β"
            else:
                label = f"α - {abs(current):.1f} β"
        
            self.canvas.create_text(
                center_x + line_length/2 + 20, y_pos,
                text=label, fill='blue', anchor='w',
                font=('DejaVu Sans', 10), state="disabled"
            )
        
            current -= step  # on descend en énergie

    
        # === 🟦 DESCRIPTEURS SOUS LE DIAGRAMME ===
        try:
            total_energy = self.df_descriptors.loc['Total π-electron energy'].values[0]
            atomization = self.df_descriptors.loc['Atomization energy per π atom'].values[0]
            gap = self.df_descriptors.loc['HOMO-LUMO gap'].values[0]
            hardness = self.df_descriptors.loc['η (hardness)'].values[0]
        except Exception as e:
            print(f"Erreur lecture descriptors: {e}")
            total_energy = atomization = gap = hardness = 'N/A'
    
        y_text = HMOViewerParameters.DIAG_Y_BOTTOM + 60
        x_center = center_x
    
        self.canvas.create_text(
            x_center, y_text,
            text=f"E: {total_energy}",
            font=('DejaVu Sans', 11), anchor='center', state="disabled"
        )
        self.canvas.create_text(
            x_center, y_text + 20,
            text=f"E_atomization per atom: {atomization}β",
            font=('DejaVu Sans', 11), anchor='center', state="disabled"
        )
        self.canvas.create_text(
            x_center, y_text + 40,
            text=f"HOMO-LUMO gap: {gap}|β|",
            font=('DejaVu Sans', 11), anchor='center', state="disabled"
        )
        self.canvas.create_text(
            x_center, y_text + 60,
            text=f"Hardness: {hardness}|β|",
            font=('DejaVu Sans', 11), anchor='center', state="disabled"
        )

    def draw_energy_levels(self):
        """
        Draws all energy levels (molecular orbitals) in the diagram.
    
        Each energy level:
        - Is represented as an image (usually a line or bar).
        - Is grouped if degenerate (multiple orbitals with same energy).
        - Shows electron occupancy via additional graphical indicators (2e or 1e overlays).
    
        Notes
        -----
        Binds a click event to each energy level image so that clicking an MO
        updates the orbital display (via `on_level_click`).
        """
        energies_sorted = sorted(self.energy_groups.keys(), reverse=True)
        # print("Sorted energies (β values):", energies_sorted)
        min_e, max_e = max(energies_sorted), min(energies_sorted)
        # print(f"draw_energy_levels {max_e=}")
        # print(f"draw_energy_levels {min_e=}")        
        scale = (HMOViewerParameters.DIAG_Y_BOTTOM - HMOViewerParameters.DIAG_Y_TOP) / (max_e - min_e + 1e-6)

        self.draw_energy_scale_and_descriptors(min_e, max_e, scale, energies_sorted)

        for energy in energies_sorted:
            group = self.energy_groups[energy]
            n_degenerate = len(group)
            total_width = n_degenerate * (self.img_energy_level.width() + 10)
            start_x = HMOViewerParameters.DIAG_X + (self.img_energy_level.width() // 2) - total_width // 2 + (self.img_energy_level.width() // 2)

            y_pos = HMOViewerParameters.DIAG_Y_BOTTOM - (energy - min_e) * scale

            for i, idx in enumerate(group):
                x_pos = start_x + i * (self.img_energy_level.width() + 10)

                img = self.canvas.create_image(x_pos, y_pos, image=self.img_energy_level)

                occ = self.occupations[idx]
                if occ == 2:
                    self.canvas.create_image(x_pos, y_pos, image=self.img_2e, state="disabled")
                elif occ == 1:
                    self.canvas.create_image(x_pos, y_pos, image=self.img_1e, state="disabled")

                self.canvas.tag_bind(img, '<Button-1>', lambda event, idx=idx: self.on_level_click(idx))
        
    def draw_skeleton_overview(self, x0, y0, width, height):
        """
        Draws a miniaturized version of the full molecule in a dedicated area under the MO diagram.
    
        Parameters
        ----------
        x0 : int
            The x-coordinate of the top-left corner of the drawing area.
        y0 : int
            The y-coordinate of the top-left corner.
        width : int
            The width of the available area.
        height : int
            The height of the available area.
    
        Notes
        -----
        - Automatically scales the molecule to fit the provided area.
        - Displays bonds and atoms; optional atom labels are shown if `self.show_atom_labels` is True.
        """
        coords = {row['Atom']: (row['X (grid units)'], row['Y (grid units)']) for _, row in self.df_atoms.iterrows()}
        colors = {row['Atom']: row.get('Color', 'gray') for _, row in self.df_atoms.iterrows()}  # <== récupère la couleur
        xs = [pos[0] for pos in coords.values()]
        ys = [pos[1] for pos in coords.values()]
        min_x, max_x = min(xs), max(xs)
        min_y, max_y = min(ys), max(ys)
    
        mol_width = max_x - min_x
        mol_height = max_y - min_y
        margin = 20
    
        scale_x = (width - margin) / mol_width if mol_width != 0 else 1
        scale_y = (height - margin) / mol_height if mol_height != 0 else 1
        scale = min(scale_x, scale_y)
    
        offset_x = x0 + width/2 - ((min_x + max_x) / 2) * scale
        offset_y = y0 + height/2 - ((min_y + max_y) / 2) * scale
    
        # Stocker les éléments créés
        items = []
    
        # tracer liaisons
        for _, row in self.df_bonds.iterrows():
            a1, a2 = row['Atom 1'], row['Atom 2']
            x1, y1 = coords[a1]
            x2, y2 = coords[a2]
            line = self.canvas.create_line(
                offset_x + x1 * scale, offset_y + y1 * scale,
                offset_x + x2 * scale, offset_y + y2 * scale,
                width=2
            )
            items.append(line)
    
        if self.show_atom_labels:
            radius_oval = 15
            for atom, (x, y) in coords.items():
                color = colors.get(atom, 'gray')
                oval = self.canvas.create_oval(
                    offset_x + x * scale - radius_oval, offset_y + y * scale - radius_oval,
                    offset_x + x * scale + radius_oval, offset_y + y * scale + radius_oval,
                    fill=color
                )
                label = self.canvas.create_text(
                    offset_x + x * scale, offset_y + y * scale,
                    text=atom, font=('DejaVu Sans', 10)
                )
                items.append(oval)
                items.append(label)
    
        # Sauvegarde pour pouvoir effacer ensuite
        self.skeleton_items.extend(items)

    def display_om(self, idx, occupied=True):
        """
        Displays a specific molecular orbital (MO) in its corresponding frame (occupied or virtual).
        
        This method draws the selected MO, including bonds and orbital lobes, in the appropriate frame (HOMO or LUMO area), scaled and centered within the available space. The visualization responds dynamically to user-controlled scaling factors from the sliders:
        - 'Molecule Scale' slider: adjusts the overall size of the molecule,
        - 'Lobe Scale' slider: adjusts the relative size of the orbital lobes.
        
        Parameters
        ----------
        idx : int
            The index of the MO to display.
        occupied : bool, default=True
            Whether the MO is an occupied orbital (True) or a virtual one (False).
        
        Notes
        -----
        - Centers and scales the molecule to fit inside its frame.
        - Draws bonds and orbital lobes, with size and color depending on MO coefficients.
        - Adds atom labels if 'Skeleton Only' mode is off.
        - The user-defined scaling factors (sliders) affect both molecule size and lobe size in real time.
        - Handles rendering of both positive (red) and negative (blue) lobes with layered drawing.
        """


        scale_user = self.scale_slider.get()
        lobe_user = self.lobe_slider.get()

        # === 🖼️ Cadre d'affichage ===
        frame_x = HMOViewerParameters.FRAME_HOMO_X if occupied else HMOViewerParameters.FRAME_LUMO_X
        frame_y = HMOViewerParameters.FRAME_HOMO_Y if occupied else HMOViewerParameters.FRAME_LUMO_Y
    
        # Effacer & redessiner le fond
        self.canvas.create_rectangle(frame_x, frame_y, frame_x + HMOViewerParameters.FRAME_WIDTH, frame_y + HMOViewerParameters.FRAME_HEIGHT, fill='white', outline='black')
    
        # === 🔖 Titre avec étiquette de l'OM ===
        nrj = f"α - {np.abs(self.energies[idx]):.2f}β" if self.energies[idx] < 0 else f"α + {np.abs(self.energies[idx]):.2f}β"
        om_label = f"OM #{idx+1} | {nrj} | {self.occupations[idx]}e"
        self.canvas.create_text(frame_x + HMOViewerParameters.FRAME_WIDTH/2, frame_y + 15, text=om_label, font=('DejaVu Sans', 10, 'italic'))
    
        # === 📈 Coordonnées des atomes ===
        coords = {row['Atom']: (row['X (grid units)'], row['Y (grid units)']) for _, row in self.df_atoms.iterrows()}
    
        xs = [pos[0] for pos in coords.values()]
        ys = [pos[1] for pos in coords.values()]
        min_x, max_x = min(xs), max(xs)
        min_y, max_y = min(ys), max(ys)
        # print(f"[DEBUG] Coordonnées min/max: X ({min_x}, {max_x}), Y ({min_y}, {max_y})")
    
        # === 📏 Calcul du scale intelligent pour viser TARGET_BOND_PX ===
        margin = 40
        
        mol_width = max_x - min_x
        mol_height = max_y - min_y
        
        # Le scale qui donnerait des liaisons de TARGET_BOND_PX
        scale_bond_based = HMOViewerParameters.TARGET_BOND_PX / self.mean_bond_length
        
        # Est-ce que ça passe dans le frame ?
        required_width = mol_width * scale_bond_based + margin
        required_height = mol_height * scale_bond_based + margin
        
        if required_width > HMOViewerParameters.FRAME_WIDTH or required_height > HMOViewerParameters.FRAME_HEIGHT:
            # Downscale forcé
            scale_x = (HMOViewerParameters.FRAME_WIDTH - margin) / mol_width if mol_width != 0 else 1
            scale_y = (HMOViewerParameters.FRAME_HEIGHT - margin) / mol_height if mol_height != 0 else 1
            scale = min(scale_x, scale_y)
            # print(f"[DEBUG] Downscaling: target scale {scale_bond_based:.2f} too big, reduced to {scale:.2f}")
        else:
            scale = scale_bond_based
            # print(f"[DEBUG] Using target scale: {scale:.2f}")
        scale *= scale_user  # <<<< ajoute ce facteur d'utilisateur

        # Offset pour bien centrer la molécule
        offset_x = frame_x + HMOViewerParameters.FRAME_WIDTH/2 - ((min_x + max_x) / 2) * scale
        offset_y = frame_y + HMOViewerParameters.FRAME_HEIGHT/2 - ((min_y + max_y) / 2) * scale
        
        # Taille max des lobes adaptée à la taille réelle après scaling
        # Trouver le coefficient max global pour normaliser
        max_coef_global = np.nanmax(np.abs(self.df_MOs.values))
        
        # Fixer la taille max du lobe : quand coef == max_coef_global, le rayon est 90% de longueur liaison/2
        desired_radius = 0.90 * self.mean_bond_length * scale / 2
        scale_lobe_factor = desired_radius / max_coef_global
        scale_lobe_factor *= lobe_user  # <<<< ajoute ce facteur d'utilisateur
        # print(f"[DEBUG] max_coef_global = {max_coef_global:.3f}, desired_radius = {desired_radius:.1f}, scale_lobe_factor = {scale_lobe_factor:.1f}")

        # === Coeffs
        coeffs = self.df_MOs.iloc[:, idx]

        # === ➿ Tracer les liaisons ===
        for _, row in self.df_bonds.iterrows():
            a1, a2 = row['Atom 1'], row['Atom 2']
            # print(f"Bond: {a1} <-> {a2}")
            x1, y1 = coords[a1]
            x2, y2 = coords[a2]
            self.canvas.create_line(
                offset_x + x1 * scale, offset_y + y1 * scale,
                offset_x + x2 * scale, offset_y + y2 * scale,
                width=2
            )
    
        # print(f"\nAffichage OM {idx} | Occupied={occupied} | Coeffs = {list(coeffs.items())}")
        if getattr(self, 'show_atom_labels', True):
            # === 🏷️ Tous les atomes : dessiner un cercle + label
            for atom, (x, y) in coords.items():
                # self.canvas.create_oval(
                #     offset_x + x * scale - 10, offset_y + y * scale - 10,
                #     offset_x + x * scale + 10, offset_y + y * scale + 10,
                #     fill='gray'
                # )
                self.canvas.create_text(offset_x + x * scale, offset_y + y * scale, text=atom)

        dx = HMOViewerParameters.LOBE_OFFSET * scale  # petit décalage droite
        dy = -HMOViewerParameters.LOBE_OFFSET * scale / 15  # petit décalage haut
        # === 🔵 Premier passage : les lobes "arrières" ===
        for atom, (x, y) in coords.items():
            coef = coeffs.get(atom, 0)
            if pd.isna(coef) or coef == 0:
                continue
            # print(f"Atome {atom} | Pos=({x},{y}) | Coef={coef}")
            size = abs(coef) * scale_lobe_factor
            color_front = 'red' if coef > 0 else 'blue'
            color_back = 'blue' if coef > 0 else 'red'
    
            # Lobe arrière (sous la liaison)
            self.canvas.create_oval(
                offset_x + x * scale + dx - size, offset_y + y * scale + dy - size,
                offset_x + x * scale + dx + size, offset_y + y * scale + dy + size,
                fill=color_back, stipple='gray50', outline=''
            )
    
        # === 🔴 Deuxième passage : les lobes "devant" ===
        for atom, (x, y) in coords.items():
            coef = coeffs.get(atom, 0)
            if pd.isna(coef) or coef == 0:
                continue
            size = abs(coef) * scale_lobe_factor
            dy = HMOViewerParameters.LOBE_OFFSET * scale / 10
            color_front = 'red' if coef > 0 else 'blue'
    
            # Lobe devant (au-dessus de la liaison)
            self.canvas.create_oval(
                offset_x + x * scale - size, offset_y + y * scale + dy - size,
                offset_x + x * scale + size, offset_y + y * scale + dy + size,
                fill=color_front, stipple='gray50', outline=''
            )

    
    def display_default_homo_lumo(self):
        """
        Automatically displays the HOMO and LUMO orbitals by detecting their indices.
    
        This method:
        - Identifies the highest occupied molecular orbital (HOMO) and the lowest unoccupied molecular orbital (LUMO).
        - Displays the HOMO in the occupied frame and the LUMO in the virtual frame.
    
        Notes
        -----
        If no HOMO or LUMO is found, nothing is displayed for that type.
        """        
        # print("Looking for HOMO/LUMO...")
    
        homo_idx = None
        lumo_idx = None
    
        homo_idx = None
        for idx in range(len(self.occupations)-1, -1, -1):
            if self.occupations[idx] > 0:
                homo_idx = idx
                break
        lumo_idx = None
        for idx in range(len(self.occupations)):
            if self.occupations[idx] == 0:
                lumo_idx = idx
                break

        # print(f"HOMO: idx {homo_idx+1}, energy {self.energies[homo_idx]}")
        # print(f"LUMO: idx {lumo_idx+1}, energy {self.energies[lumo_idx]}")
    
        if homo_idx is not None:
            self.current_occ_idx = homo_idx
            self.display_om(homo_idx, occupied=True)
        if lumo_idx is not None:
            self.current_virt_idx = lumo_idx
            self.display_om(lumo_idx, occupied=False)


    def on_level_click(self, idx):
        """
        Callback when a molecular orbital energy level is clicked in the diagram.
    
        Parameters
        ----------
        idx : int
            The index of the clicked MO.
    
        Notes
        -----
        - Determines if the MO is occupied or virtual based on its occupancy.
        - Updates the corresponding frame to display the clicked orbital using `display_om`.
        """
        occ = self.occupations[idx]
        occupied = occ > 0
        if occupied:
            self.current_occ_idx = idx
        else:
            self.current_virt_idx = idx

        self.display_om(idx, occupied)

# =============================================================================================================================================

class ChargeNode:
    """
    Represents a charge element on the canvas that cannot form bonds.
    It has a visual indicator (colored disk with sign), and affects the total number of π electrons.
    """
    charge_radius = 15  # display and interaction radius in pixels

    def __init__(self, x, y, charge='-1'):
        """
        Initialize a charge at (x, y) with the given charge label.
        Default is a single negative charge.
        """
        self.x = x
        self.y = y
        self.charge = charge
        self.selected = False

    def draw(self, canvas, x=None, y=None):
        """
        Draws the charge on the canvas. Uses scaled x/y if provided.
        """
        def ajouter_plus_si_necessaire(texte):
            if texte.startswith('-') or texte.startswith('+'):
                return texte
            try:
                # Vérifie que c'est bien un entier
                int(texte)
                return f'+{texte}'
            except ValueError:
                return texte  # Ne change rien si ce n'est pas un entier

        r = self.charge_radius
        draw_x = x if x is not None else self.x
        draw_y = y if y is not None else self.y
    
        if self.selected:
            canvas.create_oval(draw_x - r - 4, draw_y - r - 4, draw_x + r + 4, draw_y + r + 4,
                               outline='yellow', width=3)
        color = 'red' if '-' in self.charge else 'blue'
        text_color = 'white'
        # text_color = 'white' if '-' in self.charge else 'black'
        canvas.create_oval(draw_x - r, draw_y - r, draw_x + r, draw_y + r, fill=color)
        canvas.create_text(draw_x, draw_y, text=ajouter_plus_si_necessaire(self.charge), fill=text_color, font=('Helvetica', 12, 'bold'))

    def contains(self, x, y):
        """
        Checks if a point (x, y) lies within the charge disk (used for click detection).
        """
        r = self.charge_radius
        return (self.x - x)**2 + (self.y - y)**2 <= r**2

# =============================================================================================================================================

class Node:
    """
    A class representing an atom (node) in the molecular graph.

    Each Node object stores information about its position in the molecular diagram, 
    its chemical element type (e.g., C, N, O), and its π-charge (from Hückel analysis).

    Parameters
    ----------
    canvas : tk.Canvas
        The canvas on which the node (atom) will be drawn.
    x : int
        The x-coordinate of the node in the canvas (in pixels).
    y : int
        The y-coordinate of the node in the canvas (in pixels).
    radius : int, optional
        Radius of the circle representing the atom (default is 20 pixels).
    atom_type : str, optional
        Chemical element symbol for the node (default is 'C').
    color : str, optional
        Fill color for the node (default is 'white').

    Attributes
    ----------
    x : int
        X-coordinate of the node center.
    y : int
        Y-coordinate of the node center.
    radius : int
        Radius of the circle representing the atom.
    atom_type : str
        Chemical element symbol (e.g., 'C', 'N', 'O').
    color : str
        Fill color of the atom.
    circle : int
        ID of the oval shape on the canvas (for drawing reference).
    label : int
        ID of the text label on the canvas (for element symbol).
    pi_charge : float or None
        The π-electron charge on the atom, computed after Hückel analysis.

    Methods
    -------
    draw()
        Draws the node (circle + label) on the canvas.
    update_label()
        Updates the displayed label based on the atom type and charge.
    is_within(x, y)
        Checks if a given (x, y) coordinate is within the node's radius.
    move_to(new_x, new_y)
        Moves the node to a new position on the canvas.
    set_pi_charge(charge)
        Sets the π-charge for the atom and updates the label accordingly.
    """

    def __init__(self, x, y, atom_type='C·'):
        self.x = x
        self.y = y
        self.atom_type = atom_type

# =============================================================================================================================================

class MoleculeDrawer:
    """
    GUI application for drawing molecules and performing Hückel Molecular Orbital (HMO) calculations and analysis.

    The MoleculeDrawer class provides an interactive graphical interface that allows the user 
    to build a molecular structure by placing atoms and bonds on a grid. It supports the 
    customization of atom types, and includes functionalities to save/load molecules, undo/redo 
    actions, and erase atoms or bonds.

    Once a molecule is built, the tool can run a Hückel analysis to compute molecular orbital 
    energies, coefficients, π-charges, and π-bond orders. Results are displayed numerically and 
    graphically, with options to save data as Excel files or images.

    Parameters
    ----------
    master : tk.Tk or tk.Toplevel
        The root window or parent window that holds the GUI.

    Attributes
    ----------
    nodes : list of Node
        The list of atoms (nodes) in the current molecule.
    bonds : list of tuple
        The list of bonds, each defined by a tuple of two node indices.
    undo_stack : list
        A stack to store the history of molecule states for undo functionality.
    redo_stack : list
        A stack to store undone actions for redo functionality.
    df : pd.DataFrame or None
        The main DataFrame containing molecular orbital coefficients.
    df_atoms : pd.DataFrame or None
        DataFrame containing atom types, coordinates, and π charges.
    df_bonds : pd.DataFrame or None
        DataFrame listing bonds between atoms.
    df_summary : pd.DataFrame or None
        DataFrame summarizing molecular descriptors (energy, hardness, gap, etc.).
    project_name : str or None
        The user-defined name of the current project.
    safe_project_name : str or None
        A sanitized version of the project name (safe for filenames).
    mo_viewer_window : tk.Toplevel or None
        The window displaying molecular orbital visualizations.
    om_window : tk.Toplevel or None
        The window showing numerical molecular orbital data.

    Notes
    -----
    - The interface includes a toolbar with buttons for running the Hückel analysis, saving/loading 
      molecules, exporting results, and accessing help/about information.
    - Huckel parameters (α and β) are set by default to -11.0 and -2.7 but can be modified in the code.
    - The program supports exporting results as `.xlsx` files (multi-sheet Excel) and visual 
      representations as `.png` or `.pdf`.
    - The tool can gracefully handle molecule resizing and grid alignment, making drawing easy.

    Examples
    --------
    >>> root = tk.Tk()
    >>> app = MoleculeDrawer(root)
    >>> root.mainloop()
    """

    def __init__(self, master):
        self.master = master
        self.master.title("HMO Molecule Drawer")
        self.master.minsize(1200, 800)

        self.frame_main = tk.Frame(master)
        self.frame_main.pack(fill=tk.BOTH, expand=True)

        self.canvas = tk.Canvas(self.frame_main, width=1200, height=800, bg='white')
        self.canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        self.toolbar = tk.Frame(self.frame_main, width=120)
        self.toolbar.pack(side=tk.RIGHT, fill=tk.Y)

        self.nodes, self.bonds = [], []
        self.selected_node = None
        self.dragging = False
        self.eraser_mode = False
        self.current_mouse_pos = (0, 0)
        self.scale_x = 1.0
        self.scale_y = 1.0
        self.undo_stack, self.redo_stack = [], []

        # === drag an existing atom on a new position
        # Parameters for long-click selection (in milliseconds)
        self.long_click_duration = 1000  # 1 second threshold
        self.click_timer = None
        self.long_click_node = None
        # Binding mouse events for long-click detection
        self.canvas.bind("<ButtonPress-1>", self.on_button_press)
        self.canvas.bind("<ButtonRelease-1>", self.on_button_release_combined)


        self.load_icons()
        self.create_toolbar()
        self.bind_shortcuts()

        self.canvas.bind("<Button-3>", self.right_click)
        self.canvas.bind("<Motion>", self.mouse_motion)
        self.canvas.bind("<B1-Motion>", self.mouse_drag)
        self.canvas.bind("<Configure>", self.resize_canvas)

        self.om_window = None  # stocke la fenêtre DataFrame

        self.default_project_name = 'my_molecule'
        self.project_name = None
        self.safe_project_name = None
        self.molecule_is_new = True  # Nouvelle molécule par défaut
        
        self.df = None
        self.summary_data = None

        self.formal_charges = []  # list of ChargeNode instances
        self.charge_mol = 0
        self.mode = None

        self.alpha_value_num = -11.0 # eV
        self.beta_value_num = -2.7 # eV
        self.alpha_value = 0 # reduced unit
        self.beta_value = -1 # reduced unit



    def sanitize_filename(self, name):
        return re.sub(r'[\\/*?:"<>|]', "_", name)
        
    def create_toolbar(self):
        self.btn_run = self.create_button(self.icons['run'], self.on_run_huckel, "Run")
        self.create_button(self.icons['charge'], self.set_mode_add_charge, "Add a charge")
        self.create_button(self.icons['matrix'], self.show_numerical_data, "Show numerical results")
        self.create_button(self.icons['save'], self.save_molecule, "Save sigma skeleton")
        self.create_button(self.icons['load'], self.load_molecule, "Load sigma skeleton")
        self.create_button(self.icons['undo'], self.undo, "Undo last action")
        self.create_button(self.icons['redo'], self.redo, "Redo last undone action")
        self.btn_erase = self.create_button(self.icons['eraser'], self.toggle_eraser, "Erase atom or bond")
        self.create_button(self.icons['clear'], self.clear, "Delete molecule")
        tk.Frame(self.toolbar).pack(expand=True, fill=tk.BOTH)
        self.create_button(self.icons['savedata'], self.save_data, "Save data in a spreadsheet")
        self.create_button(self.icons['quit'], self.quit_program, "Quit")
        self.create_button(self.icons['about'], self.show_about, "About HMO")
        
    def load_icons(self):
        self.icons = {}
        for name in ['run', 'matrix', 'charge', 'save', 'load', 'undo', 'redo', 'eraser', 'clear', 'savedata', 'quit', 'about']:
            icon = resource_path(f"icons-logos-banner/{name}.png")
            self.icons[name] = ImageTk.PhotoImage(Image.open(icon))

    def create_button(self, icon, command, tooltip):
        btn = tk.Button(self.toolbar, image=icon, command=command, bg='lightgray')
        btn.pack(pady=2)
        btn.bind("<Enter>", lambda e: self.master.title(tooltip))
        btn.bind("<Leave>", lambda e: self.master.title("HMO Molecule Drawer"))
        ToolTip(btn, tooltip)
        return btn

    def bind_shortcuts(self):
        self.master.bind('<Control-r>', lambda e: self.on_run_huckel())
        self.master.bind('<Control-s>', lambda e: self.save_molecule())
        self.master.bind('<Control-o>', lambda e: self.load_molecule())
        self.master.bind('<Control-z>', lambda e: self.undo())
        self.master.bind('<Control-y>', lambda e: self.redo())
        self.master.bind('<Control-d>', lambda e: self.clear())
        self.master.bind('<Control-s>', lambda e: self.save_data())
        self.master.bind('<Control-q>', lambda e: self.quit_program())
        
    def set_mode_add_charge(self):
        """Sets the current interaction mode to adding a new charge."""
        self.mode = 'add_charge'
        
    def toggle_eraser(self):
        self.eraser_mode = not self.eraser_mode
        self.btn_erase.config(bg='#ffb8b9' if self.eraser_mode else 'lightgray')

    def quit_program(self):
        self.master.quit()
        
    def show_about(self):
        """
        Display the 'About' window for the HMO application.
    
        This window shows:
        - A banner/logo loaded from 'icons-logos-banner/HMO_Banner.png',
        - The author's address and credits,
        - The current version of the application,
        - A 'Close' button to dismiss the window.
    
        The window is non-resizable and can also be closed using the Escape key.
        """    
        # Create a new Toplevel window
        about_win = tk.Toplevel()
        about_win.title("About HMO")
        about_win.resizable(False, False)  # 🚫 Prevent resizing
    
        # === Load and resize the image ===
        img_path = resource_path(os.path.join('icons-logos-banner', 'HMO_Banner.png'))
        try:
            img = Image.open(img_path)
    
            # Resize: 50% of original size, with max width 500 px
            width, height = img.size
            new_width = min(width // 2, 500)
            new_height = int((new_width / width) * height)
            img_resized = img.resize((new_width, new_height), Image.LANCZOS)
    
            img_tk = ImageTk.PhotoImage(img_resized)
    
            img_label = tk.Label(about_win, image=img_tk)
            img_label.image = img_tk  # keep a reference
            img_label.pack(pady=(10, 5))
        except Exception as e:
            tk.Label(about_win, text=f"[Image not found: {e}]").pack(pady=(10, 5))
    
        # === Author text ===
        address_text = "Toulouse, France"
        address_label = tk.Label(about_win, text=address_text, font=("DejaVu Sans", 10))
        address_label.pack(pady=(5, 10))

        # === Author text ===
        author_text = "Author: Romuald Poteau, with the valuable help of ChatGPT (2025)"
        author_label = tk.Label(about_win, text=author_text, font=("DejaVu Sans", 8))
        author_label.pack(pady=(5, 5))
    
        # === Version ===
        version_label = tk.Label(about_win, text="Version 0.2.0", font=("DejaVu Sans", 10, "bold"))
        version_label.pack(pady=(5, 10))
        
        # Add Escape key binding to close the window
        about_win.bind('<Escape>', lambda event: about_win.destroy())
        
        # Close button
        close_btn = tk.Button(about_win, text="Close", command=about_win.destroy, font=("DejaVu Sans", 10))
        close_btn.pack(pady=(0, 10))
    
        # Make the window modal (optional)
        about_win.grab_set()

    def resize_canvas(self, event):
        if event.width > 0 and event.height > 0:
            self.scale_x = event.width / 800
            self.scale_y = event.height / 600
            self.redraw()

    def apply_scale(self, x, y):
        return x * self.scale_x, y * self.scale_y

    def draw_grid(self):
        self.canvas.delete("grid")
        w, h = self.canvas.winfo_width(), self.canvas.winfo_height()
        step_x = int(HuckelParameters.GRID_SIZE * self.scale_x)
        step_y = int(HuckelParameters.GRID_SIZE * self.scale_y)
        for x in range(0, w, step_x):
            for y in range(0, h, step_y):
                self.canvas.create_oval(x-1, y-1, x+1, y+1, fill='gray', outline='gray', tags="grid")

    def snap_to_grid(self, x, y):
        return (round(x / (HuckelParameters.GRID_SIZE * self.scale_x)) * HuckelParameters.GRID_SIZE * self.scale_x,
                round(y / (HuckelParameters.GRID_SIZE * self.scale_y)) * HuckelParameters.GRID_SIZE * self.scale_y)

    def find_node(self, x, y):
        for idx, node in enumerate(self.nodes):
            node_x, node_y = self.apply_scale(node.x, node.y)
            if abs(node_x - x) < HuckelParameters.ATOM_RADIUS and abs(node_y - y) < HuckelParameters.ATOM_RADIUS:
                return idx
        return None
        
    def node_exists_at(self, x, y):
        for node in self.nodes:
            if abs(node.x - x) < 1e-3 and abs(node.y - y) < 1e-3:
                return True
        return False

    def find_bond(self, x, y):
        for idx, (i, j) in enumerate(self.bonds):
            n1, n2 = self.nodes[i], self.nodes[j]
            mid_x, mid_y = self.apply_scale((n1.x + n2.x)/2, (n1.y + n2.y)/2)
            if abs(mid_x - x) < HuckelParameters.ATOM_RADIUS and abs(mid_y - y) < HuckelParameters.ATOM_RADIUS:
                return idx
        return None
        
    def on_button_release_combined(self, event):
        """
        Handles all left-button release logic, including both short and long click behaviors.
    
        This method combines and replaces the previous behaviors of `left_click` and `mouse_release`.
        It supports:
          - Long click on an atom: moves the atom to a new empty grid position (highlighted in green).
          - Short click + drag: creates a bond between the selected atom and another existing or new atom.
          - Short click on empty space: adds a new atom.
        
        It also works in conjunction with:
          - Eraser mode: activated during `on_button_press`.
          - Charge placement mode: managed during `on_button_press`.
    
        Parameters
        ----------
        event : tk.Event
            Mouse event triggered when the left mouse button is released.
        """
        if self.click_timer:
            self.master.after_cancel(self.click_timer)
            self.click_timer = None
    
        if self.long_click_node is not None:
            # Gestion clic long : déplacement d'atome existant vers un nœud vide
            x_snap, y_snap = self.snap_to_grid(event.x, event.y)
            x_grid = x_snap / self.scale_x
            y_grid = y_snap / self.scale_y
    
            if not self.node_exists_at(x_grid, y_grid) and not self.charge_exists_at(x_grid, y_grid):
                self.move_node_to_empty(self.long_click_node, x_grid, y_grid)
            else:
                messagebox.showwarning(
                    "Occupied position",
                    "The target position is already occupied. Please select an empty node."
                )
    
            self.long_click_node = None
            self.redraw()
    
        elif self.dragging and self.selected_node is not None:
            # Gestion classique : création liaison vers un nouveau nœud ou existant
            x_snap, y_snap = self.snap_to_grid(event.x, event.y)
            x_grid = x_snap / self.scale_x
            y_grid = y_snap / self.scale_y
    
            idx_target = self.find_node_by_coords(x_grid, y_grid)
            if idx_target is None:
                idx_target = self.find_node(event.x, event.y)
    
            if idx_target is None:
                self.save_state()
                self.nodes.append(Node(x_grid, y_grid))
                idx_target = len(self.nodes) - 1
    
            if (self.selected_node, idx_target) not in self.bonds and \
               (idx_target, self.selected_node) not in self.bonds and \
               self.selected_node != idx_target:
                self.save_state()
                self.bonds.append((self.selected_node, idx_target))
    
            self.selected_node = None
            self.dragging = False
            self.redraw()
    
        else:
            # Simple clic court sur un emplacement vide : création d'atome
            x_snap, y_snap = self.snap_to_grid(event.x, event.y)
            x_grid = x_snap / self.scale_x
            y_grid = y_snap / self.scale_y
            if not self.node_exists_at(x_grid, y_grid):
                self.save_state()
                self.nodes.append(Node(x_grid, y_grid))
                self.redraw()
    
        self.molecule_is_new = True


    
    def on_button_press(self, event):
        """
        Detects mouse button press and starts a timer if the user clicked on an existing atom.
    
        Parameters:
        -----------
        event : tk.Event
            Mouse event triggered by pressing button-1.
        """
        x_canvas, y_canvas = event.x, event.y
        x, y = self.snap_to_grid(x_canvas, y_canvas)
        x /= self.scale_x
        y /= self.scale_y
        
        # 1️⃣ Gomme : suppression de charges, d'atomes ou de liaisons
        if self.eraser_mode:
            for i, c in enumerate(self.formal_charges):
                if c.contains(x, y):
                    self.save_state()
                    del self.formal_charges[i]
                    self.redraw()
                    return
            idx = self.find_node(event.x, event.y)
            if idx is not None:
                self.save_state()
                self.delete_node(idx)
                self.redraw()
                return
            idx = self.find_bond(event.x, event.y)
            if idx is not None:
                self.save_state()
                del self.bonds[idx]
                self.redraw()
                return
            return
        
        # 2️⃣ Mode ajout de charge
        if self.mode == 'add_charge':
            if self.node_exists_at(x, y):
                messagebox.showwarning("Invalid location", "You cannot place a charge on top of an atom.")
                return
            for other in self.formal_charges:
                other.selected = False
            self.formal_charges.append(ChargeNode(x, y))
            self.mode = None  # exit charge mode after placing one
            self.charge_mol = -1
            self.redraw()
            return

        idx = self.find_node(event.x, event.y)
        if idx is not None:
            # Start timer to detect long-click on existing atom
            self.click_timer = self.master.after(
                self.long_click_duration,
                lambda: self.start_long_click(idx)
            )
            self.selected_node = idx
            self.dragging = True

    def start_long_click(self, idx):
        """
        Activates the long-click state by highlighting the selected atom.
    
        Parameters:
        -----------
        idx : int
            Index of the atom in self.nodes that was clicked.
        """
        self.long_click_node = idx
        self.highlight_node(idx)

    def highlight_node(self, idx):
        """
        Visually highlights the atom to indicate selection for moving.
    
        Parameters:
        -----------
        idx : int
            Index of the atom to highlight.
        """
        node = self.nodes[idx]
        x, y = self.apply_scale(node.x, node.y)
        r = HuckelParameters.ATOM_RADIUS + 6
        self.canvas.create_oval(
            x - r, y - r, x + r, y + r,
            outline='green', width=3, tag='highlight'
        )

    def move_node_to_empty(self, node_idx, new_x, new_y):
        """
        Moves the selected atom to a new empty location, preserving its existing bonds.
    
        Parameters:
        -----------
        node_idx : int
            Index of the atom to move.
        new_x : float
            X-coordinate (grid units) of the new empty location.
        new_y : float
            Y-coordinate (grid units) of the new empty location.
        """
        node = self.nodes[node_idx]
        node.x = new_x
        node.y = new_y

    def charge_exists_at(self, x, y):
        """
        Checks if a charge node already exists at the given grid position.
    
        Parameters:
        -----------
        x : float
            X-coordinate (grid units) to check.
        y : float
            Y-coordinate (grid units) to check.
    
        Returns:
        --------
        bool
            True if a charge exists at the location, False otherwise.
        """
        for charge in self.formal_charges:
            if abs(charge.x - x) < 1e-3 and abs(charge.y - y) < 1e-3:
                return True
        return False

    
    def save_state(self):
        self.undo_stack.append((self.nodes.copy(), self.bonds.copy()))
        self.redo_stack.clear()

    def undo(self):
        if self.undo_stack:
            self.redo_stack.append((self.nodes.copy(), self.bonds.copy()))
            self.nodes, self.bonds = self.undo_stack.pop()
            self.redraw()

    def redo(self):
        if self.redo_stack:
            self.undo_stack.append((self.nodes.copy(), self.bonds.copy()))
            self.nodes, self.bonds = self.redo_stack.pop()
            self.redraw()


    def right_click(self, event):
        x_canvas, y_canvas = event.x, event.y
        x, y = self.snap_to_grid(x_canvas, y_canvas)
        x /= self.scale_x
        y /= self.scale_y

        # Vérifie si une charge a été cliquée
        for c in self.formal_charges:
            if c.contains(x_canvas / self.scale_x, y_canvas / self.scale_y):
                for other in self.formal_charges:
                    other.selected = False
                c.selected = True
                new_charge = simpledialog.askstring(
                    "Edit charge", 
                    "Enter charge (e.g., -1, +2, -1, -2):", 
                    initialvalue=c.charge
                )
                if new_charge:
                    # Suppression des espaces
                    cleaned = new_charge.replace(" ", "")
                    # Vérification avec une expression régulière stricte
                    if re.fullmatch(r"[+-]?\d+", cleaned):
                        try:
                            val = int(cleaned)
                            c.charge = cleaned
                            self.charge_mol = val
                        except ValueError:
                            messagebox.showerror("Invalid charge", f"Could not convert '{cleaned}' to integer.")
                    else:
                        messagebox.showerror("Invalid charge", f"Invalid charge format: '{new_charge}'")
                        
                self.redraw()
                return  # ne continue pas vers le menu des atomes
        
        idx = self.find_node(x_canvas, y_canvas)
        if idx is not None:
            menu = tk.Menu(self.master, tearoff=0)
            menu.configure(font=font.Font(family="Arial", size=12))
            bond_count = self.count_bonds(idx)
            for atom_type in HuckelParameters.ATOM_OPTIONS:
                if atom_type == 'Me' and bond_count > 1:
                    menu.add_command(label=atom_type + " (blocked)", state='disabled')
                else:
                    menu.add_command(label=atom_type,
                        command=lambda at=atom_type: self.change_atom_type(idx, at))
            try:
                menu.tk_popup(event.x_root, event.y_root)
            finally:
                menu.grab_release()
            self.molecule_is_new = True

    def change_atom_type(self, idx, new_type):
        if new_type in HuckelParameters.ATOM_COLORS:
            self.save_state()
            self.nodes[idx].atom_type = new_type
            self.redraw()

    def mouse_drag(self, event):
        if self.dragging:
            self.current_mouse_pos = (event.x, event.y)
            self.redraw()

    def mouse_motion(self, event):
        self.current_mouse_pos = (event.x, event.y)
        self.redraw()

    def find_node_by_coords(self, x_grid, y_grid):
        for idx, node in enumerate(self.nodes):
            if abs(node.x - x_grid) < 1e-5 and abs(node.y - y_grid) < 1e-5:
                return idx
        return None

    
    def count_bonds(self, idx):
        return sum(1 for i, j in self.bonds if i == idx or j == idx)

    def delete_node(self, idx):
        del self.nodes[idx]
        self.bonds = [(i, j) for i, j in self.bonds if i != idx and j != idx]
        self.bonds = [(i - (i > idx), j - (j > idx)) for i, j in self.bonds]
        self.molecule_is_new = True
        if not self.nodes:
            self.project_name = self.default_project_name
            print(f"[INFO] Molecule is now empty; project name reset to: {self.project_name}")

    def clear(self):
        self.nodes.clear()
        self.bonds.clear()
        self.undo_stack.clear()
        self.redo_stack.clear()
        self.molecule_is_new = True
        self.project_name = self.default_project_name
        print(f"[INFO] Project name reset to default: {self.project_name}")
        
        self.df = None  # Optionnel si tu veux réinitialiser la dernière analyse Hückel
        self.redraw()

        if self.mo_viewer_window is not None and self.mo_viewer_window.winfo_exists():
            self.mo_viewer_window.destroy()
            self.mo_viewer_window = None

        if self.om_window is not None and self.om_window.winfo_exists():
            self.om_window.destroy()
            self.om_window = None

        
    def save_molecule(self):
        """
        Save the current molecular structure to a .hmo file.
    
        This method opens a file dialog prompting the user to select the save location.
        The molecule is saved in a simple text format:
        - A "Nodes:" section listing each atom with its type and (x, y) grid coordinates.
        - A "Bonds:" section listing each bond as pairs of atom indices.
    
        If a project name is available (`self.safe_project_name`), it is used as the default filename.
    
        Notes
        -----
        The saved file allows reloading the same molecular structure later using `load_molecule`.
        """
        if self.safe_project_name is None:
            path = filedialog.asksaveasfilename(
                defaultextension=".hmo",
                filetypes=[("Hückel molecule files", "*.hmo"), ("All Files", "*.*")]
            )
        else:
            path = filedialog.asksaveasfilename(
                defaultextension=".hmo",
                filetypes=[("Hückel molecule files", "*.hmo"), ("All Files", "*.*")],
                initialfile=f"{self.safe_project_name}.hmo"
            )

        if path:
            with open(path, 'w') as f:
                f.write("Nodes:\n")
                for node in self.nodes:
                    f.write(f"{node.atom_type} {node.x} {node.y}\n")
                f.write("Bonds:\n")
                for i, j in self.bonds:
                    f.write(f"{i} {j}\n")

    def load_molecule(self):
        """
        Load a molecular structure from a .hmo file.
    
        This method opens a file dialog prompting the user to select a .hmo file.
        It reads the file content and reconstructs the molecule by:
        - Parsing the "Nodes:" section to recreate the list of atoms (`self.nodes`).
        - Parsing the "Bonds:" section to recreate the list of bonds (`self.bonds`).
    
        The current state is first saved (for undo purposes), and the canvas is redrawn
        after loading. If an error occurs (e.g., corrupted file), an error message is shown.
    
        Notes
        -----
        The expected format of the file is:
        Nodes:
        AtomType x y
        ...
        Bonds:
        i j
        ...
        """    
        path = filedialog.askopenfilename(
            defaultextension=".hmo",
            filetypes=[("Hückel molecule files", "*.hmo"), ("All Files", "*.*")]
        )
        if path:
            try:
                with open(path) as f:
                    lines = f.readlines()
                self.save_state()
                self.nodes.clear()
                self.bonds.clear()
                mode = None
                for line in lines:
                    line = line.strip()
                    if line == "Nodes:":
                        mode = "nodes"
                    elif line == "Bonds:":
                        mode = "bonds"
                    elif mode == "nodes":
                        parts = line.split()
                        self.nodes.append(Node(float(parts[1]), float(parts[2]), parts[0]))
                    elif mode == "bonds":
                        i, j = map(int, line.split())
                        self.bonds.append((i, j))
                    # print("Line:", line, "| Mode:", mode)
                # print(self.nodes)
                # print(self.bonds)
                self.redraw()
                # Update project_name based on the file name (without extension)
                filename = os.path.basename(path)
                project_name_no_ext = os.path.splitext(filename)[0]
                self.project_name = project_name_no_ext
                print(f"[INFO] Project name set to: {self.project_name}")
                self.molecule_is_new = True
            except Exception as e:
                messagebox.showerror("Error", f"Failed to load molecule: {e}")

    def redraw(self):
        """
        Redraws the full molecular canvas, including grid, atoms, and bonds.
    
        This method is responsible for updating the Tkinter canvas whenever the molecule's structure changes or when the user interacts (e.g., dragging, erasing).
    
        It handles:
        
        1️⃣ Redrawing the base grid (for alignment).
        
        2️⃣ If dragging:
            - Displays a highlighted grid cell (orange) under the current mouse position.
            - Draws a temporary dashed bond between the selected node and the cursor.
            
        3️⃣ Highlights:
            - If in eraser mode, highlights the node or bond under the mouse cursor in red.
            
        4️⃣ Bonds:
            - Draws all bonds as lines between atoms, with special highlighting for bonds targeted by the eraser.
            
        5️⃣ Atoms:
            - Draws each atom as a colored circle, using the atom type's color code.
            - Highlights atoms under the eraser with a red border.
            - Displays the atom label (type) centered in each circle.
    
        Notes
        -----
        - This method fully clears the canvas before redrawing.
        - Visual cues (colors, highlights) guide user interaction for precision editing.
        - Grid snapping ensures precise atom placement during drawing.
        """

        self.canvas.delete("all")
        self.draw_grid()
        self.canvas.delete('highlight')  # remove existing highlight if any
    
        if self.dragging:
            snap_x, snap_y = self.snap_to_grid(*self.current_mouse_pos)
            self.canvas.create_oval(
                snap_x - HuckelParameters.HIGHLIGHT_RADIUS, snap_y - HuckelParameters.HIGHLIGHT_RADIUS,
                snap_x + HuckelParameters.HIGHLIGHT_RADIUS, snap_y + HuckelParameters.HIGHLIGHT_RADIUS,
                fill='orange', outline='orange'
            )

        highlight_node = self.find_node(*self.current_mouse_pos) if self.eraser_mode else None
        highlight_bond = self.find_bond(*self.current_mouse_pos) if self.eraser_mode else None

        for i, j in self.bonds:
            x1, y1 = self.apply_scale(self.nodes[i].x, self.nodes[i].y)
            x2, y2 = self.apply_scale(self.nodes[j].x, self.nodes[j].y)
            if highlight_bond is not None and (i, j) == self.bonds[highlight_bond]:
                self.canvas.create_line(x1, y1, x2, y2, fill='red', width=4)
            else:
                self.canvas.create_line(x1, y1, x2, y2, width=2)

        if self.dragging and self.selected_node is not None:
            x0, y0 = self.apply_scale(self.nodes[self.selected_node].x, self.nodes[self.selected_node].y)
            self.canvas.create_line(x0, y0, self.current_mouse_pos[0], self.current_mouse_pos[1], dash=(4,2))

        for idx, node in enumerate(self.nodes):
            x, y = self.apply_scale(node.x, node.y)
            color = HuckelParameters.ATOM_COLORS.get(node.atom_type, 'black')
            if highlight_node == idx:
                self.canvas.create_oval(
                    x - HuckelParameters.ATOM_RADIUS - 4, y - HuckelParameters.ATOM_RADIUS - 4,
                    x + HuckelParameters.ATOM_RADIUS + 4, y + HuckelParameters.ATOM_RADIUS + 4,
                    outline='red', width=2
                )
            self.canvas.create_oval(
                x - HuckelParameters.ATOM_RADIUS, y - HuckelParameters.ATOM_RADIUS,
                x + HuckelParameters.ATOM_RADIUS, y + HuckelParameters.ATOM_RADIUS,
                fill=color
            )
            self.canvas.create_text(x, y, text=node.atom_type, fill='white', font=('Arial', 8, 'bold'))
            
        for c in self.formal_charges:
            x, y = self.apply_scale(c.x, c.y)
            c.draw(self.canvas, x, y)

    
    def evaluate(self,expr, alpha, beta):
        return eval(expr, {"alpha": alpha, "beta": beta})
         
    def build_huckel_matrix(self, alpha=0.0, beta=-1.0):
        """
        Constructs the Hückel matrix (Hamiltonian) for the current molecular graph.
    
        This method builds the square matrix H of size n×n, where n is the number of atoms (nodes) 
        in the molecule. The diagonal elements (alpha) represent Coulomb integrals, while the 
        off-diagonal elements (beta) represent resonance integrals between bonded atoms.
    
        Atomic parameters are retrieved from the `Huckel_atomic_parameters` dictionary using the 
        atom's type (e.g., 'C', 'N', etc.), and bond parameters are retrieved from the 
        `Huckel_kXY_parameters` dictionary to handle heteroatomic systems with scaling factors.
    
        Parameters
        ----------
        alpha : float, optional
            The Coulomb integral (diagonal value), default is 0.
        beta : float, optional
            The resonance integral (off-diagonal value), default is -1.
            
        Note:
        -----
        The Hückel matrix H is transformed using the standard variable change:
            x = (alpha - epsilon) / beta
        where x = 0 is chosen as the energy origin.
        The entire matrix is divided by beta, so all energies (eigenvalues) are expressed in units of beta.
        Alpha is conceptually shifted to zero during this transformation.
    
        The total π-electron energy is separated into:
        - alpha_part = n_pi_total * alpha
        - beta_part  = sum of occupied eigenvalues (in beta units)    
        
        Returns
        -------
        np.ndarray
            The n×n Hückel matrix, where n is the number of atoms in the molecule.
    
        Error Handling
        --------------
        - If an atom type is not recognized or has no defined alpha parameter, an error message 
          is shown and the program exits.
        - If a bond type is not recognized or has no defined beta parameter, an error message 
          is shown and the program exits.
    
        Notes
        -----
        - The `Huckel_atomic_parameters` dictionary must provide a valid `alpha_expr` for each 
          atom type involved in the molecule.
        - The `Huckel_kXY_parameters` dictionary must provide scaling factors for bonds between 
          different atom types (e.g., C-N, C-O).
        - The method evaluates alpha expressions using the `self.evaluate()` method, which 
          substitutes the global alpha and beta values.
    
        Example
        -------
        >>> H = drawer.build_huckel_matrix(alpha=-11.0, beta=-2.7) or H = drawer.build_huckel_matrix(alpha=0, beta=-1) in reduced units
        """

        n = len(self.nodes)
        H = np.zeros((n, n))
        ### alpha (diagonal value)
        for i, node in enumerate(self.nodes):
            atom_type = node.atom_type
            try:
                param = HuckelParameters.Huckel_atomic_parameters[atom_type]
                if param is None:
                    raise ValueError(f"The parameter for atom '{atom_type}' is None.")
            except (KeyError, ValueError) as e:
                print(f"[ERROR] Unknown alpha parameter for atom '{atom_type}': {e}")
                messagebox.showerror(
                    "Error: Unknown atomic parameter",
                    f"Error: Unknown alpha parameter for atom '{atom_type}'.\n\n"
                    "Please check your molecule.\n\n"
                    "The program will now close."
                )
                self.master.quit()
                return
            H[i, i] = self.evaluate(param["alpha_expr"], self.alpha_value, self.beta_value)
        ### beta (non-diagonal value)
        for i, j in self.bonds:
            a = self.nodes[i].atom_type
            b = self.nodes[j].atom_type
        
            k_ij = HuckelParameters.Huckel_kXY_parameters.get(a, {}).get(b, None)
            if k_ij is None:
                k_ij = HuckelParameters.Huckel_kXY_parameters.get(b, {}).get(a, None)
            if k_ij is None:
                print(f"[ERROR] Unknown beta parameter for bond '{a} - {b}'.")
                messagebox.showerror(
                    "Error: Unknown bond parameter",
                    f"Error: Unknown beta parameter for bond '{a} - {b}'.\n\n"
                    "Please check your molecule.\n\n"
                    "The program will now close."
                )
                self.master.quit()
                return
            H[i, j] = H[j, i] = k_ij * self.beta_value
        return H

    def props(self, eigvals, occupation_dict):
        """
        Compute global energetic and chemical descriptors from Hückel eigenvalues.
    
        This method calculates the total π-electron energy of the molecule, decomposed into
        its alpha (on-site) and beta (interaction) parts, as well as the atomization energy
        (difference between the total π-energy and the isolated atom contributions).
    
        It also computes key global descriptors:
        - HOMO-LUMO gap
        - chemical potential (mu)
        - chemical hardness (eta)
        - softness (1/eta)
        - electrophilicity index (omega)
    
        Parameters
        ----------
        eigvals : array-like
            The eigenvalues (orbital energies) from the Hückel matrix diagonalization.
        occupation_dict : dict
            A dictionary mapping orbital indices to their occupation numbers (typically 2 for
            occupied, 0 for virtual).
    
        Notes
        -----
        - This method updates instance attributes:
          `total_energy`, `alpha_part`, `beta_part`, `atomization_energy`,
          `atomization_energy_per_atom`, `homo_lumo_gap`, `mu`, `eta`, `softness`, `omega`.
        - Includes debug print statements for tracing the calculation steps.
        """
        # Total π-electron energy = beta part, if reduced unit
        # print(eigvals)
        total_energy = sum(eigvals[j] * occ for j, occ in occupation_dict.items())
        beta_part = total_energy
        
        # for j, occ in occupation_dict.items():
        #     print(j,eigvals[j], occ)
    
        # Partie alpha réelle (somme sur chaque atome : n_pi * alpha_effectif)
        alpha_part = 0
        alpha_atoms = 0
        for i,n in enumerate(self.nodes):
            params = HuckelParameters.Huckel_atomic_parameters.get(n.atom_type)
            n_pi = params["n_pi"]
            try:
                alpha_atom = eval(params["alpha_expr"], {"alpha": self.alpha_value, "beta": self.beta_value})
            except Exception as e:
                print(f"Erreur pour {n.atom_type} : {e}")
                alpha_atom = self.alpha_value  # fallback sécurité
            alpha_atoms += n_pi * alpha_atom
            # print(i,alpha_atom,n_pi, alpha_atoms)
        alpha_part = self.total_pi_electrons
        
        # Determine if atomization energy is valid
        atom_types = set(n.atom_type for n in self.nodes)
        is_all_carbons = atom_types == {"C"}
        is_neutral = (self.charge_mol == 0)
        if is_neutral or is_all_carbons:
            atomization_energy = total_energy - alpha_atoms
            atomization_energy_per_atom = atomization_energy / len(self.nodes) if self.nodes else 0
        else:
            print("[INFO] Atomization energy skipped: molecule is charged and contains heteroatoms.")
            atomization_energy = None
            atomization_energy_per_atom = None
    
        # Atomization energy (comme avant)
        atomization_energy = total_energy - alpha_atoms
        atomization_energy_per_atom = atomization_energy / len(self.nodes) if self.nodes else 0
    
        # Trouver HOMO et LUMO
        occupied_indices = [j for j, occ in occupation_dict.items() if occ > 0]
        virtual_indices = [j for j, occ in occupation_dict.items() if occ == 0]
    
        E_HOMO = eigvals[max(occupied_indices)] if occupied_indices else None
        E_LUMO = eigvals[min(virtual_indices)] if virtual_indices else None
    
        # Descripteurs chimiques
        if E_HOMO is not None and E_LUMO is not None:
            mu = (E_HOMO + E_LUMO) / 2
            eta = (E_LUMO - E_HOMO) / 2
            softness = 1 / eta if eta != 0 else np.inf
            omega = 0.5 * mu ** 2 / eta if eta != 0 else np.inf
            homo_lumo_gap = E_LUMO - E_HOMO
        else:
            mu = eta = softness = omega = homo_lumo_gap = None
    
        # Stocker dans self
        self.total_energy = total_energy
        self.alpha_part = alpha_part
        self.beta_part = beta_part
        self.atomization_energy = atomization_energy
        self.atomization_energy_per_atom = atomization_energy_per_atom
        self.homo_lumo_gap = homo_lumo_gap
        self.mu = mu
        self.eta = eta
        self.softness = softness
        self.omega = omega
        
    def compute_charges_and_bond_orders(self, eigvecs, occupation_dict):
        """
        Compute atomic π-charges and π-bond orders from Hückel eigenvectors.
    
        This method evaluates:
        - The π-charge on each atom, as the difference between its formal π-electrons and the Mulliken-like population derived from the eigenvector coefficients.
        - The π-bond order between each pair of bonded atoms, as a weighted sum of the products of their orbital coefficients across all occupied molecular orbitals.
    
        Parameters
        ----------
        eigvecs : ndarray (n_atoms, n_orbitals)
            Matrix of eigenvectors (molecular orbital coefficients) from the Hückel diagonalization.
            Each column corresponds to an orbital; each row to an atom.
        occupation_dict : dict
            Dictionary mapping orbital indices to their occupation numbers (e.g., 2 or 0).
    
        Returns
        -------
        charges : list of float
            The computed π-charge for each atom, in the same order as `self.nodes`.
        bond_orders : dict
            A dictionary with keys as (i, j) tuples (atom indices of bonded pairs) and values as the computed π-bond order between those atoms.
    
        Notes
        -----
        - The method assumes that each atom's expected number of π-electrons (`n_pi`) is defined in `HuckelParameters.Huckel_atomic_parameters`.
        - The bond list `self.bonds` should contain tuples of bonded atom indices.
        """
        n_atoms = len(self.nodes)
        n_om = eigvecs.shape[1]
    
        # Charges : initialise à nombre π d'électrons
        charges = []
        for i, node in enumerate(self.nodes):
            n_pi = HuckelParameters.Huckel_atomic_parameters[node.atom_type]["n_pi"]
            total = 0
            for j in range(n_om):
                occ_j = occupation_dict.get(j, 0)
                coeff = eigvecs[i, j]
                total += occ_j * (coeff ** 2)
            q_i = n_pi - total
            charges.append(q_i)

        # Bond orders : dictionnaire {(i, j): valeur}
        bond_orders = {}
        for (i, j) in self.bonds:
            total = 0
            for k in range(n_om):
                occ_k = occupation_dict.get(k, 0)
                ci = eigvecs[i, k]
                cj = eigvecs[j, k]
                total += occ_k * ci * cj
            bond_orders[(i, j)] = total
    
        return charges, bond_orders


    def run_huckel_analysis(self):
        """
        Executes the full Hückel Molecular Orbital (HMO) analysis for the current molecule.
    
        This method performs the following steps:
        1. Builds the Hückel matrix using `build_huckel_matrix`.
        2. Solves for molecular orbital energies and coefficients (eigenvalues and eigenvectors).
        3. Counts the total number of π-electrons based on atomic types.
        4. Computes the occupation numbers for each molecular orbital, handling degeneracies.
        5. Calculates π-electron charges and π-bond orders using the eigenvectors and occupations.
        6. Updates molecular properties and builds a DataFrame of orbital coefficients for display.
    
        The analysis assumes default Hückel parameters (alpha = -11.0, beta = -2.7) unless specified 
        otherwise when building the matrix.
    
        Notes
        -----
        - Atoms and bonds are validated to ensure parameters exist in the Huckel parameter dictionaries.
        - If any atom or bond type is unknown, an error message is displayed and the program exits.
        - Degenerate orbitals are handled by grouping them within a numerical tolerance (1e-5 by default).
    
        Attributes Updated
        ------------------
        df : pd.DataFrame
            A DataFrame showing molecular orbital coefficients, sorted by energy levels.
        charges : list of float
            List of π-electron charges on each atom.
        bond_orders : dict
            Dictionary with keys as (i, j) tuples (atom indices) and values as bond orders.
        total_pi_electrons : int
            The total number of π-electrons in the molecule.
        alpha_value : float
            The alpha value used for the current Hückel analysis.
        beta_value : float
            The beta value used for the current Hückel analysis.
    
        Error Handling
        --------------
        - Unknown atom types or missing parameters in Huckel_atomic_parameters trigger an error.
        - Unknown bond types or missing parameters in Huckel_kXY_parameters also trigger an error.
        - In both cases, a messagebox error is shown, and the program exits gracefully.
    
        Example
        -------
        >>> MoleculeDrawer.run_huckel_analysis()
        """

        def compute_occupations(eigvals, tol=1e-5):
            """
            Computes the occupation numbers of molecular orbitals based on the total π-electron count.
        
            The function sorts eigenvalues (molecular orbital energies), identifies degenerate orbitals 
            within a specified tolerance, and distributes electrons accordingly:
            - Doubly occupies non-degenerate orbitals (max 2 electrons).
            - For degenerate groups: first assigns 1 electron to each, then adds a second electron if 
              electrons remain.
        
            Parameters
            ----------
            eigvals : array-like
                The list or array of eigenvalues (orbital energies) from Hückel analysis.
            tol : float, optional
                Tolerance used to group degenerate orbitals (default is 1e-5).
        
            Returns
            -------
            dict
                A dictionary mapping orbital indices to their occupation numbers (0, 1, or 2).
        
            Notes
            -----
            - The method handles degeneracies by grouping orbitals whose energies differ by less than `tol`.
            - Electrons are filled starting from the lowest energy orbitals upward.
            - The total number of electrons is taken from `self.total_pi_electrons` in the parent scope.
            - For degenerate orbitals, the algorithm first places 1 electron in each, then fills up to 2 
              electrons where possible, to respect electron pairing rules.
        
            Example
            -------
            >>> occupations = compute_occupations(eigvals, tol=1e-5)
            # Output: {0: 2, 1: 2, 2: 2, 3: 0, ...}
            """

            # eigvals_sorted = np.sort(eigvals)
            eigvals_sorted = eigvals.copy()
            n_electrons = self.total_pi_electrons
            
            # === Initialisation ===
            remaining_electrons = self.total_pi_electrons
            occupation_dict = {i: 0 for i in range(len(eigvals_sorted))}
            
            # Regrouper les OM dégénérées par INDICE
            tol = 1e-5
            groups = []
            current_group = [0]
            
            for j in range(1, len(eigvals_sorted)):
                if abs(eigvals_sorted[j] - eigvals_sorted[current_group[-1]]) < tol:
                    current_group.append(j)
                else:
                    groups.append(current_group)
                    current_group = [j]
            groups.append(current_group)  # ajouter le dernier groupe
            
            # === Traitement groupe par groupe ===
            for group in groups:
                group_size = len(group)
            
                if group_size == 1:
                    idx = group[0]
                    to_add = min(2, remaining_electrons)
                    occupation_dict[idx] += to_add
                    remaining_electrons -= to_add
            
                else:
                    # Premier passage : 1 électron dans chaque OM du groupe
                    for idx in group:
                        if remaining_electrons == 0:
                            break
                        occupation_dict[idx] += 1
                        remaining_electrons -= 1
            
                    # Deuxième passage : ajouter 2e électron si encore dispo
                    for idx in group:
                        if remaining_electrons == 0:
                            break
                        if occupation_dict[idx] < 2:
                            occupation_dict[idx] += 1
                            remaining_electrons -= 1
            
                # Debug utile :
                # print(f"Groupe indices: {group} ➔ Occ: {[occupation_dict[idx] for idx in group]}, Remaining: {remaining_electrons}")
            
                if remaining_electrons == 0:
                    break

            return occupation_dict

        H = self.build_huckel_matrix(self.alpha_value, self.beta_value)
        if H is None:
            print("[INFO] Hückel analysis stopped due to a matrix-building error.")
            return
        eigvals, eigvecs = np.linalg.eigh(H)
        eigvals = eigvals/ self.beta_value

        self.total_pi_electrons = 0
        for idx, n in enumerate(self.nodes):
            atom_type = n.atom_type
            try:
                param = HuckelParameters.Huckel_atomic_parameters[atom_type]
                n_pi = param["n_pi"]
                if n_pi is None:
                    raise ValueError(f"The 'n_pi' parameter is None for atom '{atom_type}'.")
                alpha_value = self.evaluate(param["alpha_expr"], self.alpha_value, self.beta_value)
                # print(f"[INFO] Atom #{idx+1} ('{atom_type}'): n_pi = {n_pi}, alpha = {alpha_value:.2f}")  # ✅ log for each atom
            except KeyError:
                print(f"[ERROR] Unknown atom type '{atom_type}' when counting π electrons.")
                messagebox.showerror(
                    "Error: Unknown atomic parameter",
                    f"Error: Atom type '{atom_type}' is unknown when counting π electrons.\n\n"
                    "Please check your molecule.\n\n"
                    "The program will now close."
                )
                self.master.quit()
                return
            except ValueError as e:
                print(f"[ERROR] {e}")
                messagebox.showerror(
                    "Error: Invalid atomic parameter",
                    f"Error: The 'n_pi' parameter is invalid for atom '{atom_type}'.\n\n"
                    "Please check your molecule.\n\n"
                    "The program will now close."
                )
                self.master.quit()
                return
            self.total_pi_electrons += n_pi
        # ➕ Correction par la charge globale
        # print(f"[INFO] Total π-electrons before charge adjustment: {self.total_pi_electrons}")
        # print(f"[INFO] Charge correction applied: {self.charge_mol}")
        self.total_pi_electrons -= self.charge_mol
        # print(f"[INFO] Final π-electron count: {self.total_pi_electrons}")
        
        occupation_dict = compute_occupations(eigvals, self.total_pi_electrons)
        # sorted_indices = np.argsort(eigvals)[::-1]
        sorted_indices = np.arange(0, len(occupation_dict))
        
        self.charges, self.bond_orders = self.compute_charges_and_bond_orders(eigvecs, occupation_dict)

        self.props(eigvals, occupation_dict)

        labels = {}
        counts = {}
        for i, n in enumerate(self.nodes):
            symbol = n.atom_type
            labels[i] = f"{symbol}{i+1}"
        data = []
        for i in range(len(self.nodes)):
            row = []
            for j in sorted_indices:
                row.append(round(eigvecs[i, j], 3))
            data.append(row)
        columns = []
        for idx, j in enumerate(sorted_indices):
            e_rounded = round(eigvals[j], 5)
            occ = occupation_dict[j]
            energy_coeff = eigvals[j]
            energy_coeff_rounded = round(energy_coeff, 2)
            
            # Évite les + -0.00 ou -0.00 avec formatage intelligent
            if abs(energy_coeff_rounded) < 1e-6:
                energy_label = "E = α"
            elif energy_coeff_rounded > 0:
                energy_label = f"E = α + {energy_coeff_rounded:.2f}β"
            else:
                energy_label = f"E = α - {abs(energy_coeff_rounded):.2f}β"
            occ_label = f"{occ}e"
            header = f"{energy_label}\n{occ_label}"
            columns.append(header)
        index_labels = [labels[i] for i in range(len(self.nodes))]
        df = pd.DataFrame(data, index=index_labels, columns=columns)
        # df = df.iloc[:, ::-1]
        self.df = df
        del df
        
    def build_dataframes(self):
        """
        Builds and updates key DataFrames summarizing the molecule's structural and electronic properties.
    
        This method generates and assigns four main DataFrames:
        1. `df_bond_orders_matrix`: A symmetric matrix showing π-bond orders between all atom pairs.
        2. `df_atoms`: A table containing atom indices, types, positions (grid units), π charges, and colors.
        3. `df_bonds`: A table listing all bonds between atoms with atom labels.
        4. `df_summary`: A summary table of molecular descriptors, such as total π-electron energy, HOMO-LUMO gap, chemical potential, hardness, softness, and electrophilicity.

        The resulting DataFrames are stored as attributes for export or display in other parts of the application.
    
        Returns
        -------
        None
        (The method updates internal attributes: `df_bond_orders_matrix`, `df_atoms`, `df_bonds`, and `df_summary`)
    
        Notes
        -----
        - Atom labels are formatted as '<AtomType><Index>' (e.g., 'C·1', 'O:2').
        - Bond orders are rounded to 3 decimal places.
        - π charges are computed from Hückel analysis results and rounded to 3 decimal places.
        - The summary DataFrame expresses all energy-related descriptors in units of β (or multiples/fractions thereof).
        - The method assumes that the Hückel analysis (and charge/bond order computation) has already been performed.
    
        Example
        -------
        >>> mol = MoleculeDrawer(root)
        >>> mol.build_dataframes()
        
        The following attributes are now available:
        - mol.df_bond_orders_matrix
        - mol.df_atoms
        - mol.df_bonds
        - mol.df_summary
        """

        # DataFrame pour les indices de liaison (π-bond orders)
        atom_labels = [f"{node.atom_type}{idx + 1}" for idx, node in enumerate(self.nodes)]
        n_atoms = len(self.nodes)
        bond_order_matrix = np.zeros((n_atoms, n_atoms))
        for (i, j), val in self.bond_orders.items():
            bond_order_matrix[i, j] = val
            bond_order_matrix[j, i] = val  # rendre symétrique
        self.df_bond_orders_matrix = pd.DataFrame(
            bond_order_matrix,
            index=atom_labels,
            columns=atom_labels
        ).round(3)
    
        # DataFrame pour les coordonnées des atomes
        atom_data = []
        for idx, node in enumerate(self.nodes):
            atom_data.append({
                "Atom": f"{node.atom_type}{idx + 1}",
                "Type": node.atom_type,
                "X (grid units)": node.x,
                "Y (grid units)": node.y,
                "π Charge": round(self.charges[idx], 3),
                "Color": HuckelParameters.ATOM_COLORS.get(node.atom_type, 'gray')                
            })
        self.df_atoms = pd.DataFrame(atom_data)
    
        # DataFrame pour les liaisons
        bond_data = []
        for idx, (i, j) in enumerate(self.bonds):
            atom1 = f"{self.nodes[i].atom_type}{i + 1}"
            atom2 = f"{self.nodes[j].atom_type}{j + 1}"
            bond_data.append({
                "Bond #": idx + 1,
                "Atom 1": atom1,
                "Atom 2": atom2
            })
        self.df_bonds = pd.DataFrame(bond_data) 
    
        # Résumé des propriétés dans une 4e feuille
        alpha = self.alpha_value
        beta = self.beta_value

        if self.homo_lumo_gap is not None:
            HLgap = f"{abs(self.homo_lumo_gap):.2f}"
            muprt = f"{self.mu:.2f}"
            etaprt = f"{abs(self.eta):.2f}"
            softnessprt = f"{abs(self.softness):.2f}"
            omegaprt = f"{self.omega:.2f}"
        else:
            HLgap = self.homo_lumo_gap
            muprt = self.mu
            etaprt = self.eta
            softnessprt = self.softness
            omegaprt = self.omega

        
        summary_data = {
            "Descriptor": [
                "Total π-electron energy",
                "Total number of π electrons",
                "Atomization energy",
                "Atomization energy per π atom",
                "HOMO-LUMO gap",
                "μ (chemical potential)",
                "η (hardness)",
                "S (softness)",
                "ω (electrophilicity)"
            ],
            "Value": [
                f"{self.alpha_part:.0f}α + {self.beta_part:.2f}β",
                f"{self.total_pi_electrons}",
                f"{self.atomization_energy:.2f}",
                f"{self.atomization_energy_per_atom:.2f}",
                f"{HLgap}",
                f"{muprt}",
                f"{etaprt}",
                f"{softnessprt}",
                f"{omegaprt}",
            ],
            "Unit": [
                f"",
                f"",
                f"β",
                f"β",
                f"|β|",
                f"β",
                f"|β|",
                f"|β|^-1",
                f"β^2",
            ]
        }
        self.df_summary = pd.DataFrame(summary_data).set_index("Descriptor")

    def save_data(self):
        if self.df is None:
            messagebox.showerror("Error", "You must run Huckel before saving the data.")
            return
        success = self.save_dataframe_as_xlsx(self.master)
        if success:
            messagebox.showinfo("Success", "The data have been successfully saved.")
        else:
            messagebox.showinfo("Cancelled", "The data have not been saved.")
    
    def save_dataframe_as_xlsx(self, win):
        """
        Export multiple pandas DataFrames to an Excel workbook with styled headers.
    
        This method saves the following DataFrames to separate sheets in a single Excel file:
        - `self.df`: Molecular orbital coefficients (sheet 'MO Coefficients')
        - `self.df_bond_orders_matrix`: π-bond order matrix (sheet 'π-bond orders')
        - `self.df_atoms`: Atom-level data (sheet 'Atoms')
        - `self.df_bonds`: Bond list (sheet 'Bond List')
        - `self.df_summary`: Molecular descriptors (sheet 'Descriptors')
    
        The user is prompted to select a save location via a file dialog. After saving,
        the header row of the 'MO Coefficients' sheet is styled with text wrapping.
    
        Parameters
        ----------
        win : tkinter.Tk or tkinter.Toplevel
            The parent window for the file dialog.
    
        Returns
        -------
        bool
            True if the file was successfully saved, False otherwise (e.g., if the user cancels).
    
        Notes
        -----
        - If no data is available (DataFrames are None), a warning dialog is shown.
        - The Excel file is saved with the '.xlsx' extension, and the default filename is
          based on `self.safe_project_name`.
        """
        from tkinter import filedialog, messagebox
        import openpyxl
    
        if self.df is None or self.df_atoms is None:
            messagebox.showwarning("Warning", "No DataFrame available to save.")
            return
    
        path = filedialog.asksaveasfilename(
            defaultextension=".xlsx",
            filetypes=[("Excel files", "*.xlsx"), ("All Files", "*.*")],    
            initialfile=f"{self.safe_project_name}.xlsx"
        )
    
        if not path:
            return False
    
        with pd.ExcelWriter(path, engine='openpyxl') as writer:
            self.df.to_excel(writer, sheet_name='MO Coefficients')
            self.df_bond_orders_matrix.to_excel(writer, sheet_name='π-bond orders')
            self.df_atoms.to_excel(writer, sheet_name='Atoms', index=False)
            self.df_bonds.to_excel(writer, sheet_name='Bond List', index=False)
            self.df_summary.to_excel(writer, sheet_name='Descriptors', index=True)
    
        # Ouvre le fichier avec openpyxl pour modifier les styles
        wb = openpyxl.load_workbook(path)
        ws = wb['MO Coefficients']
    
        # Applique wrap_text=True pour la première ligne (en-têtes)
        for cell in ws[1]:
            cell.alignment = openpyxl.styles.Alignment(wrap_text=True)
    
        wb.save(path)
        return True

    def save_dataframe_as_pdf(self, win):
    
        path = filedialog.asksaveasfilename(
            defaultextension=".pdf",
            filetypes=[
                ("PDF files", "*.pdf"),
                ("All Files", "*.*")
            ],
            parent=win,
            initialfile=f"{self.safe_project_name}.pdf"
        )
        
        if path:
            fig_size_cm = 14
            dpi = 300
            fig_size_inches = fig_size_cm / 2.54
        
            fig, ax = plt.subplots(figsize=(fig_size_inches, fig_size_inches), dpi=dpi)
            ax.axis('off')
        
            dpi = 300
            TARGET_BOND_LENGTH_PX = dpi / 2.54  # ~118.11 pixels par cm à 300 dpi
        
            if self.bonds:
                lengths = []
                for i, j in self.bonds:
                    xi, yi = self.nodes[i].x, self.nodes[i].y
                    xj, yj = self.nodes[j].x, self.nodes[j].y
                    dist = ((xi - xj)**2 + (yi - yj)**2)**0.5
                    lengths.append(dist)
                avg_length = sum(lengths) / len(lengths)
            else:
                avg_length = 100

            scale_factor_base = TARGET_BOND_LENGTH_PX / avg_length if avg_length != 0 else 1

            xs = [node.x for node in self.nodes]
            ys = [node.y for node in self.nodes]
            min_x, max_x = min(xs), max(xs)
            min_y, max_y = min(ys), max(ys)
            width_molecule = (max_x - min_x) * scale_factor_base
            height_molecule = (max_y - min_y) * scale_factor_base
            
            # print(f"[DEBUG] avg_length: {avg_length:.2f} units")
            # print(f"[DEBUG] Molecule bounds: x=({min_x:.1f}, {max_x:.1f}), y=({min_y:.1f}, {max_y:.1f})")
        
            HMOViewerParameters.CANVAS_WIDTH_px = fig_size_inches * dpi
        
            scale_limit = min(
                (HMOViewerParameters.CANVAS_WIDTH_px * 0.8) / width_molecule,
                (HMOViewerParameters.CANVAS_WIDTH_px * 0.8) / height_molecule
            ) if width_molecule > 0 and height_molecule > 0 else 1
        
            final_scale = scale_factor_base * min(scale_limit, 1)
        
            # print(f"[DEBUG] Molecule raw size (after base scale): {width_molecule:.1f} x {height_molecule:.1f} px")
            # print(f"[DEBUG] Scaling factor (liaison target ~1 cm): {scale_factor_base:.2f}")
            # print(f"[DEBUG] Downscale factor (max 14x14 cm): {min(scale_limit, 1):.2f}")
            # print(f"[DEBUG] Final scale: {final_scale:.2f}")
        
            center_x = (max_x + min_x) / 2
            center_y = (max_y + min_y) / 2
        
            # === 2️⃣ Dessiner la molécule (liaisons + atomes + labels C1, O2, etc.) ===
            # Liaisons
            for bond in self.bonds:
                i, j = bond
                xi, yi = self.nodes[i].x, self.nodes[i].y
                xj, yj = self.nodes[j].x, self.nodes[j].y
                xi_scaled = (xi - center_x) * final_scale
                yi_scaled = (yi - center_y) * final_scale
                xj_scaled = (xj - center_x) * final_scale
                yj_scaled = (yj - center_y) * final_scale
                ax.plot([xi_scaled, xj_scaled], [yi_scaled, yj_scaled], color='black')
    
            # Atomes
            for idx, node in enumerate(self.nodes):
                x_scaled = (node.x - center_x) * final_scale
                y_scaled = (node.y - center_y) * final_scale
                color = HuckelParameters.ATOM_COLORS.get(node.atom_type, 'gray')
                ax.plot(x_scaled, y_scaled, 'o', color=color, markersize=20)
                label = f"{node.atom_type}{idx+1}"
                ax.text(x_scaled, y_scaled, label, ha='center', va='center', color='white', fontsize=8)
        
            ax.set_xlim(-HMOViewerParameters.CANVAS_WIDTH_px/2, HMOViewerParameters.CANVAS_WIDTH_px/2)
            ax.set_ylim(-HMOViewerParameters.CANVAS_WIDTH_px/2, HMOViewerParameters.CANVAS_WIDTH_px/2)


            # Ajoute les propriétés en haut du PDF
            coordTxt_x = 0.1
            coordTxt_y = 0.95
            nameOfMol = Path(path).stem
            print(nameOfMol)
            fig.text(coordTxt_x, coordTxt_y, nameOfMol, fontsize=12, color = "#128d85", fontweight="bold", va='top', ha='left', fontname="DejaVu Sans")
            fig.text(coordTxt_x + 0.8, coordTxt_y, self.summary_text, fontsize=10, va='top', ha='left', fontname="DejaVu Sans")
    
            # === 3️⃣ Ajouter le DataFrame sous forme de texte monospace ===
            
            # Générer le texte propre du DataFrame
            df_text = self.df.to_string(index=False)
            
            # Afficher le texte sous la molécule
            ax.text(
                0, -HMOViewerParameters.CANVAS_WIDTH_px/2 + 50,  # placé sous la molécule ; ajuste +50 si besoin
                df_text,
                fontsize=10,
                fontfamily='monospace',
                ha='center',
                va='bottom',
                wrap=True
            )
            plt.savefig(path, bbox_inches='tight')
            plt.close()
            
    def show_numerical_data(self):
        if self.df is None:
            messagebox.showerror("Error", "You must run Huckel before viewing the data.")
            return
        self.show_dataframe_in_window()
    
    def show_dataframe_in_window(self):
        """
        Display molecular orbital coefficients, a molecular diagram, and key descriptors in a new window.
    
        This method creates a Tkinter `Toplevel` window that presents:
        
        [1] A scaled 2D visualization of the molecule (atoms and bonds), centered in a canvas.
        
        [2] A summary section with computed descriptors:
            - Total π-electron energy,
            - Number of π electrons,
            - Atomization energy (absolute and per atom),
            - HOMO-LUMO gap,
            - Chemical potential (μ),
            - Hardness (η),
            - Softness (S),
            - Electrophilicity index (ω).
            
        [3] A scrollable table that displays the molecular orbital (MO) coefficients.
            - The table is split into two Treeviews: one for the index column, one for the data.
            - Scrollbars (vertical and horizontal) are shared between the two views for smooth scrolling.
    
        The visualization includes:
        - Atoms drawn as colored circles (color based on `ATOM_COLORS`),
        - Atom labels with their type and index,
        - Bonds drawn as black lines.
    
        Features:
        ----------
        - The molecular diagram is auto-scaled to fit the canvas, with debug prints for scaling details.
        - The summary text is stored in `self.summary_text` for potential later use (e.g., export).
        - If a previous window (`self.om_window`) is still open, it is closed before opening the new one.
        - An Escape key binding allows closing the window quickly.
    
        Notes:
        ------
        - The descriptors use `self.alpha_value` and `self.beta_value` for normalization.
        - Treeviews are laid out side-by-side: the index column remains fixed, while data columns can expand.
        - This window is meant as an interactive, visual companion to numerical data (e.g., after Hückel analysis).
    
        Shortcut:
        ---------
        - Pressing `Esc` closes the window.
        """        
        
        alpha = self.alpha_value
        beta = self.beta_value

        win = tk.Toplevel()
        win.title("Hückel Molecular Orbital Coefficients")
        win.geometry("1000x900")

        # === Frame principal (vertical) ===
        frame_main = ttk.Frame(win)
        frame_main.pack(fill=tk.BOTH, expand=True)
        
        if self.om_window is not None and self.om_window.winfo_exists():
            self.om_window.destroy()
        self.om_window = win

        w = 400
        h = 400
        canvas = tk.Canvas(frame_main, width=w, height=h, bg='white')
        canvas.pack(side=tk.TOP, pady=10)
    
        # === 1️⃣ Calcul de l'échelle ===
        TARGET_BOND_LENGTH = 60  # pixels par liaison
    
        # Longueur moyenne des liaisons
        if self.bonds:
            lengths = []
            for i, j in self.bonds:
                xi, yi = self.nodes[i].x, self.nodes[i].y
                xj, yj = self.nodes[j].x, self.nodes[j].y
                dist = ((xi - xj)**2 + (yi - yj)**2)**0.5
                lengths.append(dist)
            avg_length = sum(lengths) / len(lengths)
        else:
            avg_length = 100  # fallback si aucune liaison
    
        # Premier facteur d'échelle : 1 liaison ≈ 60 px
        scale_factor_base = TARGET_BOND_LENGTH / avg_length if avg_length != 0 else 1
    
        # Dimensions de la molécule
        xs = [node.x for node in self.nodes]
        ys = [node.y for node in self.nodes]
        min_x, max_x = min(xs), max(xs)
        min_y, max_y = min(ys), max(ys)
        width_molecule_raw = (max_x - min_x)
        height_molecule_raw = (max_y - min_y)
    
        # Adaptation pour tenir dans le canvas
        scale_limit = min(
            (w * 0.8) / (width_molecule_raw * scale_factor_base),
            (h * 0.8) / (height_molecule_raw * scale_factor_base)
        ) if width_molecule_raw > 0 and height_molecule_raw > 0 else 1

        if scale_limit < 1:
            final_scale = scale_factor_base * scale_limit
        else:
            final_scale = scale_factor_base

        # print(f"[DEBUG] avg_length: {avg_length:.2f} px")
        # print(f"[DEBUG] width_raw: {width_molecule_raw:.1f}, height_raw: {height_molecule_raw:.1f}")
        # print(f"[DEBUG] scale_factor_base: {scale_factor_base:.2f}")
        # print(f"[DEBUG] scale_limit: {scale_limit:.2f}")
        # print(f"[DEBUG] final_scale: {final_scale:.2f}")
    
        # Centrage
        center_x = (max_x + min_x) / 2
        center_y = (max_y + min_y) / 2
        offset_x = w / 2
        offset_y = h / 2
        
        # === Dessiner les liaisons ===
        for i, j in self.bonds:
            xi = (self.nodes[i].x - center_x) * final_scale + offset_x
            yi = (self.nodes[i].y - center_y) * final_scale + offset_y
            xj = (self.nodes[j].x - center_x) * final_scale + offset_x
            yj = (self.nodes[j].y - center_y) * final_scale + offset_y
            canvas.create_line(xi, yi, xj, yj, fill='black', width=2)
    
        # === 2️⃣ Dessiner les atomes + numéros ===
        for idx, node in enumerate(self.nodes):
            x = (node.x - center_x) * final_scale + offset_x
            y = (node.y - center_y) * final_scale + offset_y
            atom_label = f"{node.atom_type}{idx+1}"
            color = HuckelParameters.ATOM_COLORS.get(node.atom_type, 'gray')
            canvas.create_oval(x-15, y-15, x+15, y+15, fill=color)
            canvas.create_text(x, y, text=atom_label, fill='white')

        # === Résumé des propriétés ===

        summary_title = tk.Label(
            frame_main,
            text="🔬 Energies and various descriptors",
            font=("DejaVu Sans", 10, "bold"),
            fg="#8f0000",
            justify="left"
        )
        summary_title.pack(pady=(5, 0))

        if self.homo_lumo_gap is not None:
            HLgap = f"{abs(self.homo_lumo_gap):.2f}|β|"
            muprt = f"{self.mu:.2f}β"
            etaprt = f"{abs(self.eta):.2f}|β|"
            softnessprt = f"{abs(self.softness):.2f}/|β|"
            omegaprt = f"{self.omega:.2f}β"
        else:
            HLgap = self.homo_lumo_gap
            muprt = self.mu
            etaprt = self.eta
            softnessprt = self.softness
            omegaprt = self.omega

        print(f"[DEBUG]. {self.alpha_part=}   {self.beta_part=}")

        summary_text = (
            f"Total π-electron energy: {self.alpha_part:.0f}α + {self.beta_part:.2f}β\n"
            f"Number of π electrons: {self.total_pi_electrons}\n"
            f"Atomization energy: {self.atomization_energy:.2f}β\n"
            f"Atomization energy per π atom: {self.atomization_energy_per_atom:.2f}β\n"
            f"HOMO-LUMO gap: {HLgap}\n"
            f"Electronic potential μ: {muprt}\n"
            f"Chemical hardness η: {etaprt}\n"
            f"Chemical softness S: {softnessprt}\n"
            f"Electrophilicity index ω: {omegaprt}\n"
        )
        self.summary_text = summary_text  # sauvegarde pour le PDF
        summary_label = tk.Label(frame_main, text=summary_text, font=("DejaVu Sans", 10), justify="left")
        summary_label.pack(pady=5)    
        
        # === 3️⃣ Table (Index figé + données avec scrollbars) ===
        frame_table = tk.Frame(frame_main)
        frame_table.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        columns_index = ['Index']
        columns_data = list(self.df.columns)
        col_ids_data = [str(i) for i in range(len(columns_data))]
        
        # --- Treeview pour l'Index ---
        tree_index = ttk.Treeview(frame_table, columns=columns_index, show='headings', height=15)
        tree_index.heading('Index', text='Index')
        tree_index.column('Index', width=50, anchor='center')
        
        # --- Treeview pour les données ---
        tree_data = ttk.Treeview(frame_table, columns=col_ids_data, show='headings', height=15)
        
        # Configurer les colonnes des données
        for i, col in enumerate(columns_data):
            tree_data.heading(str(i), text=col)
            tree_data.column(str(i), width=100, anchor='center')
        
        # Remplir les deux Treeviews
        for idx in self.df.index:
            # Insert dans Index
            tree_index.insert('', 'end', values=[idx])
            # Insert dans data
            tree_data.insert('', 'end', values=list(self.df.loc[idx]))
        
        # --- Scrollbars ---
        # Scrollbar verticale partagée
        scrollbar_y = ttk.Scrollbar(frame_table, orient=tk.VERTICAL)
        scrollbar_y.config(command=lambda *args: (tree_index.yview(*args), tree_data.yview(*args)))
        tree_index.configure(yscrollcommand=scrollbar_y.set)
        tree_data.configure(yscrollcommand=scrollbar_y.set)
        
        # Scrollbar horizontale pour tree_data seulement
        scrollbar_x = ttk.Scrollbar(frame_table, orient=tk.HORIZONTAL, command=tree_data.xview)
        tree_data.configure(xscrollcommand=scrollbar_x.set)
        
        # --- Grid layout ---
        tree_index.grid(row=0, column=0, sticky='nsew')
        tree_data.grid(row=0, column=1, sticky='nsew')
        scrollbar_y.grid(row=0, column=2, sticky='ns')
        scrollbar_x.grid(row=1, column=1, sticky='ew')
        
        # --- Configure la grille pour l'expansion ---
        frame_table.grid_rowconfigure(0, weight=1)
        frame_table.grid_columnconfigure(0, weight=0)  # Index column fixe
        frame_table.grid_columnconfigure(1, weight=1)  # Data colonnes extensibles

    
        # === 4️⃣ Bouton Close ===
        frame_buttons = tk.Frame(win)
        frame_buttons.pack(pady=10)
        
        def on_escape(event):
            # print("Escape key pressed. Closing the app.")
            win.destroy()
        win.bind('<Escape>', on_escape)
        
        btn_close = tk.Button(frame_buttons, text="Close", command=win.destroy)
        btn_close.pack()

    
    # Hook to existing GUI
    def on_run_huckel(self):
        """
        Executes the full Hückel analysis workflow.
    
        This method performs the following steps:

        1️⃣ Runs the Hückel analysis on the current molecular structure.

        2️⃣ Prompts the user for a project name (default: "my_molecule") to label results and files.

        3️⃣ Sanitizes the project name for safe file usage.
        
        4️⃣ Builds internal DataFrames:
            - Molecular orbital coefficients,
            - Atom-level data,
            - Bond list and bond orders,
            - Summary descriptors (energies, HOMO-LUMO gap, etc.).
            
        5️⃣ (Optional) Displays debug printouts for DataFrame shapes.
        
        6️⃣ Opens a new Tkinter window (`HMOViewer`) to visualize the molecular orbitals, atoms, bonds, and descriptors graphically.
    
        Notes
        -----
        - The visualizer window (`HMOViewer`) is initialized with all computed data for interactive exploration.
        - The project name is reused as a base name for saving exports (Excel, PDF, etc.).
        - Debug information about DataFrame shapes is printed to the console for traceability.
        """

        self.run_huckel_analysis()
        from tkinter import simpledialog

        # Determine initial value for the project name
        initial_value = self.project_name if self.project_name else self.default_project_name
    
        self.project_name = simpledialog.askstring(
            "Project name",
            "Enter a name for your project:",
            initialvalue=initial_value
        )
        if not self.project_name:
            self.project_name = self.default_project_name
        # print(f"[INFO] Project name set to: {self.project_name}")
    
        self.safe_project_name = self.sanitize_filename(self.project_name)
        self.build_dataframes()
        self.molecule_is_new = False  # Reset after analysis
        
        # self.show_dataframe_in_window()
        # # === Ouvre la visualisation des OM ===
        # print("[DEBUG] df_MOs shape:", self.df.shape)
        # print("[DEBUG] df_atoms shape:", self.df_atoms.shape)
        # print("[DEBUG] df_bonds shape:", self.df_bonds.shape)
        # print("[DEBUG] df_summary shape:", self.df_summary.shape)
        viewer_win = tk.Toplevel(self.master)
        self.mo_viewer_window = viewer_win
        viewer = HMOViewer(
            viewer_win,
            df_MOs=self.df,
            df_atoms=self.df_atoms,
            df_bonds=self.df_bonds,
            df_descriptors=self.df_summary,
            project_name=self.safe_project_name
        )

# =============================================================================================================================================

class ToolTip(object):
    """
    A class to create and manage tooltips for Tkinter widgets.

    This class attaches a tooltip (a small pop-up window displaying text) to any Tkinter widget. 
    The tooltip appears when the mouse hovers over the widget and disappears when the mouse leaves.

    The tooltip text is rendered using a custom Open Sans Regular font loaded from a TrueType Font (TTF) 
    file located at 'Fonts/OpenSans/static/OpenSans-Regular.ttf'. This ensures a consistent and modern 
    look across different platforms and display environments.

    Parameters
    ----------
    widget : tk.Widget
        The Tkinter widget to which the tooltip will be attached.
    text : str
        The text content to display in the tooltip.

    Attributes
    ----------
    widget : tk.Widget
        The widget associated with the tooltip.
    text : str
        The message displayed when hovering over the widget.
    tooltip_window : tk.Toplevel or None
        The pop-up window that displays the tooltip text (created when the mouse enters the widget).
    OpenSansReg_font : tkFont.Font
        The custom font object used to render the tooltip text, loaded from OpenSans-Regular.ttf.

    Main Methods
    -------
    show_tooltip(event=None)
        Creates and displays the tooltip near the mouse cursor using the Open Sans font.
    hide_tooltip(event=None)
        Destroys the tooltip window when the mouse leaves the widget.
    """

    def __init__(self, widget, text='widget info'):
        font_path = resource_path("Fonts/OpenSans/static/OpenSans-Regular.ttf")
        self.OpenSansReg_font = font.Font(family=font_path)
        self.widget = widget
        self.text = text
        self.tipwindow = None
        self.id = None
        self.x = self.y = 0
        widget.bind("<Enter>", self.enter)
        widget.bind("<Leave>", self.leave)
        widget.bind("<Motion>", self.motion)

    def enter(self, event=None):
        self.schedule()

    def leave(self, event=None):
        self.unschedule()
        self.hidetip()

    def motion(self, event=None):
        if self.tipwindow:
            x, y, _, _ = self.widget.bbox("insert")
            x = event.x_root + 20
            y = event.y_root + 10
            self.tipwindow.wm_geometry(f"+{x}+{y}")

    def schedule(self):
        self.unschedule()
        self.id = self.widget.after(500, self.showtip)

    def unschedule(self):
        id_ = self.id
        self.id = None
        if id_:
            self.widget.after_cancel(id_)

    def showtip(self, event=None):
        font10 = self.OpenSansReg_font.copy()
        font10.configure(size=10)
        if self.tipwindow or not self.text:
            return
        x, y, _, _ = self.widget.bbox("insert")
        x = self.widget.winfo_pointerx() + 20
        y = self.widget.winfo_pointery() + 10
        self.tipwindow = tw = tk.Toplevel(self.widget)
        tw.wm_overrideredirect(True)
        tw.wm_geometry(f"+{x}+{y}")
        label = tk.Label(
            tw, text=self.text, justify=tk.LEFT,
            background="#c5c5c5", relief=tk.SOLID, borderwidth=1,
            font=font10
        )
        label.pack(ipadx=4)

    def hidetip(self):
        tw = self.tipwindow
        self.tipwindow = None
        if tw:
            tw.destroy()

# =============================================================================================================================================

if __name__ == "__main__":
    root = tk.Tk()
    app = MoleculeDrawer(root)
    root.mainloop()

