import tkinter as tk
from tkinter import filedialog, messagebox, font, ttk
from PIL import Image, ImageTk
import numpy as np
import pandas as pd
from pathlib import Path
import re

GRID_SIZE = 30
ATOM_RADIUS = 14
HIGHLIGHT_RADIUS = 12

ATOM_COLORS = {
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

ATOM_OPTIONS = [
    'B▯', 'C·', 'N·', 'N:', 'N+·', 'O·', 'O:',  'O+·', 'F:',
    'Si·', 'P·', 'P:', 'S·', 'S:', 'Cl:', 'Br:',
    'Me:'
]

Huckel_atomic_parameters = {
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

Huckel_kXY_parameters = {
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
#### formerly independent visu_v7.py

# Frames for MO visualization
FRAME_WIDTH = 700
FRAME_HEIGHT = 450

FRAME_LUMO_X = 350
FRAME_LUMO_Y = 30

FRAME_HOMO_X = FRAME_LUMO_X
FRAME_HOMO_Y = FRAME_LUMO_Y + FRAME_HEIGHT + 40  # espacement vertical de 20 px

BOTTOM_OF_HOMO = FRAME_HOMO_Y + FRAME_HEIGHT
MARGIN_BOTTOM = 50
CANVAS_HEIGHT = BOTTOM_OF_HOMO + MARGIN_BOTTOM

CANVAS_WIDTH = FRAME_LUMO_X + FRAME_WIDTH + 40

# MOD Diagram on the left
DIAG_X = 120
DIAG_Y_TOP = 50
DIAG_Y_BOTTOM = 600
DIAG_WIDTH_UNIT = 30

SHRINK_FACTOR = 0.85
LOBE_OFFSET = 5
TARGET_BOND_PX = 60

# print(f"[DEBUG] {CANVAS_HEIGHT=}")
# print(f"[DEBUG] {CANVAS_WIDTH=}")
# print(f"[DEBUG]")

Path2Imgs = Path("DesignMOdiagram")

class HMOViewer:
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
        self.canvas = tk.Canvas(master, width=CANVAS_WIDTH, height=CANVAS_HEIGHT, bg='white')
        self.canvas.pack()

        self.master.bind('<Escape>', self.on_escape)

        # === 🟦 Ajout des boutons sous le diagramme ===
        button_frame = tk.Frame(self.master)
        button_frame.pack(pady=10)  # petit espacement vertical
        
        # Bouton Skeleton Only
        self.skeleton_button = tk.Button(button_frame, text="Skeleton Only", command=self.toggle_skeleton)
        self.skeleton_button.pack(side=tk.LEFT, padx=5)
        
        # Bouton Save as PNG
        self.save_button = tk.Button(button_frame, text="Save as PNG", command=self.save_canvas_as_png)
        self.save_button.pack(side=tk.LEFT, padx=5)

        # Bouton Close
        self.close_button = tk.Button(button_frame, text="Close", command=self.master.destroy)
        self.close_button.pack(side=tk.LEFT, padx=5)

        self.draw_layout()
        self.draw_energy_levels()
        self.display_default_homo_lumo()
        
    def on_escape(self, event):
        print("Escape key pressed. Closing the app.")
        self.master.destroy()

    def load_images(self):
        self.img_1e = ImageTk.PhotoImage(Image.open(Path2Imgs / "1e.png").resize((8, 40)))
        self.img_2e = ImageTk.PhotoImage(Image.open(Path2Imgs / "2e.png").resize((23, 41)))
        self.img_energy_level = ImageTk.PhotoImage(Image.open(Path2Imgs / "energy_level.png").resize((40, 5)))

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
        print(f"[DEBUG] Skeleton only mode: {self.show_atom_labels}")
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
        self.canvas.create_text(
            3*FRAME_LUMO_X/4,
            10,
            text=self.project_name,
            font=('DejaVu Sans', 12, 'bold'),
            fill='#128d85',
            anchor='center'
        )

        self.canvas.create_rectangle(FRAME_LUMO_X, FRAME_LUMO_Y, FRAME_LUMO_X + FRAME_WIDTH, FRAME_LUMO_Y + FRAME_HEIGHT, outline='black')
        self.canvas.create_text(FRAME_LUMO_X + FRAME_WIDTH/2, FRAME_LUMO_Y - 20, text='Virtual MOs', font=('DejaVu Sans', 14, 'bold'))

        self.canvas.create_rectangle(FRAME_HOMO_X, FRAME_HOMO_Y, FRAME_HOMO_X + FRAME_WIDTH, FRAME_HOMO_Y + FRAME_HEIGHT, outline='black')
        self.canvas.create_text(FRAME_HOMO_X + FRAME_WIDTH/2, FRAME_HOMO_Y - 20, text='Occupied MOs', font=('DejaVu Sans', 14, 'bold'))
        descriptors_height = 140
        # print(f"DEBUG y0 = {CANVAS_HEIGHT - 10=}")
        # print(f"DEBUG width = {FRAME_LUMO_X - 10=}")
        # print(f"DEBUG height = {CANVAS_HEIGHT - 10 - DIAG_Y_BOTTOM - descriptors_height =}")
        self.skeleton_x0 = 10  # ou l'endroit où tu dessines le squelette
        self.skeleton_y0 = DIAG_Y_BOTTOM + descriptors_height  # sous le diagramme HOMO
        self.skeleton_width = FRAME_LUMO_X - 10
        self.skeleton_height = CANVAS_HEIGHT - 10 - DIAG_Y_BOTTOM - descriptors_height
        self.draw_skeleton_overview(
                                        x0 = self.skeleton_x0,  # ajuste en fonction de la largeur que tu veux
                                        y0 = self.skeleton_y0,  # un peu sous le diagramme
                                        width = self.skeleton_width,
                                        height = self.skeleton_height,
                                    )

    def draw_energy_scale_and_descriptors(self, min_e, max_e, scale, energies_sorted):
        """Dessine l'échelle énergétique, les repères horizontaux et les descripteurs sous le diagramme."""
        print(f"draw_energy_scale_and_descriptors {max_e=}")
        print(f"draw_energy_scale_and_descriptors {min_e=}")
        # === 🟦 ÉCHELLE ÉNERGÉTIQUE VERTICALE ===
        start_y = DIAG_Y_BOTTOM - (0.3) * scale
        end_y = DIAG_Y_BOTTOM - ( (max_e - 0.2) - min_e ) * scale

        # Prendre un des groupes (le premier trouvé)
        for energy in energies_sorted:
            group = self.energy_groups[energy]
            n_degenerate = len(group)
            total_width = n_degenerate * (self.img_energy_level.width() + 10)
            start_x = DIAG_X + (self.img_energy_level.width() // 2) - total_width // 2 + (self.img_energy_level.width() // 2)
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
            y_pos = DIAG_Y_BOTTOM - ( current - min_e ) * scale
        
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
    
        y_text = DIAG_Y_BOTTOM + 60
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
        energies_sorted = sorted(self.energy_groups.keys(), reverse=True)
        print("Sorted energies (β values):", energies_sorted)
        min_e, max_e = max(energies_sorted), min(energies_sorted)
        print(f"draw_energy_levels {max_e=}")
        print(f"draw_energy_levels {min_e=}")        
        scale = (DIAG_Y_BOTTOM - DIAG_Y_TOP) / (max_e - min_e + 1e-6)

        self.draw_energy_scale_and_descriptors(min_e, max_e, scale, energies_sorted)

        for energy in energies_sorted:
            group = self.energy_groups[energy]
            n_degenerate = len(group)
            total_width = n_degenerate * (self.img_energy_level.width() + 10)
            start_x = DIAG_X + (self.img_energy_level.width() // 2) - total_width // 2 + (self.img_energy_level.width() // 2)

            y_pos = DIAG_Y_BOTTOM - (energy - min_e) * scale

            for i, idx in enumerate(group):
                x_pos = start_x + i * (self.img_energy_level.width() + 10)
                print(f"Energy {energy} | idx {idx} | Degeneracy {n_degenerate} | x_pos {x_pos}")

                img = self.canvas.create_image(x_pos, y_pos, image=self.img_energy_level)

                occ = self.occupations[idx]
                if occ == 2:
                    self.canvas.create_image(x_pos, y_pos, image=self.img_2e, state="disabled")
                elif occ == 1:
                    self.canvas.create_image(x_pos, y_pos, image=self.img_1e, state="disabled")

                self.canvas.tag_bind(img, '<Button-1>', lambda event, idx=idx: self.on_level_click(idx))
        
    def draw_skeleton_overview(self, x0, y0, width, height):
        """Dessine la molécule complète dans une petite fenêtre sous le diagramme d'OM."""
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
        """Affiche la représentation d'une OM dans le cadre approprié."""
    
        # === 🖼️ Cadre d'affichage ===
        frame_x = FRAME_HOMO_X if occupied else FRAME_LUMO_X
        frame_y = FRAME_HOMO_Y if occupied else FRAME_LUMO_Y
    
        # Effacer & redessiner le fond
        self.canvas.create_rectangle(frame_x, frame_y, frame_x + FRAME_WIDTH, frame_y + FRAME_HEIGHT, fill='white', outline='black')
    
        # === 🔖 Titre avec étiquette de l'OM ===
        nrj = f"α - {np.abs(self.energies[idx]):.2f}β" if self.energies[idx] < 0 else f"α + {np.abs(self.energies[idx]):.2f}β"
        om_label = f"OM #{idx+1} | {nrj} | {self.occupations[idx]}e"
        self.canvas.create_text(frame_x + FRAME_WIDTH/2, frame_y + 15, text=om_label, font=('DejaVu Sans', 10, 'italic'))
    
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
        
        # Le scale qui donnerait des liaisons de ~TARGET_BOND_PX
        scale_bond_based = TARGET_BOND_PX / self.mean_bond_length
        
        # Est-ce que ça passe dans le frame ?
        required_width = mol_width * scale_bond_based + margin
        required_height = mol_height * scale_bond_based + margin
        
        if required_width > FRAME_WIDTH or required_height > FRAME_HEIGHT:
            # Downscale forcé
            scale_x = (FRAME_WIDTH - margin) / mol_width if mol_width != 0 else 1
            scale_y = (FRAME_HEIGHT - margin) / mol_height if mol_height != 0 else 1
            scale = min(scale_x, scale_y)
            # print(f"[DEBUG] Downscaling: target scale {scale_bond_based:.2f} too big, reduced to {scale:.2f}")
        else:
            scale = scale_bond_based
            # print(f"[DEBUG] Using target scale: {scale:.2f}")
        
        # Offset pour bien centrer la molécule
        offset_x = frame_x + FRAME_WIDTH/2 - ((min_x + max_x) / 2) * scale
        offset_y = frame_y + FRAME_HEIGHT/2 - ((min_y + max_y) / 2) * scale
        
        # Taille max des lobes adaptée à la taille réelle après scaling
        # Trouver le coefficient max global pour normaliser
        max_coef_global = np.nanmax(np.abs(self.df_MOs.values))
        
        # Fixer la taille max du lobe : quand coef == max_coef_global, le rayon est 90% de longueur liaison/2
        desired_radius = 0.90 * self.mean_bond_length * scale / 2
        scale_lobe_factor = desired_radius / max_coef_global
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

        dx = LOBE_OFFSET * scale  # petit décalage droite
        dy = -LOBE_OFFSET * scale / 15  # petit décalage haut
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
            dy = LOBE_OFFSET * scale / 10
            color_front = 'red' if coef > 0 else 'blue'
    
            # Lobe devant (au-dessus de la liaison)
            self.canvas.create_oval(
                offset_x + x * scale - size, offset_y + y * scale + dy - size,
                offset_x + x * scale + size, offset_y + y * scale + dy + size,
                fill=color_front, stipple='gray50', outline=''
            )

    
    def display_default_homo_lumo(self):
        print("Looking for HOMO/LUMO...")
    
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

        print(f"HOMO: idx {homo_idx+1}, energy {self.energies[homo_idx]}")
        print(f"LUMO: idx {lumo_idx+1}, energy {self.energies[lumo_idx]}")
    
        if homo_idx is not None:
            self.current_occ_idx = homo_idx
            self.display_om(homo_idx, occupied=True)
        if lumo_idx is not None:
            self.current_virt_idx = lumo_idx
            self.display_om(lumo_idx, occupied=False)


    def on_level_click(self, idx):
        occ = self.occupations[idx]
        occupied = occ > 0
        if occupied:
            self.current_occ_idx = idx
        else:
            self.current_virt_idx = idx

        self.display_om(idx, occupied)


# =============================================================================================================================================

class Node:
    def __init__(self, x, y, atom_type='C·'):
        self.x = x
        self.y = y
        self.atom_type = atom_type

# =============================================================================================================================================

class MoleculeDrawer:
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
        
        self.load_icons()
        self.create_toolbar()
        self.bind_shortcuts()

        self.canvas.bind("<Button-1>", self.left_click)
        self.canvas.bind("<Button-3>", self.right_click)
        self.canvas.bind("<Motion>", self.mouse_motion)
        self.canvas.bind("<B1-Motion>", self.mouse_drag)
        self.canvas.bind("<ButtonRelease-1>", self.mouse_release)
        self.canvas.bind("<Configure>", self.resize_canvas)

        self.om_window = None  # stocke la fenêtre DataFrame

        self.project_name = None
        self.safe_project_name = None

        self.df = None
        self.summary_data = None

    def sanitize_filename(self, name):
        return re.sub(r'[\\/*?:"<>|]', "_", name)
        
    def create_toolbar(self):
        self.btn_run = self.create_button(self.icons['run'], self.on_run_huckel, "Run")
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
        for name in ['run', 'matrix', 'save', 'load', 'undo', 'redo', 'eraser', 'clear', 'savedata', 'quit', 'about']:
            self.icons[name] = ImageTk.PhotoImage(Image.open(f"icons/{name}.png"))

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
        
    def toggle_eraser(self):
        self.eraser_mode = not self.eraser_mode
        self.btn_erase.config(bg='#ffb8b9' if self.eraser_mode else 'lightgray')

    def quit_program(self):
        self.master.quit()

    def show_about(self):
        messagebox.showinfo("About", "HMO Molecule Drawer\nVersion 0.9")

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
        step_x = int(GRID_SIZE * self.scale_x)
        step_y = int(GRID_SIZE * self.scale_y)
        for x in range(0, w, step_x):
            for y in range(0, h, step_y):
                self.canvas.create_oval(x-1, y-1, x+1, y+1, fill='gray', outline='gray', tags="grid")

    def snap_to_grid(self, x, y):
        return (round(x / (GRID_SIZE * self.scale_x)) * GRID_SIZE * self.scale_x,
                round(y / (GRID_SIZE * self.scale_y)) * GRID_SIZE * self.scale_y)

    def find_node(self, x, y):
        for idx, node in enumerate(self.nodes):
            node_x, node_y = self.apply_scale(node.x, node.y)
            if abs(node_x - x) < ATOM_RADIUS and abs(node_y - y) < ATOM_RADIUS:
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
            if abs(mid_x - x) < ATOM_RADIUS and abs(mid_y - y) < ATOM_RADIUS:
                return idx
        return None

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

    def left_click(self, event):
        x, y = self.snap_to_grid(event.x, event.y)
        x /= self.scale_x
        y /= self.scale_y
        if self.eraser_mode:
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
        idx = self.find_node(event.x, event.y)
        if idx is None:
            if not self.node_exists_at(x, y):
                self.save_state()
                self.nodes.append(Node(x, y))
                self.redraw()
        else:
            self.selected_node = idx
            self.dragging = True

    def right_click(self, event):
        idx = self.find_node(event.x, event.y)
        if idx is not None:
            menu = tk.Menu(self.master, tearoff=0)
            menu.configure(font=font.Font(family="Arial", size=12))
            bond_count = self.count_bonds(idx)
            for atom_type in ATOM_OPTIONS:
                if atom_type == 'Me' and bond_count > 1:
                    menu.add_command(label=atom_type + " (blocked)", state='disabled')
                else:
                    menu.add_command(label=atom_type,
                        command=lambda at=atom_type: self.change_atom_type(idx, at))
            try:
                menu.tk_popup(event.x_root, event.y_root)
            finally:
                menu.grab_release()

    def change_atom_type(self, idx, new_type):
        if new_type in ATOM_COLORS:
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

    def mouse_release(self, event):

        if self.dragging and self.selected_node is not None:
            x_snap, y_snap = self.snap_to_grid(event.x, event.y)
            x_grid = x_snap / self.scale_x
            y_grid = y_snap / self.scale_y
    
            # 1️⃣ Test précis sur grille
            idx_target = self.find_node_by_coords(x_grid, y_grid)
    
            if idx_target is None:
                # 2️⃣ Test tolérant (pixels) si l'utilisateur n'était pas parfaitement sur un noeud
                idx_target = self.find_node(event.x, event.y)
    
            if idx_target is None:
                # Aucun atome existant : crée un nouveau noeud + liaison
                self.save_state()
                self.nodes.append(Node(x_grid, y_grid))
                idx_target = len(self.nodes) - 1
    
            # Crée la liaison si elle n'existe pas encore
            if (self.selected_node, idx_target) not in self.bonds and \
               (idx_target, self.selected_node) not in self.bonds and \
               self.selected_node != idx_target:
                self.save_state()
                self.bonds.append((self.selected_node, idx_target))
    
            self.selected_node = None
            self.dragging = False
            self.redraw()
            # print(f"Release: target found by coords={idx_target is not None}")
    
    def count_bonds(self, idx):
        return sum(1 for i, j in self.bonds if i == idx or j == idx)

    def delete_node(self, idx):
        del self.nodes[idx]
        self.bonds = [(i, j) for i, j in self.bonds if i != idx and j != idx]
        self.bonds = [(i - (i > idx), j - (j > idx)) for i, j in self.bonds]

    def clear(self):
        self.nodes.clear()
        self.bonds.clear()
        self.undo_stack.clear()
        self.redo_stack.clear()
        self.df = None  # Optionnel si tu veux réinitialiser la dernière analyse Hückel
        self.redraw()

        if self.mo_viewer_window is not None and self.mo_viewer_window.winfo_exists():
            self.mo_viewer_window.destroy()
            self.mo_viewer_window = None

        if self.om_window is not None and self.om_window.winfo_exists():
            self.om_window.destroy()
            self.om_window = None

    def save_molecule(self):

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
            except Exception as e:
                messagebox.showerror("Error", f"Failed to load molecule: {e}")

    def redraw(self):
        self.canvas.delete("all")
        self.draw_grid()

        if self.dragging:
            snap_x, snap_y = self.snap_to_grid(*self.current_mouse_pos)
            self.canvas.create_oval(
                snap_x - HIGHLIGHT_RADIUS, snap_y - HIGHLIGHT_RADIUS,
                snap_x + HIGHLIGHT_RADIUS, snap_y + HIGHLIGHT_RADIUS,
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
            color = ATOM_COLORS.get(node.atom_type, 'black')
            if highlight_node == idx:
                self.canvas.create_oval(
                    x - ATOM_RADIUS - 4, y - ATOM_RADIUS - 4,
                    x + ATOM_RADIUS + 4, y + ATOM_RADIUS + 4,
                    outline='red', width=2
                )
            self.canvas.create_oval(
                x - ATOM_RADIUS, y - ATOM_RADIUS,
                x + ATOM_RADIUS, y + ATOM_RADIUS,
                fill=color
            )
            self.canvas.create_text(x, y, text=node.atom_type, fill='white', font=('Arial', 8, 'bold'))
    
    def evaluate(self,expr, alpha, beta):
        return eval(expr, {"alpha": alpha, "beta": beta})
    
    def build_huckel_matrix(self, alpha=-11.0, beta=-2.7):
        n = len(self.nodes)
        H = np.zeros((n, n))
        ### alpha (diagonal value)
        for i, node in enumerate(self.nodes):
            atom_type = node.atom_type
            try:
                param = Huckel_atomic_parameters[atom_type]
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
            H[i, i] = self.evaluate(param["alpha_expr"], alpha, beta)
        ### beta (non-diagonal value)
        for i, j in self.bonds:
            a = self.nodes[i].atom_type
            b = self.nodes[j].atom_type
        
            k_ij = Huckel_kXY_parameters.get(a, {}).get(b, None)
            if k_ij is None:
                k_ij = Huckel_kXY_parameters.get(b, {}).get(a, None)
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
            H[i, j] = H[j, i] = k_ij * beta
        return H

    def props(self, eigvals, occupation_dict, sorted_indices, alpha, beta):
        # Total π-electron energy
        total_energy = sum(eigvals[j] * occ for j, occ in occupation_dict.items())
        for j, occ in occupation_dict.items():
            print(j,eigvals[j], occ)
    
        # Partie alpha réelle (somme sur chaque atome : n_pi * alpha_effectif)
        alpha_part = 0
        alpha_atoms = 0
        print("alpha",alpha,beta)
        for i,n in enumerate(self.nodes):
            params = Huckel_atomic_parameters.get(n.atom_type)
            n_pi = params["n_pi"]
            try:
                alpha_atom = eval(params["alpha_expr"], {"alpha": alpha, "beta": beta})
            except Exception as e:
                print(f"Erreur pour {n.atom_type} : {e}")
                alpha_atom = alpha  # fallback sécurité

            alpha_atoms += n_pi * alpha_atom
            alpha_part += n_pi * alpha
            print(i,alpha_atom,n_pi, alpha_part, alpha_atoms)
    
        # Partie beta : ce qui reste
        beta_part = total_energy - alpha_part
    
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
        n_atoms = len(self.nodes)
        n_om = eigvecs.shape[1]
    
        # Charges : initialise à nombre π d'électrons
        charges = []
        for i, node in enumerate(self.nodes):
            n_pi = Huckel_atomic_parameters[node.atom_type]["n_pi"]
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
        
        def compute_occupations(eigvals, tol=1e-5):
            eigvals_sorted = np.sort(eigvals)
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
                print(f"Groupe indices: {group} ➔ Occ: {[occupation_dict[idx] for idx in group]}, Remaining: {remaining_electrons}")
            
                if remaining_electrons == 0:
                    break

            return occupation_dict

        alpha = -11.0
        beta = -2.7
        self.alpha_value = alpha
        self.beta_value = beta

        H = self.build_huckel_matrix(alpha, beta)
        if H is None:
            print("[INFO] Hückel analysis stopped due to a matrix-building error.")
            return
        eigvals, eigvecs = np.linalg.eigh(H)

        self.total_pi_electrons = 0
        for idx, n in enumerate(self.nodes):
            atom_type = n.atom_type
            try:
                param = Huckel_atomic_parameters[atom_type]
                n_pi = param["n_pi"]
                if n_pi is None:
                    raise ValueError(f"The 'n_pi' parameter is None for atom '{atom_type}'.")
                alpha_value = self.evaluate(param["alpha_expr"], alpha, beta)
                print(f"[INFO] Atom #{idx+1} ('{atom_type}'): n_pi = {n_pi}, alpha = {alpha_value:.2f}")  # ✅ log for each atom
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
        
        print(f"[INFO] Total π-electron count: {self.total_pi_electrons}")
        
        occupation_dict = compute_occupations(eigvals, self.total_pi_electrons)
        sorted_indices = np.argsort(eigvals)[::-1]

        self.charges, self.bond_orders = self.compute_charges_and_bond_orders(eigvecs, occupation_dict)

        self.props(eigvals, occupation_dict, sorted_indices, alpha, beta)

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
            energy_coeff = (eigvals[j] - alpha) / beta
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
        df = df.iloc[:, ::-1]
        self.df = df
        del df
        
    def build_dataframes(self):
        """Construit tous les DataFrames internes (sans sauvegarde)."""
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
                "Color": ATOM_COLORS.get(node.atom_type, 'gray')                
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
                f"{self.alpha_part/alpha:.0f}α + {self.beta_part/beta:.2f}β",
                f"{self.total_pi_electrons}",
                f"{self.atomization_energy / beta:.2f}",
                f"{self.atomization_energy_per_atom / beta:.2f}",
                f"{abs(self.homo_lumo_gap / beta):.2f}",
                f"{self.mu / beta:.2f}",
                f"{abs(self.eta / beta):.2f}",
                f"{abs(self.softness * beta):.2f}",
                f"{self.omega / (beta**2):.2f}",
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
        return True
    
        # Ouvre le fichier avec openpyxl pour modifier les styles
        wb = openpyxl.load_workbook(path)
        ws = wb['MO Coefficients']
    
        # Applique wrap_text=True pour la première ligne (en-têtes)
        for cell in ws[1]:
            cell.alignment = openpyxl.styles.Alignment(wrap_text=True)
    
        wb.save(path)

    def save_dataframe_as_pdf(self, win):
        from tkinter import filedialog, messagebox
        import matplotlib.pyplot as plt
        from pandas.plotting import table
        from pathlib import Path
    
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
            
            print(f"[DEBUG] avg_length: {avg_length:.2f} units")
            print(f"[DEBUG] Molecule bounds: x=({min_x:.1f}, {max_x:.1f}), y=({min_y:.1f}, {max_y:.1f})")
        
            canvas_width_px = fig_size_inches * dpi
        
            scale_limit = min(
                (canvas_width_px * 0.8) / width_molecule,
                (canvas_width_px * 0.8) / height_molecule
            ) if width_molecule > 0 and height_molecule > 0 else 1
        
            final_scale = scale_factor_base * min(scale_limit, 1)
        
            print(f"[DEBUG] Molecule raw size (after base scale): {width_molecule:.1f} x {height_molecule:.1f} px")
            print(f"[DEBUG] Scaling factor (liaison target ~1 cm): {scale_factor_base:.2f}")
            print(f"[DEBUG] Downscale factor (max 14x14 cm): {min(scale_limit, 1):.2f}")
            print(f"[DEBUG] Final scale: {final_scale:.2f}")
        
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
                color = ATOM_COLORS.get(node.atom_type, 'gray')
                ax.plot(x_scaled, y_scaled, 'o', color=color, markersize=20)
                label = f"{node.atom_type}{idx+1}"
                ax.text(x_scaled, y_scaled, label, ha='center', va='center', color='white', fontsize=8)
        
            ax.set_xlim(-canvas_width_px/2, canvas_width_px/2)
            ax.set_ylim(-canvas_width_px/2, canvas_width_px/2)


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
                0, -canvas_width_px/2 + 50,  # placé sous la molécule ; ajuste +50 si besoin
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

        print(f"[DEBUG] avg_length: {avg_length:.2f} px")
        print(f"[DEBUG] width_raw: {width_molecule_raw:.1f}, height_raw: {height_molecule_raw:.1f}")
        print(f"[DEBUG] scale_factor_base: {scale_factor_base:.2f}")
        print(f"[DEBUG] scale_limit: {scale_limit:.2f}")
        print(f"[DEBUG] final_scale: {final_scale:.2f}")
    
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
            color = ATOM_COLORS.get(node.atom_type, 'gray')
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

        summary_text = (
            f"Total π-electron energy: {self.alpha_part/alpha:.0f}α + {self.beta_part/beta:.2f}β\n"
            f"Number of π electrons: {self.total_pi_electrons}\n"
            f"Atomization energy: {self.atomization_energy / beta:.2f}β\n"
            f"Atomization energy per π atom: {self.atomization_energy_per_atom / beta:.2f}β\n"
            f"HOMO-LUMO gap: {abs(self.homo_lumo_gap / beta):.2f}|β|\n"
            f"Electronic potential μ: {self.mu / beta:.2f}β\n"
            f"Chemical hardness η: {abs(self.eta / beta):.2f}|β|\n"
            f"Chemical softness S: {abs(self.softness * beta):.2f}/|β|\n"
            f"Electrophilicity index ω: {self.omega / (beta**2):.2f}β\n"
        )
        self.summary_text = summary_text  # sauvegarde pour le PDF
        summary_label = tk.Label(frame_main, text=summary_text, font=("DejaVu Sans", 10), justify="left")
        summary_label.pack(pady=5)    
        
        # === 3️⃣ Table (Treeview + scrollbars) ===
        # frame_table = tk.Frame(frame_main)
        # frame_table.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
    
        # columns_with_index = ["Index"] + list(self.df.columns)
        # col_ids = [str(i) for i in range(len(columns_with_index))]
    
        # tree = ttk.Treeview(frame_table, columns=col_ids, show='headings')
    
        # # Header Index
        # tree.heading('0', text='Index')
        # tree.column('0', width=50, anchor='center')
    
        # # Autres colonnes
        # for i, col in enumerate(self.df.columns):
        #     tree.heading(str(i + 1), text=col)
        #     tree.column(str(i + 1), width=100, anchor='center')
    
        # for idx in self.df.index:
        #     row_values = [idx] + list(self.df.loc[idx])
        #     tree.insert('', 'end', values=row_values)
    
        # tree.grid(row=0, column=0, sticky='nsew')
        # scrollbar_y = ttk.Scrollbar(frame_table, orient=tk.VERTICAL, command=tree.yview)
        # scrollbar_y.grid(row=0, column=1, sticky='ns')
        # scrollbar_x = ttk.Scrollbar(frame_table, orient=tk.HORIZONTAL, command=tree.xview)
        # scrollbar_x.grid(row=1, column=0, sticky='ew')
    
        # frame_table.grid_rowconfigure(0, weight=1)
        # frame_table.grid_columnconfigure(0, weight=1)
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
            print("Escape key pressed. Closing the app.")
            win.destroy()
        win.bind('<Escape>', on_escape)
        
        btn_close = tk.Button(frame_buttons, text="Close", command=win.destroy)
        btn_close.pack()

    
    # Hook to existing GUI
    def on_run_huckel(self):
        self.run_huckel_analysis()
        from tkinter import simpledialog

        self.project_name = simpledialog.askstring(
            "Project name",
            "Enter a name for your project:",
            initialvalue="my_molecule"
        )
        
        if not self.project_name:
            self.project_name = "my_molecule"

        self.safe_project_name = self.sanitize_filename(self.project_name)

        self.build_dataframes()
        # self.show_dataframe_in_window()
        # === Ouvre la visualisation des OM ===
        print("[DEBUG] df_MOs shape:", self.df.shape)
        print("[DEBUG] df_atoms shape:", self.df_atoms.shape)
        print("[DEBUG] df_bonds shape:", self.df_bonds.shape)
        print("[DEBUG] df_summary shape:", self.df_summary.shape)
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
    def __init__(self, widget, text='widget info'):
        
        self.OpenSansReg_font = font.Font(family="Fonts/OpenSans/static/OpenSans-Regular.ttf")
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

