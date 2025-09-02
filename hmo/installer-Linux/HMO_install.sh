#!/bin/bash

set -e

echo "Installing HMO..."

BIN_DEST="/usr/local/bin/HMO"
ICON_DEST="/usr/local/share/icons/hmo.png"
DESKTOP_DEST="/usr/local/share/applications/HMO.desktop"

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

# Copier le binaire
sudo install -Dm755 "$SCRIPT_DIR/HMO" "$BIN_DEST"

# Copier l'icône
sudo install -Dm644 "$SCRIPT_DIR/HMOicon.png" "$ICON_DEST"

# Copier/modifier le .desktop
sed "s|Exec=.*|Exec=$BIN_DEST|; s|Icon=.*|Icon=$ICON_DEST|" "$SCRIPT_DIR/HMO.desktop" | sudo tee "$DESKTOP_DEST" > /dev/null

echo "Installation completed. HMO should now appear in your app menu."

