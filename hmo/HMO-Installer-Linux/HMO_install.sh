#!/bin/bash

set -e

echo "🔧 Starting installation of HMO..."

BIN_DEST="/usr/local/bin/HMO"
ICON_DEST="/usr/local/share/icons/hmo.png"
DESKTOP_DEST="/usr/local/share/applications/HMO.desktop"

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

echo "📂 Script directory detected: $SCRIPT_DIR"

# Copy the binary
echo "⚙️ Copying binary to $BIN_DEST..."
sudo install -Dm755 "$SCRIPT_DIR/HMO" "$BIN_DEST"

# Copy the icon
echo "🖼️ Copying icon to $ICON_DEST..."
sudo install -Dm644 "$SCRIPT_DIR/HMOicon.png" "$ICON_DEST"

# Copy and modify the .desktop file
echo "📝 Installing desktop shortcut to $DESKTOP_DEST..."
sed "s|Exec=.*|Exec=$BIN_DEST|; s|Icon=.*|Icon=$ICON_DEST|" "$SCRIPT_DIR/HMO.desktop" | sudo tee "$DESKTOP_DEST" > /dev/null

echo "✅ Installation completed successfully!"
echo "👉 You should now see HMO in your application menu."

