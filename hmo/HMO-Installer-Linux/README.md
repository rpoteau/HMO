# 🚀 **HMO Installer Guide**

**HMO** is a desktop application for Hückel Molecular Orbital diagrams.

---

## 🔧 How to Install

1. **Download** the `HMO-Installer.tar.gz` archive.

2. **Extract it:**

   ```bash
   tar xzvf HMO-Installer.tar.gz
   cd HMO-Installer
   ```

3. **Make sure the installer script is executable:**

   ```bash
   chmod +x HMO_install.sh
   ```

4. **Run the installer as root:**

   ```bash
   sudo ./HMO_install.sh
   ```

---

✅ **After installation:**

- The HMO application will be installed to:

  `/usr/local/bin/HMO`

- The icon will be installed to:

  `/usr/local/share/icons/hmo.png`

- A desktop shortcut will appear in your application menu (under **Science** or **Education** categories).

You can launch HMO from your application menu or by running:

```bash
HMO
```

from a terminal.

---

## 🔄 How to Uninstall

1. **Make the uninstaller script executable (if needed):**

   ```bash
   chmod +x HMO_uninstall.sh
   ```

2. **Run the uninstaller as root:**

   ```bash
   sudo ./HMO_uninstall.sh
   ```

This will remove:

- The binary (`/usr/local/bin/HMO`)
- The icon (`/usr/local/share/icons/hmo.png`)
- The desktop shortcut (`/usr/share/applications/HMO.desktop`)

---

## ❗ Notes for Linux Users

- If the icon or desktop entry doesn’t appear immediately, try logging out and logging back in, or run:

  ```bash
  sudo update-desktop-database /usr/share/applications
  ```

- The installation requires `sudo` because it writes to system directories.

---

# 👤 Author

This application was developed by **Romuald Poteau**.

---

# 📄 License

This software is distributed under the terms of the **GNU General Public License v3.0 (GPL-3.0)**.  
See the `LICENSE` file included in this package for full details.

---

# ✅ That’s it – enjoy using **HMO**! 🎉