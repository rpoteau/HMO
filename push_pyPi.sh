#!/bin/bash

# Couleurs ANSI
GREEN='\033[1;32m'
YELLOW='\033[1;33m'
CYAN='\033[1;36m'
RED='\033[1;31m'
WHITE_BG_BLACK_TEXT='\033[30;47m'
RESET='\033[0m'

PYPROJECT="pyproject.toml"
DIST_DIR="dist"

# Ligne séparatrice claire
SEPARATOR="${WHITE_BG_BLACK_TEXT}--------------------------------------------------------------------------------${RESET}"

# Lire la version actuelle dans pyproject.toml
CURRENT_VERSION=$(grep "^version" $PYPROJECT | head -n1 | cut -d '"' -f2)
echo -e "$SEPARATOR"
echo -e "${WHITE_BG_BLACK_TEXT}  Version actuelle dans pyproject.toml : $CURRENT_VERSION  ${RESET}"

# Vérifier s'il y a un .tar.gz existant dans dist/
if [ -d "$DIST_DIR" ]; then
    echo -e "$SEPARATOR"
    echo -e "${WHITE_BG_BLACK_TEXT}  Archives trouvées dans dist/:  ${RESET}"
    ls dist/*.tar.gz 2>/dev/null || echo -e "${YELLOW}Aucune archive tar.gz trouvée.${RESET}"
else
    echo -e "$SEPARATOR"
    echo -e "${WHITE_BG_BLACK_TEXT}  Pas de répertoire dist/.  ${RESET}"
fi

# Récupérer la dernière version publiée sur PyPI (optionnel)
PACKAGE_NAME=$(grep "^name" $PYPROJECT | head -n1 | cut -d '"' -f2)
echo -e "$SEPARATOR"
echo -e "${WHITE_BG_BLACK_TEXT}  Interrogation de PyPI pour $PACKAGE_NAME...  ${RESET}"
LATEST_PYPI=$(curl -s https://pypi.org/pypi/$PACKAGE_NAME/json | jq -r '.info.version')

if [ "$LATEST_PYPI" != "null" ]; then
    echo -e "${CYAN}Dernière version publiée sur PyPI :${RESET} ${YELLOW}$LATEST_PYPI${RESET}"
else
    echo -e "${RED}Le package n'existe pas sur PyPI (ou erreur PyPI).${RESET}"
fi

# Demander si on veut incrémenter
echo -e "$SEPARATOR"
echo -e "${WHITE_BG_BLACK_TEXT}  Souhaitez-vous incrémenter la version ? (y/n)  ${RESET}"
read -r REPLY

if [[ "$REPLY" =~ ^[Yy]$ ]]; then
    echo -e "${CYAN}Quel niveau ? (patch / minor / major)${RESET}"
    read -r LEVEL

    IFS='.' read -r MAJOR MINOR PATCH <<< "$CURRENT_VERSION"

    case $LEVEL in
        patch)
            PATCH=$((PATCH + 1))
            ;;
        minor)
            MINOR=$((MINOR + 1))
            PATCH=0
            ;;
        major)
            MAJOR=$((MAJOR + 1))
            MINOR=0
            PATCH=0
            ;;
        *)
            echo -e "${RED}Type inconnu, version non modifiée.${RESET}"
            ;;
    esac

    NEW_VERSION="$MAJOR.$MINOR.$PATCH"
    echo -e "${GREEN}Mise à jour vers : $NEW_VERSION${RESET}"

    # Modifier le pyproject.toml en place
    sed -i "s/^version = \".*\"/version = \"$NEW_VERSION\"/" $PYPROJECT
else
    echo -e "${YELLOW}Version conservée.${RESET}"
fi

# Nettoyer les anciens builds
echo -e "$SEPARATOR"
echo -e "${WHITE_BG_BLACK_TEXT}  Suppression des anciens builds : rm -rf build dist *.egg-info  ${RESET}"
rm -rf build dist *.egg-info

# Build du package
echo -e "$SEPARATOR"
echo -e "${WHITE_BG_BLACK_TEXT}  Construction du package : python -m build  ${RESET}"
python -m build

# Upload vers PyPI
echo -e "$SEPARATOR"
echo -e "${WHITE_BG_BLACK_TEXT}  Upload vers PyPI : twine upload dist/*  ${RESET}"
twine upload dist/*

# Réinstallation en mode editable
echo -e "$SEPARATOR"
echo -e "${WHITE_BG_BLACK_TEXT}  Réinstallation en mode editable : pip install -e .  ${RESET}"
pip install -e .

echo -e "$SEPARATOR"
echo -e "${GREEN}🎉 Processus terminé !${RESET}"

