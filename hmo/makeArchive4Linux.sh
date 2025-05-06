#!/bin/bash
cp dist/HMO HMO-Installer-Linux/
cp icons/HMOicon.png HMO-Installer-Linux/
cp installer/* HMO-Installer-Linux/
tar zcvf HMO-Installer.tar.gz HMO-Installer-Linux/
