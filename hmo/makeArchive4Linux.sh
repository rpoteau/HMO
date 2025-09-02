#!/bin/bash
cp dist/HMO HMO-Installer-Linux/
cp icons/HMOicon.png HMO-Installer-Linux/
cp installer-Linux/* HMO-Installer-Linux/
tar zcvf HMO-Installer-Linux.tar.gz HMO-Installer-Linux/
