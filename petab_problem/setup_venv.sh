#!/usr/bin/env bash
python3 -m venv ./venv

source ./venv/bin/activate
export AMICI_CXXFLAGS="-Xpreprocessor -fopenmp"
export AMICI_LDFLAGS="-Xpreprocessor -fopenmp -lomp"
pip3 install -e git+https://github.com/ICB-DCM/AMICI@develop#egg=amici\&subdirectory=python/sdist
pip3 install -e git+https://github.com/ICB-DCM/pyPESTO@feature_207_petab#egg=pypesto
pip3 install -e git+https://github.com/PEtab-dev/PEtab@develop#egg=petab
pip3 install pysb