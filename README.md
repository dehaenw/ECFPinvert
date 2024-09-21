# ECFP invert
## what's this
This is the algorithmic ECFP inversion approach developed in our lab (https://github.com/lich-uct). It is fully refactored from the previous version and has superior performance now. If you would like to use the old version anyway, you can find it here: https://github.com/dehaenw/ECFPinvert/tree/da6301b511d076b341da67d50bb9edd7691307dc

ECFP is a type of circular fingerprint used to encode the presence of substructures in molecules. In this prototype implementation we use RDKit's MorganFingerprint as our ECFP. By a fragment per fragment build up of the molecule and checking it at every step, a good degree of reconstruction is possible. An alternative NN based method for inversion was developed earlier and can be found here: https://github.com/bayer-science-for-a-better-life/neuraldecipher Their paper also contains some good background on why it is interesting to show some ECFP can be inverted.

## installation
### PIP environment setup
requires PIP within python 3.11 or lower. Then just:
```bash
conda create -n [environment_name] python=3.11  # or any other way to get python env 
git clone https://github.com/dehaenw/ECFPinvert
cd ECFPinvert
pip install .
```
### manual environment setup

Alternatively, just make sure the conda env or venv has rdkit=2022.09 and numpy=1.23.4.
Yes, it is important to have these exact versions (for now).

## i want to try it
here is a minimally working code block to invert the ECFP4(2048) of strychnine, a rather complex molecule:
```python
from rdkit import Chem
from search import ECFPInvert
import utils
utils.initialize_atomtypes("CHEMBL")
utils.set_fp_settings(2,2048)
strychnine = Chem.MolFromSmiles("C1CN2CC3=CCO[C@H]4CC(=O)N5[C@H]6[C@H]4[C@H]3C[C@H]2[C@@]61C7=CC=CC=C75")
strychnine_fp = utils.get_fp(strychnine)
inv = ECFPInvert()
s, info = inv.run_search(strychnine_fp)
if s:
    print(f'inverted strychnine in {info["time"]} seconds')
else:
    print("failure")

```
if that worked it should print something like `inverted strychnine in 9.421 seconds`

I will add example notebooks in a few days. Stay tuned!


## what are the limitations and plans
Exact reconstruct only. Can't be used to reconstruct non existent fingerprints that do not correspond to a real structure.
This version had ca 97% success inverting ChEMBL like ECFP6(4096).
