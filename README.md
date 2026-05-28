# Miraidon

This repository holds the scripts involved in [Miraidon](https://eprint.iacr.org/2026/997) by Freeman Slaughter and Ryann Cartor, a novel MinRank identification scheme used to design a signature, ring signature, and linkable ring signature scheme.

- [Miraidon_ZKP.py](./Miraidon_ZKP.py) holds the proof-of-concept for the underlying zero-knowledge proof
- [min_secure_w.py](./min_secure_w.py) is the script which determines the minimal secure weight of the challenge vector $b$, with regards to forgery attacks
- [sig_pubkey_size.py](./sig_pubkey_size.py) computes the signature (and ring signature, and linkable ring signature) and public key sizes for a user-input parameter set

The above scripts are in Python, while the latter two are acompanied by sage script accomplishing the same task. Please use at your own risk.
