import math

def sig_and_pubkey_size(lam, q, m, n, r, k, t, w):
    """
    Computes the signature and public key size for a given security level
    
    Parameters:
    lam   : Security parameter (128, 192, 256)
    m, n  : Matrix dimensions (rows/cols)
    r     : Rank of secret matrix
    q     : Field size (2, 127)
    k     : Number of matrices in MinRank instance
    t     : Number of rounds repeated under Fiat-Shamir
    w     : Number of rounds with challenge b=1
    """
    
    # Precompute log_2(q) and weight of t when written in binary
    log2_q = math.log2(q)
    wt_t = bin(t).count('1')
            
    # Objective function components
    c0_c1_y_salt = 8 * lam
            
    # SeedTree + MerkleProof cost
    tree_inner = (t - w) * math.log2(t / (t - w)) + wt_t - 1
    seed_tree = 3 * lam * math.floor(tree_inner)
            
    # b=0 challenge (i not in I)
    resp_notin_I = (t - w) * (2 * lam + (m*n + (m + n)*r) * log2_q)
            
    # b=1 challenge (i in I)
    resp_in_I = w * (k * log2_q)
            
    # Total Size
    sig_size = c0_c1_y_salt + seed_tree + resp_notin_I + resp_in_I
    pubkey_size = lam + (m*n - k)*log2_q

    print(f"Signature Size = {sig_size / 8000:.2f} kB, Public Key Size = {pubkey_size / 8:.2f} B")


def ring_signature_size(lam, q, m, n, r, k, t, w, ell):
    """
    Computes the signature size for a ring signature
    Note that functionally, this only changes the length of v, so just the b=1 response is different
    
    Parameters:
    lam   : Security parameter (128, 192, 256)
    m, n  : Matrix dimensions (rows/cols)
    r     : Rank of secret matrix
    q     : Field size (2, 127)
    k     : Number of matrices in MinRank instance
    t     : Number of rounds repeated under Fiat-Shamir
    w     : Number of rounds with challenge b=1
    ell   : Number of members in the ring signature
    """
    
    # Precompute log_2(q) and weight of t when written in binary
    log2_q = math.log2(q)
    wt_t = bin(t).count('1')
                
    # Objective function components
    c0_c1_y_salt = 8 * lam
            
    # SeedTree + MerkleProof cost
    tree_inner = (t - w) * math.log2(t / (t - w)) + wt_t - 1
    seed_tree = 3 * lam * math.floor(tree_inner)
            
    # b=0 (i not in I)
    resp_notin_I = (t - w) * (2 * lam + (m*n + (m + n)*r) * log2_q)
            
    # b=1 (i in I)
    resp_in_I = w * ((k + ell) * log2_q)
            
    # Total Size
    sig_size = c0_c1_y_salt + seed_tree + resp_notin_I + resp_in_I

    print(f"Ring Signature Size = {sig_size / 8000:.2f} kB")


# ==========================================

sig_and_pubkey_size(lam=256, q=509, m=8, n=16, r=8, k=30, t=512, w=435)

# Note: Python uses ** for exponentiation instead of ^
# ring_signature_size(lam=128, q=127, m=19, n=20, r=2, k=60, t=176, w=115, ell=2**3)
