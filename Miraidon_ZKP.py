import os
import hashlib
import numpy as np
import galois

# ==========================================
# Primitives
# ==========================================

def commit(*args):
    """SHA-256 commmitments"""
    h = hashlib.sha256()
    for arg in args:
        if isinstance(arg, galois.FieldArray):
            h.update(arg.tobytes())
        elif isinstance(arg, bytes):
            h.update(arg)
        elif isinstance(arg, int):
            h.update(arg.to_bytes(4, 'big'))
    return h.digest()

def expand_seed(seed, m, n, GF):
    """
    Expands a seed deterministically into two invertible matrices S and T
    Reject if not in GL_m and GL_n
    """
    attempt = 0
    S, T = None, None
    
    while S is None or T is None:
        state = int.from_bytes(hashlib.sha256(seed + attempt.to_bytes(4, 'big')).digest()[:8], 'big')
        rng = np.random.default_rng(state)
        
        if S is None:
            candidate_S = GF(rng.integers(0, GF.order, size=(m, m)))
            if np.linalg.det(candidate_S) != 0:
                S = candidate_S
                
        if T is None:
            candidate_T = GF(rng.integers(0, GF.order, size=(n, n)))
            if np.linalg.det(candidate_T) != 0:
                T = candidate_T
                
        attempt += 1
        
    return S, T






# ==========================================
# Miraidon Prover & Verifier
# ==========================================

class MiraidonProver:
    def __init__(self, GF, M, a, L_base, R_base, lambda_sec=16):
        self.GF = GF
        self.M = M
        self.k = len(M)
        self.m, self.n = M[0].shape
        self.a = a
        self.L_base = L_base
        self.R_base = R_base
        self.lambda_sec = lambda_sec
        self.state = {} # Holds internal round state

    def round_1(self):
        """Generate initial commitments."""
        # Pick seed_{ST}
        self.state['seed_ST'] = os.urandom(2 * self.lambda_sec)
        
        # Expand seed to S, T
        S, T = expand_seed(self.state['seed_ST'], self.m, self.n, self.GF)
        self.state['S'] = S
        self.state['T'] = T
        
        # Select random u
        u = self.GF.Random(self.k)
        self.state['u'] = u
        
        # Set X = sum(u_i * (S * M_i * T))
        X = self.GF.Zeros((self.m, self.n))
        for i in range(self.k):
            X += u[i] * (S @ self.M[i] @ T)
        self.state['X'] = X
        
        # Solve Q = LR
        # Since E = L_base * R_base, we have Q = S*E*T = (S*L_base) * (R_base*T)
        L = S @ self.L_base
        R = self.R_base @ T
        Q = L @ R
        self.state['L'] = L
        self.state['R'] = R
        self.state['Q'] = Q
        
        # Commitments
        c0 = commit(X, Q)
        c1 = commit(S, T)
        
        return c0, c1

    def round_2(self, xi):
        """Prover receives challenge xi and sends response y."""
        self.state['xi'] = xi
        u = self.state['u']
        S = self.state['S']
        T = self.state['T']
        
        # Define v = u + xi * a
        v = u + xi * self.a
        self.state['v'] = v
        
        # Set y = commit( sum(v_i * (S * M_i * T)) )
        val = self.GF.Zeros((self.m, self.n))
        for i in range(self.k):
            val += v[i] * (S @ self.M[i] @ T)
            
        y = commit(val)
        return y

    def round_3(self, b):
        """Prover receives challenge b and reveals the corresponding secrets."""
        if b == 0:
            return {"X": self.state['X'], "L": self.state['L'], "R": self.state['R']}
        elif b == 1:
            return {"seed_ST": self.state['seed_ST'], "v": self.state['v']}
        else:
            raise ValueError("b must be 0 or 1")

class MiraidonVerifier:
    def __init__(self, GF, M):
        self.GF = GF
        self.M = M
        self.k = len(M)
        self.m, self.n = M[0].shape
        self.state = {}

    def receive_round_1(self, c0, c1):
        """Verifier stores c0, c1 and returns challenge xi."""
        self.state['c0'] = c0
        self.state['c1'] = c1
        
        # Sample xi from Fq* 
        xi = self.GF(0)
        while xi == 0:
            xi = self.GF.Random()
            
        self.state['xi'] = xi
        return xi

    def receive_round_2(self, y):
        """Verifier stores y and returns challenge b."""
        self.state['y'] = y
        
        # Pick b randomly from {0, 1}
        b = int(os.urandom(1)[0] % 2)
        self.state['b'] = b
        return b

    def verify(self, resp):
        """Verifier checks the final response based on the challenge b."""
        b = self.state['b']
        c0 = self.state['c0']
        c1 = self.state['c1']
        xi = self.state['xi']
        y = self.state['y']
        
        if b == 0:
            X, L, R = resp['X'], resp['L'], resp['R']
            
            # Reconstruct Q
            Q = L @ R
            
            # Check c0
            if c0 != commit(X, Q):
                print("Reject: c0 commitment mismatch.")
                return False
                
            # Check y = commit(X + xi * Q)
            if y != commit(X + xi * Q):
                print("Reject: y commitment mismatch (b=0).")
                return False
                
            return True
            
        elif b == 1:
            seed_ST, v = resp['seed_ST'], resp['v']
            
            # Expand seed
            S, T = expand_seed(seed_ST, self.m, self.n, self.GF)
            
            # Check c1
            if c1 != commit(S, T):
                print("Reject: c1 commitment mismatch.")
                return False
                
            # Check y = commit( sum(v_i * (S * M_i * T)) )
            val = self.GF.Zeros((self.m, self.n))
            for i in range(self.k):
                val += v[i] * (S @ self.M[i] @ T)
                
            if y != commit(val):
                print("Reject: y commitment mismatch (b=1).")
                return False
                
            return True







# ==========================================
# Execution
# ==========================================
if __name__ == "__main__":
    # 1. Setup params
    q = 11
    m, n, k, r = 10, 10, 5, 3
    GF = galois.GF(q)
    
    print(f"Constructing MinRank instance (q={q}, m={m}, n={n}, k={k}, r={r})")
    
    # 2. Generate keys
    a = GF.Random(k)
    a[-1] = 1
    
    L_base = GF.Random((m, r))
    R_base = GF.Random((r, n))
    E = L_base @ R_base
    
    M = [GF.Random((m, n)) for _ in range(k - 1)]
    
    M_k = E.copy()
    for i in range(k - 1):
        M_k -= a[i] * M[i]
    M.append(M_k)

    # 3. Instantiate ZKP
    prover = MiraidonProver(GF, M, a, L_base, R_base)
    verifier = MiraidonVerifier(GF, M)


    
    # Round 1: Commitments
    c0, c1 = prover.round_1()
    print(f"Prover -> Verifier: c0 = {c0.hex()}")
    print(f"Prover -> Verifier: c1 = {c1.hex()}\n")
    # print("Prover  -> Verifier: ({c0.hex()}, {c1.hex()})")
    
    # Round 2: Challenge 1
    xi = verifier.receive_round_1(c0, c1)
    print(f"Verifier -> Prover: xi = {xi}")
    
    # Round 3: Response 1
    y = prover.round_2(xi)
    # print("Prover  -> Verifier: {y.hex()}}")
    print(f"Prover -> Verifier: y = {y.hex()}\n")
    
    # Round 4: Challenge 2
    b = verifier.receive_round_2(y)
    print(f"Verifier -> Prover: b = {b}")
    
    # Round 5: Response 2
    resp = prover.round_3(b)
    # keys = ", ".join(resp.keys())
    # print(f"Prover  -> Verifier: resp = {{{keys}}}")


    # Format and print the dictionary contents dynamically based on 'b'
    print("Prover -> Verifier: resp = {")
    for key, value in resp.items():
        if isinstance(value, bytes):
            print(f"    '{key}': {value.hex()}")
        else:
            array_str = str(value).replace('\n', '\n          ')
            print(f"    '{key}': {array_str}")
    print("}\n")
    

    if b == 0:
        X, L, R = resp['X'], resp['L'], resp['R']
        
        # Manually reconstruct
        Q_calc = L @ R
        c0_calc = commit(X, Q_calc)
        y_calc = commit(X + xi * Q_calc)
        
        print(f"Recomputed c0 : {c0_calc.hex()}")
        print(f"Recomputed y : {y_calc.hex()}")

    elif b == 1:
        seed_ST, v = resp['seed_ST'], resp['v']
        
        # Manually reconstruct
        S_calc, T_calc = expand_seed(seed_ST, m, n, GF)
        c1_calc = commit(S_calc, T_calc)
        
        val = GF.Zeros((m, n))
        for i in range(k):
            val += v[i] * (S_calc @ M[i] @ T_calc)
        y_calc = commit(val)

        print(f"Recomputed c1 : {c1_calc.hex()}")
        print(f"Recomputed y : {y_calc.hex()}")

    # Verification
    is_valid = verifier.verify(resp)
    print(f"\n Signature {'is accepted!' if is_valid else 'IS NOT VALID'}")
