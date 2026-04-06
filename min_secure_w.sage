import math
from decimal import Decimal, getcontext

# Set high precision 
getcontext().prec = int(150)



def find_minimal_secure_w(lam, q, t, w_start, w_step):
    """
    Computes the minimal weight w to achieve the given security level against advanced forgery attack
    Starts from w_start, then iterates down in steps of size w_step
    
    Parameters:
    lam     : Security parameter (128, 192, 256)
    q       : Field size (2, 127)
    t       : Number of rounds repeated under Fiat-Shamir
    w_start : Number of b=1 challenges; code starts here (guess high)
    w_step  : Iterates down from w_start this many steps (pick 5 for first run)
    """
    
    print(f"Searching for minimum secure w (t={t}, q={q}, target={lam}-bit)...")
    
    # Precompute B(j) 
    B = []
    p_succ = Decimal(int(1)) / Decimal(int(q - 1))
    p_fail = Decimal(int(1)) - p_succ
    for j in range(t + 1):
        B.append(Decimal(math.comb(t, j)) * (p_succ ** j) * (p_fail ** (t - j)))

    # Iterate w and break early once lambda is reached
    for w in range(w_start, 0, -w_step):
        # Precompute S[alpha][j] for the current w
        S = {}
        denom_w = Decimal(math.comb(t, w))
        
        for alpha in range(w, t + 1):
            S[alpha] = []
            denom_alpha = Decimal(math.comb(t, alpha))
            denom_total = denom_w * denom_alpha
            
            for j in range(t + 1):
                s_val = Decimal(int(0))
                lower_bound = max(0, alpha - j)
                upper_bound = min(t - j, w)
                
                for w_star in range(lower_bound, upper_bound + 1):
                    if alpha - w_star >= int(0) and w - w_star >= int(0):
                        num = math.comb(t - j, w_star) * math.comb(j, alpha - w_star) * math.comb(j, w - w_star)
                        s_val += Decimal(num)
                        
                S[alpha].append(s_val / denom_total)

        # Optimize over t* and alpha
        min_obj_for_w = None
        best_t_star_for_w = None
        best_alpha_for_w = None
        
        for t_star in range(t + 1):
            P_beta = sum(B[t_star:])
            if P_beta == Decimal(int(0)):
                continue
                
            max_Pb_num = Decimal(int(0))
            best_alpha_for_t_star = None
            
            for alpha in range(w, t + 1):
                N_b = Decimal(int(0))
                for j in range(t_star, t + 1):
                    N_b += B[j] * S[alpha][j]
                    
                if N_b > max_Pb_num:
                    max_Pb_num = N_b
                    best_alpha_for_t_star = alpha
                    
            if max_Pb_num > Decimal(int(0)):
                current_obj = (Decimal(int(1)) / P_beta) + (P_beta / max_Pb_num)
                
                if min_obj_for_w is None or current_obj < min_obj_for_w:
                    min_obj_for_w = current_obj
                    best_t_star_for_w = t_star
                    best_alpha_for_w = best_alpha_for_t_star
                    
        # Check if this w meets the security target
        if min_obj_for_w is not None:
            sec_bits = float(min_obj_for_w.ln() / Decimal(int(2)).ln())
            
            if sec_bits >= lam:
                print("\n" + "="*45)
                print("TARGET SECURITY REACHED")
                print("="*45)
                print(f"First Secure w : {w}")
                print(f"Optimal t* : {best_t_star_for_w}")
                print(f"Optimal alpha  : {best_alpha_for_w}")
                print(f"Security Level : {sec_bits:.2f} bits")
                print("="*45 + "\n")
                
                return # w, best_t_star_for_w, best_alpha_for_w, sec_bits
                
        # Print occasional progress report 
        if w % w_step == w_start % w_step:
            print(f"Tested up to w={w}... (Current Security: {sec_bits:.2f} bits)")

    print("Target security not reached within the bounds of t.")
    return None



find_minimal_secure_w(lam=128, q=127, t=179, w_start=125, w_step=5)
