    """
    Computes minimal secure w with the optimal values of t_star and alpha
    
    Parameters:
    lam            : Security parameter (128, 192, 256)
    q              : Field size (2, 127)
    t              : Number of rounds repeated under Fiat-Shamir
    w_start        : Starts here then iterates down, checking forgery cost
    w_step         : Change in each iterate (so w_step = 5 and w_start = 100 checks 100, then 95, then 90, ...)
    t_star_min     : Minimal value of t_star to check
    t_star_max     : Maxmal value of t_star to check
    """


import math

def log_choose(n, k):
    if k < 0 or k > n:
        return -float('inf')
    return math.lgamma(float(n + 1)) - math.lgamma(float(k + 1)) - math.lgamma(float(n - k + 1))

def log_sum_exp(log_vals):
    valid_vals = [v for v in log_vals if v != -float('inf')]
    if not valid_vals:
        return -float('inf')
    max_v = max(valid_vals)
    return max_v + math.log(sum(math.exp(float(v - max_v)) for v in valid_vals))

def find_minimal_secure_w_fast(lam, q, t, w_start, w_step, t_star_min, t_star_max):
    lam, q, t = int(lam), int(q), int(t)
    w_start, w_step = int(w_start), int(w_step)
    
        
    print(f"Minimal secure w (t={t}, q={q}, target={lam}-bit)...")
    print(f"Restricting t* search space to [{t_star_min}, {t_star_max}]")
    
    log_B = []
    log_p_succ = math.log(float(1) / float(q - 1))
    log_p_fail = math.log(float(1) - float(1) / float(q - 1))
    
    for j in range(t + 1):
        val = log_choose(t, j) + float(j) * log_p_succ + float(t - j) * log_p_fail
        log_B.append(val)

    for w in range(w_start, 0, -w_step):
        log_S = {}
        log_denom_w = log_choose(t, w)
        
        for alpha in range(w, t + 1):
            log_S[alpha] = []
            log_denom_alpha = log_choose(t, alpha)
            log_denom_total = log_denom_w + log_denom_alpha
            
            for j in range(t + 1):
                lower_bound = max(0, alpha - j)
                upper_bound = min(t - j, w)
                
                log_terms = []
                for w_star in range(lower_bound, upper_bound + 1):
                    if alpha - w_star >= 0 and w - w_star >= 0:
                        term = log_choose(t - j, w_star) + log_choose(j, alpha - w_star) + log_choose(j, w - w_star)
                        log_terms.append(term)
                
                if not log_terms:
                    log_S[alpha].append(-float('inf'))
                else:
                    log_s_val = log_sum_exp(log_terms) - log_denom_total
                    log_S[alpha].append(log_s_val)

        min_obj_log = float('inf')
        best_t_star = None
        best_alpha = None
        
        for t_star in range(t_star_min, t_star_max + 1):
            log_P_beta = log_sum_exp(log_B[t_star:])
            if log_P_beta == -float('inf'):
                continue
                
            max_log_Pb = -float('inf')
            best_alpha_for_t_star = None
            
            for alpha in range(w, t + 1):
                log_Nb_terms = [log_B[j] + log_S[alpha][j] for j in range(t_star, t + 1)]
                log_Nb = log_sum_exp(log_Nb_terms)
                
                log_Pb_alpha = log_Nb - log_P_beta
                
                if log_Pb_alpha > max_log_Pb:
                    max_log_Pb = log_Pb_alpha
                    best_alpha_for_t_star = alpha
                    
            if max_log_Pb > -float('inf'):
                current_obj_log = log_sum_exp([-log_P_beta, -max_log_Pb])
                
                if current_obj_log < min_obj_log:
                    min_obj_log = current_obj_log
                    best_t_star = t_star
                    best_alpha = best_alpha_for_t_star
                    
        if min_obj_log != float('inf'):
            sec_bits = min_obj_log / math.log(float(2))
            
            if w % w_step == w_start % w_step:
                print(f"Tested w={w:3d} | Sec: {sec_bits:6.2f} bits | Best t*: {best_t_star}")
            
            if sec_bits >= lam:
                print("\n" + "="*45)
                print("TARGET SECURITY REACHED")
                print("="*45)
                print(f"First Secure w : {w}")
                print(f"Optimal t* : {best_t_star}")
                print(f"Optimal alpha  : {best_alpha}")
                print(f"Security Level : {sec_bits:.2f} bits")
                print("="*45 + "\n")
                
                return w, best_t_star, best_alpha, sec_bits

    print("Target security not reached within the bounds.")
    return None

# ==========================================

find_minimal_secure_w_fast(lam=256, q=509, t=512, w_start=436, w_step=1, t_star_min=50, t_star_max=60)
