import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.mixture import GaussianMixture
from scipy.stats import norm

def main():
    print("--- Generating BDP-FS v2 (GMM-Soft) Masterpiece Visualization ---")
    
    # 1. High-Fidelity Simulation (Methodological Clarity)
    np.random.seed(42)
    n_genes = 10000
    
    # Native Component (Genes with baseline technical variance)
    native = np.random.normal(loc=0.8, scale=0.3, size=int(n_genes * 0.8))
    # Distortion/Noise Component (Genes with extreme batch effects)
    noise = np.random.normal(loc=2.5, scale=0.8, size=int(n_genes * 0.2))
    
    d_scores = np.concatenate([native, noise])
    d_scores = d_scores[d_scores > 0] # Non-negative
    
    # 2. Fit GMM
    gmm = GaussianMixture(n_components=2, random_state=42)
    gmm.fit(d_scores.reshape(-1, 1))
    
    # Identify Native Component (lower mean)
    native_idx = np.argmin(gmm.means_.flatten())
    mu_native = gmm.means_.flatten()[native_idx]
    sigma_native = np.sqrt(gmm.covariances_.flatten()[native_idx])
    
    # Calculate tau_0 (95th percentile of Native)
    tau_0 = mu_native + 1.645 * sigma_native
    
    # 3. Compute Soft Weights (alpha=1.0)
    alpha = 1.0
    weights = np.exp(-alpha * np.maximum(0, d_scores - tau_0))
    
    # 4. Multi-Panel Plot
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10), sharex=True, gridspec_kw={'height_ratios': [2, 1]})
    sns.set_style("whitegrid")
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Inter', 'Roboto', 'Arial']

    # PANEL A: GMM DENSITY & ANCHOR
    x_range = np.linspace(0, 5, 1000).reshape(-1, 1)
    responsibilities = gmm.predict_proba(x_range)
    log_weights = gmm.score_samples(x_range)
    pdf = np.exp(log_weights)
    
    # Hist
    ax1.hist(d_scores, bins=100, density=True, alpha=0.15, color='gray', label='Observed Distortion Scores ($D_g$)')
    
    # GMM Components
    ax1.plot(x_range, pdf, color='black', lw=2.5, label='BDP-FS GMM Hybrid Fit')
    ax1.plot(x_range, responsibilities[:, native_idx] * norm.pdf(x_range, mu_native, sigma_native).flatten() * gmm.weights_[native_idx], 
             color='#1f77b4', lw=2, linestyle='--', label='Component 1: Native')
    ax1.plot(x_range, responsibilities[:, 1-native_idx] * norm.pdf(x_range, gmm.means_.flatten()[1-native_idx], np.sqrt(gmm.covariances_.flatten()[1-native_idx])).flatten() * gmm.weights_[1-native_idx], 
             color='#7f7f7f', lw=2, linestyle='--', label='Component 2: Technical Noise')

    # Anchor Marker
    ax1.axvline(x=tau_0, color='#d62728', lw=2, linestyle='-', label=f'τ_0 Anchor (95th Pctl Native)')
    ax1.fill_between(x_range.flatten(), 0, pdf, where=(x_range.flatten() > tau_0), color='#d62728', alpha=0.1, label='Penalization Zone')

    ax1.set_ylabel("Density", fontsize=12)
    ax1.set_title("BDP-FS v2 'Masterpiece' Logic: GMM-Anchored Soft Weighting", fontsize=15, fontweight='bold', pad=20)
    ax1.legend(frameon=True, loc='upper right', fontsize=10)

    # PANEL B: WEIGHT DECAY CURVE
    curve_x = np.linspace(0, 5, 1000)
    curve_y = np.exp(-alpha * np.maximum(0, curve_x - tau_0))
    
    ax2.plot(curve_x, curve_y, color='#2ca02c', lw=3, label='Exponential Weight Decay ($w_g$)')
    ax2.fill_between(curve_x, 0, curve_y, color='#2ca02c', alpha=0.1)
    ax2.axvline(x=tau_0, color='#d62728', lw=2, linestyle='-')
    
    ax2.set_xlabel("Distortion Score ($D_g$)", fontsize=12)
    ax2.set_ylabel("Feature Weight ($w_g$)", fontsize=12)
    ax2.set_ylim(-0.05, 1.1)
    
    # Annotations
    ax2.text(0.5, 0.8, "SAFE ZONE ($w_g=1$)", color='green', fontweight='bold', fontsize=11)
    ax2.text(tau_0 + 0.2, 0.5, "SOFT PENALTY (Decay)", color='#d62728', fontweight='bold', fontsize=11)

    plt.tight_layout()
    
    save_path = "outputs/plots/gmm_soft_weights.png"
    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    plt.savefig(save_path, bbox_inches='tight')
    plt.close()
    print(f"--- Visualization saved to {save_path} ---")

if __name__ == "__main__":
    main()
