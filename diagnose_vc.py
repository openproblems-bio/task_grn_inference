import numpy as np
import pandas as pd
import torch

# Simple diagnostic script to understand VC metric issues
def diagnose_vc_issues():
    print("=== VC Metric Diagnostic Analysis ===\n")
    
    # 1. Test GRNLayer stability
    print("1. Testing GRNLayer mathematical stability:")
    
    # Create a simple test GRN
    n_genes = 100
    A = np.random.randn(n_genes, n_genes) * 0.1  # Small random weights
    A_signs = torch.from_numpy(np.sign(A).astype(np.float32))
    A_weights = torch.nn.Parameter(torch.from_numpy(A.astype(np.float32)))
    
    # Test matrix inversion stability
    I = torch.eye(n_genes, dtype=torch.float32)
    alpha = 1.0
    
    try:
        ia = I - alpha * A_weights.t()
        eigenvals = torch.linalg.eigvals(ia).real
        min_eigval = torch.min(eigenvals)
        cond_number = torch.linalg.cond(ia)
        
        print(f"   Minimum eigenvalue: {min_eigval:.6f}")
        print(f"   Condition number: {cond_number:.2e}")
        print(f"   Matrix is {'invertible' if min_eigval > 1e-10 else 'near-singular'}")
        
        if min_eigval > 1e-10:
            ia_inv = torch.linalg.inv(ia)
            print(f"   Inversion successful")
        else:
            print(f"   ⚠️  Matrix inversion unstable!")
    except Exception as e:
        print(f"   ❌ Matrix inversion failed: {e}")
    
    print()
    
    # 2. Test typical GRN properties
    print("2. Analyzing typical GRN properties:")
    
    # Simulate different GRN densities
    for density in [0.01, 0.05, 0.1]:
        A_sparse = np.random.randn(n_genes, n_genes) * 0.1
        mask = np.random.random((n_genes, n_genes)) > density
        A_sparse[mask] = 0
        
        n_edges = np.sum(A_sparse != 0)
        print(f"   Density {density:.2%}: {n_edges} edges")
        
        # Test stability
        A_tensor = torch.from_numpy(A_sparse.astype(np.float32))
        ia = I - A_tensor.t()
        try:
            min_eigval = torch.min(torch.linalg.eigvals(ia).real)
            print(f"   -> Min eigenvalue: {min_eigval:.6f}")
        except:
            print(f"   -> Eigenvalue computation failed")
    
    print()
    
    # 3. Test baseline comparison logic
    print("3. Testing R² calculation logic:")
    
    # Simulate some losses
    ss_res_scenarios = [100, 200, 150]  # Actual GRN losses
    ss_tot_scenarios = [100, 150, 100]  # Baseline losses
    
    for ss_res, ss_tot in zip(ss_res_scenarios, ss_tot_scenarios):
        r2 = 1 - ss_res / ss_tot
        print(f"   GRN loss: {ss_res}, Baseline loss: {ss_tot} -> R²: {r2:.3f}")
        
        if r2 < 0:
            print(f"   ⚠️  Negative R²: GRN performs worse than baseline")
        elif abs(r2) < 0.01:
            print(f"   ⚠️  Near-zero R²: GRN performs similar to baseline")
    
    print()
    
    # 4. Recommendations
    print("4. Potential Solutions:")
    print("   a) Add regularization to prevent singular matrices")
    print("   b) Use smaller alpha values (e.g., 0.1 instead of 1.0)")
    print("   c) Add numerical stabilization (ridge regularization)")
    print("   d) Use different baseline construction")
    print("   e) Simplify model architecture")
    print("   f) Add early stopping based on validation R²")

if __name__ == "__main__":
    diagnose_vc_issues()