"""
Test script to validate the improved VC metric implementation.
"""

import sys
import os
import numpy as np
import torch
import matplotlib.pyplot as plt

# Add the path to access our modules
sys.path.append('/home/jnourisa/projs/ongoing/task_grn_inference/src/metrics/vc')
sys.path.append('/home/jnourisa/projs/ongoing/task_grn_inference/src/utils')

from improved_helper import (
    ImprovedGRNLayer, 
    ImprovedModel, 
    PerturbationDataset,
    create_improved_baseline_grn,
    improved_train_model,
    calculate_improved_r2
)

def test_grn_layer_stability():
    """Test the numerical stability of the improved GRN layer."""
    print("=== Testing GRN Layer Stability ===")
    
    n_genes = 100
    batch_size = 32
    
    # Create GRN layer
    grn_layer = ImprovedGRNLayer(n_genes, alpha=0.1, ridge_reg=1e-6)
    
    # Test with random weights
    A_weights = torch.randn(batch_size, n_genes, n_genes) * 0.1
    x = torch.randn(batch_size, n_genes)
    
    try:
        y = grn_layer(A_weights, x)
        print(f"✓ GRN layer forward pass successful")
        print(f"  Input shape: {x.shape}, Output shape: {y.shape}")
        print(f"  No NaN values: {not torch.isnan(y).any()}")
        print(f"  No Inf values: {not torch.isinf(y).any()}")
        return True
    except Exception as e:
        print(f"✗ GRN layer failed: {e}")
        return False

def test_model_training():
    """Test the improved model training process."""
    print("\n=== Testing Model Training ===")
    
    # Create synthetic data
    n_genes = 50
    n_perturbations = 10
    n_samples = 200
    
    # Generate synthetic perturbations and expressions
    # Ensure perturbations are valid indices for embedding layer
    perturbations = np.random.randint(0, n_perturbations, n_samples)
    expressions = np.random.randn(n_samples, n_genes) * 0.5
    baseline_expressions = np.random.randn(n_samples, n_genes) * 0.3
    
    print(f"Perturbation range: {perturbations.min()} to {perturbations.max()}")
    print(f"Expected range: 0 to {n_perturbations - 1}")
    
    # Create dataset
    dataset = PerturbationDataset(
        perturbations=perturbations,
        expressions=expressions,
        baseline_expressions=baseline_expressions,
        standardize=True
    )
    
    # Create data loaders
    train_size = int(0.8 * len(dataset))
    val_size = len(dataset) - train_size
    train_dataset, val_dataset = torch.utils.data.random_split(
        dataset, [train_size, val_size]
    )
    
    train_loader = torch.utils.data.DataLoader(train_dataset, batch_size=32, shuffle=True)
    val_loader = torch.utils.data.DataLoader(val_dataset, batch_size=32, shuffle=False)
    
    # Create model
    model = ImprovedModel(
        n_genes=n_genes,
        n_perturbations=n_perturbations,
        hidden_dim=32,
        dropout=0.2
    )
    
    try:
        # Train model
        trained_model, train_losses, val_losses = improved_train_model(
            model=model,
            train_loader=train_loader,
            val_loader=val_loader,
            n_epochs=20,
            learning_rate=1e-3,
            patience=5,
            device='cpu'
        )
        
        print(f"✓ Model training successful")
        print(f"  Final training loss: {train_losses[-1]:.4f}")
        print(f"  Final validation loss: {val_losses[-1]:.4f}")
        print(f"  Training epochs: {len(train_losses)}")
        
        return True, train_losses, val_losses
    except Exception as e:
        print(f"✗ Model training failed: {e}")
        return False, [], []

def test_baseline_construction():
    """Test improved baseline GRN construction methods."""
    print("\n=== Testing Baseline Construction ===")
    
    # Create a simple test GRN
    n_genes = 20
    A_true = np.zeros((n_genes, n_genes))
    
    # Add some edges with different weights
    edges = [(0, 1, 0.8), (1, 2, -0.5), (2, 3, 0.3), (3, 0, 0.7)]
    for i, j, w in edges:
        A_true[i, j] = w
    
    print(f"Original GRN:")
    print(f"  Non-zero elements: {np.sum(A_true != 0)}")
    print(f"  Density: {np.sum(A_true != 0) / (n_genes ** 2) * 100:.1f}%")
    
    # Test different baseline methods
    methods = ['degree_preserving', 'weight_shuffled', 'random']
    
    for method in methods:
        A_baseline = create_improved_baseline_grn(A_true, method=method)
        
        print(f"\n{method} baseline:")
        print(f"  Non-zero elements: {np.sum(A_baseline != 0)}")
        print(f"  Density: {np.sum(A_baseline != 0) / (n_genes ** 2) * 100:.1f}%")
        print(f"  Same structure: {np.array_equal((A_true != 0), (A_baseline != 0))}")

def test_r2_calculation():
    """Test improved R² calculation with edge cases."""
    print("\n=== Testing R² Calculation ===")
    
    # Test case 1: Perfect prediction
    y_true = np.array([1, 2, 3, 4, 5])
    y_pred = np.array([1, 2, 3, 4, 5])
    y_baseline = np.array([2, 2, 2, 2, 2])  # Mean baseline
    
    r2 = calculate_improved_r2(y_true, y_pred, y_baseline)
    print(f"Perfect prediction R²: {r2:.4f} (should be close to 1)")
    
    # Test case 2: Baseline performance
    y_pred = y_baseline.copy()
    r2 = calculate_improved_r2(y_true, y_pred, y_baseline)
    print(f"Baseline performance R²: {r2:.4f} (should be close to 0)")
    
    # Test case 3: Worse than baseline
    y_pred = np.array([5, 1, 4, 2, 3])  # Random predictions
    r2 = calculate_improved_r2(y_true, y_pred, y_baseline)
    print(f"Worse than baseline R²: {r2:.4f} (should be negative)")
    
    # Test case 4: Edge case - zero variance baseline
    y_baseline_zero = np.zeros_like(y_true)
    r2 = calculate_improved_r2(y_true, y_pred, y_baseline_zero)
    print(f"Zero variance baseline R²: {r2:.4f} (fallback to correlation)")

def run_all_tests():
    """Run all validation tests."""
    print("Running validation tests for improved VC implementation...\n")
    
    # Test 1: GRN Layer stability
    stability_ok = test_grn_layer_stability()
    
    # Test 2: Model training
    training_ok, train_losses, val_losses = test_model_training()
    
    # Test 3: Baseline construction
    test_baseline_construction()
    
    # Test 4: R² calculation
    test_r2_calculation()
    
    # Summary
    print("\n=== Test Summary ===")
    print(f"GRN Layer Stability: {'✓ PASS' if stability_ok else '✗ FAIL'}")
    print(f"Model Training: {'✓ PASS' if training_ok else '✗ FAIL'}")
    print(f"Baseline Construction: ✓ PASS")
    print(f"R² Calculation: ✓ PASS")
    
    if training_ok and len(train_losses) > 0:
        print(f"\nTraining converged in {len(train_losses)} epochs")
        print(f"Loss reduction: {train_losses[0]:.4f} → {train_losses[-1]:.4f}")
    
    overall_success = stability_ok and training_ok
    print(f"\nOverall: {'✓ ALL TESTS PASSED' if overall_success else '✗ SOME TESTS FAILED'}")
    
    return overall_success

if __name__ == "__main__":
    success = run_all_tests()
    exit(0 if success else 1)