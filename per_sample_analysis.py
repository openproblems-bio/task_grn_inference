"""
Analysis of Per-Sample Prediction VC Results
Comparing with previous ultra-fast implementation
"""

print("=== Per-Sample Prediction VC Analysis ===")
print()

# Current results (per-sample prediction)
current_results = {
    'positive_control': -0.0008,
    'negative_control': -0.0023,
    'pearson_corr': -0.0029
}

# Previous results (ultra-fast per-gene prediction)
previous_results = {
    'positive_control': -0.026,
    'negative_control': +0.012,
    'pearson_corr': +0.0001
}

print("1. SCORE COMPARISON:")
print("Method               | Per-Sample | Per-Gene  | Change")
print("-" * 55)
for method in current_results:
    curr = current_results[method]
    prev = previous_results[method]
    change = curr - prev
    print(f"{method:19} | {curr:+8.4f}   | {prev:+8.4f}  | {change:+7.4f}")

print()
print("2. RANKING ANALYSIS:")
print("Current (Per-Sample):")
print(f"  Positive control: {current_results['positive_control']:+.4f}")
print(f"  Negative control: {current_results['negative_control']:+.4f}")
print(f"  Difference: {current_results['positive_control'] - current_results['negative_control']:+.4f}")
print(f"  ✓ Positive > Negative: {current_results['positive_control'] > current_results['negative_control']}")

print()
print("Previous (Per-Gene):")
print(f"  Positive control: {previous_results['positive_control']:+.4f}")
print(f"  Negative control: {previous_results['negative_control']:+.4f}")
print(f"  Difference: {previous_results['positive_control'] - previous_results['negative_control']:+.4f}")
print(f"  ✗ Positive > Negative: {previous_results['positive_control'] > previous_results['negative_control']}")

print()
print("3. KEY IMPROVEMENTS:")
print("✓ FIXED RANKING: Positive control now scores higher than negative control")
print("✓ CONSISTENT SIGNS: All scores are negative (closer to expected behavior)")
print("✓ MAINTAINED SPEED: Still ~2-3 minutes per method (vs hours originally)")
print("✓ ARCHITECTURE: Single forward pass per sample (all genes simultaneously)")

print()
print("4. REMAINING CONCERNS:")
print("⚠️  Scores very close to zero (-0.0008 to -0.0029)")
print("⚠️  All negative R² values (worse than mean prediction)")
print("⚠️  Small differences between methods")

print()
print("5. TECHNICAL ANALYSIS:")
print("The per-sample prediction approach:")
print("- Predicts all 800 genes simultaneously for each sample")
print("- Uses matrix operations: X_pred = (I - αA^T)^(-1) @ X_pert")
print("- Computes MSE loss across all genes per sample")
print("- Maintains same mathematical framework but different granularity")

print()
print("6. NEXT STEPS:")
print("- Investigate why R² values are so close to zero")
print("- Check if this is expected behavior for VC metric")
print("- Consider if small negative values indicate model limitations")
print("- Validate against known good/bad GRNs with larger differences")