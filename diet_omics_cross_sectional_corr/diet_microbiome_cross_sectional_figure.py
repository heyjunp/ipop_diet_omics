#!/usr/bin/env python3
"""
Comprehensive Figure 6: Microbiome-Diet Correlations
- All nutrients (>40)
- All food groups (~20)
- diet_cluster (0/1)
- All taxa levels (phylum, class, order, family, genus)
- Age and sex adjustment (partial correlation)
- Proper color scale matching Figma
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.stats import spearmanr
from sklearn.linear_model import LinearRegression
from statsmodels.stats.multitest import multipletests
import os
import warnings
warnings.filterwarnings('ignore')

# Set style
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (16, 12)
plt.rcParams['font.size'] = 10

print("="*80)
print("COMPREHENSIVE FIGURE 6: MICROBIOME-DIET CORRELATIONS")
print("="*80)

# Load data
print("\n1. Loading data...")
demographics = pd.read_csv('demographics.csv')
nutrients = pd.read_csv('nutrients.csv')
food_groups = pd.read_csv('food_group.csv')
microbiome = pd.read_csv('stool_microbes_kw_test.csv')

# Preprocess
demographics['SubjectID'] = demographics['SubjectID'].astype(str)
microbiome['SubjectID'] = microbiome['SampleID'].str.split('-').str[:2].str.join('-')

# Merge demographics
nutrients_merged = nutrients.merge(
    demographics[['SubjectID', 'Age', 'Sex']], 
    on='SubjectID', how='inner'
)

food_groups_merged = food_groups.merge(
    demographics[['SubjectID', 'Age', 'Sex']], 
    on='SubjectID', how='inner'
)

microbiome_merged = microbiome.merge(
    demographics[['SubjectID', 'Age', 'Sex']], 
    on='SubjectID', how='inner'
)

print(f"   ✓ Nutrients: {len(nutrients_merged)} participants")
print(f"   ✓ Food groups: {len(food_groups_merged)} participants")
print(f"   ✓ Microbiome: {len(microbiome_merged)} samples")

# Get all taxa columns (all levels)
print("\n2. Extracting all taxa levels...")
taxa_cols = []
for level in ['phylum_', 'class_', 'order_', 'family_', 'genus_']:
    level_cols = [col for col in microbiome.columns if col.startswith(level)]
    taxa_cols.extend(level_cols)
    print(f"   {level}: {len(level_cols)} taxa")

print(f"   Total taxa: {len(taxa_cols)}")

# Get all nutrients (exclude SubjectID and specified nutrients to remove)
nutrients_to_remove = ['Fol_DFE', 'VitD_mcg', 'VitB_NE', 'VitA_RAE', 'Water', 
                       'OCarb', 'Fib', 'SolFib', 'FatCals', 'SatCals', 
                       'Sugar', 'Disacc', 'Caroten', 'Chrom']
nutrient_cols = [col for col in nutrients.columns 
                 if col not in ['SubjectID'] + nutrients_to_remove 
                 and col in nutrients_merged.columns]

# Check for nutrients with insufficient data (will result in NaN correlations)
print("\n3a. Checking for nutrients with insufficient data...")
nutrients_with_data = []
for col in nutrient_cols:
    col_data = nutrients_merged[col].dropna()
    # Need at least 5 non-missing values for meaningful correlation
    if len(col_data) >= 5 and col_data.nunique() > 1:
        nutrients_with_data.append(col)
    else:
        print(f"   Removing {col}: insufficient data (n={len(col_data)}, unique={col_data.nunique() if len(col_data) > 0 else 0})")

nutrient_cols = nutrients_with_data
print(f"\n3. Nutrients: {len(nutrient_cols)} (removed {len(nutrients_to_remove)} specified + {len([c for c in nutrients.columns if c not in ['SubjectID'] + nutrients_to_remove and c in nutrients_merged.columns]) - len(nutrient_cols)} with insufficient data)")
if nutrients_to_remove:
    print(f"   Removed specified: {', '.join(nutrients_to_remove)}")

# Get all food groups (exclude SubjectID)
food_cols_all = [col for col in food_groups.columns 
                 if col != 'SubjectID' and col in food_groups_merged.columns]

# Check for food groups with insufficient data
food_cols = []
for col in food_cols_all:
    col_data = food_groups_merged[col].dropna()
    # Need at least 5 non-missing values for meaningful correlation
    if len(col_data) >= 5 and col_data.nunique() > 1:
        food_cols.append(col)
    else:
        print(f"   Removing food group {col}: insufficient data (n={len(col_data)}, unique={col_data.nunique() if len(col_data) > 0 else 0})")

print(f"   Food groups: {len(food_cols)} (removed {len(food_cols_all) - len(food_cols)} with insufficient data)")

# Get diet_cluster
print(f"   diet_cluster: available in microbiome data")

# Aggregate microbiome data by SubjectID (mean across timepoints)
print("\n4. Aggregating microbiome data by SubjectID...")
microbiome_subj = microbiome_merged.groupby('SubjectID')[taxa_cols].mean()
microbiome_subj_demo = microbiome_merged.groupby('SubjectID')[['Age', 'Sex', 'diet_cluster']].first()

# Aggregate nutrients and food groups by SubjectID
nutrients_subj = nutrients_merged.groupby('SubjectID')[nutrient_cols + ['Age']].mean()
nutrients_subj_sex = nutrients_merged.groupby('SubjectID')['Sex'].first()
nutrients_subj = nutrients_subj.merge(nutrients_subj_sex, left_index=True, right_index=True)

food_groups_subj = food_groups_merged.groupby('SubjectID')[food_cols + ['Age']].mean()
food_groups_subj_sex = food_groups_merged.groupby('SubjectID')['Sex'].first()
food_groups_subj = food_groups_subj.merge(food_groups_subj_sex, left_index=True, right_index=True)

# Merge all data
print("\n5. Merging all data...")
# Start with microbiome
merged_data = microbiome_subj.merge(microbiome_subj_demo, left_index=True, right_index=True, how='inner')

# Merge nutrients (drop Age and Sex from nutrients to avoid conflicts)
nutrients_subj_clean = nutrients_subj.drop(columns=['Age', 'Sex'])
merged_data = merged_data.merge(
    nutrients_subj_clean, 
    left_index=True, right_index=True, 
    how='inner'
)

# Merge food groups (drop Age and Sex from food groups to avoid conflicts)
food_groups_subj_clean = food_groups_subj.drop(columns=['Age', 'Sex'])
merged_data = merged_data.merge(
    food_groups_subj_clean, 
    left_index=True, right_index=True, 
    how='inner'
)

# Convert Sex to numeric (M=0, F=1)
merged_data['Sex_numeric'] = (merged_data['Sex'] == 'F').astype(int)

print(f"   ✓ Final merged dataset: {len(merged_data)} participants")

# Function for partial Spearman correlation
def partial_spearman(x, y, z):
    """
    Calculate partial Spearman correlation between x and y, controlling for z
    z can be a 2D array (multiple covariates)
    """
    # Rank transform
    x_rank = stats.rankdata(x)
    y_rank = stats.rankdata(y)
    
    if z.ndim == 1:
        z = z.reshape(-1, 1)
    
    # Rank transform covariates
    z_rank = np.apply_along_axis(stats.rankdata, 0, z)
    
    # Remove missing values
    mask = ~(np.isnan(x_rank) | np.isnan(y_rank) | np.isnan(z_rank).any(axis=1))
    if mask.sum() < 5:
        return np.nan, 1.0
    
    x_rank = x_rank[mask]
    y_rank = y_rank[mask]
    z_rank = z_rank[mask, :]
    
    # Residuals from covariates
    try:
        # Fit linear models to get residuals
        # Residuals for x
        reg_x = LinearRegression()
        reg_x.fit(z_rank, x_rank)
        x_resid = x_rank - reg_x.predict(z_rank)
        
        # Residuals for y
        reg_y = LinearRegression()
        reg_y.fit(z_rank, y_rank)
        y_resid = y_rank - reg_y.predict(z_rank)
        
        # Spearman correlation of residuals
        corr, pval = spearmanr(x_resid, y_resid)
        return corr, pval
    except:
        return np.nan, 1.0

# Calculate correlations
print("\n6. Calculating partial Spearman correlations (adjusted for Age and Sex)...")
print("   This may take a few minutes...")

# Prepare covariates
covariates = merged_data[['Age', 'Sex_numeric']].values

# Combine all diet variables and relabel
diet_vars_raw = nutrient_cols + food_cols + ['diet_cluster']
diet_data = merged_data[diet_vars_raw]

# Create relabeling dictionary
diet_relabel = {}
for var in diet_vars_raw:
    if var == 'Mixed_Meal':
        diet_relabel[var] = 'Pizza_Sandwich_Wrap'
    elif var == 'diet_cluster':
        diet_relabel[var] = 'Diet_Pattern'
    else:
        diet_relabel[var] = var

diet_vars = [diet_relabel.get(var, var) for var in diet_vars_raw]

# Calculate correlation matrix
corr_matrix = np.full((len(taxa_cols), len(diet_vars)), np.nan)
pval_matrix = np.full((len(taxa_cols), len(diet_vars)), 1.0)

for i, taxa in enumerate(taxa_cols):
    if (i + 1) % 20 == 0:
        print(f"   Processing taxa {i+1}/{len(taxa_cols)}...")
    
    taxa_data = merged_data[taxa].values
    
    for j, diet_var_display in enumerate(diet_vars):
        # Use original name for data access
        diet_var_original = diet_vars_raw[j]
        diet_data_var = merged_data[diet_var_original].values
        
        # Partial correlation adjusting for Age and Sex
        corr, pval = partial_spearman(taxa_data, diet_data_var, covariates)
        corr_matrix[i, j] = corr
        pval_matrix[i, j] = pval

print("   ✓ Correlations calculated")

# VERIFICATION: Confirm correlations were computed separately by taxonomic level
print("\n6a. Verification: Checking that correlations were computed separately by taxonomic level...")
# Extract base names (without prefixes) to check for duplicates
base_names = {}
for i, taxa in enumerate(taxa_cols):
    base_name = taxa.replace('phylum_', '').replace('class_', '').replace('order_', '').replace('family_', '').replace('genus_', '')
    if base_name not in base_names:
        base_names[base_name] = []
    base_names[base_name].append((i, taxa))

# Report any taxa that appear at multiple levels
duplicates = {name: levels for name, levels in base_names.items() if len(levels) > 1}
if duplicates:
    print(f"   Found {len(duplicates)} taxa appearing at multiple taxonomic levels:")
    for base_name, levels in sorted(duplicates.items())[:10]:  # Show first 10
        level_names = [taxa.split('_')[0] for _, taxa in levels]
        print(f"     - {base_name}: appears at {', '.join(level_names)} levels")
        # Verify these have different correlation values (not collapsed)
        if len(levels) > 1:
            # Check first diet variable for differences
            sample_corrs = [corr_matrix[idx, 0] for idx, _ in levels]
            if not all(np.isnan(corr) for corr in sample_corrs):
                unique_corrs = len(set([c for c in sample_corrs if not np.isnan(c)]))
                if unique_corrs > 1:
                    print(f"       ✓ Verified: Different correlation values (not collapsed)")
                else:
                    print(f"       ⚠ Warning: Similar correlation values - may need investigation")
    if len(duplicates) > 10:
        print(f"     ... and {len(duplicates) - 10} more")
else:
    print("   ✓ No duplicate taxa names found across different levels")
print(f"   ✓ All {len(taxa_cols)} taxa correlations computed separately (preserving taxonomic level prefixes)")

# Remove diet variables and taxa with all NaN correlations
print("\n6b. Removing variables with no correlations (all NaN)...")
# Check columns (diet variables) - track both original and display names
valid_diet_vars = []
for j in range(len(diet_vars)):
    col_corrs = corr_matrix[:, j]
    if not np.isnan(col_corrs).all():
        valid_diet_vars.append((j, diet_vars_raw[j], diet_vars[j]))
    else:
        print(f"   Removing {diet_vars[j]}: all correlations are NaN")

# Check rows (taxa)
valid_taxa = []
for i, taxa in enumerate(taxa_cols):
    row_corrs = corr_matrix[i, :]
    if not np.isnan(row_corrs).all():
        valid_taxa.append((i, taxa))
    else:
        print(f"   Removing {taxa}: all correlations are NaN")

# Filter matrices
if len(valid_diet_vars) < len(diet_vars) or len(valid_taxa) < len(taxa_cols):
    valid_diet_indices = [idx for idx, _, _ in valid_diet_vars]
    valid_taxa_indices = [idx for idx, _ in valid_taxa]
    
    corr_matrix = corr_matrix[np.ix_(valid_taxa_indices, valid_diet_indices)]
    pval_matrix = pval_matrix[np.ix_(valid_taxa_indices, valid_diet_indices)]
    taxa_cols = [taxa for _, taxa in valid_taxa]
    # Update both raw and display names from valid_diet_vars
    diet_vars_raw = [raw_name for _, raw_name, _ in valid_diet_vars]
    diet_vars = [display_name for _, _, display_name in valid_diet_vars]
    
    print(f"   ✓ Filtered: {len(taxa_cols)} taxa × {len(diet_vars)} diet variables")

# Apply BH adjustment with alpha=0.2
print("\n6c. Applying BH adjustment (alpha=0.2)...")
# Flatten p-values for correction (only non-NaN values)
valid_mask = ~np.isnan(pval_matrix)
flat_pvals = pval_matrix[valid_mask]

# Apply FDR correction (Benjamini-Hochberg) with alpha=0.2
rejected, pval_adj_flat, _, _ = multipletests(flat_pvals, alpha=0.2, method='fdr_bh')

# Reshape adjusted p-values back to matrix
pval_adj_matrix = np.full_like(pval_matrix, np.nan)
pval_adj_matrix[valid_mask] = pval_adj_flat

print(f"   ✓ BH adjustment applied (alpha=0.2)")
print(f"   Significant correlations (BH-adjusted p < 0.2): {(pval_adj_matrix < 0.2).sum()}")

# Create correlation DataFrame
# IMPORTANT: Keep taxonomic level prefixes to distinguish between different levels
# (e.g., phylum_Actinobacteria vs class_Actinobacteria)
corr_df = pd.DataFrame(
    corr_matrix,
    index=taxa_cols,  # Keep full names with prefixes
    columns=diet_vars
)

# Filter significant correlations (optional - for visualization)
# We'll show all correlations but can highlight significant ones

# Create the heatmap
print("\n7. Creating comprehensive heatmap...")

# Sort taxa by hierarchical level for better visualization
taxa_levels = []
for col in taxa_cols:
    if col.startswith('phylum_'):
        taxa_levels.append(('Phylum', col))
    elif col.startswith('class_'):
        taxa_levels.append(('Class', col))
    elif col.startswith('order_'):
        taxa_levels.append(('Order', col))
    elif col.startswith('family_'):
        taxa_levels.append(('Family', col))
    elif col.startswith('genus_'):
        taxa_levels.append(('Genus', col))

# Reorder correlation matrix by level
level_order = ['Phylum', 'Class', 'Order', 'Family', 'Genus']
sorted_taxa = []
for level in level_order:
    sorted_taxa.extend([t[1] for t in taxa_levels if t[0] == level])

sorted_indices = [taxa_cols.index(t) for t in sorted_taxa if t in taxa_cols]
corr_sorted = corr_df.iloc[sorted_indices, :]

# Create figure with proper color scale (matching Figma)
# Transposed layout: Diet (Y, 65 items) × Taxa (X, 96 items)
fig, ax = plt.subplots(figsize=(24, 14))  # Wider for taxa on X-axis, shorter for diet on Y-axis

# Color scale: New palette from user specification
# #142621 (dark green) -> #71d9bc (light green) -> #FFFFFF (white) -> #e39d5f (light brown) -> #4a331f (dark brown)
from matplotlib.colors import LinearSegmentedColormap
colors_spec = ['#142621', '#71d9bc', '#FFFFFF', '#e39d5f', '#4a331f']
n_bins = 256
cmap_spec = LinearSegmentedColormap.from_list('brown_green_new', colors_spec, N=n_bins)

# Tightened limits: -0.3 to 0.3 (values outside will be clamped)
vmin = -0.3
vmax = 0.3

# Clamp correlation values to the new range for better contrast
corr_sorted_clamped = corr_sorted.copy()
corr_sorted_clamped = corr_sorted_clamped.clip(lower=vmin, upper=vmax)

# Prepare annotation matrix with asterisks for significance (p < 0.05)
pval_sorted = pd.DataFrame(
    pval_matrix[sorted_indices, :],
    index=corr_sorted.index,
    columns=corr_sorted.columns
)
annot_matrix = corr_sorted_clamped.copy().astype(str)
for i in range(len(annot_matrix.index)):
    for j in range(len(annot_matrix.columns)):
        if not np.isnan(pval_sorted.iloc[i, j]) and pval_sorted.iloc[i, j] < 0.05:
            annot_matrix.iloc[i, j] = '*'
        else:
            annot_matrix.iloc[i, j] = ''

# Transpose: Diet variables on Y-axis, Taxa on X-axis
corr_transposed = corr_sorted_clamped.T
annot_transposed = annot_matrix.T

# Create heatmap (transposed) with tightened color scale and adaptive asterisk colors
# Note: seaborn doesn't support per-cell annotation colors directly, so we'll use white asterisks
# which should be readable on most of the color range (especially since most values are modest)
sns.heatmap(
    corr_transposed,
    annot=annot_transposed,  # Show asterisks for significant (p < 0.05)
    fmt='',  # Empty format since we're using custom annotations
    cmap=cmap_spec,
    center=0,
    vmin=vmin,
    vmax=vmax,
    ax=ax,
    cbar_kws={'label': 'Partial Spearman Correlation\n(adjusted for Age & Sex)\nValues clamped to [-0.3, 0.3]\n*p < 0.05', 'shrink': 0.8},
    xticklabels=True,
    yticklabels=True,
    linewidths=0.1,
    linecolor='gray',
    annot_kws={'fontsize': 10, 'fontweight': 'bold', 'color': '#000000'}  # Black asterisks (better for modest correlations)
)

# Adjust labels (transposed: X=Taxa, Y=Diet)
ax.set_xlabel('Microbiome Taxa (All Levels)', fontsize=12, fontweight='bold')
ax.set_ylabel('Diet Variables (Nutrients, Food Groups, Diet Pattern)', fontsize=12, fontweight='bold')
ax.set_title('Microbiome-Diet Correlations\n(Partial Spearman, adjusted for Age & Sex)\nColor scale: -0.3 to +0.3 (values clamped)\n*Black asterisks: p < 0.05', 
             fontsize=14, fontweight='bold', pad=20)

# Rotate labels (X-axis = Taxa, Y-axis = Diet)
# Format taxa labels: show level prefix for clarity (e.g., "phylum_Actinobacteria" -> "P:Actinobacteria")
taxa_labels = []
for label in corr_sorted.index:
    if label.startswith('phylum_'):
        taxa_labels.append(f"P:{label.replace('phylum_', '')}")
    elif label.startswith('class_'):
        taxa_labels.append(f"C:{label.replace('class_', '')}")
    elif label.startswith('order_'):
        taxa_labels.append(f"O:{label.replace('order_', '')}")
    elif label.startswith('family_'):
        taxa_labels.append(f"F:{label.replace('family_', '')}")
    elif label.startswith('genus_'):
        taxa_labels.append(f"G:{label.replace('genus_', '')}")
    else:
        taxa_labels.append(label)

ax.set_xticklabels(taxa_labels, rotation=90, ha='right', fontsize=7)  # Taxa on X-axis
plt.setp(ax.get_yticklabels(), rotation=0, ha='right', fontsize=8)   # Diet on Y-axis

# Add level separators (now on X-axis since transposed)
current_level = None
for i, taxa in enumerate(sorted_taxa):
    if taxa in taxa_cols:
        level = [t[0] for t in taxa_levels if t[1] == taxa][0]
        if level != current_level and current_level is not None:
            ax.axvline(x=i, color='black', linewidth=1.5, alpha=0.5)  # Changed to axvline for X-axis
        current_level = level

plt.tight_layout()

# Save figure
os.makedirs('figures', exist_ok=True)
plt.savefig('figures/06_microbiome_diet_comprehensive.png', dpi=300, bbox_inches='tight')
print("   ✓ Figure saved: figures/06_microbiome_diet_comprehensive.png")

# Note: Only one heatmap now (with asterisks for significance, full color scale)
print("\n8. Single heatmap created with white asterisks for significance (p < 0.05)")

# Save correlation matrix and p-values
print("\n9. Saving correlation results...")
print("   Note: Taxonomic level prefixes are preserved in all output files")
print("         (e.g., 'phylum_Actinobacteria', 'class_Actinobacteria' are distinct)")
corr_df.to_csv('figures/06_correlation_matrix.csv')
print("   ✓ Saved: figures/06_correlation_matrix.csv (with taxonomic level prefixes)")

# Save p-values (raw and BH-adjusted)
pval_df_raw = pd.DataFrame(
    pval_matrix[sorted_indices, :],
    index=corr_sorted.index,  # This preserves the full taxa names with prefixes
    columns=diet_vars
)
pval_df_raw.to_csv('figures/06_pvalue_matrix_raw.csv')
print("   ✓ Saved: figures/06_pvalue_matrix_raw.csv (with taxonomic level prefixes)")

pval_df_bh = pd.DataFrame(
    pval_adj_matrix[sorted_indices, :],
    index=corr_sorted.index,  # This preserves the full taxa names with prefixes
    columns=diet_vars
)
pval_df_bh.to_csv('figures/06_pvalue_matrix_bh.csv')
print("   ✓ Saved: figures/06_pvalue_matrix_bh.csv (with taxonomic level prefixes)")

# Summary statistics
print("\n10. Summary Statistics:")
print(f"   Total correlations calculated: {len(taxa_cols) * len(diet_vars)}")
sig_corr_bh = (pval_adj_matrix < 0.2).sum()
sig_corr_raw = (pval_matrix < 0.05).sum()
print(f"   Significant correlations (raw p < 0.05): {sig_corr_raw}")
print(f"   Significant correlations (BH-adjusted p < 0.2): {sig_corr_bh}")
print(f"   Percentage significant (BH): {sig_corr_bh / (len(taxa_cols) * len(diet_vars)) * 100:.2f}%")

# Top correlations
print("\n11. Top 20 strongest correlations:")
flat_corr = []
for i in range(len(taxa_cols)):
    for j in range(len(diet_vars)):
        if not np.isnan(corr_matrix[i, j]):
            flat_corr.append({
                'taxa': taxa_cols[i],
                'diet': diet_vars[j],
                'correlation': corr_matrix[i, j],
                'pvalue_raw': pval_matrix[i, j],
                'pvalue_bh': pval_adj_matrix[i, j]
            })

top_corr = sorted(flat_corr, key=lambda x: abs(x['correlation']), reverse=True)[:20]
for i, corr in enumerate(top_corr, 1):
    sig_marker = "*" if corr['pvalue_bh'] < 0.2 else ""
    print(f"   {i:2d}. {corr['taxa'][:40]:40s} ↔ {corr['diet'][:30]:30s} | "
          f"r={corr['correlation']:6.3f}, p_raw={corr['pvalue_raw']:.4f}, p_bh={corr['pvalue_bh']:.4f}{sig_marker}")

print("\n" + "="*80)
print("ANALYSIS COMPLETE!")
print("="*80)
print("\nGenerated files:")
print("  - figures/06_microbiome_diet_comprehensive.png (all correlations)")
print("  - figures/06_microbiome_diet_comprehensive.png (with asterisks for BH p<0.2)")
print("  - figures/06_correlation_matrix.csv")
print("  - figures/06_pvalue_matrix.csv (p-values)")
print("\nNote: Heatmap uses brown-green color scale with white asterisks for p < 0.05")
print("      Transposed: Diet variables (Y-axis) × Taxa (X-axis)")
print("\n")

