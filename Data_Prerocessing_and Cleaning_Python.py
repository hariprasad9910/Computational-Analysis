"""
Data Preprocessing and Cleaning Pipeline
----------------------------------------
This script provides a comprehensive data preprocessing workflow for CSV datasets.
It performs:
1. Data loading and inspection
2. Missing value analysis and treatment
3. Outlier detection and handling
4. Feature engineering basics
5. Data transformation and scaling
6. Exploratory data visualization
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import MinMaxScaler, StandardScaler
import warnings

warnings.filterwarnings('ignore')
plt.style.use('ggplot')


def load_data(file_path):
    """Load CSV data and display basic information."""
    print("=" * 80)
    print("LOADING DATA")
    print("=" * 80)
    
    try:
        df = pd.read_csv(file_path)
        print(f"Successfully loaded {file_path}")
        print(f"Dataset shape: {df.shape}")
        print("\nFirst 5 rows:")
        print(df.head())
        return df
    except Exception as e:
        print(f"Error loading data: {e}")
        return None


def inspect_data(df):
    """Perform basic data inspection."""
    print("\n", "=" * 80)
    print("DATA INSPECTION")
    print("=" * 80)
    
    # Data information
    print("\nData Information:")
    print("-" * 50)
    df.info()
    
    # Check for duplicates
    duplicates = df.duplicated().sum()
    print(f"\nDuplicate rows: {duplicates}")
    
    # Identify column types
    cat_cols = [col for col in df.columns if df[col].dtype == 'object']
    num_cols = [col for col in df.columns if df[col].dtype != 'object']
    
    print(f"\nCategorical columns ({len(cat_cols)}): {cat_cols}")
    print(f"Numerical columns ({len(num_cols)}): {num_cols}")
    
    # Unique values in categorical columns
    if cat_cols:
        print("\nUnique values in categorical columns:")
        for col in cat_cols:
            unique_vals = df[col].nunique()
            print(f"{col}: {unique_vals} unique values")
    
    # Basic statistics for numerical columns
    print("\nStatistical Summary:")
    print("-" * 50)
    print(df.describe())
    
    return cat_cols, num_cols


def analyze_missing_values(df):
    """Analyze and handle missing values."""
    print("\n", "=" * 80)
    print("MISSING VALUE ANALYSIS")
    print("=" * 80)
    
    # Check for missing values
    missing = df.isnull().sum()
    missing_percent = (df.isnull().sum() / df.shape[0]) * 100
    
    # Create a dataframe with missing value information
    missing_df = pd.DataFrame({
        'Missing Values': missing,
        'Percentage': missing_percent.round(2)
    })
    
    missing_df = missing_df[missing_df['Missing Values'] > 0].sort_values(
        'Percentage', ascending=False)
    
    if missing_df.empty:
        print("No missing values found in the dataset.")
    else:
        print("Missing values by column:")
        print(missing_df)
    
    return missing_df


def handle_missing_values(df, missing_df, strategy='auto'):
    """Handle missing values using specified strategies."""
    print("\n", "=" * 80)
    print("HANDLING MISSING VALUES")
    print("=" * 80)
    
    if missing_df.empty:
        print("No missing values to handle.")
        return df
    
    df_clean = df.copy()
    
    # For each column with missing values
    for col in missing_df.index:
        # Get percentage of missing values
        missing_percent = missing_df.loc[col, 'Percentage']
        
        # If strategy is auto, choose based on missing percentage and data type
        if strategy == 'auto':
            # If more than 50% missing, consider dropping the column
            if missing_percent > 50:
                print(f"Dropping column {col} (missing: {missing_percent}%)")
                df_clean = df_clean.drop(columns=[col])
            
            # For less than 50% missing
            else:
                # For numerical columns
                if df[col].dtype != 'object':
                    # Check for outliers to decide between mean and median
                    q1 = df[col].dropna().quantile(0.25)
                    q3 = df[col].dropna().quantile(0.75)
                    iqr = q3 - q1
                    if (df[col].dropna() > (q3 + 1.5 * iqr)).sum() > 0 or \
                       (df[col].dropna() < (q1 - 1.5 * iqr)).sum() > 0:
                        # Use median for data with outliers
                        fill_value = df[col].median()
                        print(f"Filling {col} with median value: {fill_value}")
                    else:
                        # Use mean for normally distributed data
                        fill_value = df[col].mean()
                        print(f"Filling {col} with mean value: {fill_value}")
                    df_clean[col] = df_clean[col].fillna(fill_value)
                
                # For categorical columns
                else:
                    # Fill with mode (most frequent value)
                    fill_value = df[col].mode()[0]
                    print(f"Filling {col} with mode value: {fill_value}")
                    df_clean[col] = df_clean[col].fillna(fill_value)
        
        # If specific strategy provided
        else:
            if strategy == 'drop_columns':
                print(f"Dropping column {col}")
                df_clean = df_clean.drop(columns=[col])
            elif strategy == 'drop_rows':
                print(f"Dropping rows with missing values in {col}")
                df_clean = df_clean.dropna(subset=[col])
            elif strategy == 'mean':
                if df[col].dtype != 'object':
                    fill_value = df[col].mean()
                    print(f"Filling {col} with mean value: {fill_value}")
                    df_clean[col] = df_clean[col].fillna(fill_value)
                else:
                    print(f"Cannot apply mean to categorical column {col}")
            elif strategy == 'median':
                if df[col].dtype != 'object':
                    fill_value = df[col].median()
                    print(f"Filling {col} with median value: {fill_value}")
                    df_clean[col] = df_clean[col].fillna(fill_value)
                else:
                    print(f"Cannot apply median to categorical column {col}")
            elif strategy == 'mode':
                fill_value = df[col].mode()[0]
                print(f"Filling {col} with mode value: {fill_value}")
                df_clean[col] = df_clean[col].fillna(fill_value)
            elif strategy == 'zero':
                print(f"Filling {col} with zero")
                df_clean[col] = df_clean[col].fillna(0)
            else:
                print(f"Unknown strategy {strategy} for {col}, skipping")
    
    # Show the missing values after handling
    missing_after = df_clean.isnull().sum().sum()
    print(f"\nTotal missing values remaining: {missing_after}")
    
    return df_clean


def detect_outliers(df, num_cols):
    """Identify and visualize outliers in numerical columns."""
    print("\n", "=" * 80)
    print("OUTLIER DETECTION")
    print("=" * 80)
    
    if not num_cols:
        print("No numerical columns to analyze for outliers.")
        return {}
    
    outlier_info = {}
    
    # Set up the figure size based on number of columns
    num_plots = len(num_cols)
    fig_height = max(6, 2 * num_plots)
    
    # Create boxplots for each numerical column
    plt.figure(figsize=(10, fig_height))
    
    for i, col in enumerate(num_cols, 1):
        # Calculate outlier boundaries using IQR method
        q1 = df[col].quantile(0.25)
        q3 = df[col].quantile(0.75)
        iqr = q3 - q1
        lower_bound = q1 - 1.5 * iqr
        upper_bound = q3 + 1.5 * iqr
        
        # Count outliers
        outliers = df[(df[col] < lower_bound) | (df[col] > upper_bound)]
        outlier_count = len(outliers)
        outlier_percent = (outlier_count / len(df)) * 100
        
        outlier_info[col] = {
            'lower_bound': lower_bound,
            'upper_bound': upper_bound,
            'count': outlier_count,
            'percentage': outlier_percent
        }
        
        print(f"{col}: {outlier_count} outliers ({outlier_percent:.2f}%)")
        print(f"  Lower bound: {lower_bound}, Upper bound: {upper_bound}")
        
        # Create subplot for boxplot
        plt.subplot(num_plots, 1, i)
        sns.boxplot(x=df[col])
        plt.title(f'Boxplot of {col}')
        plt.tight_layout()
    
    plt.suptitle('Outlier Detection with Boxplots', y=1.02)
    plt.tight_layout()
    plt.show()
    
    return outlier_info


def handle_outliers(df, outlier_info, strategy='auto'):
    """Handle outliers using specified strategy."""
    print("\n", "=" * 80)
    print("HANDLING OUTLIERS")
    print("=" * 80)
    
    if not outlier_info:
        print("No outliers to handle.")
        return df
    
    df_clean = df.copy()
    
    for col, info in outlier_info.items():
        # Skip if no outliers
        if info['count'] == 0:
            continue
            
        # Auto-select strategy based on percentage of outliers
        if strategy == 'auto':
            if info['percentage'] < 1:
                # Very few outliers - can be capped
                method = 'cap'
            elif info['percentage'] < 5:
                # Moderate outliers - can use winsorization
                method = 'winsorize'
            else:
                # Many outliers - consider transformation
                method = 'log_transform'
        else:
            method = strategy
        
        # Apply the selected method
        if method == 'cap':
            print(f"Capping outliers in {col}")
            df_clean[col] = df_clean[col].clip(
                lower=info['lower_bound'], upper=info['upper_bound'])
        
        elif method == 'winsorize':
            print(f"Winsorizing {col} at 5th and 95th percentiles")
            p05 = df_clean[col].quantile(0.05)
            p95 = df_clean[col].quantile(0.95)
            df_clean[col] = df_clean[col].clip(lower=p05, upper=p95)
        
        elif method == 'log_transform':
            print(f"Applying log transformation to {col}")
            # Handle negative values if present
            min_val = df_clean[col].min()
            if min_val <= 0:
                shift = abs(min_val) + 1
                print(f"  Shifting values by {shift} before log transform")
                df_clean[col] = np.log(df_clean[col] + shift)
            else:
                df_clean[col] = np.log(df_clean[col])
        
        elif method == 'remove':
            print(f"Removing rows with outliers in {col}")
            df_clean = df_clean[
                (df_clean[col] >= info['lower_bound']) & 
                (df_clean[col] <= info['upper_bound'])
            ]
        
        # You can add more methods as needed
    
    return df_clean


def encode_categorical_features(df, cat_cols):
    """Encode categorical features."""
    print("\n", "=" * 80)
    print("ENCODING CATEGORICAL FEATURES")
    print("=" * 80)
    
    if not cat_cols:
        print("No categorical columns to encode.")
        return df
    
    df_encoded = df.copy()
    
    for col in cat_cols:
        n_unique = df[col].nunique()
        
        # For binary categorical variables
        if n_unique == 2:
            print(f"Label encoding {col} (binary)")
            labels = {val: idx for idx, val in enumerate(df[col].unique())}
            df_encoded[col] = df[col].map(labels)
            print(f"  Mapping: {labels}")
        
        # For nominal categorical variables with few categories
        elif n_unique <= 10:
            print(f"One-hot encoding {col} ({n_unique} categories)")
            # Get dummies and drop one category to avoid the dummy variable trap
            dummies = pd.get_dummies(df[col], prefix=col, drop_first=True)
            # Drop original column and join the dummy variables
            df_encoded = df_encoded.drop(col, axis=1)
            df_encoded = pd.concat([df_encoded, dummies], axis=1)
        
        # For high cardinality categorical variables
        else:
            print(f"Frequency encoding {col} (high cardinality: {n_unique} categories)")
            # Calculate frequency of each category
            freq_map = df[col].value_counts(normalize=True).to_dict()
            df_encoded[col] = df[col].map(freq_map)
    
    print(f"\nShape after encoding: {df_encoded.shape}")
    print(f"New columns: {list(set(df_encoded.columns) - set(df.columns))}")
    
    return df_encoded


def normalize_features(df, num_cols, method='standard'):
    """Normalize numerical features."""
    print("\n", "=" * 80)
    print("NORMALIZING FEATURES")
    print("=" * 80)
    
    if not num_cols:
        print("No numerical columns to normalize.")
        return df
    
    df_scaled = df.copy()
    
    # Identify columns to scale (only numerical columns that exist in the dataframe)
    cols_to_scale = [col for col in num_cols if col in df.columns]
    
    if method == 'minmax':
        print("Applying MinMax scaling (0-1 range)")
        scaler = MinMaxScaler()
        df_scaled[cols_to_scale] = scaler.fit_transform(df[cols_to_scale])
    
    elif method == 'standard':
        print("Applying StandardScaler (mean=0, std=1)")
        scaler = StandardScaler()
        df_scaled[cols_to_scale] = scaler.fit_transform(df[cols_to_scale])
    
    elif method == 'robust':
        print("Applying Robust scaling (based on quantiles)")
        from sklearn.preprocessing import RobustScaler
        scaler = RobustScaler()
        df_scaled[cols_to_scale] = scaler.fit_transform(df[cols_to_scale])
    
    return df_scaled


def visualize_data(df, cat_cols, num_cols):
    """Create basic visualizations of the data."""
    print("\n", "=" * 80)
    print("DATA VISUALIZATION")
    print("=" * 80)
    
    # Distribution of numerical features
    if num_cols:
        print("\nVisualizing distributions of numerical features...")
        n_cols = min(len(num_cols), 3)
        n_rows = (len(num_cols) + n_cols - 1) // n_cols
        
        plt.figure(figsize=(15, n_rows * 4))
        
        for i, col in enumerate(num_cols, 1):
            if col in df.columns:  # Check if column exists
                plt.subplot(n_rows, n_cols, i)
                sns.histplot(df[col], kde=True)
                plt.title(f'Distribution of {col}')
                plt.tight_layout()
        
        plt.show()
    
    # Categorical feature distributions
    if cat_cols:
        print("\nVisualizing distributions of categorical features...")
        valid_cat_cols = [col for col in cat_cols if col in df.columns]
        
        for col in valid_cat_cols:
            if df[col].nunique() <= 10:  # Only plot if there are 10 or fewer categories
                plt.figure(figsize=(12, 6))
                order = df[col].value_counts().index
                ax = sns.countplot(data=df, x=col, order=order)
                plt.title(f'Distribution of {col}')
                
                # Rotate x-axis labels for better visibility
                plt.xticks(rotation=45, ha='right')
                
                # Add count labels on top of bars
                for p in ax.patches:
                    ax.annotate(f'{int(p.get_height())}', 
                                (p.get_x() + p.get_width() / 2., p.get_height()),
                                ha='center', va='bottom')
                
                plt.tight_layout()
                plt.show()
    
    # Correlation matrix for numerical features
    if len(num_cols) > 1:
        print("\nCorrelation matrix for numerical features...")
        valid_num_cols = [col for col in num_cols if col in df.columns]
        
        plt.figure(figsize=(12, 10))
        corr = df[valid_num_cols].corr()
        mask = np.triu(np.ones_like(corr, dtype=bool))
        sns.heatmap(corr, mask=mask, cmap='coolwarm', annot=True, fmt='.2f',
                   linewidths=0.5, cbar_kws={'shrink': .8})
        plt.title('Correlation Matrix')
        plt.tight_layout()
        plt.show()


def preprocess_data(file_path, missing_strategy='auto', outlier_strategy='auto', 
                   scaling_method='standard'):
    """Run the complete preprocessing pipeline."""
    print("\n", "=" * 80)
    print(f"PREPROCESSING PIPELINE FOR: {file_path}")
    print("=" * 80)
    
    # Step 1: Load data
    df = load_data(file_path)
    if df is None:
        return None
    
    # Step 2: Inspect data
    cat_cols, num_cols = inspect_data(df)
    
    # Step 3: Analyze missing values
    missing_df = analyze_missing_values(df)
    
    # Step 4: Handle missing values
    df_clean = handle_missing_values(df, missing_df, strategy=missing_strategy)
    
    # Step 5: Detect outliers
    outlier_info = detect_outliers(df_clean, num_cols)
    
    # Step 6: Handle outliers
    df_clean = handle_outliers(df_clean, outlier_info, strategy=outlier_strategy)
    
    # Step 7: Encode categorical features
    valid_cat_cols = [col for col in cat_cols if col in df_clean.columns]
    df_encoded = encode_categorical_features(df_clean, valid_cat_cols)
    
    # Step 8: Normalize features
    valid_num_cols = [col for col in num_cols if col in df_encoded.columns]
    df_final = normalize_features(df_encoded, valid_num_cols, method=scaling_method)
    
    # Step 9: Visualize the processed data
    new_cat_cols = [col for col in df_final.columns if df_final[col].dtype == 'object']
    new_num_cols = [col for col in df_final.columns if df_final[col].dtype != 'object']
    visualize_data(df_final, new_cat_cols, new_num_cols)
    
    print("\n", "=" * 80)
    print("PREPROCESSING COMPLETE")
    print("=" * 80)
    print(f"Original data shape: {df.shape}")
    print(f"Processed data shape: {df_final.shape}")
    
    return df_final


if __name__ == "__main__":
    # Example usage
    # Replace 'your_dataset.csv' with your actual file path
    file_path = 'your_dataset.csv'
    
    # Run with default parameters
    # processed_df = preprocess_data(file_path)
    
    # Or run with custom parameters
    processed_df = preprocess_data(
        file_path,
        missing_strategy='auto',  # Options: 'auto', 'mean', 'median', 'mode', 'zero', 'drop_rows', 'drop_columns'
        outlier_strategy='auto',  # Options: 'auto', 'cap', 'winsorize', 'log_transform', 'remove'
        scaling_method='standard'  # Options: 'standard', 'minmax', 'robust'
    )
    
    # Save the processed data if needed
    if processed_df is not None:
        processed_df.to_csv('processed_data.csv', index=False)
        print("Processed data saved to 'processed_data.csv'")