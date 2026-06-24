import pandas as pd
import numpy as np
import argparse
import torch
from sentence_transformers import SentenceTransformer
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import pairwise_distances

def main():
    # Set up argument parser for command-line execution
    parser = argparse.ArgumentParser(description="Cluster scientific abstracts and find distinct topics.")
    parser.add_argument("--input", required=True, help="Path to the input TSV file (columns must be ID, title, abstract, without a header row)")
    parser.add_argument("--output", required=True, help="Path to the output TSV file with cluster annotations")
    parser.add_argument("--threshold", type=float, default=0.15, help="Cosine distance threshold (0.0 to 1.0). Lower = more topics, Higher = fewer topics.")
    parser.add_argument("--model", default="all-MiniLM-L6-v2", help="Name or path of the Hugging Face sentence-transformer model to use (e.g., 'all-MiniLM-L6-v2', 'pritlab/sci-bert-nli')")
    args = parser.parse_args()

    print(f"Loading data from {args.input}...")
    try:
        # Load the TSV file WITHOUT a header row and assign column names manually
        df = pd.read_csv(args.input, sep='\t', header=None, names=['ID', 'title', 'abstract'])
    except Exception as e:
        print(f"Error reading file: {e}")
        return

    # Combine title and abstract for better semantic representation
    df['text'] = df['title'].fillna('') + ". " + df['abstract'].fillna('')
    texts = df['text'].tolist()
    n_abstracts = len(texts)
    print(f"Loaded {n_abstracts} abstracts.")

    # --- DYNAMIC GPU/CPU HANDLING ---
    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"Hardware detected: {device.upper()}")
    batch_size = 64 if device == "cuda" else 32

    print(f"Generating embeddings using model: {args.model}...")
    # Load the user-specified model
    model = SentenceTransformer(args.model, device=device)
    embeddings = model.encode(
        texts, 
        batch_size=batch_size, 
        show_progress_bar=True,
        device=device
    )

    # DIAGNOSTIC: Calculate overall pairwise distances to help choose a threshold
    print("Calculating overall similarity metrics...")
    overall_distances = pairwise_distances(embeddings, metric='cosine')
    # Remove the diagonal (distance to itself is 0)
    mean_dist = np.mean(overall_distances[overall_distances > 0])
    max_dist = np.max(overall_distances)
    print(f"--- DIAGNOSTIC INFO ---")
    print(f"Average distance between any two abstracts: {mean_dist:.4f}")
    print(f"Maximum distance between two abstracts: {max_dist:.4f}")
    print(f"If clustering merges everything into 1 cluster, LOWER your threshold below {mean_dist:.4f}.")
    print(f"-----------------------")

    print(f"Clustering with Agglomerative Clustering (threshold={args.threshold})...")
    
    # Agglomerative Clustering directly on embeddings
    clusterer = AgglomerativeClustering(
        n_clusters=None,
        metric='cosine',         
        linkage='average',       
        distance_threshold=args.threshold
    )
    cluster_labels = clusterer.fit_predict(embeddings)

    # Add cluster assignments to the dataframe
    df['cluster_id'] = cluster_labels
    df['is_representative'] = False
    df['distance_to_center'] = 0.0

    unique_clusters = df['cluster_id'].unique()
    print(f"Found {len(unique_clusters)} distinct natural topics/experimental conditions. Identifying representatives...")

    # Extract Representative IDs (Medoids) and calculate distances
    for cluster_id in unique_clusters:
        cluster_indices = df[df['cluster_id'] == cluster_id].index.tolist()
        
        # If a cluster has only 1 article, it is its own representative, distance is 0.0
        if len(cluster_indices) == 1:
            df.loc[cluster_indices[0], 'is_representative'] = True
            continue
            
        # If >1 article, calculate pairwise distances within the cluster
        cluster_vectors = embeddings[cluster_indices]
        distances = pairwise_distances(cluster_vectors, metric='cosine')
        
        # The medoid is the article with the smallest average distance to all others
        mean_distances = np.mean(distances, axis=1)
        medoid_idx_in_cluster = np.argmin(mean_distances)
        
        # Annotate the distance_to_center for all items in this cluster
        for i, idx in enumerate(cluster_indices):
            df.loc[idx, 'distance_to_center'] = mean_distances[i]
        
        # Mark the representative as True
        original_df_index = cluster_indices[medoid_idx_in_cluster]
        df.loc[original_df_index, 'is_representative'] = True

    # Drop the temporary 'text' column used for processing
    df = df.drop(columns=['text'])

    print(f"Writing annotated data to {args.output}...")
    # Save the dataframe to a TSV file. 
    # We use header=False to match the input format, but you can change it to True if you prefer headers in the output.
    df.to_csv(args.output, sep='\t', index=False, header=True)
            
    print("Done.")

if __name__ == "__main__":
    main()
