#!/usr/bin/env python3
"""
TES Finder: Identify alternative transcription end sites (polyA sites)
using Gaussian Mixture Models
"""

import logging
import os

# One BLAS thread per process; parallelism is over transcripts.
for _variable in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS'):
    os.environ.setdefault(_variable, '1')

import numpy as np
from typing import List, Tuple, Dict
from collections import Counter
from sklearn.mixture import GaussianMixture
from sklearn.metrics import silhouette_score
from dataclasses import dataclass
from joblib import Parallel, delayed
from joblib.externals.loky import get_reusable_executor
from tqdm import tqdm

logger = logging.getLogger(__name__)



@dataclass
class APASite:
    """Represents an alternative polyadenylation site"""
    position: int       # 1-based genomic position
    read_count: int     # Number of supporting reads


@dataclass
class TranscriptAPA:
    """APA information for a transcript"""
    site: List[int]          # List of positions (1-based)
    count: List[int]         # List of read counts per site
    abundance: List[float]   # List of relative abundances (count/total)


class TESAnalyzer:
    """Analyze termination sites using Gaussian Mixture Models"""
    
    def __init__(self,
                 min_cluster_size: int = 10,
                 max_k: int = 5,
                 min_distance: int = 50,
                 min_relative_dominance: float = 0.3,
                 min_sharpness: float = 0.5,
                 max_subsample: int = 50000,
                 random_seed: int = 42):
        """
        Initialize TES analyzer
        
        Args:
            min_cluster_size: Minimum reads in a TES cluster
            max_k: Maximum number of GMM components to test
            min_distance: Minimum distance between TES peaks (bp)
            min_relative_dominance: Minimum relative size for alternative TES
            min_sharpness: Minimum sharpness score for peaks
            max_subsample: Maximum reads used for GMM fitting and k selection;
                           all reads are still used for final cluster assignment
            random_seed: Random seed for reproducibility
        """
        self.min_cluster_size = min_cluster_size
        self.max_k = max_k
        self.min_distance = min_distance
        self.min_relative_dominance = min_relative_dominance
        self.min_sharpness = min_sharpness
        self.max_subsample = max_subsample
        self.random_seed = random_seed
    
    def find_tes_peaks(self, read_end_sites: np.ndarray) -> List[Tuple[int, int]]:
        """
        Find TES peaks from read end sites using all reads
        
        Args:
            read_end_sites: Array of 1-based read end positions
            
        Returns:
            List of (position, read_count) tuples for significant peaks
        """
        read_end_sites = self._validate_read_end_sites(read_end_sites)
        if len(read_end_sites) < self.min_cluster_size:
            return []
        
        # Sort read ends.
        read_end_sites = np.sort(read_end_sites)

        # Subsample for GMM fitting / silhouette scoring when n is very large.
        # All reads are still used for final cluster assignment and counting.
        if len(read_end_sites) > self.max_subsample:
            rng = np.random.default_rng(self.random_seed)
            idx = np.sort(rng.choice(len(read_end_sites), size=self.max_subsample, replace=False))
            sub_sites = read_end_sites[idx]
        else:
            sub_sites = read_end_sites

        return self._find_peaks_single(read_end_sites, sub_sites)
    
    def _find_peaks_single(self, read_end_sites: np.ndarray, sub_sites: np.ndarray) -> List[Tuple[int, int]]:
        """
        Single-pass GMM peak detection.

        Args:
            read_end_sites: Full sorted array of read end positions (used for
                            final cluster assignment, counts, and sharpness).
            sub_sites: Subsample of read_end_sites (or identical to it when
                       n <= max_subsample) used for GMM fitting and silhouette
                       scoring to keep runtime bounded.
        """
        # Find optimal number of clusters using the (possibly subsampled) data.
        # Returns the already-fitted best GMM to avoid a second fit.
        k_optimal, best_gmm = self._find_optimal_k(sub_sites)

        if k_optimal == 1:
            # Single cluster - return mode from full data
            # Handle ties deterministically by choosing smallest position
            position_counts = Counter(read_end_sites)
            max_count = max(position_counts.values())
            mode_pos = min([pos for pos, cnt in position_counts.items() if cnt == max_count])
            return [(mode_pos, len(read_end_sites))]

        # Reuse the GMM fitted during k selection; predict on all reads for accurate counts
        X_full = read_end_sites.reshape(-1, 1)
        labels = best_gmm.predict(X_full)

        # Extract significant clusters using full read set for true counts
        return self._extract_peaks(read_end_sites, labels)
    
    @staticmethod
    def _validate_read_end_sites(read_end_sites: np.ndarray) -> np.ndarray:
        """Validate the read-end array: one-dimensional, numeric and finite."""
        sites = np.asarray(read_end_sites)
        if (sites.ndim != 1 or not np.issubdtype(sites.dtype, np.number)
                or np.iscomplexobj(sites) or not np.all(np.isfinite(sites))):
            raise ValueError('Read-end positions must be a one-dimensional array of finite real numbers')
        return sites

    def _find_optimal_k(self, read_end_sites: np.ndarray) -> Tuple[int, object]:
        """Find optimal number of clusters using silhouette score.

        Returns:
            (k_optimal, best_gmm): the best k and its already-fitted GMM,
            so the caller can call predict() without refitting.
            Returns (1, None) when no multi-component candidate is feasible or
            every converged model predicts a single group; the caller then
            reports the modal site. Raises if every candidate fails numerically
            or does not converge.
        """
        read_end_sites = self._validate_read_end_sites(read_end_sites)
        max_sil = -np.inf
        k_optimal = 1
        best_gmm = None
        n_samples = len(read_end_sites)
        if n_samples == 0:
            raise ValueError('Cannot select a GMM for an empty read-end distribution')
        candidates = range(2, min(self.max_k + 1, n_samples))
        if not candidates or np.unique(read_end_sites).size == 1:
            return 1, None
        X = read_end_sites.reshape(-1, 1)
        collapsed_models = 0
        failures = []

        for k in candidates:
            try:
                gmm = GaussianMixture(n_components=k, random_state=self.random_seed,
                                      max_iter=200, n_init=3)
                gmm.fit(X)
                if not gmm.converged_:
                    reason = 'did not converge'
                    failures.append(f'K={k}: {reason}')
                    logger.warning('Skipping GMM K=%d for %d read ends: %s', k, n_samples, reason)
                    continue
                labels = gmm.predict(X)
                n_labels = len(np.unique(labels))
                if n_labels < 2:
                    collapsed_models += 1
                    continue
                if n_labels >= n_samples:
                    raise ValueError('Silhouette requires fewer predicted groups than observations')

                sil_score = silhouette_score(X, labels, metric='euclidean')
                if not np.isfinite(sil_score):
                    raise ValueError('Non-finite silhouette score')
                if sil_score > max_sil:
                    max_sil = sil_score
                    k_optimal = k
                    best_gmm = gmm
            except (ValueError, np.linalg.LinAlgError, FloatingPointError) as error:
                failures.append(f'K={k}: {type(error).__name__}: {error}')
                logger.warning('Skipping GMM K=%d for %d read ends: %s: %s',
                               k, n_samples, type(error).__name__, error)
                continue

        if best_gmm is None and not collapsed_models:
            raise RuntimeError(
                f'No evaluable GMM candidate for {n_samples} read ends; '
                'all candidates failed or did not converge. ' + '; '.join(failures)
            )
        if best_gmm is None:
            logger.info('Using modal fallback: %d converged candidate(s) predicted one group; '
                        '%d other candidate(s) failed', collapsed_models, len(failures))

        return k_optimal, best_gmm
    
    def _is_peak_sharp(self, position: int, read_end_sites: np.ndarray) -> bool:
        """Check if peak is sharp using IQR method"""
        # Get reads near this peak
        window_reads = read_end_sites[
            (read_end_sites >= position - self.min_distance / 2) &
            (read_end_sites <= position + self.min_distance / 2)
        ]
        
        if len(window_reads) < self.min_cluster_size:
            return False
        
        # Calculate sharpness using IQR
        q75, q25 = np.percentile(window_reads, [75, 25])
        iqr = q75 - q25
        
        sharpness = 1.0 - (iqr / self.min_distance)
        return bool(sharpness >= self.min_sharpness)
    
    def _extract_peaks(self, read_end_sites: np.ndarray, labels: np.ndarray) -> List[Tuple[int, int]]:
        """Extract significant peaks from clustered data"""
        peaks = []
        
        # Count reads per cluster
        cluster_counts = Counter(labels)
        
        # Keep components with enough reads; order by size, then component id.
        valid_clusters = sorted(
            [(k, v) for k, v in cluster_counts.items() if v >= self.min_cluster_size],
            key=lambda x: (x[1], x[0]), reverse=True
        )
        
        if not valid_clusters:
            return []
        
        dominant_size = valid_clusters[0][1]
        used_positions = set()
        
        for cluster_id, cluster_size in valid_clusters:
            # Get positions for this cluster
            cluster_positions = read_end_sites[labels == cluster_id]
            
            # Mode = most frequent position; ties resolve to the smallest coordinate.
            position_counts = Counter(cluster_positions)
            max_count = max(position_counts.values())
            mode_position = min([pos for pos, cnt in position_counts.items() if cnt == max_count])
            
            # Check sharpness
            if not self._is_peak_sharp(mode_position, read_end_sites):
                continue
            
            # Check relative dominance
            if cluster_size < self.min_relative_dominance * dominant_size:
                continue
            
            # Check distance from other peaks
            too_close = any(abs(mode_position - pos) < self.min_distance
                           for pos in used_positions)
            
            if not too_close:
                peaks.append((mode_position, cluster_size))
                used_positions.add(mode_position)
        
        return sorted(peaks, key=lambda x: x[0])


class TESFinder:
    """Parallel TES Finder for identifying alternative polyA sites"""
    
    def __init__(self,
                 min_reads: int = 10,
                 min_cluster_size: int = 10,
                 max_k: int = 5,
                 min_distance: int = 50,
                 min_relative_dominance: float = 0.3,
                 min_sharpness: float = 0.5,
                 max_subsample: int = 50_000,
                 n_jobs: int = -1,
                 random_seed: int = 42):
        """
        Initialize TES Finder
        
        Args:
            min_reads: Minimum reads to analyze a transcript
            min_cluster_size: Minimum reads in a TES cluster
            max_k: Maximum GMM components
            min_distance: Minimum distance between peaks (bp)
            min_relative_dominance: Minimum relative size for alternative site
            min_sharpness: Minimum peak sharpness
            max_subsample: Maximum reads used for GMM fitting and k selection
                           (all reads still used for final assignment)
            n_jobs: Number of parallel jobs (-1 for all CPUs, 1 for serial/reproducible)
            random_seed: Random seed for reproducibility
        """
        self.min_reads = min_reads
        self.n_jobs = n_jobs
        self.random_seed = random_seed
        
        self.analyzer_params = {
            'min_cluster_size': min_cluster_size,
            'max_k': max_k,
            'min_distance': min_distance,
            'min_relative_dominance': min_relative_dominance,
            'min_sharpness': min_sharpness,
            'max_subsample': max_subsample,
            'random_seed': random_seed
        }
        
        logger.info(f"Initialized TES Finder with {n_jobs} workers (random_seed={random_seed})")
    
    def find_apa_sites(self,
                       transcript_reads: Dict[str, List[int]]) -> Dict[str, TranscriptAPA]:
        """
        Find alternative polyA sites for all transcripts
        
        Args:
            transcript_reads: Dict mapping transcript_id to list of read_end_positions
            
        Returns:
            Dictionary of transcript_id -> TranscriptAPA
        """
        logger.info(f"Finding APA sites for {len(transcript_reads)} transcripts")
        
        # Prepare data for parallel processing
        # Sort positions within each transcript for reproducibility
        valid_transcripts = [
            (tid, np.sort(np.array(positions)))
            for tid, positions in transcript_reads.items()
            if len(positions) >= self.min_reads
        ]
        
        logger.info(f"Analyzing {len(valid_transcripts)} transcripts with sufficient reads")
        
        # Results are deterministic with either backend and are collected in submission order.
        backend = 'loky' if self.n_jobs != 1 else 'sequential'
        if self.n_jobs == 1:
            logger.info("Running in sequential mode")
        
        # Progress is reported per completed transcript.
        results = list(tqdm(
            Parallel(n_jobs=self.n_jobs, backend=backend, return_as='generator')(
                delayed(self._process_transcript)(
                    transcript_id,
                    end_positions,
                    self.analyzer_params
                ) for transcript_id, end_positions in valid_transcripts
            ),
            total=len(valid_transcripts),
            desc="TES Analysis"
        ))

        # Release the worker pool.
        if self.n_jobs != 1:
            get_reusable_executor().shutdown(wait=True)

        # Collect results
        apa_results = {}
        n_apa = 0
        
        for transcript_id, apa_info in results:
            apa_results[transcript_id] = apa_info
            if len(apa_info.site) > 1:
                n_apa += 1
        
        logger.info(f"Found alternative TES in {n_apa} transcripts")
        
        return apa_results
    
    @staticmethod
    def _process_transcript(transcript_id: str,
                           end_positions: np.ndarray,
                           analyzer_params: Dict) -> Tuple[str, TranscriptAPA]:
        """Process a single transcript (for parallel execution)"""
        
        analyzer = TESAnalyzer(**analyzer_params)
        
        # Find TES peaks
        peaks = analyzer.find_tes_peaks(end_positions)
        
        # If no peaks found, use mode
        if not peaks:
            # Handle ties deterministically by choosing smallest position
            position_counts = Counter(end_positions)
            max_count = max(position_counts.values())
            mode_pos = min([pos for pos, cnt in position_counts.items() if cnt == max_count])
            
            return transcript_id, TranscriptAPA(
                site=[int(mode_pos)],
                count=[len(end_positions)],
                abundance=[1.0]
            )
        
        # Sort peaks by read count (descending), then by position (ascending) for determinism
        peaks = sorted(peaks, key=lambda x: (-x[1], x[0]))
        
        # Split into separate lists and calculate abundance
        sites = [pos for pos, _ in peaks]
        counts = [count for _, count in peaks]
        total_count = sum(counts)
        abundances = [count / total_count for count in counts]
        
        return transcript_id, TranscriptAPA(
            site=sites,
            count=counts,
            abundance=abundances
        )


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
