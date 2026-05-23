"""Shared NJ tree construction utilities."""

from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor


def build_nj_tree(labels, dist_matrix):
    """
    Construct a Neighbor-Joining tree from a pairwise distance matrix.

    Args:
        labels     : list of str, one per leaf
        dist_matrix: 2D indexable (numpy array or list of lists), symmetric

    Returns:
        Bio.Phylo tree object
    """
    n = len(labels)
    lower = [[float(dist_matrix[i][j]) for j in range(i + 1)] for i in range(n)]
    dm = DistanceMatrix(labels, lower)
    return DistanceTreeConstructor().nj(dm)
