"""Functions for parsing various file formats."""

import collections
import string
from typing import Sequence, Tuple

import dataclasses

DeletionMatrix = Sequence[Sequence[int]]

def parse_fasta(fasta_string: str) -> Tuple[Sequence[str], Sequence[str]]:
  """Parses FASTA string and returns list of strings with amino-acid sequences.

  Arguments:
    fasta_string: The string contents of a FASTA file.

  Returns:
    A tuple of two lists:
    * A list of sequences.
    * A list of sequence descriptions taken from the comment lines. In the
      same order as the sequences.
  """
  sequences = []
  descriptions = []
  index = -1
  for line in fasta_string.splitlines():
    line = line.strip()
    if line.startswith('>'):
      index += 1
      descriptions.append(line[1:])  # Remove the '>' at the beginning.
      sequences.append('')
      continue
    elif not line:
      continue  # Skip blank lines.
    sequences[index] += line

  return sequences, descriptions


def parse_stockholm(
    stockholm_string: str) -> Tuple[Sequence[str], DeletionMatrix]:
  """Parses sequences and deletion matrix from stockholm format alignment.

  Args:
    stockholm_string: The string contents of a stockholm file. The first
      sequence in the file should be the query sequence.

  Returns:
    A tuple of:
      * A list of sequences that have been aligned to the query. These
        might contain duplicates.
      * The deletion matrix for the alignment as a list of lists. The element
        at `deletion_matrix[i][j]` is the number of residues deleted from
        the aligned sequence i at residue position j.
  """
  name_to_sequence = collections.OrderedDict()
  for line in stockholm_string.splitlines():
    line = line.strip()
    if not line or line.startswith(('#', '//')):
      continue
    name, sequence = line.split()
    if name not in name_to_sequence:
      name_to_sequence[name] = ''
    name_to_sequence[name] += sequence

  msa = []
  deletion_matrix = []

  query = ''
  keep_columns = []
  for seq_index, sequence in enumerate(name_to_sequence.values()):
    if seq_index == 0:
      # Gather the columns with gaps from the query
      query = sequence
      keep_columns = [i for i, res in enumerate(query) if res != '-']

    # Remove the columns with gaps in the query from all sequences.
    aligned_sequence = ''.join([sequence[c] for c in keep_columns])

    msa.append(aligned_sequence)

    # Count the number of deletions w.r.t. query.
    deletion_vec = []
    deletion_count = 0
    for seq_res, query_res in zip(sequence, query):
      if seq_res != '-' or query_res != '-':
        if query_res == '-':
          deletion_count += 1
        else:
          deletion_vec.append(deletion_count)
          deletion_count = 0
    deletion_matrix.append(deletion_vec)

  return msa, deletion_matrix


def parse_a3m(a3m_string: str) -> Tuple[Sequence[str], DeletionMatrix]:
  """Parses sequences and deletion matrix from a3m format alignment.

  Args:
    a3m_string: The string contents of a a3m file. The first sequence in the
      file should be the query sequence.

  Returns:
    A tuple of:
      * A list of sequences that have been aligned to the query. These
        might contain duplicates.
      * The deletion matrix for the alignment as a list of lists. The element
        at `deletion_matrix[i][j]` is the number of residues deleted from
        the aligned sequence i at residue position j.
  """
  sequences, _ = parse_fasta(a3m_string)
  deletion_matrix = []
  for msa_sequence in sequences:
    deletion_vec = []
    deletion_count = 0
    for j in msa_sequence:
      if j.islower():
        deletion_count += 1
      else:
        deletion_vec.append(deletion_count)
        deletion_count = 0
    deletion_matrix.append(deletion_vec)

  # Make the MSA matrix out of aligned (deletion-free) sequences.
  deletion_table = str.maketrans('', '', string.ascii_lowercase)
  aligned_sequences = [s.translate(deletion_table) for s in sequences]
  return aligned_sequences, deletion_matrix
