
import numpy as np

TEMPLATE_FEATURES = {
    'template_aatype': np.float32,
    'template_all_atom_masks': np.float32,
    'template_all_atom_positions': np.float32,
    'template_domain_names': np.object,
    'template_e_value': np.float32,
    'template_neff': np.float32,
    'template_prob_true': np.float32,
    'template_release_date': np.object,
    'template_score': np.float32,
    'template_similarity': np.float32,
    'template_sequence': np.object,
    'template_sum_probs': np.float32,
    'template_confidence_scores': np.int64
}
