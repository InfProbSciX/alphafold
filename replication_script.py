
from alphafold.model import data
from alphafold.model import model
from alphafold.model import config
from alphafold.relax import relax
from alphafold.data.pipeline import DataPipeline

model_names = [f'model_{i+1}' for i in range(5)]
relaxed_pdbs = {}; plddts = {}; model_runners = {}

feature_dict = DataPipeline().process(data_path='data/6y4f/')

for model_name in model_names:
    model_config = config.model_config(model_name)
    model_config.data.eval.num_ensemble = 8 # casp14
    model_params = data.get_model_haiku_params(
        model_name=model_name, data_dir='data')
    model_runner = model.RunModel(model_config, model_params)
    model_runners[model_name] = model_runner

amber_relaxer = relax.AmberRelaxation(
    max_iterations=0,
    tolerance=2.39,
    stiffness=10.0,
    exclude_residues=[],
    max_outer_iterations=20)

for model_name, model_runner in model_runners.items():
    processed_feature_dict = model_runner.process_features(
        feature_dict, random_seed=42)

    prediction_result = model_runner.predict(processed_feature_dict)
    plddts[model_name] = np.mean(prediction_result['plddt'])

    unrelaxed_protein = protein.from_prediction(processed_feature_dict,
                                                prediction_result)

    relaxed_pdb_str, _, _ = amber_relaxer.process(prot=unrelaxed_protein)
    relaxed_pdbs[model_name] = relaxed_pdb_str

best_model_name = sorted(plddts.items(), key=lambda x: x[1], reverse=True)[0][0]
with open('data/6y4f/ranked_0.pdb', 'w') as f:
    f.write(relaxed_pdbs[model_name])

