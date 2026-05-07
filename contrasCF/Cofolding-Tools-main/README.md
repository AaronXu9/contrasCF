# Cofolding-Tools

This repository serves to reference the data and code of co-folding models used in the study "Do Deep Learning Models for Co-Folding Learn the Physics of Protein-Ligand Interactions?"

### AlphaFold3

AlphaFold3 is available via the repository here: [https://github.com/google-deepmind/alphafold3](https://github.com/google-deepmind/alphafold3). Version 3.0.0 was used in the study. Following their [installation instructions](https://github.com/google-deepmind/alphafold3/blob/main/docs/installation.md), you can run the test example by running the following Docker command:
```
docker run -it \
    --volume examples/af_input:/root/af_input \
    --volume examples/af_output:/root/af_output \
    --volume <MODEL_PARAMETERS_DIR>:/root/models \
    --volume <DATABASES_DIR>:/root/public_databases \
    --gpus all \
    alphafold3 \
    python run_alphafold.py \
    --json_path=/root/af_input/fold_input.json \
    --model_dir=/root/models \
    --output_dir=/root/af_output
```

### RosettaFold All-Atom

RFAA is available via the repository here: [https://github.com/baker-laboratory/RoseTTAFold-All-Atom](https://github.com/baker-laboratory/RoseTTAFold-All-Atom). Version 6c85140 was used in the study. Following their [installation instructions](https://github.com/baker-laboratory/RoseTTAFold-All-Atom?tab=readme-ov-file#set-up), you can run the test example by running the following commands:
```
cp examples/rfaa_input/example.yaml YOUR_CONFIG_PATH/.
python -m rf2aa.run_inference --config-name example
```

### Chai-1

Chai-1 is available via the repository here: [https://github.com/chaidiscovery/chai-lab](https://github.com/chaidiscovery/chai-lab). Version 0.5.0 was used in the study. Following their [installation instructions](https://github.com/chaidiscovery/chai-lab?tab=readme-ov-file#installation), you can run the test example by running the following command:
```
chai-lab fold --use-msa-server --use-templates-server examples/chai_input/example.fasta output_folder
```

### Boltz-1

The latest version of Boltz-1 is available via the repository here: [https://github.com/jwohlwend/boltz](https://github.com/jwohlwend/boltz). Version 0.4.1 was used in the study. Following their [installation instructions](https://github.com/jwohlwend/boltz?tab=readme-ov-file#installation), you can run the test example by running the following command:
```
boltz predict examples/boltz_input/example.yaml --use_msa_server
```

### Details
All models were run on NVIDIA GPU hardware with at least 12GB of RAM.
