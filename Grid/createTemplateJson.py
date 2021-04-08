import json

# Options:
# id (int): identify the key we are submiting.
# enable (bool): to activate/deactivate the crab jobs submission for a given dataset
# localpath: user local area for submission
# eospath: user eos area for submission
# sample: name of the dataset (for mode data_analysis or mc_analysis). For mc_private_production is the name of the datset will be created.
# mode: data_analysis (for Lumimask file), mc_analysis (for eventbased submission) and mc_private_production (for MC production)
# lumimask: only used for data_analysis. Irrelevant for other modes. 
# config: CMSSW config file.
# parameters: parameters of the CMSSW config file.
#      i.e: 'parameters':["--Mode=Muons", "--Era=B"] or 'parameters':[""]
# output: name of the output file.
# events: number of events per job (for mc_analysis) or number of events per NJOBS (for mc_private_production). NJOBS is hardcoded in gridtool library.
# site: 'T2_BR_SPRACE', T2_IT_Pisa, T2_CH_CERNBOX, T2_US_Wisconsin, T2_BR_UERJ

data = {}
data['datasets'] = []

data['datasets'].append({
    'id': 0,
    'enable': 0,
    'localpath': "PPSFramework/Grid/",
    'eospath': "store/user/",
    'sample': "PYTHIA8-SD-TOP-GEN",
    'mode': "mc_private_production",
    'lumimask': "",
    'config': "config.py",
    'parameters':[""],
    'output': "output.root",
    'events': "10",
    'site': "T2_US_Wisconsin",
})

print(json.dumps(data, indent=4))

with open("samples.json", "w") as write_file:
    json.dump(data, write_file, indent = 4)

