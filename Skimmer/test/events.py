import subprocess
import json

#command = "dasgoclient -query=\"summary dataset=/ZMM/Summer11-DESIGN42_V11_428_SLHC1-v1/GEN-SIM\""
command = "dasgoclient -query=\"summary dataset=/HXToyMC_NONE_PostTS2_850_120/dmf-HXToyMC-13TeV-miniAOD_2021-01-05_UTC12-27-38-b439bc17361eb503f053ea70fabf4631/USER instance=prod/phys03\""
result = subprocess.check_output(command, shell=True)
result = json.loads(result)

print result[0]['nevents']
