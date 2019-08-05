library(CARNIVAL)

##
# Example with known inputs
inputFile = "AREG.txt"
measFile = "class_AREG.txt"
netFile = "Omnipath_signed_Uniprot_0615.txt"
Result_dir = "AREG_Result"

results = runCARNIVAL(solverPath = "~/Desktop/cbc", netFile = netFile, measFile = measFile, inputFile = inputFile, CARNIVAL_example = NULL, timelimit = 3600, mipGAP = 0, solver = "cbc", Result_dir = Result_dir)

##
# Example with un-known (inverse CARNIVAL -- inverseCR = TRUE) inputs and with PROGENy weights
measFile = "Ex3_measurement_APAP_TGG_24hrHighDose.txt"
netFile = "Omnipath_signed_Uniprot_0615.txt"
weightFile = "Ex3_weights_APAP_TGG_24hrHighDose.txt"
Result_dir = "APAP_Result"

results = runCARNIVAL(solverPath = "~/Desktop/cbc", netFile = netFile, measFile = measFile, weightFile = weightFile, CARNIVAL_example = NULL, timelimit = 3600, mipGAP = 0, Result_dir = Result_dir, solver = "cbc", inverseCR = TRUE)
