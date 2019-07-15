library(CARNIVAL) # load CARNIVAL library

solverPath="/Users/davidhenriques/Desktop/carnival_cbc/cbc"
netFile="Ex1_network_Toy.sif"
measFile="Ex1_measurements_Toy.txt"
inputFile="Ex1_inputs_Toy.txt"
weightFile=NULL
CARNIVAL_example=NULL
Result_dir="Results_CARNIVAL"
inverseCR=F
# parallelCR=F,
parallelIdx1=1
parallelIdx2=1
nodeID="uniprot"
UP2GS=T
DOTfig=T
Export_all=F
timelimit=30
mipGAP=0.05
poolrelGAP=0.0001
limitPop=500
poolCap=100
poolIntensity=4
poolReplace=2
alphaWeight=1
betaWeight=0.2 
solver = "cbc"

CARNIVAL_Result <- runCARNIVAL(solverPath = "/Users/davidhenriques/Desktop/carnival_cbc/cbc",Result_dir="Results_CARNIVAL_Ex1", netFile = netFile, 
                               measFile = measFile, inputFile = inputFile, solver = "cbc", CARNIVAL_example=NULL, UP2GS=F)
