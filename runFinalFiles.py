import os

#os.system(r"g++ -g -o ControlPlotsSignalShapes ControlPlotsSignalShapes.cxx MakeRDF.cxx RDFObjects.cxx RDFevent.cxx -O2 $(root-config --libs --cflags) -std=c++17")
#os.system(r"g++ -g -o CutFlow CutFlow.cxx MakeRDF.cxx RDFObjects.cxx RDFevent.cxx -O2 $(root-config --libs --cflags) -std=c++17")
#os.system(r"g++ -g -o DataBackgroundComparison DataBackgroundComparison.cxx MakeRDF.cxx RDFObjects.cxx RDFevent.cxx -O2 $(root-config --libs --cflags) -std=c++17")
#os.system(r"g++ -g -o Categorization Categorization.cxx MakeRDF.cxx RDFObjects.cxx RDFevent.cxx -O2 $(root-config --libs --cflags) -std=c++17")
#os.system(r"g++ -g -o Figs_34_52 Figs_34_52.cxx MakeRDF.cxx RDFObjects.cxx RDFevent.cxx -O2 $(root-config --libs --cflags) -std=c++17")
os.system(r"g++ -g -o Section7_TablesPlots Section7_TablesPlots.cxx MakeRDF.cxx RDFObjects.cxx RDFevent.cxx -O2 $(root-config --libs --cflags) -std=c++17")

print("Success!")

#os.system(r"./ControlPlotsSignalShapes")
#os.system(r"./CutFlow")
#os.system(r"./DataBackgroundComparison")
#os.system(r"./Categorization")
#os.system(r"./Figs_34_52")
os.system(r"./Section7_TablesPlots")
