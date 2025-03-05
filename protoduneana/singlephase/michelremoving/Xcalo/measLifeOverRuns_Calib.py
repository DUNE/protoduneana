import os
from argparse import ArgumentParser as ap
#runArray=[5141, 5143, 5145, 5146, 5152, 5158, 5174, 5181, 5185, 5190, 5194, 5199, 5203, 5204, 5205, 5209, 5211, 5212, 5213, 5216,5219, 5225, 5235, 5240, 5244, 5249, 5250, 5254, 5257, 5258, 5259, 5260, 5261, 5267, 5276, 5282, 5283, 5284, 5287, 5290, 5293, 5298, 5301, 5303, 5304, 5308, 5311, 5313, 5315, 5338, 5341,5423, 5424, 5426, 5429, 5430, 5431, 5432, 5433, 5434, 5437, 5438, 5439, 5441, 5442, 5449, 5450, 5451, 5452, 5455, 5456, 5457, 5458, 5460, 5758, 5759, 5760, 5762, 5765, 5766, 5768, 5769, 5770, 5771, 5772, 5773,5774, 5775, 5776, 5777, 5778, 5779, 5780, 5783, 5784, 5785, 5786, 5788, 5791, 5792, 5794, 5796, 5797, 5809, 5810, 5814, 5815, 5816, 5817, 5818, 5819, 5824, 5825, 5826, 5827, 5831, 5833, 5834, 5835, 5836, 5837, 5838, 5839, 5840, 5841, 5842, 5843, 5844]
# From me just queing samweb

#runArray=[5834, 5835, 5838, 5839, 5840, 5841,5825, 5826, 5827, 5831, 5833, 5836, 5837,5219, 5225, 5235, 5240, 5244, 5308, 5311, 5315, 5338, 5387, 5423, 5424, 5426, 5455, 5456, 5457, 5458, 5460, 5809, 5810, 5814, 5816, 5817, 5842, 5843, 5844,5429, 5430, 5431, 5432, 5433, 5434, 5437, 5438, 5439, 5441, 5442, 5449, 5450, 5451, 5452, 5818, 5819, 5824,5777, 5778, 5779, 5780, 5783, 5784, 5785, 5786, 5788, 5791, 5792, 5794, 5796, 5797,5758, 5759, 5760, 5762, 5765, 5766, 5768, 5769, 5770, 5771, 5772, 5773, 5774, 5775, 5776,5141, 5143, 5145, 5152, 5158, 5174, 5181, 5185, 5190, 5194, 5199, 5203, 5204, 5815]

runArray=[5141, 5143, 5145, 5152, 5158, 5174, 5181, 5185, 5190, 5194, 5199, 5203, 5204, 5219, 5225, 5235, 5240, 5244, 5308, 5311, 5315, 5338, 5387, 5423, 5424, 5426, 5429, 5430, 5431, 5432, 5433, 5434, 5437, 5438, 5439, 5441, 5442, 5449, 5450, 5451, 5452, 5455, 5456, 5457, 5458, 5460, 5758, 5759, 5760, 5762, 5765, 5766, 5768, 5769, 5770, 5771, 5772, 5773, 5774, 5775, 5776, 5777, 5778, 5779, 5780, 5783, 5784, 5785, 5786, 5788, 5791, 5792, 5794, 5796, 5797, 5809, 5810, 5814, 5815, 5816, 5817, 5818, 5819, 5824, 5825, 5826, 5827, 5831, 5833, 5834, 5835, 5836, 5837, 5838, 5839, 5840, 5841, 5842, 5843, 5844]
# Official from Tingjun (unsorted and sorted)

#runArray=runArray[-10:]
print(len(runArray))


#runArray=[5770]
#runArray=[5844]
#for i in runArray:
    #os.system("samweb prestage-dataset --defname=protodune-sp_runset_"+str(i)+"_michelremoving_merged_v09_09_01_v0  --parallel=5")


for i in runArray:
    a=1
    
    os.system("python filelisting.py -r "+str(i))

    os.system("make_yz_correction file_list_"+str(i)+".txt 1")
    
    os.system("make_x_correction file_list_"+str(i)+".txt 1 1 1 ")
    os.system("make_x_correctionAlt file_list_"+str(i)+".txt 1 1 1 ")
    
     
    os.system("root -l -b -q 'plotXCaloLifetime.C(\""+str(i)+"\",\"All\",\"high\",\"second\")'")
    os.system("root -l -b -q 'plotXCaloLifetime.C(\""+str(i)+"\",\"All\",\"mid\",\"second\")'")
    os.system("root -l -b -q 'plotXCaloLifetime.C(\""+str(i)+"\",\"All\",\"low\",\"second\")'")

    os.system("root -l -b -q 'plotXCaloLifetime.C(\""+str(i)+"\",\"All\",\"high\",\"first\")'")
    os.system("root -l -b -q 'plotXCaloLifetime.C(\""+str(i)+"\",\"All\",\"mid\",\"first\")'")
    os.system("root -l -b -q 'plotXCaloLifetime.C(\""+str(i)+"\",\"All\",\"low\",\"first\")'")

    os.system("root -l -b -q 'plotXCaloLifetime.C(\""+str(i)+"\",\"All\",\"high\",\"third\")'")
    os.system("root -l -b -q 'plotXCaloLifetime.C(\""+str(i)+"\",\"All\",\"mid\",\"third\")'")
    os.system("root -l -b -q 'plotXCaloLifetime.C(\""+str(i)+"\",\"All\",\"low\",\"third\")'")
    

    os.system("root -l -b -q 'plotXCaloLifetime.C(\""+str(i)+"\",\"All\",\"high\",\"second\",0)'")
    os.system("root -l -b -q 'plotXCaloLifetime.C(\""+str(i)+"\",\"All\",\"mid\",\"second\",0)'")
    os.system("root -l -b -q 'plotXCaloLifetime.C(\""+str(i)+"\",\"All\",\"low\",\"second\",0)'")

    os.system("root -l -b -q 'plotXCaloLifetime.C(\""+str(i)+"\",\"All\",\"high\",\"first\",0)'")
    os.system("root -l -b -q 'plotXCaloLifetime.C(\""+str(i)+"\",\"All\",\"mid\",\"first\",0)'")
    os.system("root -l -b -q 'plotXCaloLifetime.C(\""+str(i)+"\",\"All\",\"low\",\"first\",0)'")

    os.system("root -l -b -q 'plotXCaloLifetime.C(\""+str(i)+"\",\"All\",\"high\",\"third\",0)'")
    os.system("root -l -b -q 'plotXCaloLifetime.C(\""+str(i)+"\",\"All\",\"mid\",\"third\",0)'")
    os.system("root -l -b -q 'plotXCaloLifetime.C(\""+str(i)+"\",\"All\",\"low\",\"third\",0)'")

    os.system("root -l -b -q 'plotXCaloLifetime.C(\""+str(i)+"\",\"All\",\"high\",\"second\",1)'")
    os.system("root -l -b -q 'plotXCaloLifetime.C(\""+str(i)+"\",\"All\",\"mid\",\"second\",1)'")
    os.system("root -l -b -q 'plotXCaloLifetime.C(\""+str(i)+"\",\"All\",\"low\",\"second\",1)'")

    os.system("root -l -b -q 'plotXCaloLifetime.C(\""+str(i)+"\",\"All\",\"high\",\"first\",1)'")
    os.system("root -l -b -q 'plotXCaloLifetime.C(\""+str(i)+"\",\"All\",\"mid\",\"first\",1)'")
    os.system("root -l -b -q 'plotXCaloLifetime.C(\""+str(i)+"\",\"All\",\"low\",\"first\",1)'")

    os.system("root -l -b -q 'plotXCaloLifetime.C(\""+str(i)+"\",\"All\",\"high\",\"third\",1)'")
    os.system("root -l -b -q 'plotXCaloLifetime.C(\""+str(i)+"\",\"All\",\"mid\",\"third\",1)'")
    os.system("root -l -b -q 'plotXCaloLifetime.C(\""+str(i)+"\",\"All\",\"low\",\"third\",1)'")
    
    
    os.system("root -l -b -q 'plotXCaloLifetimeAlt.C(\""+str(i)+"\",\"All\",\"high\",\"second\")'")
    os.system("root -l -b -q 'plotXCaloLifetimeAlt.C(\""+str(i)+"\",\"All\",\"mid\",\"second\")'")
    os.system("root -l -b -q 'plotXCaloLifetimeAlt.C(\""+str(i)+"\",\"All\",\"low\",\"second\")'")

    os.system("root -l -b -q 'plotXCaloLifetimeAlt.C(\""+str(i)+"\",\"All\",\"high\",\"first\")'")
    os.system("root -l -b -q 'plotXCaloLifetimeAlt.C(\""+str(i)+"\",\"All\",\"mid\",\"first\")'")
    os.system("root -l -b -q 'plotXCaloLifetimeAlt.C(\""+str(i)+"\",\"All\",\"low\",\"first\")'")

    os.system("root -l -b -q 'plotXCaloLifetimeAlt.C(\""+str(i)+"\",\"All\",\"high\",\"third\")'")
    os.system("root -l -b -q 'plotXCaloLifetimeAlt.C(\""+str(i)+"\",\"All\",\"mid\",\"third\")'")
    os.system("root -l -b -q 'plotXCaloLifetimeAlt.C(\""+str(i)+"\",\"All\",\"low\",\"third\")'")


    os.system("root -l -b -q 'plotXCaloLifetimeAlt.C(\""+str(i)+"\",\"All\",\"high\",\"second\",0)'")
    os.system("root -l -b -q 'plotXCaloLifetimeAlt.C(\""+str(i)+"\",\"All\",\"mid\",\"second\",0)'")
    os.system("root -l -b -q 'plotXCaloLifetimeAlt.C(\""+str(i)+"\",\"All\",\"low\",\"second\",0)'")

    os.system("root -l -b -q 'plotXCaloLifetimeAlt.C(\""+str(i)+"\",\"All\",\"high\",\"first\",0)'")
    os.system("root -l -b -q 'plotXCaloLifetimeAlt.C(\""+str(i)+"\",\"All\",\"mid\",\"first\",0)'")
    os.system("root -l -b -q 'plotXCaloLifetimeAlt.C(\""+str(i)+"\",\"All\",\"low\",\"first\",0)'")

    os.system("root -l -b -q 'plotXCaloLifetimeAlt.C(\""+str(i)+"\",\"All\",\"high\",\"third\",0)'")
    os.system("root -l -b -q 'plotXCaloLifetimeAlt.C(\""+str(i)+"\",\"All\",\"mid\",\"third\",0)'")
    os.system("root -l -b -q 'plotXCaloLifetimeAlt.C(\""+str(i)+"\",\"All\",\"low\",\"third\",0)'")

    os.system("root -l -b -q 'plotXCaloLifetimeAlt.C(\""+str(i)+"\",\"All\",\"high\",\"second\",1)'")
    os.system("root -l -b -q 'plotXCaloLifetimeAlt.C(\""+str(i)+"\",\"All\",\"mid\",\"second\",1)'")
    os.system("root -l -b -q 'plotXCaloLifetimeAlt.C(\""+str(i)+"\",\"All\",\"low\",\"second\",1)'")

    os.system("root -l -b -q 'plotXCaloLifetimeAlt.C(\""+str(i)+"\",\"All\",\"high\",\"first\",1)'")
    os.system("root -l -b -q 'plotXCaloLifetimeAlt.C(\""+str(i)+"\",\"All\",\"mid\",\"first\",1)'")
    os.system("root -l -b -q 'plotXCaloLifetimeAlt.C(\""+str(i)+"\",\"All\",\"low\",\"first\",1)'")

    os.system("root -l -b -q 'plotXCaloLifetimeAlt.C(\""+str(i)+"\",\"All\",\"high\",\"third\",1)'")
    os.system("root -l -b -q 'plotXCaloLifetimeAlt.C(\""+str(i)+"\",\"All\",\"mid\",\"third\",1)'")
    os.system("root -l -b -q 'plotXCaloLifetimeAlt.C(\""+str(i)+"\",\"All\",\"low\",\"third\",1)'")
    




















