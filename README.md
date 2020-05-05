# EMBA
Software for modeling and simulation of our dynamic pathway regulation strategy for cell factories: The Extended metabolic TF-based biosensor antithetic controller (EMBA).

This software was used to obtain the dataset for the article “Extended metabolic biosensor design for dynamic pathway regulation of cell factories” Yadira Boadaa, Alejandro Vignoni, Jesús Picó, Pablo Carbonell, 2020, iScience.
All the data obtained from the implemented mathematical model has been publish as a Mendeley dataset and can be accesed and refereced as: Boada,588Yadira; Vignoni, Alejandro; Picó, Jesús; Carbonell, Pablo (2020), “Extended metabolic589biosensor design for dynamic pathway regulation of cell factories”, Mendeley Data, V1,590doi: 10.17632/hpxhkyvctb.1

This software includes 4 matlab files:

model.m implements the ODEs of the whole system with the Antithetic controller.

model_cI_repressor.m implements the ODEs of the whole system with the Direct controller (this is to perform the comparison between the direct controller and the proposed one with the Antithetic)

parameters.m is a function to define the structure p with the model marameters.

Main_antithetic_marker.m is the main execution script running secuential simulation of both models.

