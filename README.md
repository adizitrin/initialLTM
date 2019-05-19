# initialLTM
Matlab code to construct a LTM model (without multiple images as input).
The library includes: (1) The Matlab code initial_LTM_model.m -- the code to build the model. (2) An example, input galaxy-cluster member catalog. This is the same catalog as used to construct the preliminary model for MACS0329 in Carrasco, Zitrin, Seidel & Bartelmann 2019. To run the function, copy both files to your computer to the same location. In a Matlab session, run the code with the given input catalog -- by typing: initial_LTM_model('MemberCatM0329.txt')

%Written originally by Adi Zitrin and Tom Broadhurst (Zitrin et al., 2009, MNRAS, 396, 1985)
%Updated significantly since. Last updated: May 2019.
%Suggested and implemented improvements by other contributors is highly
%appreciated especially Matthias Bartelmann, Gregor Seidel, Irene Sendra, Keiichi Umetsu
%and much advice from CLASH, Hubble Frontier Fields, and RELICS collaborators.
