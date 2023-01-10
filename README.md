# Entropy Maxima

"The Entropy Maximum principle has a corollary that is the Energy minimum principle" David Chandler
".. So there are multiple local maxima, just like there are multiple local minima." Noel Carrascal

# Installatio
From the directory where setup.py is loacated, type:

pip install -e . --user

# Third party software
To maintain all functionality in one package, Local Maxima is design to run third party software.
This is to avoid implementing functionalities that are already present in other programs. It also allows
for a smooth flow of procedures by genertating input files for the other programs and running them
internally. Local Maxima only requires CHARMM and Reduce programs to run. They most be already
installed for functionality such as prepr to work. This adds complexity to the installation for
a fully functional program, but it allows Local Maxima to do more while managing increased complexity
with simple command lines. In time, some of the functionality borrowed from third party software
will be coded into Local Maxima.

The templates to generate inputs for running third party software is in LocalMaxima/ThirdPartySoftware
which is in the .gitignore so that other people do not run into problems running it, to keep
Local Maxima purely a python software, and to allow people to customize their Local Maxima to their needs.
(To be updated soon. It is not yet decided how third party software will be included with the distribution)
# LocalMaxima
