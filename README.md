# Kick-Fukui Algorithm

![alt text](https://github.com/HumanOsv/Logos/blob/master/Gihub.jpg)

A novel program for the search of global minimum structures of atomic clusters and molecules in the gas phase, using as a Coulomb interaction between Fukui functions-based method, is introduced in this work. The proposed strategy involves different steps. First, the calculation of the Fukui functions of the different fragments identifying them as electron- acceptor or donor species, according to their chemical potential value. Second, a population of individuals, created by the interaction between molecular fragments, is generated by the stochastic distribution of these fragments in size selected spatial box. Thirstily, the best structures are selected according to their Coulombic interaction between Fukui functions values.  Subsequently, the best candidates are relaxed to the nearest local minimum using a gradient method, i.e., ab-initio or DFT methods.


# Getting Started

**1)	Step Zero**

Before starting the installation, it is important to know that Kick-Fukui uses as input file *(.frag)* the information of *Condensed values of the Fukui functions called attractors*, it needs an assisting program to calculate *Topological Analysis of the Fukui Function* like TAFF, Multiwfn, Dgrid, TopChem2, AIMALL, TopMod, AIM-UC, AIM2000, etc.

File format with extension **.frag**

    • First column         = Atom symbol for attractor (X) and element (example Si).

    • Second-Fourth column = Cartesian coordinates (Å) x-coord; y-coord; z-coord.
  
    • Fifth column         = Condensed values of the Fukui functions (attractors).

**Note: For more information of Condensed values of the Fukui functions, https://pubs.acs.org/doi/10.1021/ct100022w**

       X	0.02053684	-1.29029342	-1.36079835	0.1900536003
       X	0.02053684	 1.31835992	-1.36079835	0.1926297737
       X	0.02053684	-1.81202408	 0.20439365	0.1526948452
      Si	0.00000000	 1.57169558	-0.53352285
      Si	0.00000000	 0.00000000	 1.06704517
      Si	0.00000000	-1.57169558	-0.53352285

**1. Softwares for Topology Analysis**

  •	TAFF (https://github.com/HumanOsv/TAFF)

  •	Multiwfn (http://sobereva.com/multiwfn/)

  •	Dgrid (http://www2.cpfs.mpg.de/~kohout/dgrid.html)

  •	TopChem2 (http://www.lct.jussieu.fr/pagesperso/pilme/topchempage.html)
  
  •	AIMALL (http://aim.tkgristmill.com/)
   
  •	TopMod (http://www.lct.jussieu.fr/pagesperso/silvi/topmod_english.html)

  •	AIM-UC (https://facyt-quimicomp.neocities.org/aim_uc/manual/manual_en.html)

  •	AIM2000 (http://www.aim2000.de)
   
   
**2. Installing Perl environment.**

Once the programs for calculating Attractors using the Topological Analysis of the Fukui Function is functioning properly, the Perl environment needs to be installed as well. Perl is a highly capable, feature-rich programming language that runs on many platforms from portable to mainframes.
It can be installed from:
- https://www.perl.org/get.html

There are some additional libraries and softwares that must also be installed to allow Kick-Fukui to work:

-Install CPAN modules (http://www.cpan.org/modules/INSTALL.html or https://egoleo.wordpress.com/2008/05/19/how-to-install-perl-modules-through-cpan-on-ubuntu-hardy-server/)

    user$ sudo cpan Parallel::ForkManager
      
    user$ sudo cpan Math::Matrix

**2)	Downloading and Installing Kick-Fukui**

Kick-Fukui can be directly downloaded as a zip file from the page:

-https://github.com/HumanOsv/Kick-Fukui

Alternatively, it can be downloaded using the Git tools using the following command:

    user$ git clone https://github.com/HumanOsv/Kick-Fukui.git

    user$ cd ./Kick-Fukui

**Note: before downloading using Git tools, make sure to be in your final installation path.**

We recommend to install using Git tools to update future Kick-Fukui software easily. To update the program, use the following command:

    user$ git pull master
	
Alternatively, Kick-Fukui could be installed as follows: choose a final installation path, and then extract the ZIP file (containing the software). Provide all the basic permissions for use and, optionally, set Kick-Fukui_Algorithm.pl file as a system call.

**3)	Running Kick-Fukui**

To run Kick-Fukui the following files are necessary in the working directory:

    • f+_inputfile.frag         : File with the nucleophilic Fukui function.
    
    • f-_inputfile.frag         : File with the electrophilic Fukui function.
    
    • Config.in                 : The Kick-Fukui input file, see below for more information.
    
    • DupGrigoryanSpringborg.pl : The executable file for duplicate molecular fragments.

    • Kick-Fukui_Algorithm.pl   : The executable file for structure prediction.
    
**Note: Kick-Fukui_Algorithm.pl can be called from another path if correctly set**

Now, use the following commands to execute this program:

    user$  perl Kick-Fukui_Algorithm.pl Config.in > Output.log

Alternatively, the user can set Kick-Fukui to run in the background using one of the following methods:

    user$ nohup perl Kick-Fukui_Algorithm.pl Config.in > Output.log
    user$ setsid perl Kick-Fukui_Algorithm.pl Config.in > Output.log

**4)	Input File**

The main input file, known as Config.in, contains all the necessary parameters for a correct calculation. Each variable is explained below.

The Initial Population size (1000N, N = Atoms number).

    initial_species = 5000

Number of final geometries.

    final_species = 50

**NOTE: Species to send for quantum calculations.**

Name of molecular species in the simulation.

*File with the nucleophilic Fukui function*

    Mol_1 = f+koop-Si3.frag

*File with the electrophilic Fukui function* 

    Mol_2 = f-koop-Si9.frag

Box size the sum of the sides from both molecules multiply for a factor (Default 1). The factor size or size of the box (in Angstroms) length, width, and height.

    box_size_factor = 0.8

Order Coulombic-Interaction Value (J) in a descending (YES) or scholastic (NO) manner. 

    energy_order = YES

Search for duplicate molecular species (YES/NO).

    duplicate_species = NO

Output xyz file name. 

    output_file_name = Si12

**General Note:** Respect the spaces of separation between the symbol "=".

    Correct : initial_species = 5000
    Wrong   : initial_species=5000

**5) Kick-Fukui outputs**

After a successful run of the program, several output files will be generated in your working directory.

    Outputfile_sort.xyz       : Final coordinates XYZ file format of each species sort by descending order value J.
    Outputfile_stochastic.xyz : Final coordinates XYZ file format of each species stochastic order value J.
    Output.log                : Print summary information Kick-Fukui.
    DuplicatesCoords_GS.xyz   : Duplicates Grigoryan-Springborg Algorithm.



